#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "algos/statistics_float.h"
#include "algos/sorting.h"
#include "registration/registration.h"
#include "stacking.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

#define DEBUG_NORM

typedef struct {
	size_t Nij;
	float mij, mji;
	float sij, sji;
} overlap_stats_t;

static int compute_normalization(struct stacking_args *args);
static int compute_normalization_overlaps(struct stacking_args *args);

/* normalization: reading all images and making stats on their background level.
 * That's very long if not cached. */
int do_normalization(struct stacking_args *args) {
	if (args->normalize == NO_NORM) return ST_OK;

	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	int nb_frames = args->nb_images_to_stack;
	int nb_layers = args->seq->nb_layers;

	args->coeff.offset = malloc(nb_layers * nb_frames * sizeof(double));
	args->coeff.mul = malloc(nb_layers * nb_frames * sizeof(double));
	args->coeff.scale = malloc(nb_layers * nb_frames * sizeof(double));

	if (!args->coeff.offset || !args->coeff.mul || !args->coeff.scale) {
		PRINT_ALLOC_ERR;
		args->retval = ST_ALLOC_ERROR;
		return args->retval;
	}

	if (args->overlap_norm && compute_normalization_overlaps(args)) {
		args->retval = ST_GENERIC_ERROR;
		return args->retval;
	} 
	if (!args->overlap_norm && compute_normalization(args)) {
		args->retval = ST_GENERIC_ERROR;
		return args->retval;
	}

	if (args->seq->needs_saving)	// if we had to compute new stats
		writeseqfile(args->seq);

	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, _("Normalization computation time"));
	return ST_OK;
}

static void init_coeffs(struct stacking_args *args) {
	norm_coeff *coeff = &args->coeff;
	int nb_layers = args->seq->nb_layers;
	int nb_frames = args->nb_images_to_stack;
	for (int layer = 0; layer < nb_layers; ++layer) {
		coeff->poffset[layer] = coeff->offset + layer * nb_frames;
		coeff->pmul[layer] = coeff->mul + layer * nb_frames;
		coeff->pscale[layer] = coeff->scale + layer * nb_frames;

		for (int i = 0; i < args->nb_images_to_stack; ++i) {
			coeff->poffset[layer][i] = 0.0;
			coeff->pmul[layer][i] = 1.0;
			coeff->pscale[layer][i] = 1.0;
		}
	}
}

static int _compute_estimators_for_image(struct stacking_args *args, int i,
		threading_type threading, int image_thread_id) {
	imstats *stats[3] = { NULL };
	int nb_layers = args->seq->nb_layers;
	int retval = ST_OK;
	g_assert(nb_layers <= 3);
	g_assert(threading > 0);

	siril_debug_print("computing stats for image %d, %d threads, lite: %d\n", i, threading, args->lite_norm);
	retval = compute_all_channels_statistics_seqimage(args->seq, args->image_indices[i], NULL, (args->lite_norm) ? STATS_LITENORM : STATS_NORM, threading, image_thread_id, stats);

	for (int layer = 0; layer < args->seq->nb_layers; ++layer) {
		imstats *stat = stats[layer];
		if (!stat) {
			retval = ST_GENERIC_ERROR;
			continue;
		}
		switch (args->normalize) {
		default:
		/* In order to compute quick loc and  scale estimates
		 * if lite_norm is true
		 * location is replaced by median
		 * scale is replaced by 1.5*mad which is a good approximation of
		 * sqrt(bwmv) for a gaussian distribution (tested of large gaussian samples of 10M values)
		 * sqrt(bwmv) itself being a good enough fit for IKSS scale estimator
		 */
		case ADDITIVE_SCALING:
			args->coeff.pscale[layer][i] = (args->lite_norm) ? 1.5 * stat->mad : stat->scale;
			/* no break */
		case ADDITIVE:
			args->coeff.poffset[layer][i] = (args->lite_norm) ? stat->median : stat->location;
			break;
		case MULTIPLICATIVE_SCALING:
			args->coeff.pscale[layer][i] = (args->lite_norm) ? 1.5 * stat->mad : stat->scale;
			/* no break */
		case MULTIPLICATIVE:
			args->coeff.pmul[layer][i] = (args->lite_norm) ? stat->median : stat->location;
			break;
		}
		free_stats(stat);
	}
	return retval;
}

/* computes normalization factors based on the reference image from the
 * estimators of each image and channel. Overwrites the estimators */
static void compute_factors_from_estimators(struct stacking_args *args, int ref_index) {
	double **poffset = args->coeff.poffset, **pmul = args->coeff.pmul, **pscale = args->coeff.pscale;
	int nb_layers = args->seq->nb_layers;
	double offset0[3], mul0[3], scale0[3];

	for (int layer = 0; layer < nb_layers; ++layer) {
		offset0[layer] = poffset[layer][ref_index];
		mul0[layer] = pmul[layer][ref_index];
		scale0[layer] = pscale[layer][ref_index];
	}
	int reflayer;
	for (int layer = 0; layer < nb_layers; ++layer) {
		reflayer = (args->equalizeRGB) ? args->reglayer : layer;
		for (int i = 0; i < args->nb_images_to_stack; ++i) {
			switch (args->normalize) {
				default:
				case ADDITIVE_SCALING:
					pscale[layer][i] = (pscale[layer][i] == 0) ? 1 : scale0[reflayer] / pscale[layer][i];
					/* no break */
				case ADDITIVE:
					poffset[layer][i] = pscale[layer][i] * poffset[layer][i] - offset0[reflayer];
					break;
				case MULTIPLICATIVE_SCALING:
					pscale[layer][i] = (pscale[layer][i] == 0) ? 1 : scale0[reflayer]  / pscale[layer][i];
					/* no break */
				case MULTIPLICATIVE:
					pmul[layer][i] = (pmul[layer][i] == 0) ? 1 : mul0[reflayer] / pmul[layer][i];
					break;
			}
#ifdef DEBUG_NORM
			siril_debug_print("%2d %d %.5f %5f %.5f\n", args->image_indices[i] + 1, layer, poffset[layer][i], pscale[layer][i], pmul[layer][i]);
#endif
		}
	}
}

static int normalization_get_max_number_of_threads(sequence *seq) {
	int max_memory_MB = get_max_memory_in_MB();
	/* The normalization memory consumption, n is image size and m channel size.
	 * It uses IKSS computation in stats, which can be done only on float data.
	 * IKSS requires computing the MAD, which requires its own copy of the data.
	 * The stats are computed successively on each channel.
	 * For DATA_USHORT, we have: the image ON, rewrite of the channel without
	 * zeros O(m), a copy to float for IKSS O(2m), a copy of that for the MAD in
	 * IKSS O(2m).
	 * For DATA_FLOAT, we have: the image ON, rewrite without zeros O(m),
	 * used directly for IKSS and a copy for MAD O(m).
	 */
	guint64 memory_per_image = (guint64) seq->rx * seq->ry;
	if (get_data_type(seq->bitpix) == DATA_FLOAT)
		memory_per_image *= (seq->nb_layers + 2) * sizeof(float);
	else memory_per_image *= (seq->nb_layers + 1) * sizeof(WORD) + 2 * sizeof(float);
	unsigned int memory_per_image_MB = memory_per_image / BYTES_IN_A_MB;
	if (memory_per_image_MB == 0)
		memory_per_image_MB = 1;

	fprintf(stdout, "Memory per image: %u MB. Max memory: %d MB\n", memory_per_image_MB, max_memory_MB);

	if (memory_per_image_MB > max_memory_MB) {
		siril_log_color_message(_("Your system does not have enough memory to normalize images for stacking operation (%d MB free for %d MB required)\n"), "red", max_memory_MB, memory_per_image_MB);
		return 0;
	}

	int nb_threads = max_memory_MB / memory_per_image_MB;
	if (nb_threads > com.max_thread)
		nb_threads = com.max_thread;
	siril_log_message(_("With the current memory and thread (%d) limits, up to %d thread(s) can be used for sequence normalization\n"), com.max_thread, nb_threads);
	return nb_threads;
}

static int compute_normalization(struct stacking_args *args) {
	int ref_image_filtred_idx = -1, retval = 0, cur_nb = 0;
	int nb_layers = args->seq->nb_layers;

	init_coeffs(args);
	if (args->force_norm) { /* We empty the cache if needed (force to recompute) */
		for (int layer = 0; layer < nb_layers; layer++) {
			clear_stats(args->seq, layer);
			args->seq->needs_saving = TRUE;
		}
	}

	if (args->normalize == NO_NORM)	// should never happen here
		return 0;

	char *tmpmsg = siril_log_message(_("Computing normalization...\n"));
	tmpmsg[strlen(tmpmsg) - 1] = '\0';
	set_progress_bar_data(tmpmsg, PROGRESS_RESET);

	// first, find the index of the ref image in the filtered image list
	ref_image_filtred_idx = find_refimage_in_indices(args->image_indices,
			args->nb_images_to_stack, args->ref_image);
	if (ref_image_filtred_idx == -1) {
		siril_log_color_message(_("The reference image is not in the selected set of images. "
				"Please choose another reference image.\n"), "red");
		return ST_GENERIC_ERROR;
	}

	const char *error_msg = (_("Normalization failed."));

	// check memory first
	int nb_threads = normalization_get_max_number_of_threads(args->seq);
	if (nb_threads <= 0) {
		set_progress_bar_data(error_msg, PROGRESS_NONE);
		return ST_GENERIC_ERROR;
	}
	if (nb_threads > args->nb_images_to_stack)
		nb_threads = args->nb_images_to_stack;
	int *threads_per_thread = compute_thread_distribution(nb_threads, com.max_thread);

	set_progress_bar_data(NULL, 1.0 / (double)args->nb_images_to_stack);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nb_threads) schedule(guided) \
	if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ) && fits_is_reentrant()))
#endif

	for (int i = 0; i < args->nb_images_to_stack; ++i) {
		if (!retval) {
			if (!get_thread_run()) {
				retval = 1;
				continue;
			}
			int thread_id = -1;
			threading_type threads = SINGLE_THREADED;
#ifdef _OPENMP
			thread_id = omp_get_thread_num();
			threads = threads_per_thread[thread_id];
#endif
			if (_compute_estimators_for_image(args, i, threads, thread_id)) {
				siril_log_color_message(_("%s Check image %d first.\n"), "red",
						error_msg, args->image_indices[i] + 1);
				set_progress_bar_data(error_msg, PROGRESS_NONE);
				retval = 1;
				continue;
			}
			g_atomic_int_inc(&cur_nb);	// only used for progress bar
			set_progress_bar_data(NULL, cur_nb / (double)args->nb_images_to_stack);
		}
	}

	free(threads_per_thread);
	if (!retval)
		compute_factors_from_estimators(args, ref_image_filtred_idx);

	set_progress_bar_data(NULL, PROGRESS_DONE);
	return retval;
}

static void solve_overlap_coeffs(int nb_frames, int *index, int index_ref, size_t **Nij, double **Mij, gboolean additive, double *coeffs) {
	int c = 0;
	int N = nb_frames - 1;
	double *A = calloc(N * N, sizeof(double));
	double *B = calloc(N, sizeof(double));
	for (int i = 0; i < N; i ++) {
		int ii = index[i];
		B[i] = (additive) ? (double)Nij[ii][index_ref] * (Mij[index_ref][ii] - Mij[ii][index_ref]) :
							(double)Nij[ii][index_ref] * Mij[index_ref][ii] * Mij[ii][index_ref];
		for (int j = 0; j < N; j ++) {
			int ij = index[j];
			if (ii == ij) {
				for (int k = 0; k < nb_frames; k++) {
					if (k != ii)
						A[c] += (additive) ? (double)Nij[ii][k] :
											 (double)Nij[ii][k] * Mij[ii][k] * Mij[ii][k];
				}
			} else {
				A[c] = (additive) ? -(double)Nij[ii][ij] :
									-(double)Nij[ii][ij] * Mij[ii][ij] * Mij[ij][ii];
				if (additive)
					B[i] += Nij[ii][ij] * (Mij[ij][ii] - Mij[ii][ij]);
			}
			c++;
		}
	}
	// and solve the system
#ifdef DEBUG_NORM
	if (additive) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				siril_debug_print("%10d ", (int)A[j + i * N]);
			}
			siril_debug_print("; %16.6f\n", B[i]);
		} 
	} else {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				siril_debug_print("%+10.3f ", A[j + i * N]);
			}
			siril_debug_print("; %10.6f\n", B[i]);
		}
	}
#endif
	gsl_matrix_view m = gsl_matrix_view_array(A, N, N);
	gsl_vector_view b = gsl_vector_view_array(B, N);
	gsl_vector *x = gsl_vector_alloc(N);
	int s;
	gsl_permutation * p = gsl_permutation_alloc(N);
	gsl_linalg_LU_decomp(&m.matrix, p, &s);
	gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);
#ifdef DEBUG_NORM
	gsl_vector_fprintf(stdout, x, "%g");
#endif
	gsl_permutation_free(p);
	free(A);
	free(B);
	memcpy(coeffs, x->data, N * sizeof(double));
	gsl_vector_free(x);
}

// Computes the overlap areas on 2 images with regdata
// areai and areaj can be NULL if not required (to compute only overlap size for instance)
// returns the number of pixels in the common area
static size_t compute_overlap(struct stacking_args *args, int i, int j, rectangle *areai, rectangle *areaj) {
	sequence *seq = args->seq;
	double dxi = 0., dyi = 0., dxj = 0., dyj = 0.;
	int dx = 0, dy = 0;
	translation_from_H(seq->regparam[args->reglayer][i].H, &dxi, &dyi);
	translation_from_H(seq->regparam[args->reglayer][j].H, &dxj, &dyj);
	dx = round_to_int(dxj - dxi);
	dy = round_to_int(dyi - dyj);
	int rxi = (seq->is_variable) ? seq->imgparam[i].rx : seq->rx;
	int ryi = (seq->is_variable) ? seq->imgparam[i].ry : seq->ry;
	int rxj = (seq->is_variable) ? seq->imgparam[j].rx : seq->rx;
	int ryj = (seq->is_variable) ? seq->imgparam[j].ry : seq->ry;
	int x_tli = max(0, dx);
	int y_tli = max(0, dy);
	int x_bri = min(rxi, dx + rxj);
	int y_bri = min(ryi, dy + ryj);
	int x_tlj = max(0, -dx);
	int y_tlj = max(0, -dy);
	int x_brj = min(rxj, -dx + rxi);
	int y_brj = min(ryj, -dy + ryi);
	if (x_tli < x_bri && y_tli < y_bri) {
		if (areai)
			*areai = (rectangle) { x_tli, y_tli, x_bri - x_tli, y_bri - y_tli };
		if (areaj)
			*areaj = (rectangle) { x_tlj, y_tlj, x_brj - x_tlj, y_brj - y_tlj };
		return areai->w * areai->h;
	}
	return 0;
}

static int _compute_estimators_for_images(struct stacking_args *args, int i, int j, overlap_stats_t *stats,
		threading_type threading, int image_thread_id) {
	int nb_layers = args->seq->nb_layers;
	g_assert(nb_layers <= 3);
	g_assert(threading > 0);
	rectangle areai = { 0 };
	rectangle areaj = { 0 };
	size_t Nij = 0;
	float invnorm = 1. / USHRT_MAX;
	fits fiti = { 0 };
	fits fitj = { 0 };
	char file_i[256], file_j[256];
	sequence *seq = args->seq;

	siril_debug_print("computing stats for overlap between image %d and image %d, lite: %d\n", i + 1, j + 1, args->lite_norm);
	size_t nbdata = compute_overlap(args, i, j, &areai, &areaj);

	if (nbdata > 0) {
		siril_debug_print("image %d: boxselect %d %d %d %d\n", i + 1, areai.x, areai.y, areai.w, areai.h);
		siril_debug_print("image %d: boxselect %d %d %d %d\n", j + 1, areaj.x, areaj.y, areaj.w, areaj.h);
	} else {
		siril_debug_print("No overlap between image %d and %d\n", i + 1, j + 1);
		return 0; // no overlap
	}

	float *datai = malloc(nbdata * sizeof(float));
	float *dataj = malloc(nbdata * sizeof(float));
	fit_sequence_get_image_filename(seq, i, file_i, TRUE);
	fit_sequence_get_image_filename(seq, j, file_j, TRUE);
	if (readfits_partial_all_layers(file_i, &fiti, &areai) ||
		readfits_partial_all_layers(file_j, &fitj, &areaj)) {
		siril_log_color_message(_("Could not read overlap data between image %d and %d\n"), "red", i + 1, j + 1);
		clearfits(&fiti);
		clearfits(&fitj);
		free(datai);
		free(dataj);
		return -1;
	}
	for (int n = 0; n < nb_layers; n++) {
		Nij = 0;
		for (size_t k = 0; k < nbdata; k++) {
			if (fiti.type == DATA_FLOAT) {
				if (!fiti.fpdata[n][k] || !fitj.fpdata[n][k])
					continue;
				datai[Nij] = fiti.fpdata[n][k];
				dataj[Nij] = fitj.fpdata[n][k];
			} else if (fiti.type == DATA_USHORT) {
				if (!fiti.pdata[n][k] || !fitj.pdata[n][k])
					continue;
				datai[Nij] = (float)fiti.pdata[n][k] * invnorm;
				dataj[Nij] = (float)fitj.pdata[n][k] * invnorm;
			}
			Nij++;
		}
		siril_debug_print("%lu pixels for image %d and %d on layer %d\n", Nij, i + 1, j + 1, n);
		stats[n].Nij = Nij;
		stats[n].mij = (float)histogram_median_float(datai, Nij, SINGLE_THREADED);
		stats[n].mji = (float)histogram_median_float(dataj, Nij, SINGLE_THREADED);
		stats[n].sij = (float)siril_stats_float_mad(datai, Nij, stats[n].mij, SINGLE_THREADED, NULL);
		stats[n].sji = (float)siril_stats_float_mad(dataj, Nij, stats[n].mji, SINGLE_THREADED, NULL);
		if (!args->lite_norm) {
			double li = 0., lj = 0., si = 1., sj = 1.;
			IKSSlite(datai, Nij, stats[n].mij, stats[n].sij, &li, &si, SINGLE_THREADED);
			IKSSlite(dataj, Nij, stats[n].mji, stats[n].sji, &lj, &sj, SINGLE_THREADED);
			stats[n].mij = (float)li;
			stats[n].mji = (float)lj;
			stats[n].sij = (float)si;
			stats[n].sji = (float)sj;
		}
	}
	free(datai);
	free(dataj);
	clearfits(&fiti);
	clearfits(&fitj);
	return ST_OK;
}

static int compute_normalization_overlaps(struct stacking_args *args) {
	int index_ref = -1, retval = 0, cur_nb = 0, c = 0;
	norm_coeff *coeff = &args->coeff;
	int nb_layers = args->seq->nb_layers;
	int nb_frames = args->nb_images_to_stack;
	int N = nb_frames * (nb_frames - 1) / 2;

	if (args->normalize == NO_NORM || !args->maximize_framing)	// should never happen here
		return 0;

	init_coeffs(args);

	double ***Mij = malloc(nb_layers * sizeof(double **));
	double ***Sij = malloc(nb_layers * sizeof(double **));
	size_t ***Nij = malloc(nb_frames * sizeof(size_t **));
	for (int n = 0; n < nb_layers; n++) {
		Mij[n] = malloc(nb_frames * sizeof(double *));
		Sij[n] = malloc(nb_frames * sizeof(double *));
		Nij[n] = malloc(nb_frames * sizeof(size_t *));
		for (int i = 0; i < nb_frames; i++) {
			Mij[n][i] = calloc(nb_frames, sizeof(double));
			Sij[n][i] = calloc(nb_frames, sizeof(double));
			Nij[n][i] = calloc(nb_frames, sizeof(size_t));
		}
	}
	int *index = malloc((nb_frames - 1) * sizeof(int));

	char *tmpmsg = siril_log_message(_("Computing overlaps normalization...\n"));
	tmpmsg[strlen(tmpmsg) - 1] = '\0';
	set_progress_bar_data(tmpmsg, PROGRESS_RESET);

	// first, find the index of the ref image in the filtered image list
	index_ref = find_refimage_in_indices(args->image_indices,
			args->nb_images_to_stack, args->ref_image);
	if (index_ref == -1) {
		siril_log_color_message(_("The reference image is not in the selected set of images. "
				"Please choose another reference image.\n"), "red");
		free(Mij);
		return ST_GENERIC_ERROR;
	}

	const char *error_msg = (_("Normalization failed."));
	// check memory first
	// int nb_threads = normalization_get_max_number_of_threads(args->seq);
	// if (nb_threads <= 0) {
	// 	set_progress_bar_data(error_msg, PROGRESS_NONE);
	// 	return ST_GENERIC_ERROR;
	// }
	// if (nb_threads > args->nb_images_to_stack)
	// 	nb_threads = args->nb_images_to_stack;
	// int *threads_per_thread = compute_thread_distribution(nb_threads, com.max_thread);

	set_progress_bar_data(NULL, 0.);

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(nb_threads) schedule(guided) if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ) && fits_is_reentrant()))
// #endif
	for (int i = 0; i < nb_frames; ++i) {
		if (i != index_ref)
			index[c++] = i; // getting the filtered indexes of nonref images
		int ii = args->image_indices[i];
		for (int j = i + 1; j < nb_frames; ++j) {
			if (!retval) {
				if (!get_thread_run()) {
					retval = 1;
					continue;
				}
				int ij = args->image_indices[j];
				int thread_id = -1;
				threading_type threads = SINGLE_THREADED;
// #ifdef _OPENMP
// 			thread_id = omp_get_thread_num();
// 			threads = threads_per_thread[thread_id];
// #endif
				overlap_stats_t *ostats = calloc(nb_layers, sizeof(overlap_stats_t));
				if (_compute_estimators_for_images(args, ii, ij, ostats, threads, thread_id)) {
					siril_log_color_message(_("%s Check image %d first.\n"), "red",
							error_msg, args->image_indices[i] + 1);
					set_progress_bar_data(error_msg, PROGRESS_NONE);
					retval = 1;
					continue;
				}
				for (int n = 0; n < nb_layers; n++) {
					Mij[n][i][j] = ostats[n].mij;
					Mij[n][j][i] = ostats[n].mji;
					Sij[n][i][j] = ostats[n].sij;
					Sij[n][j][i] = ostats[n].sji;
					Nij[n][i][j] = ostats[n].Nij;
					Nij[n][j][i] = ostats[n].Nij;
				}
				g_atomic_int_inc(&cur_nb);	// only used for progress bar
				set_progress_bar_data(NULL, cur_nb / (double)N);
				free(ostats);
			}
		}
	}

#ifdef DEBUG_NORM
	for (int n = 0; n < nb_layers; n++) {
		for (int i = 0; i < nb_frames; i ++) {
			for (int j = i + 1; j < nb_frames; ++j) {
				siril_debug_print("%d;%d;%lu;%.6f;%.6f;%.6f;%.6f\n", i + 1, j + 1, Nij[n][i][j], Mij[n][i][j], Mij[n][j][i], Sij[n][i][j], Sij[n][j][i]);
			}
		}
	}
	for (int n = 0; n < nb_layers; n++) {
		siril_debug_print("Nij-%d\n", n);
		for (int i = 0; i < nb_frames; i ++) {
			for (int j = 0; j < nb_frames; ++j) {
				siril_debug_print("%10d ", (int)Nij[n][i][j]);
			}
			siril_debug_print("\n");
		}
		siril_debug_print("\n");
	}

	for (int n = 0; n < nb_layers; n++) {
		siril_debug_print("Mij-%d\n", n);
		for (int i = 0; i < nb_frames; i ++) {
			for (int j = 0; j < nb_frames; ++j) {
				siril_debug_print("%.6f ", Mij[n][i][j]);
			}
			siril_debug_print("\n");
		}
		siril_debug_print("\n");
	}

	for (int n = 0; n < nb_layers; n++) {
		siril_debug_print("Sij-%d\n", n);
		for (int i = 0; i < nb_frames; i ++) {
			for (int j = 0; j < nb_frames; ++j) {
				siril_debug_print("%.6f ", Sij[n][i][j]);
			}
			siril_debug_print("\n");
		}
		siril_debug_print("\n");
	}
#endif

	N = nb_frames - 1;
	double *coeffs = malloc(N * sizeof(double));

	if (args->normalize == ADDITIVE || args->normalize == ADDITIVE_SCALING) {
		for (int n = 0; n < nb_layers; n++) {
			solve_overlap_coeffs(nb_frames, index, index_ref, Nij[n], Mij[n], TRUE, coeffs);
			for (int i = 0; i < N; i ++) {
				coeff->poffset[n][index[i]] = -coeffs[i];
			}
		}
	}

	if (args->normalize == MULTIPLICATIVE_SCALING || args->normalize == ADDITIVE_SCALING) {
		for (int n = 0; n < nb_layers; n++) {
			solve_overlap_coeffs(nb_frames, index, index_ref, Nij[n], Sij[n], FALSE, coeffs);
			for (int i = 0; i < N; i ++) {
				coeff->pscale[n][index[i]] = 1./ coeffs[i];
			}
		}
	}

	if (args->normalize == MULTIPLICATIVE) {
		for (int n = 0; n < nb_layers; n++) {
			solve_overlap_coeffs(nb_frames, index, index_ref, Nij[n], Mij[n], FALSE, coeffs);
			for (int i = 0; i < N; i ++) {
				coeff->pmul[n][index[i]] = 1. / coeffs[i];
			}
		}
	}

	// free(threads_per_thread);
	if (!retval)
		compute_factors_from_estimators(args, index_ref);

	for (int n = 0; n < nb_layers; n++) {
		for (int i = 0; i < nb_frames; i++) {
			free(Mij[n][i]);
			free(Sij[n][i]);
			free(Nij[n][i]);
		}
		free(Mij[n]);
		free(Sij[n]);
		free(Nij[n]);
	}
	free(Mij);
	free(Sij);
	free(Nij);
	free(index);

	set_progress_bar_data(NULL, PROGRESS_DONE);
	return retval;
}

