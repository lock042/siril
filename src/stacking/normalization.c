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

// #define DEBUG_NORM

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
#ifdef DEBUG_NORM
	siril_debug_print("Normalization coeeficients\n");
#endif
	int reglayer = (args->reglayer > -1) ? args->reglayer : 1;
	for (int layer = 0; layer < nb_layers; ++layer) {
		int reflayer = (args->equalizeRGB) ? reglayer : layer;
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
			siril_debug_print("%2d %d %+8.5f %5f %.5f\n", args->image_indices[i] + 1, layer, poffset[layer][i], pscale[layer][i], pmul[layer][i]);
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

void free_ostats(overlap_stats_t **ostats, int nb_layers) {
	for (int n = 0; n < nb_layers; n++) {
		free(ostats[n]);
	}
	free(ostats);
}

overlap_stats_t **alloc_ostats(int nb_layers, int nb_frames) {
	g_assert(nb_layers > 0);
	int Npairs = nb_frames * (nb_frames - 1) / 2;
	overlap_stats_t **seq_ostats = calloc(nb_layers, sizeof(overlap_stats_t *));
	if (!seq_ostats) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	for (int n = 0; n < nb_layers; n++) {
		seq_ostats[n] = calloc(Npairs, sizeof(overlap_stats_t));
		if (!seq_ostats[n]) {
			PRINT_ALLOC_ERR;
			free_ostats(seq_ostats, nb_layers);
			return NULL;
		}
		for (int j = 0; j < Npairs; j++) {
			seq_ostats[n][j].i = -1;
			seq_ostats[n][j].j = -1;
			seq_ostats[n][j].Nij = 0;
			seq_ostats[n][j].locij = NULL_STATS;
			seq_ostats[n][j].locji = NULL_STATS;
			seq_ostats[n][j].scaij = NULL_STATS;
			seq_ostats[n][j].scaji = NULL_STATS;
			seq_ostats[n][j].medij = NULL_STATS;
			seq_ostats[n][j].medji = NULL_STATS;
			seq_ostats[n][j].madij = NULL_STATS;
			seq_ostats[n][j].madji = NULL_STATS;
			seq_ostats[n][j].areai = (rectangle){ 0, 0, 0, 0 };
			seq_ostats[n][j].areaj = (rectangle){ 0, 0, 0, 0 };
		}
	}
	return seq_ostats;
}

// returns the index in a list of (i,j) zero-based
// with i incremented first, then j (j > i to avoid permutations)
// with i in [0...N-1] and j in [i+1...N-1]
// Say for N = 4
// +---+---+------+
// | i | j | ijth |
// +---+---+------+
// | 0 | 1 |    0 |
// | 0 | 2 |    1 |
// | 0 | 3 |    2 |
// | 1 | 2 |    3 |
// | 1 | 3 |    4 |
// | 2 | 3 |    5 |
// +---+---+------+
int get_ijth_pair_index(int N, int i, int j) {
	return i * (2 * N - i -1) / 2 + j - i - 1;
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
		return (x_bri - x_tli) * (y_bri - y_tli);
	}
	return 0;
}

static int _compute_estimators_for_images(struct stacking_args *args, int i, int j) {
	int nb_layers = args->seq->nb_layers;
	g_assert(nb_layers <= 3);
	rectangle areai = { 0 };
	rectangle areaj = { 0 };
	size_t Nij = 0;
	float invnorm = 1. / USHRT_MAX;
	fits fiti = { 0 };
	fits fitj = { 0 };
	char file_i[256], file_j[256];
	sequence *seq = args->seq;
	size_t nbdata = 0;
	gboolean was_cached = FALSE;

	int ijth = get_ijth_pair_index(seq->number, i, j);
	if (seq->ostats[args->reglayer][ijth].i == i && seq->ostats[args->reglayer][ijth].j == j) { // we have some overlap data for this pair
		nbdata = seq->ostats[args->reglayer][ijth].areai.w * seq->ostats[args->reglayer][ijth].areai.h;
		areai = seq->ostats[args->reglayer][ijth].areai;
		areaj = seq->ostats[args->reglayer][ijth].areaj;
		was_cached = TRUE;
	} else {
		args->seq->needs_saving = TRUE;
		siril_debug_print("computing stats for overlap between image %d and image %d, lite: %d\n", i + 1, j + 1, args->lite_norm);
		nbdata = compute_overlap(args, i, j, &areai, &areaj);
		// we cache it for all layers
		// Normally, we should have regdata for only one layer, but what if we have for more (can't see that happening but better be safe)
		// In that case, we will assume that the differences in overlaps should be minimal (1 or 2 lines) and that
		// the first ever cached overlap stats are valid irrespective of the regdata which created them
		for (int n = 0; n < nb_layers; n++) {
			seq->ostats[n][ijth].i = i;
			seq->ostats[n][ijth].j = j;
			if (nbdata > 0) {
				seq->ostats[n][ijth].areai = areai;
				seq->ostats[n][ijth].areaj = areaj;
			}
		}
	}

	if (nbdata > 0) {
		siril_debug_print("image %d: boxselect %d %d %d %d\n", i + 1, areai.x, areai.y, areai.w, areai.h);
		siril_debug_print("image %d: boxselect %d %d %d %d\n", j + 1, areaj.x, areaj.y, areaj.w, areaj.h);
	} else {
		siril_debug_print("No overlap between image %d and %d\n", i + 1, j + 1);
		return ST_OK; // no overlap
	}

	// we'll now determine if we need to compute anything, and if yes, which data
	// we'll check on each layer even though some cases are obvious
	// For instance, with litenorm, we should have data on all 3 layers if we have reached this step
	gboolean needs_recalc_lite[3] = { 0 };
	gboolean needs_recalc_detailed[3] = { 0 };
	gboolean needs_recalc = FALSE;
	if (was_cached) {
		for (int n = 0; n < nb_layers; n++) {
			needs_recalc_lite[n] =	seq->ostats[n][ijth].madij == NULL_STATS || 
									seq->ostats[n][ijth].medij == NULL_STATS || 
									seq->ostats[n][ijth].madji == NULL_STATS || 
									seq->ostats[n][ijth].medji == NULL_STATS;
			needs_recalc_detailed[n] = !args->lite_norm && (
									seq->ostats[n][ijth].locij == NULL_STATS || 
									seq->ostats[n][ijth].scaij == NULL_STATS || 
									seq->ostats[n][ijth].locji == NULL_STATS || 
									seq->ostats[n][ijth].scaji == NULL_STATS );
			needs_recalc |= needs_recalc_lite[n] || needs_recalc_detailed[n];
		}
	} else {
		needs_recalc = TRUE;
		for (int n = 0; n < nb_layers; n++) {
			needs_recalc_lite[n] = TRUE;
			needs_recalc_detailed[n] = (!args->lite_norm);
		}
	}

	if (!needs_recalc) {
		siril_debug_print("Data for %d and %d were cached, re-using\n", i + 1, j + 1);
		return ST_OK;
	}
	args->seq->needs_saving = TRUE;
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
		if (!needs_recalc_lite[n] && !needs_recalc_detailed[n])
			continue;
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
		if (Nij > 3) { // we want at least 3 pixels with non-zero data to compute stats
			seq->ostats[n][ijth].Nij = Nij;
			if (needs_recalc_lite[n]) {
				siril_debug_print("%lu pixels for images %d and %d on layer %d\n", Nij, i + 1, j + 1, n);
				seq->ostats[n][ijth].medij = (float)histogram_median_float(datai, Nij, SINGLE_THREADED);
				seq->ostats[n][ijth].medji = (float)histogram_median_float(dataj, Nij, SINGLE_THREADED);
				seq->ostats[n][ijth].madij = (float)siril_stats_float_mad(datai, Nij, seq->ostats[n][ijth].medij, SINGLE_THREADED, NULL);
				seq->ostats[n][ijth].madji = (float)siril_stats_float_mad(dataj, Nij, seq->ostats[n][ijth].medji, SINGLE_THREADED, NULL);
			} else {
				siril_debug_print("Lite overlap stats for images %d and %d on layer %d read from cache\n", i + 1, j + 1, n);
			}
			if (needs_recalc_detailed[n]) {
				double li = 0., lj = 0., si = 1., sj = 1.;
				IKSSlite(datai, Nij, seq->ostats[n][ijth].medij, seq->ostats[n][ijth].madij, &li, &si, SINGLE_THREADED);
				IKSSlite(dataj, Nij, seq->ostats[n][ijth].medji, seq->ostats[n][ijth].madji, &lj, &sj, SINGLE_THREADED);
				seq->ostats[n][ijth].locij = (float)li;
				seq->ostats[n][ijth].locji = (float)lj;
				seq->ostats[n][ijth].scaij = (float)si;
				seq->ostats[n][ijth].scaji = (float)sj;
			}  else {
				siril_debug_print("Detailed overlap stats for images %d and %d on layer %d read from cache\n", i + 1, j + 1, n);
				siril_debug_print("Should not happen\n");
			}
		} else {
			siril_debug_print("No overlap data between image %d and %d on layer %d\n", i + 1, j + 1, n);
		}
	}
	free(datai);
	free(dataj);
	clearfits(&fiti);
	clearfits(&fitj);
	return ST_OK;
}

static int normalization_overlap_get_max_number_of_threads(struct stacking_args *args) {
	int max_memory_MB = get_max_memory_in_MB();
	int nb_frames = args->nb_images_to_stack;
	int nb_layers = args->seq->nb_layers;
	sequence *seq = args->seq;
	/* The overlap normalization memory consumption assumes:
		- n is max overlap size accross all pairs of images
		- l in the number of layers
		We need:
		- 2 * overlap parial image stored (so 2.n.l) when we read the images
		- 2 * data copies as floats (so 2.n as float) that are re-used for all layers
		- 1 * data additional copy to compute MAD (not 2 because the calcs are 
		sequential and memory freed in between calls)
	*/
	size_t nbdatamax = 0;
	for (int i = 0; i < nb_frames; ++i) {
		int ii = args->image_indices[i];
		for (int j = i + 1; j < nb_frames; ++j) {
			int ij = args->image_indices[j];
			int ijth = get_ijth_pair_index(seq->number, ii, ij);
			size_t nbdata = 0;
			if (seq->ostats[args->reglayer][ijth].i == ii && seq->ostats[args->reglayer][ijth].j == ij) { // we have some overlap data for this pair
				nbdata = seq->ostats[args->reglayer][ijth].areai.w * seq->ostats[args->reglayer][ijth].areai.h;
			} else {
				args->seq->needs_saving = TRUE;
				rectangle areai, areaj;
				siril_debug_print("computing stats for overlap between image %d and image %d, lite: %d\n", ii + 1, ij + 1, args->lite_norm);
				nbdata = compute_overlap(args, ii, ij, &areai, &areaj);
				// we cache it for all layers
				// Normally, we should have regdata for only one layer, but what if we have for more (can't see that happening but better be safe)
				// In that case, we will assume that the differences in overlaps should be minimal (1 or 2 lines) and that
				// the first ever cached overlap stats are valid irrespective of the regdata which created them
				for (int n = 0; n < nb_layers; n++) {
					seq->ostats[n][ijth].i = ii;
					seq->ostats[n][ijth].j = ij;
					if (nbdata > 0) {
						seq->ostats[n][ijth].areai = areai;
						seq->ostats[n][ijth].areaj = areaj;
					}
				}
			}
			if (nbdata > nbdatamax)
				nbdatamax = nbdata;
		}
	}
	size_t memory_per_pair = 2 * seq->nb_layers * nbdatamax * (get_data_type(seq->bitpix) == DATA_FLOAT ? sizeof(float) : sizeof(WORD)); // we size memory consumption with the largest possible overlap
	size_t memory_per_stats = 2 * nbdatamax * sizeof(float);
	unsigned int memory_per_pair_MB = ( memory_per_pair + memory_per_stats) / BYTES_IN_A_MB;
	if (memory_per_pair_MB == 0)
		memory_per_pair_MB = 1;

	fprintf(stdout, "Memory per pair: %u MB. Max memory: %d MB\n", memory_per_pair_MB, max_memory_MB);

	if (memory_per_pair_MB > max_memory_MB) {
		siril_log_color_message(_("Your system does not have enough memory to normalize images overlaps for stacking operation (%d MB free for %d MB required)\n"), "red", max_memory_MB, memory_per_pair_MB);
		return 0;
	}

	int nb_threads = max_memory_MB / memory_per_pair_MB;
	if (nb_threads > com.max_thread)
		nb_threads = com.max_thread;
	siril_log_message(_("With the current memory and thread (%d) limits, up to %d thread(s) can be used for sequence overlaps normalization\n"), com.max_thread, nb_threads);
	return nb_threads;
}

static int compute_normalization_overlaps(struct stacking_args *args) {
	int index_ref = -1, retval = 0, cur_nb = 0, c = 0;
	norm_coeff *coeff = &args->coeff;
	int nb_layers = args->seq->nb_layers;
	int nb_frames = args->nb_images_to_stack;
	int Npairs = nb_frames * (nb_frames - 1) / 2;
	int N = nb_frames - 1;
	double *coeffs = NULL;
	double ***Mij = NULL, ***Sij = NULL;
	size_t ***Nij = NULL;
	int *index = NULL;
	// imstats *refstats[3] = { NULL };

	if (args->normalize == NO_NORM || !args->maximize_framing)  // should never happen here
		return 0;

	if (args->force_norm) { /* We empty the cache if needed (force to recompute) */
		free_ostats(args->seq->ostats, args->seq->nb_layers);
		args->seq->ostats = NULL;
		args->seq->needs_saving = TRUE;
	}

	// first, find the index of the ref image in the filtered image list
	// and compute its statistics
	index_ref = find_refimage_in_indices(args->image_indices,
										 args->nb_images_to_stack, args->ref_image);
	if (index_ref == -1) {
		siril_log_color_message(_("The reference image is not in the selected set of images. "
		"Please choose another reference image.\n"), "red");
		return ST_GENERIC_ERROR;
	}

	//TODO: check if we can still equalize RGB

	//  if (compute_all_channels_statistics_seqimage(args->seq, index_ref, NULL, (args->lite_norm) ? STATS_LITENORM : STATS_NORM, SINGLE_THREADED, -1, refstats)) {
	// #ifdef DEBUG_NORM
	//      for (int n = 0; n < nb_layers; n++) {
	//          if (args->lite_norm)
	//              siril_debug_print("Reference %.6f %.6f\n", refstats[n]->median, refstats[n]->mad);
	//          else
	//              siril_debug_print("Reference %.6f %.6f\n", refstats[n]->location, refstats[n]->scale);
	//      }
	// #endif
	//      siril_log_color_message(_("Could not compute statistics of reference image"), "red");
	//      retval = 1;
	//      return ST_GENERIC_ERROR;
	//  }

	init_coeffs(args);
	// if the overlap stats have never been cached or have been cleared, we allocate there for all
	// the images of the sequence (without filtering)
	if (!args->seq->ostats) {
		args->seq->ostats = alloc_ostats(nb_layers, args->seq->number);
		if (!args->seq->ostats) {
			retval  = 1;
			goto cleanup;
		}
		args->seq->needs_saving = TRUE;
	}

	Mij = calloc(nb_layers, sizeof(double **));
	Sij = calloc(nb_layers, sizeof(double **));
	Nij = calloc(nb_layers, sizeof(size_t **));
	index = calloc(nb_frames, sizeof(int));
	coeffs = calloc(N, sizeof(double));
	if (!Mij || !Sij || !Nij || !index || !coeffs) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto cleanup2;
	}
	for (int n = 0; n < nb_layers; n++) {
		Mij[n] = malloc(nb_frames * sizeof(double *));
		Sij[n] = malloc(nb_frames * sizeof(double *));
		Nij[n] = malloc(nb_frames * sizeof(size_t *));
		if (!Mij[n] || !Sij[n] || !Nij[n]) {
			PRINT_ALLOC_ERR;
			retval = 1;
			goto cleanup2;
		}
		for (int i = 0; i < nb_frames; i++) {
			Mij[n][i] = calloc(nb_frames, sizeof(double));
			Sij[n][i] = calloc(nb_frames, sizeof(double));
			Nij[n][i] = calloc(nb_frames, sizeof(size_t));
			if (!Mij[n][i] || !Sij[n][i] || !Nij[n][i]) {
				PRINT_ALLOC_ERR;
				retval = 1;
				goto cleanup2;
			}
		}
	}

	char *tmpmsg = siril_log_message(_("Computing normalization on overlaps...\n"));
	tmpmsg[strlen(tmpmsg) - 1] = '\0';
	set_progress_bar_data(tmpmsg, PROGRESS_RESET);

	const char *error_msg = (_("Normalization failed."));
	// check memory first
	int nb_threads = normalization_overlap_get_max_number_of_threads(args);
	if (nb_threads <= 0) {
		set_progress_bar_data(error_msg, PROGRESS_NONE);
		retval = ST_GENERIC_ERROR;
		goto cleanup;
	}
	if (nb_threads > args->nb_images_to_stack)
		nb_threads = args->nb_images_to_stack;

	set_progress_bar_data(NULL, 0.);

	for (int i = 0; i < nb_frames; ++i) {
		if (i != index_ref)
			index[c++] = i; // getting the filtered indexes of nonref images
	}

	#ifdef _OPENMP
	#pragma omp parallel for num_threads(nb_threads) schedule(guided) if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ) && fits_is_reentrant()))
	#endif
	for (int i = 0; i < nb_frames; ++i) {
		int ii = args->image_indices[i];
		for (int j = i + 1; j < nb_frames; ++j) {
			if (!retval) {
				if (!get_thread_run()) {
					retval = 1;
					continue;
				}
				int ij = args->image_indices[j];
				int ijth = get_ijth_pair_index(args->seq->number, ii, ij);
				if (_compute_estimators_for_images(args, ii, ij)) {
					siril_log_color_message(_("%s Check image %d first.\n"), "red",
											error_msg, args->image_indices[i] + 1);
					set_progress_bar_data(error_msg, PROGRESS_NONE);
					#ifdef _OPENMP
					#pragma omp critical
					#endif
					{
						retval = 1;
					}
					continue;
				}
				for (int n = 0; n < nb_layers; n++) {
					if (args->seq->ostats[n][ijth].Nij == 0)
						continue;
					if (args->lite_norm) {
						Mij[n][i][j] = args->seq->ostats[n][ijth].medij;
						Mij[n][j][i] = args->seq->ostats[n][ijth].medji;
						Sij[n][i][j] = args->seq->ostats[n][ijth].madij;
						Sij[n][j][i] = args->seq->ostats[n][ijth].madji;
					} else {
						Mij[n][i][j] = args->seq->ostats[n][ijth].locij;
						Mij[n][j][i] = args->seq->ostats[n][ijth].locji;
						Sij[n][i][j] = args->seq->ostats[n][ijth].scaij;
						Sij[n][j][i] = args->seq->ostats[n][ijth].scaji;
					}

					Nij[n][i][j] = args->seq->ostats[n][ijth].Nij;
					Nij[n][j][i] = args->seq->ostats[n][ijth].Nij;
				}
				#ifdef _OPENMP
				#pragma omp atomic
				#endif
				cur_nb++;  // only used for progress bar
				set_progress_bar_data(NULL, cur_nb / (double)Npairs);
			}
		}
	}
	if (retval)
		goto cleanup;

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

	if (args->normalize == MULTIPLICATIVE_SCALING || args->normalize == ADDITIVE_SCALING) {
		for (int n = 0; n < nb_layers; n++) {
			solve_overlap_coeffs(nb_frames, index, index_ref, Nij[n], Sij[n], FALSE, coeffs);
			for (int i = 0; i < N; i ++) { // we set the coeffs for nb_frames - 1, the ref has already been init
				coeff->pscale[n][index[i]] = coeffs[i];
			}
			// we re-normalize the Mij matrix by the scales found at this step
			for (int ii = 0; ii < nb_frames; ii++) {
				for (int jj = 0; jj < nb_frames; jj++) {
					Mij[n][ii][jj] *= coeff->pscale[n][ii];
				}
			}
		}
	}

	if (args->normalize == ADDITIVE || args->normalize == ADDITIVE_SCALING) {
		for (int n = 0; n < nb_layers; n++) {
			solve_overlap_coeffs(nb_frames, index, index_ref, Nij[n], Mij[n], TRUE, coeffs);
			for (int i = 0; i < N; i ++) { // we set the coeffs for nb_frames - 1, the ref has already been init
				coeff->poffset[n][index[i]] = -coeffs[i];
			}
		}
	}

	if (args->normalize == MULTIPLICATIVE) {
		for (int n = 0; n < nb_layers; n++) {
			solve_overlap_coeffs(nb_frames, index, index_ref, Nij[n], Mij[n], FALSE, coeffs);
			for (int i = 0; i < N; i ++) {
				coeff->pmul[n][index[i]] = coeffs[i];
			}
		}
	}

	cleanup:
	for (int n = 0; n < nb_layers; n++) {
		if (Mij && Mij[n]) {
			for (int i = 0; i < nb_frames; i++) {
				free(Mij[n][i]);
			}
			free(Mij[n]);
		}
		if (Sij && Sij[n]) {
			for (int i = 0; i < nb_frames; i++) {
				free(Sij[n][i]);
			}
			free(Sij[n]);
		}
		if (Nij && Nij[n]) {
			for (int i = 0; i < nb_frames; i++) {
				free(Nij[n][i]);
			}
			free(Nij[n]);
		}
	}
	cleanup2:
	free(Mij);
	free(Sij);
	free(Nij);
	free(index);
	free(coeffs);

	set_progress_bar_data(NULL, PROGRESS_DONE);
	return retval;
}
