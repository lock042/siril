#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "stacking.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"

static int compute_normalization(struct stacking_args *args);

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

	if (compute_normalization(args)) {
		args->retval = ST_GENERIC_ERROR;
		return args->retval;
	}

	if (args->seq->needs_saving)	// if we had to compute new stats
		writeseqfile(args->seq);

	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, "Normalization computation time");
	return ST_OK;
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
		}
	}
}

static int normalization_get_max_number_of_threads(sequence *seq) {
	int max_memory_MB = get_max_memory_in_MB();
	/* The normalization memory consumption, n is image size and m channel size.
	 * It uses IKSS computation in stats, which can be done only on float data.
	 * IKSS requires computing the MAD, which requires its own copy of the data.
	 * The stats are computed successively on each channel.
	 * For DATA_USHORT, we have: the image O(n), rewrite of the channel without
	 * zeros O(m), a copy to float for IKSS O(2m), a copy of that for the MAD in
	 * IKSS O(2m).
	 * For DATA_FLOAT, we have: the image O(n), rewrite without zeros O(m),
	 * used directly for IKSS and a copy for MAD O(m).
	 */
	guint64 memory_per_image = seq->rx * seq->ry;
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

		/* We empty the cache if needed (force to recompute) */
		if (args->force_norm)
			clear_stats(args->seq, layer);
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
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb++;	// only used for progress bar
			set_progress_bar_data(NULL, cur_nb / (double)args->nb_images_to_stack);
		}
	}

	free(threads_per_thread);
	if (!retval)
		compute_factors_from_estimators(args, ref_image_filtred_idx);

	set_progress_bar_data(NULL, PROGRESS_DONE);
	return retval;
}

