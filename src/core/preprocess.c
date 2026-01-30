/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "algos/fix_xtrans_af.h"
#include "filters/cosmetic_correction.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "algos/demosaicing.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/ser.h"

#include "preprocess.h"

#define SQUARE_SIZE 512

static float evaluateNoiseOfCalibratedImage(fits *fit, fits *dark,
		float k, gboolean allow_32bit_output) {
	long chan;
	int ret = 0;
	float noise = 0.f;
	fits dark_tmp = { 0 }, fit_tmp = { 0 };
	rectangle area = { 0 };

	/* square of 512x512 in the center of the image */
	int size = SQUARE_SIZE;
	area.x = (fit->rx - size) / 2;
	area.y = (fit->ry - size) / 2;
	area.w = size;
	area.h = size;

	ret = copyfits(dark, &dark_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	if (!ret) ret = copyfits(fit, &fit_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);

	if (!ret) ret = soper(&dark_tmp, k, OPER_MUL, allow_32bit_output);
	if (!ret) ret = imoper(&fit_tmp, &dark_tmp, OPER_SUB, allow_32bit_output);
	if (ret) {
		clearfits(&dark_tmp);
		clearfits(&fit_tmp);
		return -1.0;
	}

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		/* STATS_SIGMEAN computes mean and normvalue */
		imstats *stat = statistics(NULL, -1, &fit_tmp, chan, &area, STATS_SIGMEAN, SINGLE_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return -1.0;
		}
		noise += stat->sigma / stat->normValue;
		free_stats(stat);
	}
	clearfits(&dark_tmp);
	clearfits(&fit_tmp);

	return noise;
}

#undef GR
#define GR ((sqrtf(5.f) - 1.f) / 2.f)

static float goldenSectionSearch(fits *raw, fits *dark, float a, float b,
		float tol, gboolean allow_32bits) {
	float c, d;
	float fc, fd;
	int iter = 0;

	c = b - GR * (b - a);
	d = a + GR * (b - a);
	fc = evaluateNoiseOfCalibratedImage(raw, dark, c, allow_32bits);
	fd = evaluateNoiseOfCalibratedImage(raw, dark, d, allow_32bits);
	if (fc == fd) return 1.f;
	do {
		siril_debug_print("Iter: %d (%1.2f, %1.2f)\n", ++iter, c, d);
		if (fc < 0.f || fd < 0.f)
			return -1.f;
		if (fc < fd) {
			b = d;
			d = c;
			fd = fc;
			c = b - GR * (b - a);
			fc = evaluateNoiseOfCalibratedImage(raw, dark, c, allow_32bits);
		} else {
			a = c;
			c = d;
			fc = fd;
			d = a + GR * (b - a);
			fd = evaluateNoiseOfCalibratedImage(raw, dark, d, allow_32bits);

		}
	} while (fabsf(c - d) > tol);
	return ((b + a) * 0.5f);
}

static int preprocess(fits *raw, struct preprocessing_data *args) {
	int ret = 0;

	if (args->use_bias) {
		if (args->bias_level < FLT_MAX) {
			// an offset level has been defined
			ret = soper(raw, args->bias_level, OPER_SUB, args->allow_32bit_output);
		}
		else ret = imoper(raw, args->bias, OPER_SUB, args->allow_32bit_output);
	}

	/* if dark optimization, the master-dark has already been subtracted */
	if (!ret && args->use_dark && !args->use_dark_optim) {
		ret = imoper(raw, args->dark, OPER_SUB, args->allow_32bit_output);
		if (raw->neg_ratio > 0.2f)
			siril_log_message(_("After dark subtraction, the image contains many negative pixels (%d%%), calibration frames are probably incorrect\n"), (int)(100.f*raw->neg_ratio));
#ifdef SIRIL_OUTPUT_DEBUG
		image_find_minmax(raw);
		fprintf(stdout, "after dark: min=%f, max=%f\n", raw->mini, raw->maxi);
		invalidate_stats_from_fit(raw);
#endif
	}

	if (!ret && args->use_flat) {
		// return value is an error if there is an overflow, but it is usual
		// for now, so don't treat as error
		/*ret =*/ siril_fdiv(raw, args->flat, args->normalisation, args->allow_32bit_output);
#ifdef SIRIL_OUTPUT_DEBUG
		image_find_minmax(raw);
		fprintf(stdout, "after flat: min=%f, max=%f\n", raw->mini, raw->maxi);
		invalidate_stats_from_fit(raw);
#endif
	}

	return ret;
}

static int darkOptimization(fits *raw, struct preprocessing_data *args, int in_index, GSList *hist) {
	float k0;
	float lo = 0.f, up = 2.f;
	int ret = 0;
	fits *dark = args->dark;
	fits dark_tmp = { 0 };

	if (memcmp(raw->naxes, dark->naxes, sizeof raw->naxes)) {
		siril_log_color_message(_("Images must have same dimensions\n"), "red");
		return 1;
	}

	// TODO: avoid duplicating with soper+imoper together
	// https://gitlab.com/free-astro/siril/-/issues/671
	if (copyfits(dark, &dark_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, 0))
		return 1;

	if (args->use_exposure) {
		if (dark->keywords.exposure <= 0.0) {
			siril_log_color_message(_("The dark frame contains no exposure data or incorrect exposure data.\n"), "red");
			clearfits(&dark_tmp);
			return 1;
		}
		if (raw->keywords.exposure <= 0.0) {
			siril_log_color_message(_("The light frame (%d) contains no exposure data or incorrect exposure data.\n"), "red", in_index + 1);
			clearfits(&dark_tmp);
			return 1;
		}
		/* linear scale with time */
		k0 = raw->keywords.exposure / dark->keywords.exposure;
		if (k0 > 1.f) {
			siril_log_color_message(_("Warning: master dark is shorter than lights. It is "
						"recommended that the master dark be at least as long as the lights.\n"), "salmon");
		}
	} else {
		/* Minimization of background noise to find better k */
		k0 = goldenSectionSearch(raw, &dark_tmp, lo, up, 0.001f, args->allow_32bit_output);
	}
	if (k0 < 0.f) {
		siril_log_message(_("Dark optimization of image %d failed\n"), in_index);
		ret = -1;
	} else {
		/* Multiply coefficient to master-dark */
		ret = soper(&dark_tmp, k0, OPER_MUL, args->allow_32bit_output);
		if (!ret)
			ret = imoper(raw, &dark_tmp, OPER_SUB, args->allow_32bit_output);
		if (!ret) {
			siril_log_message(_("Dark optimization of image %d: k0=%.3f\n"), in_index, k0);
			if (hist)
				hist = g_slist_append(hist, g_strdup_printf("Calibrated with an optimized master dark (factor: %.3f)", k0));
		}
		else siril_log_message(_("Dark optimization of image %d failed\n"), in_index);
	}
	clearfits(&dark_tmp);
	return ret;
}

static gint64 prepro_compute_size_hook(struct generic_seq_args *args, int nb_images) {
	struct preprocessing_data *prepro = args->user;
	gint64 size = seq_compute_size(args->seq, nb_images, args->output_type);
	if (prepro->debayer)
		size *= 3;
	return size;
}

/* we need to account for the possible 3 master frames loaded at the start and
 * the memory it takes to calibrate the images */
static int prepro_compute_mem_hook(struct generic_seq_args *args, gboolean for_writer) {
	int nb_masters = 0;
	struct preprocessing_data *prepro = args->user;
	if (prepro->use_flat && prepro->flat) nb_masters++;
	if (prepro->use_dark && prepro->dark) nb_masters++;
	if (prepro->use_bias && prepro->bias) nb_masters++;
	int input_is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;
	int output_is_float = input_is_float || !com.pref.force_16bit;
	/*
	 * are the masters float? yes if !com.pref.force_16bit, raw images too
	 *
	 * if (prepro->use_dark_optim && prepro->use_dark)
	 *	image gets three copies
	 *
	 * if (prepro->debayer)
	 * 	O(9n) for ushort
	 * 	O(3n) for float, same as dark optim
	 */
	unsigned int MB_per_input_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_input_image, NULL, &MB_avail);
	if (!input_is_float && output_is_float)
		MB_per_input_image *= 2;
	MB_avail -= nb_masters * MB_per_input_image;

	unsigned int required = MB_per_input_image;
	unsigned int MB_per_output_image = MB_per_input_image;
	if (prepro->debayer) {
		if (input_is_float || output_is_float)
			MB_per_output_image *= 3;
		else MB_per_output_image *= 9;
		required = MB_per_input_image + MB_per_output_image;
	}
	else if (prepro->use_dark_optim && prepro->use_dark) {
		required = 4 * MB_per_input_image;
	}

	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
                        thread_limit = com.max_thread;

		if (for_writer) {
                        /* we allow the already allocated thread_limit images,
                         * plus how many images can be stored in what remains
                         * unused by the main processing */
                        limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_output_image;
			siril_debug_print("%u MB avail for writer\n", MB_avail - required * thread_limit);
                } else limit = thread_limit;
	}
	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_debug_print("Memory required per thread: %u MB, per input image: %u MB, per output image %u MB, limiting to %d %s\n",
				required, MB_per_input_image, MB_per_output_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

int prepro_prepare_hook(struct generic_seq_args *args) {
	struct preprocessing_data *prepro = args->user;
	gboolean set_hist = prepro->output_seqtype != SEQ_SER;

	if (prepro->seq) {
		// handling SER and FITSEQ
		if (seq_prepare_hook(args))
			return 1;
	}

	if (set_hist) {
		if (prepro->use_dark && !prepro->use_dark_optim)
			prepro->history = g_slist_prepend(prepro->history, g_strdup("Calibrated with a master dark"));
		if (prepro->use_bias) {
			if (prepro->bias_level < FLT_MAX) {
				int bitpix = prepro->seq ? prepro->seq->bitpix : gfit.orig_bitpix;
				prepro->history = g_slist_prepend(prepro->history, g_strdup_printf("Calibrated with a synthetic bias of %d", roundf_to_int(prepro->bias_level * ((bitpix == BYTE_IMG) ? UCHAR_MAX_SINGLE : USHRT_MAX_SINGLE))));
			} else {
				prepro->history = g_slist_prepend(prepro->history, g_strdup("Calibrated with a master bias"));
			}
		}
	}

	// precompute flat levels
	if (prepro->use_flat) {
		if (prepro->equalize_cfa) {
			compute_grey_flat(prepro->flat);
		}
		if (prepro->autolevel) {
			const unsigned int width = prepro->flat->rx;
			const unsigned int height = prepro->flat->ry;

			/* due to vignetting it is better to take an area in the
			 * center of the flat image
			 */
			const unsigned int startx = width / 3;
			const unsigned int starty = height / 3;

			rectangle selection = { startx, starty, startx, starty };

			imstats *stat = statistics(NULL, -1, prepro->flat, RLAYER, &selection, STATS_BASIC, MULTI_THREADED);
			if (!stat) {
				siril_log_message(_("Error: statistics computation failed.\n"));
				return 1;
			}
			prepro->normalisation = stat->mean;
			free_stats(stat);

			// special cases for division factor, caused by how imoper works with many cases
			if (DATA_USHORT == prepro->flat->type && prepro->allow_32bit_output)
				prepro->normalisation *= INV_USHRT_MAX_SINGLE;
			else if (DATA_FLOAT == prepro->flat->type && !prepro->allow_32bit_output) {
				if (prepro->seq) { // sequence case
					printf("normalizing for seq\n");
					prepro->normalisation *= (prepro->seq->bitpix == BYTE_IMG) ? UCHAR_MAX_SINGLE : USHRT_MAX_SINGLE; // imoper deals with bitdepth through norm value
				} else { // single-image case
					printf("normalizing for single image\n");
					prepro->normalisation *= (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX_SINGLE : USHRT_MAX_SINGLE; // imoper deals with bitdepth through norm value
				}
			}

			siril_log_message(_("Normalisation value auto evaluated: %.3f\n"),
					prepro->normalisation);
		}
		if (set_hist) {
			prepro->history = g_slist_append(prepro->history, g_strdup_printf("Calibrated with a master flat, normalization of %.3f", prepro->normalisation));
		}
	}

	/** FIX XTRANS AC ISSUE **/
	if (prepro->fix_xtrans && prepro->use_dark) {
		fix_xtrans_ac(prepro->dark);
		if (set_hist)
			prepro->history = g_slist_append(prepro->history, g_strdup("Fixed X-Trans sensor artifact on the master dark"));
	}

	if (prepro->fix_xtrans && prepro->use_bias && prepro->bias) { // check for existence of bias in case of synth offset
		fix_xtrans_ac(prepro->bias);
		if (set_hist)
			prepro->history = g_slist_append(prepro->history, g_strdup("Fixed X-Trans sensor artifact on the master bias"));
	}

	// proceed to cosmetic correction
	if (prepro->use_cosmetic_correction && prepro->use_dark && prepro->cc_from_dark) {
		if (strlen(prepro->dark->keywords.bayer_pattern) > 4) {
			siril_log_color_message(_("Cosmetic correction cannot be applied on X-Trans files.\n"), "red");
			prepro->use_cosmetic_correction = FALSE;
		} else {
			if (prepro->dark->naxes[2] == 1) {
				prepro->dev = find_deviant_pixels(prepro->dark, prepro->sigma,
						&(prepro->icold), &(prepro->ihot), FALSE);
				gchar *str = ngettext("%ld corrected pixel (%ld + %ld)\n", "%ld corrected pixels (%ld + %ld)\n", prepro->icold + prepro->ihot);
				str = g_strdup_printf(str, prepro->icold + prepro->ihot, prepro->icold, prepro->ihot);
				siril_log_message(str);
				g_free(str);
				if (set_hist)
					prepro->history = g_slist_append(prepro->history, g_strdup_printf("Cosmetic correction of %ld cold pixels and %ld hot pixels", prepro->icold, prepro->ihot));
			} else
				siril_log_message(_("Darkmap cosmetic correction is only supported with single channel images\n"));
		}
	}

	return 0;
}

int prepro_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct preprocessing_data *prepro = args->user;
	GSList *history = g_slist_copy_deep(prepro->history, (GCopyFunc)g_strdup, NULL);

	if (prepro->use_dark_optim && prepro->use_dark) {
		if (darkOptimization(fit, prepro, in_index, history)) {
			g_slist_free_full(history, g_free);
			return 1;
		}
	}

	if (preprocess(fit, prepro)) {
		g_slist_free_full(history, g_free);
		return 1;
	}

	if (prepro->use_cosmetic_correction && prepro->use_dark
			&& prepro->dark->naxes[2] == 1 && prepro->cc_from_dark) {
		cosmeticCorrection(fit, prepro->dev, prepro->icold + prepro->ihot, prepro->is_cfa);
#ifdef SIRIL_OUTPUT_DEBUG
		image_find_minmax(fit);
		fprintf(stdout, "after cosmetic correction: min=%f, max=%f\n",
				fit->mini, fit->maxi);
		invalidate_stats_from_fit(fit);
#endif
	}

	if (prepro->use_cosmetic_correction && !prepro->cc_from_dark && prepro->bad_pixel_map_file) {
		apply_cosme_to_image(fit, prepro->bad_pixel_map_file, prepro->is_cfa);
	}

	if (prepro->debayer) {
		// debayering SER is now allowed - https://gitlab.com/free-astro/siril/-/issues/549
		if (!prepro->seq || prepro->seq->type == SEQ_REGULAR || prepro->seq->type == SEQ_FITSEQ || prepro->seq->type == SEQ_SER ) {
			debayer_if_needed(TYPEFITS, fit, TRUE);
		}

#ifdef SIRIL_OUTPUT_DEBUG
		image_find_minmax(fit);
		fprintf(stdout, "after debayer: min=%f, max=%f\n", fit->mini, fit->maxi);
		invalidate_stats_from_fit(fit);
#endif
	}

	/* we need to do something special here because it's a 1-channel sequence and
	 * image that will become 3-channel after this call, so stats caching will
	 * not be correct. Preprocessing does not compute new stats so we don't need
	 * to save them into cache now.
	 * Destroying the stats from the fit will also prevent saving them in the
	 * sequence, which may be done automatically by the caller.
	 */
	full_stats_invalidation_from_fit(fit);
	fit->history = g_slist_concat(fit->history, history);
	fit->keywords.lo = 0;
	return 0;
}

void clear_preprocessing_data(struct preprocessing_data *prepro) {
	if (prepro->bad_pixel_map_file)
		g_object_unref(prepro->bad_pixel_map_file);
	if (prepro->use_bias && prepro->bias)
		clearfits(prepro->bias);
	if (prepro->use_dark && prepro->dark)
		clearfits(prepro->dark);
	if (prepro->use_flat && prepro->flat)
		clearfits(prepro->flat);
	if (prepro->dev)
		free(prepro->dev);
}

static int prepro_finalize_hook(struct generic_seq_args *args) {
	int retval = seq_finalize_hook(args);
	clear_preprocessing_data((struct preprocessing_data *)args->user);
	free(args->user);
	return retval;
}

void start_sequence_preprocessing(struct preprocessing_data *prepro) {
	struct generic_seq_args *args = create_default_seqargs(prepro->seq);
	args->force_float = !com.pref.force_16bit && prepro->output_seqtype != SEQ_SER;
	if (!prepro->ignore_exclusion) {
		siril_debug_print("Ignoring exclusions: frames marked as excluded will still be processed.\n");
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = prepro->seq->selnum;
	}
	args->compute_size_hook = prepro_compute_size_hook;
	args->compute_mem_limits_hook = prepro_compute_mem_hook;
	args->prepare_hook = prepro_prepare_hook;
	args->image_hook = prepro_image_hook;
	args->finalize_hook = prepro_finalize_hook;
	args->description = _("Preprocessing");
	args->has_output = TRUE;
	args->new_seq_prefix = strdup(prepro->ppprefix);
	args->load_new_sequence = TRUE;
	args->force_ser_output = prepro->seq->type != SEQ_SER && prepro->output_seqtype == SEQ_SER;
	args->force_fitseq_output = prepro->seq->type != SEQ_FITSEQ && prepro->output_seqtype == SEQ_FITSEQ;
	args->user = prepro;
	// float output is always used in case of FITS sequence
	args->output_type = (args->force_ser_output || com.pref.force_16bit ||
			(args->seq->type == SEQ_SER && !args->force_fitseq_output)) ?
		DATA_USHORT : DATA_FLOAT;

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(prepro->ppprefix);
		free_generic_seq_args(args, TRUE);
	}
}

/********** SINGLE IMAGE (from com.uniq) ************/
int calibrate_single_image(struct preprocessing_data *args) {
	gchar *msg;
	fits fit = { 0 };
	int ret = 0;

	msg = g_strdup_printf(_("Pre-processing image %s"), com.uniq->filename);
	set_progress_bar_data(msg, 0.5);
	g_free(msg);

	copyfits(com.uniq->fit, &fit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
	copy_fits_metadata(com.uniq->fit, &fit);

	struct generic_seq_args generic = { .user = args };
	ret = prepro_prepare_hook(&generic);
	if (!ret)
		ret = prepro_image_hook(&generic, 0, 0, &fit, NULL, com.max_thread);
	clear_preprocessing_data(args);

	if (!ret) {
		gchar *filename = g_path_get_basename(com.uniq->filename);
		char *filename_noext = remove_ext_from_filename(filename);
		g_free(filename);
		gchar *dest_filename = g_strdup_printf("%s%s%s", args->ppprefix, filename_noext, com.pref.ext);
		msg = g_strdup_printf(_("Saving image %s"), filename_noext);
		set_progress_bar_data(msg, PROGRESS_NONE);
		ret = savefits(dest_filename, &fit);

		if (!ret) {
			// make the result the open image
			clearfits(com.uniq->fit);
			if (com.uniq->filename)
				free(com.uniq->filename);

			memcpy(com.uniq->fit, &fit, sizeof(fits));
			com.uniq->nb_layers = fit.naxes[2];
			com.uniq->filename = strdup(dest_filename);
			// this way of opening it will not create gfit.header
		}
		else clearfits(&fit);

		free(filename_noext);
		g_free(dest_filename);
		g_free(msg);
	} else {
		clearfits(&fit);
	}

	return ret;
}

/* single image preprocess, headless */
int preprocess_given_image(char *file, struct preprocessing_data *args) {
	fits fit = { 0 };
	int ret = 0;

	siril_log_message(_("Pre-processing image %s\n"), file);
	struct generic_seq_args generic = { .user = args };

	if (readfits(file, &fit, NULL, !com.pref.force_16bit)) {
		siril_log_message(_("Could not load the image, aborting.\n"));
		return 1;
	}

	ret = prepro_prepare_hook(&generic);
	if (!ret)
		ret = prepro_image_hook(&generic, 0, 0, &fit, NULL, com.max_thread);
	clear_preprocessing_data(args);

	if (!ret) {
		gchar *filename = g_path_get_basename(file);
		char *filename_noext = remove_ext_from_filename(filename);
		g_free(filename);
		gchar *dest_filename = g_strdup_printf("%s%s%s", args->ppprefix, filename_noext, com.pref.ext);
		siril_log_message(_("Saving image %s\n"), filename_noext);
		ret = savefits(dest_filename, &fit);
		free(filename_noext);
		g_free(dest_filename);
		clearfits(&fit);
	}
	return ret;
}


int evaluateoffsetlevel(const char* expression, fits *fit) {
	// try to find an occurence of *
	// If none -> the level is just an integer to evaluate
	// If found -> Try to find $ sign to read the offset value and its multiplier

	gchar *expressioncpy = g_strdup(expression);
	if (!expressioncpy)
		return 0;
	gchar *end = NULL;
	remove_spaces_from_str(expressioncpy);
	gchar *mulsignpos = g_strrstr(expressioncpy, "*");
	int offsetlevel, multiplier;
	if (!mulsignpos) {
		offsetlevel = g_ascii_strtoull(expressioncpy, NULL, 10);
		if (!offsetlevel) goto free_on_error;
		if (expressioncpy) g_free(expressioncpy);
		return offsetlevel;
	}
	mulsignpos[0] = '\0';
	mulsignpos += 1;
	if (!((expressioncpy[0] == '$') || (mulsignpos[0] == '$'))) goto free_on_error; //found a * char but none of the words start with a $
	if (expressioncpy[0] == '$') {
		multiplier = g_ascii_strtoull(mulsignpos, &end, 10);
		if (!g_str_equal(expressioncpy,"$OFFSET")) goto free_on_error;
	} else {
		multiplier = g_ascii_strtoull(expressioncpy, &end, 10);
		if (!g_str_equal(mulsignpos,"$OFFSET")) goto free_on_error;
	}
	if (!multiplier) goto free_on_error; // multiplier not parsed
	if (end[0] != '\0') goto free_on_error; // some characters were found after the multiplier
	offsetlevel = (int)(multiplier * fit->keywords.key_offset);
	if (expressioncpy) g_free(expressioncpy);
	return offsetlevel;
free_on_error:
	if (expressioncpy) g_free(expressioncpy);
	return 0;
}

gboolean check_for_cosme_file_sanity(GFile *file) {
	GError *error = NULL;
	gchar *line;

	GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, &error);

	if (input_stream == NULL) {
		if (error != NULL) {
			g_clear_error(&error);
			siril_log_message(_("File [%s] does not exist\n"), g_file_peek_path(file));
		}
		g_object_unref(file);
		return FALSE;
	}

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL,
				NULL, NULL))) {
		if (!g_str_has_prefix(line, "P ") && !g_str_has_prefix(line, "C ")
				&& !g_str_has_prefix(line, "C ")
				&& !g_str_has_prefix(line, "#")) {
			g_free(line);
			g_object_unref(data_input);
			g_object_unref(input_stream);
			return FALSE;
		}
		g_free(line);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);
	return TRUE;
}

static gboolean test_for_master_files(struct preprocessing_data *args) {
	GtkToggleButton *tbutton;
	GtkEntry *entry;
	gboolean has_error = FALSE;
	args->bias_level = FLT_MAX;
	args->bad_pixel_map_file = NULL;
	fits reffit = { 0 };
	gboolean isseq = (args->seq != NULL);
	if (isseq) {
		// loading the sequence reference image's metadata in case it's needed
		int image_to_load = sequence_find_refimage(args->seq);
		if (seq_read_frame_metadata(args->seq, image_to_load, &reffit)) {
			siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
			has_error = TRUE;
		}
	} else {
		reffit = gfit;
	}

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("useoffset_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("offsetname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			if (filename[0] == '=') { // offset is specified as a level not a file
				set_progress_bar_data(_("Checking offset level..."), PROGRESS_NONE);
				int offsetlevel = evaluateoffsetlevel(filename + 1, &gfit);
				if (!offsetlevel) {
					error = _("NOT USING OFFSET: the offset value could not be parsed");
					args->use_bias = FALSE;
				} else {
					siril_log_message(_("Synthetic offset: Level = %d\n"),offsetlevel);
					int maxlevel = (gfit.orig_bitpix == BYTE_IMG) ? UCHAR_MAX : USHRT_MAX;
					if ((offsetlevel > maxlevel) || (offsetlevel < -maxlevel) ) {   // not excluding all neg values here to allow defining a pedestal
						error = _("NOT USING OFFSET: the offset value is not consistent with image bitdepth");
						args->use_bias = FALSE;
					} else {
						args->bias_level = (float)offsetlevel;
						args->bias_level *= (gfit.orig_bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE; //converting to [0 1] to use with soper
						args->use_bias = TRUE;
					}
				}
			} else {
				int status;
				gchar *expression = path_parse(&reffit, filename, PATHPARSE_MODE_READ, &status);
				if (status) {
					error = _("NOT USING OFFSET: could not parse the expression\n");
					has_error = TRUE;
				} else {
					args->bias = calloc(1, sizeof(fits));
					set_progress_bar_data(_("Opening offset image..."), PROGRESS_NONE);
					if (!readfits(expression, args->bias, NULL, !com.pref.force_16bit)) {
						if (args->bias->naxes[2] != gfit.naxes[2]) {
							error = _("NOT USING OFFSET: number of channels is different");
						} else if (args->bias->naxes[0] != gfit.naxes[0] ||
								args->bias->naxes[1] != gfit.naxes[1]) {
							error = _("NOT USING OFFSET: image dimensions are different");
						} else {
							args->use_bias = TRUE;
							// if input is 8b, we assume 32b master needs to be rescaled
							if ((args->bias->type == DATA_FLOAT) && (gfit.orig_bitpix == BYTE_IMG)) {
								soper(args->bias, USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE, OPER_MUL, TRUE);
							}
						}
					} else
						error = _("NOT USING OFFSET: cannot open the file");
				}
				g_free(expression);
			}
			if (error) {
				siril_log_color_message("%s\n", "red", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				if (args->bias)
					free(args->bias);
				args->use_bias = FALSE;
				has_error = TRUE;
			}
		}
	}

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("usedark_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("darkname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			int status;
			gchar *expression = path_parse(&reffit, filename, PATHPARSE_MODE_READ, &status);
			if (status) {
				error = _("NOT USING DARK: could not parse the expression");
				has_error = TRUE;
			} else {
				set_progress_bar_data(_("Opening dark image..."), PROGRESS_NONE);
				args->dark = calloc(1, sizeof(fits));
				if (!readfits(expression, args->dark, NULL, !com.pref.force_16bit)) {
					if (args->dark->naxes[2] != gfit.naxes[2]) {
						error = _("NOT USING DARK: number of channels is different");
					} else if (args->dark->naxes[0] != gfit.naxes[0] ||
							args->dark->naxes[1] != gfit.naxes[1]) {
						error = _("NOT USING DARK: image dimensions are different");
					} else {
						args->use_dark = TRUE;
						// if input is 8b, we assume 32b master needs to be rescaled
						if ((args->dark->type == DATA_FLOAT) && (gfit.orig_bitpix == BYTE_IMG)) {
							soper(args->dark, USHRT_MAX_SINGLE / UCHAR_MAX_SINGLE, OPER_MUL, TRUE);
						}
					}

				} else error = _("NOT USING DARK: cannot open the file");
			}
			if (error) {
				siril_log_color_message("%s\n", "red", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				if (args->dark)
					clearfits(args->dark);
				args->dark = NULL; // in order to be sure it is freed
				args->use_dark = FALSE;
				has_error = TRUE;
			}
			g_free(expression);

		}

		if (args->use_dark) {
			// dark optimization
			int optim = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboDarkOptimize")));
			args->use_dark_optim = optim != 0;
			args->use_exposure = optim == 2;
			const char *error = NULL;
			if ((gfit.orig_bitpix == BYTE_IMG) && args->use_dark_optim) {
				error = _("Dark optimization: This process cannot be applied to 8b images");
			}
			if (error) {
				siril_log_color_message("%s\n", "red", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				clearfits(args->dark);
				gtk_entry_set_text(entry, "");
				args->use_dark_optim = FALSE;
				has_error = TRUE;
			}

			// cosmetic correction
			tbutton = GTK_TOGGLE_BUTTON(lookup_widget("cosmEnabledCheck"));
			args->use_cosmetic_correction = gtk_toggle_button_get_active(tbutton);

			if (args->use_cosmetic_correction) {
				tbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkSigCold"));
				if (gtk_toggle_button_get_active(tbutton)) {
					GtkSpinButton *sigCold = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeColdBox"));
					args->sigma[0] = gtk_spin_button_get_value(sigCold);
				} else args->sigma[0] = -1.0;

				tbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkSigHot"));
				if (gtk_toggle_button_get_active(tbutton)) {
					GtkSpinButton *sigHot = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHotBox"));
					args->sigma[1] = gtk_spin_button_get_value(sigHot);
				} else args->sigma[1] = -1.0;

				/* Using Bad Pixel Map ? */
				const gchar *bad_pixel_f;
				entry = GTK_ENTRY(lookup_widget("pixelmap_entry"));
				bad_pixel_f = gtk_entry_get_text(entry);
				/* test for file */
				if (bad_pixel_f[0] != '\0') {
					args->bad_pixel_map_file = g_file_new_for_path(bad_pixel_f);
					if (!check_for_cosme_file_sanity(args->bad_pixel_map_file)) {
						siril_debug_print("cosme file sanity check failed...\n");
						has_error = TRUE;
					}
				}
			}
		}
	}

	/* now we want to know which cosmetic correction we choose */
	GtkStack *stack = GTK_STACK(lookup_widget("stack_cc"));
	GtkWidget *w = gtk_stack_get_visible_child(stack);
	args->cc_from_dark = TRUE; // default value
	if (w) {
		GValue value = G_VALUE_INIT;
		g_value_init(&value, G_TYPE_INT);
		gtk_container_child_get_property(GTK_CONTAINER(stack), w, "position", &value);
		gint position = g_value_get_int(&value);
		args->cc_from_dark = position == 0 ? TRUE : FALSE;
		g_value_unset(&value);
	}

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("useflat_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			int status;
			gchar *expression = path_parse(&reffit, filename, PATHPARSE_MODE_READ, &status);
			if (status) {
				error = _("NOT USING FLAT: could not parse the expression");
				has_error = TRUE;
			} else {
				set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
				args->flat = calloc(1, sizeof(fits));
				if (!readfits(expression, args->flat, NULL, !com.pref.force_16bit)) {
					if (args->flat->naxes[2] != gfit.naxes[2]) {
						error = _("NOT USING FLAT: number of channels is different");
					} else if (args->flat->naxes[0] != gfit.naxes[0] ||
							args->flat->naxes[1] != gfit.naxes[1]) {
						error = _("NOT USING FLAT: image dimensions are different");
					} else {
						args->use_flat = TRUE;
						// no need to deal with bitdepth conversion as flat is just a division (unlike darks which need to be on same scale)
					}

				} else error = _("NOT USING FLAT: cannot open the file");
			}
			g_free(expression); // expression not used again after here, free before it falls out of scope
			if (error) {
				siril_log_color_message("%s\n", "red", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				if (args->flat)
					free(args->flat);
				args->use_flat = FALSE;
				has_error = TRUE;
			}

			if (args->use_flat) {
				GtkToggleButton *autobutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto_evaluate"));
				args->autolevel = gtk_toggle_button_get_active(autobutton);
				if (!args->autolevel) {
					GtkEntry *norm_entry = GTK_ENTRY(lookup_widget("entry_flat_norm"));
					args->normalisation = g_ascii_strtod(gtk_entry_get_text(norm_entry), NULL);
				}
			}
		}
	}
	if (isseq)
		clearfits(&reffit);
	return has_error;
}

void on_prepro_button_clicked(GtkButton *button, gpointer user_data) {
	if (!single_image_is_loaded() && !sequence_is_loaded())
		return;
	if (!single_image_is_loaded() && get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	GtkEntry *entry = GTK_ENTRY(lookup_widget("preproseqname_entry"));
	GtkToggleButton *fix_xtrans = GTK_TOGGLE_BUTTON(lookup_widget("fix_xtrans_af"));
	GtkToggleButton *CFA = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
	GtkToggleButton *debayer = GTK_TOGGLE_BUTTON(lookup_widget("checkButton_pp_dem"));
	GtkToggleButton *equalize_cfa = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
	GtkToggleButton *preprocess_excluded = GTK_TOGGLE_BUTTON(lookup_widget("toggle_preprocess_excluded"));
	GtkComboBox *output_type = GTK_COMBO_BOX(lookup_widget("prepro_output_type_combo"));

	struct preprocessing_data *args = calloc(1, sizeof(struct preprocessing_data));
	if (test_for_master_files(args)) {
		siril_log_color_message(_("Some errors have been detected, Please check the logs.\n"), "red");
		free(args);
		return;
	}

	siril_log_color_message(_("Preprocessing...\n"), "green");

	// set output filename (preprocessed file name prefix)
	args->ppprefix = strdup(gtk_entry_get_text(entry));

	args->ignore_exclusion = gtk_toggle_button_get_active(preprocess_excluded);
	args->is_cfa = gtk_toggle_button_get_active(CFA);
	args->debayer = gtk_toggle_button_get_active(debayer);
	args->equalize_cfa = gtk_toggle_button_get_active(equalize_cfa);
	args->fix_xtrans = gtk_toggle_button_get_active(fix_xtrans);

	/****/

	if (sequence_is_loaded()) {
		args->is_sequence = TRUE;
		args->seq = &com.seq;
		args->output_seqtype = gtk_combo_box_get_active(output_type);
		if (args->output_seqtype < 0 || args->output_seqtype > SEQ_FITSEQ)
			args->output_seqtype = SEQ_REGULAR;
		args->allow_32bit_output = !com.pref.force_16bit && args->output_seqtype != SEQ_SER;
		set_cursor_waiting(TRUE);
		control_window_switch_to_tab(OUTPUT_LOGS);
		start_sequence_preprocessing(args);
	} else {
		int retval;
		args->is_sequence = FALSE;
		args->output_seqtype = SEQ_REGULAR;
		args->allow_32bit_output = !com.pref.force_16bit;
		set_cursor_waiting(TRUE);
		control_window_switch_to_tab(OUTPUT_LOGS);

		retval = calibrate_single_image(args);

		free(args);

		if (retval)
			set_progress_bar_data(_("Error in preprocessing."), PROGRESS_NONE);
		else {
			set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
			invalidate_gfit_histogram();
			gui_function(open_single_image_from_gfit, NULL);
		}
		set_cursor_waiting(FALSE);
	}
}

void on_GtkButtonEvaluateCC_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;
	GtkLabel *label[2];
	GtkWidget *widget[2];
	const char *filename;
	gchar *str[2];
	double sig[2];
	long icold = 0L, ihot = 0L;
	double rate, total;
	fits fit = { 0 };
	fits reffit = { 0 };
	int status;
	gboolean isseq = FALSE;
	if (sequence_is_loaded()) {
		// loading the sequence reference image's metadata
		int image_to_load = sequence_find_refimage(&com.seq);
		if (seq_read_frame_metadata(&com.seq, image_to_load, &reffit)) {
			siril_log_message(_("Could not load the reference image of the sequence, aborting.\n"));
			return;
		}
		isseq = TRUE;
	} else {
		reffit = gfit;
	}

	set_cursor_waiting(TRUE);
	sig[0] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeCold")));
	sig[1] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHot")));
	widget[0] = lookup_widget("GtkLabelColdCC");
	widget[1] = lookup_widget("GtkLabelHotCC");
	label[0] = GTK_LABEL(lookup_widget("GtkLabelColdCC"));
	label[1] = GTK_LABEL(lookup_widget("GtkLabelHotCC"));
	entry = GTK_ENTRY(lookup_widget("darkname_entry"));
	filename = gtk_entry_get_text(entry);
	gchar *expression = path_parse(&reffit, filename, PATHPARSE_MODE_READ, &status);
	if (status || readfits(expression, &fit, NULL, !com.pref.force_16bit)) {
		str[0] = g_markup_printf_escaped(_("<span foreground=\"red\">ERROR</span>"));
		str[1] = g_markup_printf_escaped(_("<span foreground=\"red\">ERROR</span>"));
		gtk_label_set_markup(label[0], str[0]);
		gtk_label_set_markup(label[1], str[1]);
		set_cursor_waiting(FALSE);
		g_free(expression);
		return;
	}
	find_deviant_pixels(&fit, sig, &icold, &ihot, TRUE);
	printf("sig[0]=%lf\n", sig[0]);
	printf("sig[1]=%lf\n", sig[1]);
	total = fit.rx * fit.ry;
	clearfits(&fit);
	rate = (double)icold / total;
	/* 1% of cold pixels seems to be a reasonable limit */
	if (rate > 0.01) {
		str[0] = g_markup_printf_escaped("<span foreground=\"red\">%ld px</span>", icold);
		gtk_widget_set_tooltip_text(widget[0], _("This value may be too high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[0] = g_markup_printf_escaped("%ld px", icold);
		gtk_widget_set_tooltip_text(widget[0], "");
	}
	gtk_label_set_markup(label[0], str[0]);

	rate = (double)ihot / total;
	/* 1% of hot pixels seems to be a reasonable limit */
	if (rate > 0.01) {
		str[1] = g_markup_printf_escaped("<span foreground=\"red\">%ld px</span>", ihot);
		gtk_widget_set_tooltip_text(widget[1], _("This value may be too high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[1] = g_markup_printf_escaped("%ld px", ihot);
		gtk_widget_set_tooltip_text(widget[1], "");
	}
	gtk_label_set_markup(label[1], str[1]);
	g_free(str[0]);
	g_free(str[1]);
	set_cursor_waiting(FALSE);
	g_free(expression);
	if (isseq)
		clearfits(&reffit);
}
