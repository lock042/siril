/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_statistics_ushort.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "core/optimize_utils.h"
#include "core/siril_log.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "algos/median_fast.h"
#include "filters/median.h"
#include "opencv/opencv.h"

#include "cosmetic_correction.h"
static int autoDetect(fits *fit, int layer, const double sig[2], long *icold, long *ihot,
		double amount, gboolean is_cfa, threading_type threads);

/* see also getMedian3x3 in algos/PSF.c */
static float getMedian5x5_float(const float *buf, const int xx, const int yy, const int w,
		const int h, gboolean is_cfa) {

	const int step = is_cfa ? 2 : 1;
	const int radius = 2 * step;

	int n = 0;
	float value[24];
	for (int y = yy - radius; y <= yy + radius; y += step) {
		if (y >= 0 && y < h) {
			for (int x = xx - radius; x <= xx + radius; x += step) {
				if (x >= 0 && x < w) {
					// ^ limit to image bounds ^
					// exclude centre pixel v
					if (x != xx || y != yy) {
						value[n++] = buf[x + y * w];
					}
				}
			}
		}
	}
	return quickmedian_float(value, n);
}

static WORD* getAverage3x3Line(WORD *buf, const int yy, const int w,
		const int h, gboolean is_cfa) {
	g_assert(w > 0 && h > 0);
	int step, radius, x, xx, y;
	WORD *cpyline;

	if (is_cfa)
		step = radius = 2;
	else
		step = radius = 1;

	cpyline = calloc(w, sizeof(WORD));
	for (xx = 0; xx < w; ++xx) {
		int n = 0;
		float value = 0.f;
		for (y = yy - radius; y <= yy + radius; y += step) {
			if (y != yy) {	// we skip the line
				for (x = xx - radius; x <= xx + radius; x += step) {
					if (y >= 0 && y < h && x >= 0 && x < w) {
						value += (float) buf[x + y * w];
						n++;
					}
				}
			}
		}
		cpyline[xx] = n == 0 ? 0 : round_to_WORD(value / n);
	}
	return cpyline;
}

static float* getAverage3x3Line_float(const float *buf, const int yy, const int w,
		const int h, gboolean is_cfa) {
	g_assert(w > 0 && h > 0);
	int step, radius, x, xx, y;
	float *cpyline;

	if (is_cfa)
		step = radius = 2;
	else
		step = radius = 1;

	cpyline = calloc(w, sizeof(float));
	for (xx = 0; xx < w; ++xx) {
		int n = 0;
		float value = 0.f;
		for (y = yy - radius; y <= yy + radius; y += step) {
			if (y != yy) {	// we skip the line
				for (x = xx - radius; x <= xx + radius; x += step) {
					if (y >= 0 && y < h && x >= 0 && x < w) {
						value += buf[x + y * w];
						n++;
					}
				}
			}
		}
		cpyline[xx] = n == 0 ? 0.f : (value / n);
	}
	return cpyline;
}

static float getAverage3x3_float(const float *buf, const int xx, const int yy,
		const int w, const int h, gboolean is_cfa) {
	g_assert(w > 0 && h > 0);

    const int step = is_cfa ? 2 : 1;
    const int radius = step;

    int n = -1;
    float value = -buf[xx + yy * w];
    for (int y = yy - radius; y <= yy + radius; y += step) {
        if (y >= 0 && y < h) {
            for (int x = xx - radius; x <= xx + radius; x += step) {
                if (x >= 0 && x < w) {
                    value += buf[x + y * w];
                    n++;
                }
            }
        }
    }
    return n == 0 ? 0.f : value / n;
}

static float getAverage3x3_ushort(WORD *buf, const int xx, const int yy,
		const int w, const int h, gboolean is_cfa) {
	g_assert(w > 0 && h > 0);
	int step, radius, x, y;
	float value = 0.f;

	if (is_cfa)
		step = radius = 2;
	else
		step = radius = 1;

	int n = 0;
	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h) {
				if (x >= 0 && x < w) {
					if ((x != xx) || (y != yy)) {
						value += (float) buf[x + y * w];
						n++;
					}
				}
			}
		}
	}
	return n == 0 ? 0 :value / n;
}

/* Gives a list of point p containing deviant pixel coordinates, to be freed by
 * caller.
 * If eval_only is true, the function only counts the deviant pixels and
 * returns NULL. It also returns NULL when no deviant pixel is found.
 * If cold == -1 or hot == -1, this is a flag to not compute cold or hot
 */
deviant_pixel* find_deviant_pixels(fits *fit, const double sig[2], long *icold,
		long *ihot, gboolean eval_only) {
	int x, y;
	WORD *buf;
	float *fbuf;
	imstats *stat;
	double sigma, median;
	float thresHot, thresCold;
	deviant_pixel *dev;

	if (sig[0] == -1.0 && sig[1] == -1.0)
		return NULL;

	stat = statistics(NULL, -1, fit, RLAYER, NULL, STATS_BASIC, SINGLE_THREADED);
	if (!stat) {
		siril_log_message(_("Error: statistics computation failed.\n"));
		return NULL;
	}
	sigma = stat->sigma;
	median = stat->median;

	if (sig[0] == -1.0) {	// flag for no cold detection
		thresCold = -1.f;
	} else {
		double val = median - (sig[0] * sigma);
		thresCold = max((float) val, 0.f);
	}
	if (sig[1] == -1.0) {	// flag for no hot detection
		thresHot = USHRT_MAX_SINGLE + 1.f;
	} else {
		double val = median + (sig[1] * sigma);
		thresHot = min((float) val, fit->type == DATA_FLOAT ? 1.f : USHRT_MAX_SINGLE);
	}

	free_stats(stat);

	/** First we count deviant pixels **/
	*icold = 0;
	*ihot = 0;
	buf = fit->pdata[RLAYER];
	fbuf = fit->fpdata[RLAYER];
	size_t i, nbpix = fit->naxes[0] * fit->naxes[1];
	for (i = 0; i < nbpix; i++) {
		float pixel = fit->type == DATA_FLOAT ? fbuf[i] : (float) buf[i];
		if (pixel >= thresHot)
			(*ihot)++;
		else if (pixel <= thresCold)
			(*icold)++;
	}

	if (eval_only) return NULL;

	/** Second we store deviant pixels in dev */
	int n = (*icold) + (*ihot);
	if (n <= 0)
		return NULL;
	dev = malloc(n * sizeof(deviant_pixel));
	if (!dev) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	i = 0;
	for (y = 0; y < fit->ry; y++) {
		for (x = 0; x < fit->rx; x++) {
			float pixel = fit->type == DATA_FLOAT ?
							fbuf[x + y * fit->rx] : (float) buf[x + y * fit->rx];
			if (pixel >= thresHot) {
				dev[i].p.x = x;
				dev[i].p.y = y;
				dev[i].type = HOT_PIXEL;
				i++;
			} else if (pixel <= thresCold) {
				dev[i].p.x = x;
				dev[i].p.y = y;
				dev[i].type = COLD_PIXEL;
				i++;
			}
		}
	}
	return dev;
}

int cosmeticCorrOnePoint(fits *fit, deviant_pixel dev, gboolean is_cfa) {
	// Cosmetic correction, as developed here, is only used on 1-channel images
	int width = fit->rx;
	int height = fit->ry;
	int x = (int) dev.p.x;
	int y = (int) dev.p.y;

	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->pdata[RLAYER];
		WORD newpixel;
		if (dev.type == COLD_PIXEL)
			newpixel = get_median_ushort(buf, x, y, width, height, 2, is_cfa,
					FALSE);
		else
			newpixel = round_to_WORD(getAverage3x3_ushort(buf, x, y, width, height, is_cfa));
		buf[x + y * fit->rx] = newpixel;
	} else if (fit->type == DATA_FLOAT) {
		float *buf = fit->fpdata[RLAYER];
		float newpixel;
		if (dev.type == COLD_PIXEL)
			newpixel = get_median_float(buf, x, y, width, height, 2, is_cfa,
					FALSE);
		else
			newpixel = getAverage3x3_float(buf, x, y, width, height, is_cfa);
		buf[x + y * fit->rx] = newpixel;
	}

	// the caller should call invalidate_stats_from_fit(fit);
	return 0;
}

int cosmeticCorrOneLine(fits *fit, deviant_pixel dev, gboolean is_cfa) {
	if (fit->type == DATA_FLOAT) {
		float *buf = fit->fpdata[RLAYER];
		float *line, *newline;
		int width = fit->rx;
		int height = fit->ry;
		int row = (int) dev.p.y;

		line = buf + row * width;
		newline = getAverage3x3Line_float(buf, row, width, height, is_cfa);
		memcpy(line, newline, width * sizeof(float));

		free(newline);
		//invalidate_stats_from_fit(fit);
		return 0;
	} else if (fit->type == DATA_USHORT) {
		WORD *buf = fit->pdata[RLAYER];
		WORD *line, *newline;
		int width = fit->rx;
		int height = fit->ry;
		int row = (int) dev.p.y;

		line = buf + row * width;
		newline = getAverage3x3Line(buf, row, width, height, is_cfa);
		memcpy(line, newline, width * sizeof(WORD));

		free(newline);
		//invalidate_stats_from_fit(fit);
		return 0;
	}
	return 1;
}

int cosmeticCorrection(fits *fit, deviant_pixel *dev, int size, gboolean is_cfa) {
	for (int i = 0; i < size; i++) {
		cosmeticCorrOnePoint(fit, dev[i], is_cfa);
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/**** Autodetect *****/
static int cosmetic_finalize_hook(struct generic_seq_args *args) {
	int retval = seq_finalize_hook(args);
	free(args->user);
	return retval;
}

int cosmetic_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct cosmetic_data *c_args = (struct cosmetic_data*) args->user;
	int chan;
	/* Count variables, icold and ihot, need to be local in order to be parallelized */
	long icold, ihot;

	icold = ihot = 0L;
	for (chan = 0; chan < fit->naxes[2]; chan++) {
		int retval = autoDetect(fit, chan, c_args->sigma, &icold, &ihot,
				c_args->amount, c_args->is_cfa, c_args->threading);
		if (retval)
			return retval;
	}
	siril_log_color_message(_("Image %d: %ld pixel corrected (%ld + %ld)\n"),
			"bold", i, icold + ihot, icold, ihot);
	return 0;
}

static int cosmetic_mem_limits_hook(struct generic_seq_args *args, gboolean for_writer) {
	/* stats O(m)
	 * working on image channel as float O(m as float)
	 */
	unsigned int MB_per_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_image, NULL, &MB_avail);
	unsigned int required = MB_per_image;
	if (limit > 0) {
		int is_color = args->seq->nb_layers == 3;
		int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;
		unsigned int float_multiplier = is_float ? 1 : 2;
		unsigned int MB_per_float_image = MB_per_image * float_multiplier;
		unsigned int MB_per_float_channel = is_color ? MB_per_float_image / 3 : MB_per_float_image;
		required += MB_per_float_channel;
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
                        thread_limit = com.max_thread;

		if (for_writer) {
                        /* we allow the already allocated thread_limit images,
                         * plus how many images can be stored in what remains
                         * unused by the main processing */
                        limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_image;
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
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

void apply_cosmetic_to_sequence(struct cosmetic_data *cosme_args) {
	struct generic_seq_args *args = create_default_seqargs(cosme_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = cosme_args->seq->selnum;
	args->compute_mem_limits_hook = cosmetic_mem_limits_hook;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = cosmetic_finalize_hook;
	args->image_hook = cosmetic_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Cosmetic Correction");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = strdup(cosme_args->seqEntry);
	args->load_new_sequence = TRUE;
	args->user = cosme_args;

	cosme_args->fit = NULL;	// not used here

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(cosme_args->seqEntry);
		free(cosme_args);
		free_generic_seq_args(args, TRUE);
	}
}

// idle function executed at the end of the Cosmetic Correction processing
gboolean end_autoDetect(gpointer p) {
	stop_processing_thread();
	notify_gfit_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);

	return FALSE;
}

gpointer autoDetectThreaded(gpointer p) {
	struct cosmetic_data *args = (struct cosmetic_data*) p;
	struct timeval t_start, t_end;
	int retval = 0, chan;
	long icold, ihot;

	siril_log_color_message(_("Cosmetic Correction: processing...\n"), "green");
	gettimeofday(&t_start, NULL);

	icold = ihot = 0L;
	for (chan = 0; chan < args->fit->naxes[2]; chan++) {
		retval = autoDetect(args->fit, chan, args->sigma, &icold, &ihot,
				args->amount, args->is_cfa, args->threading);
		if (retval)
			break;
	}
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	gchar *str = ngettext("%ld corrected pixel (%ld + %ld)\n", "%ld corrected pixels (%ld + %ld)\n", icold + ihot);
	str = g_strdup_printf(str, icold + ihot, icold, ihot);
	siril_log_message(str);
	g_free(str);

	free(args);
	siril_add_idle(end_autoDetect, NULL);
	return GINT_TO_POINTER(retval);
}

int apply_cosme_to_image(fits *fit, GFile *file, int is_cfa) {
	int i = 0;
	gchar *line;
	deviant_pixel dev;
	int nb_tokens, retval = 0;
	double dirty;
	char type;
	GError *error = NULL;
	GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, &error);

	if (input_stream == NULL) {
		if (error != NULL) {
			g_clear_error(&error);
			siril_log_message(_("File [%s] does not exist\n"), g_file_peek_path(file));
		}
		return 1;
	}

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL,
				NULL, NULL))) {
		++i;
		switch (line[0]) {
		case '#': // comments.
			g_free(line);
			continue;
			break;
		case 'P':
			nb_tokens = sscanf(line + 2, "%lf %lf %c", &dev.p.x, &dev.p.y, &type);
			if (nb_tokens != 2 && nb_tokens != 3) {
				fprintf(stderr, "cosmetic correction: "
						"cosme file format error at line %d: %s", i, line);
				retval = 1;
				g_free(line);
				continue;
			}
			if (nb_tokens == 2)	{
				type = 'H';
			}
			if (type == 'H')
				dev.type = HOT_PIXEL;
			else
				dev.type = COLD_PIXEL;
			dev.p.y = fit->ry - dev.p.y - 1;  /* FITS are stored bottom to top */
			cosmeticCorrOnePoint(fit, dev, is_cfa);
			break;
		case 'L':
			nb_tokens = sscanf(line + 2, "%lf %lf %c", &dev.p.y, &dirty, &type);
			if (nb_tokens != 2 && nb_tokens != 3) {
				fprintf(stderr, "cosmetic correction: "
						"cosme file format error at line %d: %s\n", i, line);
				retval = 1;
				g_free(line);
				continue;
			}
			dev.type = HOT_PIXEL; // we force it
			dev.p.y = fit->ry - dev.p.y - 1; /* FITS are stored bottom to top */
			cosmeticCorrOneLine(fit, dev, is_cfa);
			break;
		case 'C':
			nb_tokens = sscanf(line + 2, "%lf %lf %c", &dev.p.y, &dirty, &type);
			if (nb_tokens != 2 && nb_tokens != 3) {
				fprintf(stderr, "cosmetic correction: "
						"cosme file format error at line %d: %s\n", i, line);
				retval = 1;
				g_free(line);
				continue;
			}
			dev.type = HOT_PIXEL; // we force it
			dev.p.y = fit->rx - dev.p.y - 1; /* FITS are stored bottom to top */
			if (cvRotateImage(fit, 90)) {
				retval = 1;
				break;
			}
			cosmeticCorrOneLine(fit, dev, is_cfa);
			if (cvRotateImage(fit, -90)) {
				retval = 1;
				break;
			}
			break;
		default:
			fprintf(stderr, _("cosmetic correction: "
					"cosme file format error at line %d: %s\n"), i, line);
			retval = 1;
		}
		g_free(line);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);

	return retval;
}

int cosme_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct cosme_data *c_args = (struct cosme_data*) args->user;

	int retval = apply_cosme_to_image(fit, c_args->file, c_args->is_cfa);
	g_object_unref(c_args->file);
	return retval;
}

static int cosme_finalize_hook(struct generic_seq_args *args) {
	int retval = seq_finalize_hook(args);
	struct cosme_data *c_args = (struct cosme_data*) args->user;
	free(c_args->prefix);
	g_clear_object(&c_args->file);
	free(args->user);
	return retval;
}

void apply_cosme_to_sequence(struct cosme_data *cosme_args) {
	struct generic_seq_args *args = create_default_seqargs(cosme_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = cosme_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = cosme_finalize_hook;
	args->image_hook = cosme_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Cosmetic Correction");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = strdup(cosme_args->prefix);
	args->load_new_sequence = TRUE;
	args->user = cosme_args;

	cosme_args->fit = NULL;	// not used here

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		g_clear_object(&cosme_args->file);
		free(cosme_args->prefix);
		free(cosme_args);
		free_generic_seq_args(args, TRUE);
	}
}

/* this is an autodetect algorithm. Cold and hot pixels
 *  are corrected in the same time */
static int autoDetect(fits *fit, int layer, const double sig[2], long *icold, long *ihot,
		double amount, gboolean is_cfa, threading_type threading) {

	/* XXX: if cfa, stats are irrelevant. We should compute them taking
	 * into account the Bayer pattern */
	imstats *stat = statistics(NULL, -1, fit, layer, NULL, STATS_BASIC | STATS_AVGDEV, threading);
	if (!stat) {
		siril_log_message(_("Error: statistics computation failed.\n"));
		return 1;
	}
	int threads = check_threading(&threading);
	const float bkg = stat->median;
	const float avgDev = stat->avgDev;
	free_stats(stat);
	const int width = fit->rx;
	const int height = fit->ry;
	const float f0 = amount;
	const float f1 = 1 - f0;

	const gboolean isFloat = fit->type == DATA_FLOAT;
	WORD *buf = fit->pdata[layer];
	float *fbuf = fit->fpdata[layer];
	const float k1 = avgDev;
	const float k2 = k1 / 2;
	const float k3 = sig[1] * k1;
	const float k4 = max(k1, k3);
	const float k = avgDev * sig[0];
	const gboolean doHot = sig[1] != -1.0;
	const gboolean doCold = sig[0] != -1.0;
	const float coldVal = doCold ? bkg - k : 0.0;
	const float hotVal = doHot ? bkg + k1 : isFloat ? 1.f : 65535.f;
	size_t n = fit->naxes[0] * fit->naxes[1];
	float *temp = malloc(n * sizeof(float));
	if (!temp) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	if (isFloat) {
		if (threads > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
			for (size_t i = 0; i < n; i++) {
				temp[i] = fbuf[i];
			}
		} else {
			// this should be faster in a single-threaded case
			memcpy(temp, fbuf, width * height * sizeof(float));
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) if(threads > 1)
#endif
		for (size_t i = 0; i < n; i++) {
			temp[i] = (float)buf[i];
		}
	}
	const int step = is_cfa ? 2 : 1;
	const int radius = 2 * step;

	long icoldL = *icold;
	long ihotL = *ihot;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) if(threads > 1)
#endif
	{
		float medianin[24];
#ifdef _OPENMP
#pragma omp for reduction(+:icoldL, ihotL) schedule(dynamic,16)
#endif
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				const float pixel = temp[x + y * width];
				if (!inInterval(pixel, coldVal, hotVal)) {
					float m;
					if (y >= radius && y < height - radius && x >= radius && x < width - radius) {
						// use fast median24 network
						int nbm = 0;
						for (int i = -radius; i <= radius; i += step) {
							for (int j = -radius; j <= radius; j += step) {
								if (i != 0 || j != 0) {
									medianin[nbm++] = temp[(y + i) * width + x + j];
								}
							}
						}
						m = median24(medianin);
					} else {
						m = getMedian5x5_float(temp, x, y, width, height, is_cfa);
					}

					/* Hot autodetect */
					if (doHot && pixel > hotVal && pixel > m + k4) {
						const float a = getAverage3x3_float(temp, x, y, width, height, is_cfa);
						if (a < m + k2) {
							ihotL++;
							if (isFloat) {
								fbuf[x + y * width] = a * f0 + pixel * f1;
							} else {
								buf[x + y * width] = a * f0 + pixel * f1;
							}
						}
					} else if (doCold && pixel < coldVal && pixel + k < m) {
						/* Cold autodetect */
						icoldL++;
						if (isFloat) {
							fbuf[x + y * width] = m * f0 + pixel * f1;
						} else {
							buf[x + y * width] = m * f0 + pixel * f1;
						}
					}
				}
			}
		}
	}
	(*icold) = icoldL;
	(*ihot) = ihotL;
	free(temp);

	invalidate_stats_from_fit(fit);
	return 0;
}

void on_button_cosmetic_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("cosmetic_dialog");
}

gboolean cosmetic_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("cosmetic_dialog");
	return TRUE;
}

void on_checkSigCosme_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	static GtkWidget *cosmeticApply = NULL;
	static GtkToggleButton *checkCosmeSigCold = NULL;
	static GtkToggleButton *checkCosmeSigHot = NULL;
	gboolean checkCold, checkHot;

	if (cosmeticApply == NULL) {
		cosmeticApply = lookup_widget("button_cosmetic_ok");
		checkCosmeSigCold = GTK_TOGGLE_BUTTON(lookup_widget("checkSigColdBox"));
		checkCosmeSigHot = GTK_TOGGLE_BUTTON(lookup_widget("checkSigHotBox"));
	}
	checkCold = gtk_toggle_button_get_active(checkCosmeSigCold);
	checkHot = gtk_toggle_button_get_active(checkCosmeSigHot);
	gtk_widget_set_sensitive(cosmeticApply, checkCold || checkHot);
}

void on_button_cosmetic_ok_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *cosmeticSeqEntry;
	GtkToggleButton *CFA, *seq;
	GtkSpinButton *sigma[2];
	GtkAdjustment *adjCosmeAmount;

	CFA = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheckBox"));
	sigma[0] = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeColdBox"));
	sigma[1] = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHotBox"));
	seq = GTK_TOGGLE_BUTTON(lookup_widget("checkCosmeticSeq"));
	cosmeticSeqEntry = GTK_ENTRY(lookup_widget("entryCosmeticSeq"));
	adjCosmeAmount = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjCosmeAmount"));

	struct cosmetic_data *args = calloc(1, sizeof(struct cosmetic_data));

	if (gtk_toggle_button_get_active(
				GTK_TOGGLE_BUTTON(lookup_widget("checkSigColdBox"))))
		args->sigma[0] = gtk_spin_button_get_value(sigma[0]);
	else
		args->sigma[0] = -1.0;

	if (gtk_toggle_button_get_active(
				GTK_TOGGLE_BUTTON(lookup_widget("checkSigHotBox"))))
		args->sigma[1] = gtk_spin_button_get_value(sigma[1]);
	else
		args->sigma[1] = -1.0;

	args->is_cfa = gtk_toggle_button_get_active(CFA);
	args->amount = gtk_adjustment_get_value(adjCosmeAmount);

	args->fit = &gfit;
	args->seqEntry = strdup(gtk_entry_get_text(cosmeticSeqEntry));
	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		if (args->seqEntry && args->seqEntry[0] == '\0') {
			free(args->seqEntry);
			args->seqEntry = strdup("cc_");
		}
		args->seq = &com.seq;
		args->threading = SINGLE_THREADED;
		gtk_toggle_button_set_active(seq, FALSE);
		apply_cosmetic_to_sequence(args);
	} else {
		args->threading = MULTI_THREADED;
		undo_save_state(&gfit, _("Cosmetic Correction"));
		if (!start_in_new_thread(autoDetectThreaded, args)) {
			free(args->seqEntry);
			free(args);
		}
	}
}

int denoise_hook_cosmetic(fits *fit) {
	long icold = 0, ihot = 0;
	double sig[2] = { 3.0, 3.0 };
	for (size_t layer = 0; layer < fit->naxes[2] ; layer++) {
		autoDetect(fit, layer, sig, &icold, &ihot, 1.0, FALSE, MULTI_THREADED);
	}
	return 0;
}
