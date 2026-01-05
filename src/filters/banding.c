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

/* ** Reduces Banding in Canon DSLR images. **
 * This code originates from CanonBandingReduction.js v0.9.1, a script
 * of PixInsight, originally written by Georg Viehoever and
 * distributed under the terms of the GNU General Public License
 */

#include <float.h>
#include <string.h>
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "opencv/opencv.h"

#include "banding.h"

static int BandingEngine(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation, threading_type threading);

/*****************************************************************************
 *      B A N D I N G      A L L O C A T O R   A N D   D E S T R U C T O R  *
 ****************************************************************************/

/* Allocator for banding_data */
struct banding_data *new_banding_data() {
	struct banding_data *args = calloc(1, sizeof(struct banding_data));
	if (args) {
		args->destroy_fn = free_banding_data;
	}
	return args;
}

/* Destructor for banding_data */
void free_banding_data(void *ptr) {
	struct banding_data *args = (struct banding_data *)ptr;
	if (!args)
		return;

	if (args->seqEntry) {
		free(args->seqEntry);
		args->seqEntry = NULL;
	}
	if (args->seq) {
		free_sequence(args->seq, TRUE);
		args->seq = NULL;
	}
	if (args->fit) {
		clearfits(args->fit);
		free(args->fit);
		args->fit = NULL;
	}
	free(ptr);
}

/*****************************************************************************
 *      B A N D I N G      R E D U C T I O N      M A N A G E M E N T        *
 ****************************************************************************/

/* Hook for sequence processing - uses generic_seq_args */
int banding_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	struct banding_data *banding_args = (struct banding_data *)args->user;
	return BandingEngine(fit, banding_args->sigma, banding_args->amount,
			banding_args->protect_highlights, banding_args->applyRotation, SINGLE_THREADED);
}

/* Hook for single image processing - uses generic_img_args */
int banding_single_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct banding_data *params = (struct banding_data *)args->user;
	if (!params)
		return 1;

	return BandingEngine(fit, params->sigma, params->amount,
			params->protect_highlights, params->applyRotation, MULTI_THREADED);
}

/* Idle function for single image processing */
static gboolean banding_single_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	if (args->retval == 0) {
		notify_gfit_modified();
	}

	// Free using the generic cleanup which will call the destructor
	free_generic_img_args(args);

	stop_processing_thread();
	return FALSE;
}

gchar *banding_log_hook(gpointer p, log_hook_detail detail) {
	struct banding_data *params = (struct banding_data*) p;
	gchar *message = NULL;
	if (!params->protect_highlights) {
		message=g_strdup_printf(_("Canon Banding Reduction (amount=%.2lf, invsigma=%.2lf)"), params->amount, params->sigma);
	} else {
		message=g_strdup_printf(_("Canon Banding Reduction (amount=%.2lf, Protect=TRUE, invsigma=%.2lf)"),
				params->amount, params->sigma);
	}
	return message;
}

static int banding_mem_limits_hook(struct generic_seq_args *args, gboolean for_writer) {
	/* [rotation => O(2n)]
	 * new image -> O(2n)
	 * + stats MAD per channel -> O(1m)
	 * [rotation => O(2n)]
	 */
	unsigned int MB_per_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_image, NULL, &MB_avail);
	unsigned int required = MB_per_image;
	if (limit > 0) {
		int is_color = args->seq->nb_layers == 3;
		unsigned int MB_per_channel = is_color ? MB_per_image / 3 : MB_per_image;
		required = 2 * MB_per_image + MB_per_channel;
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

int banding_finalize_hook(struct generic_seq_args *args) {
	struct banding_data *data = (struct banding_data *) args->user;
	int retval = seq_finalize_hook(args);
	// Note: we only free seqEntry here, not the whole struct
	// The struct itself will be freed by the sequence worker
	if (data->seqEntry) {
		free(data->seqEntry);
		data->seqEntry = NULL;
	}
	free(data);
	return retval;
}

void apply_banding_to_sequence(struct banding_data *banding_args) {
	struct generic_seq_args *args = create_default_seqargs(banding_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = args->seq->selnum;
	args->compute_mem_limits_hook = banding_mem_limits_hook;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = banding_finalize_hook;
	args->image_hook = banding_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Banding Reduction");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = strdup(banding_args->seqEntry);
	args->load_new_sequence = TRUE;
	args->user = banding_args;

	banding_args->fit = NULL;	// not used here

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(banding_args->seqEntry);
		free(banding_args);
		free_generic_seq_args(args, TRUE);
	}
}

static int fmul_layer_ushort(fits *a, int layer, float coeff) {
	WORD *buf;
	size_t i, n = a->naxes[0] * a->naxes[1];

	if (coeff < 0.0)
		return 1;
	buf = a->pdata[layer];
	for (i = 0; i < n; ++i) {
		buf[i] = round_to_WORD(buf[i] * coeff);
	}
	invalidate_stats_from_fit(a);
	return 0;
}

static int fmul_layer_float(fits *a, int layer, float coeff) {
	float *buf;
	size_t i, n = a->naxes[0] * a->naxes[1];

	if (coeff < 0.0)
		return 1;
	buf = a->fpdata[layer];
	for (i = 0; i < n; ++i) {
		buf[i] = buf[i] * coeff;
	}
	invalidate_stats_from_fit(a);
	return 0;
}

static int BandingEngine_ushort(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation, threading_type threads) {
	int chan, row, i, ret = 0;
	WORD *line, *fixline;
	double minimum = DBL_MAX, globalsigma = 0.0;
	fits *fiximage = NULL;
	double invsigma = 1.0 / sigma;

	if (applyRotation) {
		if (cvRotateImage(fit, 90)) return 1;
	}

	if (new_fit_image(&fiximage, fit->rx, fit->ry, fit->naxes[2], DATA_USHORT))
		return 1;

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, NULL, STATS_BASIC | STATS_MAD, threads);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			clearfits(fiximage);
			return 1;
		}
		double background = stat->median;
		double *rowvalue = calloc(fit->ry, sizeof(double));
		if (rowvalue == NULL) {
			PRINT_ALLOC_ERR;
			clearfits(fiximage);
			free_stats(stat);
			return 1;
		}
		if (protect_highlights) {
			globalsigma = stat->mad * MAD_NORM;
		}
		free_stats(stat);
		for (row = 0; row < fit->ry; row++) {
			line = fit->pdata[chan] + row * fit->rx;
			WORD *cpyline = calloc(fit->rx, sizeof(WORD));
			if (cpyline == NULL) {
				PRINT_ALLOC_ERR;
				free(rowvalue);
				clearfits(fiximage);
				return 1;
			}
			memcpy(cpyline, line, fit->rx * sizeof(WORD));
			int n = fit->rx;
			double median;
			if (protect_highlights) {
				quicksort_s(cpyline, n);
				WORD reject = round_to_WORD(
						background + invsigma * globalsigma);
				for (i = fit->rx - 1; i >= 0; i--) {
					if (cpyline[i] < reject)
						break;
					n--;
				}
				median = gsl_stats_ushort_median_from_sorted_data(cpyline, 1, n);
			} else {
				median = round_to_WORD(quickmedian(cpyline, n));
			}

			rowvalue[row] = background - median;
			minimum = min(minimum, rowvalue[row]);
			free(cpyline);
		}
		for (row = 0; row < fit->ry; row++) {
			fixline = fiximage->pdata[chan] + row * fiximage->rx;
			for (i = 0; i < fit->rx; i++)
				fixline[i] = round_to_WORD(rowvalue[row] - minimum);
		}
		free(rowvalue);
	}
	for (chan = 0; chan < fit->naxes[2]; chan++)
		fmul_layer_ushort(fiximage, chan, amount);
	ret = imoper(fit, fiximage, OPER_ADD, FALSE);

	invalidate_stats_from_fit(fit);
	clearfits(fiximage);
	if ((!ret) && applyRotation) {
		if (cvRotateImage(fit, -90)) return 1;
	}

	return ret;
}

static int BandingEngine_float(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation, threading_type threads) {
	int chan, row, i, ret = 0;
	float *line, *fixline;
	double minimum = DBL_MAX, globalsigma = 0.0;
	fits *fiximage = NULL;
	double invsigma = 1.0 / sigma;

	if (applyRotation) {
		if (cvRotateImage(fit, 90)) return 1;
	}

	if (new_fit_image(&fiximage, fit->rx, fit->ry, fit->naxes[2], DATA_FLOAT))
		return 1;

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, NULL, STATS_BASIC | STATS_MAD, threads);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		double background = stat->median;
		double *rowvalue = calloc(fit->ry, sizeof(double));
		if (rowvalue == NULL) {
			PRINT_ALLOC_ERR;
			free_stats(stat);
			return 1;
		}
		if (protect_highlights) {
			globalsigma = stat->mad * MAD_NORM;
		}
		free_stats(stat);
		for (row = 0; row < fit->ry; row++) {
			line = fit->fpdata[chan] + row * fit->rx;
			float *cpyline = calloc(fit->rx, sizeof(float));
			if (cpyline == NULL) {
				PRINT_ALLOC_ERR;
				free(rowvalue);
				return 1;
			}
			memcpy(cpyline, line, fit->rx * sizeof(float));
			int n = fit->rx;
			double median;
			if (protect_highlights) {
				quicksort_f(cpyline, n);
				float reject = background + invsigma * globalsigma;
				for (i = fit->rx - 1; i >= 0; i--) {
					if (cpyline[i] < reject)
						break;
					n--;
				}
				median = gsl_stats_float_median_from_sorted_data(cpyline, 1, n);
			} else {
				median = quickmedian_float(cpyline, n);
			}

			rowvalue[row] = background - median;
			minimum = min(minimum, rowvalue[row]);
			free(cpyline);
		}
		for (row = 0; row < fit->ry; row++) {
			fixline = fiximage->fpdata[chan] + row * fiximage->rx;
			for (i = 0; i < fit->rx; i++)
				fixline[i] = rowvalue[row] - minimum;
		}
		free(rowvalue);
	}
	for (chan = 0; chan < fit->naxes[2]; chan++)
		fmul_layer_float(fiximage, chan, amount);
	ret = imoper(fit, fiximage, OPER_ADD, TRUE);

	invalidate_stats_from_fit(fit);
	clearfits(fiximage);
	free(fiximage);
	if ((!ret) && applyRotation) {
		if (cvRotateImage(fit, -90)) return 1;
	}

	return ret;
}

static int BandingEngine(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation, threading_type threading) {
	int threads = check_threading(&threading);

	if (fit->type == DATA_FLOAT)
		return BandingEngine_float(fit, sigma, amount, protect_highlights, applyRotation, threads);
	if (fit->type == DATA_USHORT)
		return BandingEngine_ushort(fit, sigma, amount, protect_highlights, applyRotation, threads);
	return -1;
}

/***************** GUI for Canon Banding Reduction ********************/

void on_button_ok_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("canon_fixbanding_dialog");
}

gboolean banding_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("canon_fixbanding_dialog");
	return TRUE;
}

void on_button_apply_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	static GtkRange *range_amount = NULL;
	static GtkRange *range_invsigma = NULL;
	static GtkToggleButton *toggle_protect_highlights_banding = NULL,
		*vertical = NULL, *seq = NULL;
	static GtkEntry *bandingSeqEntry = NULL;
	double amount, invsigma;
	gboolean protect_highlights;

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (range_amount == NULL) {
		range_amount = GTK_RANGE(lookup_widget("scale_fixbanding_amount"));
		range_invsigma = GTK_RANGE(lookup_widget("scale_fixbanding_invsigma"));
		toggle_protect_highlights_banding = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_fixbanding"));
		vertical = GTK_TOGGLE_BUTTON(lookup_widget("checkBandingVertical"));
		seq = GTK_TOGGLE_BUTTON(lookup_widget("checkBandingSeq"));
		bandingSeqEntry = GTK_ENTRY(lookup_widget("entryBandingSeq"));
	}
	amount = gtk_range_get_value(range_amount);
	invsigma = gtk_range_get_value(range_invsigma);
	protect_highlights = gtk_toggle_button_get_active(
			toggle_protect_highlights_banding);
	gboolean applyRotation = gtk_toggle_button_get_active(vertical);

	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		// Sequence processing
		struct banding_data *seq_args = new_banding_data();
		if (!seq_args) {
			PRINT_ALLOC_ERR;
			set_cursor_waiting(FALSE);
			return;
		}

		const char *entry_text = gtk_entry_get_text(bandingSeqEntry);
		seq_args->seqEntry = strdup((entry_text && entry_text[0] != '\0') ? entry_text : "unband_");
		seq_args->protect_highlights = protect_highlights;
		seq_args->amount = amount;
		seq_args->sigma = invsigma;
		seq_args->applyRotation = applyRotation;
		seq_args->seq = &com.seq;
		seq_args->fit = NULL;

		gtk_toggle_button_set_active(seq, FALSE);
		apply_banding_to_sequence(seq_args);
	} else {
		// Single image processing - use generic_image_worker
		struct banding_data *params = new_banding_data();
		if (!params) {
			PRINT_ALLOC_ERR;
			set_cursor_waiting(FALSE);
			return;
		}

		params->protect_highlights = protect_highlights;
		params->amount = amount;
		params->sigma = invsigma;
		params->applyRotation = applyRotation;
		params->seqEntry = NULL;
		params->seq = NULL;
		params->fit = NULL;

		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			PRINT_ALLOC_ERR;
			free_banding_data(params);
			free(params);
			set_cursor_waiting(FALSE);
			return;
		}

		args->fit = gfit;
		args->mem_ratio = 2.0f; // Banding needs ~2x memory (rotation + fix image)
		args->image_hook = banding_single_image_hook;
		args->idle_function = banding_single_idle;
		args->description = _("Canon Banding Reduction");
		args->verbose = TRUE;
		args->user = params;
		args->log_hook = banding_log_hook;
		args->max_threads = com.max_thread;
		args->for_preview = FALSE;
		args->for_roi = FALSE;

		if (!start_in_new_thread(generic_image_worker, args)) {
			free_banding_data(params);
			free(params);
			free(args);
			set_cursor_waiting(FALSE);
		}
	}
}

void on_checkbutton_fixbanding_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	static GtkWidget *scalebandingHighlightBox = NULL;
	static GtkWidget *spinbandingHighlightBox = NULL;
	gboolean is_active;

	if (scalebandingHighlightBox == NULL) {
		scalebandingHighlightBox = lookup_widget("scale_fixbanding_invsigma");
		spinbandingHighlightBox = lookup_widget("spin_fixbanding_invsigma");
	}

	is_active = gtk_toggle_button_get_active(togglebutton);
	gtk_widget_set_sensitive(scalebandingHighlightBox, is_active);
	gtk_widget_set_sensitive(spinbandingHighlightBox, is_active);
}
