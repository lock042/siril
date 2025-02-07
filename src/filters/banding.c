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
 *      B A N D I N G      R E D U C T I O N      M A N A G E M E N T        *
 ****************************************************************************/

int banding_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	struct banding_data *banding_args = (struct banding_data *)args->user;
	return BandingEngine(fit, banding_args->sigma, banding_args->amount,
			banding_args->protect_highlights, banding_args->applyRotation, SINGLE_THREADED);
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
	args->new_seq_prefix = banding_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->user = banding_args;

	banding_args->fit = NULL;	// not used here

	if (start_in_new_thread(generic_sequence_worker, args)) {
		free(banding_args->seqEntry);
		free(banding_args);
		free_generic_seq_args(args);
	}
}

// idle function executed at the end of the BandingEngine processing
gboolean end_BandingEngine(gpointer p) {
	struct banding_data *args = (struct banding_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	notify_gfit_modified();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);

	free(args);
	return FALSE;
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

/*** Reduces Banding in Canon DSLR images.
 * This code come from CanonBandingReduction.js v0.9.1, a script of
 * PixInsight, originally written by Georg Viehoever and
 * distributed under the terms of the GNU General Public License ******/
gpointer BandingEngineThreaded(gpointer p) {
	struct banding_data *args = (struct banding_data *) p;
	struct timeval t_start, t_end;

	siril_log_color_message(_("Banding Reducing: processing...\n"), "green");
	gettimeofday(&t_start, NULL);

	int retval = BandingEngine(args->fit, args->sigma, args->amount, args->protect_highlights, args->applyRotation, MULTI_THREADED);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	siril_add_idle(end_BandingEngine, args);

	return GINT_TO_POINTER(retval);
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

	struct banding_data *args = malloc(sizeof(struct banding_data));

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

	if (!protect_highlights)
		undo_save_state(&gfit, _("Canon Banding Reduction (amount=%.2lf)"), amount);
	else
		undo_save_state(&gfit, _("Canon Banding Reduction (amount=%.2lf, Protect=TRUE, invsigma=%.2lf)"),
				amount, invsigma);

	args->fit = &gfit;
	args->protect_highlights = protect_highlights;
	args->amount = amount;
	args->sigma = invsigma;
	args->applyRotation = gtk_toggle_button_get_active(vertical);
	args->seqEntry = strdup(gtk_entry_get_text(bandingSeqEntry));
	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = strdup("unband_");
		gtk_toggle_button_set_active(seq, FALSE);
		args->seq = &com.seq;
		apply_banding_to_sequence(args);
	} else {
		if (!start_in_new_thread(BandingEngineThreaded, args))
			free(args);
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
