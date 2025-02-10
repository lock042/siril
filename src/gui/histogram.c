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

#include <gsl/gsl_histogram.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/arithm.h"
#include "algos/statistics.h"
#include "algos/colors.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/utils.h"	// for lookup_widget()
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/siril_preview.h"
#include "gui/registration_preview.h"
#include "core/undo.h"
#include "histogram.h"

#define GRADIENT_HEIGHT 12

#undef HISTO_DEBUG

/* The gsl_histogram, documented here:
 * http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_21.html
 * is able to group values into bins, it does not need to handle all values
 * separately. That is useful for display purposes, but not currently used.
 */

// Type of invocation - HISTO_STRETCH or GHT_STRETCH
static int invocation = NO_STRETCH_SET_YET;

static gboolean closing = FALSE;
static gboolean sequence_working = FALSE;
// Parameters for use in calculations
static float _B = 0.5f, _D = 0.0f, _BP = 0.0f, _LP = 0.0f, _SP = 0.0f, _HP = 1.0f;
static clip_mode_t _clip_mode = RGBBLEND;
static gboolean do_channel[3];
static int _stretchtype = STRETCH_PAYNE_NORMAL;
static int _payne_colourstretchmodel = COL_INDEP;
static ght_compute_params compute_params;

static fits* fit = &gfit;

static gboolean auto_display_compensation = FALSE;

static double invxpos = -1.0;
// colors of layers histograms		R	G	B	RGB
static double histo_color_r[] = { 1.0, 0.0, 0.0, 0.0 };
static double histo_color_g[] = { 0.0, 1.0, 0.0, 0.0 };
static double histo_color_b[] = { 0.0, 0.0, 1.0, 0.0 };
// static float graph_height = 0.f;	// the max value of all bins
static guint64 clipped[] = { 0, 0 };

// The 4th toggle pointer remains NULL, but it makes the logic easier in redraw_histo()
static GtkToggleToolButton *toggles[MAXVPORT] = { NULL };
static GtkToggleToolButton *toggleGrid = NULL, *toggleCurve = NULL, *toggleOrig = NULL;

/* the original histogram, used as starting point of each computation */
static gsl_histogram *hist_backup[MAXVPORT] = { NULL, NULL, NULL, NULL };
static gsl_histogram *hist_sat_backup = NULL;

static float _midtones, _shadows, _highlights;

static gboolean _click_on_histo = FALSE;
static ScaleType _type_of_scale;
static gboolean lp_warning_given, hp_warning_given = FALSE;

void *huebuf = NULL, *satbuf_orig = NULL, *satbuf_working = NULL, *lumbuf = NULL;

static void set_histogram(gsl_histogram *histo, int layer);

static void setup_hsl();
static void clear_hsl();
static void set_sat_histogram(gsl_histogram *histo);

static int get_width_of_histo() {
	return gtk_widget_get_allocated_width(lookup_widget("drawingarea_histograms"));
}

static int get_height_of_histo() {
	return gtk_widget_get_allocated_height(lookup_widget("drawingarea_histograms"));
}

static void clear_hist_backup() {
	if (hist_backup[0]) {
		for (int i = 0; i < fit->naxes[2]; i++) {
			gsl_histogram_free(hist_backup[i]);
			hist_backup[i] = NULL;
		}
	}
	if (hist_sat_backup) {
		gsl_histogram_free(hist_sat_backup);
		hist_sat_backup = NULL;
	}
}

static void init_toggles() {
	if (!toggles[0]) {
		toggles[0] = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolRed"));
		toggles[1] = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolGreen"));
		toggles[2] = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolBlue"));
		do_channel[0] = gtk_toggle_tool_button_get_active(toggles[0]);
		do_channel[1] = gtk_toggle_tool_button_get_active(toggles[1]);
		do_channel[2] = gtk_toggle_tool_button_get_active(toggles[2]);
		toggleGrid = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolGrid"));
		toggleCurve = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolCurve"));
		toggleOrig = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolOrig"));
	}
}

static void histo_startup() {
	add_roi_callback(histo_change_between_roi_and_image);
	roi_supported(TRUE);
	copy_gfit_to_backup();
	if (fit->naxes[2] == 3 && _payne_colourstretchmodel == COL_SAT)
		setup_hsl();

	init_toggles();
	do_channel[0] = gtk_toggle_tool_button_get_active(toggles[0]);
	do_channel[1] = gtk_toggle_tool_button_get_active(toggles[1]);
	do_channel[2] = gtk_toggle_tool_button_get_active(toggles[2]);
	if (fit->naxes[2] == 3 && !(do_channel[0] && do_channel[1] && do_channel[2])) {
		if (_payne_colourstretchmodel == COL_HUMANLUM && invocation == GHT_STRETCH) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Human luminance colour model cannot be used: setting even weighted luminance colour model."));
			_payne_colourstretchmodel = COL_EVENLUM;
			gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_payne_colour_stretch_model")), COL_EVENLUM);
		}
		if (_payne_colourstretchmodel == COL_SAT && invocation == GHT_STRETCH) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Saturation stretch cannot be used: setting independent channels colour model."));
			_payne_colourstretchmodel = COL_INDEP;
			gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_payne_colour_stretch_model")), COL_INDEP);
		}
	}
	// also get the backup histogram
	compute_histo_for_gfit();
	for (int i = 0; i < fit->naxes[2]; i++)
		hist_backup[i] = gsl_histogram_clone(com.layers_hist[i]);
	if (com.sat_hist) {
		if (hist_sat_backup)
			gsl_histogram_free(hist_sat_backup);
		hist_sat_backup = gsl_histogram_clone(com.sat_hist);
	}
}

static void histo_close(gboolean revert, gboolean update_image_if_needed) {
	if (revert) {

		for (int i = 0; i < fit->naxes[2]; i++) {
			set_histogram(hist_backup[i], i);
			hist_backup[i] = NULL;
		}
		set_sat_histogram(hist_sat_backup);
		hist_sat_backup = NULL;
		if (is_preview_active() && !copy_backup_to_gfit() && update_image_if_needed) {
			set_cursor_waiting(TRUE);
			notify_gfit_modified();
		}
	}
	// free data
	if(!sequence_working && !revert) {
		clear_hsl();
	}
	clear_backup();
	clear_hist_backup();
	roi_supported(FALSE);
	remove_roi_callback(histo_change_between_roi_and_image);
}

static void hsl_to_fit (void* h, void* s, void* l) {
	size_t npixels = fit->rx * fit->ry;
	if (fit->type == DATA_FLOAT) {
		float *hf = (float*) h;
		float* sf = (float*) s;
		float* lf = (float*) l;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0 ; i < npixels ; i++) {
			hsl_to_rgbf(hf[i], sf[i], lf[i], &fit->fpdata[0][i], &fit->fpdata[1][i], &fit->fpdata[2][i]);
		}
	} else {
		WORD *hw = (WORD*) h;
		WORD* sw = (WORD*) s;
		WORD* lw = (WORD*) l;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0 ; i < npixels ; i++) {
			hslw_to_rgbw(hw[i], sw[i], lw[i], &fit->pdata[0][i], &fit->pdata[1][i], &fit->pdata[2][i]);
		}
	}
}

static void fit_to_hsl() {
	size_t npixels = fit->rx * fit->ry;
	if (fit->type == DATA_FLOAT) {
		float* hf = (float*) huebuf;
		float* sf = (float*) satbuf_orig;
		float* lf = (float*) lumbuf;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0 ; i < npixels ; i++) {
			rgb_to_hslf(fit->fpdata[0][i], fit->fpdata[1][i], fit->fpdata[2][i], &hf[i], &sf[i], &lf[i]);
		}
		memcpy(satbuf_working, satbuf_orig, npixels * sizeof(float));
	} else {
		WORD *hw = (WORD*) huebuf;
		WORD* sw = (WORD*) satbuf_orig;
		WORD* lw = (WORD*) lumbuf;
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0 ; i < npixels ; i++) {
			rgbw_to_hslw(fit->pdata[0][i], fit->pdata[1][i], fit->pdata[2][i], &hw[i], &sw[i], &lw[i]);
		}
		memcpy(satbuf_working, satbuf_orig, npixels * sizeof(WORD));
	}
}

static void histo_recompute() {
	set_cursor("progress");
	gboolean auto_comp = auto_display_compensation;
	auto_display_compensation = FALSE;
	copy_backup_to_gfit();

	if (invocation == HISTO_STRETCH) {
		struct mtf_params params = { .shadows = _shadows, .midtones = _midtones, .highlights = _highlights, .do_red = do_channel[0], .do_green = do_channel[1], .do_blue = do_channel[2] };
		apply_linked_mtf_to_fits(fit, fit, params, TRUE);
	// com.layers_hist should be good, update_histo_mtf() is always called before
	} else if (invocation == GHT_STRETCH) {
		struct ght_params params_ght = { .B = _B, .D = _D, .LP = (float) _LP, .SP = (float) _SP, .HP = (float) _HP, .BP = _BP, .stretchtype = _stretchtype, .payne_colourstretchmodel = _payne_colourstretchmodel, .do_red = do_channel[0], .do_green = do_channel[1], .do_blue = do_channel[2], .clip_mode = _clip_mode };
		if (_payne_colourstretchmodel == COL_SAT) {
			if (fit->type == DATA_FLOAT)
				apply_linked_ght_to_fbuf_indep((float*) satbuf_orig, (float*) satbuf_working, fit->rx * fit->ry, 1, &params_ght, TRUE);
			else
				apply_linked_ght_to_Wbuf_indep((WORD*) satbuf_orig, (WORD*) satbuf_working, fit->rx * fit->ry, 1, &params_ght, TRUE);
			hsl_to_fit(huebuf, satbuf_working, lumbuf);
		} else {
			apply_linked_ght_to_fits(fit, fit, &params_ght, TRUE);
		}
	}
	if (auto_comp) {
		int depth = fit->naxes[2];
		if (depth == 1) {
			fits_change_depth(fit, 3);
			if (fit->type == DATA_FLOAT) {
				memcpy(fit->fpdata[1], fit->fdata, fit->rx * fit->ry * sizeof(float));
				memcpy(fit->fpdata[2], fit->fdata, fit->rx * fit->ry * sizeof(float));
			} else {
				memcpy(fit->pdata[1], fit->data, fit->rx * fit->ry * sizeof(WORD));
				memcpy(fit->pdata[2], fit->data, fit->rx * fit->ry * sizeof(WORD));
			}
		}
		// WARNING: the following section is *ONLY* applicable to autostretch. It should
		// not be copied into any other code as it will mess things up if the image
		// ICC profile is not the same as the monitor ICC profile.
		cmsHPROFILE temp = copyICCProfile(fit->icc_profile);
		cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = copyICCProfile(gui.icc.monitor);
		siril_colorspace_transform(fit, temp);
		cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = copyICCProfile(temp);
		cmsCloseProfile(temp);
		////////////////////////////////////////////////////////////////////////////////
		if (depth == 1) {
			size_t npixels = fit->rx * fit->ry;
			if (fit->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
				for (size_t i = 0 ; i < npixels; i++) {
					fit->fdata[i] = (fit->fpdata[0][i] + fit->fpdata[1][i] + fit->fpdata[2][i])/3.f;
				}
			} else {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
				for (size_t i = 0 ; i < npixels; i++) {
					fit->data[i] = (fit->pdata[0][i] + fit->pdata[1][i] + fit->pdata[2][i])/3.f;
				}
			}
			fits_change_depth(fit, 1);
		}
	}
	notify_gfit_modified();
}

static void _init_clipped_pixels() {
	clipped[0] = 0;
	clipped[1] = 0;
}

static void _initialize_clip_text() {
	static GtkEntry *clip_high = NULL, *clip_low = NULL;

	if (clip_high == NULL) {
		clip_high = GTK_ENTRY(lookup_widget("clip_highlights"));
		clip_low = GTK_ENTRY(lookup_widget("clip_shadows"));
	}
	gtk_entry_set_text(clip_low, "0.000%");
	gtk_entry_set_text(clip_high, "0.000%");
}

static void _update_entry_text() {
	GtkEntry *histoMidEntry = GTK_ENTRY(lookup_widget("histoMidEntry"));
	GtkEntry *histoShadEntry = GTK_ENTRY(lookup_widget("histoShadEntry"));
	GtkEntry *histoHighEntry = GTK_ENTRY(lookup_widget("histoHighEntry"));
	gchar *buffer;

	buffer = g_strdup_printf("%.7f", _shadows);
	gtk_entry_set_text(histoShadEntry, buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", _highlights);
	gtk_entry_set_text(histoHighEntry, buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", _midtones);
	gtk_entry_set_text(histoMidEntry, buffer);
	g_free(buffer);
}

static void _update_clipped_pixels(size_t data) {
	static GtkEntry *clip_high = NULL, *clip_low = NULL;
	double tmp;
	char buffer[16];

	if (clip_high == NULL) {
		clip_high = GTK_ENTRY(lookup_widget("clip_highlights"));
		clip_low = GTK_ENTRY(lookup_widget("clip_shadows"));
	}
	tmp = (double)clipped[1] * 100.0 / (double)data;
	g_snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_entry_set_text(clip_high, buffer);
	tmp = (double)clipped[0] * 100.0 / (double)data;
	g_snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_entry_set_text(clip_low, buffer);
}

static int is_histogram_visible() {
	GtkWidget *window = lookup_widget("histogram_dialog");
	return gtk_widget_get_visible(window);
}

// sets the channel names of the toggle buttons in the histogram window, based on
// the number of layers of gfit
static void set_histo_toggles_names() {
	init_toggles();

	if (fit->naxis == 2) {
		gtk_widget_set_tooltip_text(GTK_WIDGET(toggles[0]), _("Toggles whether to apply the stretch to the monochrome channel"));
		GtkWidget *w;
		if (com.pref.gui.combo_theme == 0) {
			w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/monochrome_dark.svg");
		} else {
			w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/monochrome.svg");
		}
		gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(toggles[0]), w);
		gtk_widget_show(w);
		gtk_widget_set_visible(GTK_WIDGET(toggles[1]), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[2]), FALSE);
		/* visible has no effect in GTK+ 3.12, trying sensitive too
		 * Yes it does. The solution is to call the window (widget)
		 * with gtk_widget_show and not gtk_widget_show_all */
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[1]), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[2]), FALSE);

	} else {
		gtk_widget_set_tooltip_text(GTK_WIDGET(toggles[0]), _("Red channel"));
		GtkWidget* w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/r.svg");
		gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(toggles[0]), w);
		gtk_widget_show(w);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[1]), TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[2]), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[1]), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[2]), TRUE);
	}
}

static double get_histoZoomValueH() {
	static GtkAdjustment *histoAdjZoomH = NULL;
	if (!histoAdjZoomH)
		histoAdjZoomH = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "histoAdjZoomH"));

	return gtk_adjustment_get_value(histoAdjZoomH);
}

static double get_histoZoomValueV() {
	static GtkAdjustment *histoAdjZoomV = NULL;
	if (!histoAdjZoomV)
		histoAdjZoomV = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "histoAdjZoomV"));

	return gtk_adjustment_get_value(histoAdjZoomV);
}

static void adjust_histogram_vport_size() {
	GtkWidget *drawarea, *vport;
	int targetW = 339, targetH = 281, cur_width = 339, cur_height = 281;
	double zoomH = 1.0;
	zoomH = get_histoZoomValueH();
	double zoomV = 1.0;
	zoomV = get_histoZoomValueV();

	drawarea = lookup_widget("drawingarea_histograms");
	vport = lookup_widget("viewport1");

	cur_width = gtk_widget_get_allocated_width(vport);
	cur_height = gtk_widget_get_allocated_height(vport);
	targetW = (int) (((double)cur_width) * zoomH);
	targetH = (int) (((double)cur_height) * zoomV);
	gtk_widget_set_size_request(drawarea, targetW, targetH);
#ifdef HISTO_DEBUG
	fprintf(stdout, "Histo vport size (%d, %d)\n", targetW, targetH);
#endif
}

size_t get_histo_size(fits *fit) {
	if (fit->type == DATA_USHORT) {
		if (fit->orig_bitpix == BYTE_IMG)
			return UCHAR_MAX;
	}
	return (size_t)USHRT_MAX;
}

// create a new histogram object for the passed fit and layer
gsl_histogram* computeHisto(fits *fit, int layer) {
	g_assert(layer < 3);
	size_t i, ndata, size;

	size = get_histo_size(fit);
	gsl_histogram *histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, fit->type == DATA_FLOAT ? 1.0 + 1.0 / size : size + 1);
	ndata = fit->naxes[0] * fit->naxes[1];

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
	{
		gsl_histogram *histo_thr = gsl_histogram_alloc(size + 1);
		gsl_histogram_set_ranges_uniform(histo_thr, 0, fit->type == DATA_FLOAT ? 1.0 + 1.0 / size : size + 1);

		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[layer];
#ifdef _OPENMP
#pragma omp for private(i) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				if (buf[i] == 0)
					continue;
				gsl_histogram_increment(histo_thr, (double) buf[i]);
			}
		} else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[layer];
#ifdef _OPENMP
#pragma omp for private(i) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				if (buf[i] == 0.f)
					continue;
				gsl_histogram_increment(histo_thr, (double) buf[i]);
			}
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			gsl_histogram_add(histo, histo_thr);
		}
		gsl_histogram_free(histo_thr);
	}

	return histo;
}

// create a new histogram object for the passed float buffer (used for sat)
gsl_histogram* computeHistoSat(void* buf) {
	size_t i, ndata, size;

	size = get_histo_size(&gfit);
	gsl_histogram *histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, 1.0 + 1.0 / size);
	ndata = fit->naxes[0] * fit->naxes[1];

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
	{
		gsl_histogram *histo_thr = gsl_histogram_alloc(size + 1);
		gsl_histogram_set_ranges_uniform(histo_thr, 0, 1.0 + 1.0 / size);
		if (fit->type == DATA_FLOAT) {
			float *fbuf = (float*) buf;
#ifdef _OPENMP
#pragma omp for private(i) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				gsl_histogram_increment(histo_thr, (double) fbuf[i]);
			}
		} else {
			double invnorm = 1.0 / USHRT_MAX_DOUBLE;
			WORD * Wbuf = (WORD*) buf;
#ifdef _OPENMP
#pragma omp for private(i) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				gsl_histogram_increment(histo_thr, Wbuf[i] * invnorm);
			}
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			gsl_histogram_add(histo, histo_thr);
		}
		gsl_histogram_free(histo_thr);
	}
	return histo;
}

static gboolean is_log_scale() {
	static GtkToggleButton *HistoCheckLogButton = NULL;

	if (HistoCheckLogButton == NULL)
		HistoCheckLogButton = GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckLogButton"));
	return (gtk_toggle_button_get_active(HistoCheckLogButton));
}

static void draw_curve(cairo_t *cr, int width, int height) {
	// draw curve
	int k;
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_source_rgb(cr, 0.98, 0.5, 0.45);

	// Draw x=y
	cairo_set_line_width(cr,0.25);
	cairo_move_to(cr,0,height);
	cairo_line_to(cr, width, 0);
	cairo_stroke(cr);

	cairo_set_line_width(cr, 1.0);

	if (invocation == HISTO_STRETCH) {
		for (k = 0; k < width + 1; k++) {
			float x = k / (float) width;
			float y = MTF(x, _midtones, _shadows, _highlights);
			cairo_line_to(cr, k, height * (1 - y));
		}
	} else if (invocation == GHT_STRETCH) {
		GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, _stretchtype);
		for (k = 0; k < width + 1; k++) {
			float x = k / (float) width;
			float y = (float) GHT((float) x, _B, _D, _LP, _SP, _HP, _BP, _stretchtype, &compute_params);
			cairo_line_to(cr, k, height * (1 - y));
		}
	}
	cairo_stroke(cr);
}

void draw_grid(cairo_t *cr, int width, int height) {
	double dash_format[] = { 1.0, 1.0 };

	cairo_set_line_width(cr, 1.0);
	cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
	// quarters in solid, eights in dashed line
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_move_to(cr, width * 0.25, 0);
	cairo_line_to(cr, width * 0.25, height);
	cairo_move_to(cr, width * 0.5, 0);
	cairo_line_to(cr, width * 0.5, height);
	cairo_move_to(cr, width * 0.75, 0);
	cairo_line_to(cr, width * 0.75, height);

	cairo_set_line_width(cr, 1.0);
	cairo_move_to(cr, 0, height * 0.25);
	cairo_line_to(cr, width, height * 0.25);
	cairo_move_to(cr, 0, height * 0.5);
	cairo_line_to(cr, width, height * 0.5);
	cairo_move_to(cr, 0, height * 0.75);
	cairo_line_to(cr, width, height * 0.75);

	cairo_stroke(cr);

	cairo_set_line_width(cr, 1.0);
	cairo_set_dash(cr, dash_format, 2, 0);
	cairo_move_to(cr, width * 0.125, 0);
	cairo_line_to(cr, width * 0.125, height);
	cairo_move_to(cr, width * 0.375, 0);
	cairo_line_to(cr, width * 0.375, height);
	cairo_move_to(cr, width * 0.625, 0);
	cairo_line_to(cr, width * 0.625, height);
	cairo_move_to(cr, width * 0.875, 0);
	cairo_line_to(cr, width * 0.875, height);

	cairo_set_line_width(cr, 1.0);
	cairo_move_to(cr, 0, height * 0.125);
	cairo_line_to(cr, width, height * 0.125);
	cairo_move_to(cr, 0, height * 0.375);
	cairo_line_to(cr, width, height * 0.375);
	cairo_move_to(cr, 0, height * 0.625);
	cairo_line_to(cr, width, height * 0.625);
	cairo_move_to(cr, 0, height * 0.875);
	cairo_line_to(cr, width, height * 0.875);
	cairo_stroke(cr);

}

// erase image and redraw the background color and grid
void erase_histo_display(cairo_t *cr, int width, int height) {
	gboolean drawGrid, drawCurve;
	// clear all with background color
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);

	// draw grid
	init_toggles();
	drawGrid = gtk_toggle_tool_button_get_active(toggleGrid);
	drawCurve = gtk_toggle_tool_button_get_active(toggleCurve);

	if (drawGrid)
		draw_grid(cr, width, height);
	if (drawCurve)
		draw_curve(cr, width, height);
}

static void draw_gradient(cairo_t *cr, int width, int height) {
	cairo_pattern_t *pat;

	pat = cairo_pattern_create_linear(0.0, 0.0, width, 0.0);
	cairo_pattern_add_color_stop_rgb(pat, 0, 0, 0, 0);
	cairo_pattern_add_color_stop_rgb(pat, width, 1, 1, 1);
	cairo_rectangle(cr, 0, height - GRADIENT_HEIGHT, width, GRADIENT_HEIGHT);
	cairo_set_source(cr, pat);
	cairo_fill(cr);

	cairo_pattern_destroy(pat);
}

static void draw_slider(cairo_t *cr, int width, int height, int xpos) {
	if (xpos > width / 2) {
		cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	} else {
		cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
	}
	cairo_move_to(cr, -10 + xpos, height);
	cairo_line_to(cr, 10 + xpos, height);
	cairo_line_to(cr, xpos, height - GRADIENT_HEIGHT);
	cairo_line_to(cr, -10 + xpos, height);
	cairo_stroke(cr);
}

void display_scale(cairo_t *cr, int width, int height) {
	draw_gradient(cr, width, height);
	if (invocation == HISTO_STRETCH) {
		float delta = ((_highlights - _shadows) * _midtones) + _shadows;
		draw_slider(cr, width, height, _shadows * width);
		draw_slider(cr, width, height, (delta) * width);
		draw_slider(cr, width, height, _highlights * width);
	}
}

void display_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width,
		int height, double zoomH, double zoomV, gboolean isOrig, gboolean is_log) {
	if (!histo) return;
	if (width <= 0) return;
	int current_bin;
	size_t norm = gsl_histogram_bins(histo) - 1;

	float vals_per_px = (float)norm / (float)width;	// size of a bin
	size_t i, nb_orig_bins = gsl_histogram_bins(histo);

	// We need to store the binned histogram in order to find the binned maximum
	static gfloat *displayed_values = NULL;
	static int nb_bins_allocated = 0;
	/* we create a bin for each pixel in the displayed width.
	 * nb_bins_allocated is thus equal to the width of the image */
	if (nb_bins_allocated != width) {
		gfloat *tmp;
		nb_bins_allocated = width;
		tmp = realloc(displayed_values, nb_bins_allocated * sizeof(gfloat));
		if (!tmp) {
			if (displayed_values != NULL) {
				g_free(displayed_values);
				displayed_values = NULL;
			}
			PRINT_ALLOC_ERR;
			histo_close(TRUE, TRUE);
			return;
		}
		displayed_values = tmp;
		memset(displayed_values, 0, nb_bins_allocated);
	}
	if (fit->naxis == 2)
		cairo_set_source_rgb(cr, 255.0, 255.0, 255.0);
	else if (layer < 0)
		cairo_set_source_rgb(cr, 255.0, 255.0, 0.0);
	else
		cairo_set_source_rgb(cr, histo_color_r[layer], histo_color_g[layer],
				histo_color_b[layer]);
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1.5);
	if (isOrig || layer == -2)
		cairo_set_line_width(cr, 0.5);

	// first loop builds the bins and finds the maximum
	i = 0;
	double graph_height = 0.f;
	current_bin = 0;
	do {
		double bin_val = 0.f;
		while (i < nb_orig_bins
				&& (double)i / vals_per_px <= current_bin + 0.5f) {
			bin_val += gsl_histogram_get(histo, i);
			i++;
		}
		if (is_log && bin_val != 0.f) {
			bin_val = logf(bin_val);
		}
		displayed_values[current_bin] = bin_val;
		if (bin_val > graph_height)	// check for maximum
			graph_height = bin_val;
		current_bin++;
	} while (i < nb_orig_bins && current_bin < nb_bins_allocated);
	for (i = 0; i < nb_bins_allocated; i++) {
		double bin_height = height - height * displayed_values[i] / graph_height;
		cairo_line_to(cr, i, bin_height);
	}
	cairo_stroke(cr);
}

static void apply_mtf_to_histo(gsl_histogram *histo, float norm,
		float m, float lo, float hi) {

	size_t int_norm = (size_t)norm;
	gsl_histogram *mtf_histo = gsl_histogram_alloc(int_norm + 1);

	gsl_histogram_set_ranges_uniform(mtf_histo, 0, norm);
// disabled because of ISSUE #136 (https://free-astro.org/bugs/view.php?id=136)
// Update by A Knagg-Baugh - the issue only relates to MacOS so I've limited the
// disabling to that OS with a #ifndef. Similarly for the other instance below
#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
#endif
	for (size_t i = 0; i < int_norm + 1; i++) {
		WORD mtf;
		float binval = gsl_histogram_get(histo, i);
		float pxl = ((float)i / norm);
		guint64 clip[2] = { 0, 0 };

		if (i < roundf_to_WORD(lo * norm)) {
			pxl = lo;
			clip[0] += binval;
		} else if (i > roundf_to_WORD(hi * norm)) {
			pxl = hi;
			clip[1] += binval;
		}
		mtf = roundf_to_WORD(MTF(pxl, m, lo, hi) * norm);
		gsl_histogram_accumulate(mtf_histo, (double)mtf, (double)binval);
#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp critical
#endif
#endif
		{
			clipped[0] += clip[0];
			clipped[1] += clip[1];
		}
	}
	gsl_histogram_memcpy(histo, mtf_histo);
	gsl_histogram_free(mtf_histo);
}

static void apply_ght_to_histo(gsl_histogram *histo, float norm,
		float m, float lo, float hi) {
	size_t int_norm = (size_t)norm;
	gsl_histogram *mtf_histo = gsl_histogram_alloc(int_norm + 1);

	gsl_histogram_set_ranges_uniform(mtf_histo, 0, norm);

	GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, _stretchtype);

#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
#endif
	for (size_t i = 0; i < int_norm + 1; i++) {
		float ght;
		float binval = gsl_histogram_get(histo, i);
		float pxl = ((float)i / norm);
		guint64 clip[2] = { 0, 0 };

		if (i < roundf_to_WORD(lo * norm)) {
			pxl = lo;
			clip[0] += binval;
		} else if (i > roundf_to_WORD(hi * norm)) {
			pxl = hi;
			clip[1] += binval;
		}
		ght = roundf_to_WORD(GHT(pxl, _B, _D, _LP, _SP, _HP, _BP, _stretchtype, &compute_params) * norm);
		gsl_histogram_accumulate(mtf_histo, (double) ght, (double) binval);
#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp critical
#endif
#endif
		{
			clipped[0] += clip[0];
			clipped[1] += clip[1];
		}
	}
	gsl_histogram_memcpy(histo, mtf_histo);
	gsl_histogram_free(mtf_histo);
}

void on_histoZoom100_clicked(GtkButton *button, gpointer user_data) {
	GtkAdjustment *histoAdjZoomH, *histoAdjZoomV;

	histoAdjZoomH = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "histoAdjZoomH"));
	histoAdjZoomV = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "histoAdjZoomV"));

	gtk_adjustment_set_value(histoAdjZoomH, 1.0);
	gtk_adjustment_set_value(histoAdjZoomV, 1.0);
}

static void reset_cursors_and_values(gboolean full_reset) {
	gtk_toggle_tool_button_set_active(toggleGrid, TRUE);
	gtk_toggle_tool_button_set_active(toggleCurve, TRUE);
	gtk_toggle_tool_button_set_active(toggleOrig, BOOL_TRUE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkMTFSeq")), FALSE);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entryMTFSeq")), "stretch_");
//	on_histoZoom100_clicked(NULL, NULL);
	if (full_reset) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckLogButton")),
									(com.pref.gui.display_histogram_mode == LOG_DISPLAY ? TRUE : FALSE));
		for (int i = 0; i < 3; ++i)
			gtk_toggle_tool_button_set_active(toggles[i], TRUE);
	}

	if (invocation == HISTO_STRETCH) {
		_shadows = 0.f;
		_midtones = 0.5f;
		_highlights = 1.0f;
	} else if (invocation == GHT_STRETCH) {
		_D = 0.0f;
		_B = 0.0f;
		_SP = 0.0f;
		_LP = 0.0f;
		_HP = 1.0f;
		_BP = 0.0f;
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtD")), _D);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtB")), _B);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtSP")), _SP);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtLP")), _LP);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtHP")), _HP);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtBP")), _BP);
		if (full_reset) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_payne_colour_stretch_model")), 0);
			gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_payneTyp")), 0);
			gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("histo_clip_mode")), 2);
			gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")), TRUE);
		}
	}
	_init_clipped_pixels();
	_initialize_clip_text();
	_update_entry_text();
	update_gfit_histogram_if_needed();
}

static void queue_window_redraw() {
	static GtkWidget *drawarea = NULL;
	if (!drawarea)
		drawarea = lookup_widget("drawingarea_histograms");
	gtk_widget_queue_draw(drawarea);
}

static void update_histo_mtf() {
	if (!com.layers_hist[0]) return;
	if (!hist_backup[0]) return;
	float norm = (float)gsl_histogram_bins(com.layers_hist[0]) - 1;

	_init_clipped_pixels();
	for (long i = 0; i < fit->naxes[2]; i++) {
		gsl_histogram_memcpy(com.layers_hist[i], hist_backup[i]);
		if (do_channel[i]) {
			if (invocation == HISTO_STRETCH) {
				apply_mtf_to_histo(com.layers_hist[i], norm, _midtones, _shadows, _highlights);
			} else if (invocation == GHT_STRETCH) {
				if (_payne_colourstretchmodel != COL_SAT)
					apply_ght_to_histo(com.layers_hist[i], norm, _SP, _BP, 1.0f);
			}
		}
	}
	if (invocation == GHT_STRETCH && _payne_colourstretchmodel == COL_SAT) {
		gsl_histogram_memcpy(com.sat_hist, hist_sat_backup);
		apply_ght_to_histo(com.sat_hist, norm, _SP, _BP, 1.0f);
	}
	size_t data = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	_update_clipped_pixels(data);
	queue_window_redraw();
}

static int histo_update_preview() {
	fit = gui.roi.active ? &gui.roi.fit : &gfit;
	if (!closing)
		histo_recompute();
	return 0;
}

static void set_histogram(gsl_histogram *histo, int layer) {
	g_assert(layer >= 0 && layer < MAXVPORT);
	if (com.layers_hist[layer])
		gsl_histogram_free(com.layers_hist[layer]);
	com.layers_hist[layer] = histo;
}

static void set_sat_histogram(gsl_histogram *histo) {
	if (com.sat_hist)
		gsl_histogram_free(com.sat_hist);
	com.sat_hist = histo;
}

static gboolean on_gradient(GdkEvent *event, int width, int height) {
	return ((GdkEventButton*) event)->x > 0
			&& ((GdkEventButton*) event)->x < width
			&& ((GdkEventButton*) event)->y > height - GRADIENT_HEIGHT
			&& ((GdkEventButton*) event)->y < height;
}

/*
 * Public functions
 */

gsl_histogram* computeHisto_Selection(fits* fit, int layer,
		rectangle *selection) {
	g_assert(layer < 3);

	size_t size = get_histo_size(fit);
	gsl_histogram* histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, fit->type == DATA_FLOAT ? 1.0 : size);
	size_t stridefrom = fit->rx - selection->w;

	if (fit->type == DATA_USHORT) {
		WORD *from = fit->pdata[layer] + (fit->ry - selection->y - selection->h) * fit->rx
			+ selection->x;
		for (size_t i = 0; i < selection->h; i++) {
			for (size_t j = 0; j < selection->w; j++) {
				gsl_histogram_increment(histo, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	else if (fit->type == DATA_FLOAT) {
		float *from = fit->fpdata[layer] + (fit->ry - selection->y - selection->h) * fit->rx
			+ selection->x;
		for (size_t i = 0; i < selection->h; i++) {
			for (size_t j = 0; j < selection->w; j++) {
				gsl_histogram_increment(histo, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	return histo;
}

/* call from main thread */
void compute_histo_for_gfit() {
	int nb_layers = 3;
	if (fit->naxis == 2)
		nb_layers = 1;
	for (int i = 0; i < nb_layers; i++) {
		if (!com.layers_hist[i])
			set_histogram(computeHisto(&gfit, i), i);
	}
	if (nb_layers == 3 && satbuf_working)
		set_sat_histogram(computeHistoSat(satbuf_working));
	set_histo_toggles_names();
}

/* call from any thread */
void invalidate_gfit_histogram() {
	for (int layer = 0; layer < MAXVPORT; layer++) {
		set_histogram(NULL, layer);
	}
	set_sat_histogram(NULL);
}

/* call from main thread */
void update_gfit_histogram_if_needed() {
	invalidate_gfit_histogram();
	if (is_histogram_visible()) {
		compute_histo_for_gfit();
		queue_window_redraw();
	}
}

static int mtf_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct mtf_data *m_args = (struct mtf_data*) args->user;
	apply_linked_mtf_to_fits(fit, fit, m_args->params, FALSE);
	return 0;
}

static int ght_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_, int threads) {
	struct ght_data *m_args = (struct ght_data*) args->user;
	if (m_args->params_ght->payne_colourstretchmodel == COL_SAT)
		apply_sat_ght_to_fits(fit, m_args->params_ght, FALSE);
	else
		apply_linked_ght_to_fits(fit, fit, m_args->params_ght, FALSE);
	return 0;
}

static void setup_hsl() {
	if (huebuf)
		free(huebuf);
	if (fit->type == DATA_FLOAT) {
		huebuf = malloc(fit->rx * fit->ry * sizeof(float));
		if (satbuf_orig)
			free(satbuf_orig);
		satbuf_orig = malloc(fit->rx * fit->ry * sizeof(float));
		if (satbuf_working)
			free(satbuf_working);
		satbuf_working = malloc(fit->rx * fit->ry * sizeof(float));
		if (lumbuf)
			free(lumbuf);
		lumbuf = malloc(fit->rx * fit->ry * sizeof(float));
	} else {
		huebuf = malloc(fit->rx * fit->ry * sizeof(WORD));
		if (satbuf_orig)
			free(satbuf_orig);
		satbuf_orig = malloc(fit->rx * fit->ry * sizeof(WORD));
		if (satbuf_working)
			free(satbuf_working);
		satbuf_working = malloc(fit->rx * fit->ry * sizeof(WORD));
		if (lumbuf)
			free(lumbuf);
		lumbuf = malloc(fit->rx * fit->ry * sizeof(WORD));
	}
	fit_to_hsl();
	set_sat_histogram(computeHistoSat(satbuf_working));
	if (hist_sat_backup)
		gsl_histogram_free(hist_sat_backup);
	hist_sat_backup = gsl_histogram_clone(com.sat_hist);
}

static void clear_hsl() {
	free(huebuf);
	huebuf = NULL;
	free(satbuf_orig);
	satbuf_orig = NULL;
	free(satbuf_working);
	satbuf_working = NULL;
	free(lumbuf);
	lumbuf = NULL;
	gsl_histogram_free(com.sat_hist);
	com.sat_hist = NULL;
	gsl_histogram_free(hist_sat_backup);
	hist_sat_backup = NULL;
}

/* Callback functions */

void histo_change_between_roi_and_image() {
	// This should be called if the ROI is set, changed or cleared to ensure the
	// histogram dialog continues to process the right data. Especially important
	// in GHT saturation stretch mode.
	fit = gui.roi.active ? &gui.roi.fit : &gfit;
	if (_payne_colourstretchmodel == COL_SAT) {
		clear_hsl();
		setup_hsl();
	}
	update_histo_mtf();
	gui.roi.operation_supports_roi = TRUE;
	// If we are showing the preview, update it after the ROI change.
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

gboolean redraw_histo(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width = 339, height = 281;
	double zoomH = 1, zoomV = 1;

	init_toggles();
	width = get_width_of_histo();
	height = get_height_of_histo();
#ifdef HISTO_DEBUG
	fprintf(stdout, "histogram redraw\n");
	fprintf(stdout, "w = %d and h = %d\n", width, height);
#endif
	zoomH = get_histoZoomValueH();
	zoomV = get_histoZoomValueV();

	if (height == 1)
		return FALSE;
	erase_histo_display(cr, width, height - GRADIENT_HEIGHT);

	for (i = 0; i < MAXVPORT; i++) {
		if (com.layers_hist[i]) {
			if (gtk_toggle_tool_button_get_active(toggleOrig)) {
				display_histo(hist_backup[i], cr, i, width, height - GRADIENT_HEIGHT, zoomH, zoomV, TRUE, is_log_scale());
			}
			if (!toggles[i] || gtk_toggle_tool_button_get_active(toggles[i])) {
				display_histo(com.layers_hist[i], cr, i, width, height - GRADIENT_HEIGHT, zoomH, zoomV, FALSE, is_log_scale());
			}
		}
	}
	if (invocation == GHT_STRETCH && _payne_colourstretchmodel == COL_SAT && fit->naxes[2] == 3) {
		display_histo(com.sat_hist, cr, -1, width, height - GRADIENT_HEIGHT, zoomH, zoomV, FALSE, is_log_scale());
		if (gtk_toggle_tool_button_get_active(toggleOrig))
			display_histo(hist_sat_backup, cr, -2, width, height - GRADIENT_HEIGHT, zoomH, zoomV, FALSE, is_log_scale());
	}
	display_scale(cr, width, height);
	return FALSE;
}

void on_histo_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	do_channel[0] = gtk_toggle_tool_button_get_active(toggles[0]);
	do_channel[1] = gtk_toggle_tool_button_get_active(toggles[1]);
	do_channel[2] = gtk_toggle_tool_button_get_active(toggles[2]);
	if (fit->naxes[2] == 3 && !(do_channel[0] && do_channel[1] && do_channel[2])) {
		if (_payne_colourstretchmodel == COL_HUMANLUM && invocation == GHT_STRETCH) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Human luminance colour model cannot be used: setting even weighted luminance colour model."));
			_payne_colourstretchmodel = COL_EVENLUM;
			gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_payne_colour_stretch_model")), COL_EVENLUM);
		}
		if (_payne_colourstretchmodel == COL_SAT && invocation == GHT_STRETCH) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Saturation stretch cannot be used: setting independent channels colour model."));
			_payne_colourstretchmodel = COL_INDEP;
			gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_payne_colour_stretch_model")), COL_INDEP);
		}
	}

	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_histogram_window_show(GtkWidget *object, gpointer user_data) {
	closing = FALSE;
	histo_startup();
	_initialize_clip_text();
	reset_cursors_and_values(TRUE);
	compute_histo_for_gfit();
}

void on_button_histo_close_clicked(GtkButton *button, gpointer user_data) {
	closing = TRUE;
	set_cursor_waiting(TRUE);
	histo_close(TRUE, TRUE);
	set_cursor_waiting(FALSE);
	siril_close_dialog("histogram_dialog");
}

void on_button_histo_reset_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	reset_cursors_and_values(FALSE);
	histo_close(TRUE, TRUE);
	histo_startup();
	set_cursor_waiting(FALSE);
}

gboolean on_scale_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
	return FALSE;
}

void on_button_histo_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (invocation == HISTO_STRETCH) {
		if ((_midtones == 0.5f) && (_shadows == 0.0f) && (_highlights == 1.0f)) {
			return;
		}
	} else if (invocation == GHT_STRETCH) {
		if ((_D == 0.0f) && (_BP == 0.0f)) {
			return;
		}
	}

	GtkToggleButton *seq_button = GTK_TOGGLE_BUTTON(lookup_widget("checkMTFSeq"));
	if (gtk_toggle_button_get_active(seq_button)
			&& sequence_is_loaded()) {
		/* Apply to the whole sequence */
		sequence_working = TRUE;
		if (invocation == HISTO_STRETCH) {
			struct mtf_data *args = malloc(sizeof(struct mtf_data));
			struct mtf_params params = { .shadows = _shadows, .midtones = _midtones, .highlights = _highlights, .do_red = do_channel[0], .do_green = do_channel[1], .do_blue = do_channel[2] };
			fprintf(stdout, "%d %d %d\n", params.do_red, params.do_green, params.do_blue);
			args->params = params;
			args->seqEntry = strdup( gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryMTFSeq"))));
			if (args->seqEntry && args->seqEntry[0] == '\0')
				args->seqEntry = strdup("stretch_");
			args->seq = &com.seq;
		/* here it is a bit tricky.
		 * It is better to first close the window as it is a liveview tool
		 * TODO: could we improve this behavior?
		 */
			histo_close(TRUE, FALSE);
			siril_close_dialog("histogram_dialog");

		/* apply the process */
			gtk_toggle_button_set_active(seq_button, FALSE);
			apply_mtf_to_sequence(args);
		} else if (invocation == GHT_STRETCH) {
			struct ght_data *args = malloc(sizeof(struct ght_data));
			struct ght_params *params = malloc(sizeof(struct ght_params));
			params->B = _B;
			params->D = _D;
			params->LP = _LP;
			params->SP = _SP;
			params->HP = _HP;
			params->BP = _BP;
			params->stretchtype = _stretchtype;
			params->payne_colourstretchmodel = _payne_colourstretchmodel;
			params->do_red = do_channel[0];
			params->do_green = do_channel[1];
			params->do_blue = do_channel[2];
			args->params_ght = params;
			const gchar* temp = gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryMTFSeq")));
			args->seqEntry = strdup(temp);
			if (args->seqEntry && args->seqEntry[0] == '\0')
				args->seqEntry = strdup("stretch_");
			args->seq = &com.seq;
		/* here it is a bit tricky.
		 * It is better to first close the window as it is a liveview tool
		 * TODO: could we improve this behavior?
		 */
			histo_close(TRUE, FALSE);
			siril_close_dialog("histogram_dialog");

		/* apply the process */
			gtk_toggle_button_set_active(seq_button, FALSE);
			apply_ght_to_sequence(args);
		}
	} else {
		// the apply button resets everything after recomputing with the current values
		fit = &gfit;
		if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview"))) || gui.roi.active) {
			copy_backup_to_gfit();
			if (_payne_colourstretchmodel == COL_SAT) {
				clear_hsl();
				setup_hsl();
			}
			histo_recompute();
		}
		populate_roi();
		// partial cleanup
		if (invocation == HISTO_STRETCH) {
			siril_log_message(_("Applying MTF with values %f, %f, %f\n"),
				_shadows, _midtones, _highlights);
			siril_debug_print("Applying histogram (mid=%.3f, lo=%.3f, hi=%.3f)\n",
				_midtones, _shadows, _highlights);
			undo_save_state(get_preview_gfit_backup(),
				_("Histogram Transf. (mid=%.3f, lo=%.3f, hi=%.3f)"),
				_midtones, _shadows, _highlights);
		} else if (invocation == GHT_STRETCH) {
			siril_debug_print("Applying generalised hyperbolic stretch (D=%2.3f, B=%2.3f, LP=%2.3f, SP=%2.3f, HP=%2.3f", _D, _B, _LP, _SP, _HP);
			if (_payne_colourstretchmodel != COL_SAT) {
				switch (_stretchtype) {
					case STRETCH_PAYNE_NORMAL:
						undo_save_state(get_preview_gfit_backup(), _("GHS pivot: %.3f, amount: %.2f, local: %.2f [%.2f %.2f]"), _SP, _D, _B, _LP, _HP);
						break;
					case STRETCH_PAYNE_INVERSE:
						undo_save_state(get_preview_gfit_backup(), _("GHS INV pivot: %.3f, amount: %.2f, local: %.2f [%.2f %.2f]"), _SP, _D, _B, _LP, _HP);
						break;
					case STRETCH_ASINH:
						undo_save_state(get_preview_gfit_backup(), _("GHS ASINH pivot: %.3f, amount: %.2f [%.2f %.2f]"), _SP, _D, _LP, _HP);
						break;
					case STRETCH_INVASINH:
						undo_save_state(get_preview_gfit_backup(), _("GHS ASINH INV pivot: %.3f, amount: %.2f [%.2f %.2f]"), _SP, _D, _LP, _HP);
						break;
					case STRETCH_LINEAR:
						undo_save_state(get_preview_gfit_backup(), _("GHS LINEAR BP: %.2f"), _BP);
						break;
				}
			} else {
				switch (_stretchtype) {
					case STRETCH_PAYNE_NORMAL:
						undo_save_state(get_preview_gfit_backup(), _("GHS SAT pivot: %.3f, amount: %.2f, local: %.2f [%.2f %.2f]"), _SP, _D, _B, _LP, _HP);
						break;
					case STRETCH_PAYNE_INVERSE:
						undo_save_state(get_preview_gfit_backup(), _("GHS INV SAT pivot: %.3f, amount: %.2f, local: %.2f [%.2f %.2f]"), _SP, _D, _B, _LP, _HP);
						break;
					case STRETCH_ASINH:
						undo_save_state(get_preview_gfit_backup(), _("GHS ASINH SAT pivot: %.3f, amount: %.2f [%.2f %.2f]"), _SP, _D, _LP, _HP);
						break;
					case STRETCH_INVASINH:
						undo_save_state(get_preview_gfit_backup(), _("GHS ASINH INV SAT pivot: %.3f, amount: %.2f [%.2f %.2f]"), _SP, _D, _LP, _HP);
						break;
					case STRETCH_LINEAR:
						undo_save_state(get_preview_gfit_backup(), _("GHS LINEAR SAT BP: %.2f"), _BP);
						break;
				}
			}
		}
		clear_backup();
		clear_hist_backup();
		// reinit
		histo_startup();
		reset_cursors_and_values(FALSE);

		set_cursor("default");
	}
}

void apply_histo_cancel() {
	set_cursor_waiting(TRUE);
	histo_close(TRUE, TRUE);
	set_cursor_waiting(FALSE);
}

void on_histoSpinZoom_value_changed(GtkRange *range, gpointer user_data) {
	adjust_histogram_vport_size();
	queue_window_redraw();
}

void on_histoToolAutoStretch_clicked(GtkToolButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	/* we always apply this function on original data */
	struct mtf_params params = { .shadows = _shadows, .midtones = _midtones, .highlights = _highlights, .do_red = TRUE, .do_green = TRUE, .do_blue = TRUE };
	if (!find_linked_midtones_balance_default(get_preview_gfit_backup(), &params)) {
		_shadows = params.shadows;
		_midtones = params.midtones;
		_highlights = 1.0f;
		_update_entry_text();
		update_histo_mtf();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = histo_update_preview;
		param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
		if (fit->color_managed && !profiles_identical(fit->icc_profile, gui.icc.monitor))
			auto_display_compensation = TRUE;
		notify_update((gpointer) param);
	} else {
		siril_log_color_message(_("Could not compute autostretch parameters, using default values\n"), "salmon");
	}
	set_cursor_waiting(FALSE);
}

void setup_histo_dialog() {
			gtk_window_set_title(GTK_WINDOW(lookup_widget("histogram_dialog")),_("Histogram Transformation"));

			// Hide UI elements not required by histogram stretch
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("box_ghtcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_button")), FALSE);
			// Make visible the UI elements required by histogram stretch
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("histoToolAutoStretch")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("grid32")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("histo_clip_settings")), FALSE); // MTF linked stretch cannot clip anyway

			//Force dialog size to size of visible widgets only
			gtk_window_resize(GTK_WINDOW(lookup_widget("histogram_dialog")),1,1);
}

void setup_ght_dialog() {
			gtk_window_set_title(GTK_WINDOW(lookup_widget("histogram_dialog")),_("Generalised Hyperbolic Stretch Transformations"));
			// Hide the histogram controls
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("grid32")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("histoToolAutoStretch")), FALSE);
			// Make visible and configure the GHT controls
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("box_ghtcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_button")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("histo_clip_settings")), (gfit.naxes[2] == 3));
			GtkWidget *w;
			if (com.pref.gui.combo_theme == 0) {
				w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/eyedropper_dark.svg");
			} else {
				w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/eyedropper.svg");
			}
			gtk_button_set_image(GTK_BUTTON(lookup_widget("eyedropper_button")), w);
			gtk_widget_show(w);
			// Set default parameters
			_D = 0.0f;
			_B = 0.0f;
			_BP = 0.0f;
			_LP = 0.0f;
			_SP = 0.0f;
			_HP = 1.0f;
			lp_warning_given = FALSE;
			hp_warning_given = FALSE;
			gtk_widget_set_tooltip_text(GTK_WIDGET(lookup_widget("drawingarea_histograms")), _("Clicking on the histogram sets SP"));
			gtk_window_resize(GTK_WINDOW(lookup_widget("histogram_dialog")),1,1);
}

void updateGHTcontrols() {
	if (fit->naxis == 2) {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtColourmodelcontrols")), FALSE);
	} else {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtColourmodelcontrols")), TRUE);
	}
	switch (_stretchtype) {
		case STRETCH_LINEAR:
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtDcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtSPcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtHPcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_button")), FALSE);

			break;
		case STRETCH_PAYNE_NORMAL:
		case STRETCH_PAYNE_INVERSE:
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtDcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtSPcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtHPcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_button")), TRUE);

			break;
		case STRETCH_ASINH:
		case STRETCH_INVASINH:
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtDcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtSPcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtLPcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtHPcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("ghtBPcontrols")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_button")), TRUE);
			break;

			// This resizes the window to the smallest size that fits the visible widgets
			gtk_window_resize(GTK_WINDOW(lookup_widget("histogram_dialog")),1,1);
	}
}

void toggle_histogram_window_visibility(int _invocation) {
	siril_close_preview_dialogs();

	invocation = _invocation;
	for (int i=0;i<3;i++) {
		do_channel[i] = TRUE;
	}
	icc_auto_assign_or_convert(&gfit, ICC_ASSIGN_ON_STRETCH);

	if (gtk_widget_get_visible(lookup_widget("histogram_dialog"))) {
		set_cursor_waiting(TRUE);
		histo_close(TRUE, TRUE);
		set_cursor_waiting(FALSE);
		siril_close_dialog("histogram_dialog");
	} else {
		reset_cursors_and_values(TRUE);
		copy_gfit_to_backup();
		//Ensure the colour stretch model is initialised at startup
		_payne_colourstretchmodel = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_payne_colour_stretch_model")));
		if (invocation == HISTO_STRETCH) {
			setup_histo_dialog();
		} else if (invocation == GHT_STRETCH) {
			setup_ght_dialog();
			updateGHTcontrols();
		}
		if (gui.rendering_mode == LINEAR_DISPLAY)
			setup_stretch_sliders(); // In linear mode, set sliders to 0 / 65535
		siril_open_dialog("histogram_dialog");
	}
}

gboolean on_drawingarea_histograms_motion_notify_event(GtkWidget *widget, GdkEventMotion *event,
		gpointer user_data) {
	int width = get_width_of_histo();
	int height = get_height_of_histo();

	if (on_gradient((GdkEvent *) event, width, height)) {
		if (invocation == HISTO_STRETCH)
			set_cursor("grab");
	} else {
		set_cursor("default");
	}

	if (_click_on_histo) {
			gdouble xpos = ((GdkEventButton*) event)->x / (gdouble) width;
			if (xpos < 0.0)
				xpos = 0.0;
			if (xpos > 1.0)
				xpos = 1.0;
		if (invocation == HISTO_STRETCH) {
			gchar *buffer = NULL;
			GtkEntry *histoMidEntry = GTK_ENTRY(lookup_widget("histoMidEntry"));
			GtkEntry *histoShadEntry = GTK_ENTRY(lookup_widget("histoShadEntry"));
			GtkEntry *histoHighEntry = GTK_ENTRY(lookup_widget("histoHighEntry"));

			switch (_type_of_scale) {
			case SCALE_LOW:
				if ((float)xpos > _highlights) {
					_shadows = _highlights;
				} else {
					_shadows = (float)xpos;
				}
				buffer = g_strdup_printf("%.7f", _shadows);
				gtk_entry_set_text(histoShadEntry, buffer);
				break;

			case SCALE_MID:
				if (_highlights == _shadows) {
					_midtones = _highlights;
				} else {
					_midtones = ((float)xpos - _shadows) / (_highlights - _shadows);
				}
				if (_midtones > 1.f) _midtones = 1.f;
				if (_midtones < 0.f) _midtones = 0.f;
				buffer = g_strdup_printf("%.7f", _midtones);
				gtk_entry_set_text(histoMidEntry, buffer);
				break;

			case SCALE_HI:
				if ((float)xpos < _shadows) {
					_shadows =_highlights;
				} else {
					_highlights = (float)xpos;
				}
				buffer = g_strdup_printf("%.7f", _highlights);
				gtk_entry_set_text(histoHighEntry, buffer);
				break;
			}
			set_cursor("grabbing");
			update_histo_mtf();
			g_free(buffer);
		}
	}
	return FALSE;
}

void on_drawingarea_histograms_leave_notify_event(GtkWidget *widget,
		GdkEvent *event, gpointer user_data) {
	set_cursor("default");
}

gboolean on_drawingarea_histograms_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	int width = 339;
	width = get_width_of_histo();
	int height = 281;
	height = get_height_of_histo();
	if (invocation == HISTO_STRETCH) {
		if (on_gradient((GdkEvent *) event, width, height)) {
			float delta = ((_highlights - _shadows) * _midtones) + _shadows;

			_click_on_histo = TRUE;
			gdouble xpos = ((GdkEventButton*) event)->x / (gdouble) width;
			if (fabsf((float)xpos - _highlights) < fabsf((float)xpos - _shadows) && fabsf((float)xpos - _highlights) < fabsf((float)xpos - delta)) {
				_type_of_scale = SCALE_HI;
			} else if (fabsf((float)xpos - _shadows) < fabsf((float)xpos - delta) && fabsf((float)xpos - _shadows) < fabsf((float)xpos - _highlights)) {
				_type_of_scale = SCALE_LOW;
			} else if ((_shadows == _highlights) && (_shadows > 0.f)) {
				_type_of_scale = SCALE_LOW;
			} else if ((_shadows == _highlights) && (_shadows == 0.f)) {
				_type_of_scale = SCALE_HI;
			} else {
				_type_of_scale = SCALE_MID;
			}
			set_cursor("grabbing");
		}
	} else if (invocation == GHT_STRETCH) {
		gdouble xpos = ((GdkEventButton*) event)->x / (gdouble) width;
		// Set the symmetry point by clicking on the histogram
		// Note that we need to carry out the reverse transform on xpos in order
		// to set SP to the value in unstretched space that when stretched will
		// give the clicked position. So the calculation depends on stretchtype.
		switch (_stretchtype) {
			case STRETCH_LINEAR:
				invxpos = xpos; // For linear mode SP is irrelevant so it is more useful
								// to click to set BP
				break;
			case STRETCH_PAYNE_NORMAL:
				GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, STRETCH_PAYNE_INVERSE);
				invxpos = GHT(xpos, _B, _D, _LP, _SP, _HP, _BP, STRETCH_PAYNE_INVERSE, &compute_params);
				GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, STRETCH_PAYNE_NORMAL);
				break;
			case STRETCH_PAYNE_INVERSE:
				GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, STRETCH_PAYNE_NORMAL);
				invxpos = GHT(xpos, _B, _D, _LP, _SP, _HP, _BP, STRETCH_PAYNE_NORMAL, &compute_params);
				GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, STRETCH_PAYNE_INVERSE);
				break;
			case STRETCH_ASINH:
				GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, STRETCH_INVASINH);
				invxpos = GHT(xpos, _B, _D, _LP, _SP, _HP, _BP, STRETCH_INVASINH, &compute_params);
				GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, STRETCH_ASINH);
				break;
			case STRETCH_INVASINH:
				GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, STRETCH_ASINH);
				invxpos = GHT(xpos, _B, _D, _LP, _SP, _HP, _BP, STRETCH_ASINH, &compute_params);
				GHTsetup(&compute_params, _B, _D, _LP, _SP, _HP, STRETCH_INVASINH);
				break;
		}
	}
	return FALSE;
}

gboolean on_drawingarea_histograms_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	set_cursor("default");
	if ((invocation == GHT_STRETCH) && (invxpos != -1)) {
		if (_stretchtype == STRETCH_LINEAR)
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtBP")),invxpos);
		else
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtSP")),invxpos);
		invxpos = -1.0;
	}
	if (_click_on_histo) {
		_click_on_histo = FALSE;
		set_cursor_waiting(TRUE);
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = histo_update_preview;
		param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
		notify_update((gpointer) param);
		set_cursor_waiting(FALSE);
	}
	return FALSE;
}

void on_spin_ghtD_value_changed(GtkSpinButton *button, gpointer user_data) {
	_D = (float) expm1(gtk_spin_button_get_value(button));
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_spin_ghtB_value_changed(GtkSpinButton *button, gpointer user_data) {
	_B = (float) gtk_spin_button_get_value(button);
	if (fabsf(_B) < 1.e-3f) {
		_B = 0.f;
		gtk_spin_button_set_value(button, 0.f);
	}
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_spin_ghtLP_value_changed(GtkSpinButton *button, gpointer user_data) {
	_LP = (float) gtk_spin_button_get_value(button);
	if (_LP > _SP) {
		gtk_spin_button_set_value(button,_SP);
		if (!lp_warning_given) { // Prevent spamming the warning log message
			siril_log_message(_("Shadow preservation point cannot be set higher than stretch focal point\n"));
			lp_warning_given = TRUE;
		}
	}
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_spin_ghtSP_value_changed(GtkSpinButton *button, gpointer user_data) {
	GtkSpinButton *spin_HP = GTK_SPIN_BUTTON(lookup_widget("spin_ghtHP"));
	GtkSpinButton *spin_LP = GTK_SPIN_BUTTON(lookup_widget("spin_ghtLP"));
	_SP = (float) gtk_spin_button_get_value(button);
	if (_SP < _LP) {
		gtk_spin_button_set_value(spin_LP, (double) _SP);
	}
	if (_SP > _HP) {
		gtk_spin_button_set_value(spin_HP, (double) _SP);
	}
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_eyedropper_SP_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_HP = GTK_SPIN_BUTTON(lookup_widget("spin_ghtHP"));
	GtkSpinButton *spin_LP = GTK_SPIN_BUTTON(lookup_widget("spin_ghtLP"));
	int chan;
	imstats* stats[3];
	float ref = 0.0f;
	float norm = 1.0f;
	if (get_preview_gfit_backup()->type == DATA_USHORT)
		norm = USHRT_MAX_DOUBLE;

	if (!com.selection.w || !com.selection.h) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the point or area to set SP"));
		return;
	}
	int active_channels = 0;
	for (chan = 0; chan < get_preview_gfit_backup()->naxes[2]; chan++) {
		if (do_channel[chan]) {
			stats[chan] = statistics(NULL, -1, get_preview_gfit_backup(), chan, &com.selection, STATS_BASIC, MULTI_THREADED);
			if (!stats[chan]) {
				siril_log_message(_("Error: statistics computation failed.\n"));
				return;
			}
			fprintf(stdout,"channel mean: %f\n", stats[chan]->mean);
			ref += stats[chan]->mean;
			free_stats(stats[chan]);
			active_channels++;
		}
	}
	ref /= (active_channels ? active_channels : 1.f);
	ref /= norm;
	_SP = (float) ref;
	if (_SP < _LP) {
		gtk_spin_button_set_value(spin_LP, (double) _SP);
	}
	if (_SP > _HP) {
		gtk_spin_button_set_value(spin_HP, (double) _SP);
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtSP")), (double) _SP);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_spin_ghtHP_value_changed(GtkSpinButton *button, gpointer user_data) {
	_HP = (float) gtk_spin_button_get_value(button);
	if (_HP < _SP) {
		gtk_spin_button_set_value(button, _SP);
		if (!hp_warning_given) { // Prevent spamming the warning log message
			siril_log_message(_("Highlight preservation point cannot be set lower than stretch focal point\n"));
			hp_warning_given = TRUE;
		}
	}
	update_histo_mtf();
	queue_window_redraw();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_spin_ghtBP_value_changed(GtkSpinButton *button, gpointer user_data) {
	_BP = (float) gtk_spin_button_get_value(button);
	if (_stretchtype == STRETCH_LINEAR) {
		update_histo_mtf();
		queue_window_redraw();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = histo_update_preview;
		param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
		notify_update((gpointer) param);
	}
}

void on_histo_clip_mode_changed(GtkComboBox *combo, gpointer user_data) {
	_clip_mode = (clip_mode_t) gtk_combo_box_get_active(combo);
	update_histo_mtf();
	queue_window_redraw();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_payneType_changed(GtkComboBox *combo, gpointer user_data) {
	_stretchtype = gtk_combo_box_get_active(combo);
		if (_stretchtype == STRETCH_LINEAR)
			gtk_widget_set_tooltip_text(GTK_WIDGET(lookup_widget("drawingarea_histograms")), _("Clicking on the linear stretch histogram sets BP"));
		else
			gtk_widget_set_tooltip_text(GTK_WIDGET(lookup_widget("drawingarea_histograms")), _("Clicking on the histogram sets SP"));
	updateGHTcontrols();
	update_histo_mtf();
	queue_window_redraw();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
}

void on_payne_colour_stretch_model_changed(GtkComboBox *combo, gpointer user_data) {
	int tmp = gtk_combo_box_get_active(combo);
	if (tmp == COL_SAT) {
		if (fit->naxes[2] != 3) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
			_("Not all colour channels are selected. Saturation stretch cannot be used: setting independent channels colour model."));
			gtk_combo_box_set_active(combo, COL_INDEP);
			_payne_colourstretchmodel = COL_INDEP;
		} else {
			set_cursor_waiting(TRUE);
			reset_cursors_and_values(FALSE);
			histo_close(TRUE, TRUE);
			setup_hsl();
			_payne_colourstretchmodel = tmp;
			histo_startup();
			set_cursor_waiting(FALSE);
			return;
		}
	} else {
		if (!(do_channel[0] && do_channel[1] && do_channel[2])) {
			if (tmp == COL_HUMANLUM) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Human luminance colour model cannot be used: setting even weighted luminance colour model."));
			gtk_combo_box_set_active(combo, COL_EVENLUM);
			_payne_colourstretchmodel = COL_EVENLUM;
			}
		}
		if (_payne_colourstretchmodel == COL_SAT) {
			set_cursor_waiting(TRUE);
			reset_cursors_and_values(FALSE);
			histo_close(TRUE, TRUE);
			clear_hsl();
			_payne_colourstretchmodel = tmp;
			histo_startup();
			set_cursor_waiting(FALSE);
			return;
		} else {
			_payne_colourstretchmodel = tmp;
			update_histo_mtf();
			queue_window_redraw();
			update_image *param = malloc(sizeof(update_image));
			param->update_preview_fn = histo_update_preview;
			param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
			notify_update((gpointer) param);
			return;
		}
	}
}

gboolean on_histoShadEntry_focus_out_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("histoShadEntry"));
	float lo = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
	if (lo <= 0.f) lo = 0.f;
	_shadows = lo;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", lo);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);

	return FALSE;
}

void on_histoShadEntry_activate(GtkEntry *entry, gpointer user_data) {
	float lo = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
	if (lo <= 0.f) lo = 0.f;
	_shadows = lo;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", lo);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);
}

void on_histoMidEntry_activate(GtkEntry *entry, gpointer user_data) {
	float mid = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
	_midtones = mid;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", mid);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);
}

gboolean on_histoMidEntry_focus_out_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {

	GtkEntry *entry = GTK_ENTRY(lookup_widget("histoMidEntry"));
	float mid = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
	_midtones = mid;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", mid);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);

	return FALSE;
}

gboolean on_histoHighEntry_focus_out_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("histoHighEntry"));
	float hi = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
//	if (hi <= _shadows) hi = _shadows;
	if (hi >= 1.f) hi = 1.f;
	_highlights = hi;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", hi);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);

	return FALSE;
}

void on_histoHighEntry_activate(GtkEntry *entry, gpointer user_data) {
	float hi = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
//	if (hi <= _shadows) hi = _shadows;
	if (hi >= 1.f) hi = 1.f;
	_highlights = hi;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", hi);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);
}

int mtf_finalize_hook(struct generic_seq_args *args) {
	struct mtf_data *data = (struct mtf_data *) args->user;
	int retval = seq_finalize_hook(args);
	free(data);
	return retval;
}

void apply_mtf_to_sequence(struct mtf_data *mtf_args) {
	struct generic_seq_args *args = create_default_seqargs(mtf_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = mtf_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = mtf_finalize_hook;
	args->image_hook = mtf_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Midtone Transfer Function");
	args->has_output = TRUE;
	args->new_seq_prefix = mtf_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->user = mtf_args;

	mtf_args->fit = NULL;	// not used here

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(mtf_args);
		free_generic_seq_args(args);
	}
}

int ght_finalize_hook(struct generic_seq_args *args) {
	struct ght_data *data = (struct ght_data *) args->user;
	struct ght_params *ghtp = (struct ght_params *) data->params_ght;
	int retval = seq_finalize_hook(args);
	free(ghtp);
	free(data);
	clear_hsl();
	sequence_working = FALSE;
	return retval;
}

void apply_ght_to_sequence(struct ght_data *ght_args) {
	struct generic_seq_args *args = create_default_seqargs(ght_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = ght_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = ght_finalize_hook;
	args->image_hook = ght_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Generalised Hyperbolic Transfer Function");
	args->has_output = TRUE;
	args->new_seq_prefix = ght_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->user = ght_args;

	ght_args->fit = NULL;	// not used here

	if(!start_in_new_thread(generic_sequence_worker, args)) {
		free(ght_args->seqEntry);
		free(ght_args);
		free_generic_seq_args(args);
	}
}

void on_histo_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckPreview")))) {
		/* if user click very fast */
		waiting_for_thread();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = histo_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
