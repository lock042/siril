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

#include <gsl/gsl_histogram.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "algos/colors.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/utils.h"	// for lookup_widget()
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/siril_preview.h"
#include "core/undo.h"
#include "histogram.h"
#include "histogram_utils.h"

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

static fits* fit = NULL;

static gboolean auto_display_compensation = FALSE;

static double invxpos = -1.0;
// static float graph_height = 0.f;	// the max value of all bins
static guint64 clipped[] = { 0, 0 };

// Original ICC profile, in case we don't apply a stretch and need to revert
static cmsHPROFILE original_icc = NULL;
static gboolean single_image_stretch_applied = FALSE;

// The 4th toggle pointer remains NULL, but it makes the logic easier in redraw_histo()
static GtkToggleButton *toggles[MAXVPORT] = { NULL };
static GtkToggleButton *toggleGrid = NULL, *toggleCurve = NULL, *toggleOrig = NULL;

/* the original histogram, used as starting point of each computation */
static gsl_histogram *hist_backup[MAXVPORT] = { NULL, NULL, NULL, NULL };
static gsl_histogram *hist_sat_backup = NULL;

static float _midtones, _shadows, _highlights;

static gboolean _click_on_histo = FALSE;
static ScaleType _type_of_scale;
static gboolean lp_warning_given, hp_warning_given = FALSE;

void *huebuf = NULL, *satbuf_orig = NULL, *satbuf_working = NULL, *lumbuf = NULL;

static void setup_hsl();
static void clear_hsl();
static void stretch_dialog_compute_histograms(fits *thefit);

/* create_mtf_data, destroy_mtf_data, create_ght_data, destroy_ght_data
 * moved to filters/mtf.c and filters/ght.c */

static int get_width_of_histo() {
	return gtk_widget_get_width(lookup_widget("drawingarea_histograms"));
}

static int get_height_of_histo() {
	return gtk_widget_get_height(lookup_widget("drawingarea_histograms"));
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

gboolean redraw_histo(GtkWidget *widget, cairo_t *cr, gpointer data);

/* GTK4 draw_func adapter — forwards to the legacy redraw_histo. */
static void histo_draw_cb(GtkDrawingArea *area, cairo_t *cr,
                          int width, int height, gpointer data) {
	(void) width; (void) height;
	redraw_histo(GTK_WIDGET(area), cr, data);
}

/* Forward decls so init_toggles() can reference the GTK4 event-controller
 * callbacks defined further down; their bodies depend on update_histo_mtf()
 * and on_gradient() which sit between init_toggles and the callback set. */
static void histo_press_cb(GtkGestureClick *gesture, int n_press,
                           double x, double y, gpointer user_data);
static void histo_release_cb(GtkGestureClick *gesture, int n_press,
                             double x, double y, gpointer user_data);
static void histo_motion_cb(GtkEventControllerMotion *controller,
                            double x, double y, gpointer user_data);

static void init_toggles() {
	if (!toggles[0]) {
		fit = gfit;
		toggles[0] = GTK_TOGGLE_BUTTON(lookup_widget("histoToolRed"));
		toggles[1] = GTK_TOGGLE_BUTTON(lookup_widget("histoToolGreen"));
		toggles[2] = GTK_TOGGLE_BUTTON(lookup_widget("histoToolBlue"));
		do_channel[0] = siril_toggle_get_active(GTK_WIDGET(toggles[0]));
		do_channel[1] = siril_toggle_get_active(GTK_WIDGET(toggles[1]));
		do_channel[2] = siril_toggle_get_active(GTK_WIDGET(toggles[2]));
		toggleGrid = GTK_TOGGLE_BUTTON(lookup_widget("histoToolGrid"));
		toggleCurve = GTK_TOGGLE_BUTTON(lookup_widget("histoToolCurve"));
		toggleOrig = GTK_TOGGLE_BUTTON(lookup_widget("histoToolOrig"));
		/* GTK4: GtkDrawingArea has no "draw" signal — wire the draw
		 * function explicitly.  Phase 18 stripped the .ui binding. */
		GtkWidget *drawarea = lookup_widget("drawingarea_histograms");
		if (drawarea) {
			gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawarea),
			                               histo_draw_cb, NULL, NULL);

			/* GTK4 event-controller wiring for the drawingarea sliders.
			 * GtkGestureClick replaces button-press/release-event;
			 * GtkEventControllerMotion replaces motion-notify-event. */
			GtkGesture *click = gtk_gesture_click_new();
			gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click),
			                              GDK_BUTTON_PRIMARY);
			g_signal_connect(click, "pressed",  G_CALLBACK(histo_press_cb),  NULL);
			g_signal_connect(click, "released", G_CALLBACK(histo_release_cb), NULL);
			gtk_widget_add_controller(drawarea, GTK_EVENT_CONTROLLER(click));

			GtkEventController *motion = gtk_event_controller_motion_new();
			g_signal_connect(motion, "motion", G_CALLBACK(histo_motion_cb), NULL);
			gtk_widget_add_controller(drawarea, motion);
		}
	}
}

static gboolean active_sliders = TRUE;

static void set_controls_active(gboolean state) {
	gtk_widget_set_sensitive(lookup_widget("histoShadEntry"), state);
	gtk_widget_set_sensitive(lookup_widget("histoMidEntry"), state);
	gtk_widget_set_sensitive(lookup_widget("histoHighEntry"), state);
	active_sliders = state;
}

static void histo_startup() {
	set_controls_active(TRUE);
	add_roi_callback(histo_change_between_roi_and_image);
	roi_supported(TRUE);
	copy_gfit_to_backup();
	if (fit->naxes[2] == 3 && _payne_colourstretchmodel == COL_SAT)
		setup_hsl();

	init_toggles();
	do_channel[0] = siril_toggle_get_active(GTK_WIDGET(toggles[0]));
	do_channel[1] = siril_toggle_get_active(GTK_WIDGET(toggles[1]));
	do_channel[2] = siril_toggle_get_active(GTK_WIDGET(toggles[2]));
	if (fit->naxes[2] == 3 && !(do_channel[0] && do_channel[1] && do_channel[2])) {
		if (_payne_colourstretchmodel == COL_HUMANLUM && invocation == GHT_STRETCH) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Human luminance colour model cannot be used: setting even weighted luminance colour model."));
			_payne_colourstretchmodel = COL_EVENLUM;
			gtk_drop_down_set_selected(GTK_DROP_DOWN(lookup_widget("combo_payne_colour_stretch_model")), COL_EVENLUM);
		}
		if (_payne_colourstretchmodel == COL_SAT && invocation == GHT_STRETCH) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Saturation stretch cannot be used: setting independent channels colour model."));
			_payne_colourstretchmodel = COL_INDEP;
			gtk_drop_down_set_selected(GTK_DROP_DOWN(lookup_widget("combo_payne_colour_stretch_model")), COL_INDEP);
		}
	}
	// also get the backup histogram
	stretch_dialog_compute_histograms(fit);
	for (int i = 0; i < fit->naxes[2]; i++)
		hist_backup[i] = gsl_histogram_clone(com.layers_hist[i]);
	if (com.sat_hist) {
		if (hist_sat_backup)
			gsl_histogram_free(hist_sat_backup);
		hist_sat_backup = gsl_histogram_clone(com.sat_hist);
	}
}

static void histo_close(gboolean revert, gboolean update_image_if_needed, gboolean revert_icc_profile) {
	if (revert) {

		for (int i = 0; i < fit->naxes[2]; i++) {
			set_histogram(hist_backup[i], i);
			hist_backup[i] = NULL;
		}
		set_sat_histogram(hist_sat_backup);
		hist_sat_backup = NULL;
		if (is_preview_active() && !copy_backup_to_gfit() && update_image_if_needed) {
			set_cursor_waiting(TRUE);
			notify_gfit_data_modified();
			gfit_modified_update_gui();
		}
	}
	// free data
	if(!sequence_working && !revert) {
		clear_hsl();
	}
	if (revert_icc_profile && !single_image_stretch_applied) {
		if (gfit->icc_profile)
			cmsCloseProfile(gfit->icc_profile);
		gfit->icc_profile = copyICCProfile(original_icc);
		color_manage(gfit, gfit->icc_profile != NULL);
	}
	clear_backup();
	clear_hist_backup();
	roi_supported(FALSE);
	remove_roi_callback(histo_change_between_roi_and_image);
}


static void fit_to_hsl() {
	size_t npixels = fit->rx * fit->ry;
	if (fit->type == DATA_FLOAT) {
		float* hf = (float*) huebuf;
		float* sf = (float*) satbuf_orig;
		float* lf = (float*) lumbuf;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
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
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0 ; i < npixels ; i++) {
			rgbw_to_hslw(fit->pdata[0][i], fit->pdata[1][i], fit->pdata[2][i], &hw[i], &sw[i], &lw[i]);
		}
		memcpy(satbuf_working, satbuf_orig, npixels * sizeof(WORD));
	}
}

static void histo_recompute(gboolean for_preview) {
	set_cursor("progress");
	gboolean auto_comp = auto_display_compensation;
	auto_display_compensation = FALSE;

	// Ensure we're working with the right fit
	fit = gui.roi.active ? &gui.roi.fit : gfit;

	copy_backup_to_gfit();

	// Create appropriate data structure and setup worker args
	struct generic_img_args *args = NULL;

	if (invocation == HISTO_STRETCH) {
		struct mtf_data *data = create_mtf_data();
		if (!data) {
			siril_log_color_message(_("Memory allocation failed.\n"), "red");
			return;
		}

		data->fit = fit;
		data->params.shadows = _shadows;
		data->params.midtones = _midtones;
		data->params.highlights = _highlights;
		data->params.do_red = do_channel[0];
		data->params.do_green = do_channel[1];
		data->params.do_blue = do_channel[2];
		data->auto_display_compensation = auto_comp;
		data->is_preview = TRUE;

		args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			destroy_mtf_data(data);
			siril_log_color_message(_("Memory allocation failed.\n"), "red");
			return;
		}

		args->fit = fit;
		args->mem_ratio = 1.0f;
		args->image_hook = mtf_single_image_hook;
		args->log_hook = mtf_log_hook;
		args->idle_function = mtf_single_image_idle;
		args->description = _("Midtone Transfer Function");
		args->verbose = FALSE;
		args->user = data;
		args->max_threads = com.max_thread;
		args->for_preview = TRUE;
		args->custom_undo = TRUE;
		args->mask_aware = TRUE;
		args->for_roi = gui.roi.active;

	} else if (invocation == GHT_STRETCH) {
		struct ght_data *data = create_ght_data();
		if (!data) {
			siril_log_color_message(_("Memory allocation failed.\n"), "red");
			return;
		}

		data->fit = fit;
		data->params_ght = malloc(sizeof(struct ght_params));
		if (!data->params_ght) {
			destroy_ght_data(data);
			siril_log_color_message(_("Memory allocation failed.\n"), "red");
			return;
		}

		data->params_ght->B = _B;
		data->params_ght->D = _D;
		data->params_ght->LP = _LP;
		data->params_ght->SP = _SP;
		data->params_ght->HP = _HP;
		data->params_ght->BP = _BP;
		data->params_ght->stretchtype = _stretchtype;
		data->params_ght->payne_colourstretchmodel = _payne_colourstretchmodel;
		data->params_ght->do_red = do_channel[0];
		data->params_ght->do_green = do_channel[1];
		data->params_ght->do_blue = do_channel[2];
		data->params_ght->clip_mode = _clip_mode;
		data->auto_display_compensation = auto_comp;
		data->is_preview = TRUE;

		args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			destroy_ght_data(data);
			siril_log_color_message(_("Memory allocation failed.\n"), "red");
			return;
		}

		args->fit = fit;
		args->mem_ratio = (_payne_colourstretchmodel == COL_SAT) ? 2.0f : 1.0f;
		args->image_hook = ght_single_image_hook;
		args->log_hook = ght_log_hook;
		args->idle_function = ght_single_image_idle;
		args->description = _("Generalized Hyperbolic Transform");
		args->verbose = FALSE;
		args->user = data;
		args->max_threads = com.max_thread;
		args->for_preview = TRUE;
		args->mask_aware = TRUE;
		args->for_roi = gui.roi.active;
	}

	if (args) {
		start_in_new_thread(generic_image_worker, args);
	}
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
	gtk_editable_set_text(GTK_EDITABLE(clip_low), "0.000%");
	gtk_editable_set_text(GTK_EDITABLE(clip_high), "0.000%");
}

static void _update_entry_text() {
	GtkEntry *histoMidEntry = GTK_ENTRY(lookup_widget("histoMidEntry"));
	GtkEntry *histoShadEntry = GTK_ENTRY(lookup_widget("histoShadEntry"));
	GtkEntry *histoHighEntry = GTK_ENTRY(lookup_widget("histoHighEntry"));
	gchar *buffer;

	buffer = g_strdup_printf("%.7f", _shadows);
	gtk_editable_set_text(GTK_EDITABLE(histoShadEntry), buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", _highlights);
	gtk_editable_set_text(GTK_EDITABLE(histoHighEntry), buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", _midtones);
	gtk_editable_set_text(GTK_EDITABLE(histoMidEntry), buffer);
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
	gtk_editable_set_text(GTK_EDITABLE(clip_high), buffer);
	tmp = (double)clipped[0] * 100.0 / (double)data;
	g_snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_editable_set_text(GTK_EDITABLE(clip_low), buffer);
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
		gtk_button_set_child(GTK_BUTTON(toggles[0]), w);
		gtk_widget_set_visible(w, TRUE);
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
		gtk_button_set_child(GTK_BUTTON(toggles[0]), w);
		gtk_widget_set_visible(w, TRUE);
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

	cur_width = gtk_widget_get_width(vport);
	cur_height = gtk_widget_get_height(vport);
	targetW = (int) (((double)cur_width) * zoomH);
	targetH = (int) (((double)cur_height) * zoomV);
	gtk_widget_set_size_request(drawarea, targetW, targetH);
#ifdef HISTO_DEBUG
	fprintf(stdout, "Histo vport size (%d, %d)\n", targetW, targetH);
#endif
}

// create a new histogram object for the passed float buffer (used for sat)
gsl_histogram* computeHistoSat(void* buf) {
	size_t i, ndata, size;

	size = get_histo_size(gfit);
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
	static GtkCheckButton *HistoCheckLogButton = NULL;

	if (HistoCheckLogButton == NULL)
		HistoCheckLogButton = GTK_CHECK_BUTTON(lookup_widget("HistoCheckLogButton"));
	return (siril_toggle_get_active(GTK_WIDGET(HistoCheckLogButton)));
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

// erase image and redraw the background color and grid
static void erase_histo_display(cairo_t *cr, int width, int height) {
	gboolean drawGrid, drawCurve;
	// clear all with background color
	fill_histo_background(cr, width, height);

	// draw grid
	init_toggles();
	drawGrid = siril_toggle_get_active(GTK_WIDGET(toggleGrid));
	drawCurve = siril_toggle_get_active(GTK_WIDGET(toggleCurve));

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

static void display_scale(cairo_t *cr, int width, int height) {
	draw_gradient(cr, width, height);
	if (invocation == HISTO_STRETCH) {
		float delta = ((_highlights - _shadows) * _midtones) + _shadows;
		draw_slider(cr, width, height, _shadows * width);
		draw_slider(cr, width, height, (delta) * width);
		draw_slider(cr, width, height, _highlights * width);
	}
}

static void apply_mtf_to_histo(gsl_histogram *histo, float norm,
		float m, float lo, float hi) {

	size_t int_norm = (size_t)norm;
	gsl_histogram *mtf_histo = gsl_histogram_alloc(int_norm + 1);
	gsl_histogram_set_ranges_uniform(mtf_histo, 0, norm);

	// Precompute loop-invariant values
	size_t lo_bin = (size_t)roundf_to_WORD(lo * norm);
	size_t hi_bin = (size_t)roundf_to_WORD(hi * norm);
	float inv_norm = 1.f / norm;
	float inv_range = (hi > lo) ? 1.f / (hi - lo) : 1.f;
	float m_minus_1 = m - 1.f;
	float two_m_minus_1 = 2.f * m - 1.f;

// disabled because of ISSUE #136 (https://free-astro.org/bugs/view.php?id=136)
// Update by A Knagg-Baugh - the issue only relates to MacOS so I've limited the
// disabling to that OS with a #ifndef. Similarly for the other instance below
#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
#endif
	{
		// Per-thread histogram eliminates the data race on mtf_histo
		gsl_histogram *histo_thr = gsl_histogram_alloc(int_norm + 1);
		gsl_histogram_set_ranges_uniform(histo_thr, 0, norm);
		guint64 clip0 = 0, clip1 = 0;

#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
#endif
		for (size_t i = 0; i <= int_norm; i++) {
			float binval = gsl_histogram_get(histo, i);
			if (binval == 0.0)
				continue;
			size_t bi = i;
			if (i < lo_bin) {
				bi = lo_bin;
				clip0 += (guint64)binval;
			} else if (i > hi_bin) {
				bi = hi_bin;
				clip1 += (guint64)binval;
			}
			float pxl = bi * inv_norm;
			float xp = fmaxf(0.f, fminf((pxl - lo) * inv_range, 1.f));
			WORD mtf = roundf_to_WORD(((m_minus_1 * xp) / (two_m_minus_1 * xp - m)) * norm);
			gsl_histogram_accumulate(histo_thr, (double)mtf, (double)binval);
		}

#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp critical
#endif
#endif
		{
			gsl_histogram_add(mtf_histo, histo_thr);
			clipped[0] += clip0;
			clipped[1] += clip1;
		}
		gsl_histogram_free(histo_thr);
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

	// Precompute loop-invariant values
	size_t lo_bin = (size_t)roundf_to_WORD(lo * norm);
	size_t hi_bin = (size_t)roundf_to_WORD(hi * norm);
	float inv_norm = 1.f / norm;

#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
#endif
	{
		// Per-thread histogram eliminates the data race on mtf_histo
		gsl_histogram *histo_thr = gsl_histogram_alloc(int_norm + 1);
		gsl_histogram_set_ranges_uniform(histo_thr, 0, norm);
		guint64 clip0 = 0, clip1 = 0;

#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
#endif
		for (size_t i = 0; i <= int_norm; i++) {
			float binval = gsl_histogram_get(histo, i);
			if (binval == 0.0)
				continue;
			size_t bi = i;
			if (i < lo_bin) {
				bi = lo_bin;
				clip0 += (guint64)binval;
			} else if (i > hi_bin) {
				bi = hi_bin;
				clip1 += (guint64)binval;
			}
			float pxl = bi * inv_norm;
			WORD ght = roundf_to_WORD(GHT(pxl, _B, _D, _LP, _SP, _HP, _BP, _stretchtype, &compute_params) * norm);
			gsl_histogram_accumulate(histo_thr, (double)ght, (double)binval);
		}

#ifndef __APPLE__
#ifdef _OPENMP
#pragma omp critical
#endif
#endif
		{
			gsl_histogram_add(mtf_histo, histo_thr);
			clipped[0] += clip0;
			clipped[1] += clip1;
		}
		gsl_histogram_free(histo_thr);
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
	siril_toggle_set_active(GTK_WIDGET(toggleGrid), TRUE);
	siril_toggle_set_active(GTK_WIDGET(toggleCurve), TRUE);
	siril_toggle_set_active(GTK_WIDGET(toggleOrig), BOOL_TRUE);
	siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("checkMTFSeq"))), FALSE);
	gtk_editable_set_text(GTK_EDITABLE(lookup_widget("entryMTFSeq")), "stretch_");
//	on_histoZoom100_clicked(NULL, NULL);
	if (full_reset) {
		siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckLogButton"))), (com.pref.gui.display_histogram_mode == LOG_DISPLAY ? TRUE : FALSE));
		for (int i = 0; i < 3; ++i)
			siril_toggle_set_active(GTK_WIDGET(toggles[i]), TRUE);
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
			gtk_drop_down_set_selected(GTK_DROP_DOWN(lookup_widget("combo_payne_colour_stretch_model")), 0);
			gtk_drop_down_set_selected(GTK_DROP_DOWN(lookup_widget("combo_payneTyp")), 0);
			gtk_drop_down_set_selected(GTK_DROP_DOWN(lookup_widget("histo_clip_mode")), RGBBLEND);
			siril_toggle_set_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))), TRUE);
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
	/* Lock against notify_gfit_data_modified() on the worker thread, which
	 * calls invalidate_gfit_histogram() + compute_histo_for_fit() under the
	 * same mutex.  Without this, the worker can nullify layers_hist[1] or
	 * layers_hist[2] between our NULL-check on [0] and the gsl_histogram_memcpy
	 * loop, causing a SIGSEGV. */
	g_mutex_lock(&com.histogram_mutex);
	if (!com.layers_hist[0] || !hist_backup[0]) {
		g_mutex_unlock(&com.histogram_mutex);
		return;
	}
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
		if (com.sat_hist && hist_sat_backup) {
			gsl_histogram_memcpy(com.sat_hist, hist_sat_backup);
			apply_ght_to_histo(com.sat_hist, norm, _SP, _BP, 1.0f);
		}
	}
	g_mutex_unlock(&com.histogram_mutex);
	size_t data = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	_update_clipped_pixels(data);
	queue_window_redraw();
}

static int histo_update_preview() {
	fit = gui.roi.active ? &gui.roi.fit : gfit;
	if (!closing)
		histo_recompute(TRUE);
	return 0;
}

/* Returns TRUE when the mouse is inside the gradient strip drawn at the
 * bottom of the histogram drawingarea — the band where shadow/midtone/
 * highlight markers live and where slider drags should activate. */
static gboolean on_gradient(double x, double y, int width, int height) {
	return x > 0 && x < width
	    && y > height - GRADIENT_HEIGHT
	    && y < height;
}

/*
 * Public functions
 */

static gboolean set_histo_toggles_names_idle(gpointer data G_GNUC_UNUSED) {
	set_histo_toggles_names();
	return G_SOURCE_REMOVE;
}

/* Stretch-dialog-specific histogram refresh: computes per-layer histograms via
 * the shared compute_histo_for_fit(), then also computes the saturation histogram
 * (only relevant in HSL/GHT colour modes) and updates the toggle button labels. */
static void stretch_dialog_compute_histograms(fits *thefit) {
	compute_histo_for_fit(thefit);
	if (thefit->naxes[2] == 3 && satbuf_working)
		set_sat_histogram(computeHistoSat(satbuf_working));
	siril_add_idle(set_histo_toggles_names_idle, NULL);
}

/* call from main thread */
void update_gfit_histogram_if_needed() {
	invalidate_gfit_histogram();
	if (is_histogram_visible()) {
		stretch_dialog_compute_histograms(fit);
		queue_window_redraw();
	}
}

/* Trigger a visual refresh of the histogram window if it is currently visible.
 * Call from the GTK main thread only; histogram data must already be computed. */
void refresh_histogram_if_visible() {
	if (is_histogram_visible())
		queue_window_redraw();
}


gboolean mtf_single_image_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	if (!args) {
		stop_processing_thread();
		return FALSE;
	}

	struct mtf_data *data = (struct mtf_data *)args->user;

	// Check if processing succeeded
	if (args->retval == 0) {
		// Update histograms
		stretch_dialog_compute_histograms(fit);
		queue_window_redraw();

		// Notify that gfit was modified
		gfit_modified_update_gui();

		// If this was a final apply (not preview), clear backups and reinitialize
		if (!data->is_preview) {
			clear_backup();
			clear_hist_backup();

			// Reinitialize for next use
			histo_startup();
			reset_cursors_and_values(FALSE);
		}
	} else {
		siril_log_color_message(_("MTF processing failed.\n"), "red");
	}

	// Cleanup
	free_generic_img_args(args);
	stop_processing_thread();
	return FALSE;
}

gboolean ght_single_image_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	if (!args) {
		stop_processing_thread();
		return FALSE;
	}

	struct ght_data *data = (struct ght_data *)args->user;

	// Check if processing succeeded
	if (args->retval == 0) {
		// Update histograms
		stretch_dialog_compute_histograms(fit);
		queue_window_redraw();

		// Notify that gfit was modified
		gfit_modified_update_gui();

		// If this was a final apply (not preview), clear backups and reinitialize
		if (!data->is_preview) {
			clear_backup();
			clear_hist_backup();

			// Reinitialize for next use
			histo_startup();
			reset_cursors_and_values(FALSE);
		}
	} else {
		siril_log_color_message(_("GHT processing failed.\n"), "red");
	}

	// Cleanup
	free_generic_img_args(args);
	stop_processing_thread();
	return FALSE;
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
	fit = gui.roi.active ? &gui.roi.fit : gfit;
	if (_payne_colourstretchmodel == COL_SAT) {
		clear_hsl();
		setup_hsl();
	}
	update_histo_mtf();
	gui.roi.operation_supports_roi = TRUE;
	// If we are showing the preview, update it after the ROI change.
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
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

	gboolean is_mono = (fit->naxes[2] == 1);
	for (i = 0; i < RGB_VPORT; i++) {
		if (com.layers_hist[i]) {
			if (siril_toggle_get_active(GTK_WIDGET(toggleOrig))) {
				display_histo(hist_backup[i], cr, i, width, height - GRADIENT_HEIGHT, zoomH, zoomV, TRUE, is_log_scale(), is_mono);
			}
			if (!toggles[i] || siril_toggle_get_active(GTK_WIDGET(toggles[i]))) {
				display_histo(com.layers_hist[i], cr, i, width, height - GRADIENT_HEIGHT, zoomH, zoomV, FALSE, is_log_scale(), is_mono);
			}
		}
	}
	if (invocation == GHT_STRETCH && _payne_colourstretchmodel == COL_SAT && fit->naxes[2] == 3) {
		display_histo(com.sat_hist, cr, -1, width, height - GRADIENT_HEIGHT, zoomH, zoomV, FALSE, is_log_scale(), FALSE);
		if (siril_toggle_get_active(GTK_WIDGET(toggleOrig)))
			display_histo(hist_sat_backup, cr, -2, width, height - GRADIENT_HEIGHT, zoomH, zoomV, FALSE, is_log_scale(), FALSE);
	}
	display_scale(cr, width, height);
	return FALSE;
}

void on_histo_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	do_channel[0] = siril_toggle_get_active(GTK_WIDGET(toggles[0]));
	do_channel[1] = siril_toggle_get_active(GTK_WIDGET(toggles[1]));
	do_channel[2] = siril_toggle_get_active(GTK_WIDGET(toggles[2]));
	if (fit->naxes[2] == 3 && !(do_channel[0] && do_channel[1] && do_channel[2])) {
		if (_payne_colourstretchmodel == COL_HUMANLUM && invocation == GHT_STRETCH) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Human luminance colour model cannot be used: setting even weighted luminance colour model."));
			_payne_colourstretchmodel = COL_EVENLUM;
			gtk_drop_down_set_selected(GTK_DROP_DOWN(lookup_widget("combo_payne_colour_stretch_model")), COL_EVENLUM);
		}
		if (_payne_colourstretchmodel == COL_SAT && invocation == GHT_STRETCH) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
				_("Not all colour channels are selected. Saturation stretch cannot be used: setting independent channels colour model."));
			_payne_colourstretchmodel = COL_INDEP;
			gtk_drop_down_set_selected(GTK_DROP_DOWN(lookup_widget("combo_payne_colour_stretch_model")), COL_INDEP);
		}
	}

	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
}

void on_histogram_window_show(GtkWidget *object, gpointer user_data) {
	closing = FALSE;
	histo_startup();
	_initialize_clip_text();
	reset_cursors_and_values(TRUE);
	stretch_dialog_compute_histograms(fit);
}

gboolean on_button_histo_close_clicked(GtkButton *button, gpointer user_data) {
	closing = TRUE;
	set_cursor_waiting(TRUE);
	histo_close(TRUE, TRUE, TRUE);
	set_cursor_waiting(FALSE);
	siril_close_dialog("histogram_dialog");
	return FALSE;
}

void on_button_histo_reset_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	reset_cursors_and_values(FALSE);
	histo_close(TRUE, TRUE, FALSE);
	histo_startup();
	set_cursor_waiting(FALSE);
}

gboolean on_scale_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
	return FALSE;
}


void on_button_histo_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	set_controls_active(TRUE);

	// Check if there's anything to apply
	if (invocation == HISTO_STRETCH) {
		if ((_midtones == 0.5f) && (_shadows == 0.0f) && (_highlights == 1.0f)) {
			return;
		}
	} else if (invocation == GHT_STRETCH) {
		if ((_D == 0.0f) && (_BP == 0.0f)) {
			return;
		}
	}

	GtkCheckButton *seq_button = GTK_CHECK_BUTTON(lookup_widget("checkMTFSeq"));

	if (siril_toggle_get_active(GTK_WIDGET(seq_button)) && sequence_is_loaded()) {
		/* Apply to the whole sequence */
		sequence_working = TRUE;

		if (invocation == HISTO_STRETCH) {
			struct mtf_data *args = create_mtf_data();
			if (!args) {
				siril_log_color_message(_("Memory allocation failed.\n"), "red");
				return;
			}

			args->params.shadows = _shadows;
			args->params.midtones = _midtones;
			args->params.highlights = _highlights;
			args->params.do_red = do_channel[0];
			args->params.do_green = do_channel[1];
			args->params.do_blue = do_channel[2];
			args->seqEntry = strdup(gtk_editable_get_text(GTK_EDITABLE(lookup_widget("entryMTFSeq"))));
			if (args->seqEntry && args->seqEntry[0] == '\0') {
				free(args->seqEntry);
				args->seqEntry = strdup("stretch_");
			}
			args->seq = &com.seq;

			// Close the window for sequence processing
			histo_close(TRUE, FALSE, FALSE);
			siril_close_dialog("histogram_dialog");

			// Apply the process
			siril_toggle_set_active(GTK_WIDGET(seq_button), FALSE);
			apply_mtf_to_sequence(args);

		} else if (invocation == GHT_STRETCH) {
			struct ght_data *args = create_ght_data();
			if (!args) {
				siril_log_color_message(_("Memory allocation failed.\n"), "red");
				return;
			}

			args->params_ght = malloc(sizeof(struct ght_params));
			if (!args->params_ght) {
				destroy_ght_data(args);
				siril_log_color_message(_("Memory allocation failed.\n"), "red");
				return;
			}

			args->params_ght->B = _B;
			args->params_ght->D = _D;
			args->params_ght->LP = _LP;
			args->params_ght->SP = _SP;
			args->params_ght->HP = _HP;
			args->params_ght->BP = _BP;
			args->params_ght->stretchtype = _stretchtype;
			args->params_ght->payne_colourstretchmodel = _payne_colourstretchmodel;
			args->params_ght->do_red = do_channel[0];
			args->params_ght->do_green = do_channel[1];
			args->params_ght->do_blue = do_channel[2];
			args->params_ght->clip_mode = _clip_mode;

			const gchar* temp = gtk_editable_get_text(GTK_EDITABLE(lookup_widget("entryMTFSeq")));
			args->seqEntry = strdup(temp);
			if (args->seqEntry && args->seqEntry[0] == '\0') {
				free(args->seqEntry);
				args->seqEntry = strdup("stretch_");
			}
			args->seq = &com.seq;

			// Close the window for sequence processing
			histo_close(TRUE, FALSE, TRUE);
			siril_close_dialog("histogram_dialog");

			// Apply the process
			siril_toggle_set_active(GTK_WIDGET(seq_button), FALSE);
			apply_ght_to_sequence(args);
		}

	} else {
		/* Apply to single image */
		fit = gfit;

		// Check if preview is active
		gboolean preview_active = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));

		// If preview is active and no ROI, the image is already stretched in gfit
		// We just need to save undo and finalize
		if (preview_active && !gui.roi.active) {
			// The stretch is already applied to gfit via preview
			// Just need to prepare undo with the backup

			populate_roi();

			// Prepare undo state with ICC profile handling
			// Use the backup which contains the original unstretched image
			fits *backup = get_preview_gfit_backup();
			fits undo_fit = {0};
			memcpy(&undo_fit, backup, sizeof(fits));
			undo_fit.icc_profile = original_icc;
			undo_fit.color_managed = original_icc != NULL;

			gchar *log_string = NULL;

			if (invocation == HISTO_STRETCH) {
				struct mtf_data data = {
					.destroy_fn = NULL,
					.fit = NULL,
					.seq = NULL,
					.linked = TRUE,
					.params = {
						.midtones = _midtones,
						.shadows = _shadows,
						.highlights = _highlights,
						.do_red = do_channel[0],
						.do_green = do_channel[1],
						.do_blue = do_channel[2]
					},
					.uparams = {},
					.seqEntry = NULL,
					.auto_display_compensation = FALSE,
					.is_preview = FALSE
				};
				log_string = generate_mtf_log_message(&data, SUMMARY);

			} else if (invocation == GHT_STRETCH) {
				struct ght_params params = {
					.B = _B,
					.D = _D,
					.LP = _LP,
					.SP = _SP,
					.HP = _HP,
					.BP = _BP,
					.stretchtype = _stretchtype,
					.payne_colourstretchmodel = _payne_colourstretchmodel,
					.do_red = do_channel[0],
					.do_green = do_channel[1],
					.do_blue = do_channel[2],
					.clip_mode = _clip_mode
				};
				log_string = generate_ght_log_message(&params);
			}

			if (log_string) {
				undo_save_state(&undo_fit, log_string);
				siril_log_color_message("%s\n", "green", log_string);
				g_free(log_string);
			}

			// Clear backups and reinitialize
			clear_backup();
			clear_hist_backup();
			single_image_stretch_applied = TRUE;

			// Reinitialize for next use
			histo_startup();
			reset_cursors_and_values(FALSE);

			set_cursor("default");
			return;
		}

		// If preview wasn't active or ROI is active, need to apply the stretch
		if (!preview_active || gui.roi.active) {
			copy_backup_to_gfit();
			if (_payne_colourstretchmodel == COL_SAT) {
				clear_hsl();
				setup_hsl();
			}
		}

		populate_roi();

		// Prepare undo state with ICC profile handling BEFORE applying stretch
		fits *backup = get_preview_gfit_backup();
		fits undo_fit = {0};
		memcpy(&undo_fit, backup, sizeof(fits));
		undo_fit.icc_profile = original_icc;
		undo_fit.color_managed = original_icc != NULL;

		gchar *log_string = NULL;

		if (invocation == HISTO_STRETCH) {
			struct mtf_data data = {
					.destroy_fn = NULL,
					.fit = NULL,
					.seq = NULL,
					.linked = TRUE,
					.params = {
						.midtones = _midtones,
						.shadows = _shadows,
						.highlights = _highlights,
						.do_red = do_channel[0],
						.do_green = do_channel[1],
						.do_blue = do_channel[2]
					},
					.uparams = {},
					.seqEntry = NULL,
					.auto_display_compensation = FALSE,
					.is_preview = FALSE
			};
			log_string = generate_mtf_log_message(&data, SUMMARY);

		} else if (invocation == GHT_STRETCH) {
			struct ght_params params = {
				.B = _B,
				.D = _D,
				.LP = _LP,
				.SP = _SP,
				.HP = _HP,
				.BP = _BP,
				.stretchtype = _stretchtype,
				.payne_colourstretchmodel = _payne_colourstretchmodel,
				.do_red = do_channel[0],
				.do_green = do_channel[1],
				.do_blue = do_channel[2],
				.clip_mode = _clip_mode
			};
			log_string = generate_ght_log_message(&params);
		}

		if (log_string) {
			undo_save_state(&undo_fit, log_string);
			siril_log_color_message("%s\n", "green", log_string);
			g_free(log_string);
		}

		// Now apply using generic_img_worker
		struct generic_img_args *args = NULL;

		if (invocation == HISTO_STRETCH) {
			struct mtf_data *data = create_mtf_data();
			if (!data) {
				siril_log_color_message(_("Memory allocation failed.\n"), "red");
				return;
			}

			data->fit = gfit;
			data->params.shadows = _shadows;
			data->params.midtones = _midtones;
			data->params.highlights = _highlights;
			data->params.do_red = do_channel[0];
			data->params.do_green = do_channel[1];
			data->params.do_blue = do_channel[2];
			data->auto_display_compensation = FALSE;
			data->is_preview = FALSE;

			args = calloc(1, sizeof(struct generic_img_args));
			if (!args) {
				destroy_mtf_data(data);
				siril_log_color_message(_("Memory allocation failed.\n"), "red");
				return;
			}

			args->fit = gfit;
			args->mem_ratio = 1.0f;
			args->image_hook = mtf_single_image_hook;
			args->log_hook = mtf_log_hook;
			args->idle_function = mtf_single_image_idle;
			args->description = _("Histogram Transformation");
			args->verbose = TRUE;
			args->user = data;
			args->max_threads = com.max_thread;
			args->for_preview = FALSE;
			args->mask_aware = TRUE;
			args->for_roi = gui.roi.active;

		} else if (invocation == GHT_STRETCH) {
			struct ght_data *data = create_ght_data();
			if (!data) {
				siril_log_color_message(_("Memory allocation failed.\n"), "red");
				return;
			}

			data->fit = gfit;
			data->params_ght = malloc(sizeof(struct ght_params));
			if (!data->params_ght) {
				destroy_ght_data(data);
				siril_log_color_message(_("Memory allocation failed.\n"), "red");
				return;
			}

			data->params_ght->B = _B;
			data->params_ght->D = _D;
			data->params_ght->LP = _LP;
			data->params_ght->SP = _SP;
			data->params_ght->HP = _HP;
			data->params_ght->BP = _BP;
			data->params_ght->stretchtype = _stretchtype;
			data->params_ght->payne_colourstretchmodel = _payne_colourstretchmodel;
			data->params_ght->do_red = do_channel[0];
			data->params_ght->do_green = do_channel[1];
			data->params_ght->do_blue = do_channel[2];
			data->params_ght->clip_mode = _clip_mode;
			data->auto_display_compensation = FALSE;
			data->is_preview = FALSE;

			args = calloc(1, sizeof(struct generic_img_args));
			if (!args) {
				destroy_ght_data(data);
				siril_log_color_message(_("Memory allocation failed.\n"), "red");
				return;
			}

			args->fit = gfit;
			args->mem_ratio = (_payne_colourstretchmodel == COL_SAT) ? 2.0f : 1.0f;
			args->image_hook = ght_single_image_hook;
			args->log_hook = ght_log_hook;
			args->idle_function = ght_single_image_idle;
			args->description = _("Generalised Hyperbolic Stretch");
			args->verbose = TRUE;
			args->user = data;
			args->max_threads = com.max_thread;
			args->mask_aware = TRUE;
			args->for_preview = FALSE;
			args->for_roi = gui.roi.active;
		}

		if (args) {
			single_image_stretch_applied = TRUE;
			start_in_new_thread(generic_image_worker, args);

			// Note: cleanup and reinitialization now happens in idle function
		}

		set_cursor("default");
	}
}

void apply_histo_cancel() {
	set_cursor_waiting(TRUE);
	histo_close(TRUE, TRUE, TRUE);
	set_cursor_waiting(FALSE);
	siril_close_dialog("histogram_dialog");

}

void on_histoSpinZoom_value_changed(GtkRange *range, gpointer user_data) {
	adjust_histogram_vport_size();
	queue_window_redraw();
}

void on_histoToolAutoStretch_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);

	/* We always apply this function on original data */
	struct mtf_params params = {
		.shadows = _shadows,
		.midtones = _midtones,
		.highlights = _highlights,
		.do_red = TRUE,
		.do_green = TRUE,
		.do_blue = TRUE
	};

	if (!find_linked_midtones_balance_default(get_preview_gfit_backup(), &params)) {
		_shadows = params.shadows;
		_midtones = params.midtones;
		_highlights = 1.0f;
		_update_entry_text();
		update_histo_mtf();

		// Set auto display compensation flag and trigger recompute
		if (fit->color_managed && !profiles_identical(fit->icc_profile, com.gui_icc.monitor))
			auto_display_compensation = TRUE;

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = histo_update_preview;
		param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
		notify_update((gpointer) param);
		set_controls_active(FALSE);
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
			/* GTK4: gtk_window_resize removed */;
}

void setup_ght_dialog() {
			gtk_window_set_title(GTK_WINDOW(lookup_widget("histogram_dialog")),_("Generalised Hyperbolic Stretch Transformations"));
			// Hide the histogram controls
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("grid32")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("histoToolAutoStretch")), FALSE);
			// Make visible and configure the GHT controls
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("box_ghtcontrols")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("eyedropper_button")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("histo_clip_settings")), (gfit->naxes[2] == 3));
			GtkWidget *w;
			if (com.pref.gui.combo_theme == 0) {
				w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/eyedropper_dark.svg");
			} else {
				w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/eyedropper.svg");
			}
			gtk_button_set_child(GTK_BUTTON(lookup_widget("eyedropper_button")), w);
			gtk_widget_set_visible(w, TRUE);
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
			/* GTK4: gtk_window_resize removed */;
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
			/* GTK4: gtk_window_resize removed */;
	}
}

void toggle_histogram_window_visibility(int _invocation) {
	siril_close_preview_dialogs();
	init_toggles();
	invocation = _invocation;
	for (int i=0;i<3;i++) {
		do_channel[i] = TRUE;
	}

	if (gtk_widget_get_visible(lookup_widget("histogram_dialog"))) {
		set_cursor_waiting(TRUE);
		histo_close(TRUE, TRUE, TRUE);
		set_cursor_waiting(FALSE);

		siril_close_dialog("histogram_dialog");
	} else {
		fit = gfit;
		if (original_icc)
			cmsCloseProfile(original_icc);
		original_icc = copyICCProfile(gfit->icc_profile);
		icc_auto_assign_or_convert(gfit, ICC_ASSIGN_ON_STRETCH);
		single_image_stretch_applied = FALSE;
		// When opening the dialog with a single image loaded, we cache the original ICC
		// profile (may be NULL) in case the user closes the dialog without applying a
		// stretch, in which case we will revert.
		if (single_image_is_loaded()) {
			if (original_icc) {
				cmsCloseProfile(original_icc);
				original_icc = copyICCProfile(gfit->icc_profile);
			}
			icc_auto_assign_or_convert(gfit, ICC_ASSIGN_ON_STRETCH);
		} else {
			if (original_icc) {
				cmsCloseProfile(original_icc);
				original_icc = NULL;
			}
		}
		reset_cursors_and_values(TRUE);
		copy_gfit_to_backup();
		//Ensure the colour stretch model is initialised at startup
		_payne_colourstretchmodel = gtk_drop_down_get_selected(GTK_DROP_DOWN(lookup_widget("combo_payne_colour_stretch_model")));
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

/* GTK4 motion handler.  Replaces the GTK3 button_press_event /
 * motion_notify_event pair: GtkEventControllerMotion delivers (x, y) in
 * widget coordinates directly, so we no longer extract them from a
 * GdkEvent.  Behaviour is identical to the legacy version — when
 * `_click_on_histo` is set by the press handler we drag the corresponding
 * shadow / midtone / highlight slider as the cursor moves. */
static void histo_motion_cb(GtkEventControllerMotion *controller,
                            double x, double y, gpointer user_data) {
	(void)controller; (void)user_data;
	int width = get_width_of_histo();
	int height = get_height_of_histo();

	if (on_gradient(x, y, width, height)) {
		if (invocation == HISTO_STRETCH)
			set_cursor("grab");
	} else {
		set_cursor("default");
	}

	if (_click_on_histo && invocation == HISTO_STRETCH && active_sliders) {
		gdouble xpos = x / (gdouble) width;
		if (xpos < 0.0) xpos = 0.0;
		if (xpos > 1.0) xpos = 1.0;
		gchar *buffer = NULL;
		GtkEntry *histoMidEntry  = GTK_ENTRY(lookup_widget("histoMidEntry"));
		GtkEntry *histoShadEntry = GTK_ENTRY(lookup_widget("histoShadEntry"));
		GtkEntry *histoHighEntry = GTK_ENTRY(lookup_widget("histoHighEntry"));
		switch (_type_of_scale) {
		case SCALE_LOW:
			_shadows = ((float)xpos > _highlights) ? _highlights : (float)xpos;
			buffer = g_strdup_printf("%.7f", _shadows);
			gtk_editable_set_text(GTK_EDITABLE(histoShadEntry), buffer);
			break;
		case SCALE_MID:
			if (_highlights == _shadows)
				_midtones = _highlights;
			else
				_midtones = ((float)xpos - _shadows) / (_highlights - _shadows);
			if (_midtones > 1.f) _midtones = 1.f;
			if (_midtones < 0.f) _midtones = 0.f;
			buffer = g_strdup_printf("%.7f", _midtones);
			gtk_editable_set_text(GTK_EDITABLE(histoMidEntry), buffer);
			break;
		case SCALE_HI:
			if ((float)xpos < _shadows) _shadows = _highlights;
			else _highlights = (float)xpos;
			buffer = g_strdup_printf("%.7f", _highlights);
			gtk_editable_set_text(GTK_EDITABLE(histoHighEntry), buffer);
			break;
		}
		set_cursor("grabbing");
		update_histo_mtf();
		g_free(buffer);
	}
}


void on_drawingarea_histograms_leave_notify_event(GtkWidget *widget,
		GdkEvent *event, gpointer user_data) {
	set_cursor("default");
}

/* GTK4 press handler: classifies the click into a slider grab (HISTO_STRETCH)
 * or a symmetry-point set (GHT_STRETCH) and stashes the choice in module
 * statics for the motion / release handlers to act on. */
static void histo_press_cb(GtkGestureClick *gesture, int n_press,
                           double x, double y, gpointer user_data) {
	(void)gesture; (void)n_press; (void)user_data;
	int width  = get_width_of_histo();
	int height = get_height_of_histo();
	if (invocation == HISTO_STRETCH && active_sliders) {
		if (on_gradient(x, y, width, height)) {
			float delta = ((_highlights - _shadows) * _midtones) + _shadows;
			_click_on_histo = TRUE;
			gdouble xpos = x / (gdouble) width;
			if (fabsf((float)xpos - _highlights) < fabsf((float)xpos - _shadows)
			 && fabsf((float)xpos - _highlights) < fabsf((float)xpos - delta)) {
				_type_of_scale = SCALE_HI;
			} else if (fabsf((float)xpos - _shadows) < fabsf((float)xpos - delta)
			        && fabsf((float)xpos - _shadows) < fabsf((float)xpos - _highlights)) {
				_type_of_scale = SCALE_LOW;
			} else if (_shadows == _highlights && _shadows > 0.f) {
				_type_of_scale = SCALE_LOW;
			} else if (_shadows == _highlights && _shadows == 0.f) {
				_type_of_scale = SCALE_HI;
			} else {
				_type_of_scale = SCALE_MID;
			}
			set_cursor("grabbing");
		}
	} else if (invocation == GHT_STRETCH) {
		gdouble xpos = x / (gdouble) width;
		/* Reverse-transform the click position so the user-visible click
		 * lands on the right value in *unstretched* space (which is what
		 * the BP / SP spin buttons expect). */
		switch (_stretchtype) {
		case STRETCH_LINEAR:
			/* Linear mode has no symmetry point; clicking sets BP instead. */
			invxpos = xpos;
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
}


/* GTK4 release handler.  Commits the in-progress drag: pushes the staged
 * BP/SP into the spin buttons (GHT_STRETCH) or kicks the preview pipeline
 * (HISTO_STRETCH) so the histogram preview reflects the new sliders. */
static void histo_release_cb(GtkGestureClick *gesture, int n_press,
                             double x, double y, gpointer user_data) {
	(void)gesture; (void)n_press; (void)x; (void)y; (void)user_data;
	set_cursor("default");
	if (invocation == GHT_STRETCH && invxpos != -1.0) {
		if (_stretchtype == STRETCH_LINEAR)
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtBP")), invxpos);
		else
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ghtSP")), invxpos);
		invxpos = -1.0;
	}
	if (_click_on_histo) {
		_click_on_histo = FALSE;
		set_cursor_waiting(TRUE);
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = histo_update_preview;
		param->show_preview = siril_toggle_get_active(
		    GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
		notify_update((gpointer) param);
		set_cursor_waiting(FALSE);
	}
}


void on_spin_ghtD_value_changed(GtkSpinButton *button, gpointer user_data) {
	_D = (float) expm1(gtk_spin_button_get_value(button));
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
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
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
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
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
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
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
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
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
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
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
}

void on_spin_ghtBP_value_changed(GtkSpinButton *button, gpointer user_data) {
	_BP = (float) gtk_spin_button_get_value(button);
	if (_stretchtype == STRETCH_LINEAR) {
		update_histo_mtf();
		queue_window_redraw();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = histo_update_preview;
		param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
		notify_update((gpointer) param);
	}
}

void on_histo_clip_mode_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec;
	_clip_mode = (clip_mode_t) gtk_drop_down_get_selected(combo);
	update_histo_mtf();
	queue_window_redraw();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
}

void on_payneType_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec;
	_stretchtype = gtk_drop_down_get_selected(combo);
		if (_stretchtype == STRETCH_LINEAR)
			gtk_widget_set_tooltip_text(GTK_WIDGET(lookup_widget("drawingarea_histograms")), _("Clicking on the linear stretch histogram sets BP"));
		else
			gtk_widget_set_tooltip_text(GTK_WIDGET(lookup_widget("drawingarea_histograms")), _("Clicking on the histogram sets SP"));
	updateGHTcontrols();
	update_histo_mtf();
	queue_window_redraw();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
}

void on_payne_colour_stretch_model_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec;
	int tmp = gtk_drop_down_get_selected(combo);
	if (tmp == COL_SAT) {
		if (fit->naxes[2] != 3) {
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Channels warning"),
			_("Not all colour channels are selected. Saturation stretch cannot be used: setting independent channels colour model."));
			gtk_drop_down_set_selected(combo, COL_INDEP);
			_payne_colourstretchmodel = COL_INDEP;
		} else {
			set_cursor_waiting(TRUE);
			reset_cursors_and_values(FALSE);
			histo_close(TRUE, TRUE, FALSE);
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
			gtk_drop_down_set_selected(combo, COL_EVENLUM);
			_payne_colourstretchmodel = COL_EVENLUM;
			}
		}
		if (_payne_colourstretchmodel == COL_SAT) {
			set_cursor_waiting(TRUE);
			reset_cursors_and_values(FALSE);
			histo_close(TRUE, TRUE, FALSE);
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
			param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
			notify_update((gpointer) param);
			return;
		}
	}
}

gboolean on_histoShadEntry_focus_out_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("histoShadEntry"));
	float lo = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
	if (lo <= 0.f) lo = 0.f;
	_shadows = lo;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", lo);
	gtk_editable_set_text(GTK_EDITABLE(entry), str);
	g_free(str);
	set_cursor_waiting(FALSE);

	return FALSE;
}

void on_histoShadEntry_activate(GtkEntry *entry, gpointer user_data) {
	float lo = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
	if (lo <= 0.f) lo = 0.f;
	_shadows = lo;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", lo);
	gtk_editable_set_text(GTK_EDITABLE(entry), str);
	g_free(str);
	set_cursor_waiting(FALSE);
}

void on_histoMidEntry_activate(GtkEntry *entry, gpointer user_data) {
	float mid = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
	_midtones = mid;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", mid);
	gtk_editable_set_text(GTK_EDITABLE(entry), str);
	g_free(str);
	set_cursor_waiting(FALSE);
}

gboolean on_histoMidEntry_focus_out_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {

	GtkEntry *entry = GTK_ENTRY(lookup_widget("histoMidEntry"));
	float mid = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
	_midtones = mid;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", mid);
	gtk_editable_set_text(GTK_EDITABLE(entry), str);
	g_free(str);
	set_cursor_waiting(FALSE);

	return FALSE;
}

gboolean on_histoHighEntry_focus_out_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("histoHighEntry"));
	float hi = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
//	if (hi <= _shadows) hi = _shadows;
	if (hi >= 1.f) hi = 1.f;
	_highlights = hi;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", hi);
	gtk_editable_set_text(GTK_EDITABLE(entry), str);
	g_free(str);
	set_cursor_waiting(FALSE);

	return FALSE;
}

void on_histoHighEntry_activate(GtkEntry *entry, gpointer user_data) {
	float hi = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
//	if (hi <= _shadows) hi = _shadows;
	if (hi >= 1.f) hi = 1.f;
	_highlights = hi;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = histo_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))));
	notify_update((gpointer) param);
	gchar *str = g_strdup_printf("%8.7f", hi);
	gtk_editable_set_text(GTK_EDITABLE(entry), str);
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
	args->new_seq_prefix = strdup(mtf_args->seqEntry);
	args->load_new_sequence = TRUE;
	args->user = mtf_args;

	mtf_args->fit = NULL;	// not used here

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(mtf_args->seqEntry);
		free(mtf_args);
		free_generic_seq_args(args, TRUE);
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
	args->new_seq_prefix = strdup(ght_args->seqEntry);
	args->load_new_sequence = TRUE;
	args->user = ght_args;

	ght_args->fit = NULL;	// not used here

	if(!start_in_new_thread(generic_sequence_worker, args)) {
		free(ght_args->seqEntry);
		free(ght_args);
		free_generic_seq_args(args, TRUE);
	}
}

void on_histo_preview_toggled(GtkCheckButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!siril_toggle_get_active(GTK_WIDGET(GTK_CHECK_BUTTON(lookup_widget("HistoCheckPreview"))))) {
		/* if user click very fast */
		cancel_and_wait_for_preview();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = histo_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
