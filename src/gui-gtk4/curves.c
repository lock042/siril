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

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/icc_profile.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "algos/colors.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/siril_preview.h"
#include "curves.h"
#include "histogram.h"
#include "histogram_utils.h"
#include "filters/curve_transform.h"

// ---------------------------------------------------------------------------
// FORWARD DECLARATIONS
// ---------------------------------------------------------------------------
static int curves_process_with_worker(gboolean for_preview, gboolean for_roi);
static gboolean is_curves_log_scale();
static void update_range_ui_from_state();
static int curves_update_preview();
static void curves_update_image();
static void _initialize_clip_text();
static void curves_setup_hsl();
static void curves_clear_hsl();
static void curves_compute_special_histograms();
static gsl_histogram* curves_compute_histo_from_buffer(void* buf, fits* f);
static void init_all_curves();
static void init_all_curves_and_reset_channel();
static void clear_display_histogram();

void _update_entry_text();
gboolean curve_preview_idle(gpointer p);
gboolean curve_apply_idle(gpointer p);
gboolean redraw_curves(GtkWidget *widget, cairo_t *cr, gpointer data);
void draw_curve_points(cairo_t *cr, int width, int height);

/* GTK4 draw_func adapter + event-controller callbacks (wired in C, since
 * GTK4 dropped the GtkDrawingArea "draw" signal and the legacy
 * button-press / release / motion event signals). */
static void curves_draw_cb(GtkDrawingArea *area, cairo_t *cr,
                           int width, int height, gpointer data);
static void curves_da_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);
static void curves_da_released(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);
static void curves_da_drag_update(GtkGestureDrag *gesture,
		double offset_x, double offset_y, gpointer user_data);
static void curves_prev_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);
static void curves_next_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);

// ---------------------------------------------------------------------------
// UI CACHE
// ---------------------------------------------------------------------------

static GtkAdjustment *curves_adj_zoom = NULL;
static GtkDropDown *curves_interpolation_combo = NULL;
static GtkEntry *curves_id_entry = NULL, *curves_x_entry = NULL, *curves_y_entry = NULL, *curves_clip_low = NULL, *curves_clip_high = NULL, *curves_seq_entry = NULL;
static GtkGrid *curves_point_grid = NULL;
static GtkCheckButton *curves_sequence_check = NULL, *curves_preview_check = NULL, *curves_log_check = NULL;
static GtkToggleButton *curves_grid_toggle = NULL;
static GtkWidget *curves_drawingarea = NULL, *curves_viewport = NULL, *curves_dialog = NULL;
static GtkToggleButton *curves_chan_rgb_radio = NULL;
static GtkToggleButton *curves_chan_r_radio = NULL;
static GtkToggleButton *curves_chan_g_radio = NULL;
static GtkToggleButton *curves_chan_b_radio = NULL;
static GtkToggleButton *curves_chan_l_radio = NULL;
static GtkToggleButton *curves_chan_c_radio = NULL;
static GtkToggleButton *curves_chan_s_radio = NULL;
static GtkCheckButton *curves_range_check = NULL;
static GtkScale *curves_range_min_scale = NULL;
static GtkScale *curves_range_max_scale = NULL;
static GtkScale *curves_range_feather_scale = NULL;
static GtkToggleButton *curves_pipette_button = NULL;
static GtkButton *curves_undo_stage_button = NULL;
static GtkButton *curves_apply_stage_button = NULL;

// ---------------------------------------------------------------------------
// STATE
// ---------------------------------------------------------------------------

static curve_channel_config gui_channels[CHAN_COUNT];

static enum curve_channel current_channel = CHAN_RGB_K;
enum curve_algorithm algorithm = AKIMA_SPLINE;

static gboolean closing = FALSE;
static fits *fit = NULL;
static long clipped[] = {0, 0};
static gsl_histogram *display_histogram[MAXVPORT] = {NULL, NULL, NULL, NULL};
static gsl_histogram *display_histogram_lum = NULL;
static gsl_histogram *display_histogram_sat = NULL;
static gsl_histogram *display_histogram_chroma = NULL;
// Input (pre-curve) histogram backups for mathematical live update during drag
static gsl_histogram *hist_input[MAXVPORT] = {NULL, NULL, NULL, NULL};
static gsl_histogram *hist_input_lum = NULL;
static gsl_histogram *hist_input_sat = NULL;
static gsl_histogram *hist_input_chroma = NULL;
static gboolean is_drawingarea_pressed = FALSE;
static void *curves_huebuf = NULL, *curves_satbuf = NULL, *curves_lumbuf = NULL;
static cmsHPROFILE original_icc = NULL;
static gboolean single_image_stretch_applied = FALSE;
static point *selected_point = NULL;
static double point_size = 5.0;

// Stage stack (non-destructive iterative workflow)
static GList *stage_stack = NULL;
static fits *original_fit_copy = NULL;

// Pipette
static float pipette_value = -1.0f;

int compare_points(const void *a, const void *b) {
	point *point_a = (point *) a;
	point *point_b = (point *) b;
	if (point_a->x < point_b->x) return -1;
	if (point_a->x > point_b->x) return 1;
	return 0;
}

static point* copy_point(const point *src, gpointer user_data) {
	point *new_p = g_new(point, 1);
	new_p->x = src->x;
	new_p->y = src->y;
	return new_p;
}

// ---------------------------------------------------------------------------
// MATHEMATICAL HISTOGRAM UPDATE (like apply_mtf_to_histo in histogram.c)
// ---------------------------------------------------------------------------

static void apply_curve_to_histo(gsl_histogram *src, gsl_histogram *dst, enum curve_channel chan) {
	if (!src || !dst) return;
	size_t n = gsl_histogram_bins(src);
	gsl_histogram *result = gsl_histogram_alloc(n);
	gsl_histogram_set_ranges_uniform(result, 0.0, (double)(n - 1));

	GList *points = gui_channels[chan].points;
	int np = points ? (int)g_list_length(points) : 0;

	if (np < 2) {
		gsl_histogram_memcpy(dst, src);
		gsl_histogram_free(result);
		return;
	}

	akima_spline_data akima;
	cubic_spline_data cspline;
	double linear_slopes[MAX_POINTS];

	gboolean use_akima = (algorithm == AKIMA_SPLINE && np >= 5);
	if (algorithm == LINEAR) {
		linear_fit(points, linear_slopes);
	} else if (use_akima) {
		akima_spline_fit(points, &akima);
	} else {
		cubic_spline_fit(points, &cspline);
	}

	float inv_n = 1.0f / (float)(n - 1);
	for (size_t i = 0; i < n; i++) {
		double count = gsl_histogram_get(src, i);
		if (count == 0.0) continue;
		float x = (float)i * inv_n;
		float y;
		if (algorithm == LINEAR)
			y = linear_interpolate(x, points, linear_slopes);
		else if (use_akima)
			y = akima_spline_interpolate(x, &akima);
		else
			y = cubic_spline_interpolate(x, &cspline);
		y = fmaxf(0.0f, fminf(1.0f, y));
		size_t dest = (size_t)(y * (float)(n - 1) + 0.5f);
		if (dest >= n) dest = n - 1;
		gsl_histogram_accumulate(result, (double)dest, count);
	}
	gsl_histogram_memcpy(dst, result);
	gsl_histogram_free(result);
}

// Apply RGB/K then a per-channel curve to one histogram slot (accounts for
// any RGB/K edit made before the user switched to an individual channel).
static void apply_rgb_channel_pipeline(int c, enum curve_channel per_chan) {
	if (c >= (int)fit->naxes[2] || !hist_input[c] || !display_histogram[c]) return;
	apply_curve_to_histo(hist_input[c], display_histogram[c], CHAN_RGB_K);
	gsl_histogram *tmp = gsl_histogram_clone(display_histogram[c]);
	apply_curve_to_histo(tmp, display_histogram[c], per_chan);
	gsl_histogram_free(tmp);
}

static void update_display_histogram_from_curve() {
	if (!fit) return;
	int nc = (int)fit->naxes[2];
	switch (current_channel) {
	case CHAN_RGB_K:
		for (int c = 0; c < nc; c++)
			apply_curve_to_histo(hist_input[c], display_histogram[c], CHAN_RGB_K);
		break;
	case CHAN_R: apply_rgb_channel_pipeline(0, CHAN_R); break;
	case CHAN_G: apply_rgb_channel_pipeline(1, CHAN_G); break;
	case CHAN_B: apply_rgb_channel_pipeline(2, CHAN_B); break;
	case CHAN_L:
		apply_curve_to_histo(hist_input_lum, display_histogram_lum, CHAN_L);
		break;
	case CHAN_C:
		apply_curve_to_histo(hist_input_chroma, display_histogram_chroma, CHAN_C);
		break;
	case CHAN_S:
		apply_curve_to_histo(hist_input_sat, display_histogram_sat, CHAN_S);
		break;
	default:
		break;
	}
}

static void refresh_hist_input_from_fit() {
	if (!fit) return;
	for (int i = 0; i < (int)fit->naxes[2]; i++) {
		if (hist_input[i]) { gsl_histogram_free(hist_input[i]); hist_input[i] = NULL; }
		if (com.layers_hist[i]) hist_input[i] = gsl_histogram_clone(com.layers_hist[i]);
	}
	if (hist_input_lum) { gsl_histogram_free(hist_input_lum); hist_input_lum = NULL; }
	if (display_histogram_lum) hist_input_lum = gsl_histogram_clone(display_histogram_lum);
	if (hist_input_sat) { gsl_histogram_free(hist_input_sat); hist_input_sat = NULL; }
	if (display_histogram_sat) hist_input_sat = gsl_histogram_clone(display_histogram_sat);
	if (hist_input_chroma) { gsl_histogram_free(hist_input_chroma); hist_input_chroma = NULL; }
	if (display_histogram_chroma) hist_input_chroma = gsl_histogram_clone(display_histogram_chroma);
}

// ---------------------------------------------------------------------------
// STAGE STACK HELPERS
// ---------------------------------------------------------------------------

static void free_snapshot(gpointer data) {
	curve_state_snapshot *snap = (curve_state_snapshot *) data;
	for (int i = 0; i < CHAN_COUNT; i++) {
		if (snap->channels[i].points)
			g_list_free_full(snap->channels[i].points, g_free);
	}
	g_free(snap);
}

static void free_stage_stack() {
	if (stage_stack) {
		g_list_free_full(stage_stack, free_snapshot);
		stage_stack = NULL;
	}
}

static curve_state_snapshot *snapshot_current_state() {
	curve_state_snapshot *snap = g_new(curve_state_snapshot, 1);
	for (int i = 0; i < CHAN_COUNT; i++) {
		snap->channels[i] = gui_channels[i];
		snap->channels[i].points = gui_channels[i].points
			? g_list_copy_deep(gui_channels[i].points, (GCopyFunc) copy_point, NULL)
			: NULL;
	}
	snap->algorithm = algorithm;
	return snap;
}

static void apply_snapshot_to_fit(fits *f, curve_state_snapshot *snap) {
	struct curve_params *p = new_curve_params();
	for (int i = 0; i < CHAN_COUNT; i++) {
		p->channels[i].active       = snap->channels[i].active;
		p->channels[i].range_enabled= snap->channels[i].range_enabled;
		p->channels[i].lum_min      = snap->channels[i].lum_min;
		p->channels[i].lum_max      = snap->channels[i].lum_max;
		p->channels[i].feather      = snap->channels[i].feather;
		p->channels[i].points       = snap->channels[i].points
			? g_list_copy_deep(snap->channels[i].points, (GCopyFunc) copy_point, NULL)
			: NULL;
	}
	p->algorithm = snap->algorithm;
	p->fit = f;
	p->verbose = FALSE;
	apply_curve(f, f, p, TRUE);
	free_curve_params(p);
}

static void update_stage_buttons() {
	if (curves_undo_stage_button)
		gtk_widget_set_sensitive(GTK_WIDGET(curves_undo_stage_button), stage_stack != NULL);
}

static void curves_rebuild_from_stages() {
	if (!original_fit_copy) return;

	set_cursor_waiting(TRUE);
	copyfits(original_fit_copy, gfit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
	for (GList *l = stage_stack; l; l = l->next)
		apply_snapshot_to_fit(gfit, (curve_state_snapshot *) l->data);

	copy_gfit_to_backup();
	compute_histo_for_fit(fit);
	clear_display_histogram();
	for (int i = 0; i < (int)fit->naxes[2]; i++)
		display_histogram[i] = gsl_histogram_clone(com.layers_hist[i]);
	if (fit->naxes[2] == 3) {
		if (curves_huebuf) { free(curves_huebuf); curves_huebuf = NULL; }
		if (curves_satbuf) { free(curves_satbuf); curves_satbuf = NULL; }
		if (curves_lumbuf) { free(curves_lumbuf); curves_lumbuf = NULL; }
		curves_compute_special_histograms();
	}
	refresh_hist_input_from_fit();
	init_all_curves();
	selected_point = (point *) gui_channels[current_channel].points->data;
	update_range_ui_from_state();
	_update_entry_text();
	update_stage_buttons();
	gtk_widget_queue_draw(curves_drawingarea);
	notify_gfit_data_modified();
	gfit_modified_update_gui();
	curves_update_image();
	set_cursor_waiting(FALSE);
}

// ---------------------------------------------------------------------------
// CONFIG
// ---------------------------------------------------------------------------

static void init_single_channel_config(int ch) {
	if (gui_channels[ch].points) {
		g_list_free_full(gui_channels[ch].points, g_free);
		gui_channels[ch].points = NULL;
	}
	point *p = g_new(point, 1);
	p->x = 0.0; p->y = 0.0;
	gui_channels[ch].points = g_list_insert_sorted(gui_channels[ch].points, p, (GCompareFunc) compare_points);
	p = g_new(point, 1);
	p->x = 1.0; p->y = 1.0;
	gui_channels[ch].points = g_list_insert_sorted(gui_channels[ch].points, p, (GCompareFunc) compare_points);

	gui_channels[ch].active = FALSE;
	gui_channels[ch].range_enabled = FALSE;
	gui_channels[ch].lum_min = 0.0f;
	gui_channels[ch].lum_max = 1.0f;
	gui_channels[ch].feather = 0.25f;
}

static void init_all_curves() {
	for (int i = 0; i < CHAN_COUNT; i++)
		init_single_channel_config(i);
	selected_point = (point *) gui_channels[current_channel].points->data;
}

static void init_all_curves_and_reset_channel() {
	for (int i = 0; i < CHAN_COUNT; i++)
		init_single_channel_config(i);
	current_channel = CHAN_RGB_K;
	selected_point = (point *) gui_channels[CHAN_RGB_K].points->data;
}

static void update_range_ui_from_state() {
	if (!curves_range_check) return;
	g_signal_handlers_block_by_func(curves_range_check, on_curve_check_range_button_toggled, NULL);
	g_signal_handlers_block_by_func(curves_range_min_scale, on_curves_range_value_changed, NULL);
	g_signal_handlers_block_by_func(curves_range_max_scale, on_curves_range_value_changed, NULL);
	g_signal_handlers_block_by_func(curves_range_feather_scale, on_curves_feather_value_changed, NULL);

	curve_channel_config *cfg = &gui_channels[current_channel];
	siril_toggle_set_active(GTK_WIDGET(curves_range_check), cfg->range_enabled);
	gtk_range_set_value(GTK_RANGE(curves_range_min_scale), cfg->lum_min * 100.0);
	gtk_range_set_value(GTK_RANGE(curves_range_max_scale), cfg->lum_max * 100.0);
	gtk_range_set_value(GTK_RANGE(curves_range_feather_scale), cfg->feather * 100.0);

	g_signal_handlers_unblock_by_func(curves_range_check, on_curve_check_range_button_toggled, NULL);
	g_signal_handlers_unblock_by_func(curves_range_min_scale, on_curves_range_value_changed, NULL);
	g_signal_handlers_unblock_by_func(curves_range_max_scale, on_curves_range_value_changed, NULL);
	g_signal_handlers_unblock_by_func(curves_range_feather_scale, on_curves_feather_value_changed, NULL);
}

static void switch_channel_view(enum curve_channel new_chan) {
	current_channel = new_chan;
	if (gui_channels[new_chan].points)
		selected_point = (point *) gui_channels[new_chan].points->data;
	else {
		init_single_channel_config(new_chan);
		selected_point = (point *) gui_channels[new_chan].points->data;
	}
	update_range_ui_from_state();

	GtkToggleButton *chan_radios[CHAN_COUNT] = {
		curves_chan_rgb_radio,
		curves_chan_r_radio,
		curves_chan_g_radio,
		curves_chan_b_radio,
		curves_chan_l_radio,
		curves_chan_c_radio,
		curves_chan_s_radio,
	};
	if (new_chan >= 0 && new_chan < CHAN_COUNT && chan_radios[new_chan]) {
		GtkToggleButton *btn = chan_radios[new_chan];
		g_signal_handlers_block_matched(btn, G_SIGNAL_MATCH_DATA, 0, 0, NULL, NULL, NULL);
		gtk_toggle_button_set_active(btn, TRUE);
		g_signal_handlers_unblock_matched(btn, G_SIGNAL_MATCH_DATA, 0, 0, NULL, NULL, NULL);
	}
}

// ---------------------------------------------------------------------------
// UI SETUP
// ---------------------------------------------------------------------------

void curves_dialog_init_statics() {
	if (curves_adj_zoom == NULL) {
		curves_adj_zoom = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "curves_adj_zoom"));
		curves_interpolation_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "curves_interpolation_combo"));
		curves_dialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_dialog"));
		curves_drawingarea = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_drawingarea"));

		/* GTK4: GtkDrawingArea has no "draw" signal — wire the draw
		 * function explicitly, plus the GtkGestureClick / GtkGestureDrag
		 * controllers that replace the GTK3 button-press / release /
		 * motion-notify signals. */
		if (curves_drawingarea) {
			gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(curves_drawingarea),
			                               curves_draw_cb, NULL, NULL);

			GtkGesture *click = gtk_gesture_click_new();
			gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click), 0);
			g_signal_connect(click, "pressed",  G_CALLBACK(curves_da_pressed),  NULL);
			g_signal_connect(click, "released", G_CALLBACK(curves_da_released), NULL);
			gtk_widget_add_controller(curves_drawingarea, GTK_EVENT_CONTROLLER(click));

			GtkGesture *drag = gtk_gesture_drag_new();
			gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag), GDK_BUTTON_PRIMARY);
			g_signal_connect(drag, "drag-update",
			                 G_CALLBACK(curves_da_drag_update), NULL);
			gtk_widget_add_controller(curves_drawingarea, GTK_EVENT_CONTROLLER(drag));

			/* Group click + drag so they share sequence state — without
			 * this whichever gesture claims first DENIES the other for
			 * the rest of the sequence and drag-update never fires on
			 * macOS. Same pattern as histo_display.c / image_interactions.c. */
			gtk_gesture_group(click, drag);
		}

		/* GTK3 wired prev/next via "button-press-event" on the wrapper
		 * GtkBox of each arrow GtkImage; GTK4 uses GtkGestureClick on the
		 * same wrapper boxes (the .ui keeps the ids). */
		GtkWidget *prev_box = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_prev_event_box"));
		GtkWidget *next_box = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_next_event_box"));
		if (prev_box) {
			GtkGesture *gprev = gtk_gesture_click_new();
			gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(gprev), GDK_BUTTON_PRIMARY);
			g_signal_connect(gprev, "pressed", G_CALLBACK(curves_prev_pressed), NULL);
			gtk_widget_add_controller(prev_box, GTK_EVENT_CONTROLLER(gprev));
		}
		if (next_box) {
			GtkGesture *gnext = gtk_gesture_click_new();
			gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(gnext), GDK_BUTTON_PRIMARY);
			g_signal_connect(gnext, "pressed", G_CALLBACK(curves_next_pressed), NULL);
			gtk_widget_add_controller(next_box, GTK_EVENT_CONTROLLER(gnext));
		}

		curves_id_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_id_entry"));
		curves_x_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_x_entry"));
		curves_y_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_y_entry"));
		curves_clip_low = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_clip_low"));
		curves_clip_high = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_clip_high"));
		curves_seq_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_seq_entry"));
		curves_point_grid = GTK_GRID(gtk_builder_get_object(gui.builder, "curves_point_grid"));
		curves_sequence_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "curves_sequence_check"));
		curves_preview_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "curves_preview_check"));
		curves_log_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "curves_log_check"));
		curves_grid_toggle = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_grid_toggle"));
		curves_viewport = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_viewport"));

		curves_range_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "curve_check_range_button"));
		curves_range_min_scale = GTK_SCALE(gtk_builder_get_object(gui.builder, "curve_slider_min"));
		curves_range_max_scale = GTK_SCALE(gtk_builder_get_object(gui.builder, "curve_slider_max"));
		curves_range_feather_scale = GTK_SCALE(gtk_builder_get_object(gui.builder, "curve_slider_feather"));
		curves_chan_rgb_radio = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_chan_rgb"));
		curves_chan_r_radio = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_chan_r"));
		curves_chan_g_radio = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_chan_g"));
		curves_chan_b_radio = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_chan_b"));
		curves_chan_l_radio = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_chan_l"));
		curves_chan_c_radio = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_chan_c"));
		curves_chan_s_radio = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_chan_s"));
		curves_pipette_button = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_pipette_button"));
		curves_undo_stage_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "curves_undo_stage_button"));
		curves_apply_stage_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "curves_apply_stage_button"));
	}
}

static void draw_background_gradient(cairo_t *cr, int width, int height) {
	cairo_pattern_t *pat = cairo_pattern_create_linear(0, height, 0, 0);

	if (current_channel == CHAN_S || current_channel == CHAN_C) {
		cairo_pattern_add_color_stop_rgb(pat, 0.0, 0.1, 0.1, 0.1);

		if (current_channel == CHAN_S) {
			cairo_pattern_add_color_stop_rgb(pat, 1.0, 0.16, 0.12, 0.27);
		} else {
			cairo_pattern_add_color_stop_rgb(pat, 1.0, 0.27, 0.16, 0.12);
		}

		cairo_set_source(cr, pat);
		cairo_rectangle(cr, 0, 0, width, height);
		cairo_fill(cr);

		cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
		cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_set_font_size(cr, 10);

		const char *lbl = (current_channel == CHAN_S) ? "High Sat" : "High Chroma";
		cairo_move_to(cr, 5, 12);
		cairo_show_text(cr, lbl);

		lbl = (current_channel == CHAN_S) ? "Low Sat" : "Low Chroma";
		cairo_move_to(cr, 5, height - 5);
		cairo_show_text(cr, lbl);

	} else {

		double r=0.12, g=0.12, b=0.12;

		if (current_channel == CHAN_R)      { r=0.15; g=0.08; b=0.08; }
		else if (current_channel == CHAN_G) { r=0.08; g=0.15; b=0.08; }
		else if (current_channel == CHAN_B) { r=0.08; g=0.08; b=0.18; }
		else if (current_channel == CHAN_L) { r=0.14; g=0.14; b=0.14; }

		cairo_pattern_add_color_stop_rgb(pat, 0.0, 0.05, 0.05, 0.05);
		cairo_pattern_add_color_stop_rgb(pat, 1.0, r, g, b);

		cairo_set_source(cr, pat);
		cairo_rectangle(cr, 0, 0, width, height);
		cairo_fill(cr);
	}

	cairo_pattern_destroy(pat);

	if (siril_toggle_get_active(GTK_WIDGET(curves_grid_toggle))) {
		cairo_set_source_rgba(cr, 0.4, 0.4, 0.4, 0.3);
		cairo_set_line_width(cr, 1.0);
		double dashes[] = {2.0};
		cairo_set_dash(cr, dashes, 1, 0);

		for(int i=1; i<4; i++) {
			cairo_move_to(cr, i*width/4.0, 0); cairo_line_to(cr, i*width/4.0, height);
			cairo_move_to(cr, 0, i*height/4.0); cairo_line_to(cr, width, i*height/4.0);
		}
		cairo_stroke(cr);
		cairo_set_dash(cr, NULL, 0, 0);

		cairo_set_source_rgba(cr, 0.3, 0.3, 0.3, 0.5);
		cairo_move_to(cr, 0, height); cairo_line_to(cr, width, 0);
		cairo_stroke(cr);
	}
}

static void render_spline(cairo_t *cr, GList *points, int width, int height) {
	if (!points) return;

	if (algorithm == LINEAR) {
		double slopes[MAX_POINTS];
		linear_fit(points, slopes);
		for (double x = 0; x <= 1; x += 1.0 / width) {
			double y = linear_interpolate(x, points, slopes);
			if (x == 0) cairo_move_to(cr, x * width, height - (y * height));
			else cairo_line_to(cr, x * width, height - (y * height));
		}
	} else if (algorithm == CUBIC_SPLINE) {
		cubic_spline_data cspline_data;
		cubic_spline_fit(points, &cspline_data);
		for (double x = 0; x <= 1; x += 1.0 / width) {
			double y = cubic_spline_interpolate(x, &cspline_data);
			if (x == 0) cairo_move_to(cr, x * width, height - (y * height));
			else cairo_line_to(cr, x * width, height - (y * height));
		}
	} else if (algorithm == AKIMA_SPLINE) {
		akima_spline_data akima_data;
		akima_spline_fit(points, &akima_data);
		for (double x = 0; x <= 1; x += 1.0 / width) {
			double y = akima_spline_interpolate(x, &akima_data);
			if (x == 0) cairo_move_to(cr, x * width, height - (y * height));
			else cairo_line_to(cr, x * width, height - (y * height));
		}
	}
}

static void draw_ghost_curves(cairo_t *cr, int width, int height) {
	cairo_set_line_width(cr, 1.0);

	for(int i=0; i<CHAN_COUNT; i++) {
		if (i == current_channel) continue;

		gboolean active = FALSE;
		if (gui_channels[i].points && g_list_length(gui_channels[i].points) > 2) active = TRUE;
		else if (gui_channels[i].points) {
			point *p1 = (point*)gui_channels[i].points->data;
			point *p2 = (point*)gui_channels[i].points->next->data;
			if (fabs(p1->y) > 0.01 || fabs(p2->y - 1.0) > 0.01) active = TRUE;
		}

		if (active) {
			switch(i) {
				case CHAN_R: cairo_set_source_rgba(cr, 1.0, 0.2, 0.2, 0.3); break;
				case CHAN_G: cairo_set_source_rgba(cr, 0.2, 1.0, 0.2, 0.3); break;
				case CHAN_B: cairo_set_source_rgba(cr, 0.2, 0.2, 1.0, 0.3); break;
				case CHAN_L: cairo_set_source_rgba(cr, 0.9, 0.9, 0.9, 0.3); break;
				default:     cairo_set_source_rgba(cr, 0.6, 0.6, 0.6, 0.2); break;
			}
			render_spline(cr, gui_channels[i].points, width, height);
			cairo_stroke(cr);
		}
	}
}

static void draw_active_curve(cairo_t *cr, int width, int height) {
	double r = 1.0, g = 1.0, b = 1.0;

	if (current_channel == CHAN_R) { r=1.0; g=0.3; b=0.3; }
	else if (current_channel == CHAN_G) { r=0.3; g=1.0; b=0.3; }
	else if (current_channel == CHAN_B) { r=0.4; g=0.4; b=1.0; }
	else if (current_channel == CHAN_L) { r=0.9; g=0.9; b=0.9; } // L = White/Grey
	else if (current_channel == CHAN_C) { r=1.0; g=0.7; b=0.2; } // C = Gold/Orange
	else if (current_channel == CHAN_S) { r=0.9; g=0.2; b=0.9; } // S = Magenta/Purple
	else { r=1.0; g=1.0; b=1.0; } // RGB/K = White

	cairo_set_source_rgb(cr, r, g, b);
	cairo_set_line_width(cr, 2.0);

	render_spline(cr, gui_channels[current_channel].points, width, height);

	cairo_stroke(cr);
}

void draw_curve_points(cairo_t *cr, int width, int height) {
	GList *iter;
	point *p;

	cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
	cairo_set_line_width(cr, 1.5);

	GList *points = gui_channels[current_channel].points;

	for (iter = points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		double x = p->x * width;
		double y = height - (p->y * height);

		cairo_set_source_rgba(cr, 0.1, 0.1, 0.1, 0.8);
		cairo_arc(cr, x, y, point_size, 0, 2*M_PI);
		cairo_fill(cr);
		cairo_set_source_rgba(cr, 0.8, 0.8, 0.8, 0.8);
		cairo_arc(cr, x, y, point_size, 0, 2*M_PI);
		cairo_stroke(cr);
	}
	if (selected_point) {
		double x = selected_point->x * width;
		double y = height - (selected_point->y * height);

		cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 0.4);
		cairo_arc(cr, x, y, point_size+3, 0, 2*M_PI);
		cairo_fill(cr);

		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
		cairo_arc(cr, x, y, point_size+1, 0, 2*M_PI);
		cairo_fill(cr);
	}
}

void draw_range_overlays(cairo_t *cr, int width, int height) {
	if (gui_channels[current_channel].range_enabled) {
		float mn = gui_channels[current_channel].lum_min;
		float mx = gui_channels[current_channel].lum_max;

		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 0.6);

		if (mn > 0.0) {
			cairo_rectangle(cr, 0, 0, mn * width, height);
			cairo_fill(cr);
		}
		if (mx < 1.0) {
			cairo_rectangle(cr, mx * width, 0, (1.0 - mx) * width, height);
			cairo_fill(cr);
		}

		cairo_set_source_rgba(cr, 1.0, 0.6, 0.0, 0.8);
		cairo_set_line_width(cr, 1.5);
		double dashes[] = {4.0};
		cairo_set_dash(cr, dashes, 1, 0);

		if (mn > 0.0) {
			cairo_move_to(cr, mn * width, 0); cairo_line_to(cr, mn * width, height); cairo_stroke(cr);
		}
		if (mx < 1.0) {
			cairo_move_to(cr, mx * width, 0); cairo_line_to(cr, mx * width, height); cairo_stroke(cr);
		}
		cairo_set_dash(cr, NULL, 0, 0);
	}
}

// ---------------------------------------------------------------------------
// HSL HISTOGRAM FUNCTIONS
// ---------------------------------------------------------------------------

// Draw a histogram with a custom RGB color, matching the style of display_histo()
static void draw_histo_with_color(gsl_histogram *histo, cairo_t *cr, int width, int height,
                                   double r, double g, double b, gboolean is_log) {
	if (!histo) return;
	if (width <= 0) return;

	size_t norm = gsl_histogram_bins(histo) - 1;
	float vals_per_px = (float)norm / (float)width;
	size_t i, nb_orig_bins = gsl_histogram_bins(histo);

	// Allocate array to store binned values for display
	static gfloat *displayed_values = NULL;
	static int nb_bins_allocated = 0;

	if (nb_bins_allocated != width) {
		gfloat *tmp;
		nb_bins_allocated = width;
		tmp = realloc(displayed_values, nb_bins_allocated * sizeof(gfloat));
		if (!tmp) {
			if (displayed_values != NULL) {
				free(displayed_values);
				displayed_values = NULL;
			}
			return;
		}
		displayed_values = tmp;
		memset(displayed_values, 0, nb_bins_allocated * sizeof(gfloat));
	}

	// Set color and line style
	cairo_set_source_rgb(cr, r, g, b);
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1.5);

	// First loop: build the bins and find the maximum
	i = 0;
	double graph_height = 0.0;
	int current_bin = 0;

	do {
		double bin_val = 0.0;
		while (i < nb_orig_bins && (double)i / vals_per_px <= current_bin + 0.5) {
			bin_val += gsl_histogram_get(histo, i);
			i++;
		}
		if (is_log && bin_val != 0.0) {
			bin_val = logf(bin_val);
		}
		displayed_values[current_bin] = bin_val;
		if (bin_val > graph_height)
			graph_height = bin_val;
		current_bin++;
	} while (i < nb_orig_bins && current_bin < nb_bins_allocated);

	if (!graph_height)
		return;

	// Second loop: draw the histogram as a continuous line
	for (i = 0; i < nb_bins_allocated; i++) {
		double bin_height = height - height * displayed_values[i] / graph_height;
		cairo_line_to(cr, i, bin_height);
	}
	cairo_stroke(cr);
}

// Convert RGB image to HSL and store in buffers
static void curves_fit_to_hsl() {
	if (!fit || fit->naxes[2] != 3) return;

	size_t npixels = fit->rx * fit->ry;
	if (fit->type == DATA_FLOAT) {
		float* hf = (float*) curves_huebuf;
		float* sf = (float*) curves_satbuf;
		float* lf = (float*) curves_lumbuf;
		// rgb_to_hslf has internal branches and an early return, so the loop
		// body cannot be SIMD-vectorized; use plain thread parallelism.
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0 ; i < npixels ; i++) {
			rgb_to_hslf(fit->fpdata[0][i], fit->fpdata[1][i], fit->fpdata[2][i], &hf[i], &sf[i], &lf[i]);
		}
	} else {
		WORD *hw = (WORD*) curves_huebuf;
		WORD* sw = (WORD*) curves_satbuf;
		WORD* lw = (WORD*) curves_lumbuf;
		// rgbw_to_hslw wraps rgb_to_hslf (branches + early return), so the loop
		// body cannot be SIMD-vectorized; use plain thread parallelism.
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0 ; i < npixels ; i++) {
			rgbw_to_hslw(fit->pdata[0][i], fit->pdata[1][i], fit->pdata[2][i], &hw[i], &sw[i], &lw[i]);
		}
	}
}

// Compute histogram from a buffer (luminance, saturation, or chroma)
static gsl_histogram* curves_compute_histo_from_buffer(void* buf, fits* f) {
	if (!buf || !f) return NULL;

	size_t i, ndata, size;
	size = get_histo_size(f);
	gsl_histogram *histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, 1.0 + 1.0 / size);
	ndata = f->naxes[0] * f->naxes[1];

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
	{
		gsl_histogram *histo_thr = gsl_histogram_alloc(size + 1);
		gsl_histogram_set_ranges_uniform(histo_thr, 0, 1.0 + 1.0 / size);
		if (f->type == DATA_FLOAT) {
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

// Compute chroma histogram (chroma is computed from saturation and luminance)
static gsl_histogram* curves_compute_histo_chroma(fits* f) {
	if (!curves_satbuf || !curves_lumbuf || !f || f->naxes[2] != 3) return NULL;

	size_t ndata = f->naxes[0] * f->naxes[1];

	// Allocate temporary buffer for chroma values
	void* chroma_buf = NULL;
	if (f->type == DATA_FLOAT) {
		chroma_buf = malloc(ndata * sizeof(float));
		if (!chroma_buf) return NULL;

		float *cf = (float*) chroma_buf;
		float *sf = (float*) curves_satbuf;
		float *lf = (float*) curves_lumbuf;

		// Chroma approximation: saturation weighted by (1 - |2*L - 1|)
		// This gives higher chroma for saturated colors at mid-luminance
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < ndata; i++) {
			float l_factor = 1.0f - fabsf(2.0f * lf[i] - 1.0f);
			cf[i] = sf[i] * l_factor;
		}
	} else {
		chroma_buf = malloc(ndata * sizeof(WORD));
		if (!chroma_buf) return NULL;

		WORD *cw = (WORD*) chroma_buf;
		WORD *sw = (WORD*) curves_satbuf;
		WORD *lw = (WORD*) curves_lumbuf;
		double invnorm = 1.0 / USHRT_MAX_DOUBLE;

#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < ndata; i++) {
			double l_norm = lw[i] * invnorm;
			double s_norm = sw[i] * invnorm;
			double l_factor = 1.0 - fabs(2.0 * l_norm - 1.0);
			cw[i] = (WORD)(s_norm * l_factor * USHRT_MAX_DOUBLE);
		}
	}

	gsl_histogram* histo = curves_compute_histo_from_buffer(chroma_buf, f);
	free(chroma_buf);
	return histo;
}

static void curves_setup_hsl() {
	if (!fit || fit->naxes[2] != 3) return;

	// Free existing buffers
	if (curves_huebuf) free(curves_huebuf);
	if (curves_satbuf) free(curves_satbuf);
	if (curves_lumbuf) free(curves_lumbuf);

	// Allocate new buffers
	if (fit->type == DATA_FLOAT) {
		curves_huebuf = malloc(fit->rx * fit->ry * sizeof(float));
		curves_satbuf = malloc(fit->rx * fit->ry * sizeof(float));
		curves_lumbuf = malloc(fit->rx * fit->ry * sizeof(float));
	} else {
		curves_huebuf = malloc(fit->rx * fit->ry * sizeof(WORD));
		curves_satbuf = malloc(fit->rx * fit->ry * sizeof(WORD));
		curves_lumbuf = malloc(fit->rx * fit->ry * sizeof(WORD));
	}

	// Convert RGB to HSL
	curves_fit_to_hsl();
}

static void curves_clear_hsl() {
	if (curves_huebuf) {
		free(curves_huebuf);
		curves_huebuf = NULL;
	}
	if (curves_satbuf) {
		free(curves_satbuf);
		curves_satbuf = NULL;
	}
	if (curves_lumbuf) {
		free(curves_lumbuf);
		curves_lumbuf = NULL;
	}
	if (display_histogram_lum) {
		gsl_histogram_free(display_histogram_lum);
		display_histogram_lum = NULL;
	}
	if (display_histogram_sat) {
		gsl_histogram_free(display_histogram_sat);
		display_histogram_sat = NULL;
	}
	if (display_histogram_chroma) {
		gsl_histogram_free(display_histogram_chroma);
		display_histogram_chroma = NULL;
	}
}

static void curves_compute_special_histograms() {
	if (!fit || fit->naxes[2] != 3) return;

	// Free existing special histograms
	if (display_histogram_lum) {
		gsl_histogram_free(display_histogram_lum);
		display_histogram_lum = NULL;
	}
	if (display_histogram_sat) {
		gsl_histogram_free(display_histogram_sat);
		display_histogram_sat = NULL;
	}
	if (display_histogram_chroma) {
		gsl_histogram_free(display_histogram_chroma);
		display_histogram_chroma = NULL;
	}

	// Compute HSL if not already done
	if (!curves_lumbuf) {
		curves_setup_hsl();
	}

	// Compute histograms
	if (curves_lumbuf)
		display_histogram_lum = curves_compute_histo_from_buffer(curves_lumbuf, fit);
	if (curves_satbuf)
		display_histogram_sat = curves_compute_histo_from_buffer(curves_satbuf, fit);
	display_histogram_chroma = curves_compute_histo_chroma(fit);
}

// ---------------------------------------------------------------------------
// LOGIC
// ---------------------------------------------------------------------------

static gboolean is_curves_log_scale() {
	return (siril_toggle_get_active(GTK_WIDGET(curves_log_check)));
}

/* GTK4 draw_func adapter — forwards to redraw_curves. */
static void curves_draw_cb(GtkDrawingArea *area, cairo_t *cr,
                           int width, int height, gpointer data) {
	(void) width; (void) height;
	redraw_curves(GTK_WIDGET(area), cr, data);
}

gboolean redraw_curves(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int width = gtk_widget_get_width(curves_drawingarea);
	int height = gtk_widget_get_height(curves_drawingarea);
	if (height == 1) return FALSE;

	draw_background_gradient(cr, width, height);

	// Display appropriate histogram based on current channel
	if (current_channel == CHAN_L && display_histogram_lum) {
		// Display luminance histogram in white/gray
		draw_histo_with_color(display_histogram_lum, cr, width, height, 0.9, 0.9, 0.9, is_curves_log_scale());
	} else if (current_channel == CHAN_S && display_histogram_sat) {
		// Display saturation histogram in magenta/purple
		draw_histo_with_color(display_histogram_sat, cr, width, height, 0.9, 0.2, 0.9, is_curves_log_scale());
	} else if (current_channel == CHAN_C && display_histogram_chroma) {
		// Display chroma histogram in orange/gold
		draw_histo_with_color(display_histogram_chroma, cr, width, height, 1.0, 0.7, 0.2, is_curves_log_scale());
	} else {
		// Display RGB histograms
		gboolean is_mono = fit ? (fit->naxes[2] == 1) : FALSE;
		for (int i = 0; i < MAXVPORT; i++) {
			if (current_channel == CHAN_R && i != 0) continue;
			if (current_channel == CHAN_G && i != 1) continue;
			if (current_channel == CHAN_B && i != 2) continue;
			display_histo(display_histogram[i], cr, i, width, height, 1.0, 1.0, FALSE, is_curves_log_scale(), is_mono);
		}
	}

	draw_range_overlays(cr, width, height);
	draw_ghost_curves(cr, width, height);
	draw_active_curve(cr, width, height);
	draw_curve_points(cr, width, height);

	// Pipette marker: vertical dashed line at sampled channel value
	if (pipette_value >= 0.0f && pipette_value <= 1.0f) {
		double px = pipette_value * width;
		double dashes[] = {4.0, 3.0};
		cairo_set_dash(cr, dashes, 2, 0);
		cairo_set_source_rgba(cr, 1.0, 0.9, 0.1, 0.85);
		cairo_set_line_width(cr, 1.5);
		cairo_move_to(cr, px, 0);
		cairo_line_to(cr, px, height);
		cairo_stroke(cr);
		cairo_set_dash(cr, NULL, 0, 0);
		cairo_arc(cr, px, height * 0.5, 4.5, 0, 2 * G_PI);
		cairo_fill(cr);
	}

	return FALSE;
}

static int curves_update_preview() {
	fit = gui.roi.active ? &gui.roi.fit : gfit;
	if (!closing) {
		copy_backup_to_gfit();
		curves_process_with_worker(TRUE, gui.roi.active);
	}
	return 0;
}

static void curves_update_image() {
	set_cursor_waiting(TRUE);
	if (processing_is_job_active())
		stop_processing_thread();
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = curves_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(curves_preview_check));
	notify_update((gpointer) param);
	set_cursor_waiting(FALSE);
}

static void _initialize_clip_text() {
	gtk_editable_set_text(GTK_EDITABLE(curves_clip_low), "0.000%");
	gtk_editable_set_text(GTK_EDITABLE(curves_clip_high), "0.000%");
}

// Helper to update UI from global stats
static void update_clip_ui_from_stats() {
	if (!fit) return;
	size_t total_pixels = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	if (total_pixels == 0) total_pixels = 1; // Safety

	double tmp;
	char buffer[32];

	// High
	tmp = (double) clipped[1] * 100.0 / (double) total_pixels;
	if (tmp < 0) tmp = 0;
	snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_editable_set_text(GTK_EDITABLE(curves_clip_high), buffer);

	// Low
	tmp = (double) clipped[0] * 100.0 / (double) total_pixels;
	if (tmp < 0) tmp = 0;
	snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_editable_set_text(GTK_EDITABLE(curves_clip_low), buffer);
}

static void reset_cursors_and_values(gboolean full_reset) {
	_initialize_clip_text();
	gtk_adjustment_set_value(curves_adj_zoom, 1.0);
	gtk_editable_set_text(GTK_EDITABLE(curves_seq_entry), "curve_");
	siril_toggle_set_active(GTK_WIDGET(curves_preview_check), TRUE);
	siril_toggle_set_active(GTK_WIDGET(curves_grid_toggle), TRUE);

	if (full_reset) {
		gtk_drop_down_set_selected(GTK_DROP_DOWN(curves_interpolation_combo), AKIMA_SPLINE);
		algorithm = AKIMA_SPLINE;
		if (curves_chan_rgb_radio) {
			g_signal_handlers_block_matched(curves_chan_rgb_radio, G_SIGNAL_MATCH_DATA, 0, 0, NULL, NULL, NULL);
			gtk_toggle_button_set_active(curves_chan_rgb_radio, TRUE);
			g_signal_handlers_unblock_matched(curves_chan_rgb_radio, G_SIGNAL_MATCH_DATA, 0, 0, NULL, NULL, NULL);
			current_channel = CHAN_RGB_K;
		}
		init_all_curves();
	} else {
		init_single_channel_config(current_channel);
		selected_point = (point*) gui_channels[current_channel].points->data;
	}
	update_range_ui_from_state();
	_update_entry_text();
	update_gfit_curves_histogram_if_needed();
}

static void clear_display_histogram() {
	if (display_histogram[0]) {
		for (int i = 0; i < fit->naxes[2]; i++) {
			gsl_histogram_free(display_histogram[i]);
			display_histogram[i] = NULL;
		}
	}
	if (display_histogram_lum) {
		gsl_histogram_free(display_histogram_lum);
		display_histogram_lum = NULL;
	}
	if (display_histogram_sat) {
		gsl_histogram_free(display_histogram_sat);
		display_histogram_sat = NULL;
	}
	if (display_histogram_chroma) {
		gsl_histogram_free(display_histogram_chroma);
		display_histogram_chroma = NULL;
	}
}

static void curves_startup() {
	init_all_curves();
	add_roi_callback(curves_histogram_change_between_roi_and_image);
	roi_supported(TRUE);
	copy_gfit_to_backup();

	// Save a copy of the original image for stage-stack undo (only at first open)
	if (!original_fit_copy) {
		original_fit_copy = calloc(1, sizeof(fits));
		copyfits(gfit, original_fit_copy, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
		free_stage_stack();
	}

	compute_histo_for_fit(fit);
	for (int i = 0; i < (int)fit->naxes[2]; i++) {
		if (display_histogram[i]) { gsl_histogram_free(display_histogram[i]); display_histogram[i] = NULL; }
		display_histogram[i] = gsl_histogram_clone(com.layers_hist[i]);
	}

	// Compute special histograms for L, C, S channels
	if (fit->naxes[2] == 3) {
		curves_compute_special_histograms();
	}

	// Snapshot histograms as "input" reference for mathematical live-drag update
	refresh_hist_input_from_fit();
}

static void curves_close(gboolean update_image_if_needed, gboolean revert_icc_profile) {
	// Disable pipette mode if active
	if (curves_pipette_button && siril_toggle_get_active(GTK_WIDGET(curves_pipette_button))) {
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		siril_toggle_set_active(GTK_WIDGET(curves_pipette_button), FALSE);
	}
	pipette_value = -1.0f;

	// Clean up stage stack
	if (original_fit_copy) {
		clearfits(original_fit_copy);
		free(original_fit_copy);
		original_fit_copy = NULL;
	}
	free_stage_stack();

	for (int i = 0; i < (int)fit->naxes[2]; i++) {
		set_histogram(display_histogram[i], i);
		display_histogram[i] = NULL;
		if (hist_input[i]) { gsl_histogram_free(hist_input[i]); hist_input[i] = NULL; }
	}
	if (hist_input_lum) { gsl_histogram_free(hist_input_lum); hist_input_lum = NULL; }
	if (hist_input_sat) { gsl_histogram_free(hist_input_sat); hist_input_sat = NULL; }
	if (hist_input_chroma) { gsl_histogram_free(hist_input_chroma); hist_input_chroma = NULL; }
	if (is_preview_active() && !copy_backup_to_gfit() && update_image_if_needed) {
		set_cursor_waiting(TRUE);
		notify_gfit_data_modified();
		gfit_modified_update_gui();
	}
	if (revert_icc_profile && !single_image_stretch_applied) {
		if (gfit->icc_profile) cmsCloseProfile(gfit->icc_profile);
		gfit->icc_profile = copyICCProfile(original_icc);
		color_manage(gfit, gfit->icc_profile != NULL);
	}
	clear_backup();
	clear_display_histogram();
	curves_clear_hsl();
	roi_supported(FALSE);
	remove_roi_callback(curves_histogram_change_between_roi_and_image);
}

// Build a fully-populated curve_params from the GUI channel state. Caller owns
// the returned params and must free with free_curve_params().
static struct curve_params *build_curve_params_from_gui() {
	struct curve_params *params = new_curve_params();
	if (!params) return NULL;
	for (int i = 0; i < CHAN_COUNT; i++) {
		params->channels[i].active = gui_channels[i].active;
		params->channels[i].range_enabled = gui_channels[i].range_enabled;
		params->channels[i].lum_min = gui_channels[i].lum_min;
		params->channels[i].lum_max = gui_channels[i].lum_max;
		params->channels[i].feather = gui_channels[i].feather;
		if (gui_channels[i].points)
			params->channels[i].points = g_list_copy_deep(gui_channels[i].points, (GCopyFunc) copy_point, NULL);
		else
			params->channels[i].points = NULL;
	}
	params->algorithm = algorithm;
	params->target_channel = current_channel;
	return params;
}

static int curves_process_with_worker(gboolean for_preview, gboolean for_roi) {
	struct curve_params *params = build_curve_params_from_gui();
	if (!params) { PRINT_ALLOC_ERR; return 1; }

	// Reset global clipping and link it to the worker params
	clipped[0] = 0;
	clipped[1] = 0;
	params->clipped_count = clipped;

	params->fit = for_roi ? &gui.roi.fit : gfit;
	params->verbose = !for_preview;
	params->for_preview = for_preview;

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_curve_params(params);
		return 1;
	}

	// Set the fit based on whether ROI is active
	args->fit = params->fit;
	args->mem_ratio = 2.0f; // Curves need memory for depth conversions
	args->image_hook = curve_image_hook;
	args->idle_function = for_preview ? curve_preview_idle : curve_apply_idle;
	args->description = _("Curve Transformation");
	args->verbose = !for_preview;
	args->user = params;
	args->log_hook = curves_log_hook;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;
	args->for_roi = for_roi;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return 1;
	}
	return 0;
}

// Idle function implementation to update UI stats
gboolean curve_preview_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	// Update UI from stats gathered by worker
	update_clip_ui_from_stats();

	stop_processing_thread();
	if (args->retval == 0) gfit_modified_update_gui();
	free_generic_img_args(args);
	return FALSE;
}

gboolean curve_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	// Update UI from stats gathered by worker
	update_clip_ui_from_stats();

	stop_processing_thread();
	if (args->retval == 0) {
		gfit_modified_update_gui();

		copy_gfit_to_backup();

		// Full commit: reset stage stack and update the original reference
		free_stage_stack();
		if (original_fit_copy) { clearfits(original_fit_copy); free(original_fit_copy); original_fit_copy = NULL; }
		original_fit_copy = calloc(1, sizeof(fits));
		copyfits(gfit, original_fit_copy, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);

		compute_histo_for_fit(fit);
		for (int i = 0; i < (int)fit->naxes[2]; i++) {
			if (display_histogram[i]) { gsl_histogram_free(display_histogram[i]); display_histogram[i] = NULL; }
			display_histogram[i] = gsl_histogram_clone(com.layers_hist[i]);
		}

		// Clear HSL buffers to force recomputation from the new transformed image
		if (fit->naxes[2] == 3) {
			if (curves_huebuf) { free(curves_huebuf); curves_huebuf = NULL; }
			if (curves_satbuf) { free(curves_satbuf); curves_satbuf = NULL; }
			if (curves_lumbuf) { free(curves_lumbuf); curves_lumbuf = NULL; }
			curves_compute_special_histograms();
		}

		// Refresh input histogram backup for live drag updates
		refresh_hist_input_from_fit();

		// Reset all curves to identity - transformation is now baked into the image
		init_all_curves();
		selected_point = (point*) gui_channels[current_channel].points->data;
		update_range_ui_from_state();
		_update_entry_text();
		update_stage_buttons();

		gtk_widget_queue_draw(curves_drawingarea);
	}
	free_generic_img_args(args);
	return FALSE;
}

void _update_entry_text() {
	if (selected_point == NULL) return;
	gchar *buffer;
	buffer = g_strdup_printf("%.7f", selected_point->x);
	gtk_editable_set_text(GTK_EDITABLE(curves_x_entry), buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", selected_point->y);
	gtk_editable_set_text(GTK_EDITABLE(curves_y_entry), buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%d", g_list_index(gui_channels[current_channel].points, selected_point));
	gtk_editable_set_text(GTK_EDITABLE(curves_id_entry), buffer);
	g_free(buffer);
}

void update_gfit_curves_histogram_if_needed() {
	invalidate_gfit_histogram();
	if (gtk_widget_get_visible(curves_dialog)) {
		compute_histo_for_fit(fit);
		// Clear HSL buffers to force recomputation from current image
		if (fit->naxes[2] == 3) {
			if (curves_huebuf) { free(curves_huebuf); curves_huebuf = NULL; }
			if (curves_satbuf) { free(curves_satbuf); curves_satbuf = NULL; }
			if (curves_lumbuf) { free(curves_lumbuf); curves_lumbuf = NULL; }
			curves_compute_special_histograms();
		}
		gtk_widget_queue_draw(curves_drawingarea);
	}
}

void curves_reset_after_undo() {
	if (!gtk_widget_get_visible(curves_dialog))
		return;

	fit = gfit;
	copy_gfit_to_backup();

	compute_histo_for_fit(fit);

	clear_display_histogram();
	for (int i = 0; i < fit->naxes[2]; i++)
		display_histogram[i] = gsl_histogram_clone(com.layers_hist[i]);

	// Clear HSL buffers to force recomputation from current image
	if (fit->naxes[2] == 3) {
		if (curves_huebuf) { free(curves_huebuf); curves_huebuf = NULL; }
		if (curves_satbuf) { free(curves_satbuf); curves_satbuf = NULL; }
		if (curves_lumbuf) { free(curves_lumbuf); curves_lumbuf = NULL; }
		curves_compute_special_histograms();
	}

	init_all_curves_and_reset_channel();
	switch_channel_view(CHAN_RGB_K);
	selected_point = (point*) gui_channels[CHAN_RGB_K].points->data;
	update_range_ui_from_state();
	_update_entry_text();

	gtk_widget_queue_draw(curves_drawingarea);

	single_image_stretch_applied = FALSE;
}

void curves_histogram_change_between_roi_and_image() {
	fit = gui.roi.active ? &gui.roi.fit : gfit;
	gui.roi.operation_supports_roi = TRUE;
	curves_update_image();
}

// ---------------------------------------------------------------------------
// EVENTS
// ---------------------------------------------------------------------------

static void on_curves_channel_radio_toggled(GtkToggleButton *button, gpointer user_data) {
	if (!siril_toggle_get_active(GTK_WIDGET(button)))
		return;
	enum curve_channel new_chan = GPOINTER_TO_INT(user_data);
	if (new_chan != current_channel) {
		switch_channel_view(new_chan);
		pipette_value = -1.0f;
		update_display_histogram_from_curve();
		gtk_widget_queue_draw(curves_drawingarea);
		curves_update_image();
		_update_entry_text();
	}
}

void on_curves_chan_rgb_toggled(GtkToggleButton *button, gpointer user_data) {
	on_curves_channel_radio_toggled(button, GINT_TO_POINTER(CHAN_RGB_K));
}
void on_curves_chan_r_toggled(GtkToggleButton *button, gpointer user_data) {
	on_curves_channel_radio_toggled(button, GINT_TO_POINTER(CHAN_R));
}
void on_curves_chan_g_toggled(GtkToggleButton *button, gpointer user_data) {
	on_curves_channel_radio_toggled(button, GINT_TO_POINTER(CHAN_G));
}
void on_curves_chan_b_toggled(GtkToggleButton *button, gpointer user_data) {
	on_curves_channel_radio_toggled(button, GINT_TO_POINTER(CHAN_B));
}
void on_curves_chan_l_toggled(GtkToggleButton *button, gpointer user_data) {
	on_curves_channel_radio_toggled(button, GINT_TO_POINTER(CHAN_L));
}
void on_curves_chan_c_toggled(GtkToggleButton *button, gpointer user_data) {
	on_curves_channel_radio_toggled(button, GINT_TO_POINTER(CHAN_C));
}
void on_curves_chan_s_toggled(GtkToggleButton *button, gpointer user_data) {
	on_curves_channel_radio_toggled(button, GINT_TO_POINTER(CHAN_S));
}

void on_curve_check_range_button_toggled(GtkToggleButton *button, gpointer user_data) {
	gboolean active = siril_toggle_get_active(GTK_WIDGET(button));
	gui_channels[current_channel].range_enabled = active;

	GtkWidget *expander = lookup_widget("curve_range_expander");
	if (expander) {
		gtk_widget_set_sensitive(expander, active);
		gtk_expander_set_expanded((GtkExpander*)expander, active);
	}

	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image();
}

void on_curves_range_value_changed(GtkRange *range, gpointer user_data) {
	double min_v = gtk_range_get_value(GTK_RANGE(curves_range_min_scale));
	double max_v = gtk_range_get_value(GTK_RANGE(curves_range_max_scale));

	if (min_v >= max_v) {
		if (GTK_WIDGET(range) == GTK_WIDGET(curves_range_min_scale)) {
			max_v = min_v + 1.0; if (max_v > 100) max_v = 100;
			gtk_range_set_value(GTK_RANGE(curves_range_max_scale), max_v);
		} else {
			min_v = max_v - 1.0; if (min_v < 0) min_v = 0;
			gtk_range_set_value(GTK_RANGE(curves_range_min_scale), min_v);
		}
	}

	gui_channels[current_channel].lum_min = min_v / 100.0f;
	gui_channels[current_channel].lum_max = max_v / 100.0f;

	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image();
}

void on_curves_feather_value_changed(GtkRange *range, gpointer user_data) {
	double val = gtk_range_get_value(range);
	gui_channels[current_channel].feather = val / 100.0f;
	curves_update_image();
}

void on_curves_display_toggle(GtkToggleButton *togglebutton, gpointer user_data) {
	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image();
}

void on_curves_window_show(GtkWidget *object, gpointer user_data) {
	fit = gfit;
	closing = FALSE;
	pipette_value = -1.0f;
	curves_startup();
	_initialize_clip_text();
	reset_cursors_and_values(TRUE);
	compute_histo_for_fit(fit);
	update_stage_buttons();
}

void on_curves_close_button_clicked(GtkButton *button, gpointer user_data) {
	closing = TRUE;
	set_cursor_waiting(TRUE);
	curves_close(TRUE, TRUE);
	set_cursor_waiting(FALSE);
	siril_close_dialog("curves_dialog");
}

void on_curves_reset_button_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	init_single_channel_config(current_channel);
	selected_point = (point*) gui_channels[current_channel].points->data;
	update_range_ui_from_state();
	_initialize_clip_text();
	gtk_adjustment_set_value(curves_adj_zoom, 1.0);
	_update_entry_text();
	curves_update_image();
	gtk_widget_queue_draw(curves_drawingarea);
	init_all_curves();
	set_cursor_waiting(FALSE);
}

void on_curves_apply_button_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa()) return;

	if (siril_toggle_get_active(GTK_WIDGET(curves_sequence_check)) && sequence_is_loaded()) {
		struct curve_data *args = calloc(1, sizeof(struct curve_data));
		for(int i=0; i<CHAN_COUNT; i++) {
			args->params.channels[i].active = gui_channels[i].active;
			args->params.channels[i].range_enabled = gui_channels[i].range_enabled;
			args->params.channels[i].lum_min = gui_channels[i].lum_min;
			args->params.channels[i].lum_max = gui_channels[i].lum_max;
			args->params.channels[i].feather = gui_channels[i].feather;
			if (gui_channels[i].points)
				args->params.channels[i].points = g_list_copy_deep(gui_channels[i].points, (GCopyFunc)copy_point, NULL);
		}
		args->params.algorithm = algorithm;
		args->params.target_channel = current_channel;

		args->seq_entry = strdup(gtk_editable_get_text(GTK_EDITABLE(curves_seq_entry)));
		args->seq = &com.seq;
		if (args->seq_entry && args->seq_entry[0] == '\0') {
			free(args->seq_entry);
			args->seq_entry = strdup("curve_");
		}
		curves_close(FALSE, TRUE);
		siril_close_dialog("curves_dialog");
		siril_toggle_set_active(GTK_WIDGET(curves_sequence_check), FALSE);
		apply_curve_to_sequence(args);
	} else {
		fit = gfit;
		gboolean preview_active = siril_toggle_get_active(GTK_WIDGET(curves_preview_check));

		if (preview_active && !gui.roi.active) {
			// The curve is already applied to gfit via preview; we are
			// about to discard the backup that holds the pre-curve
			// pixels, so push them onto the undo stack first.  Without
			// this, Undo after Apply (preview-on path) is a no-op
			// because the worker that would normally call
			// undo_save_state never runs — the curve was applied
			// incrementally by the preview pipeline.
			struct curve_params *undo_params = build_curve_params_from_gui();
			gchar *summary = curves_log_hook(undo_params, SUMMARY);
			fits undo_fit = { 0 };
			memcpy(&undo_fit, get_preview_gfit_backup(), sizeof(fits));
			undo_fit.icc_profile = original_icc;
			undo_fit.color_managed = original_icc != NULL;
			undo_save_state(&undo_fit, "%s", summary);
			g_free(summary);
			free_curve_params(undo_params);

			populate_roi();
			clear_backup();
			clear_display_histogram();
			single_image_stretch_applied = TRUE;
			// Full commit: reset stage stack and let curves_startup() save the new original
			free_stage_stack();
			if (original_fit_copy) { clearfits(original_fit_copy); free(original_fit_copy); original_fit_copy = NULL; }
			curves_startup();
			reset_cursors_and_values(FALSE);
			update_stage_buttons();
			set_cursor("default");
			return;
		}

		// Preview not active or ROI active: apply the curve now, defer reinit to curve_apply_idle
		copy_backup_to_gfit();
		populate_roi();
		curves_process_with_worker(FALSE, gui.roi.active);
		single_image_stretch_applied = TRUE;
		set_cursor("default");
	}
}

int curve_finalize_hook(struct generic_seq_args *args) {
	struct curve_data *curve_data = (struct curve_data *) args->user;
	int return_value = seq_finalize_hook(args);
	for(int i=0; i<CHAN_COUNT; i++) {
		if (curve_data->params.channels[i].points)
			g_list_free_full(curve_data->params.channels[i].points, g_free);
	}
	free(curve_data);
	return return_value;
}

static int curve_image_hook_seq(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	struct curve_data *curve_args = (struct curve_data *) args->user;
	apply_curve(fit, fit, &curve_args->params, FALSE);
	return 0;
}

void apply_curve_to_sequence(struct curve_data *curve_args) {
	struct generic_seq_args *args = create_default_seqargs(curve_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = curve_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = curve_finalize_hook;
	args->image_hook = curve_image_hook_seq;
	args->stop_on_error = FALSE;
	args->description = _("Curves Transform");
	args->has_output = TRUE;
	args->new_seq_prefix = strdup(curve_args->seq_entry);
	args->load_new_sequence = TRUE;
	args->user = curve_args;
	curve_args->fit = NULL;

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(curve_args->seq_entry);
		free(curve_args);
		free_generic_seq_args(args, TRUE);
	}
}

void apply_curves_cancel() {
	set_cursor_waiting(TRUE);
	curves_close(TRUE, TRUE);
	set_cursor_waiting(FALSE);
	siril_close_dialog("curves_dialog");
}

void on_curves_reset_zoom_clicked(GtkButton *button, gpointer user_data) {
	gtk_adjustment_set_value(curves_adj_zoom, 1.0);
}

void on_curves_spin_zoom_value_changed(GtkRange *range, gpointer user_data) {
	int target_width, current_width;
	double zoom = gtk_adjustment_get_value(curves_adj_zoom);
	current_width = gtk_widget_get_width(curves_viewport);
	target_width = (int) (((double) current_width) * zoom);
	gtk_widget_set_size_request(curves_drawingarea, target_width, 281);
	gtk_widget_queue_draw(curves_drawingarea);
}

void on_curves_interpolation_combo_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *widget = GTK_DROP_DOWN(obj);
	(void)pspec;
	algorithm = gtk_drop_down_get_selected(widget);
	curves_update_image();
	gtk_widget_queue_draw(curves_drawingarea);
}

void setup_curve_dialog() {
	gtk_window_set_title(GTK_WINDOW(curves_dialog), _("Curve Transformation"));
	gtk_widget_set_visible(GTK_WIDGET(curves_point_grid), TRUE);
	/* GTK4: gtk_window_resize removed */
}

void toggle_curves_window_visibility() {
	if (gtk_widget_get_visible(GTK_WIDGET(gtk_builder_get_object(gui.builder, "histogram_dialog"))))
	siril_close_dialog("histogram_dialog");

	curves_dialog_init_statics();

	if (gtk_widget_get_visible(curves_dialog)) {
		set_cursor_waiting(TRUE);
		curves_close(TRUE, TRUE);
		set_cursor_waiting(FALSE);
		siril_close_dialog("curves_dialog");
	} else {
		single_image_stretch_applied = FALSE;
		if (single_image_is_loaded()) {
			if (original_icc) {
				cmsCloseProfile(original_icc);
				original_icc = copyICCProfile(gfit->icc_profile);
			}
			icc_auto_assign_or_convert(gfit, ICC_ASSIGN_ON_STRETCH);
		} else {
			if (original_icc) { cmsCloseProfile(original_icc); original_icc = NULL; }
		}
		reset_cursors_and_values(TRUE);
		copy_gfit_to_backup();
		setup_curve_dialog();
		if (gui.rendering_mode == LINEAR_DISPLAY) setup_stretch_sliders();
		siril_open_dialog("curves_dialog");
	}
}

/* GtkGestureDrag drag-update — fires for each motion event while a button is
 * held (including on macOS, where GtkEventControllerMotion does not).  Drags
 * the selected curve point in normalised [0,1] coordinates. */
static void curves_da_drag_update(GtkGestureDrag *gesture,
		double offset_x, double offset_y, gpointer user_data) {
	(void)user_data;
	if (!is_drawingarea_pressed || !selected_point) return;
	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
	double sx, sy;
	gtk_gesture_drag_get_start_point(gesture, &sx, &sy);
	double x = sx + offset_x;
	double y = sy + offset_y;
	gdouble xpos = x / (gdouble) gtk_widget_get_width(curves_drawingarea);
	/* y position is inverted (top-left origin -> bottom-up curve). */
	gdouble ypos = 1 - y / (gdouble) gtk_widget_get_height(curves_drawingarea);
	xpos = (xpos < 0) ? 0 : (xpos > 1) ? 1 : xpos;
	selected_point->y = (ypos < 0) ? 0 : (ypos > 1) ? 1 : ypos;
	/* Constrain x movement to the gap between the neighbouring points. */
	GList *points = gui_channels[current_channel].points;
	point *prev = NULL, *next = NULL;
	for (GList *iter = points; iter != NULL; iter = iter->next) {
		point *p = (point *) iter->data;
		if (p->x < selected_point->x) prev = p;
		else if (p->x > selected_point->x) { next = p; break; }
	}
	if (prev && prev->x >= xpos)      selected_point->x = prev->x + 0.00001;
	else if (next && next->x <= xpos) selected_point->x = next->x - 0.00001;
	else                              selected_point->x = xpos;
	_update_entry_text();
	gtk_widget_queue_draw(widget);
}

void on_curves_drawingarea_leave_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	set_cursor("default");
}

/* GTK4 press: primary-click selects an existing point or adds a new one;
 * secondary-click removes a point. */
static void curves_da_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)n_press; (void)user_data;
	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
	gdouble xpos = x / (gdouble) gtk_widget_get_width(curves_drawingarea);
	gdouble ypos = 1 - y / (gdouble) gtk_widget_get_height(curves_drawingarea);

	double click_precision = 0.02 / gtk_adjustment_get_value(curves_adj_zoom);

	GList **points_ptr = &gui_channels[current_channel].points;
	GList *points = *points_ptr;

	guint button = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	if (button == GDK_BUTTON_PRIMARY) {
		gboolean was_point_clicked = FALSE;
		for (GList *iter = points; iter != NULL; iter = iter->next) {
			point *p = (point *) iter->data;
			if (fabs(p->x - xpos) < click_precision && fabs(p->y - ypos) < click_precision) {
				selected_point = p;
				set_cursor("hand1");
				was_point_clicked = TRUE;
				break;
			}
		}
		/* No existing point under the cursor → insert a new one. */
		if (!was_point_clicked && g_list_length(points) < MAX_POINTS) {
			for (GList *iter = points; iter != NULL; iter = iter->next) {
				point *p = (point *) iter->data;
				if (p->x == xpos) {
					GtkWidget *popover = popover_new(curves_drawingarea,
					    "A point already exists at this x position");
					gtk_widget_set_visible(popover, TRUE);
					return;
				}
			}
			point *p1 = g_new(point, 1);
			p1->x = xpos;
			p1->y = ypos;
			*points_ptr = g_list_insert_sorted(points, p1, (GCompareFunc) compare_points);
			selected_point = p1;
		}
		is_drawingarea_pressed = TRUE;
		gtk_widget_queue_draw(widget);
		_update_entry_text();
	} else if (button == GDK_BUTTON_SECONDARY && !is_drawingarea_pressed) {
		/* Right-click on a point removes it (if more than two remain). */
		for (GList *iter = points; iter != NULL; iter = iter->next) {
			point *p = (point *) iter->data;
			if (fabs(p->x - xpos) < click_precision && fabs(p->y - ypos) < click_precision) {
				if (g_list_length(points) > 2) {
					GList *new_selected_point;
					if (p == g_list_nth_data(points, 0))
						new_selected_point = iter->next;
					else
						new_selected_point = iter->prev;
					*points_ptr = g_list_remove(points, p);
					g_free(p);
					selected_point = (point *) new_selected_point->data;
					gtk_widget_queue_draw(widget);
					_update_entry_text();
				}
				break;
			}
		}
	}
}

/* GTK4 release: clear the drag state, refresh the live histogram and rebuild
 * the previewed image. */
static void curves_da_released(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)gesture; (void)n_press; (void)x; (void)y; (void)user_data;
	is_drawingarea_pressed = FALSE;
	update_display_histogram_from_curve();
	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image();
}

/* GTK4 prev / next: walk the selected_point through the per-channel list. */
static void curves_prev_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)gesture; (void)x; (void)y; (void)user_data;
	if (n_press != 1) return;
	GList *points = gui_channels[current_channel].points;
	int idx = g_list_index(points, selected_point);
	if (idx > 0) {
		selected_point = g_list_nth_data(points, idx - 1);
		gtk_widget_queue_draw(curves_drawingarea);
		_update_entry_text();
	}
}

static void curves_next_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)gesture; (void)x; (void)y; (void)user_data;
	if (n_press != 1) return;
	GList *points = gui_channels[current_channel].points;
	int idx = g_list_index(points, selected_point);
	if (idx < (int)g_list_length(points) - 1) {
		selected_point = g_list_nth_data(points, idx + 1);
		gtk_widget_queue_draw(curves_drawingarea);
		_update_entry_text();
	}
}

void on_curves_id_entry_activate(GtkEntry *entry, gpointer user_data) {
	GList *points = gui_channels[current_channel].points;
	int entered_index = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
	int len = g_list_length(points);
	if (entered_index >= 0 && entered_index < len) selected_point = g_list_nth_data(points, entered_index);
	else if (entered_index >= len) selected_point = g_list_nth_data(points, len - 1);
	else selected_point = g_list_nth_data(points, 0);
	gtk_widget_queue_draw(curves_drawingarea);
	_update_entry_text();
}

gboolean on_curves_id_entry_focus_out_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	on_curves_id_entry_activate(curves_id_entry, user_data);
	return FALSE;
}

void on_curves_x_entry_activate(GtkEntry *entry, gpointer user_data) {
	float curve_x = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
	GList *points = gui_channels[current_channel].points;
	GList *iter; point *p;
	for (iter = points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		if (p != selected_point && p->x == curve_x) {
			gchar *str = g_strdup_printf("%8.7f", selected_point->x);
			gtk_editable_set_text(GTK_EDITABLE(entry), str);
			g_free(str);
			GtkWidget *popover = popover_new(curves_drawingarea, "A point already exists at this x position");
			gtk_widget_set_visible(popover, TRUE);
			return;
		}
	}
	if (curve_x < 0 || curve_x > 1) {
		gchar *str = g_strdup_printf("%8.7f", selected_point->x);
		gtk_editable_set_text(GTK_EDITABLE(entry), str);
		g_free(str);
		return;
	}
	selected_point->x = curve_x;
	gui_channels[current_channel].points = g_list_sort(points, (GCompareFunc) compare_points);
	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image(); _update_entry_text();
}

gboolean on_curves_x_entry_focus_out_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	on_curves_x_entry_activate(curves_x_entry, user_data);
	return FALSE;
}

void on_curves_y_entry_activate(GtkEntry *entry, gpointer user_data) {
	float curve_y = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
	curve_y = (curve_y < 0) ? 0 : (curve_y > 1) ? 1 : curve_y;
	selected_point->y = curve_y;
	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image();
	gchar *str = g_strdup_printf("%8.7f", curve_y);
	gtk_editable_set_text(GTK_EDITABLE(entry), str);
	g_free(str);
}

gboolean on_curves_y_entry_focus_out_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	on_curves_y_entry_activate(curves_y_entry, user_data);
	return FALSE;
}

void on_curves_preview_toggled(GtkCheckButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!siril_toggle_get_active(GTK_WIDGET(curves_preview_check))) {
		cancel_and_wait_for_preview();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();
		curves_update_image();
	}
}

void on_curves_log_check_toggled(GtkCheckButton *button, gpointer user_data) {
	gtk_widget_queue_draw(curves_drawingarea);
}

void on_curves_apply_stage_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa()) return;

	set_cursor_waiting(TRUE);

	gboolean preview_active = siril_toggle_get_active(GTK_WIDGET(curves_preview_check));

	if (preview_active && !gui.roi.active) {
		// gfit already has backup + current curves; just freeze it
		stage_stack = g_list_append(stage_stack, snapshot_current_state());
		copy_gfit_to_backup();
	} else {
		// No live preview: apply synchronously to build the new backup
		curve_state_snapshot *snap = snapshot_current_state();
		stage_stack = g_list_append(stage_stack, snap);
		copy_backup_to_gfit();
		apply_snapshot_to_fit(gfit, snap);
		copy_gfit_to_backup();
		notify_gfit_data_modified();
		gfit_modified_update_gui();
	}

	// Refresh histograms from the new backup
	compute_histo_for_fit(fit);
	clear_display_histogram();
	for (int i = 0; i < (int)fit->naxes[2]; i++)
		display_histogram[i] = gsl_histogram_clone(com.layers_hist[i]);
	if (fit->naxes[2] == 3) {
		if (curves_huebuf) { free(curves_huebuf); curves_huebuf = NULL; }
		if (curves_satbuf) { free(curves_satbuf); curves_satbuf = NULL; }
		if (curves_lumbuf) { free(curves_lumbuf); curves_lumbuf = NULL; }
		curves_compute_special_histograms();
	}
	refresh_hist_input_from_fit();

	// Reset curves to identity for the next stage
	init_all_curves();
	selected_point = (point *) gui_channels[current_channel].points->data;
	update_range_ui_from_state();
	_update_entry_text();
	update_stage_buttons();
	gtk_widget_queue_draw(curves_drawingarea);

	// Trigger preview redraw (backup + identity = backup = accumulated stages)
	curves_update_image();
	set_cursor_waiting(FALSE);
}

void on_curves_undo_stage_clicked(GtkButton *button, gpointer user_data) {
	if (!stage_stack) return;

	GList *last = g_list_last(stage_stack);
	free_snapshot((curve_state_snapshot *) last->data);
	stage_stack = g_list_delete_link(stage_stack, last);

	curves_rebuild_from_stages();
}

void on_curves_pipette_button_toggled(GtkToggleButton *button, gpointer user_data) {
	if (siril_toggle_get_active(GTK_WIDGET(button))) {
		mouse_status = MOUSE_ACTION_CURVES_PIPETTE;
	} else {
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		pipette_value = -1.0f;
		gtk_widget_queue_draw(curves_drawingarea);
	}
}

void curves_handle_pipette_click(int x, int y) {
	if (!gfit || x < 0 || y < 0 || x >= (int)gfit->rx || y >= (int)gfit->ry)
		return;

	float r, g, b;
	size_t idx = (size_t) y * gfit->rx + x;

	if (gfit->type == DATA_FLOAT) {
		if (gfit->naxes[2] == 3) {
			r = gfit->fpdata[0][idx];
			g = gfit->fpdata[1][idx];
			b = gfit->fpdata[2][idx];
		} else {
			r = g = b = gfit->fdata[idx];
		}
	} else {
		double inv = 1.0 / USHRT_MAX_DOUBLE;
		if (gfit->naxes[2] == 3) {
			r = (float)(gfit->pdata[0][idx] * inv);
			g = (float)(gfit->pdata[1][idx] * inv);
			b = (float)(gfit->pdata[2][idx] * inv);
		} else {
			r = g = b = (float)(gfit->data[idx] * inv);
		}
	}

	switch (current_channel) {
		case CHAN_RGB_K:
			pipette_value = (r + g + b) / 3.0f;
			break;
		case CHAN_R:
			pipette_value = r;
			break;
		case CHAN_G:
			pipette_value = g;
			break;
		case CHAN_B:
			pipette_value = b;
			break;
		case CHAN_L: {
			// Linear luminance approximation (close to CIE L*)
			float lum = 0.2126f * r + 0.7152f * g + 0.0722f * b;
			pipette_value = CLAMP(lum, 0.0f, 1.0f);
			break;
		}
		case CHAN_S: {
			float h, s, l;
			rgb_to_hslf(r, g, b, &h, &s, &l);
			pipette_value = s;
			break;
		}
		case CHAN_C: {
			// Chroma approximation: same as the histogram computation
			float h, s, l;
			rgb_to_hslf(r, g, b, &h, &s, &l);
			float l_factor = 1.0f - fabsf(2.0f * l - 1.0f);
			pipette_value = CLAMP(s * l_factor, 0.0f, 1.0f);
			break;
		}
		default:
			pipette_value = -1.0f;
			break;
	}

	if (curves_drawingarea)
		gtk_widget_queue_draw(curves_drawingarea);
}
