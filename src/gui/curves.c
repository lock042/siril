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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "curves.h"
#include "histogram.h"
#include "filters/curve_transform.h"

// ---------------------------------------------------------------------------
// FORWARD DECLARATIONS
// ---------------------------------------------------------------------------
static int curves_process_with_worker(gboolean for_preview, gboolean for_roi);
static gboolean is_curves_log_scale(void);
static void update_range_ui_from_state(void);
static int curves_update_preview(void);
static void curves_update_image(void);
static void _update_entry_text(void);
static void _initialize_clip_text(void);

// Forward declaration of Idle functions defined in this file
gboolean curve_preview_idle(gpointer p);
gboolean curve_apply_idle(gpointer p);

// ---------------------------------------------------------------------------
// UI CACHE
// ---------------------------------------------------------------------------

static GtkAdjustment *curves_adj_zoom = NULL;
static GtkComboBoxText *curves_interpolation_combo = NULL;
static GtkEntry *curves_id_entry = NULL, *curves_x_entry = NULL, *curves_y_entry = NULL, *curves_clip_low = NULL, *curves_clip_high = NULL, *curves_seq_entry = NULL;
static GtkGrid *curves_point_grid = NULL;
static GtkToggleButton *curves_sequence_check = NULL, *curves_preview_check = NULL, *curves_log_check = NULL;
static GtkToggleToolButton *curves_grid_toggle = NULL;
static GtkWidget *curves_drawingarea = NULL, *curves_viewport = NULL, *curves_dialog = NULL;
static GtkComboBoxText *curves_channel_combo = NULL;
static GtkCheckButton *curves_range_check = NULL;
static GtkScale *curves_range_min_scale = NULL;
static GtkScale *curves_range_max_scale = NULL;
static GtkScale *curves_range_feather_scale = NULL;
static GtkWidget *curves_undo_stage_btn = NULL;
static GtkLabel *curves_stage_lbl = NULL;

// ---------------------------------------------------------------------------
// STATE
// ---------------------------------------------------------------------------

static curve_channel_config gui_channels[CHAN_COUNT];
static GList *undo_stack = NULL; // Stack of curve_state_snapshot*
static int stage_count = 0;

static enum curve_channel current_channel = CHAN_RGB_K;
static enum curve_algorithm algorithm = AKIMA_SPLINE;

static gboolean closing = FALSE;
static gboolean do_channel[3];
static fits *fit = NULL;
static long clipped[] = {0, 0};
static gsl_histogram *display_histogram[MAXVPORT] = {NULL, NULL, NULL, NULL};
static gboolean is_drawingarea_pressed = FALSE;
static cmsHPROFILE original_icc = NULL;
static gboolean single_image_stretch_applied = FALSE;
static point *selected_point = NULL;
static double point_size = 5.0;

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
// UNDO / STAGE MANAGEMENT
// ---------------------------------------------------------------------------

static void free_snapshot(curve_state_snapshot *snap) {
	if (!snap) return;
	for(int i=0; i<CHAN_COUNT; i++) {
		if (snap->channels[i].points)
			g_list_free_full(snap->channels[i].points, g_free);
	}
	g_free(snap);
}

static void push_undo_state() {
	curve_state_snapshot *snap = g_new0(curve_state_snapshot, 1);
	snap->algorithm = algorithm;

	for(int i=0; i<CHAN_COUNT; i++) {
		snap->channels[i].active = gui_channels[i].active;
		snap->channels[i].range_enabled = gui_channels[i].range_enabled;
		snap->channels[i].lum_min = gui_channels[i].lum_min;
		snap->channels[i].lum_max = gui_channels[i].lum_max;
		snap->channels[i].feather = gui_channels[i].feather;
		if (gui_channels[i].points)
			snap->channels[i].points = g_list_copy_deep(gui_channels[i].points, (GCopyFunc)copy_point, NULL);
	}

	undo_stack = g_list_prepend(undo_stack, snap);
	stage_count++;

	char buf[32];
	snprintf(buf, 32, "Stages: %d", stage_count);
	gtk_label_set_text(curves_stage_lbl, buf);
	gtk_widget_set_sensitive(curves_undo_stage_btn, TRUE);
}

static void pop_undo_state() {
	if (!undo_stack) return;

	curve_state_snapshot *snap = (curve_state_snapshot*)undo_stack->data;
	undo_stack = g_list_remove(undo_stack, snap);

	// Restore State
	algorithm = snap->algorithm;
	gtk_combo_box_set_active(GTK_COMBO_BOX(curves_interpolation_combo), algorithm);

	for(int i=0; i<CHAN_COUNT; i++) {
		// Clear current
		if (gui_channels[i].points) g_list_free_full(gui_channels[i].points, g_free);

		// Restore
		gui_channels[i].active = snap->channels[i].active;
		gui_channels[i].range_enabled = snap->channels[i].range_enabled;
		gui_channels[i].lum_min = snap->channels[i].lum_min;
		gui_channels[i].lum_max = snap->channels[i].lum_max;
		gui_channels[i].feather = snap->channels[i].feather;

		// Move ownership
		gui_channels[i].points = snap->channels[i].points;
		snap->channels[i].points = NULL;
	}

	free_snapshot(snap);
	stage_count--;

	char buf[32];
	snprintf(buf, 32, "Stages: %d", stage_count);
	gtk_label_set_text(curves_stage_lbl, buf);
	gtk_widget_set_sensitive(curves_undo_stage_btn, (undo_stack != NULL));

	// Refresh UI
	selected_point = (point*) gui_channels[current_channel].points->data;
	update_range_ui_from_state();
	_update_entry_text();
	curves_update_image();
	gtk_widget_queue_draw(curves_drawingarea);
}

static void clear_undo_stack() {
	if (undo_stack) {
		g_list_free_full(undo_stack, (GDestroyNotify)free_snapshot);
		undo_stack = NULL;
	}
	stage_count = 0;
	if (curves_stage_lbl) gtk_label_set_text(curves_stage_lbl, "Stages: 0");
	if (curves_undo_stage_btn) gtk_widget_set_sensitive(curves_undo_stage_btn, FALSE);
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
	for(int i=0; i<CHAN_COUNT; i++) init_single_channel_config(i);
	current_channel = CHAN_RGB_K;
	selected_point = (point *) gui_channels[CHAN_RGB_K].points->data;
}

static void update_range_ui_from_state(void) {
	if (!curves_range_check) return;
	g_signal_handlers_block_by_func(curves_range_check, on_curve_check_range_button_toggled, NULL);
	g_signal_handlers_block_by_func(curves_range_min_scale, on_curves_range_value_changed, NULL);
	g_signal_handlers_block_by_func(curves_range_max_scale, on_curves_range_value_changed, NULL);
	g_signal_handlers_block_by_func(curves_range_feather_scale, on_curves_feather_value_changed, NULL);

	curve_channel_config *cfg = &gui_channels[current_channel];
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(curves_range_check), cfg->range_enabled);
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
}

// ---------------------------------------------------------------------------
// UI SETUP
// ---------------------------------------------------------------------------

void curves_dialog_init_statics() {
	if (curves_adj_zoom == NULL) {
		curves_adj_zoom = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "curves_adj_zoom"));
		curves_interpolation_combo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "curves_interpolation_combo"));
		curves_dialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_dialog"));
		curves_drawingarea = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_drawingarea"));
		curves_id_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_id_entry"));
		curves_x_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_x_entry"));
		curves_y_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_y_entry"));
		curves_clip_low = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_clip_low"));
		curves_clip_high = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_clip_high"));
		curves_seq_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_seq_entry"));
		curves_point_grid = GTK_GRID(gtk_builder_get_object(gui.builder, "curves_point_grid"));
		curves_sequence_check = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_sequence_check"));
		curves_preview_check = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_preview_check"));
		curves_log_check = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_log_check"));
		curves_grid_toggle = GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object(gui.builder, "curves_grid_toggle"));
		curves_viewport = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_viewport"));

		curves_range_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "curve_check_range_button"));
		curves_range_min_scale = GTK_SCALE(gtk_builder_get_object(gui.builder, "curve_slider_min"));
		curves_range_max_scale = GTK_SCALE(gtk_builder_get_object(gui.builder, "curve_slider_max"));
		curves_range_feather_scale = GTK_SCALE(gtk_builder_get_object(gui.builder, "curve_slider_feather"));
		curves_channel_combo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "curves_channel_combo"));

		curves_undo_stage_btn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_undo_stage_btn"));
		curves_stage_lbl = GTK_LABEL(gtk_builder_get_object(gui.builder, "curves_undo_stage_lbl"));
	}
}

// ---------------------------------------------------------------------------
// DRAWING (THE COOL STUFF)
// ---------------------------------------------------------------------------

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

	if (gtk_toggle_tool_button_get_active(curves_grid_toggle)) {
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
// LOGIC
// ---------------------------------------------------------------------------

static gboolean is_curves_log_scale(void) {
	return (gtk_toggle_button_get_active(curves_log_check));
}

gboolean redraw_curves(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int width = gtk_widget_get_allocated_width(curves_drawingarea);
	int height = gtk_widget_get_allocated_height(curves_drawingarea);
	if (height == 1) return FALSE;

	draw_background_gradient(cr, width, height);

	for (int i = 0; i < MAXVPORT; i++) {
		if (current_channel == CHAN_R && i != 0) continue;
		if (current_channel == CHAN_G && i != 1) continue;
		if (current_channel == CHAN_B && i != 2) continue;
		display_histo(display_histogram[i], cr, i, width, height, 1.0, 1.0, FALSE, is_curves_log_scale());
	}

	draw_range_overlays(cr, width, height);
	draw_ghost_curves(cr, width, height);
	draw_active_curve(cr, width, height);
	draw_curve_points(cr, width, height);

	return FALSE;
}

static int curves_update_preview(void) {
	fit = gui.roi.active ? &gui.roi.fit : gfit;
	if (!closing) {
		copy_backup_to_gfit();
		// Reset is handled in curves_process_with_worker via local param
		curves_process_with_worker(TRUE, gui.roi.active);
	}
	return 0;
}

static void curves_update_image(void) {
	set_cursor_waiting(TRUE);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = curves_update_preview;
	param->show_preview = gtk_toggle_button_get_active(curves_preview_check);
	notify_update((gpointer) param);
	set_cursor_waiting(FALSE);
}

static void _initialize_clip_text(void) {
	gtk_entry_set_text(curves_clip_low, "0.000%");
	gtk_entry_set_text(curves_clip_high, "0.000%");
}

// Helper to update UI from global stats
static void update_clip_ui_from_stats(void) {
	if (!fit) return;
	size_t total_pixels = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	if (total_pixels == 0) total_pixels = 1; // Safety

	double tmp;
	char buffer[32];

	// High
	tmp = (double) clipped[1] * 100.0 / (double) total_pixels;
	if (tmp < 0) tmp = 0;
	snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_entry_set_text(curves_clip_high, buffer);

	// Low
	tmp = (double) clipped[0] * 100.0 / (double) total_pixels;
	if (tmp < 0) tmp = 0;
	snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_entry_set_text(curves_clip_low, buffer);
}

static void reset_cursors_and_values(gboolean full_reset) {
	_initialize_clip_text();
	gtk_adjustment_set_value(curves_adj_zoom, 1.0);
	gtk_entry_set_text(curves_seq_entry, "curve_");
	gtk_toggle_button_set_active(curves_preview_check, TRUE);

	if (full_reset) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(curves_interpolation_combo), AKIMA_SPLINE);
		algorithm = AKIMA_SPLINE;
		if (curves_channel_combo) gtk_combo_box_set_active(GTK_COMBO_BOX(curves_channel_combo), 0);
		init_all_curves();
		clear_undo_stack();
	} else {
		init_single_channel_config(current_channel);
		selected_point = (point*) gui_channels[current_channel].points->data;
	}
	update_range_ui_from_state();
	_update_entry_text();
	update_gfit_curves_histogram_if_needed();
}

static void set_histogram(gsl_histogram *histo, int layer) {
	g_assert(layer >= 0 && layer < MAXVPORT);
	if (com.layers_hist[layer]) gsl_histogram_free(com.layers_hist[layer]);
	com.layers_hist[layer] = histo;
}

static void clear_display_histogram(void) {
	if (display_histogram[0]) {
		for (int i = 0; i < fit->naxes[2]; i++) {
			gsl_histogram_free(display_histogram[i]);
			display_histogram[i] = NULL;
		}
	}
}

static void curves_startup(void) {
	init_all_curves();
	add_roi_callback(curves_histogram_change_between_roi_and_image);
	roi_supported(TRUE);
	copy_gfit_to_backup();
	compute_histo_for_fit(fit);
	for (int i = 0; i < fit->naxes[2]; i++)
		display_histogram[i] = gsl_histogram_clone(com.layers_hist[i]);
}

static void curves_close(gboolean update_image_if_needed, gboolean revert_icc_profile) {
	clear_undo_stack();

	for (int i = 0; i < fit->naxes[2]; i++) {
		set_histogram(display_histogram[i], i);
		display_histogram[i] = NULL;
	}
	if (is_preview_active() && !copy_backup_to_gfit() && update_image_if_needed) {
		set_cursor_waiting(TRUE);
		notify_gfit_modified();
	}
	if (revert_icc_profile && !single_image_stretch_applied) {
		if (gfit->icc_profile) cmsCloseProfile(gfit->icc_profile);
		gfit->icc_profile = copyICCProfile(original_icc);
		color_manage(gfit, gfit->icc_profile != NULL);
	}
	clear_backup();
	clear_display_histogram();
	roi_supported(FALSE);
	remove_roi_callback(curves_histogram_change_between_roi_and_image);
}

static int curves_process_with_worker(gboolean for_preview, gboolean for_roi) {
	struct curve_params *params = new_curve_params();
	if (!params) { PRINT_ALLOC_ERR; return 1; }

	// Deep Copy Configuration for Thread Safety
	for(int i=0; i<CHAN_COUNT; i++) {
		params->channels[i].active = gui_channels[i].active;
		params->channels[i].range_enabled = gui_channels[i].range_enabled;
		params->channels[i].lum_min = gui_channels[i].lum_min;
		params->channels[i].lum_max = gui_channels[i].lum_max;
		params->channels[i].feather = gui_channels[i].feather;

		if (gui_channels[i].points) {
			params->channels[i].points = g_list_copy_deep(gui_channels[i].points, (GCopyFunc)copy_point, NULL);
		} else {
			params->channels[i].points = NULL;
		}
	}

	// Reset global clipping and link it to the worker params
	clipped[0] = 0;
	clipped[1] = 0;
	params->clipped_count = clipped;

	params->algorithm = algorithm;
	params->target_channel = current_channel;
	params->fit = for_roi ? &gui.roi.fit : gfit;
	params->verbose = !for_preview;
	params->for_preview = for_preview;

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = params->fit;
	args->mem_ratio = 2.0f;
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
	if (args->retval == 0) notify_gfit_modified();
	free_generic_img_args(args);
	return FALSE;
}

gboolean curve_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;

	// Update UI from stats gathered by worker
	update_clip_ui_from_stats();

	stop_processing_thread();
	if (args->retval == 0) notify_gfit_modified();
	free_generic_img_args(args);
	return FALSE;
}

void _update_entry_text(void) {
	if (selected_point == NULL) return;
	gchar *buffer;
	buffer = g_strdup_printf("%.7f", selected_point->x);
	gtk_entry_set_text(curves_x_entry, buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", selected_point->y);
	gtk_entry_set_text(curves_y_entry, buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%d", g_list_index(gui_channels[current_channel].points, selected_point));
	gtk_entry_set_text(curves_id_entry, buffer);
	g_free(buffer);
}

void update_gfit_curves_histogram_if_needed(void) {
	invalidate_gfit_histogram();
	if (gtk_widget_get_visible(curves_dialog)) {
		compute_histo_for_fit(fit);
		gtk_widget_queue_draw(curves_drawingarea);
	}
}

void curves_histogram_change_between_roi_and_image(void) {
	fit = gui.roi.active ? &gui.roi.fit : gfit;
	gui.roi.operation_supports_roi = TRUE;
	curves_update_image();
}

// ---------------------------------------------------------------------------
// EVENTS
// ---------------------------------------------------------------------------

void on_curves_channel_combo_changed(GtkComboBox *widget, gpointer user_data) {
	enum curve_channel new_chan = gtk_combo_box_get_active(widget);
	if (new_chan != current_channel) {
		switch_channel_view(new_chan);
		gtk_widget_queue_draw(curves_drawingarea);
		curves_update_image();
		_update_entry_text();
	}
}

void on_curve_check_range_button_toggled(GtkToggleButton *button, gpointer user_data) {
	gboolean active = gtk_toggle_button_get_active(button);
	gui_channels[current_channel].range_enabled = active;

	GtkWidget *expander = lookup_widget("curve_range_expander");
	gtk_expander_set_resize_toplevel((GtkExpander*)expander, TRUE);
	gtk_widget_set_sensitive(expander, active);
	gtk_expander_set_expanded((GtkExpander*)expander,  active);

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
	curves_startup();
	_initialize_clip_text();
	reset_cursors_and_values(TRUE);
	compute_histo_for_fit(fit);
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
	set_cursor_waiting(FALSE);
}

void on_curves_apply_button_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa()) return;

	if (gtk_toggle_button_get_active(curves_sequence_check) && sequence_is_loaded()) {
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

		args->seq_entry = strdup(gtk_entry_get_text(curves_seq_entry));
		args->seq = &com.seq;
		if (args->seq_entry && args->seq_entry[0] == '\0') {
			free(args->seq_entry);
			args->seq_entry = strdup("curve_");
		}
		curves_close(FALSE, TRUE);
		siril_close_dialog("curves_dialog");
		gtk_toggle_button_set_active(curves_sequence_check, FALSE);
		apply_curve_to_sequence(args);
	} else {
		push_undo_state();

		fit = gfit;
		fits undo_fit = {0};
		memcpy(&undo_fit, get_preview_gfit_backup(), sizeof(fits));
		undo_fit.icc_profile = original_icc;
		undo_fit.color_managed = original_icc != NULL;
		copy_backup_to_gfit();

		curves_process_with_worker(FALSE, FALSE);

		single_image_stretch_applied = TRUE;
		populate_roi();
		clear_backup();
		clear_display_histogram();
		curves_startup();

		reset_cursors_and_values(FALSE);
		set_cursor("default");
	}
}

void on_curves_undo_stage_clicked(GtkButton *button, gpointer user_data) {
	if (!undo_stack) return;

	set_cursor_waiting(TRUE);

	// FIX LINKER ERROR: Removed explicit call to cmd_interpreter/undo.
	// This button now acts as "Restore Previous Curve State"
	// To undo image changes, user must use the main Undo button (Ctrl+Z).
	// siril_log_message(_("Undo Stage: Restoring curve parameters... (Use Ctrl+Z to undo image changes)"));

	pop_undo_state();

	if (stage_count == 0) single_image_stretch_applied = FALSE;

	set_cursor_waiting(FALSE);
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
	current_width = gtk_widget_get_allocated_width(curves_viewport);
	target_width = (int) (((double) current_width) * zoom);
	gtk_widget_set_size_request(curves_drawingarea, target_width, 281);
	gtk_widget_queue_draw(curves_drawingarea);
}

void on_curves_interpolation_combo_changed(GtkComboBox *widget, gpointer user_data) {
	algorithm = gtk_combo_box_get_active(widget);
	curves_update_image();
	gtk_widget_queue_draw(curves_drawingarea);
}

void setup_curve_dialog() {
	gtk_window_set_title(GTK_WINDOW(curves_dialog), _("Curve Transformation")); 
	gtk_widget_set_visible(GTK_WIDGET(curves_point_grid), TRUE);
	gtk_window_resize(GTK_WINDOW(curves_dialog), 1, 1);
}

void toggle_curves_window_visibility() {
	if (gtk_widget_get_visible(lookup_widget("histogram_dialog")))
	siril_close_dialog("histogram_dialog");

	curves_dialog_init_statics();

	if (gtk_widget_get_visible(curves_dialog)) {
		set_cursor_waiting(TRUE);
		curves_close(TRUE, TRUE);
		set_cursor_waiting(FALSE);
		siril_close_dialog("curves_dialog");
	} else {
		for (int i = 0; i < 3; i++) do_channel[i] = TRUE;
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

gboolean on_curves_drawingarea_motion_notify_event(GtkWidget *widget, GdkEventMotion *event, gpointer user_data) {
	if (!is_drawingarea_pressed) return FALSE;
	gdouble xpos = ((GdkEventMotion *) event)->x / (gdouble) gtk_widget_get_allocated_width(curves_drawingarea);
	gdouble ypos = 1 - ((GdkEventMotion *) event)->y / (gdouble) gtk_widget_get_allocated_height(curves_drawingarea);
	xpos = (xpos < 0) ? 0 : (xpos > 1) ? 1 : xpos;
	selected_point->y = (ypos < 0) ? 0 : (ypos > 1) ? 1 : ypos;

	GList *points = gui_channels[current_channel].points;
	GList *iter;
	point *p;
	point *prev = NULL;
	point *next = NULL;
	for (iter = points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		if (p->x < selected_point->x) prev = p;
		else if (p->x > selected_point->x) { next = p; break; }
	}
	if (prev && prev->x >= xpos) selected_point->x = prev->x + 0.00001;
	else if (next && next->x <= xpos) selected_point->x = next->x - 0.00001;
	else selected_point->x = xpos;

	_update_entry_text();
	gtk_widget_queue_draw(widget);
	return FALSE;
}

void on_curves_drawingarea_leave_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	set_cursor("default");
}

gboolean on_curves_drawingarea_button_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	gdouble xpos = ((GdkEventButton *) event)->x / (gdouble) gtk_widget_get_allocated_width(curves_drawingarea);
	gdouble ypos = 1 - ((GdkEventButton *) event)->y / (gdouble) gtk_widget_get_allocated_height(curves_drawingarea);
	double click_precision = 0.02 / gtk_adjustment_get_value(curves_adj_zoom);

	GList **points_ptr = &gui_channels[current_channel].points;
	GList *points = *points_ptr;

	if (event->button == GDK_BUTTON_PRIMARY) {
		gboolean was_point_clicked = FALSE;
		GList *iter;
		point *p;
		for (iter = points; iter != NULL; iter = iter->next) {
			p = (point *) iter->data;
			if (fabs(p->x - xpos) < click_precision && fabs(p->y - ypos) < click_precision) {
				selected_point = p;
				set_cursor("hand1");
				was_point_clicked = TRUE;
				break;
			}
		}

		if (!was_point_clicked && g_list_length(points) < MAX_POINTS) {
			for (iter = points; iter != NULL; iter = iter->next) {
				p = (point *) iter->data;
				if (p->x == xpos) return FALSE;
			}
			point *p1 = g_new(point, 1);
			p1->x = xpos; p1->y = ypos;
			*points_ptr = g_list_insert_sorted(points, p1, (GCompareFunc) compare_points);
			selected_point = p1;
		}
		is_drawingarea_pressed = TRUE;
		gtk_widget_queue_draw(widget);
		_update_entry_text();
	} else if (event->button == GDK_BUTTON_SECONDARY && !is_drawingarea_pressed) {
		GList *iter;
		point *p;
		for (iter = points; iter != NULL; iter = iter->next) {
			p = (point *) iter->data;
			if (fabs(p->x - xpos) < click_precision && fabs(p->y - ypos) < click_precision) {
				if (g_list_length(points) > 2) {
					GList *new_selected_point;
					if	(p == g_list_nth_data(points, 0)) new_selected_point = iter->next;
					else new_selected_point = iter->prev;
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
	return FALSE;
}

gboolean on_curves_drawingarea_button_release_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	is_drawingarea_pressed = FALSE;
	curves_update_image();
	return FALSE;
}

gboolean on_curves_prev_button_clicked(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	if (event->type != GDK_BUTTON_PRESS) return TRUE;
	GList *points = gui_channels[current_channel].points;
	int selected_point_index = g_list_index(points, selected_point);
	if (selected_point_index > 0) {
		selected_point = g_list_nth_data(points, selected_point_index - 1);
		gtk_widget_queue_draw(curves_drawingarea);
		_update_entry_text();
	}
	return FALSE;
}

gboolean on_curves_next_button_clicked(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	if (event->type != GDK_BUTTON_PRESS) return TRUE;
	GList *points = gui_channels[current_channel].points;
	int selected_point_index = g_list_index(points, selected_point);
	if (selected_point_index < g_list_length(points) - 1) {
		selected_point = g_list_nth_data(points, selected_point_index + 1);
		gtk_widget_queue_draw(curves_drawingarea);
		_update_entry_text();
	}
	return FALSE;
}

void on_curves_id_entry_activate(GtkEntry *entry, gpointer user_data) {
	GList *points = gui_channels[current_channel].points;
	int entered_index = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
	int len = g_list_length(points);
	if (entered_index >= 0 && entered_index < len) selected_point = g_list_nth_data(points, entered_index);
	else if (entered_index >= len) selected_point = g_list_nth_data(points, len - 1);
	else selected_point = g_list_nth_data(points, 0);
	gtk_widget_queue_draw(curves_drawingarea);
	_update_entry_text();
}
gboolean on_curves_id_entry_focus_out_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) { on_curves_id_entry_activate(curves_id_entry, user_data); return FALSE; }

void on_curves_x_entry_activate(GtkEntry *entry, gpointer user_data) {
	float curve_x = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
	GList *points = gui_channels[current_channel].points;
	GList *iter; point *p;
	for (iter = points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		if (p != selected_point && p->x == curve_x) return;
	}
	if (curve_x < 0 || curve_x > 1) return;
	selected_point->x = curve_x;
	gui_channels[current_channel].points = g_list_sort(points, (GCompareFunc) compare_points);
	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image(); _update_entry_text();
}
gboolean on_curves_x_entry_focus_out_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) { on_curves_x_entry_activate(curves_x_entry, user_data); return FALSE; }

void on_curves_y_entry_activate(GtkEntry *entry, gpointer user_data) {
	float curve_y = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
	curve_y = (curve_y < 0) ? 0 : (curve_y > 1) ? 1 : curve_y;
	selected_point->y = curve_y;
	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image();
	gchar *str = g_strdup_printf("%8.7f", curve_y);
	gtk_entry_set_text(entry, str);
	g_free(str);
}
gboolean on_curves_y_entry_focus_out_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) { on_curves_y_entry_activate(curves_y_entry, user_data); return FALSE; }

void on_curves_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(curves_preview_check)) {
		waiting_for_thread();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();
		curves_update_image();
	}
}

void on_curves_log_check_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_widget_queue_draw(curves_drawingarea);
}