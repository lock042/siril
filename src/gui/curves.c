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

#include <string.h>
#include <math.h>
#include <float.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "algos/colors.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/siril_preview.h"
#include "gui/registration_preview.h"
#include "core/undo.h"
#include "curves.h"
#include "histogram.h"
#include "filters/curve_transform.h"

/* The gsl_histogram, documented here:
 * http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_21.html
 * is able to group values into bins, it does not need to handle all values
 * separately. That is useful for display purposes, but not currently used.
 */

// Cache UI elements
static GtkAdjustment *curves_adj_zoom = NULL;
static GtkComboBoxText *curves_interpolation_combo = NULL;
static GtkEntry *curves_id_entry = NULL, *curves_x_entry = NULL, *curves_y_entry = NULL, *curves_clip_low = NULL, *curves_clip_high = NULL, *curves_seq_entry = NULL;
static GtkGrid *curves_point_grid = NULL;
static GtkToggleButton *curves_sequence_check = NULL, *curves_preview_check = NULL, *curves_log_check = NULL;
static GtkToggleToolButton *curves_red_toggle = NULL, *curves_green_toggle = NULL, *curves_blue_toggle = NULL, *curves_grid_toggle = NULL;
static GtkWidget *curves_drawingarea = NULL, *curves_viewport = NULL, *curves_dialog = NULL;

static gboolean closing = FALSE;

static gboolean do_channel[3];

static fits *fit = &gfit;

// static float graph_height = 0.f;	// the max value of all bins
static long clipped[] = {0, 0};

// The histogram displayed on the drawingarea
static gsl_histogram *display_histogram[MAXVPORT] = {NULL, NULL, NULL, NULL};

static gboolean is_drawingarea_pressed = FALSE;

// compare function for points
int compare_points(const void *a, const void *b) {
	point *point_a = (point *) a;
	point *point_b = (point *) b;

	if (point_a->x < point_b->x) return -1;
	if (point_a->x > point_b->x) return 1;
	return 0;
}

// Curve points
static GList *curve_points = NULL;
static point *selected_point = NULL;
static double point_size = 5.0;

enum curve_algorithm algorithm = CUBIC_SPLINE;

void curves_dialog_init_statics() {
	if (curves_adj_zoom == NULL) {
		// GtkAdjustment
		curves_adj_zoom = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "curves_adj_zoom"));
		// GtkComboBoxText
		curves_interpolation_combo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "curves_interpolation_combo"));
		// GtkDialog
		curves_dialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_dialog"));
		// GtkDrawingArea
		curves_drawingarea = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_drawingarea"));
		// GtkEntry
		curves_id_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_id_entry"));
		curves_x_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_x_entry"));
		curves_y_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_y_entry"));
		curves_clip_low = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_clip_low"));
		curves_clip_high = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_clip_high"));
		curves_seq_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "curves_seq_entry"));
		// GtkGrid
		curves_point_grid = GTK_GRID(gtk_builder_get_object(gui.builder, "curves_point_grid"));
		// GtkToggleButton
		curves_sequence_check = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_sequence_check"));
		curves_preview_check = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_preview_check"));
		curves_log_check = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_log_check"));
		// GtkToggleToolButton
		curves_red_toggle = GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object(gui.builder, "curves_red_toggle"));
		curves_green_toggle = GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object(gui.builder, "curves_green_toggle"));
		curves_blue_toggle = GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object(gui.builder, "curves_blue_toggle"));
		curves_grid_toggle = GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object(gui.builder, "curves_grid_toggle"));
		// GtkViewport
		curves_viewport = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_viewport"));
	}
}

static void set_histogram(gsl_histogram *histo, int layer);

static void clear_display_histogram() {
	if (display_histogram[0]) {
		for (int i = 0; i < fit->naxes[2]; i++) {
			gsl_histogram_free(display_histogram[i]);
			display_histogram[i] = NULL;
		}
	}
}

static void update_do_channel() {
	do_channel[0] = gtk_toggle_tool_button_get_active(curves_red_toggle);
	do_channel[1] = gtk_toggle_tool_button_get_active(curves_green_toggle);
	do_channel[2] = gtk_toggle_tool_button_get_active(curves_blue_toggle);
}

static int curves_update_preview();

static void curves_update_image() {
	set_cursor_waiting(TRUE);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = curves_update_preview;
	param->show_preview = gtk_toggle_button_get_active(curves_preview_check);
	notify_update((gpointer) param);
	set_cursor_waiting(FALSE);
}

// sets the channel names of the toggle buttons in the curves window, based on the number of color channels
void set_curves_toggles_names() {
	update_do_channel();

	if (fit->naxis == 2) {
		gtk_widget_set_tooltip_text(GTK_WIDGET(curves_red_toggle),
									_("Toggles whether to apply the curve to the monochrome channel"));
		GtkWidget *w;
		if (com.pref.gui.combo_theme == 0) {
			w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/monochrome_dark.svg");
		} else {
			w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/monochrome.svg");
		}
		gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(curves_red_toggle), w);
		gtk_widget_show(w);
		gtk_widget_set_visible(GTK_WIDGET(curves_green_toggle), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(curves_blue_toggle), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(curves_green_toggle), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(curves_blue_toggle), FALSE);
	} else {
		gtk_widget_set_tooltip_text(GTK_WIDGET(curves_red_toggle),
									_("Toggles whether to apply the curve to the red channel"));
		GtkWidget *w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/r.svg");
		gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(curves_red_toggle), w);
		gtk_widget_show(w);
		gtk_widget_set_sensitive(GTK_WIDGET(curves_green_toggle), TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(curves_blue_toggle), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(curves_green_toggle), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(curves_blue_toggle), TRUE);
	}
}

static void init_curve_points() {
	curve_points = NULL;

	// Add two initial points at [0,0] and [1,1]
	point *p = g_new(point, 1);
	p->x = 0.0;
	p->y = 0.0;
	curve_points = g_list_insert_sorted(curve_points, p, (GCompareFunc) compare_points);
	p = g_new(point, 1);
	p->x = 1.0;
	p->y = 1.0;
	curve_points = g_list_insert_sorted(curve_points, p, (GCompareFunc) compare_points);

	selected_point = (point *) curve_points->data;
}

static void curves_startup() {
	init_curve_points();

	add_roi_callback(curves_histogram_change_between_roi_and_image);
	roi_supported(TRUE);
	update_do_channel();
	copy_gfit_to_backup();
	// also get the display histogram
	compute_histo_for_gfit();
	set_curves_toggles_names();
	for (int i = 0; i < fit->naxes[2]; i++)
		display_histogram[i] = gsl_histogram_clone(com.layers_hist[i]);
}

static void curves_close(gboolean update_image_if_needed) {
	for (int i = 0; i < fit->naxes[2]; i++) {
		set_histogram(display_histogram[i], i);
		display_histogram[i] = NULL;
	}
	if (is_preview_active() && !copy_backup_to_gfit() && update_image_if_needed) {
		set_cursor_waiting(TRUE);
		notify_gfit_modified();
	}

	clear_backup();
	clear_display_histogram();
	roi_supported(FALSE);
	remove_roi_callback(curves_histogram_change_between_roi_and_image);
}

static void _update_clipped_pixels(size_t data) {
	if (fit->type == DATA_USHORT) {
		for (size_t i = 0; i < fit->naxes[2]; i++) {
			for (size_t j = 0; j < fit->naxes[0] * fit->naxes[1]; j++) {
				if (fit->pdata[i][j] <= 0) clipped[0]++;
				else if (fit->pdata[i][j] >= USHRT_MAX_DOUBLE) clipped[1]++;
			}
		}
	} else if (fit->type == DATA_FLOAT) {
		for (size_t i = 0; i < fit->naxes[2]; i++) {
			for (size_t j = 0; j < fit->naxes[0] * fit->naxes[1]; j++) {
				if (fit->fpdata[i][j] <= 0) clipped[0]++;
				else if (fit->fpdata[i][j] >= 1) clipped[1]++;
			}
		}
	}

	double tmp;
	char buffer[16];

	tmp = max((double) clipped[1] * 100.0 / (double) data, 0);
	g_snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_entry_set_text(curves_clip_high, buffer);
	tmp = max((double) clipped[0] * 100.0 / (double) data, 0);
	g_snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_entry_set_text(curves_clip_low, buffer);
}

static void curves_recompute() {
	set_cursor("progress");
	copy_backup_to_gfit();

	// Set the clipped to 0, then count the clipped pixels but set them as negative
	// later on, count the clipped pixels again but that time increase the value
	clipped[0] = 0;
	clipped[1] = 0;
	if (fit->type == DATA_USHORT) {
		for (size_t i = 0; i < fit->naxes[2]; i++) {
			for (size_t j = 0; j < fit->naxes[0] * fit->naxes[1]; j++) {
				if (fit->pdata[i][j] <= 0) clipped[0]--;
				else if (fit->pdata[i][j] >= USHRT_MAX_DOUBLE) clipped[1]--;
			}
		}
	} else if (fit->type == DATA_FLOAT) {
		for (size_t i = 0; i < fit->naxes[2]; i++) {
			for (size_t j = 0; j < fit->naxes[0] * fit->naxes[1]; j++) {
				if (fit->fpdata[i][j] <= 0) clipped[0]--;
				else if (fit->fpdata[i][j] >= 1) clipped[1]--;
			}
		}
	}

	struct curve_params params = {.points = curve_points, .algorithm = algorithm, .do_channel = {do_channel[0],
																								 do_channel[1],
																								 do_channel[2]}};
	apply_curve(fit, fit, params, TRUE);


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

	if (depth == 1) {
		size_t npixels = fit->rx * fit->ry;
		float factor = 1 / 3.f;
		if (fit->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0; i < npixels; i++) {
				fit->fdata[i] = (fit->fpdata[0][i] + fit->fpdata[1][i] + fit->fpdata[2][i]) * factor;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for simd num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0; i < npixels; i++) {
				fit->data[i] = (fit->pdata[0][i] + fit->pdata[1][i] + fit->pdata[2][i]) * factor;
			}
		}
		fits_change_depth(fit, 1);
	}
	size_t data = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	_update_clipped_pixels(data);
	notify_gfit_modified();
}

static void _initialize_clip_text() {
	gtk_entry_set_text(curves_clip_low, "0.000%");
	gtk_entry_set_text(curves_clip_high, "0.000%");
}

void _update_entry_text() {
	if (selected_point == NULL)
		return;

	gchar *buffer;

	buffer = g_strdup_printf("%.7f", selected_point->x);
	gtk_entry_set_text(curves_x_entry, buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", selected_point->y);
	gtk_entry_set_text(curves_y_entry, buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%d", g_list_index(curve_points, selected_point));
	gtk_entry_set_text(curves_id_entry, buffer);
	g_free(buffer);
}

static void reset_curve_points() {
	GList *iter;
	point *p;
	for (iter = curve_points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		g_free(p);
	}
	g_list_free(curve_points);
	curve_points = NULL;
	init_curve_points();
}

static void adjust_curves_vport_size() {
	int target_width, target_height, current_width, current_height;
	double zoom = gtk_adjustment_get_value(curves_adj_zoom);

	current_width = gtk_widget_get_allocated_width(curves_viewport);
	current_height = gtk_widget_get_allocated_height(curves_viewport);

	target_width = (int) (((double) current_width) * zoom);
	target_height = (int) (((double) current_height) * 1.0); // The vertical zoom is always 1.0
	gtk_widget_set_size_request(curves_drawingarea, target_width, target_height);
}

static void draw_curve(cairo_t *cr, int width, int height) {
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_source_rgb(cr, 0.98, 0.5, 0.45);
	cairo_set_line_width(cr, 1.0);

	if (algorithm == LINEAR) {
		GList *iter;
		point *p = NULL;
		for (iter = curve_points; iter != NULL; iter = iter->next) {
			p = (point *) iter->data;
			if (iter == curve_points) {
				cairo_move_to(cr, 0, height - (p->y * height));
			}
			cairo_line_to(cr, p->x * width, height - (p->y * height));
		}
		if (p) cairo_line_to(cr, width, height - (p->y * height));
	} else if (algorithm == CUBIC_SPLINE) {
		cubic_spline_data cspline_data;
		cubic_spline_fit(curve_points, &cspline_data);
		for (double x = 0; x <= 1; x += 1.0 / width) {
			double y = cubic_spline_interpolate(x, &cspline_data);
			if (x == 0) {
				cairo_move_to(cr, x * width, height - (y * height));
			} else {
				cairo_line_to(cr, x * width, height - (y * height));
			}
		}
	}

	cairo_stroke(cr);
}

void draw_curve_points(cairo_t *cr, int width, int height) {
	GList *iter;
	point *p;
	cairo_set_source_rgb(cr, 0.0, 1.0, 0.0);
	cairo_set_line_width(cr, 1.0);

	for (iter = curve_points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		cairo_rectangle(cr, p->x * width - point_size / 2, height - (p->y * height) - point_size / 2, point_size,
						point_size);
		cairo_stroke(cr);
	}

	cairo_rectangle(cr, selected_point->x * width - point_size / 2,
					height - (selected_point->y * height) - point_size / 2, point_size, point_size);
	cairo_fill(cr);
}

// erase image and redraw the background color and grid
void erase_curves_histogram_display(cairo_t *cr, int width, int height) {
	// clear all with background color
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);

	update_do_channel();
	if (gtk_toggle_tool_button_get_active(curves_grid_toggle))
		draw_grid(cr, width, height);
}

static void reset_cursors_and_values(gboolean full_reset) {
	_initialize_clip_text();
	gtk_adjustment_set_value(curves_adj_zoom, 1.0);
	gtk_entry_set_text(curves_seq_entry, "curve_");
	gtk_toggle_button_set_active(curves_preview_check, TRUE);
	gtk_toggle_button_set_active(curves_sequence_check, FALSE);
	gtk_toggle_tool_button_set_active(curves_grid_toggle, TRUE);
	if (full_reset) {
		gtk_toggle_tool_button_set_active(curves_red_toggle, TRUE);
		gtk_toggle_tool_button_set_active(curves_green_toggle, TRUE);
		gtk_toggle_tool_button_set_active(curves_blue_toggle, TRUE);
		gtk_combo_box_set_active(GTK_COMBO_BOX(curves_interpolation_combo), CUBIC_SPLINE);
		gtk_toggle_button_set_active(curves_log_check, com.pref.gui.display_histogram_mode == LOG_DISPLAY ? TRUE : FALSE);
	}
	reset_curve_points();
	_update_entry_text();
	update_gfit_curves_histogram_if_needed();
}

static int curves_update_preview() {
	fit = gui.roi.active ? &gui.roi.fit : &gfit;
	if (!closing)
		curves_recompute();
	return 0;
}

static gboolean is_curves_log_scale() {
	return (gtk_toggle_button_get_active(curves_log_check));
}

static void set_histogram(gsl_histogram *histo, int layer) {
	g_assert(layer >= 0 && layer < MAXVPORT);
	if (com.layers_hist[layer])
		gsl_histogram_free(com.layers_hist[layer]);
	com.layers_hist[layer] = histo;
}

/* Public functions */

/* call from main thread */
void update_gfit_curves_histogram_if_needed() {
	invalidate_gfit_histogram();
	if (gtk_widget_get_visible(curves_dialog)) {
		compute_histo_for_gfit();
		gtk_widget_queue_draw(curves_drawingarea);
	}
}

/* Callback functions */

void curves_histogram_change_between_roi_and_image() {
	// This should be called if the ROI is set, changed or cleared to ensure the
	// curves dialog continues to process the right data.
	fit = gui.roi.active ? &gui.roi.fit : &gfit;

	gui.roi.operation_supports_roi = TRUE;
	curves_update_image();
}

gboolean redraw_curves(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width = 339, height = 281;
	double zoom = 1;

	update_do_channel();
	width = gtk_widget_get_allocated_width(curves_drawingarea);
	height = gtk_widget_get_allocated_height(curves_drawingarea);
	zoom = gtk_adjustment_get_value(curves_adj_zoom);

	if (height == 1)
		return FALSE;

	erase_curves_histogram_display(cr, width, height);

	for (i = 0; i < MAXVPORT; i++)
		display_histo(display_histogram[i], cr, i, width, height, zoom, 1.0, FALSE, is_curves_log_scale());

	draw_curve(cr, width, height);
	draw_curve_points(cr, width, height);
	return FALSE;
}

void on_curves_display_toggle(GtkToggleButton *togglebutton, gpointer user_data) {
	update_do_channel();
	gtk_widget_queue_draw(curves_drawingarea);
	curves_update_image();
}

void on_curves_window_show(GtkWidget *object, gpointer user_data) {
	closing = FALSE;
	curves_startup();
	_initialize_clip_text();
	reset_cursors_and_values(TRUE);
	compute_histo_for_gfit();
}

void on_curves_close_button_clicked(GtkButton *button, gpointer user_data) {
	closing = TRUE;
	set_cursor_waiting(TRUE);
	curves_close(TRUE);
	set_cursor_waiting(FALSE);
	siril_close_dialog("curves_dialog");
}

void on_curves_reset_button_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	reset_cursors_and_values(FALSE);
	curves_close(TRUE);
	curves_startup();
	gtk_widget_queue_draw(curves_drawingarea);
	set_cursor_waiting(FALSE);
}

void on_curves_apply_button_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	siril_log_message(_("Applying %s curve transformation with %d points\n"),
					  algorithm == LINEAR ? "linear" : "cubic spline", g_list_length(curve_points));

	if (gtk_toggle_button_get_active(curves_sequence_check)
		&& sequence_is_loaded()) {
		// apply to sequence
		int algoritm = gtk_combo_box_get_active(GTK_COMBO_BOX(curves_interpolation_combo));
		struct curve_params params = {.points = curve_points, .algorithm = algoritm, .do_channel = {do_channel[0],
																									do_channel[1],
																									do_channel[2]}};
		struct curve_data *args = malloc(sizeof(struct mtf_data));

		args->params = params;
		args->seq_entry = strdup(gtk_entry_get_text(curves_seq_entry));
		args->seq = &com.seq;
		// If entry text is empty, set the sequence prefix to default 'curve_'
		if (args->seq_entry && args->seq_entry[0] == '\0') {
			args->seq_entry = strdup("stretch_");
		}

		curves_close(FALSE);
		siril_close_dialog("curves_dialog");

		gtk_toggle_button_set_active(curves_sequence_check, FALSE);
		apply_curve_to_sequence(args);
	} else {
		// the apply button resets everything after recomputing with the current values
		fit = &gfit;
		if (!gtk_toggle_button_get_active(curves_preview_check) || gui.roi.active) {
			copy_backup_to_gfit();

			curves_recompute();
		}
		populate_roi();
		undo_save_state(get_preview_gfit_backup(), "Curves Transformation with %d points", g_list_length(curve_points));

		clear_backup();
		clear_display_histogram();
		// reinit
		curves_startup();
		reset_cursors_and_values(FALSE);
		set_cursor("default");
	}
}

int curve_finalize_hook(struct generic_seq_args *args) {
	struct curve_data *curve_data = (struct curve_data *) args->user;
	int return_value = seq_finalize_hook(args);
	free(curve_data);
	return return_value;
}

static int curve_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
							rectangle *_, int threads) {
	struct curve_data *curve_args = (struct curve_data *) args->user;
	apply_curve(fit, fit, curve_args->params, FALSE);
	return 0;
}

void apply_curve_to_sequence(struct curve_data *curve_args) {
	struct generic_seq_args *args = create_default_seqargs(curve_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = curve_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = curve_finalize_hook;
	args->image_hook = curve_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Curves Transform");
	args->has_output = TRUE;
	args->new_seq_prefix = curve_args->seq_entry;
	args->load_new_sequence = TRUE;
	args->user = curve_args;
	curve_args->fit = NULL;

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(curve_args->seq_entry);
		free(curve_args);
		free_generic_seq_args(args);
	}
}

void apply_curves_cancel() {
	set_cursor_waiting(TRUE);
	curves_close(TRUE);
	set_cursor_waiting(FALSE);
}

void on_curves_reset_zoom_clicked(GtkButton *button, gpointer user_data) {
	gtk_adjustment_set_value(curves_adj_zoom, 1.0);
}

void on_curves_spin_zoom_value_changed(GtkRange *range, gpointer user_data) {
	adjust_curves_vport_size();
	gtk_widget_queue_draw(curves_drawingarea);
}

void on_curves_interpolation_combo_changed(GtkComboBox *widget, gpointer user_data) {
	algorithm = gtk_combo_box_get_active(widget);
	curves_update_image();
	gtk_widget_queue_draw(curves_drawingarea);
}

void setup_curve_dialog() {
	gtk_window_set_title(GTK_WINDOW(curves_dialog), _("Curves Transformation"));
	gtk_widget_set_visible(GTK_WIDGET(curves_point_grid), TRUE);
	gtk_window_resize(GTK_WINDOW(curves_dialog), 1, 1);
}

void toggle_curves_window_visibility() {
	if (gtk_widget_get_visible(lookup_widget("histogram_dialog")))
		siril_close_dialog("histogram_dialog");

	// Initialize UI elements
	curves_dialog_init_statics();

	for (int i = 0; i < 3; i++) {
		do_channel[i] = TRUE;
	}
	icc_auto_assign_or_convert(&gfit, ICC_ASSIGN_ON_STRETCH);

	if (gtk_widget_get_visible(curves_dialog)) {
		set_cursor_waiting(TRUE);
		curves_close(TRUE);
		set_cursor_waiting(FALSE);
		siril_close_dialog("curves_dialog");
	} else {
		reset_cursors_and_values(TRUE);
		copy_gfit_to_backup();
		setup_curve_dialog();

		if (gui.rendering_mode == LINEAR_DISPLAY)
			setup_stretch_sliders(); // In linear mode, set sliders to 0 / 65535
		siril_open_dialog("curves_dialog");
	}
}

gboolean on_curves_drawingarea_motion_notify_event(GtkWidget *widget, GdkEventMotion *event, gpointer user_data) {
	if (!is_drawingarea_pressed)
		return FALSE;

	gdouble xpos = ((GdkEventMotion *) event)->x / (gdouble) gtk_widget_get_allocated_width(curves_drawingarea);
	// y position is inverted
	gdouble ypos = 1 - ((GdkEventMotion *) event)->y / (gdouble) gtk_widget_get_allocated_height(curves_drawingarea);

	// clamp the coordinates to 0 and 1
	xpos = (xpos < 0) ? 0 : (xpos > 1) ? 1 : xpos;
	selected_point->y = (ypos < 0) ? 0 : (ypos > 1) ? 1 : ypos;

	// Find the closest points to the selected point based on x position
	// Limit the range of movement to the closest points
	GList *iter;
	point *p;
	point *prev = NULL;
	point *next = NULL;
	for (iter = curve_points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		if (p->x < selected_point->x)
			prev = p;
		else if (p->x > selected_point->x) {
			next = p;
			break;
		}
	}

	if (prev && prev->x >= xpos)
		selected_point->x = prev->x + 0.00001;
	else if (next && next->x <= xpos)
		selected_point->x = next->x - 0.00001;
	else
		selected_point->x = xpos;

	_update_entry_text();
	gtk_widget_queue_draw(widget);

	return FALSE;
}

void on_curves_drawingarea_leave_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	set_cursor("default");
}

gboolean on_curves_drawingarea_button_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	gdouble xpos = ((GdkEventButton *) event)->x / (gdouble) gtk_widget_get_allocated_width(curves_drawingarea);
	// y position is inverted
	gdouble ypos = 1 - ((GdkEventButton *) event)->y / (gdouble) gtk_widget_get_allocated_height(curves_drawingarea);

	// This value dictates how close the cursor has to be to a point to select it (makes the selection smoother)
	double click_precision = 0.02 / gtk_adjustment_get_value(curves_adj_zoom);

	if (event->button == GDK_BUTTON_PRIMARY) {
		gboolean was_point_clicked = FALSE;

		// If a point is clicked, select it
		GList *iter;
		point *p;
		for (iter = curve_points; iter != NULL; iter = iter->next) {
			p = (point *) iter->data;
			if (fabs(p->x - xpos) < click_precision && fabs(p->y - ypos) < click_precision) {
				selected_point = p;
				set_cursor("hand1");
				was_point_clicked = TRUE;
				break;
			}
		}

		// If no point is clicked, add a new point
		if (!was_point_clicked && g_list_length(curve_points) < MAX_POINTS) {
			// If a point already exists at the same x position, don't add a new point
			for (iter = curve_points; iter != NULL; iter = iter->next) {
				p = (point *) iter->data;
				if (p->x == xpos) {
					GtkWidget *popover = popover_new(curves_drawingarea, "A point already exists at this x position");
					gtk_widget_show_all(popover);
					return FALSE;
				}
			}

			point *p1 = g_new(point, 1);
			p1->x = xpos;
			p1->y = ypos;
			curve_points = g_list_insert_sorted(curve_points, p1, (GCompareFunc) compare_points);
			selected_point = p1;
		}

		is_drawingarea_pressed = TRUE;
		gtk_widget_queue_draw(widget);
		_update_entry_text();
	} else if (event->button == GDK_BUTTON_SECONDARY && !is_drawingarea_pressed) {
		// If a point is right-clicked, remove it and select the previous point. You can't remove a point if there are only two points
		GList *iter;
		point *p;
		for (iter = curve_points; iter != NULL; iter = iter->next) {
			p = (point *) iter->data;
			if (fabs(p->x - xpos) < click_precision && fabs(p->y - ypos) < click_precision) {
				if (g_list_length(curve_points) > 2) {
					GList *new_selected_point;
					if	(p == g_list_nth_data(curve_points, 0))
						new_selected_point = iter->next;
					else
						new_selected_point = iter->prev;

					curve_points = g_list_remove(curve_points, p);
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

gboolean on_curves_drawingarea_button_release_event(GtkWidget *widget,
													GdkEventButton *event, gpointer user_data) {
	is_drawingarea_pressed = FALSE;
	curves_update_image();
	return FALSE;
}

gboolean on_curves_prev_button_clicked(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	// Ignore double clicks (double click would register as 3 button presses)
	if (event->type != GDK_BUTTON_PRESS)
		return TRUE;

	int selected_point_index = g_list_index(curve_points, selected_point);
	if (selected_point_index > 0) {
		selected_point = g_list_nth_data(curve_points, selected_point_index - 1);
		gtk_widget_queue_draw(curves_drawingarea);
		_update_entry_text();
	}

	return FALSE;
}

gboolean on_curves_next_button_clicked(GtkWidget *widget, GdkEventButton *event, gpointer user_data) {
	// Ignore double clicks (double click would register as 3 button presses)
	if (event->type != GDK_BUTTON_PRESS)
		return TRUE;

	int selected_point_index = g_list_index(curve_points, selected_point);
	if (selected_point_index < g_list_length(curve_points) - 1) {
		selected_point = g_list_nth_data(curve_points, selected_point_index + 1);
		gtk_widget_queue_draw(curves_drawingarea);
		_update_entry_text();
	}

	return FALSE;
}

void on_curves_id_entry_activate(GtkEntry *entry, gpointer user_data) {
	int entered_index = g_ascii_strtod(gtk_entry_get_text(entry), NULL);
	int curve_points_length = g_list_length(curve_points);

	// If the entered index is out of bounds, select the first or last point
	if (entered_index >= 0 && entered_index < curve_points_length)
		selected_point = g_list_nth_data(curve_points, entered_index);
	else if (entered_index >= curve_points_length)
		selected_point = g_list_nth_data(curve_points, curve_points_length - 1);
	else
		selected_point = g_list_nth_data(curve_points, 0);

	gtk_widget_queue_draw(curves_drawingarea);
	_update_entry_text();
}

gboolean on_curves_id_entry_focus_out_event(GtkWidget *widget, GdkEvent *event,
											gpointer user_data) {
	on_curves_id_entry_activate(curves_id_entry, user_data);
	return FALSE;
}

void on_curves_x_entry_activate(GtkEntry *entry, gpointer user_data) {
	float curve_x = g_ascii_strtod(gtk_entry_get_text(entry), NULL);

	GList *iter;
	point *p;
	for (iter = curve_points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		if (p != selected_point && p->x == curve_x) {
			// Reset the entry text to the previous value and show a popover
			gchar *str = g_strdup_printf("%8.7f", selected_point->x);
			gtk_entry_set_text(entry, str);
			GtkWidget *popover = popover_new(curves_drawingarea, "A point already exists at this x position");
			gtk_widget_show_all(popover);
			return;
		}
	}

	if (curve_x < 0 || curve_x > 1) {
		gchar *str = g_strdup_printf("%8.7f", selected_point->x);
		gtk_entry_set_text(entry, str);
		return;
	}

	selected_point->x = curve_x;
	curve_points = g_list_sort(curve_points, (GCompareFunc) compare_points);
	gtk_widget_queue_draw(curves_drawingarea);

	curves_update_image();
	_update_entry_text();
}

gboolean on_curves_x_entry_focus_out_event(GtkWidget *widget, GdkEvent *event,
										   gpointer user_data) {
	on_curves_x_entry_activate(curves_x_entry, user_data);
	return FALSE;
}

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

gboolean on_curves_y_entry_focus_out_event(GtkWidget *widget, GdkEvent *event,
										   gpointer user_data) {
	on_curves_y_entry_activate(curves_y_entry, user_data);
	return FALSE;
}

void on_curves_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (!gtk_toggle_button_get_active(curves_preview_check)) {
		copy_backup_to_gfit();
		redraw(REMAP_ALL);
	} else {
		copy_gfit_to_backup();
		curves_update_image();
	}
}

void on_curves_log_check_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_widget_queue_draw(curves_drawingarea);
}
