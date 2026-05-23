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
#include <float.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/icc_profile.h"
#include "core/siril_log.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/siril_preview.h"
#include "core/undo.h"
#include "curves.h"
#include "histogram.h"
#include "histogram_utils.h"
#include "filters/curve_transform.h"

/* The gsl_histogram, documented here:
 * http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_21.html
 * is able to group values into bins, it does not need to handle all values
 * separately. That is useful for display purposes, but not currently used.
 */

// Cache UI elements
static GtkAdjustment *curves_adj_zoom = NULL;
static GtkDropDown *curves_interpolation_combo = NULL;
static GtkEntry *curves_id_entry = NULL, *curves_x_entry = NULL, *curves_y_entry = NULL, *curves_clip_low = NULL, *curves_clip_high = NULL, *curves_seq_entry = NULL;
static GtkGrid *curves_point_grid = NULL;
static GtkCheckButton *curves_sequence_check = NULL, *curves_preview_check = NULL, *curves_log_check = NULL;
static GtkToggleButton *curves_red_toggle = NULL, *curves_green_toggle = NULL, *curves_blue_toggle = NULL, *curves_grid_toggle = NULL;
static GtkWidget *curves_drawingarea = NULL, *curves_viewport = NULL, *curves_dialog = NULL;
gboolean redraw_curves(GtkWidget *widget, cairo_t *cr, gpointer data);
static void curves_draw_cb(GtkDrawingArea *area, cairo_t *cr,
                           int width, int height, gpointer data);

static gboolean closing = FALSE;

static gboolean do_channel[3];

static fits *fit = NULL;

static long clipped[] = {0, 0};

// The histogram displayed on the drawingarea
static gsl_histogram *display_histogram[MAXVPORT] = {NULL, NULL, NULL, NULL};

/* All references to is_drawingarea_pressed live inside #if 0 blocks
 * pending Phase 8 conversion of the legacy button-press handlers
 * to GtkGestureClick.  Mark it unused so -Werror=unused-variable doesn't
 * fire while the migration is in flight. */
static gboolean is_drawingarea_pressed = FALSE;

/* Forward decls for the GTK4 event-controller callbacks below.
 * curves_dialog_init_statics() needs to attach them. */
static void curves_da_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);
static void curves_da_released(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);
static void curves_da_motion(GtkEventControllerMotion *controller,
		double x, double y, gpointer user_data);
static void curves_prev_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);
static void curves_next_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data);

// Original ICC profile, in case we don't apply a stretch and need to revert
static cmsHPROFILE original_icc = NULL;
static gboolean single_image_stretch_applied = FALSE;

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
		// GtkDropDown
		curves_interpolation_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "curves_interpolation_combo"));
		// GtkDialog
		curves_dialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_dialog"));
		// GtkDrawingArea
		curves_drawingarea = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_drawingarea"));
		/* GTK4: GtkDrawingArea has no "draw" signal — wire the draw
		 * function explicitly.  Phase 18 stripped the .ui binding. */
		if (curves_drawingarea) {
			gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(curves_drawingarea),
			                               curves_draw_cb, NULL, NULL);

			/* GTK4 event-controllers: replace the GTK3 button-press /
			 * release / motion-notify signals on curves_drawingarea. */
			GtkGesture *click = gtk_gesture_click_new();
			gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click), 0);
			g_signal_connect(click, "pressed",  G_CALLBACK(curves_da_pressed),  NULL);
			g_signal_connect(click, "released", G_CALLBACK(curves_da_released), NULL);
			gtk_widget_add_controller(curves_drawingarea, GTK_EVENT_CONTROLLER(click));

			GtkEventController *motion = gtk_event_controller_motion_new();
			g_signal_connect(motion, "motion", G_CALLBACK(curves_da_motion), NULL);
			gtk_widget_add_controller(curves_drawingarea, motion);
		}

		/* GTK3 wired prev/next via "button-press-event" on the parent
		 * GtkBox of each arrow GtkImage (curves_prev_event_box /
		 * curves_next_event_box).  GTK4 uses GtkGestureClick on the same
		 * wrapper boxes — the .ui keeps the ids so GtkBuilder lookups
		 * still resolve. */
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
		curves_sequence_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "curves_sequence_check"));
		curves_preview_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "curves_preview_check"));
		curves_log_check = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "curves_log_check"));
		// GtkToggleButton
		curves_red_toggle = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_red_toggle"));
		curves_green_toggle = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_green_toggle"));
		curves_blue_toggle = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_blue_toggle"));
		curves_grid_toggle = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "curves_grid_toggle"));
		// GtkViewport
		curves_viewport = GTK_WIDGET(gtk_builder_get_object(gui.builder, "curves_viewport"));
	}
}



static void clear_display_histogram() {
	if (display_histogram[0]) {
		for (int i = 0; i < fit->naxes[2]; i++) {
			gsl_histogram_free(display_histogram[i]);
			display_histogram[i] = NULL;
		}
	}
}

static void update_do_channel() {
	do_channel[0] = siril_toggle_get_active(GTK_WIDGET(curves_red_toggle));
	do_channel[1] = siril_toggle_get_active(GTK_WIDGET(curves_green_toggle));
	do_channel[2] = siril_toggle_get_active(GTK_WIDGET(curves_blue_toggle));
}

static int curves_update_preview();
static gboolean curve_apply_idle(gpointer p);

static void curves_update_image() {
	set_cursor_waiting(TRUE);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = curves_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(curves_preview_check));
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
		gtk_button_set_child(GTK_BUTTON(curves_red_toggle), w);
		gtk_widget_set_visible(w, TRUE);
		gtk_widget_set_visible(GTK_WIDGET(curves_green_toggle), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(curves_blue_toggle), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(curves_green_toggle), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(curves_blue_toggle), FALSE);
	} else {
		gtk_widget_set_tooltip_text(GTK_WIDGET(curves_red_toggle),
									_("Toggles whether to apply the curve to the red channel"));
		GtkWidget *w = gtk_image_new_from_resource("/org/siril/ui/pixmaps/r.svg");
		gtk_button_set_child(GTK_BUTTON(curves_red_toggle), w);
		gtk_widget_set_visible(w, TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(curves_green_toggle), TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(curves_blue_toggle), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(curves_green_toggle), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(curves_blue_toggle), TRUE);
	}
	/* GTK4: GtkToolButton's `icon-widget` is gone; the .ui's `<child>`
	 * slot of green / blue / grid is empty.  Wire the SVG resources
	 * programmatically as the button children. */
	gtk_button_set_child(GTK_BUTTON(curves_green_toggle),
	    gtk_image_new_from_resource("/org/siril/ui/pixmaps/g.svg"));
	gtk_button_set_child(GTK_BUTTON(curves_blue_toggle),
	    gtk_image_new_from_resource("/org/siril/ui/pixmaps/b.svg"));
	gtk_button_set_child(GTK_BUTTON(curves_grid_toggle),
	    gtk_image_new_from_resource("/org/siril/ui/pixmaps/grid.svg"));
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
	compute_histo_for_fit(fit);
	set_curves_toggles_names();
	for (int i = 0; i < fit->naxes[2]; i++)
		display_histogram[i] = gsl_histogram_clone(com.layers_hist[i]);
}

static void curves_close(gboolean update_image_if_needed, gboolean revert_icc_profile) {
	for (int i = 0; i < fit->naxes[2]; i++) {
		set_histogram(display_histogram[i], i);
		display_histogram[i] = NULL;
	}
	if (is_preview_active() && !copy_backup_to_gfit() && update_image_if_needed) {
		set_cursor_waiting(TRUE);
		notify_gfit_data_modified();
		gfit_modified_update_gui();
	}

	if (revert_icc_profile && !single_image_stretch_applied) {
		if (current_icc_profile())
			cmsCloseProfile(current_icc_profile());
		gfit->icc_profile = copyICCProfile(original_icc);
		color_manage(gfit, current_icc_profile() != NULL);
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
	gtk_editable_set_text(GTK_EDITABLE(curves_clip_high), buffer);
	tmp = max((double) clipped[0] * 100.0 / (double) data, 0);
	g_snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_editable_set_text(GTK_EDITABLE(curves_clip_low), buffer);
}

/* Create and launch curve processing with generic_image_worker */
static int curves_process_with_worker(gboolean for_preview, gboolean for_roi) {
	// Allocate parameters
	struct curve_params *params = new_curve_params();
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Set up curve parameters
	params->points = curve_points;
	params->algorithm = algorithm;
	params->do_channel[0] = do_channel[0];
	params->do_channel[1] = do_channel[1];
	params->do_channel[2] = do_channel[2];
	params->fit = for_roi ? &gui.roi.fit : gfit;
	params->verbose = !for_preview;
	params->for_preview = for_preview;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_curve_params(params);
		return 1;
	}

	// Set the fit based on whether ROI is active
	args->fit = for_roi ? &gui.roi.fit : gfit;
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

static void _initialize_clip_text() {
	gtk_editable_set_text(GTK_EDITABLE(curves_clip_low), "0.000%");
	gtk_editable_set_text(GTK_EDITABLE(curves_clip_high), "0.000%");
}

void _update_entry_text() {
	if (selected_point == NULL)
		return;

	gchar *buffer;

	buffer = g_strdup_printf("%.7f", selected_point->x);
	gtk_editable_set_text(GTK_EDITABLE(curves_x_entry), buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", selected_point->y);
	gtk_editable_set_text(GTK_EDITABLE(curves_y_entry), buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%d", g_list_index(curve_points, selected_point));
	gtk_editable_set_text(GTK_EDITABLE(curves_id_entry), buffer);
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

	current_width = gtk_widget_get_width(curves_viewport);
	current_height = gtk_widget_get_height(curves_viewport);

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
	fill_histo_background(cr, width, height);

	update_do_channel();
	if (siril_toggle_get_active(GTK_WIDGET(curves_grid_toggle)))
		draw_grid(cr, width, height);
}

static void reset_cursors_and_values(gboolean full_reset) {
	_initialize_clip_text();
	gtk_adjustment_set_value(curves_adj_zoom, 1.0);
	gtk_editable_set_text(GTK_EDITABLE(curves_seq_entry), "curve_");
	siril_toggle_set_active(GTK_WIDGET(curves_preview_check), TRUE);
	siril_toggle_set_active(GTK_WIDGET(curves_sequence_check), FALSE);
	siril_toggle_set_active(GTK_WIDGET(curves_grid_toggle), TRUE);
	if (full_reset) {
		siril_toggle_set_active(GTK_WIDGET(curves_red_toggle), TRUE);
		siril_toggle_set_active(GTK_WIDGET(curves_green_toggle), TRUE);
		siril_toggle_set_active(GTK_WIDGET(curves_blue_toggle), TRUE);
		gtk_drop_down_set_selected(GTK_DROP_DOWN(curves_interpolation_combo), CUBIC_SPLINE);
		siril_toggle_set_active(GTK_WIDGET(curves_log_check), com.pref.gui.display_histogram_mode == LOG_DISPLAY ? TRUE : FALSE);
	}
	reset_curve_points();
	_update_entry_text();
	update_gfit_curves_histogram_if_needed();
}

static int curves_update_preview() {
	fit = gui.roi.active ? &gui.roi.fit : gfit;
	if (!closing) {
		copy_backup_to_gfit();

		// Reset clipped pixels counter
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

		// Process with worker
		curves_process_with_worker(TRUE, gui.roi.active);

	}
	return 0;
}

/* Idle function for preview updates */
gboolean curve_preview_idle(gpointer p) {
	// Update clipped pixels after processing completes
	size_t data = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	_update_clipped_pixels(data);

	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		gfit_modified_update_gui();
	}
	free_generic_img_args(args);
	return FALSE;
}



static gboolean curve_apply_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0) {
		size_t data = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
		_update_clipped_pixels(data);
		gfit_modified_update_gui();
		clear_backup();
		clear_display_histogram();
		curves_startup();
		reset_cursors_and_values(FALSE);
	}
	free_generic_img_args(args);
	return FALSE;
}

static gboolean is_curves_log_scale() {
	return (siril_toggle_get_active(GTK_WIDGET(curves_log_check)));
}

/* Public functions */

/* call from main thread */
void update_gfit_curves_histogram_if_needed() {
	invalidate_gfit_histogram();
	if (gtk_widget_get_visible(curves_dialog)) {
		compute_histo_for_fit(fit);   // shared version: no sat/toggle-names side effects
		gtk_widget_queue_draw(curves_drawingarea);
	}
}

/* Callback functions */

void curves_histogram_change_between_roi_and_image() {
	// This should be called if the ROI is set, changed or cleared to ensure the
	// curves dialog continues to process the right data.
	fit = gui.roi.active ? &gui.roi.fit : gfit;

	gui.roi.operation_supports_roi = TRUE;
	curves_update_image();
}

/* GTK4 draw_func adapter — forwards to the legacy redraw_curves. */
static void curves_draw_cb(GtkDrawingArea *area, cairo_t *cr,
                           int width, int height, gpointer data) {
	(void) width; (void) height;
	redraw_curves(GTK_WIDGET(area), cr, data);
}

gboolean redraw_curves(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width = 339, height = 281;
	double zoom = 1;

	update_do_channel();
	width = gtk_widget_get_width(curves_drawingarea);
	height = gtk_widget_get_height(curves_drawingarea);
	zoom = gtk_adjustment_get_value(curves_adj_zoom);

	if (height == 1)
		return FALSE;

	erase_curves_histogram_display(cr, width, height);

	gboolean is_mono = fit ? (fit->naxes[2] == 1) : FALSE;
	for (i = 0; i < MAXVPORT; i++)
		display_histo(display_histogram[i], cr, i, width, height, zoom, 1.0, FALSE, is_curves_log_scale(), is_mono);

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
	reset_cursors_and_values(FALSE);
	curves_close(TRUE, FALSE);
	curves_startup();
	gtk_widget_queue_draw(curves_drawingarea);
	set_cursor_waiting(FALSE);
}

void on_curves_apply_button_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	if (siril_toggle_get_active(GTK_WIDGET(curves_sequence_check))
		&& sequence_is_loaded()) {
		// apply to sequence
		int algoritm = gtk_drop_down_get_selected(GTK_DROP_DOWN(curves_interpolation_combo));
		struct curve_params temp_params = {.points = curve_points, .algorithm = algoritm, .do_channel = {do_channel[0],
																									do_channel[1],
																									do_channel[2]}};
		struct curve_data *args = calloc(1, sizeof(struct curve_data));

		args->params = temp_params;
		args->seq_entry = strdup(gtk_editable_get_text(GTK_EDITABLE(curves_seq_entry)));
		args->seq = &com.seq;
		// If entry text is empty, set the sequence prefix to default 'curve_'
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
			int algoritm = gtk_drop_down_get_selected(GTK_DROP_DOWN(curves_interpolation_combo));
			struct curve_params undo_params = {
				.points = curve_points,
				.algorithm = algoritm,
				.do_channel = { do_channel[0], do_channel[1], do_channel[2] }
			};
			gchar *summary = curves_log_hook(&undo_params, SUMMARY);
			fits undo_fit = { 0 };
			memcpy(&undo_fit, get_preview_gfit_backup(), sizeof(fits));
			undo_fit.icc_profile = original_icc;
			undo_fit.color_managed = original_icc != NULL;
			undo_save_state(&undo_fit, "%s", summary);
			g_free(summary);

			populate_roi();
			clear_backup();
			clear_display_histogram();
			single_image_stretch_applied = TRUE;
			curves_startup();
			reset_cursors_and_values(FALSE);
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
	adjust_curves_vport_size();
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
	gtk_window_set_title(GTK_WINDOW(curves_dialog), _("Curves Transformation"));
	gtk_widget_set_visible(GTK_WIDGET(curves_point_grid), TRUE);
	/* GTK4: gtk_window_resize removed */;
}

void toggle_curves_window_visibility() {
	if (gtk_widget_get_visible(GTK_WIDGET(gtk_builder_get_object(gui.builder, "histogram_dialog"))))
	siril_close_dialog("histogram_dialog");

	// Initialize UI elements
	curves_dialog_init_statics();

	if (gtk_widget_get_visible(curves_dialog)) {
		set_cursor_waiting(TRUE);
		curves_close(TRUE, TRUE);
		set_cursor_waiting(FALSE);
		siril_close_dialog("curves_dialog");
	} else {
		for (int i = 0; i < 3; i++) {
			do_channel[i] = TRUE;
		}
		single_image_stretch_applied = FALSE;
		// When opening the dialog with a single image loaded, we cache the original ICC
		// profile (may be NULL) in case the user closes the dialog without applying a
		// stretch, in which case we will revert.
		if (single_image_is_loaded()) {
			if (original_icc) {
				cmsCloseProfile(original_icc);
				original_icc = copyICCProfile(current_icc_profile());
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
		setup_curve_dialog();

		if (gui.rendering_mode == LINEAR_DISPLAY)
			setup_stretch_sliders(); // In linear mode, set sliders to 0 / 65535
		siril_open_dialog("curves_dialog");
	}
}

/* GTK4 motion: drag the selected control point.  Replaces the GTK3
 * "motion-notify-event"; coordinates come straight from the
 * GtkEventControllerMotion. */
static void curves_da_motion(GtkEventControllerMotion *controller,
		double x, double y, gpointer user_data) {
	(void)user_data;
	if (!is_drawingarea_pressed || !selected_point) return;
	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(controller));
	gdouble xpos = x / (gdouble) gtk_widget_get_width(curves_drawingarea);
	/* y position is inverted (top-left origin -> bottom-up curve). */
	gdouble ypos = 1 - y / (gdouble) gtk_widget_get_height(curves_drawingarea);
	xpos = (xpos < 0) ? 0 : (xpos > 1) ? 1 : xpos;
	selected_point->y = (ypos < 0) ? 0 : (ypos > 1) ? 1 : ypos;
	/* Constrain x movement to the gap between the neighbouring points. */
	point *prev = NULL, *next = NULL;
	for (GList *iter = curve_points; iter != NULL; iter = iter->next) {
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

/* GTK4 press: primary-click selects an existing point or adds a new
 * one; secondary-click removes a point.  Replaces the GTK3
 * "button-press-event". */
static void curves_da_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)n_press; (void)user_data;
	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
	gdouble xpos = x / (gdouble) gtk_widget_get_width(curves_drawingarea);
	gdouble ypos = 1 - y / (gdouble) gtk_widget_get_height(curves_drawingarea);

	/* This value dictates how close the cursor has to be to a point to
	 * select it (makes the selection smoother under zoom). */
	double click_precision = 0.02 / gtk_adjustment_get_value(curves_adj_zoom);

	guint button = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	if (button == GDK_BUTTON_PRIMARY) {
		gboolean was_point_clicked = FALSE;
		for (GList *iter = curve_points; iter != NULL; iter = iter->next) {
			point *p = (point *) iter->data;
			if (fabs(p->x - xpos) < click_precision && fabs(p->y - ypos) < click_precision) {
				selected_point = p;
				set_cursor("hand1");
				was_point_clicked = TRUE;
				break;
			}
		}
		/* No existing point under the cursor → insert a new one. */
		if (!was_point_clicked && g_list_length(curve_points) < MAX_POINTS) {
			for (GList *iter = curve_points; iter != NULL; iter = iter->next) {
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
			curve_points = g_list_insert_sorted(curve_points, p1,
			                                    (GCompareFunc) compare_points);
			selected_point = p1;
		}
		is_drawingarea_pressed = TRUE;
		gtk_widget_queue_draw(widget);
		_update_entry_text();
	} else if (button == GDK_BUTTON_SECONDARY && !is_drawingarea_pressed) {
		/* Right-click on a point removes it (if more than two remain). */
		for (GList *iter = curve_points; iter != NULL; iter = iter->next) {
			point *p = (point *) iter->data;
			if (fabs(p->x - xpos) < click_precision && fabs(p->y - ypos) < click_precision) {
				if (g_list_length(curve_points) > 2) {
					GList *new_selected_point;
					if (p == g_list_nth_data(curve_points, 0))
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
}

/* GTK4 release: clear the drag state and rebuild the previewed image. */
static void curves_da_released(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)gesture; (void)n_press; (void)x; (void)y; (void)user_data;
	is_drawingarea_pressed = FALSE;
	curves_update_image();
}

/* GTK4 prev / next: walk the selected_point through the sorted list.
 * In GTK3 the handlers were bound to the GtkImage's parent GtkEventBox
 * via "button-press-event"; in GTK4 we attach a GtkGestureClick to the
 * arrow's wrapper GtkBox (curves_prev_event_box / curves_next_event_box).
 * The double-click filter from the GTK3 body is unnecessary because
 * GtkGestureClick reports clicks via n_press — n_press > 1 is delivered
 * as a separate event, so the body just runs once per single-click. */
static void curves_prev_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)gesture; (void)x; (void)y; (void)user_data;
	if (n_press != 1) return;
	int idx = g_list_index(curve_points, selected_point);
	if (idx > 0) {
		selected_point = g_list_nth_data(curve_points, idx - 1);
		gtk_widget_queue_draw(curves_drawingarea);
		_update_entry_text();
	}
}

static void curves_next_pressed(GtkGestureClick *gesture, int n_press,
		double x, double y, gpointer user_data) {
	(void)gesture; (void)x; (void)y; (void)user_data;
	if (n_press != 1) return;
	int idx = g_list_index(curve_points, selected_point);
	if (idx < (int)g_list_length(curve_points) - 1) {
		selected_point = g_list_nth_data(curve_points, idx + 1);
		gtk_widget_queue_draw(curves_drawingarea);
		_update_entry_text();
	}
}


void on_curves_id_entry_activate(GtkEntry *entry, gpointer user_data) {
	int entered_index = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);
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

gboolean on_curves_id_entry_focus_out_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	on_curves_id_entry_activate(curves_id_entry, user_data);
	return FALSE;
}

void on_curves_x_entry_activate(GtkEntry *entry, gpointer user_data) {
	float curve_x = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(entry)), NULL);

	GList *iter;
	point *p;
	for (iter = curve_points; iter != NULL; iter = iter->next) {
		p = (point *) iter->data;
		if (p != selected_point && p->x == curve_x) {
			// Reset the entry text to the previous value and show a popover
			gchar *str = g_strdup_printf("%8.7f", selected_point->x);
			gtk_editable_set_text(GTK_EDITABLE(entry), str);
			GtkWidget *popover = popover_new(curves_drawingarea, "A point already exists at this x position");
			gtk_widget_set_visible(popover, TRUE);
			return;
		}
	}

	if (curve_x < 0 || curve_x > 1) {
		gchar *str = g_strdup_printf("%8.7f", selected_point->x);
		gtk_editable_set_text(GTK_EDITABLE(entry), str);
		return;
	}

	selected_point->x = curve_x;
	curve_points = g_list_sort(curve_points, (GCompareFunc) compare_points);
	gtk_widget_queue_draw(curves_drawingarea);

	curves_update_image();
	_update_entry_text();
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



