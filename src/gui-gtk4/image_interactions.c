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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_world_cs.h"
#include "algos/demosaicing.h"
#include "algos/siril_wcs.h"
#include "algos/photometry.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "image_interactions.h"
#include "gui-gtk4/mouse_action_functions.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/histo_display.h"
#include "gui-gtk4/mpp_ap_editor.h"
#include "gui-gtk4/derotation.h"
#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "progress_and_log.h"

//#define DEBUG_SCROLL

static GtkApplicationWindow *imgint_app_win = NULL;
static GtkToggleButton *imgint_zoom_fit_btn = NULL;
static GtkLabel *imgint_label_rgb = NULL;
static GtkLabel *imgint_labels_wcs[5] = { NULL };
static GtkLabel *imgint_labels_density[5] = { NULL };
static GtkLabel *imgint_labels_zoom[5] = { NULL };

static void image_interactions_init_statics(void) {
	if (imgint_app_win) return;
	imgint_app_win = GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window"));
	imgint_zoom_fit_btn = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "zoom_to_fit_check_button"));
	imgint_label_rgb = GTK_LABEL(gtk_builder_get_object(gui.builder, "label-rgb"));
	const gchar *wcs_names[] = { "labelwcs_red", "labelwcs_green", "labelwcs_blue", "labelwcs_rgb", "labelwcs_mask" };
	const gchar *density_names[] = { "labeldensity_red", "labeldensity_green", "labeldensity_blue", "labeldensity_rgb", "labeldensity_mask" };
	const gchar *zoom_names[] = { "labelzoom_red", "labelzoom_green", "labelzoom_blue", "labelzoom_rgb", "labelzoom_mask" };
	for (int i = 0; i < 5; i++) {
		imgint_labels_wcs[i] = GTK_LABEL(gtk_builder_get_object(gui.builder, wcs_names[i]));
		imgint_labels_density[i] = GTK_LABEL(gtk_builder_get_object(gui.builder, density_names[i]));
		imgint_labels_zoom[i] = GTK_LABEL(gtk_builder_get_object(gui.builder, zoom_names[i]));
	}
}

mouse_status_enum mouse_status;
cut_method cutting;
static double margin_size = 10;
static release_action button_release = { 0, mouse_nullfunction };

mouse_status_enum get_mouse_status() {
	return mouse_status;
}

#define MAX_CALLBACKS_PER_EVENT 10
/* selection zone event management */
static selection_update_callback _registered_selection_callbacks[MAX_CALLBACKS_PER_EVENT];
static int _nb_selection_callbacks = 0;

mouse_action* create_mouse_action(guint button, SirilEventType type, SirilModifierMask state, const mouse_function_metadata *metadata) {
	mouse_action* action = (mouse_action*) malloc(sizeof(mouse_action));
	action->button = button;
	action->type = type;
	action->state = state;
	action->data = metadata;
	return action;
}

scroll_action* create_scroll_action(SirilModifierMask state, const scroll_function_metadata *metadata, SirilScrollDirection direction) {
	scroll_action* action = (scroll_action*)malloc(sizeof(scroll_action));
	action->direction = direction;
	action->state = state;
	action->data = metadata;
	return action;
}

void initialize_mouse_actions() {
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_PRIMARY, SIRIL_EVENT_BUTTON_PRESS, 0, &main_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_PRIMARY, SIRIL_EVENT_BUTTON_PRESS, SIRIL_MOD_PRIMARY, &drag_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_PRIMARY, SIRIL_EVENT_BUTTON_PRESS, SIRIL_MOD_PRIMARY | SIRIL_MOD_SHIFT, &measure_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_PRIMARY, SIRIL_EVENT_2BUTTON_PRESS, 0, &open_if_unloaded_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_MIDDLE, SIRIL_EVENT_BUTTON_PRESS, 0, &drag_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_MIDDLE, SIRIL_EVENT_BUTTON_PRESS, SIRIL_MOD_PRIMARY, &photometry_box_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_MIDDLE, SIRIL_EVENT_2BUTTON_PRESS, 0, &zoom_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_SECONDARY, SIRIL_EVENT_BUTTON_PRESS, 0, &second_action));
	gui.mouse_actions = g_slist_reverse(gui.mouse_actions);
}

void load_or_initialize_mouse_actions() {
	if (com.pref.gui.mouse_cfg.mouse_actions_array) {
		gui.mouse_actions = mouse_actions_config_to_list(com.pref.gui.mouse_cfg.mouse_actions_array);
	} else {
		initialize_mouse_actions();
	}
}

void initialize_scroll_actions() {
	gui.scroll_actions = g_slist_prepend(gui.scroll_actions, create_scroll_action(0, &scroll_zooms_action, MOUSE_VERTICAL_SCROLL));
	gui.scroll_actions = g_slist_prepend(gui.scroll_actions, create_scroll_action(SIRIL_MOD_SHIFT, &scroll_traverses_seq_action, MOUSE_VERTICAL_SCROLL));
	gui.scroll_actions = g_slist_reverse(gui.scroll_actions);
}

void load_or_initialize_scroll_actions() {
	if (com.pref.gui.mouse_cfg.scroll_actions_array) {
		gui.scroll_actions = scroll_actions_config_to_list(com.pref.gui.mouse_cfg.scroll_actions_array);
	} else {
		initialize_scroll_actions();
	}
}

/* GTK4-native helpers used by the EventController-based handlers below.
 * Callers pass the SIRIL_*-namespaced event type and modifier mask;
 * the GDK→SIRIL translation happens once at handler entry via
 * siril_event_type_from_n_press / siril_mods_from_gdk, which implicitly
 * filters NumLock/CapsLock/button-state bits — no eventmask AND needed
 * here. */
static inline gboolean mouse_action_match(guint button, SirilEventType type,
                                          SirilModifierMask state,
                                          mouse_action *action) {
	return (action->button == button
	     && action->type   == type
	     && action->state  == state);
}

static gboolean scroll_action_match(GdkScrollDirection direction_raw,
                                    SirilModifierMask state,
                                    double dx, double dy,
                                    scroll_action *action) {
	SirilScrollDirection event_dir;
	if (direction_raw == GDK_SCROLL_DOWN || direction_raw == GDK_SCROLL_UP)
		event_dir = MOUSE_VERTICAL_SCROLL;
	else if (direction_raw == GDK_SCROLL_LEFT || direction_raw == GDK_SCROLL_RIGHT)
		event_dir = MOUSE_HORIZ_SCROLL;
	else
		event_dir = MOUSE_SMOOTH_SCROLL;

#ifdef DEBUG_SCROLL
	siril_log_debug("Event direction: %u (raw direction: %u), delta x: %f, delta y: %f, modifiers: %u\n",
		event_dir, direction_raw, dx, dy, state);
	siril_log_debug("Action drection: %u, modifiers: %u\n", action->direction, action->state);
#endif

	if (action->state == state && action->direction == event_dir) {
		return TRUE;
	} else if (action->state == state && direction_raw == GDK_SCROLL_SMOOTH) {
		if (action->direction == MOUSE_HORIZ_SCROLL && dx != 0.0)
			return TRUE;
		if (action->direction == MOUSE_VERTICAL_SCROLL && dy != 0.0)
			return TRUE;
	}
	return FALSE;
}


void register_release_callback(mouse_function function, guint button) {
	button_release.function = function;
	button_release.button = button;
}

void clear_release_callback() {
	button_release.function = mouse_nullfunction;
	button_release.button = 0;
}

void register_selection_update_callback(selection_update_callback f) {
	if (_nb_selection_callbacks < MAX_CALLBACKS_PER_EVENT) {
		_registered_selection_callbacks[_nb_selection_callbacks] = f;
		_nb_selection_callbacks++;
	}
}

void unregister_selection_update_callback(const selection_update_callback f) {
	int i;
	for (i = 0; i < _nb_selection_callbacks; ++i) {
		if (_registered_selection_callbacks[i] == f) {
			_registered_selection_callbacks[i] =
				_registered_selection_callbacks[_nb_selection_callbacks];
			_registered_selection_callbacks[_nb_selection_callbacks] = NULL;
			_nb_selection_callbacks--;
			return;
		}
	}
}

// send the events
gboolean new_selection_zone(gpointer user_data) {
	int i;
	siril_log_debug("selection: %d,%d,\t%dx%d\n", com.selection.x,
			com.selection.y, com.selection.w, com.selection.h);
	for (i = 0; i < _nb_selection_callbacks; ++i) {
		if (_registered_selection_callbacks[i])
			_registered_selection_callbacks[i]();
	}
	redraw(REDRAW_OVERLAY);
	return FALSE;
}

void delete_selected_area() {
	memset(&com.selection, 0, sizeof(rectangle));
	if (!com.script)
		gui_function(new_selection_zone, NULL);
	if (gui.roi.active && com.pref.gui.roi_mode == ROI_AUTO)
		on_clear_roi();
}

void reset_display_offset() {
	siril_log_debug("resetting display offset\n");
	gui.display_offset.x = 0;
	gui.display_offset.y = 0;
}

void reset_menu_toggle_button() {
	image_interactions_init_statics();
	GtkApplicationWindow *app_win = imgint_app_win;
	GAction *action_tilt = g_action_map_lookup_action(G_ACTION_MAP(app_win), "show-tilt");
	GAction *action_disto = g_action_map_lookup_action(G_ACTION_MAP(app_win), "show-disto");

	GVariant *state = g_action_get_state(action_tilt);
	if (g_variant_get_boolean(g_action_get_state(action_tilt))) {
		g_action_change_state(action_tilt, g_variant_new_boolean(FALSE));
	}
	g_variant_unref(state);
	state = g_action_get_state(action_disto);
	if (g_variant_get_boolean(g_action_get_state(action_disto))) {
		g_action_change_state(action_disto, g_variant_new_boolean(FALSE));
	}
	g_variant_unref(state);
}

void reset_zoom_default() {
	gui.zoom_value = ZOOM_DEFAULT;
	if (gui.zoom_value == ZOOM_FIT && !com.script) {
		image_interactions_init_statics();
		/* Light the fit toggle by setting the bound action's state directly.
		 * siril_toggle_set_active() would click the action-bound button and
		 * re-enter change_zoom_fit_state(), which can clobber gui.zoom_value. */
		GAction *a = g_action_map_lookup_action(G_ACTION_MAP(imgint_app_win), "zoom-fit");
		if (a)
			g_simple_action_set_state(G_SIMPLE_ACTION(a), g_variant_new_boolean(TRUE));
	}
}

void update_zoom_fit_button() {
	image_interactions_init_statics();
	/* Drive the bound action's state rather than the button's "active":
	 * a programmatic set_active() on an action-bound GTK4 button updates only
	 * the visual and leaves the action state stale, so the next click on the
	 * fit button toggled the wrong way and did nothing (needing a second
	 * click).  change-state runs change_zoom_fit_state(), which also freezes
	 * the current zoom so panning/zooming leaves fit mode cleanly. */
	GAction *a = g_action_map_lookup_action(G_ACTION_MAP(imgint_app_win), "zoom-fit");
	if (!a)
		return;
	GVariant *st = g_action_get_state(a);
	gboolean active = st && g_variant_get_boolean(st);
	if (st)
		g_variant_unref(st);
	if (active)
		g_action_change_state(a, g_variant_new_boolean(FALSE));
}

gboolean is_over_the_left_side_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	int s = round_to_int(margin_size / zoom);
	if (zoomed.x > com.selection.x - s && zoomed.x < com.selection.x + s) {
		if (zoomed.y > com.selection.y - s
				&& zoomed.y < com.selection.y + com.selection.h + s)
			return TRUE;
	}
	return FALSE;
}

gboolean is_over_the_right_side_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	int s = round_to_int(margin_size / zoom);
	if (zoomed.x > com.selection.x + com.selection.w - s
				&& zoomed.x < com.selection.x + com.selection.w + s) {
		if (zoomed.y > com.selection.y - s
				&& zoomed.y < com.selection.y + com.selection.h + s)
			return TRUE;
	}
	return FALSE;
}

gboolean is_over_the_bottom_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	int s = round_to_int(margin_size / zoom);
	if (zoomed.y > com.selection.y + com.selection.h - s
				&& zoomed.y < com.selection.y + com.selection.h + s) {
		if (zoomed.x > com.selection.x - s
				&& zoomed.x < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

gboolean is_over_the_top_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	int s = round_to_int(margin_size / zoom);
	if (zoomed.y > com.selection.y - s && zoomed.y < com.selection.y + s) {
		if (zoomed.x > com.selection.x - s
				&& zoomed.x < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

gboolean is_inside_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	int s = round_to_int(margin_size / zoom);
	if (zoomed.x >= com.selection.x + s && zoomed.x <= com.selection.x + com.selection.w - s) {
		if (zoomed.y >= com.selection.y + s && zoomed.y <= com.selection.y + com.selection.h - s)
			return TRUE;
	}
	return FALSE;
}

/* Clamp given coordinates to image boundaries.
   Returns true if point was inside, false otherwise.
*/
static gboolean clamp2image(pointi* pt) {
	gboolean x_inside = FALSE;
	if (pt->x < 0) {
		pt->x = 0;
	} else if (pt->x > gfit->rx) {
		pt->x = gfit->rx - 1;
	} else {
		x_inside = pt->x < gfit->rx;
	}

	gboolean y_inside = FALSE;
	if (pt->y < 0) {
		pt->y = 0;
	} else if (pt->y > gfit->ry) {
		pt->y = gfit->ry - 1;
	} else {
		y_inside = pt->y < gfit->ry;
	}
	return x_inside && y_inside;
}

/*
 * Button events
 */

void init_mouse() {
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
}

void init_draw_poly() {
	mouse_status = MOUSE_ACTION_DRAW_POLY;
}

void init_add_poly_to_mask() {
	uint32_t color = 0x00FF0040;
	gui.poly_fill = FALSE;
	gui.poly_ink = uint32_to_gdk_rgba(color);
	mouse_status = MOUSE_ACTION_ADD_POLY_TO_MASK;
}

void init_clear_poly_from_mask() {
	uint32_t color = 0xFF000040;
	gui.poly_fill = FALSE;
	gui.poly_ink = uint32_to_gdk_rgba(color);
	mouse_status = MOUSE_ACTION_CLEAR_POLY_FROM_MASK;
}

GdkModifierType get_primary() {
	/* GTK4: GdkKeymap is gone.  The "primary" accelerator modifier is
	 * Cmd (META) on macOS and Ctrl elsewhere. */
#ifdef OS_OSX
	return GDK_META_MASK;
#else
	return GDK_CONTROL_MASK;
#endif
}

void enforce_ratio_and_clamp() {
	if (gui.ratio > 0.0
		&& !(gui.freezeX && gui.freezeY)) {
		// Enforce a ratio for the selection
		if (gui.freezeY) {
			com.selection.h = round_to_int(com.selection.w / gui.ratio);
		} else if (gui.freezeX) {
			com.selection.w = round_to_int(com.selection.h * gui.ratio);
		} else {
			gint delta_w = round_to_int(com.selection.h * gui.ratio) - com.selection.w;
			com.selection.w += delta_w;

			if (com.selection.x < gui.start.x) { // Changing selection from the left
				com.selection.x -= delta_w;
			}
		}

		// clamp the selection dimensions
		if (com.selection.w > gfit->rx) {
			com.selection.w = gfit->rx;
			com.selection.h = round_to_int(com.selection.w / gui.ratio);
		}
		else if (com.selection.h > gfit->ry) {
			com.selection.h = gfit->ry;
			com.selection.w = round_to_int(com.selection.h * gui.ratio);
		}
	}

	// clamp the selection inside the image (needed when enforcing a ratio or moving)
	com.selection.x = set_int_in_interval(com.selection.x, 0, gfit->rx - com.selection.w);
	com.selection.y = set_int_in_interval(com.selection.y, 0, gfit->ry - com.selection.h);

	// If the image is CFA ensure the selection is aligned to a Bayer repeat
	// This ensures CFA statistics (for Bayer patterns) will be valid
	int cfa = get_cfa_pattern_index_from_string(gfit->keywords.bayer_pattern);
	switch (cfa) {
		case BAYER_FILTER_NONE:
			break;
		case BAYER_FILTER_RGGB: // Fallthrough intentional
		case BAYER_FILTER_BGGR:
		case BAYER_FILTER_GBRG:
		case BAYER_FILTER_GRBG:
			// Bayer (2x2 guarantees at least 1 pixel per subchannel for statistics)
			if(com.selection.w < 2)
				com.selection.w = 2;
			if (com.selection.h < 2)
				com.selection.h = 2;
			break;
		default: // X-trans (3x3 guarantees at least 1 pixel per subchannel for statistics)
			if(com.selection.w < 3)
				com.selection.w = 3;
			if (com.selection.h < 3)
				com.selection.h = 3;
	}
	if (cfa != BAYER_FILTER_NONE) {
		com.selection.x &= ~1;
		com.selection.y &= ~1;
		com.selection.w &= ~1;
		com.selection.h &= ~1;
		if(com.selection.w < 2)
			com.selection.w = 2;
		if (com.selection.h < 2)
			com.selection.h = 2;
	}
}

/* GTK4-native press handler.  Bound to a GtkGestureClick "pressed" signal
 * by attach_drawingarea_event_controllers().  GDK-namespaced inputs
 * (n_press, GdkModifierType) are translated into the SIRIL_EVENT_* /
 * SIRIL_MOD_* namespaces here so the matcher and the on-disk config
 * format stay decoupled from GDK enum churn. */
static void on_drawingarea_pressed_cb(GtkGestureClick *gesture, int n_press,
                                      double x, double y, gpointer user_data) {
	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
	guint button = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	GdkModifierType gdk_state = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));
	SirilEventType type = siril_event_type_from_n_press(n_press);
	SirilModifierMask mods = siril_mods_from_gdk(gdk_state);

	cache_widgets();
	if (!gui.mouse_actions)
		load_or_initialize_mouse_actions();

	double zoom = get_zoom_val();
	point evpos = { x, y };
	cairo_matrix_transform_point(&gui.image_matrix, &evpos.x, &evpos.y);
	pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
	gboolean inside = clamp2image(&zoomed);
	/* mouse_data carries raw GDK state for the action body — handlers
	 * test it with GDK_SHIFT_MASK / get_primary() against the live event. */
	mouse_data data = {
		.widget = widget, .event = NULL,
		.button = button, .state = gdk_state, .type = GDK_BUTTON_PRESS,
		.x = x, .y = y,
		.zoom = zoom, .evpos = evpos, .zoomed = zoomed, .inside = inside,
		.mouse_status = &mouse_status, .cutting = &cutting,
	};

	for (GSList *iter = gui.mouse_actions; iter != NULL; iter = g_slist_next(iter)) {
		mouse_action* action = (mouse_action*) iter->data;
		if (mouse_action_match(button, type, mods, action)) {
			action->data->function(&data);
			break;
		}
	}
}


/* GTK4-native release handler. */
static void on_drawingarea_released_cb(GtkGestureClick *gesture, int n_press,
                                       double x, double y, gpointer user_data) {
	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
	guint button = gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(gesture));
	GdkModifierType state = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));

	if (button_release.button == button) {
		double zoom = get_zoom_val();
		cache_widgets();
		point evpos = { x, y };
		cairo_matrix_transform_point(&gui.image_matrix, &evpos.x, &evpos.y);
		pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
		gboolean inside = clamp2image(&zoomed);
		mouse_data data = {
			.widget = widget, .event = NULL,
			.button = button, .state = state, .type = GDK_BUTTON_RELEASE,
			.x = x, .y = y,
			.zoom = zoom, .evpos = evpos, .zoomed = zoomed, .inside = inside,
			.mouse_status = &mouse_status, .cutting = &cutting,
		};

		button_release.function(&data);
		clear_release_callback();
		gui_function(update_MenuItem, NULL);
	}
}


/* Per-event context shared between the hover (GtkEventControllerMotion)
 * and drag (GtkGestureDrag) handlers.  Both need the same coordinate
 * transforms; computing it once keeps the two handlers in sync without
 * duplicating the body. */
typedef struct {
	GtkWidget       *widget;
	double           x, y;          /* widget-local cursor position */
	GdkModifierType  state;
	double           zoom;
	point            evpos;         /* image-coord cursor position (float) */
	pointi           zoomed;        /* image-coord cursor position, int + clamped */
	gboolean         inside;        /* TRUE if the cursor is over the image */
	gboolean         valid;         /* FALSE if no image is loaded; caller should bail */
} drawingarea_ctx;

/* Last cursor position in drawingarea-local coords, written by the motion
 * controller and read by the scroll controller as the zoom anchor.  GTK4's
 * scroll signal does not deliver pointer x/y, so we cache it here. */
static double last_motion_widget_x = 0.0;
static double last_motion_widget_y = 0.0;

static drawingarea_ctx compute_drawingarea_ctx(GtkWidget *widget,
                                               double x, double y,
                                               GdkModifierType state) {
	drawingarea_ctx ctx = {
		.widget = widget, .x = x, .y = y, .state = state, .valid = FALSE,
	};
	cache_widgets();
	if ((!single_image_is_loaded() && !sequence_is_loaded()) ||
	    gfit->type == DATA_UNSUPPORTED)
		return ctx;
	ctx.zoom = get_zoom_val();
	ctx.evpos = (point){ x, y };
	cairo_matrix_transform_point(&gui.image_matrix, &ctx.evpos.x, &ctx.evpos.y);
	ctx.zoomed = (pointi){ (int)ctx.evpos.x, (int)ctx.evpos.y };
	ctx.inside = clamp2image(&ctx.zoomed);
	ctx.valid = TRUE;
	return ctx;
}

/* Hover-feedback work shared between the motion controller (Linux) and
 * the drag gesture (macOS, where the motion controller doesn't fire while
 * a button is held).  Reads pixel values, formats labels, and updates the
 * histogram-overlay cursor.  Pure-feedback; touches no drag state. */
static void update_pointer_feedback(const drawingarea_ctx *c) {
	if (c->inside) {
		histogram_update_cursor_value(c->zoomed.x, c->zoomed.y);
	} else {
		histogram_clear_cursor_value();
	}

	image_interactions_init_statics();
	gboolean blank_density_cvport = TRUE, blank_wcs_cvport = TRUE;
	gtk_label_set_text(imgint_label_rgb, "");

	/* The remaining feedback body reads gfit->pdata / fpdata / mask /
	 * dimensions / wcs.  Take a reader trylock so a concurrent writer
	 * (image processing) can't tear the read, and skip the feedback if
	 * we can't grab the lock right now — keeps the GUI thread responsive
	 * during long-running ops.  The histogram_update_cursor_value call
	 * above acquired-and-released its own lock, so we're not nesting. */
	gboolean feedback_locked = c->inside && g_rw_lock_reader_trylock(&gfit->rwlock);
	if (c->inside && feedback_locked) {
		const pointi zoomed = c->zoomed;
		if (gui.cvport == RGB_VPORT) {
			static gchar buffer[256] = { 0 };
			if (gfit->type == DATA_USHORT) {
				g_sprintf(buffer, "<span foreground=\"#FF0000\"><b>R=%.3lf%%</b></span>\n<span foreground=\"#00FF00\"><b>G=%.3lf%%</b></span>\n<span foreground=\"#0054FF\"><b>B=%.3lf%%</b></span>",
						gfit->pdata[RLAYER][gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x] / USHRT_MAX_DOUBLE * 100.0,
						gfit->pdata[GLAYER][gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x] / USHRT_MAX_DOUBLE * 100.0,
						gfit->pdata[BLAYER][gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x] / USHRT_MAX_DOUBLE * 100.0);
			} else if (gfit->type == DATA_FLOAT) {
				g_sprintf(buffer, "<span foreground=\"#FF0000\"><b>R=%.3lf%%</b></span>\n<span foreground=\"#00FF00\"><b>G=%.3lf%%</b></span>\n<span foreground=\"#0054FF\"><b>B=%.3lf%%</b></span>",
						gfit->fpdata[RLAYER][gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x] * 100.0,
						gfit->fpdata[GLAYER][gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x] * 100.0,
						gfit->fpdata[BLAYER][gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x] * 100.0);
			}
			gtk_label_set_markup(GTK_LABEL(imgint_label_rgb), buffer);
		}
		static gchar buffer[256] = { 0 };
		static gchar format[256] = { 0 };
		int coords_width = 3;
		int vport = select_vport(gui.cvport);

		if (gfit->rx >= 1000 || gfit->ry >= 1000)
			coords_width = 4;
		if (gfit->rx >= 10000 || gfit->ry >= 10000)
			coords_width = 5;
		if (vport < MASK_VPORT) {
			if (gfit->type == DATA_USHORT && gfit->pdata[vport] != NULL) {
				int val_width = 3;
				char *format_base_ushort;
				if (gui.cvport < RGB_VPORT) {
					format_base_ushort = "x: %%.%dd y: %%.%dd (=%%.%dd)";
				} else {
					format_base_ushort = "x: %%.%dd y: %%.%dd";
				}
				if (gfit->keywords.hi >= 1000)
					val_width = 4;
				if (gfit->keywords.hi >= 10000)
					val_width = 5;
				g_sprintf(format, format_base_ushort,
						coords_width, coords_width, val_width);
				if (gui.cvport < RGB_VPORT) {
					g_sprintf(buffer, format, zoomed.x, zoomed.y, gfit->pdata[vport][gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x]);
				} else {
					g_sprintf(buffer, format, zoomed.x, zoomed.y);
				}
			} else if (gfit->type == DATA_FLOAT && gfit->fpdata[vport] != NULL) {
				char *format_base_float;
				if (gui.cvport < RGB_VPORT) {
					format_base_float = "x: %%.%dd y: %%.%dd (=%%f)";
				} else {
					format_base_float = "x: %%.%dd y: %%.%dd";
				}
				g_sprintf(format, format_base_float, coords_width, coords_width);
				if (gui.cvport < RGB_VPORT) {
					g_sprintf(buffer, format, zoomed.x, zoomed.y, gfit->fpdata[vport][gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x]);
				} else {
					g_sprintf(buffer, format, zoomed.x, zoomed.y);
				}
			}
		} else if (vport == MASK_VPORT) {
			if (gfit->mask->bitpix < 32 && gfit->mask->data != NULL) {
				int val_width = 3;
				char *format_base_ushort;
				format_base_ushort = "x: %%.%dd y: %%.%dd (=%%.%dd)";
				if (gfit->mask->bitpix < 16)
					val_width = 3;
				else
					val_width = 5;
				g_sprintf(format, format_base_ushort,
						coords_width, coords_width, val_width);
				if (gfit->mask->bitpix < 16) {
					uint8_t* m = (uint8_t*) gfit->mask->data;
					g_sprintf(buffer, format, zoomed.x, zoomed.y, m[gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x]);
				} else {
					uint16_t* m = (uint16_t*) gfit->mask->data;
					g_sprintf(buffer, format, zoomed.x, zoomed.y, m[gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x]);
				}
			} else if (gfit->mask->bitpix == 32 && gfit->mask->data != NULL) {
				char *format_base_float;
				format_base_float = "x: %%.%dd y: %%.%dd (=%%f)";
				g_sprintf(format, format_base_float, coords_width, coords_width);
				float *m = (float*) gfit->mask->data;
				g_sprintf(buffer, format, zoomed.x, zoomed.y, m[gfit->rx * (gfit->ry - zoomed.y - 1) + zoomed.x]);
			}
		}

		if (buffer[0] != '\0') {
			gtk_label_set_text(imgint_labels_density[gui.cvport], buffer);
			blank_density_cvport = FALSE;
		}

		static gchar wcs_buffer[256] = { 0 };
		if (has_wcs(gfit)) {
			double world_x, world_y;
			pix2wcs(gfit, (double) zoomed.x, (double) (gfit->ry - zoomed.y - 1), &world_x, &world_y);
			if (world_x >= 0.0 && !isnan(world_x) && !isnan(world_y)) {
				SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(world_x, world_y);
				if (world_cs) {
					gchar *ra, *dec;
					if (com.pref.gui.show_deciasec) {
						ra = siril_world_cs_alpha_format(world_cs, "%02dh%02dm%04.1lfs");
						dec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%04.1lf\"");
					} else {
						ra = siril_world_cs_alpha_format(world_cs, "%02dh%02dm%02ds");
						dec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
					}
					g_sprintf(wcs_buffer, "α: %s δ: %s", ra, dec);

					gtk_label_set_text(imgint_labels_wcs[gui.cvport], wcs_buffer);
					blank_wcs_cvport = FALSE;

					g_free(ra);
					g_free(dec);
					siril_world_cs_unref(world_cs);
				}
			}
		}
	}

	if (feedback_locked)
		g_rw_lock_reader_unlock(&gfit->rwlock);

	if (blank_wcs_cvport && gtk_label_get_text(imgint_labels_wcs[gui.cvport])[0] != '\0')
		gtk_label_set_text(imgint_labels_wcs[gui.cvport], " ");
	if (blank_density_cvport && gtk_label_get_text(imgint_labels_density[gui.cvport])[0] != '\0')
		gtk_label_set_text(imgint_labels_density[gui.cvport], " ");

}

/* Cursor-shape update.  Runs only on hover (motion controller); during
 * a drag the cursor stays whatever was set when the drag started.
 * Skipped while a background job is running so we don't fight the
 * "progress" cursor set by set_cursor_waiting. */
static void update_drawingarea_cursor(const drawingarea_ctx *c) {
	if (processing_is_job_active() || processing_is_reserved_for_python())
		return;
	if (!c->inside) {
		set_cursor("default");
		return;
	}
	if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
		set_cursor("cell");
		return;
	}
	const pointi zoomed = c->zoomed;
	const double zoom = c->zoom;
	/* Order matters when the selection is small enough that edge
	 * detection overlaps — must match on_drawingarea_pressed_cb. */
	gboolean right  = is_over_the_right_side_of_sel(zoomed, zoom);
	gboolean left   = is_over_the_left_side_of_sel(zoomed, zoom);
	gboolean bottom = is_over_the_bottom_of_sel(zoomed, zoom);
	gboolean top    = is_over_the_top_of_sel(zoomed, zoom);

	if (!gui.drawing && !gui.translating) {
		if      (bottom && right)             set_cursor("se-resize");
		else if (top && right)                set_cursor("ne-resize");
		else if (right)                       set_cursor("e-resize");
		else if (bottom && left)              set_cursor("sw-resize");
		else if (top && left)                 set_cursor("nw-resize");
		else if (left)                        set_cursor("w-resize");
		else if (bottom)                      set_cursor("s-resize");
		else if (top)                         set_cursor("n-resize");
		else if (is_inside_of_sel(zoomed, zoom)) set_cursor("all-scroll");
		else                                  set_cursor("crosshair");
	} else {
		if ((c->state & get_primary()) || is_inside_of_sel(zoomed, zoom))
			set_cursor("all-scroll");
		else
			set_cursor("crosshair");
	}
}

/* Drag-state updates.  Driven exclusively by the GtkGestureDrag's
 * drag-update — the press handler set one of {gui.drawing_polygon,
 * gui.translating, cutting, gui.drawing, gui.measure_start.x != -1}
 * and this routine carries that drag forward.  Decoupled from the
 * hover-path motion controller so macOS works (where AppKit only
 * delivers motion-while-pressed to the drag gesture). */
static void apply_drag_motion(const drawingarea_ctx *c) {
	const pointi zoomed = c->zoomed;

	if (derotation_get_drag() > 0) {
		/* Derotation disk fit: a centre/handle was grabbed on press; track it
		 * until release. Only active during a disc drag, so a plain selection
		 * drag started elsewhere falls through to the normal handling below. */
		if (c->inside)
			derotation_drag_to(zoomed.x, zoomed.y);
		return;
	}

	if (mouse_status == MOUSE_ACTION_EDIT_APS) {
		/* AP editor: a button is held on an AP, drag it (mpp_ap_move each
		 * motion event for live tracking). The drag index was set on press
		 * in mouse_action_functions.c; if nothing is grabbed this is a no-op. */
		mpp_run_t *run = mpp_get_cached_run();
		const int drag = mpp_ap_editor_get_drag_idx();
		if (run && c->inside && drag >= 0) {
			int ap_x, ap_y;
			mpp_display_to_ap_coord(run, (int)gfit->rx, (int)gfit->ry,
			                        com.seq.current, zoomed.x, zoomed.y,
			                        &ap_x, &ap_y);
			/* One undo step per drag: coalesce all motion events for this
			 * AP. The +(1<<24) offset keeps a drag distinct from a resize
			 * of the same AP (which uses the bare AP index). */
			mpp_ap_editor_record_undo((1 << 24) + drag);
			mpp_ap_move(run, drag, ap_x, ap_y);
			redraw(REDRAW_OVERLAY);
		}
		return;
	}

	if (c->inside && gui.measure_start.x != -1) {
		gui.measure_end.x = zoomed.x;
		gui.measure_end.y = zoomed.y;
		redraw(REDRAW_OVERLAY);
	}

	if (gui.drawing_polygon) {
		point *ev = malloc(sizeof(point));
		ev->x = zoomed.x;
		ev->y = zoomed.y;
		gui.drawing_polypoints = g_slist_prepend(gui.drawing_polypoints, ev);
		redraw(REDRAW_OVERLAY);
	} else if (gui.translating) {
		update_zoom_fit_button();
		pointi ev = { (int)c->x, (int)c->y };
		point delta = { ev.x - gui.start.x, ev.y - gui.start.y };
		gui.start = ev;
		gui.display_offset.x += delta.x;
		gui.display_offset.y += delta.y;
		adjust_vport_size_to_image();
		redraw(REDRAW_OVERLAY);
	} else if (cutting) {	// button 1 down, dragging a line for the pixel profile cut
		if (c->state & GDK_SHIFT_MASK)
			cutting = CUT_VERT_OR_HORIZ;
		else
			cutting = CUT_UNCONSTRAINED;
		pointi tmp = { zoomed.x, zoomed.y };
		if (cutting == CUT_VERT_OR_HORIZ) {
			if (fabs(tmp.y - gui.cut.cut_start.y) > fabs(tmp.x - gui.cut.cut_start.x))
				tmp.x = gui.cut.cut_start.x;
			else
				tmp.y = gui.cut.cut_start.y;
		}
		gui.cut.cut_end.x = tmp.x;
		gui.cut.cut_end.y = tmp.y;
		redraw(REDRAW_OVERLAY);
	} else if (gui.drawing) {	// with button 1 down
		if (!gui.freezeX) {
			if (zoomed.x > gui.start.x) {
				com.selection.x = gui.start.x;
				com.selection.w = zoomed.x - com.selection.x;
			} else {
				com.selection.x = zoomed.x;
				com.selection.w = gui.start.x - zoomed.x;
			}
		}
		if (!gui.freezeY) {
			if (zoomed.y > gui.start.y) {
				com.selection.y = gui.start.y;
				com.selection.h = zoomed.y - com.selection.y;
			} else {
				com.selection.y = zoomed.y;
				com.selection.h = gui.start.y - zoomed.y;
			}
		}
		if (gui.freezeX && gui.freezeY) {  // Move selection
			com.selection.x = (zoomed.x - gui.start.x) + gui.origin.x;
			com.selection.y = (zoomed.y - gui.start.y) + gui.origin.y;
		}
		enforce_ratio_and_clamp();
		update_display_selection();
		redraw(REDRAW_OVERLAY);
	}
}

/* Hover-only handler.  Early-returns when any mouse button is held so
 * the drag gesture's drag-update is the single source of truth for
 * motion-while-pressed on every platform.  On macOS the motion
 * controller already wouldn't fire during press (AppKit routes
 * NSEventTypeMouseDragged exclusively to drag gestures); the filter
 * makes Linux match, which both eliminates the previously-redundant
 * feedback call during drag and keeps the cursor-shape logic from
 * fighting whatever drag-update sets. */
static void on_drawingarea_motion_cb(GtkEventControllerMotion *ctrl,
                                     double x, double y, gpointer user_data) {
	(void)user_data;
	GdkModifierType state = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(ctrl));
	if (state & (GDK_BUTTON1_MASK | GDK_BUTTON2_MASK | GDK_BUTTON3_MASK))
		return;
	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(ctrl));
	last_motion_widget_x = x;
	last_motion_widget_y = y;
	drawingarea_ctx ctx = compute_drawingarea_ctx(widget, x, y, state);
	if (!ctx.valid) return;
	update_pointer_feedback(&ctx);
	update_drawingarea_cursor(&ctx);
	if (mouse_status == MOUSE_ACTION_EDIT_APS) {
		/* AP editor hover (no button held): update the hover index so the
		 * overlay can highlight the AP under the cursor. */
		mpp_run_t *run = mpp_get_cached_run();
		if (run && ctx.inside) {
			int ap_x, ap_y;
			mpp_display_to_ap_coord(run, (int)gfit->rx, (int)gfit->ry,
			                        com.seq.current, ctx.zoomed.x, ctx.zoomed.y,
			                        &ap_x, &ap_y);
			const int new_hover = mpp_ap_hit_test(run, ap_x, ap_y);
			if (new_hover != mpp_ap_editor_get_hover_idx()) {
				mpp_ap_editor_set_hover_idx(new_hover);
				redraw(REDRAW_OVERLAY);
			}
		} else if (mpp_ap_editor_get_hover_idx() != -1) {
			/* cursor left the image — clear hover */
			mpp_ap_editor_set_hover_idx(-1);
			redraw(REDRAW_OVERLAY);
		}
	}
}

/* Drag controllers.  GtkGestureDrag.drag-update is the only event
 * source for motion-while-pressed on macOS (AppKit fires
 * NSEventTypeMouseDragged, which GDK Quartz routes to drag gestures
 * exclusively, not to GtkEventControllerMotion).  On Linux it also
 * fires, redundantly with the motion controller, but the work is
 * idempotent so the duplicate is harmless.
 *
 * drag-begin / drag-end are no-ops: the press handler (running through
 * GtkGestureClick.pressed) already set up the drag flags via the
 * mouse_action dispatch, and the release handler clears them.  We
 * could wire drag-end as a fallback in case GtkGestureClick.released
 * is suppressed by drag claiming the sequence, but that hasn't been
 * observed and would risk double-firing the release callback. */
static void on_drawingarea_drag_begin_cb(GtkGestureDrag *g,
                                         double x, double y,
                                         gpointer user_data) {
	(void)g; (void)x; (void)y; (void)user_data;
}

static void on_drawingarea_drag_update_cb(GtkGestureDrag *g,
                                          double offset_x, double offset_y,
                                          gpointer user_data) {
	(void)user_data;
	GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(g));
	GdkModifierType state = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(g));
	double sx, sy;
	gtk_gesture_drag_get_start_point(g, &sx, &sy);
	double abs_x = sx + offset_x;
	double abs_y = sy + offset_y;
	last_motion_widget_x = abs_x;
	last_motion_widget_y = abs_y;
	drawingarea_ctx ctx = compute_drawingarea_ctx(widget, abs_x, abs_y, state);
	if (!ctx.valid) return;
	update_pointer_feedback(&ctx);
	apply_drag_motion(&ctx);
}

static void on_drawingarea_drag_end_cb(GtkGestureDrag *g,
                                       double offset_x, double offset_y,
                                       gpointer user_data) {
	(void)g; (void)offset_x; (void)offset_y; (void)user_data;
}


/* GTK4-native enter / leave handlers, bound to GtkEventControllerMotion's
 * "enter" / "leave" signals.  Both signals receive (controller, x, y, ud)
 * for "enter" and (controller, ud) for "leave"; we don't need the position
 * here, only the loaded-image / job-state checks. */
static void on_drawingarea_enter_cb(GtkEventControllerMotion *ctrl,
                                    double x, double y, gpointer user_data) {
	(void) ctrl; (void) x; (void) y; (void) user_data;
	if (single_image_is_loaded() || sequence_is_loaded()) {
		if (processing_is_job_active() || processing_is_reserved_for_python()) {
			set_cursor_waiting(TRUE);
		} else {
			/* trick to get default cursor */
			set_cursor_waiting(FALSE);
		}
	}
}

static void on_drawingarea_leave_cb(GtkEventControllerMotion *ctrl,
                                    gpointer user_data) {
	(void) ctrl; (void) user_data;
	if (single_image_is_loaded() || sequence_is_loaded()) {
		if (processing_is_job_active() || processing_is_reserved_for_python()) {
			set_cursor_waiting(TRUE);
		} else {
			/* trick to get default cursor */
			set_cursor_waiting(FALSE);
		}
	}
}

static gboolean set_label_zoom_text_idle(gpointer p) {
	const gchar *txt = (const gchar *) p;
	image_interactions_init_statics();
	if (gfit->naxes[2] == 3)
		for (int i = 0; i < G_N_ELEMENTS(imgint_labels_zoom); i++)
			gtk_label_set_text(imgint_labels_zoom[i], txt);
	else gtk_label_set_text(imgint_labels_zoom[0], txt);
	return FALSE;
}

void update_zoom_label() {
	static gchar zoom_buffer[256] = { 0 };
	double zoom = gui.zoom_value;
	if ((single_image_is_loaded() || sequence_is_loaded())) {
		if (zoom < 0)
			zoom = get_zoom_val();
		g_sprintf(zoom_buffer, "%d%%", (int) (zoom * 100.0));
	} else {
		g_sprintf(zoom_buffer, " ");
	}
	g_idle_add(set_label_zoom_text_idle, zoom_buffer);
}

gboolean update_zoom_label_idle(gpointer user_data) {
	update_zoom_label();
	return FALSE;
}

gboolean update_zoom(gdouble x, gdouble y, double scale) {
	// event position in image coordinates before changing the zoom value
	point evpos = { x, y };
	cairo_matrix_transform_point(&gui.image_matrix, &evpos.x, &evpos.y);
	gdouble factor;
	gboolean zoomed = FALSE;

	update_zoom_fit_button();

	gui.zoom_value = get_zoom_val();

	factor = gui.zoom_value * scale;

	/* For very large images the fit-to-window zoom can fall below ZOOM_MIN.
	 * When zooming out, lower the floor to the actual fit zoom so the full
	 * image always remains reachable via the scroll wheel.
	 * Also clamp factor to min_zoom instead of rejecting: smooth-scroll steps
	 * are variable-sized, so the user can end up between two discrete levels
	 * with no way to land exactly on min_zoom via multiplication alone. */
	double min_zoom = ZOOM_MIN;
	if (scale < 1.0 && gfit->rx > 0 && gfit->ry > 0) {
		int ww = gtk_widget_get_width(gui.view[RED_VPORT].drawarea);
		int wh = gtk_widget_get_height(gui.view[RED_VPORT].drawarea);
		if (ww > 1 && wh > 1)
			min_zoom = min(ZOOM_MIN, min((double)ww / gfit->rx, (double)wh / gfit->ry));
		if (factor < min_zoom)
			factor = min_zoom;
	}

	if (factor >= min_zoom && factor <= ZOOM_MAX) {
		zoomed = TRUE;
		gui.zoom_value = factor;
		update_zoom_label();
		adjust_vport_size_to_image();
		cairo_matrix_transform_point(&gui.display_matrix, &evpos.x, &evpos.y);
		gui.display_offset.x += x - evpos.x;
		gui.display_offset.y += y - evpos.y;
		adjust_vport_size_to_image();
		redraw(REDRAW_IMAGE);
	}
	return zoomed;
}

/* GTK4-native scroll handler.  Bound to GtkEventControllerScroll "scroll".
 * GTK4 always emits dx/dy with the signal; the discrete direction (UP/DOWN/
 * LEFT/RIGHT/SMOOTH) lives on the underlying GdkEvent and we fetch it via
 * the controller's current-event accessor. */
static gboolean on_drawingarea_scroll_cb(GtkEventControllerScroll *ctrl,
                                         double dx, double dy, gpointer user_data) {
	/* AP editor: with an AP selected, the scroll wheel resizes it (grow on
	 * scroll-up, shrink on scroll-down) instead of zooming. Gate on the
	 * editor being open and an AP selected (same as the +/- key path). */
	if (mpp_ap_editor_is_open()) {
		const int sel = mpp_ap_editor_get_selected_idx();
		mpp_run_t *run = mpp_get_cached_run();
		if (sel >= 0 && run && run->aps && sel < run->aps->count) {
			GdkEvent *sevt = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(ctrl));
			GdkScrollDirection sdir = sevt ? gdk_scroll_event_get_direction(sevt) : GDK_SCROLL_SMOOTH;
			int step = 0;
			if (sdir == GDK_SCROLL_UP)        step = 2;
			else if (sdir == GDK_SCROLL_DOWN) step = -2;
			else if (sdir == GDK_SCROLL_SMOOTH) {
				if (dy < 0.0)      step = 2;
				else if (dy > 0.0) step = -2;
			}
			if (step != 0) {
				mpp_ap_editor_record_undo(sel);   /* coalesce the resize burst */
				mpp_ap_resize(run, sel, step);
				redraw(REDRAW_OVERLAY);
				return TRUE;   /* consume — no zoom while resizing an AP */
			}
		}
	}

	if (!gui.scroll_actions)
		load_or_initialize_scroll_actions();

	GdkEvent *evt = gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(ctrl));
	GdkScrollDirection direction_raw = evt ? gdk_scroll_event_get_direction(evt) : GDK_SCROLL_SMOOTH;
	GdkModifierType gdk_state = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(ctrl));
	SirilModifierMask mods = siril_mods_from_gdk(gdk_state);

	/* Zoom anchor: GTK4's scroll signal does not deliver pointer x/y, and
	 * gdk_event_get_position() returns surface-relative coordinates that
	 * are not what update_zoom() expects.  Use the last cursor position
	 * recorded by the motion controller — for a scroll to fire over our
	 * drawingarea the pointer must have crossed it (firing motion) first,
	 * so the cached value is always current. */
	double x_pos = last_motion_widget_x;
	double y_pos = last_motion_widget_y;

	for (GSList *iter = gui.scroll_actions; iter; iter = g_slist_next(iter)) {
		scroll_action *action = (scroll_action*) iter->data;
		if (scroll_action_match(direction_raw, mods, dx, dy, action)) {
			scroll_data data = {
				.event = evt,
				.direction_raw = direction_raw,
				.state = gdk_state,
				.x = x_pos, .y = y_pos,
				.dx = dx, .dy = dy,
				.direction = action->direction,
			};
			action->data->function(&data);
			break;
		}
	}
	return TRUE;
}

/* Public entry point: attach the GTK4 event controllers (click, drag,
 * motion, enter/leave, scroll) onto a drawingarea so the legacy
 * mouse_action / scroll_action dispatchers run again under GTK4. */
void attach_drawingarea_event_controllers(GtkWidget *area) {
	if (!area) return;

	/* Click — one controller, two signals.  Set button to 0 so it captures
	 * primary, secondary and middle clicks (GtkGestureSingle defaults to
	 * primary only). */
	GtkGesture *click = gtk_gesture_click_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click), 0);
	g_signal_connect(click, "pressed",  G_CALLBACK(on_drawingarea_pressed_cb),  NULL);
	g_signal_connect(click, "released", G_CALLBACK(on_drawingarea_released_cb), NULL);
	gtk_widget_add_controller(area, GTK_EVENT_CONTROLLER(click));

	/* Drag — fires drag-update for every motion event while a button
	 * is held.  This is the single source of truth for motion-while-
	 * pressed on every platform; the motion controller below filters
	 * itself out when a button mask is set so the two paths don't
	 * both run on Linux.  On macOS the motion controller already
	 * wouldn't have fired during press (AppKit routes
	 * NSEventTypeMouseDragged exclusively to drag gestures), so the
	 * filter is a no-op there but makes the two platforms behave
	 * identically. */
	GtkGesture *drag = gtk_gesture_drag_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag), 0);
	g_signal_connect(drag, "drag-begin",  G_CALLBACK(on_drawingarea_drag_begin_cb),  NULL);
	g_signal_connect(drag, "drag-update", G_CALLBACK(on_drawingarea_drag_update_cb), NULL);
	g_signal_connect(drag, "drag-end",    G_CALLBACK(on_drawingarea_drag_end_cb),    NULL);
	gtk_widget_add_controller(area, GTK_EVENT_CONTROLLER(drag));

	/* Group click + drag so they share sequence state.  Without this,
	 * whichever gesture claims the button sequence first puts the other
	 * into the DENIED state for the rest of the sequence — and DENIED
	 * gestures stop processing subsequent events.  On Linux GtkGestureClick
	 * defers its claim long enough that the drag fires anyway, but on
	 * macOS the Quartz backend pushes click into CLAIMED earlier and
	 * drag-update never fires (selection rectangle, pan and measure all
	 * stop responding).  Grouping keeps both controllers alive on the
	 * same sequence.  Same fix the histogram overlay needed; see
	 * histo_display.c:1392. */
	gtk_gesture_group(click, drag);

	/* Motion + enter/leave.  Motion handles hover (label updates, cursor
	 * shape).  Drag updates moved to the drag controller. */
	GtkEventController *motion = gtk_event_controller_motion_new();
	g_signal_connect(motion, "motion", G_CALLBACK(on_drawingarea_motion_cb), NULL);
	g_signal_connect(motion, "enter",  G_CALLBACK(on_drawingarea_enter_cb),  NULL);
	g_signal_connect(motion, "leave",  G_CALLBACK(on_drawingarea_leave_cb),  NULL);
	gtk_widget_add_controller(area, motion);

	/* Scroll — both axes; the action dispatcher decides which to honour. */
	GtkEventController *scroll = gtk_event_controller_scroll_new(
		GTK_EVENT_CONTROLLER_SCROLL_BOTH_AXES);
	g_signal_connect(scroll, "scroll", G_CALLBACK(on_drawingarea_scroll_cb), NULL);
	gtk_widget_add_controller(area, scroll);
}

