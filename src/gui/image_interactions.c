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
#include "gui/mouse_action_functions.h"
#include "image_display.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "progress_and_log.h"

//#define DEBUG_SCROLL

/* This filters the modifier keys that are checked (to avoid things like checks
 * failing because NumLock is pressed). If we want to support additional modifiers
 * in future (Alt perhaps) they must be added to eventmask. This variable should
 * not be modified except using set_mouse_event_mask() at GUI startup.
 */
static GdkModifierType eventmask;

void set_mouse_event_mask() {
	eventmask = GDK_SHIFT_MASK | get_primary();
}

mouse_status_enum mouse_status;
cut_method cutting;
static double margin_size = 10;
static release_action button_release = { 0, mouse_nullfunction };

#define MAX_CALLBACKS_PER_EVENT 10
/* selection zone event management */
static selection_update_callback _registered_selection_callbacks[MAX_CALLBACKS_PER_EVENT];
static int _nb_selection_callbacks = 0;

mouse_action* create_mouse_action(guint button, GdkEventType type, GdkModifierType state, const mouse_function_metadata *metadata) {
	mouse_action* action = (mouse_action*) malloc(sizeof(mouse_action));
	action->button = button;
	action->type = type;
	action->state = state;
	action->data = metadata;
	return action;
}

scroll_action* create_scroll_action(GdkModifierType state, const scroll_function_metadata *metadata, SirilScrollDirection direction) {
	scroll_action* action = (scroll_action*)malloc(sizeof(scroll_action));
	action->direction = direction;
	action->state = state;
	action->data = metadata;
	return action;
}

void initialize_mouse_actions() {
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_PRIMARY, GDK_BUTTON_PRESS, 0, &main_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_PRIMARY, GDK_BUTTON_PRESS, get_primary(), &drag_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_PRIMARY, GDK_BUTTON_PRESS, get_primary() | GDK_SHIFT_MASK, &measure_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_PRIMARY, GDK_DOUBLE_BUTTON_PRESS, 0, &open_if_unloaded_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_MIDDLE, GDK_BUTTON_PRESS, 0, &drag_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_MIDDLE, GDK_BUTTON_PRESS, get_primary(), &photometry_box_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_MIDDLE, GDK_DOUBLE_BUTTON_PRESS, 0, &zoom_action));
	gui.mouse_actions = g_slist_prepend(gui.mouse_actions, create_mouse_action(GDK_BUTTON_SECONDARY, GDK_BUTTON_PRESS, 0, &second_action));
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
	gui.scroll_actions = g_slist_prepend(gui.scroll_actions, create_scroll_action(GDK_SHIFT_MASK, &scroll_traverses_seq_action, MOUSE_VERTICAL_SCROLL));
	gui.scroll_actions = g_slist_reverse(gui.scroll_actions);
}

void load_or_initialize_scroll_actions() {
	if (com.pref.gui.mouse_cfg.scroll_actions_array) {
		gui.scroll_actions = scroll_actions_config_to_list(com.pref.gui.mouse_cfg.scroll_actions_array);
	} else {
		initialize_scroll_actions();
	}
}

static inline gboolean mouse_action_match_to_event(GdkEventButton *event, mouse_action *action) {
	return (action->button == event->button && action->type == event->type && action->state == (event->state & eventmask));
}

static gboolean scroll_action_match_to_event(GdkEventScroll *event, scroll_action *action) {
	point delta;
	gdk_event_get_scroll_deltas((GdkEvent*) event, &delta.x, &delta.y);
	GdkModifierType filtered_event_state = event->state & eventmask;
	SirilScrollDirection eventScrollDirection;
	if (event->direction == GDK_SCROLL_DOWN || event->direction == GDK_SCROLL_UP)
		eventScrollDirection = MOUSE_VERTICAL_SCROLL;
	else if (event->direction == GDK_SCROLL_LEFT || event->direction == GDK_SCROLL_RIGHT)
		eventScrollDirection = MOUSE_HORIZ_SCROLL;
	else eventScrollDirection = MOUSE_SMOOTH_SCROLL;

#ifdef DEBUG_SCROLL
	siril_debug_print("Event direction: %u (raw direction: %u), delta x: %f, delta y: %f, modifiers: %u\n", eventScrollDirection, event->direction, delta.x, delta.y, event->state);
	siril_debug_print("Action drection: %u, modifiers: %u\n", action->direction, action->state);
#endif

	if ((action->state == filtered_event_state) && (action->direction == eventScrollDirection)) {
		return TRUE;
	} else if (action->state == filtered_event_state && event->direction == GDK_SCROLL_SMOOTH) {
		if (action->direction == MOUSE_HORIZ_SCROLL && delta.x != 0.0) {
			return TRUE;
		} else if (action->direction == MOUSE_VERTICAL_SCROLL && delta.y != 0.0) {
			return TRUE;
		}
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
	siril_debug_print("selection: %d,%d,\t%dx%d\n", com.selection.x,
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
	siril_debug_print("resetting display offset\n");
	gui.display_offset.x = 0;
	gui.display_offset.y = 0;
}

void reset_menu_toggle_button() {
	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(lookup_widget("control_window"));
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
		gtk_toggle_tool_button_set_active(GTK_TOGGLE_TOOL_BUTTON(lookup_widget("zoom_to_fit_check_button")), TRUE);
	}
}

void update_zoom_fit_button() {
	GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("zoom_to_fit_check_button"));
	if (gtk_toggle_tool_button_get_active(button)) {
		gtk_toggle_tool_button_set_active(button, FALSE);
	}
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
	} else if (pt->x > gfit.rx) {
		pt->x = gfit.rx - 1;
	} else {
		x_inside = pt->x < gfit.rx;
	}

	gboolean y_inside = FALSE;
	if (pt->y < 0) {
		pt->y = 0;
	} else if (pt->y > gfit.ry) {
		pt->y = gfit.ry - 1;
	} else {
		y_inside = pt->y < gfit.ry;
	}
	return x_inside && y_inside;
}

/*
 * Button events
 */

void init_mouse() {
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
}

GdkModifierType get_primary() {
	return gdk_keymap_get_modifier_mask(
			gdk_keymap_get_for_display(gdk_display_get_default()),
			GDK_MODIFIER_INTENT_PRIMARY_ACCELERATOR);
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
		if (com.selection.w > gfit.rx) {
			com.selection.w = gfit.rx;
			com.selection.h = round_to_int(com.selection.w / gui.ratio);
		}
		else if (com.selection.h > gfit.ry) {
			com.selection.h = gfit.ry;
			com.selection.w = round_to_int(com.selection.h * gui.ratio);
		}
	}

	// clamp the selection inside the image (needed when enforcing a ratio or moving)
	com.selection.x = set_int_in_interval(com.selection.x, 0, gfit.rx - com.selection.w);
	com.selection.y = set_int_in_interval(com.selection.y, 0, gfit.ry - com.selection.h);

	// If the image is CFA ensure the selection is aligned to a Bayer repeat
	// This ensures CFA statistics (for Bayer patterns) will be valid
	int cfa = get_cfa_pattern_index_from_string(gfit.keywords.bayer_pattern);
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

gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {

	cache_widgets();
	if (!gui.mouse_actions)
		load_or_initialize_mouse_actions();

	double zoom = get_zoom_val();
	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&gui.image_matrix, &evpos.x, &evpos.y);
	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
	gboolean inside = clamp2image(&zoomed);
	//siril_debug_print("clicked at %g, %g, in image it's %d, %d (pointer is%s inside)\n",
	//		event->x, event->y, zoomed.x, zoomed.y, inside ? "" : " not");
	mouse_data data = {.widget = widget, .event = event, .zoom = zoom, .evpos = evpos, .zoomed = zoomed, .inside = inside, .mouse_status = &mouse_status, .cutting = &cutting };

	for (GSList *iter = gui.mouse_actions ; iter != NULL; iter = g_slist_next(iter)) {
		mouse_action* action = (mouse_action*) iter->data;
		if (mouse_action_match_to_event(event, action)) {
			action->data->function(&data);
			break;
		}
	}
	return FALSE;
}

gboolean on_drawingarea_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {

	if (button_release.button == event->button) {
		double zoom = get_zoom_val();
		cache_widgets();
		// evpos.x/evpos.y = cursor position in image coordinate
		point evpos = { event->x, event->y };
		cairo_matrix_transform_point(&gui.image_matrix, &evpos.x, &evpos.y);
		// same as evpos but rounded to integer and clamped to image bounds
		pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
		gboolean inside = clamp2image(&zoomed);
		mouse_data data = {.widget = widget, .event = event, .zoom = zoom, .evpos = evpos, .zoomed = zoomed, .inside = inside, .mouse_status = &mouse_status, .cutting = &cutting };

		button_release.function(&data);
		clear_release_callback();
		gui_function(update_MenuItem, NULL);
	}
	return FALSE;
}

gboolean on_drawingarea_motion_notify_event(GtkWidget *widget,
		GdkEventMotion *event, gpointer user_data) {
	cache_widgets();
	if ((!single_image_is_loaded() && !sequence_is_loaded())
			|| gfit.type == DATA_UNSUPPORTED) {
		return FALSE;
	}

	double zoom = get_zoom_val();

	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&gui.image_matrix, &evpos.x, &evpos.y);

	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
	gboolean inside = clamp2image(&zoomed);
	//siril_debug_print("pointer at %g, %g, in image it's %d, %d (pointer is%s inside)\n",
	//		event->x, event->y, zoomed.x, zoomed.y, inside ? "" : " not");

	const gchar *label_density_names[] = { "labeldensity_red", "labeldensity_green", "labeldensity_blue", "labeldensity_rgb" };
	const gchar *label_wcs_names[] = { "labelwcs_red", "labelwcs_green", "labelwcs_blue", "labelwcs_rgb" };
	static GtkLabel *labels_wcs[G_N_ELEMENTS(label_wcs_names)] = { 0 };
	static GtkLabel *labels_density[G_N_ELEMENTS(label_density_names)] = { 0 };
	if (!labels_wcs[0]) {
		for (int i = 0; i < G_N_ELEMENTS(label_wcs_names); i++)
			labels_wcs[i] = GTK_LABEL(lookup_widget(label_wcs_names[i]));
		for (int i = 0; i < G_N_ELEMENTS(label_density_names); i++)
			labels_density[i] = GTK_LABEL(lookup_widget(label_density_names[i]));
	}

	gboolean blank_density_cvport = TRUE, blank_wcs_cvport = TRUE;
	gtk_label_set_text(GTK_LABEL(lookup_widget("label-rgb")), "");

	if (inside) {
		if (gui.measure_start.x != -1) {
			gui.measure_end.x = zoomed.x;
			gui.measure_end.y = zoomed.y;
			redraw(REDRAW_OVERLAY);
		}
		if (gui.cvport == RGB_VPORT) {
			static gchar buffer[256] = { 0 };
			if (gfit.type == DATA_USHORT) {
				g_sprintf(buffer, "<span foreground=\"#FF0000\"><b>R=%.3lf%%</b></span>\n<span foreground=\"#00FF00\"><b>G=%.3lf%%</b></span>\n<span foreground=\"#0054FF\"><b>B=%.3lf%%</b></span>",
						gfit.pdata[RLAYER][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x] / USHRT_MAX_DOUBLE * 100.0,
						gfit.pdata[GLAYER][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x] / USHRT_MAX_DOUBLE * 100.0,
						gfit.pdata[BLAYER][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x] / USHRT_MAX_DOUBLE * 100.0);
			} else if (gfit.type == DATA_FLOAT) {
				g_sprintf(buffer, "<span foreground=\"#FF0000\"><b>R=%.3lf%%</b></span>\n<span foreground=\"#00FF00\"><b>G=%.3lf%%</b></span>\n<span foreground=\"#0054FF\"><b>B=%.3lf%%</b></span>",
						gfit.fpdata[RLAYER][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x] * 100.0,
						gfit.fpdata[GLAYER][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x] * 100.0,
						gfit.fpdata[BLAYER][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x] * 100.0);
			}
			gtk_label_set_markup(GTK_LABEL(lookup_widget("label-rgb")), buffer);
		}
		static gchar buffer[256] = { 0 };
		static gchar format[256] = { 0 };
		int coords_width = 3;
		int vport = select_vport(gui.cvport);

		if (gfit.rx >= 1000 || gfit.ry >= 1000)
			coords_width = 4;
		if (gfit.rx >= 10000 || gfit.ry >= 10000)
			coords_width = 5;
		if (gfit.type == DATA_USHORT && gfit.pdata[vport] != NULL) {
			int val_width = 3;
			char *format_base_ushort;
			if (gui.cvport < RGB_VPORT) {
				format_base_ushort = "x: %%.%dd y: %%.%dd (=%%.%dd)";
			} else {
				format_base_ushort = "x: %%.%dd y: %%.%dd";
			}
			if (gfit.keywords.hi >= 1000)
				val_width = 4;
			if (gfit.keywords.hi >= 10000)
				val_width = 5;
			g_sprintf(format, format_base_ushort,
					coords_width, coords_width, val_width);
			if (gui.cvport < RGB_VPORT) {
				g_sprintf(buffer, format, zoomed.x, zoomed.y, gfit.pdata[vport][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x]);
			} else {
				g_sprintf(buffer, format, zoomed.x, zoomed.y);
			}
		} else if (gfit.type == DATA_FLOAT && gfit.fpdata[vport] != NULL) {
			char *format_base_float;
			if (gui.cvport < RGB_VPORT) {
				format_base_float = "x: %%.%dd y: %%.%dd (=%%f)";
			} else {
				format_base_float = "x: %%.%dd y: %%.%dd";
			}
			g_sprintf(format, format_base_float, coords_width, coords_width);
			if (gui.cvport < RGB_VPORT) {
				g_sprintf(buffer, format, zoomed.x, zoomed.y, gfit.fpdata[vport][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x]);
			} else {
				g_sprintf(buffer, format, zoomed.x, zoomed.y);
			}
		}

		if (buffer[0] != '\0') {
			gtk_label_set_text(labels_density[gui.cvport], buffer);
			blank_density_cvport = FALSE;
		}

		static gchar wcs_buffer[256] = { 0 };
		if (has_wcs(&gfit)) {
			double world_x, world_y;
			pix2wcs(&gfit, (double) zoomed.x, (double) (gfit.ry - zoomed.y - 1), &world_x, &world_y);
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

					gtk_label_set_text(labels_wcs[gui.cvport], wcs_buffer);
					blank_wcs_cvport = FALSE;

					g_free(ra);
					g_free(dec);
					siril_world_cs_unref(world_cs);
				}
			}
		}
	}

	if (blank_wcs_cvport && gtk_label_get_text(labels_wcs[gui.cvport])[0] != '\0')
		gtk_label_set_text(labels_wcs[gui.cvport], " ");
	if (blank_density_cvport && gtk_label_get_text(labels_density[gui.cvport])[0] != '\0')
		gtk_label_set_text(labels_density[gui.cvport], " ");

	if (gui.translating) {
		update_zoom_fit_button();

		pointi ev = { (int)(event->x), (int)(event->y) };
		point delta = { ev.x - gui.start.x , ev.y - gui.start.y };
		gui.start = ev;
		gui.display_offset.x += delta.x;
		gui.display_offset.y += delta.y;
		adjust_vport_size_to_image();
		redraw(REDRAW_OVERLAY);
	} else if (cutting) {	// button 1 down, dragging a line for the pixel profile cut
		if (event->state & GDK_SHIFT_MASK)
			cutting = CUT_VERT_OR_HORIZ;
		else
			cutting = CUT_UNCONSTRAINED;
		pointi tmp;
		tmp.x = zoomed.x;
		tmp.y = zoomed.y;
		if (cutting == CUT_VERT_OR_HORIZ) {
			if (fabs(tmp.y - gui.cut.cut_start.y) > fabs(tmp.x - gui.cut.cut_start.x)) {
				tmp.x = gui.cut.cut_start.x;
			} else {
				tmp.y = gui.cut.cut_start.y;
			}
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

		if (gui.freezeX && gui.freezeY) { // Move selection
			com.selection.x = (zoomed.x - gui.start.x) + gui.origin.x;
			com.selection.y = (zoomed.y - gui.start.y) + gui.origin.y;
		}

		// Enforce a ratio and clamp selection to the image
		enforce_ratio_and_clamp();

		// Display the dimensions of the selection while drawing it
		update_display_selection();
		redraw(REDRAW_OVERLAY);
	}

	/* don't change cursor if thread is running or if Python
	 claims the thread */
	if (get_thread_run() || com.python_claims_thread) return FALSE;

	if (inside) {
		if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
			set_cursor("cell");
		} else {
			// The order matters if the selection is so small that edge detection overlaps
			// and need to be the same as in the on_drawingarea_button_press_event()
			gboolean right = is_over_the_right_side_of_sel(zoomed, zoom);
			gboolean left = is_over_the_left_side_of_sel(zoomed, zoom);
			gboolean bottom = is_over_the_bottom_of_sel(zoomed, zoom);
			gboolean top = is_over_the_top_of_sel(zoomed, zoom);

			if (!gui.drawing && !gui.translating) {
				if (bottom && right) {
					set_cursor("se-resize");
				} else if (top && right) {
					set_cursor("ne-resize");
				} else if (right) {
					set_cursor("e-resize");
				} else if (bottom && left) {
					set_cursor("sw-resize");
				} else if (top && left) {
					set_cursor("nw-resize");
				} else if (left) {
					set_cursor("w-resize");
				} else if (bottom) {
					set_cursor("s-resize");
				} else if (top) {
					set_cursor("n-resize");
				} else if (is_inside_of_sel(zoomed, zoom)) {
					set_cursor("all-scroll");
				} else {
					set_cursor("crosshair");
				}
			} else {
				if ((event->state & get_primary()) || (is_inside_of_sel(zoomed, zoom))) {
					set_cursor("all-scroll");
				} else {
					set_cursor("crosshair");
				}
			}
		}
	} else {
		set_cursor("default");
	}

	return FALSE;
}

void on_drawingarea_enter_notify_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		if (get_thread_run() || com.python_claims_thread) {
			set_cursor_waiting(TRUE);
		} else {
			/* trick to get default cursor */
			set_cursor_waiting(FALSE);
		}
	}
}

void on_drawingarea_leave_notify_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		if (get_thread_run() || com.python_claims_thread) {
			set_cursor_waiting(TRUE);
		} else {
			/* trick to get default cursor */
			set_cursor_waiting(FALSE);
		}
	}
}

static const gchar *label_zoom[] = { "labelzoom_red", "labelzoom_green", "labelzoom_blue", "labelzoom_rgb" };

static gboolean set_label_zoom_text_idle(gpointer p) {
	const gchar *txt = (const gchar *) p;
	static GtkLabel *labels[sizeof label_zoom] = { NULL };
	if (!labels[0]) {
		for (int i = 0; i < G_N_ELEMENTS(label_zoom); i++)
			labels[i] = GTK_LABEL(lookup_widget(label_zoom[i]));
	}
	if (gfit.naxes[2] == 3)
		for (int i = 0; i < G_N_ELEMENTS(label_zoom); i++)
			gtk_label_set_text(labels[i], txt);
	else gtk_label_set_text(labels[0], txt);
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
	gdk_threads_add_idle(set_label_zoom_text_idle, zoom_buffer);
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

	if (factor >= ZOOM_MIN && factor <= ZOOM_MAX) {
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

// scroll event is the wheel interception for zoom
gboolean on_drawingarea_scroll_event(GtkWidget *widget, GdkEventScroll *event, gpointer user_data) {
	if (!gui.scroll_actions)
		load_or_initialize_scroll_actions();
	point delta;
	gdk_event_get_scroll_deltas((GdkEvent *) event, &delta.x, &delta.y);
	GSList *iter = gui.scroll_actions;
	while (iter) {
		scroll_action* action = (scroll_action*) iter->data;
		if (scroll_action_match_to_event(event, action)) {
			scroll_data data = { .event = event, .direction = action->direction };
			action->data->function(&data);
			break;
		}
		iter = g_slist_next(iter);
	}

	return FALSE;
}
