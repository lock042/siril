/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include "core/command.h"
#include "gui/cut.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/siril_world_cs.h"
#include "algos/background_extraction.h"
#include "algos/siril_wcs.h"
#include "algos/photometry.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "gui/open_dialog.h"
#include "gui/dialogs.h"
#include "gui/PSF_list.h"
#include "image_interactions.h"
#include "image_display.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "progress_and_log.h"
#include "message_dialog.h"
#include "registration_preview.h"

mouse_status_enum mouse_status;
cut_method cutting;

// caching widgets
static GtkWidget *rotation_dlg = NULL;
static GtkWidget *cut_dialog = NULL, *dynpsf_dlg = NULL;
static GtkLabel *label_wn1_x = NULL, *label_wn1_y = NULL, *label_wn2_x = NULL, *label_wn2_y = NULL;

static void cache_widgets() {
	if (!rotation_dlg) {
		rotation_dlg = lookup_widget("rotation_dialog");
		label_wn1_x = GTK_LABEL(lookup_widget("label_wn1_x"));
		label_wn1_y = GTK_LABEL(lookup_widget("label_wn1_y"));
		label_wn2_x = GTK_LABEL(lookup_widget("label_wn2_x"));
		label_wn2_y = GTK_LABEL(lookup_widget("label_wn2_y"));
		cut_dialog = lookup_widget("cut_dialog");
		dynpsf_dlg = lookup_widget("stars_list_window");
	}
}

/* mouse callbacks */
static double margin_size = 10;

static gboolean is_over_the_left_side_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	int s = round_to_int(margin_size / zoom);
	if (zoomed.x > com.selection.x - s && zoomed.x < com.selection.x + s) {
		if (zoomed.y > com.selection.y - s
				&& zoomed.y < com.selection.y + com.selection.h + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_right_side_of_sel(pointi zoomed, double zoom) {
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

static gboolean is_over_the_bottom_of_sel(pointi zoomed, double zoom) {
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

static gboolean is_over_the_top_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	int s = round_to_int(margin_size / zoom);
	if (zoomed.y > com.selection.y - s && zoomed.y < com.selection.y + s) {
		if (zoomed.x > com.selection.x - s
				&& zoomed.x < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_inside_of_sel(pointi zoomed, double zoom) {
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

#define MAX_CALLBACKS_PER_EVENT 10
/* selection zone event management */
static selection_update_callback _registered_selection_callbacks[MAX_CALLBACKS_PER_EVENT];
static int _nb_selection_callbacks = 0;

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
void new_selection_zone() {
	int i;
	siril_debug_print("selection: %d,%d,\t%dx%d\n", com.selection.x,
			com.selection.y, com.selection.w, com.selection.h);
	for (i = 0; i < _nb_selection_callbacks; ++i) {
		if (_registered_selection_callbacks[i])
			_registered_selection_callbacks[i]();
	}
	redraw(REDRAW_OVERLAY);
}


void delete_selected_area() {
	memset(&com.selection, 0, sizeof(rectangle));
	if (!com.script)
		new_selection_zone();
	if (gui.roi.active && com.pref.gui.roi_mode == ROI_AUTO)
		on_clear_roi();
}

void reset_display_offset() {
	siril_debug_print("resetting display offset\n");
	gui.display_offset.x = 0;
	gui.display_offset.y = 0;
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

/*
 * Button events
 */

//static void do_popup_rgbmenu(GtkWidget *my_widget, GdkEventButton *event) {
//	static GtkMenu *menu = NULL;
//
//	if (!menu) {
//		menu = GTK_MENU(gtk_builder_get_object(gui.builder, "menurgb"));
//		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
//	}
//
//#if GTK_CHECK_VERSION(3, 22, 0)
//	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
//#else
//	int button, event_time;
//
//	if (event) {
//		button = event->button;
//		event_time = event->time;
//	} else {
//		button = 0;
//		event_time = gtk_get_current_event_time();
//	}
//
//	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button,
//			event_time);
//#endif
//}

/* Gray popup menu */

static void do_popup_graymenu(GtkWidget *my_widget, GdkEventButton *event) {
	static GtkMenu *menu = NULL;

	gboolean is_a_single_image_loaded = single_image_is_loaded() && (!sequence_is_loaded()
			|| (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE
					|| com.seq.current == SCALED_IMAGE)));

	if (!menu) {
		menu = GTK_MENU(gtk_builder_get_object(gui.builder, "menugray"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
	}

	// selection submenu
	double original_ratio = (double)gfit.rx / (double)gfit.ry;
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_free")), gui.ratio == 0.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_preserve")), gui.ratio == original_ratio);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_16_9")), gui.ratio == 16.0 / 9.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_4_3")), gui.ratio == 4.0 / 3.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_3_2")), gui.ratio == 3.0 / 2.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_1_1")), gui.ratio == 1.0 / 1.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_3_4")), gui.ratio == 3.0 / 4.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_2_3")), gui.ratio == 2.0 / 3.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_9_16")), gui.ratio == 9.0 / 16.0);
	gtk_widget_set_sensitive(lookup_widget("menuitem_selection_preserve"), is_a_single_image_loaded || sequence_is_loaded());
	gtk_widget_set_sensitive(lookup_widget("menuitem_selection_all"), is_a_single_image_loaded || sequence_is_loaded());
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_guides_0")), com.pref.gui.selection_guides == 0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_guides_2")), com.pref.gui.selection_guides == 2);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_guides_3")), com.pref.gui.selection_guides == 3);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_guides_5")), com.pref.gui.selection_guides == 5);

#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#else
	int button, event_time;

	if (event) {
		button = event->button;
		event_time = event->time;
	} else {
		button = 0;
		event_time = gtk_get_current_event_time();
	}

	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button, event_time);
#endif
}

//gboolean rgb_area_popup_menu_handler(GtkWidget *widget) {
//	do_popup_rgbmenu(widget, NULL);
//	return TRUE;
//}

void init_mouse() {
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
}

static GdkModifierType get_primary() {
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
}

gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {

	cache_widgets();
	/* when double clicking on drawing area (if no images loaded)
	 * you can load an image This feature is in GIMP and I really
	 * love it: lazy world :).
	 */
	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		if (event->button == GDK_BUTTON_PRIMARY
				&& event->type == GDK_DOUBLE_BUTTON_PRESS) {
			header_open_button_clicked();
			launch_clipboard_survey();
		}
		return FALSE;
	}

	double zoom = get_zoom_val();

	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&gui.image_matrix, &evpos.x, &evpos.y);

	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
	gboolean inside = clamp2image(&zoomed);
	//siril_debug_print("clicked at %g, %g, in image it's %d, %d (pointer is%s inside)\n",
	//		event->x, event->y, zoomed.x, zoomed.y, inside ? "" : " not");

	if (inside) {
		/* if Ctrl-Shift is pressed, prepare to measure */
		if (event->button == GDK_BUTTON_PRIMARY && (event->state & GDK_SHIFT_MASK) && (event->state & get_primary())) {
			gui.measure_start.x = zoomed.x;
			gui.measure_start.y = zoomed.y;
			gui.measure_end.x = zoomed.x;
			gui.measure_end.y = zoomed.y;
		}

		/* Ctrl click to drag */
		else if (event->state & get_primary()) {
			if (event->button == GDK_BUTTON_PRIMARY) {
				// viewport translation
				gui.translating = TRUE;
				gui.start.x = (int)(event->x);
				gui.start.y = (int)(event->y);
				return TRUE;
			}
		}

		/* else, click on gray image */
		else if (event->button == GDK_BUTTON_PRIMARY) {	// left click
			point pt;
			int radius, s;
			gboolean right, left, bottom, top;
			rectangle area;
			struct phot_config *ps = NULL;
			switch (mouse_status) {
				case MOUSE_ACTION_SELECT_REG_AREA:
					if (gui.drawing) {
						gui.drawing = FALSE;
					} else {
						gui.drawing = TRUE;
						if (is_inside_of_sel(zoomed, zoom)) {
							// Move selection
							gui.freezeX = gui.freezeY = TRUE;
							gui.start = zoomed;
							gui.origin.x = com.selection.x;
							gui.origin.y = com.selection.y;
						} else {
							// Default values
							gui.freezeX = gui.freezeY = FALSE;
							// The order matters if the selection is so small that edge detection overlaps
							// and need to be the same as in the on_drawingarea_motion_notify_event()
							right = is_over_the_right_side_of_sel(zoomed, zoom);
							left = is_over_the_left_side_of_sel(zoomed, zoom);
							bottom = is_over_the_bottom_of_sel(zoomed, zoom);
							top = is_over_the_top_of_sel(zoomed, zoom);
							if (right || left || bottom || top) {
								// Freeze one axis when grabbing an edge far enough from a corner
								if (right) {
									gui.start.x = com.selection.x;
									if (!bottom && !top)
										gui.freezeY = TRUE;
								} else if (left) {
									gui.start.x = com.selection.x + com.selection.w;
									if (!bottom && !top)
										gui.freezeY = TRUE;
								}
								if (bottom) {
									gui.start.y = com.selection.y;
									if (!left && !right)
										gui.freezeX = TRUE;
								} else if (top) {
									gui.start.y = com.selection.y + com.selection.h;
									if (!left && !right)
										gui.freezeX = TRUE;
								}
							} else {
								gui.start = zoomed;
								com.selection.h = 0;
								com.selection.w = 0;
							}
						}
					}
					redraw(REDRAW_OVERLAY);
					break;
				case MOUSE_ACTION_DRAW_SAMPLES:
					radius = get_background_sample_radius();

					pt.x = (gdouble) zoomed.x;
					pt.y = (gdouble) zoomed.y;

					if (pt.x + radius < gfit.rx && pt.y + radius < gfit.ry
							&& pt.x - radius > 0 && pt.y - radius > 0) {
						com.grad_samples = add_background_sample(com.grad_samples, &gfit, pt);

						redraw(REDRAW_OVERLAY);
						redraw_previews();
					}
					break;
				case MOUSE_ACTION_PHOTOMETRY:
					s = com.pref.phot_set.outer * 1.2;
					area.x = zoomed.x - s;
					area.y = zoomed.y - s;
					area.w = s * 2;
					area.h = s * 2;
					if (zoomed.x - area.w > 0 && zoomed.x + area.w < gfit.rx
							&& zoomed.y - area.h > 0 && zoomed.y + area.h < gfit.ry) {
						ps = phot_set_adjusted_for_image(&gfit);
						gui.qphot = psf_get_minimisation(&gfit, select_vport(gui.cvport), &area, TRUE, ps, TRUE, com.pref.starfinder_conf.profile, NULL);
						free(ps);
						if (gui.qphot) {
							gui.qphot->xpos = gui.qphot->x0 + area.x;
							if (gfit.top_down)
								gui.qphot->ypos = gui.qphot->y0 + area.y;
							else
								gui.qphot->ypos = area.y + area.h - gui.qphot->y0;
							redraw(REDRAW_OVERLAY);
							popup_psf_result(gui.qphot, &area, &gfit);
						}
					}
					break;
				case MOUSE_ACTION_CUT_SELECT:
					// Reset the cut line before setting new coords in order to avoid
					// drawing artefacts
					gui.cut.cut_start.x = -1.;
					gui.cut.cut_start.y = -1.;
					gui.cut.cut_end.x = -1.;
					gui.cut.cut_end.y = -1.;
					if (event->state & GDK_SHIFT_MASK) {
						cutting = CUT_VERT_OR_HORIZ;
					} else {
						cutting = CUT_UNCONSTRAINED;
					}
					gui.cut.cut_start.x = zoomed.x;
					gui.cut.cut_start.y = zoomed.y;
					// This is a new cut line so reset any spectroscopic wavenumber points
					gui.cut.cut_wn1.x = -1.;
					gui.cut.cut_wn1.y = -1.;
					gui.cut.cut_wn2.x = -1.;
					gui.cut.cut_wn2.y = -1.;
					gtk_label_set_text(label_wn1_x, "");
					gtk_label_set_text(label_wn1_y, "");
					gtk_label_set_text(label_wn2_x, "");
					gtk_label_set_text(label_wn2_y, "");
					break;
				default:
					break;
			}
		} else if (event->button == GDK_BUTTON_SECONDARY) {	// right click
			// Reset the cut line if one has been drawn
			gui.cut.cut_start.x = -1;
			gui.cut.cut_start.y = -1;
			gui.cut.cut_end.x = -1;
			gui.cut.cut_end.y = -1;
			if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				point pt;
				int radius = (int) (25 / 2);

				pt.x = (gdouble) zoomed.x;
				pt.y = (gdouble) zoomed.y;

				if (pt.x + radius <= gfit.rx && pt.y + radius <= gfit.ry
						&& pt.x - radius >= 0 && pt.y - radius >= 0) {
					com.grad_samples = remove_background_sample(com.grad_samples, &gfit, pt);

					redraw(REDRAW_OVERLAY);
					redraw_previews();
				}
			} else if (mouse_status == MOUSE_ACTION_PHOTOMETRY) {
				if (sequence_is_loaded()) {
					int s = com.pref.phot_set.outer * 1.2;
					rectangle area = { zoomed.x - s, zoomed.y - s, s * 2, s * 2 };
					if (zoomed.x - area.w > 0 && zoomed.x + area.w < gfit.rx
							&& zoomed.y - area.h > 0 && zoomed.y + area.h < gfit.ry) {
						memcpy(&com.selection, &area, sizeof(rectangle));
						process_seq_psf(0);
						delete_selected_area();
					}
				}
			}
		}
	}
	return FALSE;
}

gboolean on_drawingarea_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {

	cache_widgets();
	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&gui.image_matrix, &evpos.x, &evpos.y);

	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
	gboolean inside = clamp2image(&zoomed);

	if (event->button == GDK_BUTTON_PRIMARY && gui.measure_start.x != -1.) {
		gui.measure_end.x = zoomed.x;
		gui.measure_end.y = zoomed.y;
		gboolean use_arcsec = (gfit.wcsdata.pltsolvd || gui.cut.pref_as);
		measure_line(&gfit, gui.measure_start, gui.measure_end, use_arcsec);
		gui.measure_start.x = -1.;
		gui.measure_start.y = -1.;
		gui.measure_end.x = -1.;
		gui.measure_end.y = -1.;
	}

	else if (event->button == GDK_BUTTON_PRIMARY) {	// left click
		if (gui.translating) {
			gui.translating = FALSE;
		} else if (gui.drawing && mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
			gui.drawing = FALSE;
			/* finalize selection rectangle coordinates */
			if (!gui.freezeX) {
				if (zoomed.x >= gui.start.x) {
					com.selection.x = gui.start.x;
					com.selection.w = zoomed.x - com.selection.x;
				} else {
					com.selection.x = zoomed.x;
					com.selection.w = gui.start.x - zoomed.x;
				}
			}
			if (!gui.freezeY) {
				if (zoomed.y >= gui.start.y) {
					com.selection.y = gui.start.y;
					com.selection.h = zoomed.y - com.selection.y;
				} else {
					com.selection.y = zoomed.y;
					com.selection.h = gui.start.y - zoomed.y;
				}
			}
			// Clicking in displayed psf selects star in list if the DynamicPSF dialog is open
			if (((com.selection.w == 0 || com.selection.w == 1) && (com.selection.h == 0 || com.selection.h == 1)) && gtk_widget_is_visible(dynpsf_dlg)) {
				set_iter_of_clicked_psf((double) com.selection.x, (double) com.selection.y);
			}

			// never let selection be null if rotation_dlg is visible
			// reinstate full image instead
			if (!gui.freezeX && com.selection.w == 0 && gtk_widget_is_visible(rotation_dlg))
				com.selection.w = gfit.rx;
			if (!gui.freezeY && com.selection.h == 0 && gtk_widget_is_visible(rotation_dlg))
				com.selection.h = gfit.ry;

			if (gui.freezeX && gui.freezeY) { // Move selection
				com.selection.x = (zoomed.x - gui.start.x) + gui.origin.x;
				com.selection.y = (zoomed.y - gui.start.y) + gui.origin.y;
			}

			// Enforce a ratio and clamp selection to the image
			enforce_ratio_and_clamp();

			/* we have a new rectangular selection zone,
			 * or an unselection (empty zone) */
			new_selection_zone();

			// Terminate any specific selection modification mode
			gui.freezeX = gui.freezeY = FALSE;
		} else if (mouse_status == MOUSE_ACTION_SELECT_PREVIEW1) {
			set_preview_area(0, zoomed.x, zoomed.y);
			mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
			// redraw to get the position of the new preview area
			redraw(REDRAW_OVERLAY);
		} else if (mouse_status == MOUSE_ACTION_SELECT_PREVIEW2) {
			set_preview_area(1, zoomed.x, zoomed.y);
			mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
			redraw(REDRAW_OVERLAY);
		} else if (mouse_status == MOUSE_ACTION_GET_COMP_CENTER_COORDINATE) {
			if (gui.comp_layer_centering) {
				gui.comp_layer_centering->center.x = zoomed.x;
				gui.comp_layer_centering->center.y = zoomed.y;
				gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(gui.comp_layer_centering->centerbutton), FALSE);
			}
			mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
			redraw(REDRAW_OVERLAY);
		} else if (mouse_status == MOUSE_ACTION_CUT_SELECT) {
			point tmp;
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
			measure_line(&gfit, gui.cut.cut_start, gui.cut.cut_end, gui.cut.pref_as);
			cutting = CUT_NOT_CUTTING;
			redraw(REDRAW_OVERLAY);
			// Deselect the Cut button once the cut is made
			if (!gtk_widget_is_visible(cut_dialog))
				siril_open_dialog("cut_dialog");
		} else if (mouse_status == MOUSE_ACTION_CUT_WN1) {
			gui.cut.cut_wn1.x = zoomed.x;
			gui.cut.cut_wn1.y = zoomed.y;
			// Snap the selected pixel to the closest point on the line
			gui.cut.cut_wn1 = closest_point_on_line(gui.cut.cut_wn1, gui.cut.cut_start, gui.cut.cut_end);
			gchar* l1x = g_strdup_printf("%d", (int) gui.cut.cut_wn1.x);
			gchar* l1y = g_strdup_printf("%d", (int) gui.cut.cut_wn1.y);
			gtk_label_set_text(label_wn1_x, l1x);
			gtk_label_set_text(label_wn1_y, l1y);
			g_free(l1x);
			g_free(l1y);
			set_cursor("default");
			mouse_status = MOUSE_ACTION_NONE;
		} else if (mouse_status == MOUSE_ACTION_CUT_WN2) {
			gui.cut.cut_wn2.x = zoomed.x;
			gui.cut.cut_wn2.y = zoomed.y;
			// Snap the selected pixel to the closest point on the line
			gui.cut.cut_wn2 = closest_point_on_line(gui.cut.cut_wn2, gui.cut.cut_start, gui.cut.cut_end);
			gchar* l2x = g_strdup_printf("%d", (int) gui.cut.cut_wn2.x);
			gchar* l2y = g_strdup_printf("%d", (int) gui.cut.cut_wn2.y);
			gtk_label_set_text(label_wn2_x, l2x);
			gtk_label_set_text(label_wn2_y, l2y);
			g_free(l2x);
			g_free(l2y);
			set_cursor("default");
			mouse_status = MOUSE_ACTION_NONE;
		}
	} else if (event->button == GDK_BUTTON_MIDDLE) {	// middle click
		if (inside) {
			double dX = 1.5 * com.pref.phot_set.outer;
			double dY = dX;
			double w = 3 * com.pref.phot_set.outer;
			double h = w;

			if ((dX <= zoomed.x) && (dY <= zoomed.y)
					&& (zoomed.x - dX + w < gfit.rx)
					&& (zoomed.y - dY + h < gfit.ry)) {

				com.selection.x = zoomed.x - dX;
				com.selection.y = zoomed.y - dY;
				com.selection.w = w;
				com.selection.h = h;

				new_selection_zone();
			}
		}

	} else if (event->button == GDK_BUTTON_SECONDARY) {	// right click
		if (mouse_status != MOUSE_ACTION_DRAW_SAMPLES && mouse_status != MOUSE_ACTION_PHOTOMETRY) {
			do_popup_graymenu(widget, NULL);
		}
	}
	update_MenuItem();

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
			if (gfit.hi >= 1000)
				val_width = 4;
			if (gfit.hi >= 10000)
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

	/* don't change cursor if thread is running */
	if (get_thread_run()) return FALSE;

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

void on_drawingarea_leave_notify_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		if (get_thread_run()) {
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
	gboolean handled = FALSE;

	if (!single_image_is_loaded() && !sequence_is_loaded())
		return FALSE;

	if (gui.icc.iso12646)
		disable_iso12646_conditions(FALSE, TRUE, TRUE);

	if (event->state & get_primary()) {
		point delta;

		switch (event->direction) {
		case GDK_SCROLL_SMOOTH:	// what's that?
			gdk_event_get_scroll_deltas((GdkEvent*) event, &delta.x, &delta.y);
			if (delta.y < 0) {
				return update_zoom(event->x, event->y, ZOOM_IN);
			}
			if (delta.y > 0) {
				return update_zoom(event->x, event->y, ZOOM_OUT);
			}
			break;
		case GDK_SCROLL_DOWN:
			return update_zoom(event->x, event->y, ZOOM_OUT);
		case GDK_SCROLL_UP:
			return update_zoom(event->x, event->y, ZOOM_IN);
		default:
			handled = FALSE;
		}
	}
	return handled;
}
