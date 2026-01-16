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
#include "core/command_line_processor.h"
#include "core/masks.h"
#include "core/siril_log.h"
#include "gui/cut.h"
#include "gui/icc_profile.h"
#include "core/processing.h"
#include "core/undo.h"
#include "algos/background_extraction.h"
#include "algos/photometry.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "gui/open_dialog.h"
#include "gui/dialogs.h"
#include "gui/PSF_list.h"
#include "image_interactions.h"
#include "gui/mouse_action_functions.h"
#include "gui/masks_gui.h"
#include "image_display.h"
#include "gui/callbacks.h"
#include "gui/save_dialog.h"
#include "gui/utils.h"
#include "progress_and_log.h"
#include "registration_preview.h"
#include "gui/user_polygons.h"

// Tracks whether double middle click will zoom to fit or zoom to 1:1, for the
// toggle
static gboolean zoom_action_zooms_to_fit = FALSE;

// caching widgets
static GtkWidget *rotation_dlg = NULL;
static GtkWidget *cut_dialog = NULL, *dynpsf_dlg = NULL;
static GtkLabel *label_wn1_x = NULL, *label_wn1_y = NULL, *label_wn2_x = NULL, *label_wn2_y = NULL;

/* Mouse release functions should be coded as static gboolean function(mouse_data *)
 * and only be referenced from the button press functions further down this file: if
 * a button press function needs to set a mouse release action to be handled by the
 * button release event handler, it must do so using register_release_callback()
 * Mouse release functions can call other functions, some of which are declared as
 * static functions at the top of the file.
 */

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
	double original_ratio = (double)gfit->rx / (double)gfit->ry;
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

/* Mask popup menu */
static void do_popup_maskmenu(GtkWidget *my_widget, GdkEventButton *event) {
	static GtkMenu *menu = NULL;

	if (!menu) {
		menu = GTK_MENU(gtk_builder_get_object(gui.builder, "menumask_rmb"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
	}

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

void cache_widgets() {
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

// ######### Define functions that can be assigned to mouse buttons #########

/* Notes for developers *****************************************************
 * Mouse actions can be assigned to button functions by assigning to the
 * relevant function pointer in image_interactions.c. If the action requires
 * additional things to happen when the mouse button is released, it must
 * set the callback function using register_release_callback()
 * This callback is automatically removed on completion, so there is no need
 * for each release callback function to do so.
 */

// Just returns immediately. Since each mouse action calls its function pointer,
// this is required in order to avoid calling a NULL pointer if the mouse
// action does not have any other function assigned.
gboolean mouse_nullfunction(mouse_data *data) {
	return TRUE;
}

static gboolean select_reg_area_release(mouse_data *data) {
	if (gui.drawing) {
		gui.drawing = FALSE;
		/* finalize selection rectangle coordinates */
		if (!gui.freezeX) {
			if (data->zoomed.x >= gui.start.x) {
				com.selection.x = gui.start.x;
				com.selection.w = data->zoomed.x - com.selection.x;
			} else {
				com.selection.x = data->zoomed.x;
				com.selection.w = gui.start.x - data->zoomed.x;
			}
		}
		if (!gui.freezeY) {
			if (data->zoomed.y >= gui.start.y) {
				com.selection.y = gui.start.y;
				com.selection.h = data->zoomed.y - com.selection.y;
			} else {
				com.selection.y = data->zoomed.y;
				com.selection.h = gui.start.y - data->zoomed.y;
			}
		}
		// Clicking in displayed psf selects star in list if the DynamicPSF dialog is open
		if (((com.selection.w == 0 || com.selection.w == 1) && (com.selection.h == 0 || com.selection.h == 1)) && gtk_widget_is_visible(dynpsf_dlg)) {
			set_iter_of_clicked_psf((double) com.selection.x, (double) com.selection.y);
		}

		// never let selection be null if rotation_dlg is visible
		// reinstate full image instead
		if (!gui.freezeX && com.selection.w == 0 && gtk_widget_is_visible(rotation_dlg))
			com.selection.w = gfit->rx;
		if (!gui.freezeY && com.selection.h == 0 && gtk_widget_is_visible(rotation_dlg))
			com.selection.h = gfit->ry;

		if (gui.freezeX && gui.freezeY) { // Move selection
			com.selection.x = (data->zoomed.x - gui.start.x) + gui.origin.x;
			com.selection.y = (data->zoomed.y - gui.start.y) + gui.origin.y;
		}

		// Enforce a ratio and clamp selection to the image
		enforce_ratio_and_clamp();

		/* we have a new rectangular selection zone,
			* or an unselection (empty zone) */
		gui_function(new_selection_zone, NULL);

		// Terminate any specific selection modification mode
		gui.freezeX = gui.freezeY = FALSE;
	}
	return TRUE;
}

static gboolean select_preview1_release(mouse_data *data) {
	set_preview_area(0, data->zoomed.x, data->zoomed.y);
	*data->mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	// redraw to get the position of the new preview area
	redraw(REDRAW_OVERLAY);
	return TRUE;
}

static gboolean select_preview2_release(mouse_data *data) {
	set_preview_area(1, data->zoomed.x, data->zoomed.y);
	*data->mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	redraw(REDRAW_OVERLAY);
	return TRUE;
}

static gboolean get_comp_center_coordinate_release(mouse_data *data) {
	if (gui.comp_layer_centering) {
		gui.comp_layer_centering->center.x = data->zoomed.x;
		gui.comp_layer_centering->center.y = data->zoomed.y;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(gui.comp_layer_centering->centerbutton), FALSE);
	}
	*data->mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	redraw(REDRAW_OVERLAY);
	return TRUE;
}

static gboolean cut_select_release(mouse_data *data) {
	point tmp;
	tmp.x = data->zoomed.x;
	tmp.y = data->zoomed.y;
	if (*data->cutting == CUT_VERT_OR_HORIZ) {
		if (fabs(tmp.y - gui.cut.cut_start.y) > fabs(tmp.x - gui.cut.cut_start.x)) {
			tmp.x = gui.cut.cut_start.x;
		} else {
			tmp.y = gui.cut.cut_start.y;
		}
	}
	gui.cut.cut_end.x = tmp.x;
	gui.cut.cut_end.y = tmp.y;
	measure_line(gfit, gui.cut.cut_start, gui.cut.cut_end, gui.cut.pref_as);
	*data->cutting = CUT_NOT_CUTTING;
	redraw(REDRAW_OVERLAY);
	// Deselect the Cut button once the cut is made
	if (!gtk_widget_is_visible(cut_dialog))
		siril_open_dialog("cut_dialog");
	return TRUE;
}

static gboolean cut_wn1_release(mouse_data *data) {
	gui.cut.cut_wn1.x = data->zoomed.x;
	gui.cut.cut_wn1.y = data->zoomed.y;
	// Snap the selected pixel to the closest point on the line
	gui.cut.cut_wn1 = closest_point_on_line(gui.cut.cut_wn1, gui.cut.cut_start, gui.cut.cut_end);
	gchar* l1x = g_strdup_printf("%d", (int) gui.cut.cut_wn1.x);
	gchar* l1y = g_strdup_printf("%d", (int) gui.cut.cut_wn1.y);
	gtk_label_set_text(label_wn1_x, l1x);
	gtk_label_set_text(label_wn1_y, l1y);
	g_free(l1x);
	g_free(l1y);
	set_cursor("default");
	*data->mouse_status = MOUSE_ACTION_NONE;
	return TRUE;
}

static gboolean cut_wn2_release(mouse_data *data) {
	gui.cut.cut_wn2.x = data->zoomed.x;
	gui.cut.cut_wn2.y = data->zoomed.y;
	// Snap the selected pixel to the closest point on the line
	gui.cut.cut_wn2 = closest_point_on_line(gui.cut.cut_wn2, gui.cut.cut_start, gui.cut.cut_end);
	gchar* l2x = g_strdup_printf("%d", (int) gui.cut.cut_wn2.x);
	gchar* l2y = g_strdup_printf("%d", (int) gui.cut.cut_wn2.y);
	gtk_label_set_text(label_wn2_x, l2x);
	gtk_label_set_text(label_wn2_y, l2y);
	g_free(l2x);
	g_free(l2y);
	set_cursor("default");
	*data->mouse_status = MOUSE_ACTION_NONE;
	return TRUE;
}

static gboolean drag_release(mouse_data *data) {
	gui.translating = FALSE;
	return TRUE;
}

static gboolean measure_release (mouse_data *data) {
	if (gui.measure_start.x != -1.) {
		gui.measure_end.x = data->zoomed.x;
		gui.measure_end.y = data->zoomed.y;
		gboolean use_arcsec = (gfit->keywords.wcsdata.pltsolvd || gui.cut.pref_as);
		measure_line(gfit, gui.measure_start, gui.measure_end, use_arcsec);
		gui.measure_start.x = -1.;
		gui.measure_start.y = -1.;
		gui.measure_end.x = -1.;
		gui.measure_end.y = -1.;
	}
	return TRUE;
}

static gboolean draw_poly_release(mouse_data *data) {
	// Parse gui.drawing_polypoints into a Polygon and add it using add_user_polygon
	gui.drawing_polygon = FALSE;
	*data->mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	UserPolygon *poly = create_user_polygon_from_points(gui.drawing_polypoints);
	add_existing_polygon(poly, &gui.poly_ink, gui.poly_fill);
	queue_redraw(REDRAW_OVERLAY);
	// Free and NULL gui.drawing_polypoints
	g_slist_free_full(gui.drawing_polypoints, free);
	gui.drawing_polypoints = NULL;
	return TRUE;
}

static gboolean mask_add_poly_release(mouse_data *data) {
	// Parse gui.drawing_polypoints into a Polygon and add it using add_user_polygon
	gui.drawing_polygon = FALSE;
	*data->mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	UserPolygon *poly = create_user_polygon_from_points(gui.drawing_polypoints);
	// Free and NULL gui.drawing_polypoints
	g_slist_free_full(gui.drawing_polypoints, free);
	if (!gfit->mask) { // we need something to add the polygon to, so create a zeroes-like mask
		mask_create_zeroes_like(gfit, get_default_mask_bitpix());
	}
	set_poly_in_mask(poly, gfit, TRUE);
	free_user_polygon(poly);
	gui.drawing_polypoints = NULL;
	queue_redraw_mask();
	return TRUE;
}

static gboolean mask_clear_poly_release(mouse_data *data) {
	// Parse gui.drawing_polypoints into a Polygon and add it using add_user_polygon
	gui.drawing_polygon = FALSE;
	*data->mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	UserPolygon *poly = create_user_polygon_from_points(gui.drawing_polypoints);
	if (!gfit->mask) { // we need somthing to subtract the polygon from, so create a ones-like mask
		mask_create_ones_like(gfit, get_default_mask_bitpix());
	}
	set_poly_in_mask(poly, gfit, FALSE);
	// Free and NULL gui.drawing_polypoints
	g_slist_free_full(gui.drawing_polypoints, free);
	gui.drawing_polypoints = NULL;
	queue_redraw_mask();
	return TRUE;
}

static gboolean sample_mask_color_release(mouse_data *data) {
	mask_color_handle_image_click(data->zoomed.x, data->zoomed.y);
	return TRUE;
}

static gboolean show_popup_menu(mouse_data *data) {
	if (gui.cvport == MASK_VPORT) {
		do_popup_maskmenu(data->widget, NULL);
	} else if (*data->mouse_status != MOUSE_ACTION_DRAW_SAMPLES && *data->mouse_status != MOUSE_ACTION_PHOTOMETRY) {
		do_popup_graymenu(data->widget, NULL);
	}
	return TRUE;
}

mouse_function_metadata null_action = { mouse_nullfunction, NAME_NULL_ACTION, "", MOUSE_REF_NULL };

mouse_function_metadata open_if_unloaded_action = { open_if_nothing_loaded,
	NAME_OPEN_IF_UNLOADED_ACTION, N_("This action opens a new image if there is no image currently loaded"), FALSE, MOUSE_REF_OPEN_UNLOADED };

gboolean open_if_nothing_loaded(mouse_data *data) {
	/* with this action (if no images loaded)
	 * you can load an image This feature is in GIMP and I really
	 * love it: lazy world :).
	 */
	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		header_open_button_clicked();
		gui_function(launch_clipboard_survey, NULL);
	}
	return TRUE;
}

mouse_function_metadata measure_action = { measure_click, NAME_MEASURE_ACTION,
	N_("This action starts a measurement from one point in the image to another. Releasing the button ends the measurement."), TRUE, MOUSE_REF_MEASURE };

gboolean measure_click (mouse_data *data) {
	if (data->inside) {
		gui.measure_start.x = data->zoomed.x;
		gui.measure_start.y = data->zoomed.y;
		gui.measure_end.x = data->zoomed.x;
		gui.measure_end.y = data->zoomed.y;
	}
	register_release_callback(measure_release, data->event->button);
	return TRUE;
}

mouse_function_metadata zoom_action = { zoom_function, NAME_ZOOM_ACTION,
	N_("This action toggles the image zoom level as configured in Preferences / User Interface"), FALSE, MOUSE_REF_ZOOM };

gboolean zoom_function (mouse_data *data) {
	if (data->inside) {
		update_zoom_fit_button();
		switch (com.pref.gui.mmb_action) {
			case MMB_ZOOM_FIT:
				gui.zoom_value = ZOOM_FIT;
				reset_display_offset();
				break;
			case MMB_ZOOM_100:
				gui.zoom_value = ZOOM_NONE;
				double x = data->evpos.x;
				double y = data->evpos.y;
				cairo_matrix_transform_point(&gui.display_matrix, &data->evpos.x, &data->evpos.y);
				gui.display_offset.x = data->evpos.x - x;
				gui.display_offset.y = data->evpos.y - y;
				break;
			case MMB_ZOOM_TOGGLE:
				if (zoom_action_zooms_to_fit) {
					gui.zoom_value = ZOOM_FIT;
					reset_display_offset();
				} else {
					gui.zoom_value = ZOOM_NONE;
					double x = data->evpos.x;
					double y = data->evpos.y;
					cairo_matrix_transform_point(&gui.display_matrix, &data->evpos.x, &data->evpos.y);
					gui.display_offset.x = data->evpos.x - x;
					gui.display_offset.y = data->evpos.y - y;
				}
				break;
		}
		update_zoom_label();
		redraw(REDRAW_IMAGE);
		zoom_action_zooms_to_fit = !zoom_action_zooms_to_fit;
	}
	return TRUE;
}

mouse_function_metadata drag_action = { drag_click, NAME_DRAG_ACTION,
	N_("This action drags the image around the viewport"), TRUE, MOUSE_REF_DRAG };

gboolean drag_click(mouse_data *data) {
	if (data->inside) {

		// viewport translation
		gui.translating = TRUE;
		gui.start.x = (int)(data->event->x);
		gui.start.y = (int)(data->event->y);
	}
	register_release_callback(drag_release, data->event->button);
	return TRUE;
}

mouse_function_metadata main_action = { main_action_click, NAME_MAIN_ACTION,
	N_("This is the main Siril mouse action used for making selections, selecting background samples and previews, doing quick photometry etc."),
	TRUE, MOUSE_REF_MAIN_ACTION};

gboolean main_action_click(mouse_data *data) {
	if (data->inside) {
		point pt;
		int radius, s;
		gboolean right, left, bottom, top;
		rectangle area;
		struct phot_config *ps = NULL;
		switch (*data->mouse_status) {
			case MOUSE_ACTION_SELECT_REG_AREA: {
				if (gui.drawing) {
					gui.drawing = FALSE;
				} else {
					gui.drawing = TRUE;
					if (is_inside_of_sel(data->zoomed, data->zoom)) {
						// Move selection
						gui.freezeX = gui.freezeY = TRUE;
						gui.start = data->zoomed;
						gui.origin.x = com.selection.x;
						gui.origin.y = com.selection.y;
					} else {
						// Default values
						gui.freezeX = gui.freezeY = FALSE;
						// The order matters if the selection is so small that edge detection overlaps
						// and need to be the same as in the on_drawingarea_motion_notify_event()
						right = is_over_the_right_side_of_sel(data->zoomed, data->zoom);
						left = is_over_the_left_side_of_sel(data->zoomed, data->zoom);
						bottom = is_over_the_bottom_of_sel(data->zoomed, data->zoom);
						top = is_over_the_top_of_sel(data->zoomed, data->zoom);
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
							gui.start = data->zoomed;
							com.selection.h = 0;
							com.selection.w = 0;
						}
					}
				}
				redraw(REDRAW_OVERLAY);
				register_release_callback(select_reg_area_release, data->event->button);
				break;
			}
			case MOUSE_ACTION_DRAW_SAMPLES: {
				radius = get_background_sample_radius();

				pt.x = (gdouble) data->zoomed.x;
				pt.y = (gdouble) data->zoomed.y;

				if (pt.x + radius < gfit->rx && pt.y + radius < gfit->ry
						&& pt.x - radius > 0 && pt.y - radius > 0) {
					sample_mutex_lock();
					com.grad_samples = add_background_sample(com.grad_samples, gfit, pt);
					sample_mutex_unlock();

					redraw(REDRAW_OVERLAY);
					gui_function(redraw_previews, NULL);
				}
				break;
			}
			case MOUSE_ACTION_SELECT_PREVIEW1: {
				register_release_callback(select_preview1_release, data->event->button);
				break;
			}
			case MOUSE_ACTION_SELECT_PREVIEW2: {
				register_release_callback(select_preview2_release, data->event->button);
				break;
			}
			case MOUSE_ACTION_GET_COMP_CENTER_COORDINATE: {
				register_release_callback(get_comp_center_coordinate_release, data->event->button);
				break;
			}
			case MOUSE_ACTION_PHOTOMETRY: {
				s = com.pref.phot_set.outer * 1.2;
				area.x = data->zoomed.x - s;
				area.y = data->zoomed.y - s;
				area.w = s * 2;
				area.h = s * 2;
				if (data->zoomed.x - area.w > 0 && data->zoomed.x + area.w < gfit->rx
						&& data->zoomed.y - area.h > 0 && data->zoomed.y + area.h < gfit->ry) {
					ps = phot_set_adjusted_for_image(gfit);
					psf_error error = PSF_NO_ERR;
					gui.qphot = psf_get_minimisation(gfit, select_vport(gui.cvport), &area, TRUE, TRUE, ps, TRUE, com.pref.starfinder_conf.profile, &error);
					free(ps);
					if (gui.qphot) {
						if (!gui.qphot->phot_is_valid || error != PSF_NO_ERR) {
							free_psf(gui.qphot);
							gui.qphot = NULL;
							break;
						}
						gui.qphot->xpos = gui.qphot->x0 + area.x;
						if (gfit->top_down)
							gui.qphot->ypos = gui.qphot->y0 + area.y;
						else
							gui.qphot->ypos = area.y + area.h - gui.qphot->y0;
						redraw(REDRAW_OVERLAY);
						popup_psf_result(gui.qphot, &area, gfit);
					}
				}
				break;
			}
			case MOUSE_ACTION_CUT_SELECT: {
				// Reset the cut line before setting new coords in order to avoid
				// drawing artefacts
				gui.cut.cut_start.x = -1.;
				gui.cut.cut_start.y = -1.;
				gui.cut.cut_end.x = -1.;
				gui.cut.cut_end.y = -1.;
				if (data->event->state & GDK_SHIFT_MASK) {
					*data->cutting = CUT_VERT_OR_HORIZ;
				} else {
					*data->cutting = CUT_UNCONSTRAINED;
				}
				gui.cut.cut_start.x = data->zoomed.x;
				gui.cut.cut_start.y = data->zoomed.y;
				// This is a new cut line so reset any spectroscopic wavenumber points
				gui.cut.cut_wn1.x = -1.;
				gui.cut.cut_wn1.y = -1.;
				gui.cut.cut_wn2.x = -1.;
				gui.cut.cut_wn2.y = -1.;
				gtk_label_set_text(label_wn1_x, "");
				gtk_label_set_text(label_wn1_y, "");
				gtk_label_set_text(label_wn2_x, "");
				gtk_label_set_text(label_wn2_y, "");
				update_spectro_labels();
				register_release_callback(cut_select_release, data->event->button);
				break;
			}
			case MOUSE_ACTION_CUT_WN1: {
				register_release_callback(cut_wn1_release, data->event->button);
				break;
			}
			case MOUSE_ACTION_CUT_WN2: {
				register_release_callback(cut_wn2_release, data->event->button);
				break;
			}
			case MOUSE_ACTION_DRAW_POLY: {
				g_assert(gui.drawing_polypoints == NULL);
				point *ev = malloc(sizeof(point));
				ev->x = data->zoomed.x;
				ev->y = data->zoomed.y;
				gui.drawing_polypoints = g_slist_prepend(gui.drawing_polypoints, ev);
				register_release_callback(draw_poly_release, data->event->button);
				gui.drawing_polygon = TRUE;
				break;
			}
			case MOUSE_ACTION_ADD_POLY_TO_MASK: {
				g_assert(gui.drawing_polypoints == NULL);
				point *ev = malloc(sizeof(point));
				ev->x = data->zoomed.x;
				ev->y = data->zoomed.y;
				gui.drawing_polypoints = g_slist_prepend(gui.drawing_polypoints, ev);
				register_release_callback(mask_add_poly_release, data->event->button);
				gui.drawing_polygon = TRUE;
				break;
			}
			case MOUSE_ACTION_CLEAR_POLY_FROM_MASK: {
				g_assert(gui.drawing_polypoints == NULL);
				point *ev = malloc(sizeof(point));
				ev->x = data->zoomed.x;
				ev->y = data->zoomed.y;
				gui.drawing_polypoints = g_slist_prepend(gui.drawing_polypoints, ev);
				register_release_callback(mask_clear_poly_release, data->event->button);
				gui.drawing_polygon = TRUE;
				break;
			}
			case MOUSE_ACTION_SAMPLE_MASK_COLOR: {
				register_release_callback(sample_mask_color_release, data->event->button);
				break;
			}
			default: {
				break;
			}
		}
	}
	return TRUE;
}

mouse_function_metadata second_action = { second_action_click, NAME_SECOND_ACTION,
	N_("Remove background samples or photometry boxes, or show the mouse context menu"), FALSE, MOUSE_REF_SECOND_ACTION };

gboolean second_action_click(mouse_data *data) {
	if (data->inside) {
		// Reset the cut line if one has been drawn
		gui.cut.cut_start.x = -1;
		gui.cut.cut_start.y = -1;
		gui.cut.cut_end.x = -1;
		gui.cut.cut_end.y = -1;
		if (*data->mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
			point pt;
			int radius = (int) (25 / 2);

			pt.x = (gdouble) data->zoomed.x;
			pt.y = (gdouble) data->zoomed.y;

			if (pt.x + radius <= gfit->rx && pt.y + radius <= gfit->ry
					&& pt.x - radius >= 0 && pt.y - radius >= 0) {
				sample_mutex_lock();
				com.grad_samples = remove_background_sample(com.grad_samples, gfit, pt);
				sample_mutex_unlock();

				redraw(REDRAW_OVERLAY);
				gui_function(redraw_previews, NULL);
			}
		} else if (*data->mouse_status == MOUSE_ACTION_PHOTOMETRY) {
			if (sequence_is_loaded()) {
				int s = com.pref.phot_set.outer * 1.2;
				rectangle area = { data->zoomed.x - s, data->zoomed.y - s, s * 2, s * 2 };
				if (data->zoomed.x - area.w > 0 && data->zoomed.x + area.w < gfit->rx
						&& data->zoomed.y - area.h > 0 && data->zoomed.y + area.h < gfit->ry) {
					memcpy(&com.selection, &area, sizeof(rectangle));
					seq_qphot(&com.seq, select_vport(gui.cvport));
					delete_selected_area();
				}
			}
		} else {
			// If we aren't in a special action, releasing the mouse button shows the popup menu
			register_release_callback(show_popup_menu, data->event->button);
		}
	}
	return TRUE;
}

mouse_function_metadata photometry_box_action = { photometry_box_set,
	NAME_PHOTOMETRY_BOX_ACTION, N_("This action sets the photometry box"), FALSE, MOUSE_REF_PHOTOMETRY_BOX };

gboolean photometry_box_set(mouse_data *data) {
	/* Ctrl middle-click to set the photometry box */
	if (data->inside && data->event->state & get_primary()) {
		double dX = 1.5 * com.pref.phot_set.outer;
		double dY = dX;
		double w = 3 * com.pref.phot_set.outer;
		double h = w;

		if ((dX <= data->zoomed.x) && (dY <= data->zoomed.y)
				&& (data->zoomed.x - dX + w < gfit->rx)
				&& (data->zoomed.y - dY + h < gfit->ry)) {

			com.selection.x = data->zoomed.x - dX;
			com.selection.y = data->zoomed.y - dY;
			com.selection.w = w;
			com.selection.h = h;

			gui_function(new_selection_zone, NULL);
		}
	}
	return TRUE;
}

static gpointer mouse_save_as(gpointer user_data) {
	on_header_save_as_button_clicked();
	return GINT_TO_POINTER(0);
}

mouse_function_metadata save_on_click_action = { save_on_click,
	NAME_SAVE_ON_CLICK_ACTION, N_("This action saves the current image"), FALSE, MOUSE_REF_SAVE };

	gboolean save_on_click(mouse_data *data) {
	if (single_image_is_loaded() && com.uniq->fileexist) {
		savefits(com.uniq->filename, gfit);
	} else {
		start_in_new_thread(mouse_save_as, NULL);
	}
	return TRUE;
}

mouse_function_metadata save_as_on_click_action = { save_as_on_click,
	NAME_SAVE_AS_ON_CLICK_ACTION, N_("This action saves the current image with a new name"), FALSE, MOUSE_REF_SAVE_AS };

gboolean save_as_on_click(mouse_data *data) {
	start_in_new_thread(mouse_save_as, NULL);
	return TRUE;
}

mouse_function_metadata snapshot_on_click_action = { snapshot_on_click,
	NAME_SNAPSHOT_ON_CLICK_ACTION, N_("This action saves a snapshot of the current image"), FALSE, MOUSE_REF_SNAPSHOT };

gboolean snapshot_on_click(mouse_data *data) {
	on_header_snapshot_button_clicked(FALSE);
	return TRUE;
}

mouse_function_metadata undo_on_click_action = { undo_on_click,
	NAME_UNDO_ON_CLICK_ACTION, N_("This action undoes the last operation"), FALSE, MOUSE_REF_UNDO };

gboolean undo_on_click(mouse_data *data) {
	set_cursor_waiting(TRUE);
	undo_display_data(UNDO);
	set_cursor_waiting(FALSE);
	return TRUE;
}

mouse_function_metadata redo_on_click_action = { redo_on_click,
	NAME_REDO_ON_CLICK_ACTION, N_("This action redoes the last operation"), FALSE, MOUSE_REF_REDO };

gboolean redo_on_click(mouse_data *data) {
	set_cursor_waiting(TRUE);
	undo_display_data(REDO);
	set_cursor_waiting(FALSE);
	return TRUE;
}

mouse_function_metadata findstar_on_click_action = { findstar_on_click,
	NAME_FINDSTAR_ON_CLICK_ACTION, N_("This action finds stars in the current image"), FALSE, MOUSE_REF_FINDSTAR };

extern void on_process_starfinder_button_clicked(GtkButton *button, gpointer user_data);

gboolean findstar_on_click(mouse_data *data) {
	set_cursor_waiting(TRUE);
	on_process_starfinder_button_clicked(NULL, NULL);
	set_cursor_waiting(FALSE);
	return TRUE;
}

mouse_function_metadata platesolve_on_click_action = { platesolve_on_click,
	NAME_PLATESOLVE_ON_CLICK_ACTION, N_("This action plate solves the current image"), FALSE, MOUSE_REF_PLATESOLVE };

gboolean platesolve_on_click(mouse_data *data) {
	siril_open_dialog("astrometry_dialog");
	return TRUE;
}

mouse_function_metadata spcc_on_click_action = { spcc_on_click,
	NAME_SPCC_ON_CLICK_ACTION, N_("This action opens the spectrophotometric color calibration dialog"), FALSE, MOUSE_REF_SPCC };

gboolean spcc_on_click(mouse_data *data) {
	siril_open_dialog("s_pcc_dialog");
	return TRUE;
}

/*
 * Mouse Scroll functions
 */

scroll_function_metadata scroll_null_action = { scroll_nullfunction, NAME_NULL_ACTION, "", SCROLL_REF_NULL };

gboolean scroll_nullfunction(scroll_data *data) {
	return FALSE;
}

scroll_function_metadata scroll_zooms_action = { scroll_zooms, NAME_SCROLL_ZOOM_ACTION,
	N_("Scroll wheel zooms in / out") , SCROLL_REF_ZOOM};

gboolean scroll_zooms(scroll_data *data) {
	// We stored the GdkEventScroll* in a different pointer type in data
	// to save having a completely different struct for scroll events;
	// cast it back now
	GdkEventScroll *event = data->event;
	if (!single_image_is_loaded() && !sequence_is_loaded())
		return FALSE;

	if (gui.icc.iso12646)
		disable_iso12646_conditions(FALSE, TRUE, TRUE);

	point delta;

	// The handler is written to act on either horizontal or vertical scroll events
	// However the handler will only ever be called for one or the other depending
	// on the configuration.
	double speed_limit = com.pref.gui.mouse_speed_limit;
	switch (event->direction) {
		case GDK_SCROLL_SMOOTH:
			gdk_event_get_scroll_deltas((GdkEvent*) event, &delta.x, &delta.y);
			if (speed_limit) {
				if (delta.x > speed_limit)
					delta.x = speed_limit;
				else if (delta.x < -speed_limit)
					delta.x = -speed_limit;
				if (delta.y > speed_limit)
					delta.y = speed_limit;
				else if (delta.y < -speed_limit)
					delta.y = -speed_limit;
			}
			if (data->direction == MOUSE_VERTICAL_SCROLL) {
				if (delta.y < 0) {
					return update_zoom(event->x, event->y, ((ZOOM_IN - 1.0) * fabs(delta.y)) + 1.0);
				} else if (delta.y > 0) {
					return update_zoom(event->x, event->y, 1.0 / (((ZOOM_IN - 1.0) * delta.y) + 1.0));
				}
			} else if (data->direction == MOUSE_HORIZ_SCROLL) {
				if (delta.x > 0) {
					return update_zoom(event->x, event->y, ((ZOOM_IN - 1.0) * delta.x) + 1.0);
				} else if (delta.x < 0) {
					return update_zoom(event->x, event->y, 1.0 / (((ZOOM_IN - 1.0) * fabs(delta.x)) + 1.0));
				}
			}
			break;
			// Fallthrough intentional
		case GDK_SCROLL_DOWN:
		case GDK_SCROLL_LEFT:
			return update_zoom(event->x, event->y, ZOOM_IN);
			// Fallthrough intentional
		case GDK_SCROLL_UP:
		case GDK_SCROLL_RIGHT:
			return update_zoom(event->x, event->y, ZOOM_OUT);
		default:
			break;
	}
	return FALSE;
}

scroll_function_metadata scroll_traverses_seq_action = { scroll_traverses_sequence,
	NAME_SCROLL_SEQ_ACTION, N_("Scroll wheel scrolls through the frames in the current sequence"), SCROLL_REF_SEQ };

static double seq_accumulator = 0.0;

gboolean scroll_traverses_sequence(scroll_data *data) {
	GdkEventScroll *event = data->event;
	if (!sequence_is_loaded())
		return FALSE;
	gint idx = -1;
	point delta;
	switch (event->direction) {
		case GDK_SCROLL_SMOOTH:
			gdk_event_get_scroll_deltas((GdkEvent*) data->event, &delta.x, &delta.y);
			if (data->direction == MOUSE_VERTICAL_SCROLL) {
				seq_accumulator += delta.y;
			} else if (data->direction == MOUSE_HORIZ_SCROLL) {
				seq_accumulator += delta.x;
			}
			break;
		// Fallthrough intentional
		case GDK_SCROLL_DOWN:
		case GDK_SCROLL_LEFT:
			seq_accumulator = 0.0;
			idx = com.seq.current - 1;
			break;
		// Fallthrough intentional
		case GDK_SCROLL_UP:
		case GDK_SCROLL_RIGHT:
			seq_accumulator = 0.0;
			idx = com.seq.current + 1;
			break;
		default:
			break;
	}
	if (seq_accumulator >= 1.0) {
		idx = com.seq.current + (com.seq.current < com.seq.number ? 1 : 0);
		seq_accumulator = 0.0;
	} else if (seq_accumulator <= -1.0) {
		idx = com.seq.current - (com.seq.current> 0 ? 1 : 0);
		seq_accumulator = 0.0;
	}
	if (idx >= 0 && idx < com.seq.number) {
		if (seq_load_image(&com.seq, idx, TRUE)) { // if loading fails, we fall back reloading the reference image
			seq_load_image(&com.seq, com.seq.reference_image, TRUE);
			siril_log_color_message(_("Error: image load failed, reloading reference image!\n"), "red");
		}
	}
	return FALSE;
}

scroll_function_metadata scroll_changes_tab_action = { scroll_changes_tab,
	NAME_SCROLL_TAB_ACTION, N_("Scroll wheel scrolls through the Siril GUI tabs"), SCROLL_REF_TABS };

static double tab_accumulator = 0.0;

gboolean scroll_changes_tab(scroll_data *data) {
	GdkEventScroll *event = data->event;
	GtkNotebook* notebook = GTK_NOTEBOOK(lookup_widget("notebook_center_box"));
	gint tab = gtk_notebook_get_current_page(notebook);
	point delta;
	gdk_event_get_scroll_deltas((GdkEvent*) data->event, &delta.x, &delta.y);
	switch (event->direction) {
		case GDK_SCROLL_SMOOTH:
			if (data->direction == MOUSE_VERTICAL_SCROLL) {
				tab_accumulator += delta.y;
			} else if (data->direction == MOUSE_HORIZ_SCROLL) {
				tab_accumulator += delta.x;
			}
			break;
		// Fallthrough intentional
		case GDK_SCROLL_DOWN:
		case GDK_SCROLL_LEFT:
			tab = tab > 0 ? tab - 1 : tab;
			tab_accumulator = 0.0;
			break;
		// Fallthrough intentional
		case GDK_SCROLL_UP:
		case GDK_SCROLL_RIGHT:
			tab = tab < OUTPUT_LOGS ? tab + 1 : tab;
			tab_accumulator = 0.0;
			break;
		default:
			break;
	}
	if (tab_accumulator >= 1.0) {
		tab = tab < OUTPUT_LOGS ? tab + 1 : tab;
		tab_accumulator = 0.0;
	} else if (tab_accumulator <= -1.0) {
		tab = tab > 0 ? tab - 1 : tab;
		tab_accumulator = 0.0;
	}
	gtk_notebook_set_current_page(notebook, tab);
	return FALSE;
}
