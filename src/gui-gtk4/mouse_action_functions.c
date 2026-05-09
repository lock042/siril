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
#include "gui-gtk4/cut.h"
#include "gui-gtk4/icc_profile.h"
#include "core/processing.h"
#include "core/undo.h"
#include "algos/background_extraction.h"
#include "algos/photometry.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "gui-gtk4/open_dialog.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/PSF_list.h"
#include "image_interactions.h"
#include "gui-gtk4/mouse_action_functions.h"
#include "gui-gtk4/masks_gui.h"
#include "image_display.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/save_dialog.h"
#include "gui-gtk4/siril_actions.h"
#include "gui-gtk4/utils.h"
#include "progress_and_log.h"
#include "registration_preview.h"
#include "gui-gtk4/user_polygons.h"

// Tracks whether double middle click will zoom to fit or zoom to 1:1, for the
// toggle
static gboolean zoom_action_zooms_to_fit = FALSE;

// caching widgets
static GtkWidget *rotation_dlg = NULL;
static GtkWidget *cut_dialog = NULL, *dynpsf_dlg = NULL;
static GtkLabel *label_wn1_x = NULL, *label_wn1_y = NULL, *label_wn2_x = NULL, *label_wn2_y = NULL;
static GtkCheckButton *maf_mi_free = NULL, *maf_mi_preserve = NULL;
static GtkCheckButton *maf_mi_16_9 = NULL, *maf_mi_4_3 = NULL, *maf_mi_3_2 = NULL;
static GtkCheckButton *maf_mi_1_1 = NULL, *maf_mi_3_4 = NULL, *maf_mi_2_3 = NULL, *maf_mi_9_16 = NULL;
static GtkWidget *maf_mi_all = NULL;
static GtkCheckButton *maf_mi_guides_0 = NULL, *maf_mi_guides_2 = NULL;
static GtkCheckButton *maf_mi_guides_3 = NULL, *maf_mi_guides_5 = NULL;
static GtkCheckButton *maf_bkg_grad_descent = NULL;
static GtkNotebook *maf_notebook_center_box = NULL;

static void mouse_action_functions_init_statics(void) {
	if (rotation_dlg) return;
	rotation_dlg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "rotation_dialog"));
	label_wn1_x = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_wn1_x"));
	label_wn1_y = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_wn1_y"));
	label_wn2_x = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_wn2_x"));
	label_wn2_y = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_wn2_y"));
	cut_dialog = GTK_WIDGET(gtk_builder_get_object(gui.builder, "cut_dialog"));
	dynpsf_dlg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stars_list_window"));
	maf_mi_free = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_free"));
	maf_mi_preserve = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_preserve"));
	maf_mi_16_9 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_16_9"));
	maf_mi_4_3 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_4_3"));
	maf_mi_3_2 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_3_2"));
	maf_mi_1_1 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_1_1"));
	maf_mi_3_4 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_3_4"));
	maf_mi_2_3 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_2_3"));
	maf_mi_9_16 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_9_16"));
	maf_mi_all = GTK_WIDGET(gtk_builder_get_object(gui.builder, "menuitem_selection_all"));
	maf_mi_guides_0 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_guides_0"));
	maf_mi_guides_2 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_guides_2"));
	maf_mi_guides_3 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_guides_3"));
	maf_mi_guides_5 = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "menuitem_selection_guides_5"));
	maf_bkg_grad_descent = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "bkg_grad_descent_button"));
	maf_notebook_center_box = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook_center_box"));
}

/* Mouse release functions should be coded as static gboolean function(mouse_data *)
 * and only be referenced from the button press functions further down this file: if
 * a button press function needs to set a mouse release action to be handled by the
 * button release event handler, it must do so using register_release_callback()
 * Mouse release functions can call other functions, some of which are declared as
 * static functions at the top of the file.
 */

/* GTK4 right-click context menus: re-parent the existing GtkPopover
 * (menugray / menumask_rmb) to the drawingarea where the click landed,
 * point at the click position, and pop up.  We reparent on every show
 * because the user can right-click on any of the five viewport
 * drawingareas (R/G/B/RGB/Mask) and a popover can only have one parent
 * at a time. */
/* Hardcoded GMenu structure for the drawingarea right-click context
 * menus.  We tried walking the builder-loaded GtkPopover's widget tree
 * to discover sections and submenus, but the legacy GTK3
 * <child type="submenu"> markup parses inconsistently in GTK4 and the
 * walk produced misnamed/duplicated submenus.  Hardcoding is the most
 * reliable path; it also lets us drop submenu items that would need
 * stateful GAction wiring (the Selection ratio / guide toggles) which
 * the existing .ui implements via signal handlers.
 *
 * Each entry is { label, action } where:
 *   - label = NULL                   → ends the current section (separator)
 *   - label set, action = NULL       → starts a submenu (next entries until matching end)
 *   - label = "<<end>>"              → ends the current submenu
 *   - label set, action set          → regular action item
 */
typedef struct { const char *label; const char *action; } MenuEntry;
#define MENU_SEP        { NULL,                              NULL                 }
#define MENU_SUBMENU(L) { L,                                 NULL                 }
#define MENU_END        { "<<end>>",                         NULL                 }

static const MenuEntry menugray_items[] = {
	{ N_("Undo"),                                    "win.undo" },
	{ N_("Redo"),                                    "win.redo" },
	MENU_SEP,
	{ N_("PSF"),                                     "win.psf" },
	{ N_("PSF for the Sequence"),                    "win.seq-psf" },
	{ N_("Pick a Star"),                             "win.pickstar" },
	MENU_SEP,
	MENU_SUBMENU(N_("Selection")),
		/* Aspect-ratio radio group.  GMenu renders entries that target the
		 * same stateful action as a radio set; the action's current state
		 * decides which entry is checked. */
		{ N_("Free ratio"),                          "win.aspect-ratio::free" },
		{ N_("Preserve ratio"),                      "win.aspect-ratio::preserve" },
		MENU_SEP,
		{ "16:9",                                    "win.aspect-ratio::16:9" },
		{ "3:2",                                     "win.aspect-ratio::3:2" },
		{ "4:3",                                     "win.aspect-ratio::4:3" },
		{ "1:1",                                     "win.aspect-ratio::1:1" },
		{ "3:4",                                     "win.aspect-ratio::3:4" },
		{ "2:3",                                     "win.aspect-ratio::2:3" },
		{ "9:16",                                    "win.aspect-ratio::9:16" },
		MENU_SEP,
		{ N_("Select All"),                          "win.select-all" },
		MENU_SEP,
		/* Selection guides — second radio group on a separate stateful
		 * action so the two groups don't share state. */
		{ N_("No guides"),                           "win.selection-guides::0" },
		{ N_("Center lines"),                        "win.selection-guides::2" },
		{ N_("Rule of thirds"),                      "win.selection-guides::3" },
		{ N_("Rule of fifths"),                      "win.selection-guides::5" },
	MENU_END,
	MENU_SUBMENU(N_("Crop")),
		{ N_("Crop"),                                "win.crop" },
		{ N_("Rotate & Crop"),                       "win.rotation-processing" },
		{ N_("Crop Sequence..."),                    "win.seq-crop" },
	MENU_END,
	MENU_SUBMENU(N_("Mask")),
		{ N_("Mask from current image"),             "win.mask_from_image" },
		{ N_("Mask from color"),                     "win.mask_from_color" },
		{ N_("Mask from stars"),                     "win.mask_from_stars" },
		{ N_("Mask from image file"),                "win.mask_from_file" },
		MENU_SEP,
		{ N_("Add freehand to mask"),                "win.mask_add_from_poly" },
		{ N_("Subtract freehand from mask"),         "win.mask_clear_from_poly" },
		MENU_SEP,
		{ N_("Clear mask"),                          "win.clear_mask" },
	MENU_END,
	MENU_SUBMENU(N_("ROI")),
		{ N_("Set ROI to selection"),                "win.set_roi" },
		{ N_("Clear ROI"),                           "win.clear_roi" },
	MENU_END,
	MENU_SUBMENU(N_("RGB align")),
		{ N_("Global alignment"),                    "win.align-global" },
		{ N_("KOMBAT alignment (planetary/deep-sky)"), "win.align-kombat" },
		{ N_("One star registration (deep-sky)"),    "win.align-psf" },
		{ N_("Image pattern alignment (planetary/deep-sky)"), "win.align-dft" },
	MENU_END,
};

static const MenuEntry menumask_items[] = {
	MENU_SUBMENU(N_("Stretch")),
		{ N_("Autostretch Mask"),                    "win.autostretch_mask" },
		{ N_("Apply Thresholds to Mask"),            "win.threshold_mask" },
	MENU_END,
	MENU_SUBMENU(N_("Modify")),
		{ N_("Blur Mask"),                           "win.blur_mask" },
		{ N_("Feather Mask"),                        "win.feather_mask" },
		{ N_("Invert Mask"),                         "win.invert_mask" },
		{ N_("Multiply mask"),                       "win.scale_mask" },
		{ N_("Gradient of mask"),                    "win.mask_from_gradient" },
	MENU_END,
	MENU_SEP,
	{ N_("Clear Mask"),                              "win.clear_mask" },
};

/* Recursive builder.  *idx is the cursor through items[]; advances past
 * the closing MENU_END marker on return.  Sections are flushed into the
 * `current` GMenu when a separator or submenu boundary is hit. */
static void build_menu_recursive(GMenu *current, const MenuEntry *items,
                                 size_t total, size_t *idx) {
	GMenu *section = g_menu_new();
	gboolean section_has_items = FALSE;
	while (*idx < total) {
		const MenuEntry *e = &items[*idx];
		if (e->label && g_strcmp0(e->label, "<<end>>") == 0) {
			(*idx)++;
			break;
		}
		if (!e->label) {
			/* Separator: flush current section. */
			if (section_has_items) {
				g_menu_append_section(current, NULL, G_MENU_MODEL(section));
				g_object_unref(section);
				section = g_menu_new();
				section_has_items = FALSE;
			}
			(*idx)++;
			continue;
		}
		if (!e->action) {
			/* Submenu opener.  Recurse into a child GMenu. */
			GMenu *sub = g_menu_new();
			(*idx)++;
			build_menu_recursive(sub, items, total, idx);
			g_menu_append_submenu(section, _(e->label), G_MENU_MODEL(sub));
			g_object_unref(sub);
			section_has_items = TRUE;
			continue;
		}
		g_menu_append(section, _(e->label), e->action);
		section_has_items = TRUE;
		(*idx)++;
	}
	if (section_has_items)
		g_menu_append_section(current, NULL, G_MENU_MODEL(section));
	g_object_unref(section);
}

static GMenu *build_menu_model(const MenuEntry *items, size_t n) {
	GMenu *root = g_menu_new();
	size_t idx = 0;
	build_menu_recursive(root, items, n, &idx);
	return root;
}

static GMenu *build_menu_model_from_popover(const char *popover_id) {
	if (g_strcmp0(popover_id, "menugray") == 0)
		return build_menu_model(menugray_items, G_N_ELEMENTS(menugray_items));
	if (g_strcmp0(popover_id, "menumask_rmb") == 0)
		return build_menu_model(menumask_items, G_N_ELEMENTS(menumask_items));
	return NULL;
}

/* Replace the builder-loaded GtkPopover with a runtime-built
 * GtkPopoverMenu — the latter is GTK4's well-tested context-menu path
 * (it's what GtkMenuButton.set_menu_model creates internally) and
 * doesn't have the surface-clipping problems we saw when popping a
 * builder GtkPopover from inside a deeply-nested widget hierarchy. */
static void show_drawingarea_popover(const char *popover_id,
                                     GtkWidget *anchor, double x, double y) {
	if (!anchor || !GTK_IS_WIDGET(anchor)) return;

	/* Cache the GtkPopoverMenu per popover_id.  The model is built once;
	 * re-parenting moves it between drawingareas as needed. */
	static GHashTable *cache = NULL;
	if (!cache)
		cache = g_hash_table_new(g_str_hash, g_str_equal);

	GtkPopover *popover = g_hash_table_lookup(cache, popover_id);
	if (!popover) {
		GMenu *model = build_menu_model_from_popover(popover_id);
		if (!model) return;
		popover = GTK_POPOVER(gtk_popover_menu_new_from_model(G_MENU_MODEL(model)));
		g_object_unref(model);
		gtk_popover_set_has_arrow(popover, FALSE);
		gtk_popover_set_autohide(popover, TRUE);
		/* GtkPopoverMenu's internal stack uses GTK_SCROLLED_WINDOW_NEVER for
		 * its hadjustment but the natural height clamps to a small default,
		 * forcing scrollbars on short menus.  Bump the min size on a
		 * per-menu basis so each fits without scrolling. */
		int min_w = 220, min_h = 240;
		if (g_strcmp0(popover_id, "menugray") == 0) {
			min_h = 360;   /* Selection submenu has the deepest content */
		}
		gtk_widget_set_size_request(GTK_WIDGET(popover), min_w, min_h);
		/* Keep a strong ref via the cache so the popover survives
		 * unparent/reparent cycles. */
		g_object_ref_sink(popover);
		g_hash_table_insert(cache, (gpointer) popover_id, popover);
	}

	gtk_popover_popdown(popover);
	GtkWidget *cur_parent = gtk_widget_get_parent(GTK_WIDGET(popover));
	if (cur_parent != anchor) {
		if (cur_parent)
			gtk_widget_unparent(GTK_WIDGET(popover));
		gtk_widget_set_parent(GTK_WIDGET(popover), anchor);
	}

	/* Mirror the GTK3 do_popup_graymenu_unused() behaviour: refresh the
	 * radio actions for menugray's Selection submenu from the current
	 * gui.ratio / com.pref.gui.selection_guides values just before showing
	 * the menu, so the radio dot lands on the right entry even when the
	 * state was changed outside the menu. */
	if (g_strcmp0(popover_id, "menugray") == 0) {
		GtkRoot *root = gtk_widget_get_root(anchor);
		if (root && GTK_IS_APPLICATION_WINDOW(root))
			sync_selection_action_state(root);
	}

	GdkRectangle rect = { (int) x, (int) y, 1, 1 };
	gtk_popover_set_pointing_to(popover, &rect);
	gtk_popover_popup(popover);
}

static void do_popup_graymenu(GtkWidget *my_widget, double x, double y) {
	show_drawingarea_popover("menugray", my_widget, x, y);
}
static void do_popup_maskmenu(GtkWidget *my_widget, double x, double y) {
	show_drawingarea_popover("menumask_rmb", my_widget, x, y);
}

void cache_widgets() {
	mouse_action_functions_init_statics();
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
		siril_toggle_set_active(GTK_WIDGET(gui.comp_layer_centering->centerbutton), FALSE);
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
	double ink[4] = { gui.poly_ink.red, gui.poly_ink.green, gui.poly_ink.blue, gui.poly_ink.alpha };
	add_existing_polygon(poly, ink, gui.poly_fill);
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
	undo_save_state(gfit, _("Add polygon to image mask"));
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
	undo_save_state(gfit, _("Subtract polygon from image mask"));
	if (!gfit->mask) { // we need somthing to subtract the polygon from, so create a ones-like mask
		mask_create_ones_like(gfit, get_default_mask_bitpix());
	}
	set_poly_in_mask(poly, gfit, FALSE);
	free_user_polygon(poly);
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
		do_popup_maskmenu(data->widget, data->x, data->y);
	} else if (*data->mouse_status != MOUSE_ACTION_DRAW_SAMPLES && *data->mouse_status != MOUSE_ACTION_PHOTOMETRY) {
		do_popup_graymenu(data->widget, data->x, data->y);
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
	register_release_callback(measure_release, data->button);
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
		gui.start.x = (int)(data->x);
		gui.start.y = (int)(data->y);
	}
	register_release_callback(drag_release, data->button);
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
				register_release_callback(select_reg_area_release, data->button);
				break;
			}
			case MOUSE_ACTION_DRAW_SAMPLES: {
				radius = get_background_sample_radius();

				pt.x = (gdouble) data->zoomed.x;
				pt.y = (gdouble) data->zoomed.y;

				if (pt.x + radius < gfit->rx && pt.y + radius < gfit->ry
						&& pt.x - radius > 0 && pt.y - radius > 0) {
					gboolean gd = siril_toggle_get_active(GTK_WIDGET(maf_bkg_grad_descent));
					sample_mutex_lock();
					com.grad_samples = add_background_sample(com.grad_samples, gfit, pt, gd);
					sample_mutex_unlock();

					redraw(REDRAW_OVERLAY);
					gui_function(redraw_previews, NULL);
				}
				break;
			}
			case MOUSE_ACTION_SELECT_PREVIEW1: {
				register_release_callback(select_preview1_release, data->button);
				break;
			}
			case MOUSE_ACTION_SELECT_PREVIEW2: {
				register_release_callback(select_preview2_release, data->button);
				break;
			}
			case MOUSE_ACTION_GET_COMP_CENTER_COORDINATE: {
				register_release_callback(get_comp_center_coordinate_release, data->button);
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
				if (data->state & GDK_SHIFT_MASK) {
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
				register_release_callback(cut_select_release, data->button);
				break;
			}
			case MOUSE_ACTION_CUT_WN1: {
				register_release_callback(cut_wn1_release, data->button);
				break;
			}
			case MOUSE_ACTION_CUT_WN2: {
				register_release_callback(cut_wn2_release, data->button);
				break;
			}
			case MOUSE_ACTION_DRAW_POLY: {
				g_assert(gui.drawing_polypoints == NULL);
				point *ev = malloc(sizeof(point));
				ev->x = data->zoomed.x;
				ev->y = data->zoomed.y;
				gui.drawing_polypoints = g_slist_prepend(gui.drawing_polypoints, ev);
				register_release_callback(draw_poly_release, data->button);
				gui.drawing_polygon = TRUE;
				break;
			}
			case MOUSE_ACTION_ADD_POLY_TO_MASK: {
				g_assert(gui.drawing_polypoints == NULL);
				point *ev = malloc(sizeof(point));
				ev->x = data->zoomed.x;
				ev->y = data->zoomed.y;
				gui.drawing_polypoints = g_slist_prepend(gui.drawing_polypoints, ev);
				register_release_callback(mask_add_poly_release, data->button);
				gui.drawing_polygon = TRUE;
				break;
			}
			case MOUSE_ACTION_CLEAR_POLY_FROM_MASK: {
				g_assert(gui.drawing_polypoints == NULL);
				point *ev = malloc(sizeof(point));
				ev->x = data->zoomed.x;
				ev->y = data->zoomed.y;
				gui.drawing_polypoints = g_slist_prepend(gui.drawing_polypoints, ev);
				register_release_callback(mask_clear_poly_release, data->button);
				gui.drawing_polygon = TRUE;
				break;
			}
			case MOUSE_ACTION_SAMPLE_MASK_COLOR: {
				register_release_callback(sample_mask_color_release, data->button);
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
			register_release_callback(show_popup_menu, data->button);
		}
	}
	return TRUE;
}

mouse_function_metadata photometry_box_action = { photometry_box_set,
	NAME_PHOTOMETRY_BOX_ACTION, N_("This action sets the photometry box"), FALSE, MOUSE_REF_PHOTOMETRY_BOX };

gboolean photometry_box_set(mouse_data *data) {
	/* Ctrl middle-click to set the photometry box */
	if (data->inside && data->state & get_primary()) {
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

	if (com.gui_icc.iso12646)
		disable_iso12646_conditions(FALSE, TRUE, TRUE);

	point delta;

	// The handler is written to act on either horizontal or vertical scroll events
	// However the handler will only ever be called for one or the other depending
	// on the configuration.
	double speed_limit = com.pref.gui.mouse_speed_limit;
	switch (data->direction_raw) {
		case GDK_SCROLL_SMOOTH:
			gdk_scroll_event_get_deltas(event, &delta.x, &delta.y);
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
					return update_zoom(data->x, data->y, ((ZOOM_IN - 1.0) * fabs(delta.y)) + 1.0);
				} else if (delta.y > 0) {
					return update_zoom(data->x, data->y, 1.0 / (((ZOOM_IN - 1.0) * delta.y) + 1.0));
				}
			} else if (data->direction == MOUSE_HORIZ_SCROLL) {
				if (delta.x > 0) {
					return update_zoom(data->x, data->y, ((ZOOM_IN - 1.0) * delta.x) + 1.0);
				} else if (delta.x < 0) {
					return update_zoom(data->x, data->y, 1.0 / (((ZOOM_IN - 1.0) * fabs(delta.x)) + 1.0));
				}
			}
			break;
			// Fallthrough intentional
		case GDK_SCROLL_DOWN:
		case GDK_SCROLL_LEFT:
			return update_zoom(data->x, data->y, ZOOM_IN);
			// Fallthrough intentional
		case GDK_SCROLL_UP:
		case GDK_SCROLL_RIGHT:
			return update_zoom(data->x, data->y, ZOOM_OUT);
		default:
			break;
	}
	return FALSE;
}

scroll_function_metadata scroll_traverses_seq_action = { scroll_traverses_sequence,
	NAME_SCROLL_SEQ_ACTION, N_("Scroll wheel scrolls through the frames in the current sequence"), SCROLL_REF_SEQ };

static double seq_accumulator = 0.0;

gboolean scroll_traverses_sequence(scroll_data *data) {
	if (!sequence_is_loaded())
		return FALSE;
	gint idx = -1;
	point delta;
	switch (data->direction_raw) {
		case GDK_SCROLL_SMOOTH:
			gdk_scroll_event_get_deltas(data->event, &delta.x, &delta.y);
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
	mouse_action_functions_init_statics();
	GtkNotebook* notebook = maf_notebook_center_box;
	gint tab = gtk_notebook_get_current_page(notebook);
	point delta;
	gdk_scroll_event_get_deltas(data->event, &delta.x, &delta.y);
	switch (data->direction_raw) {
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
