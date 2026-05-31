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

#include "core/siril.h"
#include "core/proto.h"
#include "core/command.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/masks.h"
#include "core/siril_update.h"
#include "gui-gtk4/siril_cmd_help.h"
#include "core/siril_log.h"
#include "core/initfile.h"
#include "compositing/align_rgb.h"
#include "io/annotation_catalogues.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "algos/astrometry_solver.h"
#include "algos/noise.h"
#include "algos/photometry.h"
#include "algos/siril_wcs.h"
#include "algos/ccd-inspector.h"
#include "gui-gtk4/ccd-inspector.h"
#include "compositing/compositing.h"
#include "filters/synthstar.h"
#include "gui-gtk4/about_dialog.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/colors.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/curves.h"
#include "gui-gtk4/documentation.h"
#include "gui-gtk4/histogram.h"
#include "gui-gtk4/histo_display.h"
#include "gui-gtk4/icc_profile.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/open_dialog.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/PSF_list.h"
#include "gui-gtk4/save_dialog.h"
#include "gui-gtk4/script_menu.h"
#include "gui-gtk4/python_gui.h"
#include "gui-gtk4/sequence_list.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/photometric_cc.h"
#include "gui-gtk4/menu_gray_geometry.h"
#include "gui-gtk4/registration.h"
#include "gui-gtk4/registration_preview.h"
#include "gui-gtk4/remixer.h"
#include "gui-gtk4/user_polygons.h"
#include "registration/registration.h"
#include "io/siril_catalogues.h"

#include "siril_actions.h"

#define CHECK_FOR_OPENED_DIALOG \
    do { \
        if (is_a_dialog_opened()) { \
            gui_iface.message_dialog(SIRIL_MSG_INFO, _("Cannot process image"), _("The image can't be processed while another processing dialog is opened.")); \
            return; \
        } \
    } while (0)

void open_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	header_open_button_clicked();
	gui_function(launch_clipboard_survey, NULL);
}

void cwd_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	cwd_btton_clicked();
	gui_function(update_MenuItem, NULL);
}

void livestacking_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(gui.builder, "livestacking_player"));

	gtk_widget_set_visible(w, TRUE);
	/* GTK4: gtk_window_set_keep_above removed */;
}

void save_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_save_button_clicked();
}

void save_as_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_save_as_button_clicked();
}

void snapshot_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_snapshot_button_clicked(FALSE);
}

void clipboard_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_snapshot_button_clicked(TRUE);
}

/* Scripts menu — entries appended by initialize_script_menu() reference
 * win.script-getscripts and win.script-pythonpad.  Without these GAction
 * registrations the GMenu items render insensitive (greyed out). */
void script_getscripts_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_get_scripts_clicked(NULL);
}

void script_pythonpad_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_open_pythonpad(NULL, NULL);
}

/* Boolean toggle backing the headerbar Scripts ▸ "Enable Python debug mode"
 * check item.  It shares the python_debug flag with the script editor's
 * Script ▸ debug toggle; set_python_debug_mode() updates the flag and keeps
 * both check items showing the same state. */
void script_pythondebug_change_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	g_simple_action_set_state(action, state);
	set_python_debug_mode(g_variant_get_boolean(state));
}

void undo_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	gui_iface.set_busy(TRUE);
	undo_display_data(UNDO);
	gui_iface.set_busy(FALSE);
}

void redo_action_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	gui_iface.set_busy(TRUE);
	undo_display_data(REDO);
	gui_iface.set_busy(FALSE);
}

void quit_action_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	siril_quit();
}

void about_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_show_about_dialog();
}

void preferences_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("settings_window");
}

void close_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	process_close(0);
}

void updates_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
#if defined(HAVE_LIBCURL)
	siril_check_updates(TRUE);
#else
	siril_log_message(_("Cannot check for updates with this version, missing dependency\n"));
#endif
}

void doc_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_documentation(NULL);
}

static gboolean is_extended = FALSE;

void full_screen_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWindow *window;
	GtkWidget *toolbarbox = GTK_WIDGET(gtk_builder_get_object(gui.builder, "toolbarbox"));
	GtkRevealer *rev = GTK_REVEALER(gtk_builder_get_object(gui.builder, "center_revealer"));
	GtkButton *button = GTK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "button_paned")));
	gboolean is_control_box_visible;
	gboolean is_fullscreen;

	window = GTK_WINDOW(GTK_APPLICATION_WINDOW(user_data));

	is_fullscreen = gtk_window_is_fullscreen(window);
	is_control_box_visible = gtk_widget_get_visible(
	    gtk_paned_get_end_child(GTK_PANED(gtk_builder_get_object(gui.builder, "main_panel"))));

	if (is_fullscreen) {
		gtk_window_unfullscreen(window);
		if (is_extended)
			g_signal_emit_by_name(button, "clicked");
	} else {
		gtk_window_fullscreen(window);
		if (is_control_box_visible) {
			g_signal_emit_by_name(button, "clicked");
		}
		is_extended = is_control_box_visible;
	}
	gtk_widget_set_visible(toolbarbox, is_fullscreen);
}

void panel_animate(gboolean show);   /* defined in callbacks.c */

void panel_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkPaned *paned = GTK_PANED(gtk_builder_get_object(gui.builder, "main_panel"));
	GtkImage *image = GTK_IMAGE(gtk_button_get_child(GTK_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "button_paned")))));

	gboolean is_visible = gtk_widget_get_visible(gtk_paned_get_end_child(paned));

	panel_animate(!is_visible);
	gtk_image_set_from_icon_name(image, !is_visible ? "pan-end-symbolic" : "pan-start-symbolic");

	if (!is_visible && com.gui_icc.iso12646)
		disable_iso12646_conditions(TRUE, FALSE, TRUE);

	if (com.pref.gui.remember_windows) {
		com.pref.gui.is_extended = !is_visible;
		writeinitfile();
	}
}

void keyboard_shortcuts_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWindow *window;

	window = GTK_WINDOW(GTK_APPLICATION_WINDOW(user_data));

	siril_cmd_help_keyboard_shortcuts(window);
}

void tab_conversion_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(FILE_CONVERSION);
}

void tab_sequence_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(IMAGE_SEQ);
}

void tab_prepro_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PRE_PROC);
}

void tab_registration_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(REGISTRATION);
}

void tab_plot_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PLOT);
}

void tab_stacking_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

void tab_logs_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(OUTPUT_LOGS);
}

void toolbar_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(gui.builder, "toolbarbox"));
	gtk_widget_set_visible(w, !gtk_widget_get_visible(w));
}

void on_histogram_overlay_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	gboolean new_state = g_variant_get_boolean(state);
	set_histogram_overlay_visible(new_state);
	g_simple_action_set_state(action, state);
}

void on_histogram_overlay_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;
	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

/* Selection-aspect-ratio radio group (replaces the GTK3 .ui-driven
 * GtkRadioMenuItem set in the menugray "Selection" submenu).  The action
 * state carries a string token: "free", "preserve", "16:9", "3:2", "4:3",
 * "1:1", "3:4", "2:3", "9:16".  Each menu item references the action with
 * the corresponding ::target string and GMenu renders them as a radio
 * group automatically. */
static double aspect_ratio_for_token(const char *tok) {
	if (g_strcmp0(tok, "free") == 0)     return 0.0;
	if (g_strcmp0(tok, "preserve") == 0) {
		return gfit ? (double)gfit->rx / (double)gfit->ry : 0.0;
	}
	if (g_strcmp0(tok, "16:9") == 0) return 16.0 / 9.0;
	if (g_strcmp0(tok, "3:2") == 0)  return  3.0 / 2.0;
	if (g_strcmp0(tok, "4:3") == 0)  return  4.0 / 3.0;
	if (g_strcmp0(tok, "1:1") == 0)  return  1.0;
	if (g_strcmp0(tok, "3:4") == 0)  return  3.0 / 4.0;
	if (g_strcmp0(tok, "2:3") == 0)  return  2.0 / 3.0;
	if (g_strcmp0(tok, "9:16") == 0) return  9.0 / 16.0;
	return 0.0;
}

void aspect_ratio_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	g_action_change_state(G_ACTION(action), parameter);
}

void aspect_ratio_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	const char *tok = g_variant_get_string(state, NULL);
	double ratio = aspect_ratio_for_token(tok);
	gui.ratio = ratio;
	if (ratio != 0.0) {
		enforce_ratio_and_clamp();
		update_display_selection();
		gui_function(new_selection_zone, NULL);
		gui_iface.redraw_image(REDRAW_OVERLAY);
	}
	g_simple_action_set_state(action, state);
}

void selection_guides_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	g_action_change_state(G_ACTION(action), parameter);
}

void selection_guides_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	const char *tok = g_variant_get_string(state, NULL);
	com.pref.gui.selection_guides = atoi(tok);
	gui_iface.redraw_image(REDRAW_OVERLAY);
	g_simple_action_set_state(action, state);
}

void select_all_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!gfit) return;
	com.selection.x = 0;
	com.selection.y = 0;
	com.selection.w = gfit->rx;
	com.selection.h = gfit->ry;
	double original = (double)gfit->rx / (double)gfit->ry;
	if (gui.ratio != original) {
		gui.ratio = 0.0;
	} else {
		gui.ratio = original;
		enforce_ratio_and_clamp();
	}
	update_display_selection();
	gui_function(new_selection_zone, NULL);
	gui_iface.redraw_image(REDRAW_OVERLAY);
}

/* Public sync helpers used by mouse_action_functions.c::show_drawingarea_popover
 * to refresh the radio state from gui.ratio / com.pref.gui.selection_guides
 * just before popping menugray, mirroring the GTK3 sync that happened in
 * do_popup_graymenu_unused().  This keeps the radio in step with state that
 * external code paths (e.g. keyboard shortcuts, registration UI) may have
 * mutated since the last menu open. */
static const char *aspect_token_for_ratio(double r) {
	if (r == 0.0) return "free";
	if (gfit && r == (double)gfit->rx / (double)gfit->ry) return "preserve";
	if (r == 16.0 / 9.0) return "16:9";
	if (r ==  3.0 / 2.0) return "3:2";
	if (r ==  4.0 / 3.0) return "4:3";
	if (r ==  1.0)       return "1:1";
	if (r ==  3.0 / 4.0) return "3:4";
	if (r ==  2.0 / 3.0) return "2:3";
	if (r ==  9.0 / 16.0) return "9:16";
	return "free";
}

void sync_selection_action_state(gpointer window) {
	if (!window || !G_IS_ACTION_MAP(window)) return;
	GAction *ar = g_action_map_lookup_action(G_ACTION_MAP(window), "aspect-ratio");
	if (ar) g_simple_action_set_state(G_SIMPLE_ACTION(ar),
		g_variant_new_string(aspect_token_for_ratio(gui.ratio)));
	GAction *sg = g_action_map_lookup_action(G_ACTION_MAP(window), "selection-guides");
	if (sg) {
		char tok[8];
		g_snprintf(tok, sizeof tok, "%d", com.pref.gui.selection_guides);
		g_simple_action_set_state(G_SIMPLE_ACTION(sg), g_variant_new_string(tok));
	}
}

void change_zoom_fit_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	if (g_variant_get_boolean(state)) {
		gui.zoom_value = ZOOM_FIT;
		if (com.gui_icc.iso12646)
			disable_iso12646_conditions(FALSE, TRUE, TRUE);
		reset_display_offset();
		gui_iface.redraw_image(REDRAW_IMAGE);
	} else {
		gui.zoom_value = get_zoom_val();
	}
	update_zoom_label();
	g_simple_action_set_state(action, state);
}

void zoom_fit_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void zoom_in_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	point center = get_center_of_vport();
	update_zoom(center.x, center.y, ZOOM_IN);
	if (com.gui_icc.iso12646)
		disable_iso12646_conditions(FALSE, TRUE, TRUE);
}

void zoom_out_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	point center = get_center_of_vport();
	update_zoom(center.x, center.y, ZOOM_OUT);
	if (com.gui_icc.iso12646)
		disable_iso12646_conditions(FALSE, TRUE, TRUE);
}

void zoom_one_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	update_zoom_fit_button();
	gui.zoom_value = ZOOM_NONE;
	reset_display_offset();
	update_zoom_label();
	gui_iface.redraw_image(REDRAW_IMAGE);
}

void negative_view_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	g_simple_action_set_state(action, state);
	gui_iface.set_busy(TRUE);
	notify_gfit_data_modified(); // here the data isn't modified but we need to trigger the remap
	gui_iface.redraw_image(REDRAW_ALL);
	gui_function(redraw_previews, NULL);
	gui_iface.set_busy(FALSE);
}

void negative_view_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void photometry_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	mouse_status = g_variant_get_boolean(state) ? MOUSE_ACTION_PHOTOMETRY : MOUSE_ACTION_SELECT_REG_AREA;
	g_simple_action_set_state(action, state);
	free(gui.qphot);
	gui.qphot = NULL;
	gui_iface.redraw_image(REDRAW_OVERLAY);
}

void photometry_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void color_map_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	g_simple_action_set_state(action, state);
	gui_iface.set_busy(TRUE);
	notify_gfit_data_modified();
	gui_iface.redraw_image(REDRAW_ALL);
	gui_function(redraw_previews, NULL);
	gui_iface.set_busy(FALSE);
}

void color_map_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

static void update_chain_channels_ui(gboolean linked) {
	set_unlink_channels(!linked);

	const gchar *path = linked ? "/org/siril/ui/pixmaps/chain-linked.svg"
	                           : "/org/siril/ui/pixmaps/chain.svg";
	GtkImage *image = GTK_IMAGE(gtk_builder_get_object(gui.builder, "autostretch_linked_icon"));
	if (image) {
		/* Use the GtkIconPaintable / paintable pipeline rather than
		 * gtk_image_set_from_resource: the latter goes through GdkPixbuf
		 * for SVG which silently fails on systems without the librsvg
		 * pixbuf loader, leaving the chain icon stuck on its initial
		 * .ui-loaded paintable.  siril_paintable_from_resource uses
		 * GTK4's native SVG renderer. */
		GdkPaintable *paintable = siril_paintable_from_resource(path, 24);
		if (paintable) {
			gtk_image_set_from_paintable(image, paintable);
			g_object_unref(paintable);
		} else {
			gtk_image_set_from_resource(image, path);
		}
	}

	GtkWidget *button = GTK_WIDGET(gtk_builder_get_object(gui.builder, "linked_autostretch_button"));
	gchar *tooltip_text = g_strdup_printf(_("Link/unlink channels in autostretch viewer mode.\nCurrent state: %s."),
										linked ? _("linked") : _("unlinked"));
	gtk_widget_set_tooltip_text(button, tooltip_text);

	g_free(tooltip_text);

	/* Toggling channel-linking in autostretch mode changes the lookup
	 * table per channel, so the displayed image needs to be remapped —
	 * mirroring the path on_display_item_toggled() takes when the
	 * preview-mode dropdown changes. */
	if (single_image_is_loaded() || sequence_is_loaded()) {
		notify_gfit_data_modified();
		redraw(REDRAW_ALL);
		gui_function(redraw_previews, NULL);
	}
}

void chain_channels_state_change(GSimpleAction *action, GVariant *state, gpointer user_data) {
	gboolean linked = g_variant_get_boolean(state);
	g_simple_action_set_state(action, state);
	update_chain_channels_ui(linked);
}

gboolean chain_channels_idle_callback(gpointer user_data) {
	gboolean linked = GPOINTER_TO_INT(user_data);
	update_chain_channels_ui(linked);
	return G_SOURCE_REMOVE;
}

void chain_channels_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void astrometry_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	open_astrometry_dialog();
}

void dyn_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("stars_list_window");
}

void pick_star_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	pick_a_star();
}

void psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	psf_star *result = NULL;
	int layer = select_vport(gui.cvport);

	if (layer == -1)
		return;
	if (!(com.selection.h && com.selection.w))
		return;
	if (com.selection.w > 300 || com.selection.h > 300) {
		gui_iface.message_dialog(SIRIL_MSG_WARNING, _("Current selection is too large"),
				_("To determine the PSF, please make a selection around a star."));
		return;
	}
	struct phot_config *ps = phot_set_adjusted_for_image(gfit);
	psf_error error = PSF_NO_ERR;
	result = psf_get_minimisation(gfit, layer, &com.selection, TRUE, FALSE, ps, TRUE, com.pref.starfinder_conf.profile, &error);
	free(ps);
	if (result)
		popup_psf_result(result, &com.selection, gfit);
	free_psf(result);
}

void seq_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	if (!sequence_is_loaded()) {
		siril_log_error(_("Error: no sequence loaded.\n"));
		return;
	}

	// If we reach here, sequence is loaded
	sequence *seq = &com.seq;
	int layer = select_vport(gui.cvport);

	// Validate selection size
	if (com.selection.w > 300 || com.selection.h > 300) {
		siril_log_error(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
		return;
	}
	if (com.selection.w <= 0 || com.selection.h <= 0) {
		siril_log_error(_("Select an area first\n"));
		return;
	}

	// Determine framing mode
	framing_mode framing = REGISTERED_FRAME;
	if (!seq->regparam[layer])
		framing = ORIGINAL_FRAME;
	if (framing == ORIGINAL_FRAME) {
		// com.headless is FALSE, so we can read the toggle directly
		if (registration_get_follow_star())
			framing = FOLLOW_STAR_FRAME;
	}

	// Run PSF
	siril_log_message(_("Running the PSF on the sequence, layer %d\n"), layer);
	int retval = seqpsf(seq, layer, FALSE, FALSE, FALSE, framing, TRUE, FALSE) ? 1 : 0;

	if (retval != 0) {
		siril_log_error(_("Error running the PSF on the sequence\n"));
	}
}

void crop_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_crop();
}

void seq_crop_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (valid_rgbcomp_seq()) {
		crop_rgbcomp_seq();
	} else {
		siril_open_dialog("crop_dialog");
	}
}

void annotate_dialog_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("annotate_dialog");
}

void annotate_object_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	if (g_variant_get_boolean(state)) {
		if (has_wcs(gfit)) {
			com.found_object = find_objects_in_field(gfit);
		}
	} else {
		clear_user_polygons();
		g_slist_free(com.found_object);
		com.found_object = NULL;
		purge_user_catalogue(CAT_AN_USER_TEMP);
		refresh_annotation_visibility();
	}
	g_simple_action_set_state(action, state);
	gui_iface.redraw_image(REDRAW_OVERLAY);
}

void wcs_grid_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	gui.show_wcs_grid = g_variant_get_boolean(state);
	g_simple_action_set_state(action, state);
	gui_iface.redraw_image(REDRAW_OVERLAY);
}

void annotate_object_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void wcs_grid_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void regframe_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	GtkToggleButton *drawframe;
	drawframe = GTK_TOGGLE_BUTTON(GTK_WIDGET(gtk_builder_get_object(gui.builder, "drawframe_check")));
	siril_toggle_set_active(GTK_WIDGET(drawframe), g_variant_get_boolean(state));
	g_simple_action_set_state(action, state);
	gui_iface.redraw_image(REDRAW_OVERLAY);
}

void regframe_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;
	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void seq_list_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_widget_get_visible(GTK_WIDGET(gtk_builder_get_object(gui.builder, "seqlist_dialog")))) {
		siril_close_dialog("seqlist_dialog");
	} else {
		gboolean confirm = TRUE;
		if (com.seq.current == RESULT_IMAGE) {
			confirm = gui_iface.confirm_dialog(_("Save your changes before loading a frame of the sequence."),
					_("The image currently displayed is the result of the previous stack. "
							"If you load an image from the sequence, you might lose the entire process you performed on the image, "
							"but not the image itself. You need to save your data before doing this."),
					_("Load another image"));
		}
		if (confirm) {
			int layer = get_registration_layer(&com.seq);
			update_seqlist(layer);
			siril_open_dialog("seqlist_dialog");
		}
	}
}

void statistics_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	gui_iface.set_busy(TRUE);
	computeStat();
	siril_open_dialog("StatWindow");
	gui_iface.set_busy(FALSE);
}

void noise_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	evaluate_noise_in_image();
}

void ccd_inspector_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	compute_aberration_inspector();
}

void show_tilt_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void show_tilt_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	gui_iface.set_busy(TRUE);
	if (g_variant_get_boolean(state)) {
		draw_sensor_tilt(gfit);

	} else {
		clear_sensor_tilt();
		gui_iface.redraw_image(REDRAW_OVERLAY);
	}
	g_simple_action_set_state(action, state);
	gui_iface.set_busy(FALSE);
}

void show_disto_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void show_disto_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	if (!has_wcs(gfit) || !gfit->keywords.wcslib->lin.dispre) {
		siril_log_error(_("This command only works on plate solved images with distortions included\n"));
		return;
	}
	gui_iface.set_busy(TRUE);

	gui.show_wcs_disto = g_variant_get_boolean(state);
	gui_iface.redraw_image(REDRAW_OVERLAY);
	g_simple_action_set_state(action, state);

	gui_iface.set_busy(FALSE);
}

void image_information_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("file_information");
}

void image_fits_header_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("keywords_dialog");
}

/******* processing menu **************/

void remove_green_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("SCNR_dialog");
}

void saturation_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (value_check( gfit))
		siril_open_dialog("satu_dialog");
}

void color_calib_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	initialize_calibration_interface();
	siril_open_dialog("color_calibration");
}

void pcc_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	initialize_photometric_cc_dialog();
	siril_open_dialog("s_pcc_dialog");
}

void spcc_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	initialize_spectrophotometric_cc_dialog();
	siril_open_dialog("s_pcc_dialog");
}

void split_channel_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	siril_open_dialog("extract_channel_dialog");
}

void negative_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (value_check(gfit)) {
		CHECK_FOR_OPENED_DIALOG;
		negative_processing();
	}
}

void histo_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (value_check(gfit))
		toggle_histogram_window_visibility(1);
}

void curves_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (value_check(gfit))
		toggle_curves_window_visibility();
}

/* Forward declarations for per-dialog init_statics functions that cache
 * widget pointers / attach handlers.  Each must run before its dialog
 * is shown — previously they only ran from the dialog's Apply path,
 * leaving widget pointers NULL and (for fft / linear_match) the file
 * chooser buttons inert until the user clicked Apply. */
extern void linear_match_init_statics(void);
extern void fft_dialog_init_statics(void);
extern void banding_dialog_init_statics(void);
extern void cosmetic_dialog_init_statics(void);
extern void rgradient_dialog_init_statics(void);
extern void split_cfa_init_statics(void);

void fix_banding_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	banding_dialog_init_statics();
	siril_open_dialog("canon_fixbanding_dialog");
}

void cosmetic_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	cosmetic_dialog_init_statics();
	siril_open_dialog("cosmetic_dialog");
}

void background_extr_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("background_extraction_dialog");
}

void asinh_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (value_check(gfit))
		siril_open_dialog("asinh_dialog");
}

void epf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("epf_dialog");
}

void unpurple_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
        siril_open_dialog("unpurple_dialog");
}

void starnet_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("starnet_dialog");
}

void deconvolution_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("bdeconv_dialog");
}

void payne_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (value_check(gfit))
		toggle_histogram_window_visibility(2);
}

void binning_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("binxy_dialog");
}

void resample_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("resample_dialog");
}

void rotation_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (com.selection.w == 0 || com.selection.h == 0) {
		com.selection = (rectangle){ 0, 0, gfit->rx, gfit->ry };
	}
	siril_open_dialog("rotation_dialog");
	gui_iface.redraw_image(REDRAW_OVERLAY);
}

void rotation90_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	CHECK_FOR_OPENED_DIALOG;
	siril_rotate90();
}

void rotation270_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	CHECK_FOR_OPENED_DIALOG;
	siril_rotate270();
}

void mirrorx_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	CHECK_FOR_OPENED_DIALOG;
	mirrorx_gui(gfit);
}

void mirrory_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	CHECK_FOR_OPENED_DIALOG;
	mirrory_gui(gfit);
}

void wavelets_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	siril_open_dialog("wavelets_dialog");
}

void split_wavelets_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("extract_wavelets_layers_dialog");
}

void medianfilter_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("Median_dialog");
}

void rgradient_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	rgradient_dialog_init_statics();
	siril_open_dialog("rgradient_dialog");
}

void clahe_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("CLAHE_dialog");
}

void linearmatch_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	linear_match_init_statics();
	siril_open_dialog("linearmatch_dialog");
}

void fft_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	fft_dialog_init_statics();
	/* Seed the "initial folder" hint that the dialog will use when the
	 * button is clicked without a prior selection. */
	GtkWidget *magbutton = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filechooser_mag"));
	GtkWidget *phasebutton = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filechooser_phase"));
	siril_file_chooser_set_current_folder_path(magbutton, com.wd);
	siril_file_chooser_set_current_folder_path(phasebutton, com.wd);
	siril_open_dialog("dialog_FFT");
}

void rgb_compositing_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (single_image_is_loaded() && gui_iface.confirm_dialog(_("Close current image?"), _("Opening the RGB Composition dialog will close the current image without saving. Are you sure?"), _("Yes"))) {
		close_single_image();
	}
	close_sequence(FALSE);
	open_compositing_window();
}

void star_remix_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (single_image_is_loaded() && gui_iface.confirm_dialog(_("Close current image?"), _("Opening the Star Recomposition dialog will close the current image without saving. Are you sure?"), _("Yes"))) {
		close_single_image();
	}
	close_sequence(FALSE);
	toggle_remixer_window_visibility(CALL_FROM_MENU, NULL, NULL);
}

void pixel_math_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("pixel_math_dialog");
}

void split_cfa_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	split_cfa_init_statics();
	siril_open_dialog("split_cfa_dialog");
}

void nina_lc_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("nina_light_curve");
}

void compstars_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("compstars");
}

void catmag_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("catmag");
}

void denoise_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("denoise_dialog");
}

void merge_cfa_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("merge_cfa_dialog");
}

void star_desaturate_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	CHECK_FOR_OPENED_DIALOG;
	if (!check_ok_if_cfa()) return;
	control_window_switch_to_tab(OUTPUT_LOGS);
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = gfit;
	args->image_hook = unclip_image_hook;
	args->log_hook = unclip_log_hook;
	args->description = _("Unclip stars");
	args->verbose = TRUE;
	if (!start_in_new_thread(generic_image_worker, args))
		free_generic_img_args(args);
}

void star_synthetic_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	CHECK_FOR_OPENED_DIALOG;
	if (!check_ok_if_cfa())
		return;
	control_window_switch_to_tab(OUTPUT_LOGS);
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = gfit;
	args->image_hook = synthstar_image_hook;
	args->log_hook = synthstar_log_hook;
	args->description = _("Synthetic stars");
	args->verbose = TRUE;
	if (!start_in_new_thread(generic_image_worker, args))
		free_generic_img_args(args);
}

void align_dft_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	undo_save_state(gfit, _("RGB alignment (DFT)"));
	rgb_align(1);
}

void align_global_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	undo_save_state(gfit, _("RGB alignment (Global stars)"));
	rgb_align(2);
}

void align_kombat_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	undo_save_state(gfit, _("RGB alignment (KOMBAT)"));
	rgb_align(3);
}

void align_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	undo_save_state(gfit, _("RGB alignment (PSF)"));
	if (com.selection.w > 300 || com.selection.h > 300) {
		siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
		return;
	}
	rgb_align(0);
}

void icc_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("icc_dialog");
}

/* Stateful boolean action backing the bottom-toolbar intensity profile
 * toggle.  The old non-stateful handler read the GtkToggleButton's
 * active property to decide whether to open the dialog, but in GTK4 a
 * non-stateful action doesn't drive the button's active state — the
 * read returned FALSE on every click and the dialog never opened.
 * Following the photometry_state / photometry_activate pattern keeps
 * the button and the dialog in sync. */
void cut_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	g_simple_action_set_state(action, state);
	if (g_variant_get_boolean(state)) {
		mouse_status = MOUSE_ACTION_CUT_SELECT;
		siril_open_dialog("cut_dialog");
	} else {
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		siril_close_dialog("cut_coords_dialog");
		siril_close_dialog("cut_spectroscopy_dialog");
		siril_close_dialog("cut_dialog");
	}
}

void cut_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action),
		g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void clear_roi(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_clear_roi();
}

void set_roi(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_set_roi();
}

void ccm_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("ccm_dialog");
}

void mask_from_image_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	mask_from_image_dialog_set_file_mode(FALSE);
	siril_open_dialog("mask_from_image_dialog");
}

void mask_from_stars_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("mask_from_stars_dialog");
}

void mask_from_color_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("mask_from_color_dialog");
}

void mask_from_file_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	mask_from_image_dialog_set_file_mode(TRUE);
	siril_open_dialog("mask_from_image_dialog");
}

void clear_mask_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!gfit || !gfit->mask) {
		return;
	}

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit = gfit;
	args->mask_hook = mask_clear_hook;
	args->description = _("Clear mask");
	args->verbose = TRUE;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}

void autostretch_mask_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!gfit || !gfit->mask || !gfit->mask->data) {
		gui_iface.message_dialog(SIRIL_MSG_ERROR, _("No mask present"),
		                     _("There is no mask to autostretch."));
		return;
	}

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit = gfit;
	args->mem_ratio = 1.0f;
	args->mask_hook = mask_autostretch_hook;
	args->description = _("Autostretch mask");
	args->verbose = TRUE;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}

void blur_mask_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("mask_blur_dialog");
}

void feather_mask_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("mask_feather_dialog");
}

void mask_scale_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("mask_scale_dialog");
}

void invert_mask_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!gfit || !gfit->mask || !gfit->mask->data) {
		gui_iface.message_dialog(SIRIL_MSG_ERROR, _("No mask present"),
		                     _("There is no mask to invert."));
		return;
	}

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit = gfit;
	args->mem_ratio = 0.0f;
	args->mask_hook = mask_invert_hook;
	args->log_hook = NULL;
	args->description = _("Invert mask");
	args->verbose = TRUE;
	args->user = NULL;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}

void mask_add_from_poly_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	init_add_poly_to_mask();
}

void mask_clear_from_poly_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	init_clear_poly_from_mask();
}

void mask_from_gradient_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!gfit) {
		return;
	}

	struct generic_mask_args *args = calloc(1, sizeof(struct generic_mask_args));
	args->fit = gfit;
	args->mask_hook = mask_from_gradient_hook;
	args->description = _("Gradient of mask");
	args->verbose = TRUE;
	args->max_threads = com.max_thread;

	start_in_new_thread(generic_mask_worker, args);
}

void threshold_mask_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("mask_thresholds_dialog");
}
