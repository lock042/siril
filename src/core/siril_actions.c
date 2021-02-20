/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/command.h"
#include "core/undo.h"
#include "core/siril_update.h"
#include "core/siril_cmd_help.h"
#include "algos/colors.h"
#include "algos/noise.h"
#include "algos/siril_wcs.h"
#include "algos/plateSolver.h"
#include "gui/about_dialog.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/open_dialog.h"
#include "gui/save_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/script_menu.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/photometric_cc.h"

#include "siril_actions.h"

void open_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	header_open_button_clicked();
}

void cwd_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	cwd_btton_clicked();
}

void save_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	on_header_save_button_clicked();
}

void save_as_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	on_header_save_as_button_clicked();
}

void snapshot_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	on_header_snapshot_button_clicked();
}

void undo_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	undo_display_data(UNDO);
	set_cursor_waiting(FALSE);
}

void redo_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	undo_display_data(REDO);
	set_cursor_waiting(FALSE);
}

void quit_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_quit();
}

void about_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_show_about_dialog();
}

void preferences_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("settings_window");
}

void close_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	process_close(0);
}

void scripts_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_get_on_script_pages();
}

void updates_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_check_updates(TRUE);
}

static gboolean is_extended = FALSE;

void full_screen_activated(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	GtkWindow *window;
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	GtkWidget *control_center_box = lookup_widget("control_center_box");
	GtkButton *button = GTK_BUTTON(lookup_widget("button_paned"));
	gboolean is_control_box_visible;
	gboolean is_fullscreen;

	window = GTK_WINDOW(GTK_APPLICATION_WINDOW(user_data));

	GdkWindow *gdk_window = gtk_widget_get_window(GTK_WIDGET(window));
	is_fullscreen = gdk_window_get_state(gdk_window) & GDK_WINDOW_STATE_FULLSCREEN;
	is_control_box_visible = gtk_widget_get_visible(control_center_box);

	if (is_fullscreen) {
		gtk_window_unfullscreen(window);
		if (is_extended)
			gtk_button_clicked(button);
	} else {
		gtk_window_fullscreen(window);
		if (is_control_box_visible) {
			gtk_button_clicked(button);
		}
		is_extended = is_control_box_visible;
	}
	gtk_widget_set_visible(toolbarbox, is_fullscreen);
}

void keyboard_shortcuts_activated(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	GtkWindow *window;

	window = GTK_WINDOW(GTK_APPLICATION_WINDOW(user_data));

	siril_cmd_help_keyboard_shortcuts(window);
}

void tab_conversion_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(FILE_CONVERSION);
}

void tab_sequence_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(IMAGE_SEQ);
}

void tab_prepro_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PRE_PROC);
}

void tab_registration_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(REGISTRATION);
}

void tab_plot_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PLOT);
}

void tab_stacking_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

void tab_logs_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(OUTPUT_LOGS);
}

void toolbar_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	GtkWidget *w = lookup_widget("toolbarbox");
	gtk_widget_set_visible(w, !gtk_widget_get_visible(w));
}

void change_zoom_fit_state(GSimpleAction *action, GVariant *state,
		gpointer user_data) {
	if (g_variant_get_boolean(state)) {
		com.zoom_value = ZOOM_FIT;
		reset_display_offset();
		redraw(com.cvport, REMAP_NONE);
	} else {
		com.zoom_value = get_zoom_val();
	}
	g_simple_action_set_state(action, state);
}

void zoom_fit_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action),
			g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void zoom_in_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	point center = get_center_of_vport();
	update_zoom(center.x, center.y, ZOOM_IN);
}

void zoom_out_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	point center = get_center_of_vport();
	update_zoom(center.x, center.y, ZOOM_OUT);
}

void astrometry_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	open_astrometry_dialog();
}

void dyn_psf_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	printf("OKOK\n");
	siril_open_dialog("stars_list_window");
}

void search_object_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	if (has_wcs(&gfit))
		siril_open_dialog("search_objects");
}

void statistics_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	computeStat();
	siril_open_dialog("StatWindow");
	set_cursor_waiting(FALSE);
}

void noise_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	evaluate_noise_in_image();
}

void image_information_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("file_information");
}

/******* processing menu **************/

void remove_green_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("SCNR_dialog");
}

void saturation_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("satu_dialog");
}

void color_calib_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	initialize_calibration_interface();
	siril_open_dialog("color_calibration");
}

void pcc_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	initialize_photometric_cc_dialog();
	siril_open_dialog("ImagePlateSolver_Dial");
}

void split_channel_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("extract_channel_dialog");
}

void negative_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	negative_processing();
}

void histo_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("histogram_dialog");
}

void fix_banding_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("canon_fixbanding_dialog");
}

void cosmetic_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("cosmetic_dialog");
}

void background_extr_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("background_extraction_dialog");

}
