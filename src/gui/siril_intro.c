/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include "siril_intro.h"

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"

#define INTRO_DELAY 1000

static guint tip_index;
static gboolean go_next;

const SirilTipIntro intro_tips[] = {
		{"headerbar", N_("Welcome to the new Siril version. Please take a moment to read tips about this release"), 5},
		{"notebook1", N_("All the application windows have been merged into this window. In the left panel, you can see the image preview with the Red, Green, Blue channels and the RGB mix"), 6},
		{"label22", N_("Pre-processing steps are grouped in the right panel. You can reach each step with the F1, F2, ... F7 keys"), 5},
		{"right_panel_image", N_("Hitting this button will hide the right panel. You can also try the full screen mode (Control-F)"), 4},
		{"hamburger-menu", N_("Press F10 or click on this button to open the menu. You can find here the shortcut list and the preferences dialogs where many options are available"), 6},
		{"cwd_button", N_("You can now change your working directory by hitting this button. The working directory is shown right below the title at the center of the headerbar"), 6},
		{"header_open_button", N_("You can open a single image or FITS/SER sequence"), 3},
		{"recent_menu_button", N_("Here are listed the most recent FITS file you have opened"), 3},
		{"header_processing_button", N_("Processing algorithms are all in this single menu"), 3},
		{"header_undo_button", N_("Use this button to undo an operation"), 3},
		{"header_redo_button", N_("This one to redo an operation"), 3},
		{"header_precision_button", N_("Now Siril works in the 32-bit per chanel precision by default. You can change it in Preferences and you can change the currently loaded single image precision too with this selector"), 8},
		{"header_save_as_button", N_("Save your work as many times as needed by choosing a new name ..."), 3},
		{"header_save_button", N_("... or save the current FITS image with the same name"), 3},
		{"command", N_("As usual you can enter Siril commands. To have an overview of all commands, type \"help\""), 4},
		{"GtkToolMainBar", N_("Basic viewing operations are available in the main toolbar. Zooming is now available through Control-Wheel"), 5},
		{"drawingarear", N_("Enjoy using the new Siril"), 6}
};

static GtkWidget *intro_popover(GtkWidget *widget, const gchar *text) {
	gchar *markup_txt = g_strdup_printf("<big><b>%s</b></big>", text);
	GtkWidget *popover = popover_new(widget, markup_txt);
#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_popover_popup(GTK_POPOVER(popover));
#else
	gtk_widget_show(popover);
#endif
	g_free(markup_txt);
	return popover;
}

static gboolean intro_popover_close(gpointer user_data) {
	gtk_widget_hide((GtkWidget *) user_data);
	go_next = TRUE;
	return FALSE;
}

static gboolean update_popover(gpointer user_data) {
	GtkWidget *popover;
	if (go_next) {
		popover = intro_popover(lookup_widget(intro_tips[tip_index].widget), _(intro_tips[tip_index].tip));
		g_timeout_add(INTRO_DELAY * intro_tips[tip_index].delay, (GSourceFunc) intro_popover_close, (gpointer) popover);
		go_next = FALSE;
		tip_index++;
	}
	return tip_index < G_N_ELEMENTS(intro_tips);
}

static void intro_notify_update(gpointer user_data) {
	g_timeout_add(250, (GSourceFunc) update_popover, user_data);
}

void start_intro_script() {
	tip_index = 0;
	go_next = TRUE;
	intro_notify_update(NULL);
}
