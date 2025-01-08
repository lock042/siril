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

#include "core/siril.h"
#include "core/proto.h"
#include "gui/utils.h"
#include "gui/utils.h"

#include "siril_intro.h"

#define INTRO_DELAY 1000

static guint tip_index;
static gboolean go_next;

const SirilTipIntro intro_tips[] = {
		{"headerbar", N_("Welcome to the newest version of Siril. Please take a moment to read some tips about this release"), 8},
		{"labelRGB", N_("The RGB tab is now ready for processing. You can select and work on this tab which is open by default if the image is RGB."), 8},
		{"label22", N_("Pre-processing steps are grouped together in the right panel. You can reach each step with the F1...F7 keys"), 8},
		{"button_paned", N_("This button will hide the right panel. You can also try the full screen mode (Control - F)"), 8},
		{"hamburger-menu", N_("Press F10 or click on this button to open the menu. Here you can find the shortcut list and the preferences dialog where many of the options are available"), 9},
		{"hamburger-menu", N_("You can get more scripts by clicking on the \"Get scripts\" item in the menu"), 6},
		{"hamburger-menu", N_("The documentation set is also available via the menu \"Siril Manual\""), 6},
		{"cwd_button", N_("You can change your working directory by hitting this button. The working directory is shown right below the title at the center of the headerbar"), 9},
		{"header_open_button", N_("You can open a single image or FITS/SER sequence"), 6},
		{"recent_menu_button", N_("Here is a list of the most recent FITS files you’ve opened"), 6},
		{"header_processing_button", N_("Processing algorithms are all in this single menu"), 6},
		{"header_tools_button", N_("This new button brings together all the tools that were previously scattered throughout the application. It includes some of what was contained in the Hamburger menu, as well as tools for photometry, astrometry and much more. This menu is generally essential once pre-processing is complete."), 15},
		{"header_undo_button", N_("Use this button to undo an operation"), 5},
		{"header_redo_button", N_("Use this button to redo an operation"), 5},
		{"header_precision_button", N_("Siril works in 32-bit per channel precision by default. You can change it in Preferences and you can change the currently loaded image precision with this selector"), 11},
		{"header_save_as_button", N_("Save your work as many times as needed by choosing a new name ..."), 6},
		{"header_save_button", N_("... or save the current FITS image with the same name"), 6},
		{"header_snapshot_button", N_("Take a snapshot of your image at any time, with false color, in negative view or with annotation"), 8},
		{"command", N_("As usual you can enter Siril commands. To have an overview of all commands, type \"help\""), 7},
		{"GtkToolMainBar", N_("Basic viewing operations are available in the main toolbar. Zooming is available with Ctrl-Scroll up and down"), 8},
		{"GtkToolMainBar", N_("To pan your image you can use Ctrl-Click (and drag)"), 6},
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
	SirilUIIntro *ui = (SirilUIIntro *) user_data;

	gtk_widget_hide(ui->popover);
#ifndef OS_OSX // very slow on macOS
	gtk_style_context_remove_class(gtk_widget_get_style_context(ui->widget), "siril-intro-highlight");
#endif
	g_free(ui);
	go_next = TRUE;
	return FALSE;
}

static gboolean intro_popover_update(gpointer user_data) {
	if (go_next) {
		SirilUIIntro *ui = g_new(SirilUIIntro, 1);
		ui->widget = lookup_widget(intro_tips[tip_index].widget);
#ifndef OS_OSX // very slow on macOS
		gtk_style_context_add_class(gtk_widget_get_style_context(ui->widget), "siril-intro-highlight");
#endif
		ui->popover = intro_popover(ui->widget, _(intro_tips[tip_index].tip));
		g_timeout_add(INTRO_DELAY * intro_tips[tip_index].delay, (GSourceFunc) intro_popover_close, (gpointer) ui);
		go_next = FALSE;
		tip_index++;
	}
	return tip_index < G_N_ELEMENTS(intro_tips);
}

static void intro_notify_update(gpointer user_data) {
	g_timeout_add(250, (GSourceFunc) intro_popover_update, user_data);
}

void start_intro_script() {
	tip_index = 0;
	go_next = TRUE;
	intro_notify_update(NULL);
}
