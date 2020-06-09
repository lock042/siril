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

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"

#include "siril_intro.h"

#define PREVIEW_DELAY 5000

static guint timer_id;
static guint tip_index;
static gboolean go_next;

const SirilTipHelper help_tips[] = {
		{"headerbar", N_("Welcome to the new Siril version. Please take a moment to read some tips about this release")},
		{"hamburger-menu", N_("Press F10 or click on this button to open the menu and show the shortcut help window")},
		{"cwd_button", N_("You can now change your working directory by hitting this button")},
};

static void free_struct(gpointer user_data) {
	timer_id = 0;
}

static GtkWidget *helper_popover(GtkWidget *widget, const gchar *text) {
	GtkWidget *popover = popover_new(widget, text);
#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_popover_popup(GTK_POPOVER(popover));
#else
	gtk_widget_show(popover);
#endif
	return popover;
}

static gboolean helper_popover_close(gpointer user_data) {
	gtk_widget_hide((GtkWidget *) user_data);
	go_next = TRUE;
	return FALSE;
}

static gboolean update_popover(gpointer user_data) {
	GtkWidget *popover;
	if (go_next) {
		popover = helper_popover(lookup_widget(help_tips[tip_index].widget), _(help_tips[tip_index].tip));
		g_timeout_add(PREVIEW_DELAY, (GSourceFunc) helper_popover_close, (gpointer) popover);
		go_next = FALSE;
		tip_index++;
	}
	return tip_index < G_N_ELEMENTS(help_tips);
}

static void helper_notify_update(gpointer user_data) {
	if (timer_id != 0) {
		g_source_remove(timer_id);
	}
	timer_id = g_timeout_add_full(G_PRIORITY_DEFAULT_IDLE, 250,
			(GSourceFunc) update_popover, user_data,
			(GDestroyNotify) free_struct);
}

void start_helper_script() {
	tip_index = 0;
	go_next = TRUE;
	helper_notify_update(NULL);
}
