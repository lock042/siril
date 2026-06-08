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

/*
 * GUI-side script console: status-bar logging, command history, and the
 * end-of-script cleanup idle.  Moved here from core/command_line_processor.c
 * so that the GTK code lives exclusively in a GUI-build file.
 */

#ifndef SRC_GUI_SCRIPT_CONSOLE_H_
#define SRC_GUI_SCRIPT_CONSOLE_H_

#include <gtk/gtk.h>

/* Schedule a status-bar update from any thread (passes ownership of msg). */
void console_log_status(const gchar *msg, int line);

/* Show the command-help popup for the entry widget (GSourceFunc signature). */
gboolean show_command_help_popup(gpointer user_data);

/* Initialise command-line completion and keyboard handler. */
void init_command(void);
/* GTK callback: show help popup when ? button is clicked. */
void on_GtkCommandHelper_clicked(GtkButton *button, gpointer user_data);

/* Clear the script status bar (call from GTK main thread or via gui_function). */
void console_clear_status_bar(void);

/* GTK idle: runs all post-script GUI cleanup.  Registered via gui_function
 * or gui_iface.end_script_gui().  Returns FALSE (one-shot). */
gboolean end_script(gpointer p);

#endif /* SRC_GUI_SCRIPT_CONSOLE_H_ */
