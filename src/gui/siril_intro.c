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
static gboolean go_prev;
static GSList *key_handler_windows = NULL;
static SirilUIIntro *current_ui = NULL;
static guint current_timeout_id = 0;
static gboolean processing_key = FALSE;

/* Structure to keep track of windows and their key handlers */
typedef struct {
	GtkWidget *window;
	gulong handler_id;
} WindowKeyHandler;



const SirilTipIntro intro_tips[] = {
		{"headerbar", N_("Welcome to the newest version of Siril, %s. Please take a moment to read some tips about this release. If the rate is too slow, skip to the next tip using the right arrow of your keyboard. You can also go back to the previous tip using the left arrow and stop the animation with the escape key."), 14, SIRIL_INTRO_GTK_POPOVER},
		{"header_processing_button", N_("The Processing menu has been reorganized by theme to streamline its structure and reduce its size. Many new processing tools have been developed and integrated in this menu."), 10, SIRIL_INTRO_GTK_POPOVER},
		{"header_tools_button", N_("A new menu, Tools, has been created to centralize Siril's tools, which were previously scattered throughout the interface. It now includes statistics, FITS header, astrometry and photometry tools, as well as image analysis features."), 12, SIRIL_INTRO_GTK_POPOVER},
		{"header_scripts_button", N_("The Script menu has also evolved, now featuring both Python scripts and a script editor."), 7, SIRIL_INTRO_GTK_POPOVER},
		{"python_window", N_("The new script editor can be used to write both Siril and Python scripts. It offers syntax highlighting and many other useful features. Additionally, you can execute the script directly from the editor's interface."), 11, SIRIL_INTRO_WORKAROUND_POPOVER},
		{"hamburger-menu", N_("The hamburger menu has been streamlined and now includes quick access to the documentation, preferences, keyboard shortcuts help, and the script tab in the preferences, making it easy to add new scripts."), 11, SIRIL_INTRO_GTK_POPOVER},
		{"scripts_page", N_("Script management has been entirely revamped. You can now download new scripts directly from the preferences interface. These scripts are hosted by the Siril team, but can also be contributed by the community."), 11, SIRIL_INTRO_WORKAROUND_POPOVER},
		{"label28", N_("The Registration, Plot and Stacking tabs have been redesigned for improved clarity and better organization. This release adds true Drizzle registration and mosaic stacking."), 9, SIRIL_INTRO_GTK_POPOVER},
		{"icc_main_window_button", N_("This version of Siril includes a color management tool. A left-click on this button allows you to manage ICC profiles, while a right-click displays the image in soft proofing mode."), 10, SIRIL_INTRO_GTK_POPOVER},
		{"annotate_button", N_("Right-clicking the annotations button now allows selection of the annotation catalogues to display, and support for IAU constellation names and outlines has been added."), 9, SIRIL_INTRO_GTK_POPOVER},
		{"cut_button", N_("A new intensity profiling tool has been added. This allows creating mono or RGB profiles of pixel values along a line in the image."), 8, SIRIL_INTRO_GTK_POPOVER},
		{"cm_page", N_("All color management options are now centralized in the preferences. This allows for professional-level work with any mono or RGB ICC profile."), 8, SIRIL_INTRO_WORKAROUND_POPOVER},
		{"user_page", N_("For more precise control, such as customizing annotation colors or configuring mouse interactions, visit the User Interface tab in the preferences. Don't hesitate to use your mouse scroll to view all the available options."), 11, SIRIL_INTRO_WORKAROUND_POPOVER},
		{"drawingarear", N_("Enjoy using the new Siril. You can restart this introduction at any moment in the Miscellaneous tab of the preferences."), 7, SIRIL_INTRO_GTK_POPOVER}
};

static void hide_all_except(GtkWindow *keep_visible) {
	GList *toplevels = gtk_window_list_toplevels();

	for (GList *l = toplevels; l != NULL; l = l->next) {
		GtkWindow *window = GTK_WINDOW(l->data);

		if (window != keep_visible) {
			gtk_widget_hide(GTK_WIDGET(window));
		}
	}

	g_list_free(toplevels);
}

static void on_window_destroy(GtkWidget *widget, gpointer user_data) {
	WindowKeyHandler *handler = (WindowKeyHandler *)user_data;

	if (handler) {
		handler->window = NULL;  // invalidate window
		handler->handler_id = 0;
	}
}

/* Global key handler for any window */
static gboolean on_key_press(GtkWidget *widget, GdkEventKey *event, gpointer user_data) {
	if (!processing_key && current_ui) {
		processing_key = TRUE;  // Set lock

		if (current_timeout_id > 0) {
			g_source_remove(current_timeout_id);
			current_timeout_id = 0;
		}

		if (event->keyval == GDK_KEY_Right && tip_index < G_N_ELEMENTS(intro_tips)) {
			gtk_widget_hide(current_ui->popover);
			gtk_style_context_remove_class(gtk_widget_get_style_context(current_ui->widget), "siril-intro-highlight");
			g_free(current_ui);
			current_ui = NULL;
			hide_all_except(GTK_WINDOW(lookup_widget("control_window")));
			go_next = TRUE;
			go_prev = FALSE;
		} else if (event->keyval == GDK_KEY_Left && tip_index > 1) {
			gtk_widget_hide(current_ui->popover);
			gtk_style_context_remove_class(gtk_widget_get_style_context(current_ui->widget), "siril-intro-highlight");
			g_free(current_ui);
			current_ui = NULL;
			hide_all_except(GTK_WINDOW(lookup_widget("control_window")));
			tip_index -= 2; // Go back two steps (one for the previous tip, one to counter the increment)
			go_next = TRUE;
			go_prev = TRUE;
		} else if (event->keyval == GDK_KEY_Escape) {
			gtk_widget_hide(current_ui->popover);
			gtk_style_context_remove_class(gtk_widget_get_style_context(current_ui->widget), "siril-intro-highlight");
			g_free(current_ui);
			current_ui = NULL;
			hide_all_except(GTK_WINDOW(lookup_widget("control_window")));
			tip_index = G_N_ELEMENTS(intro_tips); // Go at the end, and close everything
		}

		processing_key = FALSE;  // Release lock
	}
	return TRUE;
}

/* Function to disconnect key handlers from all windows */
static void disconnect_key_handlers() {
	GSList *current = key_handler_windows;
	GSList *next;

	while (current != NULL) {
		WindowKeyHandler *handler = current->data;
		next = current->next;

		if (handler) {
			if (handler->window && GTK_IS_WINDOW(handler->window)) {
				g_signal_handlers_disconnect_by_func(handler->window, G_CALLBACK(on_window_destroy), handler);

				if (handler->handler_id > 0) {
					g_signal_handler_disconnect(handler->window, handler->handler_id);
				}
			}
			g_free(handler);
		}

		key_handler_windows = g_slist_delete_link(key_handler_windows, current);
		current = next;
	}
}

/* Function to add key handler to a window */
static void add_key_handler_to_widget(GtkWidget *widget) {
	if (!widget || !GTK_IS_WIDGET(widget))
		return;

	WindowKeyHandler *handler = g_new0(WindowKeyHandler, 1);
	handler->window = widget;
	handler->handler_id = g_signal_connect(widget, "key-press-event", G_CALLBACK(on_key_press), NULL);

	g_signal_connect(widget, "destroy", G_CALLBACK(on_window_destroy), handler);

	key_handler_windows = g_slist_prepend(key_handler_windows, handler);
}


static void ensure_widget_and_parents_visible(GtkWidget *widget) {
	if (!widget)
		return;

	GtkWidget *parent = widget;
	GSList *parents_to_show = NULL;

	while (parent != NULL) {
		if (GTK_IS_NOTEBOOK(parent)) {
			GtkWidget *page = gtk_widget_get_ancestor(widget, GTK_TYPE_WIDGET);
			gint page_num = gtk_notebook_page_num(GTK_NOTEBOOK(parent), page);
			if (page_num >= 0) {
				gtk_notebook_set_current_page(GTK_NOTEBOOK(parent), page_num);
			}
		} else if (GTK_IS_STACK(parent)) {
			GtkWidget *target_child = widget;

			while (target_child && gtk_widget_get_parent(target_child) != parent) {
				target_child = gtk_widget_get_parent(target_child);
			}

			if (target_child) {
				gtk_stack_set_visible_child(GTK_STACK(parent), target_child);
			}
		}

		if (!gtk_widget_get_visible(parent)) {
			parents_to_show = g_slist_prepend(parents_to_show, parent);
		}

		parent = gtk_widget_get_parent(parent);
	}

	for (GSList *l = parents_to_show; l != NULL; l = l->next) {
		GtkWidget *widget_to_show = GTK_WIDGET(l->data);
		gtk_widget_show(widget_to_show);
	}

	g_slist_free(parents_to_show);
}

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

static GtkWidget* floating_window_new(GtkWidget *widget, const gchar *text) {
	GtkWidget *window, *box, *image, *label;
	gchar *markup_txt = g_strdup_printf("<big><b>%s</b></big>", text);
	GtkWindow *parent_window = GTK_WINDOW(gtk_widget_get_toplevel(widget));

	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_transient_for(GTK_WINDOW(window), parent_window);
	gtk_window_set_decorated(GTK_WINDOW(window), FALSE);
	gtk_window_set_type_hint(GTK_WINDOW(window), GDK_WINDOW_TYPE_HINT_DIALOG);

	box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_container_set_border_width(GTK_CONTAINER(box), 12);

	image = gtk_image_new_from_icon_name("dialog-information-symbolic", GTK_ICON_SIZE_DIALOG);

	label = gtk_label_new(NULL);
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_line_wrap_mode(GTK_LABEL(label), PANGO_WRAP_WORD);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 80);
	gtk_label_set_markup(GTK_LABEL(label), markup_txt);

	gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_NONE);
	gtk_widget_set_valign(label, GTK_ALIGN_CENTER);

	gtk_box_pack_start(GTK_BOX(box), image, FALSE, FALSE, 10);
	gtk_box_pack_start(GTK_BOX(box), label, TRUE, TRUE, 10);
	gtk_container_add(GTK_CONTAINER(window), box);

	gtk_window_set_resizable(GTK_WINDOW(window), FALSE);

	gtk_widget_show_all(window);

	gint parent_x, parent_y, parent_width, parent_height;
	gtk_window_get_position(parent_window, &parent_x, &parent_y);
	gtk_window_get_size(parent_window, &parent_width, &parent_height);

	gint window_width, window_height;
	gtk_window_get_size(GTK_WINDOW(window), &window_width, &window_height);

	GdkDisplay *display = gdk_display_get_default();
	GdkMonitor *monitor = gdk_display_get_primary_monitor(display);
	GdkRectangle workarea;
	gdk_monitor_get_workarea(monitor, &workarea);

	gint x = parent_x + (parent_width - window_width) / 2;
	gint y;

	x = CLAMP(x, workarea.x, workarea.x + workarea.width - window_width);

	if (parent_y - window_height - 10 >= workarea.y) {
		y = parent_y - window_height - 10;
	} else if (parent_y + parent_height + window_height + 10 <= workarea.y + workarea.height) {
		y = parent_y + parent_height + 10;
	} else {
		y = parent_y + (parent_height - window_height) / 2;
	}

	gtk_window_move(GTK_WINDOW(window), x, y);

	g_free(markup_txt);
	return window;
}

static gboolean intro_popover_close(gpointer user_data) {
	SirilUIIntro *ui = (SirilUIIntro*) user_data;

	if (ui == current_ui) {
		gtk_widget_hide(ui->popover);
		gtk_style_context_remove_class(gtk_widget_get_style_context(ui->widget), "siril-intro-highlight");
		g_free(ui);
		current_ui = NULL;
		hide_all_except(GTK_WINDOW(lookup_widget("control_window")));
		go_next = TRUE;
	}

	current_timeout_id = 0;
	return FALSE;
}

static gboolean intro_popover_update(gpointer user_data) {
	if (go_next) {
		current_ui = g_new(SirilUIIntro, 1);
		current_ui->widget = lookup_widget(intro_tips[tip_index].widget);
		ensure_widget_and_parents_visible(current_ui->widget);
		gtk_style_context_add_class(gtk_widget_get_style_context(current_ui->widget), "siril-intro-highlight");

		gchar *tip_text;
		if (tip_index == 0) {
			gchar *translated_format = _(intro_tips[tip_index].tip);
			tip_text = g_strdup_printf(translated_format, PACKAGE_STRING);
		} else {
			tip_text = g_strdup(_(intro_tips[tip_index].tip));
		}

		if (intro_tips[tip_index].type == SIRIL_INTRO_GTK_POPOVER) {
			current_ui->popover = intro_popover(current_ui->widget, tip_text);
		} else {
			current_ui->popover = floating_window_new(current_ui->widget, tip_text);
		}
		g_free(tip_text);

		add_key_handler_to_widget(current_ui->popover);

		current_timeout_id = g_timeout_add(INTRO_DELAY * intro_tips[tip_index].delay, (GSourceFunc) intro_popover_close, (gpointer) current_ui);

		go_next = FALSE;
		go_prev = FALSE;
		tip_index++;
	}

	if (tip_index >= G_N_ELEMENTS(intro_tips)) {
		disconnect_key_handlers();
		return FALSE;
	}

	return TRUE;
}

static void intro_notify_update(gpointer user_data) {
	g_timeout_add(250, (GSourceFunc) intro_popover_update, user_data);
}

void start_intro_script() {
	tip_index = 0;
	go_next = TRUE;
	go_prev = FALSE;
	current_ui = NULL;
	current_timeout_id = 0;
	key_handler_windows = NULL;

	/* Add key handler to main window and all toplevel windows */
	GList *toplevels = gtk_window_list_toplevels();
	for (GList *l = toplevels; l != NULL; l = l->next) {
		add_key_handler_to_widget(GTK_WIDGET(l->data));
	}
	g_list_free(toplevels);

	intro_notify_update(NULL);
}
