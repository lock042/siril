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
		{"headerbar", N_("Welcome to the newest version of Siril, "PACKAGE_STRING". Please take a moment to read some tips about this release"), 8, SIRIL_INTRO_GTK_POPOVER},
		{"header_processing_button", N_("The Processing menu has been reorganized by theme to streamline its structure and reduce its size. Many new processing tools have been developed and integrated in this menu."), 9, SIRIL_INTRO_GTK_POPOVER},
		{"header_tools_button", N_("A new menu, Tools, has been created to centralize Siril's tools, which were previously scattered throughout the interface. It now includes statistics, astrometry and photometry tools, as well as image analysis features."), 12, SIRIL_INTRO_GTK_POPOVER},
		{"header_scripts_button", N_("The Script menu has also evolved, now featuring both Python scripts and a script editor."), 7, SIRIL_INTRO_GTK_POPOVER},
		{"python_window", N_("The new script editor can be used to write both Siril and Python scripts. It offers syntax highlighting and many other useful features. Additionally, you can execute the script directly from the editor's interface."), 12, SIRIL_INTRO_WORKAROUND_POPOVER},
		{"icc_main_window_button", N_("This version of Siril includes a color management tool. A left-click on this button allows you to manage ICC profiles, while a right-click displays the image in soft proofing mode."), 9, SIRIL_INTRO_GTK_POPOVER},
		{"scripts_page", N_("Script management has been entirely revamped. You can now download new scripts directly from the preferences interface. These scripts are hosted by the Siril team, but can also be contributed by the community."), 12, SIRIL_INTRO_WORKAROUND_POPOVER},
		{"drawingarear", N_("Enjoy using the new Siril. You can restart this introduction at any moment in the Miscellaneous tab of the preferences"), 8, SIRIL_INTRO_GTK_POPOVER}
};

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

//static void open_menu_if_needed(GtkWidget *widget) {
//	if (!widget)
//		return;
//
//	GtkWidget *parent = widget;
//
//	while (parent != NULL) {
//		if (GTK_IS_MENU_BUTTON(parent)) {
//			GtkPopover *popover = gtk_menu_button_get_popover(GTK_MENU_BUTTON(parent));
//			if (popover && !gtk_widget_is_visible(GTK_WIDGET(popover))) {
//				gtk_widget_show(GTK_WIDGET(popover));
//			}
//		}
//
//		parent = gtk_widget_get_parent(parent);
//	}
//}

//static void hide_all_except(GtkWindow *keep_visible) {
//	GList *toplevels = gtk_window_list_toplevels();
//
//	for (GList *l = toplevels; l != NULL; l = l->next) {
//		GtkWindow *window = GTK_WINDOW(l->data);
//
//		if (window != keep_visible) {
//			gtk_widget_hide(GTK_WIDGET(window));
//		}
//	}
//
//	g_list_free(toplevels);
//}


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

	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(gtk_widget_get_toplevel(widget)));
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
	g_free(markup_txt);
	return window;
}


static gboolean intro_popover_close(gpointer user_data) {
	SirilUIIntro *ui = (SirilUIIntro *) user_data;

	gtk_widget_hide(ui->popover);
#ifndef OS_OSX // very slow on macOS
	gtk_style_context_remove_class(gtk_widget_get_style_context(ui->widget), "siril-intro-highlight");
#endif
	g_free(ui);
//	hide_all_except(GTK_WINDOW(lookup_widget("control_window")));
	go_next = TRUE;
	return FALSE;
}

static gboolean intro_popover_update(gpointer user_data) {
	if (go_next) {
		SirilUIIntro *ui = g_new(SirilUIIntro, 1);
		ui->widget = lookup_widget(intro_tips[tip_index].widget);
		ensure_widget_and_parents_visible(ui->widget);
//		open_menu_if_needed(ui->widget);
#ifndef OS_OSX // very slow on macOS
		gtk_style_context_add_class(gtk_widget_get_style_context(ui->widget), "siril-intro-highlight");
#endif
		if (intro_tips[tip_index].type == SIRIL_INTRO_GTK_POPOVER) {
			ui->popover = intro_popover(ui->widget, _(intro_tips[tip_index].tip));
		} else {
			ui->popover = floating_window_new(ui->widget, _(intro_tips[tip_index].tip));
		}
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
