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

#include "io/sequence.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"

static void set_selection_ratio(double ratio) {
	gui.ratio = ratio;
	enforce_ratio_and_clamp();
	update_display_selection();
	gui_function(new_selection_zone, NULL);
	redraw(REDRAW_OVERLAY);
}

void on_menuitem_selection_free_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		gui.ratio = 0.0;
	}
}

void on_menuitem_selection_preserve_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio((double)gfit.rx / (double)gfit.ry);
	}
}

void on_menuitem_selection_16_9_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(16.0 / 9.0);
	}
}

void on_menuitem_selection_3_2_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(3.0 / 2.0);
	}
}

void on_menuitem_selection_4_3_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(4.0 / 3.0);
	}
}

void on_menuitem_selection_1_1_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(1.0 / 1.0);
	}
}

void on_menuitem_selection_3_4_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(3.0 / 4.0);
	}
}

void on_menuitem_selection_2_3_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(2.0 / 3.0);
	}
}

void on_menuitem_selection_9_16_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(9.0 / 16.0);
	}
}

void on_menuitem_selection_all_activate(GtkMenuItem *menuitem, gpointer user_data) {
	com.selection.x = 0;
	com.selection.y = 0;
	com.selection.w = gfit.rx;
	com.selection.h = gfit.ry;
	// "Select All" need to reset any enforced ratio that would not match the ratio of the image
	// 1. it's nice to NOT enforce a ratio when the user just want to select the whole image
	// 2. it's nice to keep the enforced ratio if it does match the image
	if (gui.ratio != ((double)gfit.rx / (double)gfit.ry)) {
		set_selection_ratio(0.0);
	} else {
		set_selection_ratio((double)gfit.rx / (double)gfit.ry); // triggers the new_selection() callbacks etc.
	}
}

void menuitem_selection_guides_0_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.pref.gui.selection_guides = 0;
	}
}

void menuitem_selection_guides_2_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.pref.gui.selection_guides = 2;
	}
}

void menuitem_selection_guides_3_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.pref.gui.selection_guides = 3;
	}
}

void menuitem_selection_guides_5_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.pref.gui.selection_guides = 5;
	}
}
