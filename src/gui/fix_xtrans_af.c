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

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/proto.h"
#include "gui/gui_state.h"
#include "gui/fix_xtrans_af.h"

static GtkEntry *xtrans_af_x = NULL, *xtrans_af_y = NULL;
static GtkEntry *xtrans_af_w = NULL, *xtrans_af_h = NULL;
static GtkEntry *xtrans_sample_x = NULL, *xtrans_sample_y = NULL;
static GtkEntry *xtrans_sample_w = NULL, *xtrans_sample_h = NULL;

static void fix_xtrans_af_init_statics(void) {
	if (xtrans_af_x) return;
	xtrans_af_x = GTK_ENTRY(gtk_builder_get_object(gui.builder, "xtrans_af_x"));
	xtrans_af_y = GTK_ENTRY(gtk_builder_get_object(gui.builder, "xtrans_af_y"));
	xtrans_af_w = GTK_ENTRY(gtk_builder_get_object(gui.builder, "xtrans_af_w"));
	xtrans_af_h = GTK_ENTRY(gtk_builder_get_object(gui.builder, "xtrans_af_h"));
	xtrans_sample_x = GTK_ENTRY(gtk_builder_get_object(gui.builder, "xtrans_sample_x"));
	xtrans_sample_y = GTK_ENTRY(gtk_builder_get_object(gui.builder, "xtrans_sample_y"));
	xtrans_sample_w = GTK_ENTRY(gtk_builder_get_object(gui.builder, "xtrans_sample_w"));
	xtrans_sample_h = GTK_ENTRY(gtk_builder_get_object(gui.builder, "xtrans_sample_h"));
}

void init_xtrans_ui_pixels() {
	char pixel_value[256];
	fix_xtrans_af_init_statics();

	if (com.pref.prepro.xtrans_af.w != 0 && com.pref.prepro.xtrans_af.h != 0) {
		g_snprintf(pixel_value, 256, "%d", com.pref.prepro.xtrans_af.x);
		gtk_entry_set_text(xtrans_af_x, pixel_value);
		g_snprintf(pixel_value, 256, "%d", com.pref.prepro.xtrans_af.y);
		gtk_entry_set_text(xtrans_af_y, pixel_value);
		g_snprintf(pixel_value, 256, "%d", com.pref.prepro.xtrans_af.w);
		gtk_entry_set_text(xtrans_af_w, pixel_value);
		g_snprintf(pixel_value, 256, "%d", com.pref.prepro.xtrans_af.h);
		gtk_entry_set_text(xtrans_af_h, pixel_value);
	}

	if (com.pref.prepro.xtrans_sample.w != 0 && com.pref.prepro.xtrans_sample.h != 0) {
		g_snprintf(pixel_value, 256, "%d", com.pref.prepro.xtrans_sample.x);
		gtk_entry_set_text(xtrans_sample_x, pixel_value);
		g_snprintf(pixel_value, 256, "%d", com.pref.prepro.xtrans_sample.y);
		gtk_entry_set_text(xtrans_sample_y, pixel_value);
		g_snprintf(pixel_value, 256, "%d", com.pref.prepro.xtrans_sample.w);
		gtk_entry_set_text(xtrans_sample_w, pixel_value);
		g_snprintf(pixel_value, 256, "%d", com.pref.prepro.xtrans_sample.h);
		gtk_entry_set_text(xtrans_sample_h, pixel_value);
	}

}
