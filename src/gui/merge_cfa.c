/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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
#include "core/command.h"
#include "algos/demosaicing.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"

void on_merge_cfa_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("merge_cfa_dialog");
}

void on_merge_cfa_apply_clicked(GtkButton *button, gpointer user_data) {
	fits *cfa0 = NULL, *cfa1 = NULL, *cfa2 = NULL, *cfa3 = NULL;
	gchar *f_cfa0 = NULL, *f_cfa1 = NULL, *f_cfa2 = NULL, *f_cfa3 = NULL;
	GtkFileChooser *filechooser0 = GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa0"));
	GtkFileChooser *filechooser1 = GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa1"));
	GtkFileChooser *filechooser2 = GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa2"));
	GtkFileChooser *filechooser3 = GTK_FILE_CHOOSER(lookup_widget("filechooser_cfa3"));
	GtkComboBox *combo_pattern = GTK_COMBO_BOX(lookup_widget("merge_cfa_pattern"));
	gint p = gtk_combo_box_get_active(combo_pattern);
	sensor_pattern pattern = (sensor_pattern) p;
	f_cfa0 = g_strdup(gtk_file_chooser_get_filename(filechooser0));
	if (readfits(f_cfa0, cfa0, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser0);
		return;
	}
	f_cfa1 = g_strdup(gtk_file_chooser_get_filename(filechooser1));
	if (readfits(f_cfa1, cfa1, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser1);
		return;
	}
	f_cfa2 = g_strdup(gtk_file_chooser_get_filename(filechooser2));
	if (readfits(f_cfa2, cfa2, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser2);
		return;
	}
	f_cfa3 = g_strdup(gtk_file_chooser_get_filename(filechooser3));
	if (readfits(f_cfa3, cfa3, NULL, FALSE)) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		gtk_file_chooser_unselect_all(filechooser3);
		return;
	}
	merge_cfa(cfa0, cfa1, cfa2, cfa3, pattern);
}

