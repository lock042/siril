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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"

#include "filters/linear_match.h"

static gchar *get_reference_filename() {
	GtkFileChooser *linearmatch_ref = GTK_FILE_CHOOSER(lookup_widget("reference_filechooser_linearmatch"));
	return siril_file_chooser_get_filename(linearmatch_ref);
}

static gdouble get_high_rejection() {
	GtkSpinButton *button = GTK_SPIN_BUTTON(lookup_widget("spin_linearmatch_high"));

	return gtk_spin_button_get_value(button);
}

static gdouble get_low_rejection() {
	GtkSpinButton *button = GTK_SPIN_BUTTON(lookup_widget("spin_linearmatch_low"));

	return gtk_spin_button_get_value(button);
}

/*** callbacks **/

void on_linearmatch_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("linearmatch_dialog");
}

gboolean linearmatch_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("linearmatch_dialog");
	return TRUE;
}

static gboolean linearmatch_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	free_generic_img_args(args);
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_linearmatch_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (!single_image_is_loaded())
		return;
	gchar *filename = get_reference_filename();
	if (!filename)
		return;
	double low = get_low_rejection();
	double high = get_high_rejection();
	fits ref = { 0 };
	if (readfits(filename, &ref, NULL, gfit->type == DATA_FLOAT)) {
		g_free(filename);
		return;
	}
	g_free(filename);

	struct linear_match_data *data = new_linear_match_data(&ref, low, high);
	if (!data) {
		clearfits(&ref);
		return;
	}

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		free_linear_match_data(data);
		return;
	}
	args->fit = gfit;
	args->image_hook = linear_match_image_hook;
	args->log_hook = linear_match_log_hook;
	args->idle_function = linearmatch_idle;
	args->description = _("Linear Match");
	args->verbose = TRUE;
	args->user = data;
	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		set_cursor_waiting(FALSE);
	}
}
