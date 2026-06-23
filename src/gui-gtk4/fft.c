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

/* GTK callbacks and idle/thread functions for the Fourier Transform dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "filters/fft.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/file_browser.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/registration_preview.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"
#include "io/sequence.h"

static GtkNotebook *fft_notebook = NULL;
static GtkCheckButton *fft_centered = NULL;
static GtkEntry *fft_mag_entry = NULL, *fft_phase_entry = NULL;
static GtkWidget *fft_filechooser_mag = NULL, *fft_filechooser_phase = NULL;

/* Exposed so the activate handler can wire the file-chooser buttons
 * before the dialog is shown.  Without this the magnitude/phase buttons
 * have no click handler attached until the user hits Apply, by which
 * point they're useless. */
void fft_dialog_init_statics(void);
void fft_dialog_init_statics(void) {
	if (fft_notebook) return;
	fft_notebook = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook_fft"));
	fft_centered = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "fft_centered"));
	fft_mag_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "fftd_mag_entry"));
	fft_phase_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "fftd_phase_entry"));
	fft_filechooser_mag = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filechooser_mag"));
	fft_filechooser_phase = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filechooser_phase"));
	const gchar *fits_pattern =
		"*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS;*.fit.fz;*.FIT.fz;*.fits.fz;*.FITS.fz";
	siril_image_button_init(fft_filechooser_mag,
		_("Select magnitude image"), _("FITS files"), fits_pattern, NULL, NULL);
	siril_image_button_init(fft_filechooser_phase,
		_("Select phase image"), _("FITS files"), fits_pattern, NULL, NULL);
}

/* Idle function for the generic_image_worker path */
static gboolean fft_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (args->retval == 0)
		gfit_modified_update_gui();
	gui_function(redraw_previews, NULL);
	free_generic_img_args(args);
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_button_fft_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	gchar *mag, *phase;
	char *type = NULL, page;
	int type_order = -1;

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	fft_dialog_init_statics();
	page = gtk_notebook_get_current_page(fft_notebook);

	if (page == 0) {
		if (sequence_is_loaded()) {
			char *msg = siril_log_message(_("FFT does not work with sequences !\n"));
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
			set_cursor_waiting(FALSE);
			return;
		}
		if (!single_image_is_loaded()) {
			char *msg = siril_log_message(_("Open an image first !\n"));
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
			set_cursor_waiting(FALSE);
			return;
		}

		type_order = !siril_toggle_get_active(GTK_WIDGET(fft_centered));
		type = strdup("fftd");
		mag = g_strdup(gtk_editable_get_text(GTK_EDITABLE(fft_mag_entry)));
		phase = g_strdup(gtk_editable_get_text(GTK_EDITABLE(fft_phase_entry)));
	} else {
		type = strdup("ffti");
		mag = siril_file_chooser_get_filename(fft_filechooser_mag);
		phase = siril_file_chooser_get_filename(fft_filechooser_phase);

		if (mag == NULL || phase == NULL) {
			char *msg = siril_log_message(_("Select magnitude and phase before !\n"));
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
			set_cursor_waiting(FALSE);
			free(type);
			g_free(mag);
			g_free(phase);
			return;
		}
		open_single_image(mag);
	}

	if ((mag != NULL) && (phase != NULL)) {
		set_cursor_waiting(TRUE);
		struct fft_data *fft_args = calloc(1, sizeof(struct fft_data));
		fft_args->destroy_fn = free_fft_data;
		fft_args->fit = gfit;
		fft_args->type = type;
		fft_args->modulus = mag;
		fft_args->phase = phase;
		fft_args->type_order = type_order;

		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			PRINT_ALLOC_ERR;
			free_fft_data(fft_args);
			return;
		}
		args->fit = gfit;
		args->mem_ratio = 2.0f;
		args->image_hook = fft_image_hook;
		args->log_hook = fft_log_hook;
		args->idle_function = fft_idle;
		args->description = _("Fourier Transform");
		args->verbose = TRUE;
		args->user = fft_args;

		if (!start_in_new_thread(generic_image_worker, args)) {
			free_generic_img_args(args);
		}
	} else {
		free(type);
		g_free(mag);
		g_free(phase);
	}
}

void on_button_fft_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("dialog_FFT");
}

gboolean fft_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("dialog_FFT");
	return TRUE;
}
