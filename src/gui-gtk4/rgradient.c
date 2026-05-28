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

/* GTK callbacks for the Rotational Gradient dialog. */

#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "algos/PSF.h"
#include "filters/rgradient.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"

static GtkEntry *rgradient_xc = NULL, *rgradient_yc = NULL;
static GtkRange *rgradient_scale_radial = NULL, *rgradient_scale_rot = NULL;

/* Exposed so the activate handler can cache the widget pointers before
 * the dialog is shown — previously init only ran from Apply / the
 * "use selection" button handler. */
void rgradient_dialog_init_statics(void);
void rgradient_dialog_init_statics(void) {
	if (rgradient_xc) return;
	rgradient_xc = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_rgradient_xc"));
	rgradient_yc = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entry_rgradient_yc"));
	rgradient_scale_radial = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_radial_rgradient"));
	rgradient_scale_rot = GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_rot_rgradient"));
}

static double get_xc(void) {
	return g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(rgradient_xc)), NULL);
}

static double get_yc(void) {
	return g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(rgradient_yc)), NULL);
}

static double get_dR(void) {
	return gtk_range_get_value(rgradient_scale_radial);
}

static double get_da(void) {
	return gtk_range_get_value(rgradient_scale_rot);
}

static int rgradient_process_with_worker(void) {
	if (!single_image_is_loaded())
		return 1;

	if (gfit->orig_bitpix == BYTE_IMG) {
		siril_log_error(_("Error: this process cannot be applied to 8b images\n"));
		return 1;
	}

	double xc = get_xc();
	double yc = get_yc();
	double dR = get_dR();
	double da = get_da();

	if ((xc >= gfit->rx) || (yc >= gfit->ry)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong center coordinates"),
			_("The coordinates cannot be greater than the size of the image. Please change their values and retry."));
		return 1;
	}

	struct rgradient_data *params = new_rgradient_data();
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	params->xc = xc;
	params->yc = yc;
	params->dR = dR;
	params->da = da;
	params->fit = gfit;
	params->verbose = TRUE;

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_rgradient_data(params);
		return 1;
	}

	args->fit = gfit;
	args->mem_ratio = 3.0f;
	args->image_hook = rgradient_image_hook;
	args->log_hook = rgradient_log_hook;
	args->idle_function = NULL;
	args->description = _("Rotational Gradient");
	args->verbose = TRUE;
	args->user = params;
	args->mask_aware = TRUE;
	args->max_threads = com.max_thread;
	args->for_preview = FALSE;
	args->for_roi = FALSE;

	set_cursor_waiting(TRUE);

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return 1;
	}

	return 0;
}

void on_rgradient_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("rgradient_dialog");
}

gboolean rgradient_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("rgradient_dialog");
	return TRUE;
}

void on_rgradient_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	rgradient_dialog_init_statics();
	rgradient_process_with_worker();
}

void on_button_rgradient_selection_clicked(GtkButton *button, gpointer user_data) {
	rgradient_dialog_init_statics();
	if (com.selection.h && com.selection.w) {
		psf_error error = PSF_NO_ERR;
		psf_star *result = psf_get_minimisation(gfit, 0, &com.selection, FALSE, FALSE, NULL, TRUE, PSF_GAUSSIAN, &error);
		if (result && error == PSF_NO_ERR) {
			gchar *x0 = g_strdup_printf("%.3lf", result->x0 + com.selection.x);
			gtk_editable_set_text(GTK_EDITABLE(rgradient_xc), x0);
			gchar *y0 = g_strdup_printf("%.3lf", com.selection.y + com.selection.h - result->y0);
			gtk_editable_set_text(GTK_EDITABLE(rgradient_yc), y0);
			g_free(x0);
			g_free(y0);
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Center coordinate selection error"),
				_("No valid PSF found within selection."));
		}
		free_psf(result);
	}
}
