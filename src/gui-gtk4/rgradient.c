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
#include "core/op_descriptors.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/nde_replay.h"
#include "algos/PSF.h"
#include "filters/rgradient.h"
#include "gui-gtk4/rgradient.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/utils.h"
#include "io/single_image.h"

static GtkEntry *rgradient_xc = NULL, *rgradient_yc = NULL;
static GtkRange *rgradient_scale_radial = NULL, *rgradient_scale_rot = NULL;
static GtkWidget *rgradient_amend_note = NULL;

/* Amend mode (convergence C5b): the dialog edits a filters.rgradient history
 * record.  gfit holds the pre-record state installed by
 * nde_amend_preview_start.  rgradient has no preview machinery (dialogs.c
 * has_preview=FALSE), so there is no backup to arm or clear: the display shows
 * the installed pre-record state while the form is open, and Apply/Cancel route
 * through nde_amend_preview_end. */
static gboolean rgradient_amend_mode = FALSE;

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
	rgradient_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "rgradient_amend_note"));
}

/* Fill the widgets from the amended record's current parameters through the
 * op's own deserializer.  Centre coordinates are plain text entries. */
static void rgradient_prefill_from_amend(void) {
	struct rgradient_data *p = op_desc_rgradient.deserialize(nde_amend_preview_params(),
	                                                        nde_amend_preview_op_version());
	if (!p)
		return;
	gchar buf[G_ASCII_DTOSTR_BUF_SIZE];
	gtk_editable_set_text(GTK_EDITABLE(rgradient_xc),
			g_ascii_dtostr(buf, sizeof buf, p->xc));
	gtk_editable_set_text(GTK_EDITABLE(rgradient_yc),
			g_ascii_dtostr(buf, sizeof buf, p->yc));
	gtk_range_set_value(rgradient_scale_radial, p->dR);
	gtk_range_set_value(rgradient_scale_rot, p->da);
	if (p->destroy_fn)
		p->destroy_fn(p);
	else
		free(p);
}

void on_rgradient_dialog_show(GtkWidget *widget, gpointer user_data) {
	rgradient_dialog_init_statics();
	gtk_widget_set_visible(rgradient_amend_note, rgradient_amend_mode);
	if (rgradient_amend_mode) {
		/* gfit already shows the pre-record state.  No ROI. */
		if (gui.roi.active)
			on_clear_roi();
		set_notify_block(TRUE);
		rgradient_prefill_from_amend();
		set_notify_block(FALSE);
	}
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
	args->op = &op_desc_rgradient;
	args->idle_function = NULL;
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
	if (rgradient_amend_mode) {
		/* Cancel: leave amend mode without committing (restore true pixels). */
		nde_amend_preview_end(FALSE, NULL);
		rgradient_amend_mode = FALSE;
	}
	siril_close_dialog("rgradient_dialog");
}

gboolean rgradient_hide_on_delete(GtkWidget *widget) {
	if (rgradient_amend_mode) {
		nde_amend_preview_end(FALSE, NULL);
		rgradient_amend_mode = FALSE;
	}
	siril_close_dialog("rgradient_dialog");
	return TRUE;
}

void on_rgradient_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (rgradient_amend_mode) {
		/* Serialize the widget state through the SAME struct and serializer
		 * the normal apply uses, then route it through the amend path.  The
		 * centre bounds check still applies. */
		rgradient_dialog_init_statics();
		struct rgradient_data applied = { 0 };
		applied.xc = get_xc();
		applied.yc = get_yc();
		applied.dR = get_dR();
		applied.da = get_da();
		if ((applied.xc >= gfit->rx) || (applied.yc >= gfit->ry)) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong center coordinates"),
				_("The coordinates cannot be greater than the size of the image. Please change their values and retry."));
			return;
		}
		gchar *blob = op_desc_rgradient.serialize(&applied);
		nde_amend_preview_end(TRUE, blob);
		g_free(blob);
		rgradient_amend_mode = FALSE;
		siril_close_dialog("rgradient_dialog");
		return;
	}

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

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void rgradient_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	rgradient_amend_mode = TRUE;
	siril_open_dialog("rgradient_dialog");
}

gboolean rgradient_open_amend(gint64 record_id) {
	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed. */
	nde_amend_preview_start(record_id, rgradient_amend_ready, NULL);
	return TRUE;
}
