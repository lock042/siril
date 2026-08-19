/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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

#include <math.h>

#include "core/siril.h"
#include "core/op_descriptors.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "filters/epf.h"
#include "gui-gtk4/epf.h"
#include "gui-gtk4/nde_editors.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/file_browser.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/siril_preview.h"
#include "core/undo.h"
#include "core/nde_replay.h"
#include "core/nde_history.h"
#include "opencv/opencv.h"

// Statics declarations
GtkButton *epf_undo = NULL, *epf_cancel = NULL, *epf_apply = NULL;
GtkDropDown *ep_filter_type = NULL;
GtkWindow *epf_dialog = NULL;
GtkWidget *guided_filter_guideimage = NULL;
GtkGrid *guide_image_widgets = NULL, *epf_sigma_spatial_settings = NULL, *epf_mod_settings = NULL;
GtkLabel *label176 = NULL, *label1 = NULL, *label177 = NULL;
GtkScale *scale_epf_d = NULL;
GtkSpinButton *spin_epf_d = NULL, *spin_epf_sigma_spatial = NULL, *spin_epf_sigma_col = NULL, *spin_epf_mod = NULL;
GtkCheckButton *guided_filter_selfguide = NULL, *epf_preview = NULL;
static GtkWidget *epf_amend_note = NULL;

/* Amend mode (convergence C5b): the dialog edits an existing history record.
 * gfit holds the pre-record state installed by nde_amend_preview_start;
 * previews run against it through the ordinary backup machinery, and Apply/
 * Cancel route through nde_amend_preview_end instead of a worker run + undo
 * save. */
static gboolean epf_amend_mode = FALSE;

// Static for loaded guide image
static fits loaded_fit = { 0 };
// Path of the loaded guide image, pinned into the record (phase 4.5 Convention 1)
static gchar *loaded_fit_path = NULL;

static void on_guided_filter_guideimage_picked(GtkWidget *button, const gchar *path, gpointer user_data);

// Statics init
void epf_dialog_init_statics() {
	if (epf_undo == NULL) {
		// GtkButton
		epf_undo = GTK_BUTTON(gtk_builder_get_object(gui.builder, "epf_undo"));
		epf_cancel = GTK_BUTTON(gtk_builder_get_object(gui.builder, "epf_cancel"));
		epf_apply = GTK_BUTTON(gtk_builder_get_object(gui.builder, "epf_apply"));
		// GtkDropDown
		ep_filter_type = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "ep_filter_type"));
		// GtkDialog
		epf_dialog = GTK_WINDOW(gtk_builder_get_object(gui.builder, "epf_dialog"));
		// Picker button for the guide image.
		guided_filter_guideimage = GTK_WIDGET(gtk_builder_get_object(gui.builder, "guided_filter_guideimage"));
		siril_image_button_init(guided_filter_guideimage,
			_("Select guide image"), _("FITS files"),
			"*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS;*.fit.fz;*.FIT.fz;*.fits.fz;*.FITS.fz",
			on_guided_filter_guideimage_picked, NULL);
		// GtkGrid
		guide_image_widgets = GTK_GRID(gtk_builder_get_object(gui.builder, "guide_image_widgets"));
		epf_sigma_spatial_settings = GTK_GRID(gtk_builder_get_object(gui.builder, "epf_sigma_spatial_settings"));
		epf_mod_settings = GTK_GRID(gtk_builder_get_object(gui.builder, "epf_mod_settings"));
		// GtkLabel
		label176 = GTK_LABEL(gtk_builder_get_object(gui.builder, "label176"));
		label1 = GTK_LABEL(gtk_builder_get_object(gui.builder, "label1"));
		label177 = GTK_LABEL(gtk_builder_get_object(gui.builder, "label177"));
		// GtkScale
		scale_epf_d = GTK_SCALE(gtk_builder_get_object(gui.builder, "scale_epf_d"));
		// GtkSpinButton
		spin_epf_d = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_epf_d"));
		spin_epf_sigma_spatial = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_epf_sigma_spatial"));
		spin_epf_sigma_col = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_epf_sigma_col"));
		spin_epf_mod = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_epf_mod"));
		// GtkToggleButton
		guided_filter_selfguide = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "guided_filter_selfguide"));
		epf_preview = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "epf_preview"));
		epf_amend_note = GTK_WIDGET(gtk_builder_get_object(gui.builder, "epf_amend_note"));
	}
}

/* Fill the widgets from the amended record's current parameters through the
 * op's own deserializer — the same struct the normal apply builds.  A file-
 * guide record is loaded back into loaded_fit so the guide pin survives the
 * round trip; a self-guide re-arms the self-guide toggle. */
static void epf_prefill_from_amend(void) {
	struct epfargs *p = op_desc_epf.deserialize(nde_amend_preview_params(),
	                                            nde_amend_preview_op_version());
	if (!p)
		return;
	gtk_spin_button_set_value(spin_epf_d, p->d);
	gtk_spin_button_set_value(spin_epf_sigma_col, p->sigma_col);
	gtk_spin_button_set_value(spin_epf_sigma_spatial, p->sigma_space);
	gtk_spin_button_set_value(spin_epf_mod, p->mod);
	gtk_drop_down_set_selected(ep_filter_type, p->filter);

	clearfits(&loaded_fit);
	g_clear_pointer(&loaded_fit_path, g_free);
	if (p->filter == EP_GUIDED) {
		if (p->guide_path && *p->guide_path) {
			/* file guide: reload the pinned guide image so Apply can re-pin
			 * it exactly (dimensions already matched at capture). */
			if (!readfits(p->guide_path, &loaded_fit, NULL, FALSE)) {
				loaded_fit_path = g_strdup(p->guide_path);
				gtk_widget_set_sensitive(GTK_WIDGET(guided_filter_selfguide), TRUE);
				siril_toggle_set_active(GTK_WIDGET(guided_filter_selfguide), FALSE);
			} else {
				/* guide file gone: fall back to self-guide so the form is
				 * still usable (the veto below should prevent reaching here). */
				clearfits(&loaded_fit);
				siril_toggle_set_active(GTK_WIDGET(guided_filter_selfguide), TRUE);
			}
		} else {
			siril_toggle_set_active(GTK_WIDGET(guided_filter_selfguide), TRUE);
		}
	}
	gtk_widget_set_visible(GTK_WIDGET(guide_image_widgets), p->filter != EP_BILATERAL);
	gtk_widget_set_visible(GTK_WIDGET(epf_sigma_spatial_settings), p->filter == EP_BILATERAL);

	if (p->destroy_fn)
		p->destroy_fn(p);
	else
		free(p);
}

/* Helper function to get current widget values */
static void get_epf_values(double *d, double *sigma_col, double *sigma_space, double *mod, ep_filter_t *filter_type) {
	if (d)
		*d = gtk_spin_button_get_value(spin_epf_d);
	if (sigma_col)
		*sigma_col = gtk_spin_button_get_value(spin_epf_sigma_col);
	if (sigma_space)
		*sigma_space = gtk_spin_button_get_value(spin_epf_sigma_spatial);
	if (mod)
		*mod = gtk_spin_button_get_value(spin_epf_mod);
	if (filter_type)
		*filter_type = (ep_filter_t)gtk_drop_down_get_selected(ep_filter_type);
}

/* Create and launch EPF processing */
static int epf_process_with_worker(gboolean for_preview) {
	// Allocate parameters
	struct epfargs *params = new_epf_args();
	if (!params) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	// Get current values from widgets
	get_epf_values(&params->d, &params->sigma_col, &params->sigma_space, &params->mod, &params->filter);
	params->fit = gfit;
	// Set up guide image
	params->guide_needs_freeing = FALSE;
	if (params->filter == EP_GUIDED) {
		if (siril_toggle_get_active(GTK_WIDGET(guided_filter_selfguide))) {
			params->guidefit = params->fit;
		} else {
			if (loaded_fit.rx != 0) {
				params->guidefit = &loaded_fit;
				/* pin the separate guide image for NDE replay (Convention 1) */
				params->guide_path = g_strdup(loaded_fit_path);
			} else {
				free_epf_args(params);
				return 1;
			}
		}
	} else {
		params->guidefit = NULL;
	}

	params->verbose = !for_preview;
	params->applying = !for_preview;

	// Allocate worker args
	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		free_epf_args(params);
		return 1;
	}

	// Set the fit based on whether ROI is active
	args->fit = gfit;
	args->op = &op_desc_epf;
	args->idle_function = NULL;
	args->verbose = !for_preview;
	args->user = params;
	args->max_threads = com.max_thread;
	args->for_preview = for_preview;

	if (!start_in_new_thread(generic_image_worker, args)) {
		free_generic_img_args(args);
		return 1;
	}
	return 0;
}

/* Update preview using the worker */
static int epf_update_preview() {
	if (siril_toggle_get_active(GTK_WIDGET(epf_preview))) {
		copy_backup_to_gfit();
		return epf_process_with_worker(TRUE);
	}
	return 0;
}

void epf_change_between_roi_and_image() {
	// If we are showing the preview, update it after the ROI change.
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(epf_preview));
	notify_update((gpointer) param);
}

static void epf_startup() {
	copy_gfit_to_backup();
	add_roi_callback(epf_change_between_roi_and_image);
	roi_declare_op(&op_desc_epf);
}

static void epf_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		copy_backup_to_gfit();
		gfit_modified_update_gui();
	} else {
		invalidate_stats_from_fit(gfit);
		double d, sigma_col, sigma_space;
		get_epf_values(&d, &sigma_col, &sigma_space, NULL, NULL);
	}

	roi_declare_op(NULL);
	remove_roi_callback(epf_change_between_roi_and_image);
	clearfits(&loaded_fit);
	clear_backup();
	set_cursor_waiting(FALSE);
}

static void apply_epf_changes() {
	double sigma_col, sigma_space;
	get_epf_values(NULL, &sigma_col, &sigma_space, NULL, NULL);
	gboolean status = (sigma_col != 0.0) || (sigma_space != 0.0);
	epf_close(!status);
}

void apply_epf_cancel() {
	if (epf_amend_mode) {
		/* Leave amend mode: drop the pre-record backup and restore the true
		 * pixels (nothing was committed). */
		cancel_pending_update();
		cancel_and_wait_for_preview();
		clearfits(&loaded_fit);
		g_clear_pointer(&loaded_fit_path, g_free);
		roi_declare_op(NULL);
		remove_roi_callback(epf_change_between_roi_and_image);
		clear_backup();
		nde_amend_preview_end(FALSE, NULL);
		epf_amend_mode = FALSE;
		siril_close_dialog("epf_dialog");
		return;
	}
	epf_close(TRUE);
	siril_close_dialog("epf_dialog");
}

/*** callbacks **/

void on_epf_dialog_show(GtkWidget *widget, gpointer user_data) {
	epf_dialog_init_statics();
	nde_amend_note_update(epf_amend_note, epf_amend_mode,
	                      &op_desc_epf);

	if (epf_amend_mode) {
		/* gfit already shows the pre-record state.  No ICC juggling — an
		 * amend only changes pixels.  The ROI stays armed: the backup
		 * epf_startup() takes IS that state, so a region preview crops the
		 * right pixels and the worker replays the record's successors inside
		 * the rectangle (nde_replay.h). */
		clearfits(&loaded_fit);
		g_clear_pointer(&loaded_fit_path, g_free);
		epf_startup();   /* ROI callback + descriptor + backup := pre-K */

		set_notify_block(TRUE);
		siril_toggle_set_active(GTK_WIDGET(epf_preview), TRUE);
		epf_prefill_from_amend();
		set_notify_block(FALSE);

		/* One preview tick so the dialog opens showing the record applied
		 * with its current parameters. */
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = epf_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
		return;
	}

	epf_startup();
	clearfits(&loaded_fit);

	set_notify_block(TRUE);
	gtk_widget_set_visible(GTK_WIDGET(guide_image_widgets), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(epf_sigma_spatial_settings), TRUE);
	/* GTK4: gtk_file_chooser_unselect_all removed */;
	siril_toggle_set_active(GTK_WIDGET(guided_filter_selfguide), TRUE);
	siril_toggle_set_active(GTK_WIDGET(epf_preview), FALSE);
	gtk_widget_set_sensitive(GTK_WIDGET(guided_filter_selfguide), FALSE);
	gtk_drop_down_set_selected(ep_filter_type, EP_BILATERAL);
	gtk_spin_button_set_value(spin_epf_d, 0.0);
	gtk_spin_button_set_value(spin_epf_sigma_spatial, 11.0);
	gtk_spin_button_set_value(spin_epf_sigma_col, 11.0);
	gtk_spin_button_set_value(spin_epf_mod, 1.0);
	set_notify_block(FALSE);
}

void on_epf_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_epf_cancel();
}

void on_epf_apply_clicked(GtkButton *button, gpointer user_data) {
	if (epf_amend_mode) {
		/* Serialize the widget state through the SAME epfargs struct and
		 * serializer the normal apply builds, then route it through the
		 * amend path.  guidefit + fit must be set so epf_serialize emits the
		 * guide token (self / file + pinned path) just as a capture does. */
		cancel_pending_update();
		cancel_and_wait_for_preview();
		struct epfargs *applied = new_epf_args();
		if (!applied) {
			PRINT_ALLOC_ERR;
			return;
		}
		get_epf_values(&applied->d, &applied->sigma_col, &applied->sigma_space,
				&applied->mod, &applied->filter);
		applied->fit = gfit;
		if (applied->filter == EP_GUIDED) {
			if (siril_toggle_get_active(GTK_WIDGET(guided_filter_selfguide))) {
				applied->guidefit = applied->fit;
			} else if (loaded_fit.rx != 0) {
				applied->guidefit = &loaded_fit;
				applied->guide_path = g_strdup(loaded_fit_path);
			} else {
				/* guided with no guide image: cannot serialize faithfully */
				free_epf_args(applied);
				siril_message_dialog(GTK_MESSAGE_ERROR,
						_("No guide image"),
						_("Select a guide image or use self-guide before applying."));
				return;
			}
		}
		gchar *blob = op_desc_epf.serialize(applied);
		/* free_epf_args would clearfits our loaded_fit guide; detach it */
		if (applied->guidefit == &loaded_fit)
			applied->guidefit = NULL;
		applied->guide_needs_freeing = FALSE;
		free_epf_args(applied);
		clearfits(&loaded_fit);
		g_clear_pointer(&loaded_fit_path, g_free);
		roi_declare_op(NULL);
		remove_roi_callback(epf_change_between_roi_and_image);
		clear_backup();
		nde_amend_preview_end(TRUE, blob);
		g_free(blob);
		epf_amend_mode = FALSE;
		siril_close_dialog("epf_dialog");
		return;
	}

	if (!check_ok_if_cfa())
		return;

	// If preview is on, need to copy backup to gfit first
	if (siril_toggle_get_active(GTK_WIDGET(epf_preview))) {
		copy_backup_to_gfit();
	}

	// Always process full image when Apply is clicked
	epf_process_with_worker(FALSE);

	apply_epf_changes();
	siril_close_dialog("epf_dialog");
}

void on_epf_dialog_close(GtkWindow *dialog, gpointer user_data) {
	apply_epf_changes();
}

void on_epf_undo_clicked(GtkButton *button, gpointer user_data) {
	if (epf_amend_mode) {
		/* Reset means "back to the record's recorded parameters" here. */
		set_notify_block(TRUE);
		epf_prefill_from_amend();
		siril_toggle_set_active(GTK_WIDGET(epf_preview), TRUE);
		set_notify_block(FALSE);
		copy_backup_to_gfit();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = epf_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
		return;
	}
	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_epf_d, 0);
	gtk_spin_button_set_value(spin_epf_sigma_col, 11);
	gtk_spin_button_set_value(spin_epf_sigma_spatial, 11);
	gtk_spin_button_set_value(spin_epf_mod, 1);
	gtk_drop_down_set_selected(ep_filter_type, EP_BILATERAL);
	siril_toggle_set_active(GTK_WIDGET(epf_preview), TRUE);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(epf_preview));
	notify_update((gpointer) param);
}

/*** adjusters **/
void on_ep_filter_type_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec;
	ep_filter_t filter_type = gtk_drop_down_get_selected(combo);
	gtk_widget_set_visible(GTK_WIDGET(guide_image_widgets), filter_type != EP_BILATERAL);
	gtk_widget_set_visible(GTK_WIDGET(epf_sigma_spatial_settings), filter_type == EP_BILATERAL);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(epf_preview));
	notify_update((gpointer) param);
}

void on_guided_filter_selfguide_toggled(GtkCheckButton *button, gpointer user_data) {
	gboolean active = siril_toggle_get_active(GTK_WIDGET(button));
	if (active) {
		clearfits(&loaded_fit);
		g_clear_pointer(&loaded_fit_path, g_free);
		/* GTK4: gtk_file_chooser_unselect_all removed */;
		gtk_widget_set_sensitive(GTK_WIDGET(guided_filter_selfguide), FALSE);
	}
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(epf_preview));
	notify_update((gpointer) param);
}

static void on_guided_filter_guideimage_picked(GtkWidget *button, const gchar *path, gpointer user_data) {
	(void)button; (void)user_data;
	clearfits(&loaded_fit);
	g_clear_pointer(&loaded_fit_path, g_free);
	if (!path) return;
	if (readfits(path, &loaded_fit, NULL, FALSE)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error: image could not be loaded"),
			_("Image loading failed"));
		clearfits(&loaded_fit);
		return;
	}
	if (loaded_fit.rx != gfit->rx || loaded_fit.ry != gfit->ry) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error: image dimensions do not match"),
			_("Image loading failed"));
		clearfits(&loaded_fit);
		return;
	}
	loaded_fit_path = g_strdup(path);
	gtk_widget_set_sensitive(GTK_WIDGET(guided_filter_selfguide), TRUE);
	siril_toggle_set_active(GTK_WIDGET(guided_filter_selfguide), FALSE);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(epf_preview));
	notify_update((gpointer) param);
}

void on_epf_parameter_changed(GtkWidget *button, gpointer user_data) {
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = epf_update_preview;
	param->show_preview = siril_toggle_get_active(GTK_WIDGET(epf_preview));
	notify_update((gpointer) param);
}

void on_epf_preview_toggled(GtkCheckButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!siril_toggle_get_active(GTK_WIDGET(epf_preview))) {
		/* if user click very fast */
		cancel_and_wait_for_preview();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = epf_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}

/* ---- amend-mode entry (nde_editors registry) --------------------------- */

static void epf_amend_ready(gboolean ok, gpointer user) {
	(void)user;
	if (!ok)
		return;   /* the core logged the reason; nothing was changed */
	epf_amend_mode = TRUE;
	siril_open_dialog("epf_dialog");
}

gboolean epf_open_amend(gint64 record_id) {
	/* A guided record pinning a SEPARATE guide file can only be represented
	 * here if that file is still loadable at matching dimensions.  If it is
	 * gone or resized, veto so the Edit button falls back to the kv-grid
	 * editor rather than silently collapsing the record to a self-guide. */
	GPtrArray *snap = nde_history_snapshot(NULL);
	gboolean ok = TRUE;
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *rec = g_ptr_array_index(snap, i);
		if (rec->record_id != record_id)
			continue;
		GHashTable *kv = rec->params ? nde_kv_parse(rec->params) : NULL;
		if (kv) {
			const char *guide = nde_kv_get_str(kv, "guide");
			if (guide && !strcmp(guide, "file")) {
				const char *path = nde_kv_get_str(kv, "operand_path");
				fits probe = { 0 };
				if (!path || !*path || readfits(path, &probe, NULL, FALSE)) {
					ok = FALSE;
				} else {
					if (probe.rx != gfit->rx || probe.ry != gfit->ry)
						ok = FALSE;
					clearfits(&probe);
				}
			}
			g_hash_table_unref(kv);
		}
		break;
	}
	if (snap)
		g_ptr_array_unref(snap);
	if (!ok)
		return FALSE;

	/* The dialog opens from the ready callback, once the pre-record state
	 * has been synthesized and installed. */
	nde_amend_preview_start(record_id, epf_amend_ready, NULL);
	return TRUE;
}
