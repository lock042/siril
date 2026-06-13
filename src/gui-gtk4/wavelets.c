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

/* GTK callbacks, preview logic and idle functions for the Wavelets dialogs. */

#include <math.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/undo.h"
#include "algos/Def_Wavelet.h"
#include "filters/wavelets.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/siril_preview.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/wavelets.h"
#include "io/single_image.h"

static float wavelet_value[6];
/* Guards the per-layer preview toggle handler against re-entrancy while we
 * set the checkbutton states programmatically (e.g. from reset_scale_w). */
static gboolean block_preview_toggle = FALSE;

static GtkRange *wavelets_scale_w[6] = { NULL };
static GtkSpinButton *wavelets_plans_spin = NULL;
static GtkCheckButton *wavelets_preview_btn[6] = { NULL };
static GtkWidget *wavelets_frame = NULL, *wavelets_reset_btn = NULL;
static GtkDropDown *wavelets_type_combo = NULL;
static GtkSpinButton *wavelets_extract_spin = NULL;
static GtkDropDown *wavelets_extract_interp = NULL;
static GtkWidget *wavelets_box_w[6] = { NULL };

/* A layer takes part in the live preview when its box is visible (i.e. within
 * the current number of layers) and its per-layer preview toggle is active. */
static gboolean layer_preview_active(int i) {
	return gtk_widget_get_visible(wavelets_box_w[i])
			&& siril_toggle_get_active(GTK_WIDGET(wavelets_preview_btn[i]));
}

static gboolean any_layer_preview_active(void) {
	for (int i = 0; i < 6; i++)
		if (layer_preview_active(i))
			return TRUE;
	return FALSE;
}

static void wavelets_dialog_init_statics(void) {
	if (wavelets_plans_spin) return;
	char name[20];
	for (int i = 0; i < 6; i++) {
		snprintf(name, sizeof(name), "scale_w%d", i);
		wavelets_scale_w[i] = GTK_RANGE(gtk_builder_get_object(gui.builder, name));
		snprintf(name, sizeof(name), "box_w%d", i);
		wavelets_box_w[i] = GTK_WIDGET(gtk_builder_get_object(gui.builder, name));
	}
	for (int i = 0; i < 6; i++) {
		snprintf(name, sizeof(name), "preview_w%d", i);
		wavelets_preview_btn[i] = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, name));
	}
	wavelets_plans_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbutton_plans_w"));
	wavelets_frame = GTK_WIDGET(gtk_builder_get_object(gui.builder, "frame_wavelets"));
	wavelets_reset_btn = GTK_WIDGET(gtk_builder_get_object(gui.builder, "button_reset_w"));
	wavelets_type_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combobox_type_w"));
	wavelets_extract_spin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spinbutton_extract_w"));
	wavelets_extract_interp = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_interpolation_extract_w"));
}

static void reset_scale_w() {
	set_notify_block(TRUE);
	block_preview_toggle = TRUE;
	for (int i = 0; i < 6; i++)
		gtk_range_set_value(wavelets_scale_w[i], 1.f);
	gtk_spin_button_set_value(wavelets_plans_spin, 6);
	for (int i = 0; i < 6; i++)
		siril_toggle_set_active(GTK_WIDGET(wavelets_preview_btn[i]), TRUE);
	block_preview_toggle = FALSE;
	set_notify_block(FALSE);
}

/* Reconstruct gfit from the wavelet transform files using the given per-layer
 * coefficients. */
static int reconstruct_gfit(float *coef) {
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;

	tmpdir = g_get_tmp_dir();

	set_cursor_waiting(TRUE);

	for (int i = 0; i < gfit->naxes[2]; i++) {
		dir[i] = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		if (gfit->type == DATA_USHORT) {
			wavelet_reconstruct_file(dir[i], coef, gfit->pdata[i]);
		} else {
			wavelet_reconstruct_file_float(dir[i], coef, gfit->fpdata[i]);
		}
		g_free(dir[i]);
	}
	notify_gfit_data_modified();
	redraw(REDRAW_ALL);
	return 0;
}

/* Preview reconstruction: only the layers whose preview toggle is active use
 * their adjusted strength; every other layer is held at the neutral value of
 * 1.0 so the preview isolates the effect of the previewed layers. */
static int update_wavelets() {
	float coef[6];
	for (int i = 0; i < 6; i++)
		coef[i] = layer_preview_active(i) ? wavelet_value[i] : 1.f;
	return reconstruct_gfit(coef);
}

static void wavelets_startup() {
	for (int i = 0; i < 6; i++) {
		wavelet_value[i] = 1.f;
	}
	copy_gfit_to_backup();
}

/* Idle function: called after generic_image_worker finishes the wrecons apply */
static gboolean wrecons_idle(gpointer p) {
	struct generic_img_args *args = (struct generic_img_args *)p;
	stop_processing_thread();
	if (!args->retval)
		gfit_modified_update_gui();
	free_generic_img_args(args);
	set_cursor_waiting(FALSE);
	siril_close_dialog("wavelets_dialog");
	return FALSE;
}

static gboolean wavelet_compute_idle(gpointer p) {
	stop_processing_thread();
	gtk_widget_set_sensitive(wavelets_frame, TRUE);
	gtk_widget_set_sensitive(wavelets_reset_btn, TRUE);
	reset_scale_w();
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_wavelets_dialog_show(GtkWidget *widget, gpointer user_data) {
	wavelets_dialog_init_statics();
	wavelets_startup();
}

void on_wavelets_dialog_hide(GtkWidget *widget, gpointer user_data) {
	gtk_widget_set_sensitive(wavelets_frame, FALSE);
	gtk_widget_set_sensitive(wavelets_reset_btn, FALSE);
	clear_backup();
}

void on_button_reset_w_clicked(GtkButton *button, gpointer user_data) {
	gtk_drop_down_set_selected(wavelets_type_combo, 1);
	reset_scale_w();
	update_wavelets();
}

void apply_wavelets_cancel(void) {
	if (gtk_widget_get_sensitive(wavelets_frame) == TRUE) {
		reset_scale_w();
		update_wavelets();
	}
	siril_close_dialog("wavelets_dialog");
}

void on_button_ok_w_clicked(GtkButton *button, gpointer user_data) {
	if (gtk_widget_get_sensitive(wavelets_frame) != TRUE) {
		/* Wavelets were never computed: nothing to apply. */
		siril_close_dialog("wavelets_dialog");
		return;
	}
	/* The live preview may have been showing an isolated reconstruction (only
	 * some layers at their adjusted strength). Whatever was previewed, the
	 * applied result must use every layer's real strength. If a preview is
	 * active, gfit holds reconstructed pixels and the backup holds the
	 * original: restore the original first so the worker snapshots the correct
	 * undo "before" state before applying the full reconstruction. */
	if (is_preview_active())
		copy_backup_to_gfit();

	set_cursor_waiting(TRUE);
	struct wrecons_data *wrecons_args = calloc(1, sizeof(struct wrecons_data));
	wrecons_args->destroy_fn = free_wrecons_data;
	wrecons_args->nb_chan = gfit->naxes[2];
	for (int i = 0; i < 6; i++)
		wrecons_args->coef[i] = wavelet_value[i];

	struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
	args->fit = gfit;
	args->image_hook = wrecons_image_hook;
	args->log_hook = wrecons_log_hook;
	args->idle_function = wrecons_idle;
	args->description = _("Wavelets Transformation");
	args->verbose = TRUE;
	args->user = wrecons_args;
	if (!start_in_new_thread(generic_image_worker, args))
		free_generic_img_args(args);
}

gboolean on_button_cancel_w_clicked(GtkButton *button, gpointer user_data) {
	apply_wavelets_cancel();
	return FALSE;
}

void on_button_compute_w_clicked(GtkButton *button, gpointer user_data) {
	if (is_preview_active()) {
		copy_backup_to_gfit();
	}

	int Nbr_Plan = gtk_spin_button_get_value(wavelets_plans_spin);
	int Type_Transform = gtk_drop_down_get_selected(wavelets_type_combo) + 1;

	int mins = min(gfit->rx, gfit->ry);
	int maxplan = log(mins) / log(2) - 2;

	if (Nbr_Plan > maxplan) {
		char *msg = siril_log_message(
				_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
		Nbr_Plan = maxplan;
		gtk_spin_button_set_value(wavelets_plans_spin, Nbr_Plan);
	}

	if (Type_Transform != TO_PAVE_LINEAR && Type_Transform != TO_PAVE_BSPLINE) {
		char *msg = siril_log_message(_("Wavelet: type must be %d or %d\n"),
		TO_PAVE_LINEAR, TO_PAVE_BSPLINE);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
	}

	struct wavelet_transform_data *args = calloc(1, sizeof(struct wavelet_transform_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}
	args->Nbr_Plan = Nbr_Plan;
	args->Type_Transform = Type_Transform;

	set_cursor_waiting(TRUE);
	if (!start_in_new_thread(wavelet_transform_worker, args)) {
		free(args);
		set_cursor_waiting(FALSE);
		return;
	}
	/* In GUI mode start_in_new_thread is fire-and-forget; the idle below
	 * re-enables the widgets once the worker completes.  In script mode it
	 * blocks, but on_button_compute_w_clicked is never called from scripts. */
	siril_add_idle(wavelet_compute_idle, NULL);
}

/****************** GUI for Wavelet Layers Extraction *****************/

void on_button_extract_w_ok_clicked(GtkButton *button, gpointer user_data) {
	int Nbr_Plan, Type, maxplan, mins;

	Nbr_Plan = gtk_spin_button_get_value(wavelets_extract_spin);
	Type = gtk_drop_down_get_selected(wavelets_extract_interp) + 1;// 1: linear, 2: bspline

	set_cursor_waiting(TRUE);
	mins = min(gfit->rx, gfit->ry);
	maxplan = log(mins) / log(2) - 2;

	if (Nbr_Plan > maxplan) {
		char *msg = siril_log_message(_("Wavelet: maximum number "
				"of plans for this image size is %d\n"), maxplan);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
		set_cursor_waiting(FALSE);
		return;
	}

	struct wavelets_filter_data *args = calloc(1, sizeof(struct wavelets_filter_data));

	args->Type = Type;
	args->Nbr_Plan = Nbr_Plan;
	args->fit = gfit;
	if (!start_in_new_thread(extract_plans, args))
		free(args);
}

void on_button_extract_w_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("extract_wavelets_layers_dialog");
}

void on_spinbutton_plans_w_value_changed(GtkSpinButton *button, gpointer user_data) {
	gint current_value = gtk_spin_button_get_value_as_int(button);
	for (int i = 0; i < 6; i++)
		gtk_widget_set_visible(wavelets_box_w[i], current_value > i);
}

void on_spin_w0_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[0] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = any_layer_preview_active();
	notify_update((gpointer) param);
}

void on_spin_w1_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[1] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = any_layer_preview_active();
	notify_update((gpointer) param);
}

void on_spin_w2_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[2] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = any_layer_preview_active();
	notify_update((gpointer) param);
}

void on_spin_w3_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[3] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = any_layer_preview_active();
	notify_update((gpointer) param);
}

void on_spin_w4_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[4] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = any_layer_preview_active();
	notify_update((gpointer) param);
}

void on_spin_w5_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[5] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = any_layer_preview_active();
	notify_update((gpointer) param);
}

/* Shared handler for every per-layer preview toggle. The preview shows the
 * combined effect of all currently-previewed layers, so any toggle simply
 * recomputes (or tears down) the single shared preview. */
void on_wavelet_preview_toggled(GtkCheckButton *button, gpointer user_data) {
	if (block_preview_toggle)
		return;
	cancel_pending_update();
	if (!any_layer_preview_active()) {
		/* Last previewed layer was just switched off: restore the original. */
		cancel_and_wait_for_preview();
		siril_preview_hide();
	} else {
		/* When no preview was active yet, gfit still holds the original
		 * pixels: snapshot them so the preview can be undone. Skip this when a
		 * preview is already showing, otherwise the reconstructed pixels would
		 * overwrite the original backup. */
		if (!is_preview_active())
			copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = update_wavelets;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
