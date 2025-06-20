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

#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "io/image_format_fits.h"
#include "algos/Def_Wavelet.h"
#include "wavelets.h"

static float wavelet_value[6];
static gboolean wavelet_show_preview;

static void reset_scale_w() {
	static GtkRange *range_w[6] = { NULL, NULL, NULL, NULL, NULL, NULL };

	set_notify_block(TRUE);
	for (int i = 0; i < 6; i++) {
		char widget_name[10];
		sprintf(widget_name, "scale_w%d", i);
		range_w[i] = GTK_RANGE(lookup_widget(widget_name));
		gtk_range_set_value(range_w[i], 1.f);
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_plans_w")), 6);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("wavelet_preview")), TRUE);
	set_notify_block(FALSE);
}

static int update_wavelets() {
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;

	tmpdir = g_get_tmp_dir();

	set_cursor_waiting(TRUE);

	for (int i = 0; i < gfit.naxes[2]; i++) {
		dir[i] = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		if (gfit.type == DATA_USHORT) {
			wavelet_reconstruct_file(dir[i], wavelet_value, gfit.pdata[i]);
		} else {
			wavelet_reconstruct_file_float(dir[i], wavelet_value,
					gfit.fpdata[i]);
		}
		g_free(dir[i]);
	}
	redraw(REMAP_ALL);
	return 0;
}

static void wavelets_startup() {
	for (int i = 0; i < 6; i++) {
		wavelet_value[i] = 1.f;
	}
	copy_gfit_to_backup();
}

/* This function computes wavelets with the number of Nbr_Plan and
 * extracts plan "Plan" in fit parameters */

int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer) {
	int chan, start, end, retval = 0;
	wave_transf_des wavelet[3] = { 0 };

	g_assert(fit->naxes[2] <= 3);

	float *Imag = NULL;
	size_t n = fit->naxes[0] * fit->naxes[1];
	if (fit->type == DATA_USHORT) {
		Imag = f_vector_alloc(n);
		if (!Imag) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	if (reqlayer < 0 || reqlayer >= 3) {
		start = 0;
		end = fit->naxes[2];
	}
	else {
		start = reqlayer;
		end = start + 1;
	}

	for (chan = start; chan < end; chan++) {
		int Nl, Nc;

		if (fit->type == DATA_USHORT) {
			/* float wavelet of data [0, 65535] */
			if (wavelet_transform(Imag, fit->ry, fit->rx, &wavelet[chan],
						Type, Nbr_Plan, fit->pdata[chan])) {
				retval = 1;
				break;
			}
		}
		else if (fit->type == DATA_FLOAT) {
			/* float wavelet of data [0, 1] */
			Imag = fit->fpdata[chan];
			if (wavelet_transform_float(Imag, fit->ry, fit->rx, &wavelet[chan],
						Type, Nbr_Plan)) {
				retval = 1;
				break;
			}
		} else { // Unknown fit->type
			retval = 1;
			break;
		}
		Nl = wavelet[chan].Nbr_Ligne;
		Nc = wavelet[chan].Nbr_Col;
		pave_2d_extract_plan(wavelet[chan].Pave.Data, Imag, Nl, Nc, Plan);
		if (fit->type == DATA_USHORT)
			reget_rawdata(Imag, Nl, Nc, fit->pdata[chan]);
		wave_io_free(&wavelet[chan]);
	}

	/* Free */
	if (fit->type == DATA_USHORT)
		free(Imag);
	return retval;
}

static gboolean end_wavelets_filter(gpointer p) {
	struct wavelets_filter_data *args = (struct wavelets_filter_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);

	set_cursor_waiting(FALSE);
	free(args);
	return FALSE;
}

gpointer extract_plans(gpointer p) {
	int i;
	fits fit = { 0 };
	struct wavelets_filter_data *args = (struct wavelets_filter_data *) p;

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);

	for (i = 0; i < args->Nbr_Plan; i++) {
		gchar *filename, *msg;
		if (copyfits(args->fit, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
			siril_log_message(_("Could not copy image, aborting\n"));
			siril_add_idle(end_wavelets_filter, args);
			return GINT_TO_POINTER(1);
		}
		filename = g_strdup_printf("layer%02d", i);
		msg = g_strdup_printf(_("Extracting %s..."), filename);
		set_progress_bar_data(msg, (float)i / args->Nbr_Plan);
		get_wavelet_layers(&fit, args->Nbr_Plan, i, args->Type, -1);
		savefits(filename, &fit);
		g_free(filename);
		g_free(msg);
	}
	clearfits(&fit);
	siril_add_idle(end_wavelets_filter, args);
	return GINT_TO_POINTER(0);
}

void on_wavelets_dialog_show(GtkWidget *widget, gpointer user_data) {
	wavelet_show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("wavelet_preview")));
	wavelets_startup();
}

void on_wavelets_dialog_hide(GtkWidget *widget, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("frame_wavelets"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("button_reset_w"), FALSE);
	clear_backup();
}

void on_button_reset_w_clicked(GtkButton *button, gpointer user_data) {
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combobox_type_w")), 1);
	reset_scale_w();
	update_wavelets();
}

void apply_wavelets_cancel() {
	if (gtk_widget_get_sensitive(lookup_widget("frame_wavelets")) == TRUE) {
		reset_scale_w();
		update_wavelets();
	}
	siril_close_dialog("wavelets_dialog");
}

void on_button_ok_w_clicked(GtkButton *button, gpointer user_data) {
	gboolean is_active = gtk_widget_get_sensitive(lookup_widget("frame_wavelets")) == TRUE;
	if (is_active && wavelet_show_preview == FALSE) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = update_wavelets;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
	undo_save_state(get_preview_gfit_backup(), _("Wavelets Transformation"));

	siril_close_dialog("wavelets_dialog");
}

gboolean on_button_cancel_w_clicked(GtkButton *button, gpointer user_data) {
	apply_wavelets_cancel();
	return FALSE;
}

void on_button_compute_w_clicked(GtkButton *button, gpointer user_data) {
	if (wavelet_show_preview) {
		copy_backup_to_gfit();
	}

	int Type_Transform, Nbr_Plan, maxplan, mins, i;
	int nb_chan = gfit.naxes[2];
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;

	tmpdir = g_get_tmp_dir();

	Nbr_Plan = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_plans_w")));
	Type_Transform = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combobox_type_w"))) + 1;

	mins = min(gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if (Nbr_Plan > maxplan) {
		char *msg = siril_log_message(
				_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
		Nbr_Plan = maxplan;
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_plans_w")), Nbr_Plan);
	}

	if (Type_Transform != TO_PAVE_LINEAR && Type_Transform != TO_PAVE_BSPLINE) {
		char *msg = siril_log_message(_("Wavelet: type must be %d or %d\n"),
		TO_PAVE_LINEAR, TO_PAVE_BSPLINE);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
	}

	set_cursor_waiting(TRUE);

	if (gfit.type == DATA_USHORT) {
		size_t n = gfit.naxes[0] * gfit.naxes[1] * sizeof(float);
		float *Imag = malloc(n);
		if (Imag) {
			for (i = 0; i < nb_chan; i++) {
				dir[i] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[i]) + 2);
				strcpy(dir[i], tmpdir);
				strcat(dir[i], G_DIR_SEPARATOR_S);
				strcat(dir[i], File_Name_Transform[i]);
				wavelet_transform_file(Imag, gfit.ry, gfit.rx, dir[i],
						Type_Transform, Nbr_Plan, gfit.pdata[i]);
				free(dir[i]);
			}
			free(Imag);
		}
	} else {
		for (i = 0; i < nb_chan; i++) {
			dir[i] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[i]) + 2);
			strcpy(dir[i], tmpdir);
			strcat(dir[i], G_DIR_SEPARATOR_S);
			strcat(dir[i], File_Name_Transform[i]);
			wavelet_transform_file_float(gfit.fpdata[i], gfit.ry, gfit.rx, dir[i],
					Type_Transform, Nbr_Plan);
			free(dir[i]);
		}
	}
	gtk_widget_set_sensitive(lookup_widget("frame_wavelets"), TRUE);
	gtk_widget_set_sensitive(lookup_widget("button_reset_w"), TRUE);
	reset_scale_w();
	set_cursor_waiting(FALSE);
	return;
}

/****************** GUI for Wavelet Layers Extraction *****************/

void on_button_extract_w_ok_clicked(GtkButton *button, gpointer user_data) {
	int Nbr_Plan, Type, maxplan, mins;
	static GtkSpinButton *Spin_Nbr_Plan = NULL;
	static GtkComboBox *Combo_Wavelets_Type = NULL;

	if (Spin_Nbr_Plan == NULL) {
		Spin_Nbr_Plan = GTK_SPIN_BUTTON(lookup_widget("spinbutton_extract_w"));
		Combo_Wavelets_Type = GTK_COMBO_BOX(
				lookup_widget("combo_interpolation_extract_w"));
	}

	Nbr_Plan = gtk_spin_button_get_value(Spin_Nbr_Plan);
	Type = gtk_combo_box_get_active(Combo_Wavelets_Type) + 1;// 1: linear, 2: bspline

	set_cursor_waiting(TRUE);
	mins = min(gfit.rx, gfit.ry);
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
	args->fit = &gfit;
	if (!start_in_new_thread(extract_plans, args))
		free(args);
}

void on_button_extract_w_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("extract_wavelets_layers_dialog");
}

void on_spinbutton_plans_w_value_changed(GtkSpinButton *button, gpointer user_data) {
	int i;
	gint current_value = gtk_spin_button_get_value_as_int(button);
	for (i = 0; i < 6; i++) {
		gchar *tmp = g_strdup_printf("box_w%d", i);
		gtk_widget_set_visible(lookup_widget(tmp), current_value > i);
		g_free(tmp);
	}
}

void on_spin_w0_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[0] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = wavelet_show_preview;
	notify_update((gpointer) param);
}

void on_spin_w1_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[1] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = wavelet_show_preview;
	notify_update((gpointer) param);
}

void on_spin_w2_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[2] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = wavelet_show_preview;
	notify_update((gpointer) param);
}

void on_spin_w3_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[3] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = wavelet_show_preview;
	notify_update((gpointer) param);
}

void on_spin_w4_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[4] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = wavelet_show_preview;
	notify_update((gpointer) param);
}

void on_spin_w5_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[5] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	param->show_preview = wavelet_show_preview;
	notify_update((gpointer) param);
}

void on_wavelet_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	wavelet_show_preview = gtk_toggle_button_get_active(button);
	cancel_pending_update();
	if (!wavelet_show_preview) {
	/* if user click very fast */
		waiting_for_thread();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = update_wavelets;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}

