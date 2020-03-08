/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#include "core/proto.h"
#include "algos/statistics.h"
#include "core/arithm.h"
#include "core/processing.h"
#include "core/undo.h"
#include "gui/callbacks.h"
#include "gui/siril_preview.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "ddp.h"

static gboolean ddp_show_preview;
static float ddp_stretch_value;
static float ddp_intensity_value;
static float ddp_constant_value;

static void ddp_startup() {
	/* we convert to float */
	imstats *stat = statistics(NULL, -1, &gfit, 0, NULL, STATS_BASIC, TRUE);
	ddp_constant_value = (float)stat->mean * 3.f;
	free_stats(stat);

	if (gfit.type == DATA_USHORT) {
		long ndata = gfit.naxes[0] * gfit.naxes[1] * gfit.naxes[2];
		fit_replace_buffer(&gfit, ushort_buffer_to_float(gfit.data, ndata), DATA_FLOAT);
		set_precision_switch();
	} else {
		ddp_constant_value *= USHRT_MAX_SINGLE;
	}

	copy_gfit_to_backup();
	ddp_stretch_value = 1.f;
	ddp_intensity_value = 1.f;
}

static void ddp_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				"Processing: DDP (stretch=%.2f, int=%.2f)", ddp_stretch_value,
				ddp_intensity_value);
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int ddp_update_preview() {
	copy_backup_to_gfit();

	set_cursor_waiting(TRUE);

	ddp(&gfit, ddp_constant_value, ddp_stretch_value, ddp_intensity_value);
	return 0;
}

int ddp(fits *a, int level, float coeff, float sigma) {
	fits fit = { 0 };
	if (level < 0 || level > USHRT_MAX) {
		siril_log_color_message(_("ddp level argument must be [0, 65535]\n"), "red");
		return 1;
	}
	float l = ushort_to_float_range(level);

	int ret = copyfits(a, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
	if (!ret) ret = unsharp(&fit, sigma, 0, FALSE);
	if (!ret) ret = soper(&fit, l, OPER_ADD, TRUE);
	if (!ret) ret = nozero(&fit, 1);
	if (!ret) ret = siril_fdiv(a, &fit, l, TRUE);
	if (!ret) ret = soper(a, coeff, OPER_MUL, TRUE);
	clearfits(&fit);
	invalidate_stats_from_fit(a);
	return ret;
}


void on_menuitem_ddp_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("ddp_dialog");
}

void apply_ddp_cancel() {
	ddp_close(TRUE);
	siril_close_dialog("ddp_dialog");
}

void on_cancel_ddp_clicked(GtkMenuItem *menuitem, gpointer user_data) {
	ddp_close(TRUE);
	siril_close_dialog("ddp_dialog");
}

void on_apply_ddp_clicked(GtkButton *button, gpointer user_data) {
	if (ddp_show_preview == FALSE) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = ddp_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}

	ddp_close(FALSE);
	siril_close_dialog("ddp_dialog");
}

void on_ddp_dialog_show(GtkWidget *widget, gpointer user_data) {
	ddp_startup();;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_stretch_ddp")), ddp_stretch_value);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinbutton_int_ddp")),	ddp_intensity_value);
	set_notify_block(FALSE);

	ddp_show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("ddp_preview")));

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = ddp_update_preview;
	param->show_preview = ddp_show_preview;
	notify_update((gpointer) param);
}

/** adjusters **/
void on_spinbutton_stretch_ddp_value_changed(GtkSpinButton *button, gpointer user_data) {
	ddp_stretch_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = ddp_update_preview;
	param->show_preview = ddp_show_preview;
	notify_update((gpointer) param);
}

void on_spinbutton_int_ddp_value_changed(GtkSpinButton *button, gpointer user_data) {
	ddp_intensity_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = ddp_update_preview;
	param->show_preview = ddp_show_preview;
	notify_update((gpointer) param);
}

void on_ddp_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (ddp_show_preview == TRUE) {
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = ddp_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
	ddp_show_preview = !ddp_show_preview;
}
