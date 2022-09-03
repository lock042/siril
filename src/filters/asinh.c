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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"

#include "asinh.h"

static gboolean asinh_rgb_space = FALSE;
static float asinh_stretch_value = 0.0f, asinh_black_value = 0.0f;
static gboolean asinh_show_preview;

static void asinh_startup() {
	copy_gfit_to_backup();
}

static void asinh_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				_("Asinh Transformation: (stretch=%6.1lf, bp=%7.5lf)"),
				asinh_stretch_value, asinh_black_value);
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int asinh_update_preview() {
	copy_backup_to_gfit();
	asinhlut(&gfit, asinh_stretch_value, asinh_black_value, asinh_rgb_space);
	return 0;
}

int asinhlut_ushort(fits *fit, float beta, float offset, gboolean human_luminance) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };

	float norm = get_normalized_value(fit);
	float invnorm = 1.0 / norm;
	float asinh_beta = asinh(beta);
	float factor_red = human_luminance ? 0.2126f : 0.3333f;
	float factor_green = human_luminance ? 0.7152f : 0.3333f;
	float factor_blue = human_luminance ? 0.0722f : 0.3333f;

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->rx * 16)
#endif
		for (i = 0; i < n; i++) {
			float r = buf[RLAYER][i] * invnorm;
			float g = buf[GLAYER][i] * invnorm;
			float b = buf[BLAYER][i] * invnorm;
			float rprime = max(0, (r - offset) / (1.0f - offset));
			float gprime = max(0, (g - offset) / (1.0f - offset));
			float bprime = max(0, (b - offset) / (1.0f - offset));

			float x = factor_red * rprime + factor_green * gprime + factor_blue * bprime;

			float k = (x == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinh(beta * x) / (x * asinh_beta);

			buf[RLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f,(rprime) * k)));
			buf[GLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f,(gprime) * k)));
			buf[BLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f,(bprime) * k)));
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->rx * 16)
#endif
		for (i = 0; i < n; i++) {
			float x = buf[RLAYER][i] * invnorm;
			float xprime = max(0, (x - offset) / (1.0f - offset));
			float k = (xprime == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinh(beta * xprime) / (xprime * asinh_beta);
			buf[RLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f,(xprime) * k)));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int asinhlut_float(fits *fit, float beta, float offset, gboolean human_luminance) {
	float *buf[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	float asinh_beta = asinh(beta);
	float factor_red = human_luminance ? 0.2126f : 0.3333f;
	float factor_green = human_luminance ? 0.7152f : 0.3333f;
	float factor_blue = human_luminance ? 0.0722f : 0.3333f;

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {
			float r = buf[RLAYER][i];
			float g = buf[GLAYER][i];
			float b = buf[BLAYER][i];
			float rprime = max(0.0f, (r - offset) / (1.0f - offset));
			float gprime = max(0.0f, (g - offset) / (1.0f - offset));
			float bprime = max(0.0f, (b - offset) / (1.0f - offset));

			float x = factor_red * rprime + factor_green * gprime + factor_blue * bprime;

			float k = (x == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinh(beta * x) / (x * asinh_beta);

			buf[RLAYER][i] = min(1.0f, max(0.0f, (rprime * k)));
			buf[GLAYER][i] = min(1.0f, max(0.0f, (gprime * k)));
			buf[BLAYER][i] = min(1.0f, max(0.0f, (bprime * k)));
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {
			float x, k;
			x = buf[RLAYER][i];
			float xprime = max(0.0f, (x - offset) / (1.0f - offset));
			k = (xprime == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinh(beta * xprime) / (xprime * asinh_beta);
			buf[RLAYER][i] = min(1.0f, max(0.0f,(x * k)));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int asinhlut(fits *fit, float beta, float offset, gboolean human_luminance) {
	if (fit->type == DATA_USHORT)
		return asinhlut_ushort(fit, beta, offset, human_luminance);
	if (fit->type == DATA_FLOAT)
		return asinhlut_float(fit, beta, offset, human_luminance);
	return 1;
}


static void apply_asinh_changes() {
	gboolean status = (asinh_stretch_value != 1.0f) || (asinh_black_value != 0.0f) || !asinh_rgb_space;
	asinh_close(!status);
}

void apply_asinh_cancel() {
	asinh_close(TRUE);
	siril_close_dialog("asinh_dialog");
}

/*** callbacks **/

void on_asinh_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_stretch = GTK_SPIN_BUTTON(lookup_widget("spin_asinh"));
	GtkSpinButton *spin_black_p = GTK_SPIN_BUTTON(lookup_widget("black_point_spin_asinh"));
	GtkToggleButton *toggle_rgb = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace"));

	asinh_startup();
	asinh_stretch_value = 0.0f;
	asinh_black_value = 0.0f;
	asinh_rgb_space = TRUE;

	set_notify_block(TRUE);
	gtk_toggle_button_set_active(toggle_rgb, asinh_rgb_space);
	gtk_spin_button_set_value(spin_stretch, asinh_stretch_value);
	gtk_spin_button_set_value(spin_black_p, asinh_black_value);
	gtk_spin_button_set_increments(spin_black_p, 0.001, 0.01);
	set_notify_block(FALSE);

	asinh_show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = asinh_show_preview;
	notify_update((gpointer) param);
}

void on_asinh_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_asinh_cancel();
}

void on_asinh_ok_clicked(GtkButton *button, gpointer user_data) {
	if (asinh_show_preview == FALSE) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = asinh_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}

	apply_asinh_changes();
	siril_close_dialog("asinh_dialog");
}

void on_asinh_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_asinh_changes();
}

void on_asinh_undo_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_stretch = GTK_SPIN_BUTTON(lookup_widget("spin_asinh"));
	GtkSpinButton *spin_black_p = GTK_SPIN_BUTTON(lookup_widget("black_point_spin_asinh"));
	GtkToggleButton *toggle_rgb = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace"));
	asinh_stretch_value = 0.0;
	asinh_black_value = 0.0;
	asinh_rgb_space = TRUE;

	set_notify_block(TRUE);
	gtk_toggle_button_set_active(toggle_rgb, asinh_rgb_space);
	gtk_spin_button_set_value(spin_stretch, asinh_stretch_value);
	gtk_spin_button_set_value(spin_black_p, asinh_black_value);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = asinh_show_preview;
	notify_update((gpointer) param);
}

/*** adjusters **/
void on_spin_asinh_value_changed(GtkSpinButton *button, gpointer user_data) {
	asinh_stretch_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = asinh_show_preview;
	notify_update((gpointer) param);
}

void on_black_point_spin_asinh_value_changed(GtkSpinButton *button, gpointer user_data) {
	asinh_black_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = asinh_show_preview;
	notify_update((gpointer) param);
}

void on_asinh_RGBspace_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	asinh_rgb_space = gtk_toggle_button_get_active(togglebutton);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = asinh_show_preview;
	notify_update((gpointer) param);
}

void on_asinh_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (asinh_show_preview == TRUE) {
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = asinh_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
	asinh_show_preview = !asinh_show_preview;
}
