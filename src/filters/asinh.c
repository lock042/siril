/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "algos/statistics.h"
#include "core/arithm.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"

#include "asinh.h"

static gboolean asinh_rgb_space = FALSE;
static float asinh_stretch_value = 0.0f, asinh_black_value = 0.0f;

static int asinh_update_preview() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview"))))
		copy_backup_to_gfit();
	fits *fit = gui.roi.active ? &gui.roi.fit : &gfit;
	asinhlut(fit, asinh_stretch_value, asinh_black_value, asinh_rgb_space);
	notify_gfit_modified();
	return 0;
}

void asinh_change_between_roi_and_image() {
	gui.roi.operation_supports_roi = TRUE;
	// If we are showing the preview, update it after the ROI change.
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));
	notify_update((gpointer) param);
}

static void asinh_startup() {
	add_roi_callback(asinh_change_between_roi_and_image);
	roi_supported(TRUE);
	copy_gfit_to_backup();
}

static void asinh_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		if (asinh_stretch_value != 0.0f || asinh_black_value != 0.0f) {
			siril_preview_hide();
		} else {
			clear_backup();
		}
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				_("Asinh Transformation: (stretch=%6.1lf, bp=%7.5lf)"),
				asinh_stretch_value, asinh_black_value);
	}
	roi_supported(FALSE);
	remove_roi_callback(asinh_change_between_roi_and_image);
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int asinh_process_all() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview"))))
		copy_backup_to_gfit();
	asinhlut(&gfit, asinh_stretch_value, asinh_black_value, asinh_rgb_space);
	populate_roi();
	notify_gfit_modified();
	return 0;
}

int asinhlut_ushort(fits *fit, float beta, float offset, gboolean human_luminance) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	float m_CB = 1.f;
	float norm = get_normalized_value(fit);
	float invnorm = 1.0f / norm;
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
			blend_data data;
#pragma omp simd
			for (int ii = 0; ii < 3; ii++)
				data.do_channel[ii] = TRUE;
			float r = buf[RLAYER][i] * invnorm;
			float g = buf[GLAYER][i] * invnorm;
			float b = buf[BLAYER][i] * invnorm;
			float rout, gout, bout;
			float rprime = max(0, (r - offset) / (1.0f - offset));
			float gprime = max(0, (g - offset) / (1.0f - offset));
			float bprime = max(0, (b - offset) / (1.0f - offset));

			float x = factor_red * rprime + factor_green * gprime + factor_blue * bprime;
			data.tf[RLAYER] = (rprime == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * rprime) / (rprime * asinh_beta);
			data.tf[GLAYER] = (gprime == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * gprime) / (gprime * asinh_beta);
			data.tf[BLAYER] = (bprime == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * bprime) / (bprime * asinh_beta);
			float k = (x == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * x) / (x * asinh_beta);
			data.sf[RLAYER] = min(1.0f, max(0.0f, (rprime * k)));
			data.sf[GLAYER] = min(1.0f, max(0.0f, (gprime * k)));
			data.sf[BLAYER] = min(1.0f, max(0.0f, (bprime * k)));

			rgbblend(&data, &rout, &gout, &bout, m_CB);

			buf[RLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f, rout)));
			buf[GLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f, gout)));
			buf[BLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f, bout)));
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->rx * 16)
#endif
		for (i = 0; i < n; i++) {
			float x = buf[RLAYER][i] * invnorm;
			float xprime = max(0, (x - offset) / (1.0f - offset));
			float k = (xprime == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * xprime) / (xprime * asinh_beta);
			buf[RLAYER][i] = roundf_to_WORD(norm * min(1.0f, max(0.0f,(xprime) * k)));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int asinhlut_float(fits *fit, float beta, float offset, gboolean human_luminance) {
	float *buf[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	float m_CB = 1.f;
	float asinh_beta = asinhf(beta);
	float factor_red = human_luminance ? 0.2126f : 0.3333f;
	float factor_green = human_luminance ? 0.7152f : 0.3333f;
	float factor_blue = human_luminance ? 0.0722f : 0.3333f;

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (i = 0; i < n; i++) {
			blend_data data;
#pragma omp simd
			for (int ii = 0; ii < 3; ii++)
				data.do_channel[ii] = TRUE;
			float r = buf[RLAYER][i];
			float g = buf[GLAYER][i];
			float b = buf[BLAYER][i];
			float rprime = max(0.0f, (r - offset) / (1.0f - offset));
			float gprime = max(0.0f, (g - offset) / (1.0f - offset));
			float bprime = max(0.0f, (b - offset) / (1.0f - offset));

			float x = factor_red * rprime + factor_green * gprime + factor_blue * bprime;
			data.tf[RLAYER] = (rprime == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * rprime) / (rprime * asinh_beta);
			data.tf[GLAYER] = (gprime == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * gprime) / (gprime * asinh_beta);
			data.tf[BLAYER] = (bprime == 0.f) ? 0.f : (beta == 0.f) ? 1.f : asinhf(beta * bprime) / (bprime * asinh_beta);
			float k = (x == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * x) / (x * asinh_beta);
			data.sf[RLAYER] = min(1.0f, max(0.0f, (rprime * k)));
			data.sf[GLAYER] = min(1.0f, max(0.0f, (gprime * k)));
			data.sf[BLAYER] = min(1.0f, max(0.0f, (bprime * k)));

			rgbblend(&data, &buf[RLAYER][i], &buf[GLAYER][i], &buf[BLAYER][i], m_CB);
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (i = 0; i < n; i++) {
			float x, k;
			x = buf[RLAYER][i];
			float xprime = max(0.0f, (x - offset) / (1.0f - offset));
			k = (xprime == 0.0f) ? 0.0f : (beta == 0.0f) ? 1.0f : asinhf(beta * xprime) / (xprime * asinh_beta);
			buf[RLAYER][i] = min(1.0f, max(0.0f,(xprime * k)));
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

	if (gui.rendering_mode == LINEAR_DISPLAY)
		setup_stretch_sliders(); // In linear mode, set sliders to 0 / 65535

	icc_auto_assign_or_convert(&gfit, ICC_ASSIGN_ON_STRETCH);

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

	/* default parameters do not transform image, no need to update preview */
}

void on_asinh_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_asinh_cancel();
}

void on_asinh_ok_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview"))) || gui.roi.active) {
		asinh_process_all();
	}

	apply_asinh_changes();
	siril_close_dialog("asinh_dialog");
}

void on_asinh_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_asinh_changes();
}

void on_asinh_undo_clicked(GtkButton *button, gpointer user_data) {
	double prev_stretch = asinh_stretch_value;
	double prev_bp = asinh_black_value;
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
	// Update the preview if the parameters have changed
	if (prev_stretch != 0.0 || prev_bp != 0.0) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = asinh_update_preview;
		param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));
		notify_update((gpointer) param);
	}
}

/*** adjusters **/
void on_spin_asinh_value_changed(GtkSpinButton *button, gpointer user_data) {
	asinh_stretch_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));
	notify_update((gpointer) param);
}

void on_black_point_spin_asinh_value_changed(GtkSpinButton *button, gpointer user_data) {
	asinh_black_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));
	notify_update((gpointer) param);
}

void on_asinh_RGBspace_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	asinh_rgb_space = gtk_toggle_button_get_active(togglebutton);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")));
	notify_update((gpointer) param);
}

void on_asinh_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("asinh_preview")))) {
		copy_backup_to_gfit();
		redraw(REMAP_ALL);
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = asinh_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
}
