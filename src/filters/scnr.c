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

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"

#include "scnr.h"

const char *scnr_type_to_string(scnr_type t) {
	switch (t) {
		default:
		case SCNR_AVERAGE_NEUTRAL:
			return _("average neutral");
		case SCNR_MAXIMUM_NEUTRAL:
			return _("maximum neutral");
		case SCNR_MAXIMUM_MASK:
			return _("maximum mask");
		case SCNR_ADDITIVE_MASK:
			return _("additive mask");
	}
}

/* Subtractive Chromatic Noise Reduction */
gpointer scnr(gpointer p) {
	struct scnr_data *args = (struct scnr_data *) p;
	g_assert(args->fit->type == DATA_USHORT || args->fit->type == DATA_FLOAT);
	size_t i, nbdata = args->fit->naxes[0] * args->fit->naxes[1];
	gint nb_above_1 = 0;
	struct timeval t_start, t_end;
	double norm = get_normalized_value(args->fit);
	double invnorm = 1.0 / norm;

	siril_log_color_message(_("SCNR: processing with %s algorithm%s...\n"),
			"green", scnr_type_to_string(args->type),
			args->preserve ? _(", preserving lightness") : "");
	gettimeofday(&t_start, NULL);

	int error = 0;

	cmsHPROFILE cielab_profile = NULL;
	cmsColorSpaceSignature sig;
	cmsUInt32Number src_type, dest_type;
//	cmsUInt32Number datasize, bytesperline, bytesperplane;
	cmsHTRANSFORM transform = NULL, invtransform = NULL;
	if (com.icc.available) {
		cielab_profile = cmsCreateLab4Profile(NULL);
		sig = cmsGetColorSpace(args->fit->icc_profile);
		src_type = get_planar_formatter_type(sig, args->fit->type, FALSE);
		dest_type = get_planar_formatter_type(cmsSigLabData, args->fit->type, FALSE);
		transform = cmsCreateTransform(args->fit->icc_profile, src_type, cielab_profile, dest_type, com.pref.icc.processing_intent, 0);
		invtransform = cmsCreateTransform(cielab_profile, dest_type, args->fit->icc_profile, src_type, com.pref.icc.processing_intent, 0);
		cmsCloseProfile(cielab_profile);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (i = 0; i < nbdata; i++) {
		float red, green, blue;
		switch (args->fit->type) {
			case DATA_USHORT:
				red = args->fit->pdata[RLAYER][i] * invnorm;
				green = args->fit->pdata[GLAYER][i] * invnorm;
				blue = args->fit->pdata[BLAYER][i] * invnorm;
				break;
			case DATA_FLOAT:
				red = args->fit->fpdata[RLAYER][i];
				green = args->fit->fpdata[GLAYER][i];
				blue = args->fit->fpdata[BLAYER][i];
				break;
			default: // Default needs to be included in the switch to avoid warning about omitting DATA_UNSUPPORTED
				break;
		}

		float x, y, z, L, a, b, m;
		if (args->preserve) {
			if (com.icc.available) {
				// This is horribly un-optimized! Only saving grace is it isn't worse than the original.
				// TODO: rewrite this code so that the La*b* transform is done in one hit (needs additional memory though)
				float in[3] = { red, green, blue };
				float out[3];
				cmsDoTransform(transform, (void*) &in, (void*) &out, 1);
				x = out[0];
				y = out[1];
				z = out[2];
			} else {
				rgb_to_xyzf(red, green, blue, &x, &y, &z);
				xyz_to_LABf(x, y, z, &L, &a, &b);
			}
		}

		switch (args->type) {
			case SCNR_AVERAGE_NEUTRAL:
				m = 0.5f * (red + blue);
				green = min(green, m);
				break;
			case SCNR_MAXIMUM_NEUTRAL:
				m = max(red, blue);
				green = min(green, m);
				break;
			case SCNR_MAXIMUM_MASK:
				m = max(red, blue);
				green = (green * (1.0f - args->amount) * (1.0f - m)) + (m * green);
				break;
			case SCNR_ADDITIVE_MASK:
				m = min(1.0f, red + blue);
				green = (green * (1.0f - args->amount) * (1.0f - m)) + (m * green);
		}

		if (args->preserve) {
			if (com.icc.available) {
				float in[3] = { red, green, blue };
				float out[3];
				cmsDoTransform(transform, (void*) &in, (void*) &out, 1);
				out[0] = L;
				cmsDoTransform(invtransform, (void*) &out, (void*) &in, 1);
				red = in[0];
				green = in[1];
				blue = in[2];
			} else {
				float tmp;
				rgb_to_xyzf(red, green, blue, &x, &y, &z);
				xyz_to_LABf(x, y, z, &tmp, &a, &b);
				LAB_to_xyzf(L, a, b, &x, &y, &z);
				xyz_to_rgbf(x, y, z, &red, &green, &blue);
			}
			if (red > 1.000001f || green > 1.000001f || blue > 1.000001f)
				g_atomic_int_inc(&nb_above_1);
		}

		if (args->fit->type == DATA_USHORT) {
			if (args->fit->orig_bitpix == BYTE_IMG) {
				args->fit->pdata[RLAYER][i] = round_to_BYTE(red * norm);
				args->fit->pdata[GLAYER][i] = round_to_BYTE(green * norm);
				args->fit->pdata[BLAYER][i] = round_to_BYTE(blue * norm);
			} else {
				args->fit->pdata[RLAYER][i] = round_to_WORD(red * norm);
				args->fit->pdata[GLAYER][i] = round_to_WORD(green * norm);
				args->fit->pdata[BLAYER][i] = round_to_WORD(blue * norm);
			}
		}
		else if (args->fit->type == DATA_FLOAT) {
			args->fit->fpdata[RLAYER][i] = set_float_in_interval(red, 0.0f, 1.0f);
			args->fit->fpdata[GLAYER][i] = set_float_in_interval(green, 0.0f, 1.0f);
			args->fit->fpdata[BLAYER][i] = set_float_in_interval(blue, 0.0f, 1.0f);
		}
	}

	cmsDeleteTransform(transform);
	cmsDeleteTransform(invtransform);
	/* normalize in case of preserve, it can under/overshoot */
	if (args->preserve && nb_above_1)
		siril_log_message("%d pixels were truncated to a maximum value of 1\n", nb_above_1);

	free(args);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	notify_gfit_modified();
	return GINT_TO_POINTER(error);
}

void on_SCNR_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkComboBox *comboscnr = GTK_COMBO_BOX(
			gtk_builder_get_object(gui.builder, "combo_scnr"));
	int type = gtk_combo_box_get_active(comboscnr);

	if (type == -1)
		gtk_combo_box_set_active(comboscnr, 0);
}

void on_SCNR_Apply_clicked(GtkButton *button, gpointer user_data) {
	int type = gtk_combo_box_get_active(
			GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combo_scnr")));
	GtkToggleButton *light_button = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(gui.builder, "preserve_light"));
	gboolean preserve = gtk_toggle_button_get_active(light_button);
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(gui.builder, "scale_scnr")));

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	undo_save_state(&gfit, _("SCNR (type=%s, amount=%0.2lf, preserve=%s)"),
			scnr_type_to_string(type), amount, preserve ? "true" : "false");

	args->fit = &gfit;
	args->type = type;
	args->amount = amount;
	args->preserve = preserve;
	set_cursor_waiting(TRUE);
	start_in_new_thread(scnr, args);
}

void on_SCNR_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("SCNR_dialog");
}

void on_combo_scnr_changed(GtkComboBoxText *box, gpointer user_data) {
	int type = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_scnr")));
	GtkScale *scale = GTK_SCALE(lookup_widget("scale_scnr"));
	GtkLabel *label = GTK_LABEL(lookup_widget("label56"));
	GtkSpinButton *spinButton = GTK_SPIN_BUTTON(lookup_widget("spin_scnr"));

	gtk_widget_set_sensitive(GTK_WIDGET(scale), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(label), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(spinButton), type > 1);
}
