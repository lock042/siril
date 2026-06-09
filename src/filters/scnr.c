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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "algos/colors.h"

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

/*****************************************************************************
 *      S C N R      A L L O C A T O R   A N D   D E S T R U C T O R        *
 ****************************************************************************/

struct scnr_data *new_scnr_data(void) {
	struct scnr_data *args = calloc(1, sizeof(struct scnr_data));
	if (args) {
		args->destroy_fn = free_scnr_data;
	}
	return args;
}

void free_scnr_data(void *ptr) {
	struct scnr_data *args = (struct scnr_data *)ptr;
	if (!args)
		return;
	free(ptr);
}

gchar *scnr_log_hook(gpointer p, log_hook_detail detail) {
	struct scnr_data* args = (struct scnr_data*) p;
	return g_strdup_printf(_("SCNR: %s algorithm%s..."),
				scnr_type_to_string(args->type),
				args->preserve ? _(", preserving lightness") : "");
}

/* Subtractive Chromatic Noise Reduction - core processing function */
static int scnr_process(struct scnr_data *args, fits *fit) {
	g_assert(fit->type == DATA_USHORT || fit->type == DATA_FLOAT);
	size_t i, nbdata = fit->naxes[0] * fit->naxes[1];
	gint nb_above_1 = 0;

	gchar *msg = scnr_log_hook(args, SUMMARY);
	gui_iface.set_progress(PROGRESS_PULSATE, msg);
	g_free(msg);

	double norm = get_normalized_value(fit);
	double invnorm = 1.0 / norm;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (i = 0; i < nbdata; i++) {
		double red, green, blue;
		switch (fit->type) {
			case DATA_USHORT:
				red = fit->pdata[RLAYER][i] * invnorm;
				green = fit->pdata[GLAYER][i] * invnorm;
				blue = fit->pdata[BLAYER][i] * invnorm;
				break;
			case DATA_FLOAT:
				red = (double)fit->fpdata[RLAYER][i];
				green = (double)fit->fpdata[GLAYER][i];
				blue = (double)fit->fpdata[BLAYER][i];
				break;
			default:
				break;
		}

		double x, y, z, L, a, b, m;
		if (args->preserve) {
			linrgb_to_xyz(red, green, blue, &x, &y, &z, TRUE);
			xyz_to_LAB(x, y, z, &L, &a, &b);
		}

		switch (args->type) {
			case SCNR_AVERAGE_NEUTRAL:
				m = 0.5 * (red + blue);
				green = min(green, m);
				break;
			case SCNR_MAXIMUM_NEUTRAL:
				m = max(red, blue);
				green = min(green, m);
				break;
			case SCNR_MAXIMUM_MASK:
				m = max(red, blue);
				green = (green * (1.0 - args->amount) * (1.0 - m)) + (m * green);
				break;
			case SCNR_ADDITIVE_MASK:
				m = min(1.0, red + blue);
				green = (green * (1.0 - args->amount) * (1.0 - m)) + (m * green);
		}

		if (args->preserve) {
			double tmp;
			linrgb_to_xyz(red, green, blue, &x, &y, &z, TRUE);
			xyz_to_LAB(x, y, z, &tmp, &a, &b);
			LAB_to_xyz(L, a, b, &x, &y, &z);
			xyz_to_linrgb(x, y, z, &red, &green, &blue, TRUE);
			if (red > 1.000001 || green > 1.000001 || blue > 1.000001)
				g_atomic_int_inc(&nb_above_1);
		}

		if (fit->type == DATA_USHORT) {
			if (fit->orig_bitpix == BYTE_IMG) {
				fit->pdata[RLAYER][i] = round_to_BYTE(red * norm);
				fit->pdata[GLAYER][i] = round_to_BYTE(green * norm);
				fit->pdata[BLAYER][i] = round_to_BYTE(blue * norm);
			} else {
				fit->pdata[RLAYER][i] = round_to_WORD(red * norm);
				fit->pdata[GLAYER][i] = round_to_WORD(green * norm);
				fit->pdata[BLAYER][i] = round_to_WORD(blue * norm);
			}
		}
		else if (fit->type == DATA_FLOAT) {
			fit->fpdata[RLAYER][i] = set_float_in_interval(red, 0.0f, 1.0f);
			fit->fpdata[GLAYER][i] = set_float_in_interval(green, 0.0f, 1.0f);
			fit->fpdata[BLAYER][i] = set_float_in_interval(blue, 0.0f, 1.0f);
		}
	}

	if (args->preserve && nb_above_1)
		siril_log_message("%d pixels were truncated to a maximum value of 1\n", nb_above_1);

	return 0;
}

/* The actual SCNR processing hook for generic_image_worker */
int scnr_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct scnr_data *params = (struct scnr_data *)args->user;
	if (!params)
		return 1;
	return scnr_process(params, fit);
}

/* GUI callbacks moved to src/gui/scnr.c */
