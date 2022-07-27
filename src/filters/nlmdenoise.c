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

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
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
#include "io/image_format_fits.h"
#include "filters/nlmdenoise.h"
#include "opencv/opencv.h"

/* Non-Local Means Denoising */
gpointer nlmdenoise(gpointer p) {
	NLMDenoise_data *args = (NLMDenoise_data *) p;
	fits *image = args->image;

	size_t ndata = image->naxes[0] * image->naxes[1];
	struct timeval t_start, t_end;
	gboolean orig_32bit = (image->type == DATA_FLOAT) ? TRUE : FALSE;
	int retval;

	siril_log_color_message(_("Non-local means denoise: processing...\n"), "green");
	gettimeofday(&t_start, NULL);

	/////////////////   MONO IMAGE   //////////////////////////
	if (image->naxes[2] == 1) {
		// FastNlMeansDenoising requires 8 or 16 bit input, if we have a 32bit fits we have to convert to 16bit
		if (orig_32bit) {
			siril_debug_print("Changing FIT data to USHORT\n");
			fit_replace_buffer(image, float_buffer_to_ushort(image->fdata, ndata), DATA_USHORT);
		}
		set_progress_bar_data("Running non-local means denoising...", 0.33);
		cvNLMDenoiseMono(image, image->naxes[0], image->naxes[1], args->h_lum);

		if (orig_32bit) {// If necessary, restore to 32bit to maintain precision for further calcs
			siril_debug_print("Changing FIT data to FLOAT\n");
			fit_replace_buffer(image, ushort_buffer_to_float(image->data, ndata), DATA_FLOAT);
		}

		retval = 0;
	}

	//////////////////   RGB IMAGE   //////////////////////////
	else {
		double x, y, z, *L, *A, *B, r, g, b, norm, invnorm;
		WORD *L_in, *A_in, *B_in;
		WORD *L_out, *A_out, *B_out;
		L = calloc(ndata, sizeof(double));
		A = calloc(ndata, sizeof(double));
		B = calloc(ndata, sizeof(double));
		L_in = calloc(ndata, sizeof(WORD));
		A_in = calloc(ndata, sizeof(WORD));
		B_in = calloc(ndata, sizeof(WORD));
		L_out = calloc(ndata, sizeof(WORD));
		A_out = calloc(ndata, sizeof(WORD));
		B_out = calloc(ndata, sizeof(WORD));
		norm = USHRT_MAX_DOUBLE;
		invnorm = 1.0 / norm;
		imstats *stat[image->naxes[2]];
		retval = compute_all_channels_statistics_single_image(image, STATS_BASIC | STATS_MAD, MULTI_THREADED, stat);
		gboolean seems_linear = FALSE;
		for (int i=0; i<image->naxes[2];i++) {
			if (stat[i]->median < 0.04)
				seems_linear = TRUE;
		}
		for (int i=0;i<image->naxes[2];i++)
			free_stats(stat[i]);
		if (seems_linear)
			siril_log_color_message(_("Warning: low mean pixel values in this image suggest it may still be linear. This non-local means algorithm is intended for use on stretched images. It is likely to produce bad results when applied to a linear image.\n"), "red");
		for (size_t i=0; i < ndata; i++) {
			if (args->image->type == DATA_USHORT) {
				r = image->pdata[RLAYER][i] * invnorm;
				g = image->pdata[GLAYER][i] * invnorm;
				b = image->pdata[BLAYER][i] * invnorm;
			} else if (args->image->type == DATA_FLOAT) {
				r = (double) image->fpdata[RLAYER][i];
				g = (double) image->fpdata[GLAYER][i];
				b = (double) image->fpdata[BLAYER][i];
			}
			rgb_to_xyz(r, g, b, &x, &y, &z);
			xyz_to_LAB(x, y, z, &L[i], &A[i], &B[i]);

			L_in[i] = (WORD) L[i] * norm / 100;
			A_in[i] = (WORD) (A[i] + 127) * norm / 256;
			B_in[i] = (WORD) (B[i] + 127) * norm / 256;
		}

		free(L);
		free(A);
		free(B);

		cvNLMDenoiseWORD(L_in, L_out, image->rx, image->ry, args->h_lum);
		set_progress_bar_data("Running non-local means denoising...", 0.33);
		cvNLMDenoiseWORD(A_in, A_out, image->rx, image->ry, args->h_AB);
		set_progress_bar_data("Running non-local means denoising...", 0.67);
		cvNLMDenoiseWORD(B_in, B_out, image->rx, image->ry, args->h_AB);
		set_progress_bar_data("Running non-local means denoising...", 1.00);

		free(L_in);
		free(A_in);
		free(B_in);

		for (size_t i=0; i < ndata; i++) {
			double L_final = L_out[i] * invnorm * 100;
			double A_final = (A_out[i] * invnorm * 256) - 127;
			double B_final = (B_out[i] * invnorm * 256) - 127;
			LAB_to_xyz(L_final, A_final, B_final, &x, &y, &z);
			xyz_to_rgb(x, y, z, &r, &g, &b);
			if(image->type == DATA_USHORT) {
				image->pdata[RLAYER][i] = (WORD) r * norm;
				image->pdata[GLAYER][i] = (WORD) g * norm;
				image->pdata[BLAYER][i] = (WORD) b * norm;
			} else if (image->type == DATA_FLOAT) {
				image->fpdata[0][i] = r;
				image->fpdata[1][i] = g;
				image->fpdata[2][i] = b;
			}
		}

		free(L_out);
		free(A_out);
		free(B_out);

		if (!orig_32bit)
			fit_replace_buffer(image, float_buffer_to_ushort(image->fdata, ndata), DATA_USHORT);

		retval = 0;
	}

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	notify_gfit_modified();
	set_progress_bar_data("Ready.", PROGRESS_RESET);
	return GINT_TO_POINTER(retval);
}
