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

#include <string.h>
#include <stdlib.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/gui_iface.h"
#include "io/image_format_fits.h"
#include "algos/Def_Wavelet.h"
#include "wavelets.h"

/************* wavelet transform worker (wavelet command and GUI compute path) *************/

gpointer wavelet_transform_worker(gpointer p) {
	struct wavelet_transform_data *args = (struct wavelet_transform_data *)p;
	const char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" };
	const char *tmpdir = g_get_tmp_dir();
	g_rw_lock_reader_lock(&gfit->rwlock);
	int nb_chan = gfit->naxes[2];

	if (gfit->type == DATA_USHORT) {
		float *Imag = f_vector_alloc(gfit->rx * gfit->ry);
		if (!Imag) {
			PRINT_ALLOC_ERR;
			free(args);
			g_rw_lock_reader_unlock(&gfit->rwlock);
			siril_add_idle(end_generic, NULL);
			return GINT_TO_POINTER(1);
		}
		for (int i = 0; i < nb_chan; i++) {
			gchar *dir = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
			wavelet_transform_file(Imag, gfit->ry, gfit->rx, dir,
					args->Type_Transform, args->Nbr_Plan, gfit->pdata[i], args->anscombe);
			g_free(dir);
		}
		free(Imag);
	} else if (gfit->type == DATA_FLOAT) {
		for (int i = 0; i < nb_chan; i++) {
			gchar *dir = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
			wavelet_transform_file_float(gfit->fpdata[i], gfit->ry, gfit->rx, dir,
					args->Type_Transform, args->Nbr_Plan, args->anscombe);
			g_free(dir);
		}
	} else {
		free(args);
		g_rw_lock_reader_unlock(&gfit->rwlock);
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}
	siril_log_message(_("Wavelet decomposition computed (%d plans)\n"), args->Nbr_Plan);
	free(args);
	g_rw_lock_reader_unlock(&gfit->rwlock);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(0);
}

/************* wrecons hook (shared with GUI OK path and process_wrecons) *************/

void free_wrecons_data(void *p) {
	free(p);
}

int wrecons_image_hook(struct generic_img_args *gargs, fits *fit, int threads) {
	struct wrecons_data *args = (struct wrecons_data *)gargs->user;
	const char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave", "b_rawdata.wave" };
	const char *tmpdir = g_get_tmp_dir();
	/* ROI preview: the transform covers the full image, so reconstruct into a
	 * full-size scratch buffer and copy the selection window into the smaller
	 * fit the worker passes. */
	const size_t full = (size_t) args->full_rx * (size_t) args->full_ry;
	for (int i = 0; i < args->nb_chan; i++) {
		gchar *dir = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		int ret;
		if (!args->for_roi) {
			if (fit->type == DATA_USHORT)
				ret = wavelet_reconstruct_file(dir, args->coef, &args->denoise, fit->pdata[i]);
			else if (fit->type == DATA_FLOAT)
				ret = wavelet_reconstruct_file_float(dir, args->coef, &args->denoise, fit->fpdata[i]);
			else {
				g_free(dir);
				return 1;
			}
		} else {
			const int w = fit->rx, h = fit->ry;
			if (fit->type == DATA_USHORT) {
				WORD *tmp = malloc(full * sizeof(WORD));
				if (!tmp) { g_free(dir); return 1; }
				ret = wavelet_reconstruct_file(dir, args->coef, &args->denoise, tmp);
				if (!ret)
					for (int r = 0; r < h; r++)
						memcpy(fit->pdata[i] + (size_t) r * w,
								tmp + (size_t) (args->full_ry - r - args->roi_y - 1) * args->full_rx + args->roi_x,
								(size_t) w * sizeof(WORD));
				free(tmp);
			} else if (fit->type == DATA_FLOAT) {
				float *tmp = malloc(full * sizeof(float));
				if (!tmp) { g_free(dir); return 1; }
				ret = wavelet_reconstruct_file_float(dir, args->coef, &args->denoise, tmp);
				if (!ret)
					for (int r = 0; r < h; r++)
						memcpy(fit->fpdata[i] + (size_t) r * w,
								tmp + (size_t) (args->full_ry - r - args->roi_y - 1) * args->full_rx + args->roi_x,
								(size_t) w * sizeof(float));
				free(tmp);
			} else {
				g_free(dir);
				return 1;
			}
		}
		g_free(dir);
		if (ret) return 1;
	}
	return 0;
}

gchar *wrecons_log_hook(gpointer p, log_hook_detail detail) {
	return g_strdup(_("Wavelet reconstruction"));
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
	stop_processing_thread();
	gui_iface.set_progress(PROGRESS_DONE, PROGRESS_TEXT_RESET);
	gui_iface.set_busy(FALSE);
	free(args);
	return FALSE;
}

gpointer extract_plans(gpointer p) {
	int i;
	fits fit = { 0 };
	struct wavelets_filter_data *args = (struct wavelets_filter_data *) p;

	gui_iface.set_progress(PROGRESS_RESET, NULL);

	for (i = 0; i < args->Nbr_Plan; i++) {
		gchar *filename, *msg;
		/* Extracted layers are written to disk for later recombination
		 * (e.g. via PixelMath); keep full metadata so a recombined result
		 * can retain the source's WCS/header. Geometry is unchanged, so the
		 * solution stays valid. */
		if (copyfits(args->fit, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT | CP_METADATA_HEAP, -1)) {
			siril_log_error(_("Could not copy image, aborting\n"));
			siril_add_idle(end_wavelets_filter, args);
			return GINT_TO_POINTER(1);
		}
		filename = g_strdup_printf("layer%02d", i);
		msg = g_strdup_printf(_("Extracting %s..."), filename);
		gui_iface.set_progress((double)i / args->Nbr_Plan, msg);
		get_wavelet_layers(&fit, args->Nbr_Plan, i, args->Type, -1);
		savefits(filename, &fit);
		g_free(filename);
		g_free(msg);
	}
	clearfits(&fit);
	siril_add_idle(end_wavelets_filter, args);
	return GINT_TO_POINTER(0);
}

/* GUI callbacks moved to src/gui/wavelets.c */
