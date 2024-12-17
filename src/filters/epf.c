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
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "core/arithm.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "gui/callbacks.h"

#include "filters/epf.h"

gboolean end_epf(gpointer p) {
	set_cursor_waiting(FALSE);
	stop_processing_thread();
	notify_gfit_modified();
	populate_roi();
	return FALSE;
}

gpointer epfhandler (gpointer args) {
	lock_roi_mutex();
	struct epfargs *p = (struct epfargs*) args;
	set_cursor_waiting(TRUE);
	int retval = edge_preserving_filter(p);
	unlock_roi_mutex();
	if (!com.script)
		siril_add_idle(end_epf, NULL);
	return GINT_TO_POINTER(retval);
}

gpointer epf_filter (gpointer args) {
	struct epfargs *p = (struct epfargs*) args;
	set_cursor_waiting(TRUE);
	int retval = edge_preserving_filter(p);
	if (!com.script)
		siril_add_idle(end_epf, NULL);
	return GINT_TO_POINTER(retval);
}

int match_guide_to_roi(fits *guide, fits *guide_roi) {
	int retval = 0;
	if (!gui.roi.active)
		return 1;
	uint32_t nchans = guide->naxes[2];
	size_t npixels = guide->rx * guide->ry;
	gboolean rgb = (nchans == 3);
	size_t npixels_roi = gui.roi.fit.rx * gui.roi.fit.ry;
	copyfits(guide, guide_roi, CP_FORMAT, -1);
	guide_roi->rx = guide_roi->naxes[0] = gui.roi.selection.w;
	guide_roi->ry = guide_roi->naxes[1] = gui.roi.selection.h;
	guide_roi->naxes[2] = nchans;
	guide_roi->naxis = nchans == 1 ? 2 : 3;
	if (guide->type == DATA_FLOAT) {
		guide_roi->fdata = malloc(npixels_roi * nchans * sizeof(float));
		if (!guide_roi->fdata)
			retval = 1;
		guide_roi->fpdata[0] = gui.roi.fit.fdata;
		guide_roi->fpdata[1] = rgb? guide_roi->fdata + npixels_roi : guide_roi->fdata;
		guide_roi->fpdata[2] = rgb? guide_roi->fdata + 2 * npixels_roi : guide_roi->fdata;
		for (uint32_t c = 0 ; c < nchans ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				float *srcindex = guide->fdata + (npixels * c) + ((guide->ry - y - gui.roi.selection.y) * guide->rx) + gui.roi.selection.x;
				float *destindex = guide_roi->fdata + (npixels_roi * c) + (guide_roi->rx * y);
				memcpy(destindex, srcindex, (gui.roi.selection.w) * sizeof(float));
			}
		}
	} else {
		guide_roi->data = malloc(npixels_roi * nchans * sizeof(WORD));
		if (!guide_roi->data)
			retval = 1;
		guide_roi->pdata[0] = gui.roi.fit.data;
		guide_roi->pdata[1] = rgb? guide_roi->data + npixels_roi : guide_roi->data;
		guide_roi->pdata[2] = rgb? guide_roi->data + 2 * npixels_roi : guide_roi->data;
		for (uint32_t c = 0 ; c < nchans ; c++) {
			for (uint32_t y = 0; y < gui.roi.selection.h ; y++) {
				WORD *srcindex = guide->data + (npixels * c) + ((guide->ry - y - gui.roi.selection.y) * guide->rx) + gui.roi.selection.x;
				WORD *destindex = guide_roi->data + (npixels_roi * c) + (guide_roi->rx * y);
				memcpy(destindex, srcindex, (gui.roi.selection.w) * sizeof(WORD));
			}
		}
	}
	return retval;
}

int edge_preserving_filter(struct epfargs *args) {
	fits *fit = args->fit;
	fits *guide = args->guidefit;
	double d = args->d;
	double sigma_col = args->sigma_col;
	double sigma_space = args->sigma_space;
	double mod = args->mod;
	gboolean verbose = args->verbose;
	ep_filter_t filter_type = args->filter;

	struct timeval t_start, t_end;
	if (sigma_col <= 0.0 || (sigma_space <= 0.0 && filter_type == EP_BILATERAL))
		return 1;
	if (verbose) {
		siril_log_color_message(_("Bilateral filter: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}
	sigma_col /= 100.0;

	if (fit->naxes[2] == 1) {
		// This makes the settings behave more consistently between color and mono images
		sigma_col /= 25.0;
		if (filter_type == EP_GUIDED) {
			sigma_col /= 10.0;
		}
	}
	if (filter_type == EP_GUIDED) {
		// This makes the settings behave more consistently between the two filter types
		sigma_col /= 5.0;
		if (d == 0)
			d = sigma_space;
		d /= 3.0;
	}
	if (fit->type == DATA_FLOAT) {
		// This makes the settings behave more consistently between 16-bit and 32-bit images
		sigma_col *= 2.0;
	}

	// cv::BilateralFilter() only works on 8u and 32f images, so we convert 16-bit to 32-bit
	size_t ndata = fit->rx * fit->ry * fit->naxes[2];
	data_type orig_type = fit->type;
	if (orig_type == DATA_USHORT) {
		fit_replace_buffer(fit, ushort_buffer_to_float(fit->data, ndata), DATA_FLOAT);
	}
	fits orig = { 0 }; // for use with modulation
	if (mod < (1.0 - DBL_EPSILON)) {
		copyfits(fit, &orig, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	}
	double eps = sigma_col * sigma_col;
	fits *guide_roi = NULL, *guidance = NULL;
	gboolean roi_fitting_needed;
	switch (filter_type) {
		case EP_BILATERAL:
			cvBilateralFilter(fit, d, eps, sigma_space);
			break;
		case EP_GUIDED:
			guide_roi = calloc(1, sizeof(fits));
			roi_fitting_needed = (fit == &gui.roi.fit && guide != &gui.roi.fit && gui.roi.active);
			if (roi_fitting_needed)
				match_guide_to_roi(guide, guide_roi);
			guidance = roi_fitting_needed ? guide_roi : guide;
			cvGuidedFilter(fit, guidance, d, eps);
			clearfits(guide_roi);
			free(guide_roi);
			break;
	}
	if (mod < (1.0 - DBL_EPSILON)) {
		for (size_t j = 0 ; j < fit->ry ; j++) {
			size_t offset = j * fit->rx;
			float modrem = 1.0 - mod;
			if (fit->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread) collapse(2) if (com.max_thread > 1)
#endif
				for (size_t k = 0 ; k < fit->naxes[2] ; k++) {
					for (size_t i = 0 ; i < fit->rx ; i++) {
						fit->fpdata[k][i + offset] = (float) mod * fit->fpdata[k][i + offset] + modrem * orig.fpdata[k][i + offset];
					}
				}
			} else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread) collapse(2) if (com.max_thread > 1)
#endif
				for (size_t k = 0 ; k < fit->naxes[2] ; k++) {
					for (size_t i = 0 ; i < fit->rx ; i++) {
						fit->pdata[k][i + offset] = roundf_to_WORD((float)(mod * fit->pdata[k][i + offset] + modrem * orig.pdata[k][i + offset]));
					}
				}
			}
		}
		clearfits(&orig);
	}
	if (orig_type == DATA_USHORT) {
		fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
	}

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}

	if (args->guide_needs_freeing)
		free(args->guidefit);
	free(args);
	return 0;
}
