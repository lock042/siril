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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "algos/sorting.h"
#include "algos/statistics.h"
#include "algos/siril_wcs.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "drizzle/cdrizzlebox.h"
#include "drizzle/cdrizzlemap.h"
#include "drizzle/cdrizzleutil.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "io/fits_keywords.h"
#include "registration/registration.h"
#include "algos/demosaicing.h"

#include "opencv/opencv.h"

#define MIN_RATIO 0.1 // minimum fraction of ref image dimensions to validate output sequence is worth creating

static int new_ref_index = -1;
static gboolean check_applyreg_output(struct registration_args *regargs);

static regdata *apply_reg_get_current_regdata(struct registration_args *regargs) {
	regdata *current_regdata;
	if (regargs->seq->regparam[regargs->layer]) {
		siril_log_message(
				_("Applying existing registration from layer #%d to transform the images\n"), regargs->layer);
		current_regdata = regargs->seq->regparam[regargs->layer];
	} else {
		siril_log_message(
				_("No registration data exists for this layer\n"));
		return NULL;
	}
	return current_regdata;
}

static void update_framing(regframe *framing, sequence *seq, int index) {
	if (seq->is_variable) {
		double rx2 = (double)seq->imgparam[index].rx;
		double ry2 = (double)seq->imgparam[index].ry;
		framing->pt[1].x = rx2 - 1.;
		framing->pt[2].x = rx2 - 1.;
		framing->pt[2].y = ry2 - 1.;
		framing->pt[3].y = ry2 - 1.;
	}
}

// This is the equivalent of opencv cv::detail::RotationWarper warpRoi
void compute_roi(Homography *H, int rx, int ry, framing_roi *roi) {
	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;
	Homography Href = { 0 };
	cvGetEye(&Href);
	regframe framing = (regframe){(point){0., 0.}, (point){(double)rx - 1., 0.}, (point){(double)rx - 1., (double)ry - 1.}, (point){0., (double)ry - 1.}};
	for (int j = 0; j < 4; j++) {
		cvTransfPoint(&framing.pt[j].x, &framing.pt[j].y, *H,  Href, 1.);
		if (xmin > framing.pt[j].x) xmin = framing.pt[j].x;
		if (ymin > framing.pt[j].y) ymin = framing.pt[j].y;
		if (xmax < framing.pt[j].x) xmax = framing.pt[j].x;
		if (ymax < framing.pt[j].y) ymax = framing.pt[j].y;
		// siril_debug_print("Point #%d: %3.2f %3.2f\n", j, framing.pt[j].x, framing.pt[j].y);
	}
	int x0 = (int)xmin;
	int y0 = (int)ymin;
	int w = (int)xmax - (int)xmin + 1;
	int h = (int)ymax - (int)ymin + 1;
	*roi = (framing_roi){ x0, y0, w, h };
}

static gboolean compute_framing(struct registration_args *regargs) {
	if (regargs->seq->number < 1)
		return FALSE;
	// validity of matrices has already been checked before this call
	// and null matrices have been discarded
	// Homography Href = regargs->seq->regparam[regargs->layer][regargs->reference_image].H;
	Homography Href = regargs->framingd.Htransf;
	// cvGetEye(&Href);
	Homography Hshift = { 0 };
	cvGetEye(&Hshift);
	int rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
	int ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
	int x0, y0, rx_0 = rx, ry_0 = ry, n;
	double xmin, xmax, ymin, ymax, cogx, cogy;

	// corners expressed in opencv conventions
	regframe framing = (regframe){(point){0., 0.}, (point){(double)rx - 1., 0.}, (point){(double)rx - 1., (double)ry - 1.}, (point){0., (double)ry - 1.}};
	gboolean retval = TRUE;

	switch (regargs->framing) {
		case FRAMING_CURRENT:
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i + 1);
				regframe current_framing = framing;
				update_framing(&current_framing, regargs->seq, i);
				double xs[4], ys[4];
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, Href, 1.);
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
					xs[j] = current_framing.pt[j].x;
					ys[j] = current_framing.pt[j].y;
				}
				quicksort_d(&xs[0], 4);
				quicksort_d(&ys[0], 4);
				// check we have overlap with the reference
				if (xs[3] < 0. || xs[0] > (double)rx || ys[3] < 0. || ys[0] > (double)ry) {
					siril_log_color_message(_("Image %d has no overlap with the reference\n"), "red", i + 1);
					retval = FALSE;
				}
			}
			break;
		case FRAMING_MAX:
			xmin = DBL_MAX;
			xmax = -DBL_MAX;
			ymin = DBL_MAX;
			ymax = -DBL_MAX;
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i);
				regframe current_framing = framing;
				update_framing(&current_framing, regargs->seq, i);
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, Href, 1.);
					if (xmin > current_framing.pt[j].x) xmin = current_framing.pt[j].x;
					if (ymin > current_framing.pt[j].y) ymin = current_framing.pt[j].y;
					if (xmax < current_framing.pt[j].x) xmax = current_framing.pt[j].x;
					if (ymax < current_framing.pt[j].y) ymax = current_framing.pt[j].y;
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
				}
			}
			rx_0 = (int)xmax - (int)xmin + 1;
			ry_0 = (int)ymax - (int)ymin + 1;
			x0 = (int)xmin;
			y0 = (int)ymin;
			siril_debug_print("new size: %d %d\n", rx_0, ry_0);
			siril_debug_print("new origin: %d %d\n", x0, y0);
			Hshift.h02 = (double)x0;
			Hshift.h12 = (double)y0;
			break;
		case FRAMING_MIN:
			xmin = -DBL_MAX;
			xmax = DBL_MAX;
			ymin = -DBL_MAX;
			ymax = DBL_MAX;
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i + 1);
				regframe current_framing = framing;
				update_framing(&current_framing, regargs->seq, i);
				double xs[4], ys[4];
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, Href, 1.);
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
					xs[j] = current_framing.pt[j].x;
					ys[j] = current_framing.pt[j].y;
				}
				quicksort_d(&xs[0], 4);
				quicksort_d(&ys[0], 4);
				if (xmin < xs[1]) xmin = xs[1];
				if (ymin < ys[1]) ymin = ys[1];
				if (xmax > xs[2]) xmax = xs[2];
				if (ymax > ys[2]) ymax = ys[2];
			}
			rx_0 = (int)xmax - (int)xmin + 1;
			ry_0 = (int)ymax - (int)ymin + 1;
			if (rx_0 < 0 || ry_0 < 0) {
				siril_log_color_message(_("The intersection of all images is null or negative\n"), "red");
				retval = FALSE;
			}
			x0 = (int)xmin;
			y0 = (int)ymin;
			siril_debug_print("new size: %d %d\n", rx_0, ry_0);
			siril_debug_print("new origin: %d %d\n", x0, y0);
			Hshift.h02 = (double)x0;
			Hshift.h12 = (double)y0;
			break;
		case FRAMING_COG:
			cogx = 0.;
			cogy = 0;
			n = 0;
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				int rx2 = (regargs->seq->is_variable) ? regargs->seq->imgparam[i].rx : regargs->seq->rx;
				int ry2 = (regargs->seq->is_variable) ? regargs->seq->imgparam[i].ry : regargs->seq->ry;
				double currcogx = (double)rx2 * 0.5;
				double currcogy = (double)ry2 * 0.5;
				cvTransfPoint(&currcogx, &currcogy,regargs->seq->regparam[regargs->layer][i].H, Href, 1.);
				cogx += currcogx;
				cogy += currcogy;
				n++;
			}
			if (n == 0)
				return FALSE; // If there are no frames filtered in, we have nothing to do.

			cogx /= (double)n;
			cogy /= (double)n;
			x0 = (int)(cogx - (double)rx * 0.5 - 0.5);
			y0 = (int)(cogy - (double)ry * 0.5 - 0.5);
			siril_log_message(_("Framing: Shift from reference origin: %d, %d\n"), x0, y0);
			Hshift.h02 = (double)x0;
			Hshift.h12 = (double)y0;
			Homography newHref = { 0 };
			cvMultH(Href, Hshift, &newHref);
			// we now check overlaps
			for (int i = 0; i < regargs->seq->number; i++) {
				if (!regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
					continue;
				siril_debug_print("Image #%d:\n", i + 1);
				regframe current_framing = framing;
				update_framing(&current_framing, regargs->seq, i);
				double xs[4], ys[4];
				for (int j = 0; j < 4; j++) {
					cvTransfPoint(&current_framing.pt[j].x, &current_framing.pt[j].y,regargs->seq->regparam[regargs->layer][i].H, newHref, 1.);
					siril_debug_print("Point #%d: %3.2f %3.2f\n", j, current_framing.pt[j].x, current_framing.pt[j].y);
					xs[j] = current_framing.pt[j].x;
					ys[j] = current_framing.pt[j].y;
				}
				quicksort_d(&xs[0], 4);
				quicksort_d(&ys[0], 4);
				// check we have overlap with the reference
				if (xs[3] < 0. || xs[0] > (double)rx || ys[3] < 0. || ys[0] > (double)ry) {
					siril_log_color_message(_("Image %d has no overlap with the reference\n"), "red", i + 1);
					retval = FALSE;
				}
			}
			break;
		default:
			return FALSE;
	}
	cvMultH(Href, Hshift, &regargs->framingd.Htransf);
	regargs->framingd.roi_out.w = (regargs->output_scale != 1.f) ? (int)(roundf((float)rx_0 * regargs->output_scale)) : rx_0;
	regargs->framingd.roi_out.h = (regargs->output_scale != 1.f) ? (int)(roundf((float)ry_0 * regargs->output_scale)) : ry_0;
	regargs->framingd.Hshift = Hshift;

	return retval;
}

// For framing max, we don't want to export the image with black borders
// From the homographies and the input image size
// we compute:
// - H transf, the warping transformation in place
// - H shift, the residual shift wrt to ref image
// - the size of the transformed image after applying the transformation
static void compute_Hmax(Homography *Himg, Homography *Href, int src_rx_in, int src_ry_in, double scale, Homography *H, Homography *Hshift, int *dst_rx_out, int *dst_ry_out) {
	*dst_rx_out = 0;
	*dst_ry_out = 0;
	regframe framing = { 0 };
	framing = (regframe){(point){0., 0.}, (point){(double)src_rx_in - 1., 0.}, (point){(double)src_rx_in - 1., (double)src_ry_in - 1.}, (point){0., (double)src_ry_in - 1.}};
	double xmin, xmax, ymin, ymax;
	xmin = DBL_MAX;
	xmax = -DBL_MAX;
	ymin = DBL_MAX;
	ymax = -DBL_MAX;
	// cvGetEye(Href);
	for (int j = 0; j < 4; j++) {
		cvTransfPoint(&framing.pt[j].x, &framing.pt[j].y, *Himg, *Href, scale);
		if (xmin > framing.pt[j].x) xmin = framing.pt[j].x;
		if (ymin > framing.pt[j].y) ymin = framing.pt[j].y;
		if (xmax < framing.pt[j].x) xmax = framing.pt[j].x;
		if (ymax < framing.pt[j].y) ymax = framing.pt[j].y;
		siril_debug_print("Point #%d: %3.2f %3.2f\n", j, framing.pt[j].x, framing.pt[j].y);
	}
	*dst_rx_out = (int)xmax - (int)xmin + 1;
	*dst_ry_out = (int)ymax - (int)ymin + 1;
	Hshift->h02 = round(xmin);
	Hshift->h12 = round(ymin);
	cvTransfH(Himg, Href, H);
	// the shift matrix is at the final scale, while the transformation matrix
	// is still at the orginal scale (will be upscaled in cvTransformImage)
	Homography Hcorr = { 0 };
	cvGetEye(&Hcorr);
	Hcorr.h02 = -Hshift->h02 / scale;
	Hcorr.h12 = -Hshift->h12 / scale;
	cvMultH(Hcorr, *H, H);
}

int apply_reg_prepare_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;

	fits fit = { 0 };

	/* preparing reference data from reference fit and making sanity checks*/
	sadata->current_regdata = apply_reg_get_current_regdata(regargs);
	if (!sadata->current_regdata) return -2;


	if (seq_read_frame_metadata(args->seq, regargs->reference_image, &fit)) {
		siril_log_color_message(_("Could not load reference image\n"), "red");
		args->seq->regparam[regargs->layer] = NULL;
		clearfits(&fit);
		return 1;
	}
	if (!regargs->driz && fit.naxes[2] == 1 && fit.keywords.bayer_pattern[0] != '\0')
		siril_log_color_message(_("Applying transformation on a sequence opened as CFA is a bad idea.\n"), "red");

	if (regargs->undistort) {
		siril_log_message(_("Distortion data was found in the sequence file, undistortion will be applied\n"));
	}

		// We prepare the distortion structure maps if required
	if (regargs->undistort && init_disto_map(fit.rx, fit.ry, regargs->disto)) {
		siril_log_color_message(
				_("Could not init distortion mapping\n"), "red");
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		clearfits(&fit);
		return 1;
	}
	clearfits(&fit);
	return registration_prepare_results(args);
}

/* reads the image and apply existing transformation */
int apply_reg_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	struct driz_args_t *driz = regargs->driz;
	float scale = regargs->output_scale;
	struct driz_param_t *p = NULL;
	disto_data *disto = NULL;

	Homography H = { 0 };
	Homography Himg = { 0 };
	Homography Hs = { 0 };
	cvGetEye(&Hs);
	int dst_rx = regargs->framingd.roi_out.w;
	int dst_ry = regargs->framingd.roi_out.h;
	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes

	if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
		siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
	}

	// Composing transformation wrt reference image
	Himg = regargs->seq->regparam[regargs->layer][in_index].H;
	if (guess_transform_from_H(Himg) == NULL_TRANSFORMATION)
		return 1; // in case H is null and -selected was not passed

	if (regargs->framing != FRAMING_MAX)
		cvTransfH(&Himg, &regargs->framingd.Htransf, &H);
	else {
		compute_Hmax(&Himg, &regargs->framingd.Htransf, fit->rx, fit->ry, scale, &H, &Hs, &dst_rx, &dst_ry);
	}

	if (regargs->undistort) {
		if (regargs->disto->dtype == DISTO_MAP_D2S || regargs->disto->dtype == DISTO_MAP_S2D) {
			disto = regargs->disto;
		} else {
			disto = &regargs->disto[in_index];
		}
	}
	if (has_wcs(fit)) {
		free_wcs(fit); // we remove the current solution in all cases
		if (regargs->wcsref) { // we will update it only if ref image has a solution
			fit->keywords.wcslib = wcs_deepcopy(regargs->wcsref, NULL); // we copy the reference astrometry
			Homography Hscale = { 0 }, Hshift = { 0 };
			if (regargs->distoparam.velocity.x != 0.f || regargs->distoparam.velocity.y != 0.f) {
				cvGetEye(&Hscale);
				pointf reg = { 0.f, 0.f };
				get_comet_shift(regargs->reference_date, fit->keywords.date_obs, regargs->distoparam.velocity, &reg);
				Hshift.h02 = -reg.x;
				Hshift.h12 =  reg.y;
				cvApplyFlips(&Hscale, dst_ry / scale, dst_ry);
				reframe_wcs(fit->keywords.wcslib, &Hshift);
			}
			if (regargs->framing == FRAMING_MAX) {
				cvGetEye(&Hshift);
				Hshift = Hs;
				// Hshift is the correction at full scale, we need to de-scale it
				Hshift.h02 /= scale;
				Hshift.h12 /= scale;
				Homography Href = regargs->framingd.Hshift;
				cvMultH(Href, Hshift, &Hshift);
				Hshift.h02 *= -1.;
				// we use (int)(dst_ry / scale) to find the projected area origin size unscaled
				// we also need to recast fit->ry bec. it's a uint
				Hshift.h12 += (double)(int)(dst_ry / scale) - (double)fit->ry;
				reframe_wcs(fit->keywords.wcslib, &Hshift);
			}
			if (scale != 1.f) {
				cvGetEye(&Hscale);
				Hscale.h00 = regargs->output_scale;
				Hscale.h11 = regargs->output_scale;
				cvApplyFlips(&Hscale, dst_ry / scale, dst_ry);
				reframe_wcs(fit->keywords.wcslib, &Hscale);
			}
		}
	}

	if (regargs->driz) {
		p = calloc(1, sizeof(struct driz_param_t));
		driz_param_init(p);
		p->kernel = driz->kernel;
		p->driz = driz;
		p->error = malloc(sizeof(struct driz_error_t));
		p->scale = scale;
		p->pixel_fraction = driz->pixel_fraction;
		p->cfa = driz->cfa;
		// Set bounds equal to whole image
		p->xmin = p->ymin = 0;
		p->xmax = fit->rx - 1;
		p->ymax = fit->ry - 1;
		p->pixmap = calloc(1, sizeof(imgmap_t));
		p->pixmap->rx = fit->rx;
		p->pixmap->ry = fit->ry;
		p->threads = threads;

		if (map_image_coordinates_h(fit, H, p->pixmap, dst_rx, dst_ry, scale, disto, threads)) {
			free(p->error);
			free(p->pixmap);
			free(p);
			return 1;
		}

		if (!p->pixmap->xmap) {
			siril_log_color_message(_("Error generating mapping array.\n"), "red");
			free(p->error);
			free(p->pixmap);
			free(p);
			return 1;
		}
		p->data = fit;
		// Convert fit to 32-bit float if required
		float *newbuf = NULL;
		if (fit->type == DATA_USHORT) {
			siril_debug_print("Replacing ushort buffer for drizzling\n");
			size_t ndata = fit->rx * fit->ry * fit->naxes[2];
			newbuf = malloc(ndata * sizeof(float));
			float invnorm = 1.f / USHRT_MAX_SINGLE;
			for (size_t i = 0 ; i < ndata ; i++) {
				newbuf[i] = fit->data[i] * invnorm;
			}
			fit_replace_buffer(fit, newbuf, DATA_FLOAT);
		}
		/* Set up output fits */
		fits out = { 0 };
		copyfits(fit, &out, CP_FORMAT, -1);
		// copy the DATE_OBS
		out.keywords.date_obs = g_date_time_ref(fit->keywords.date_obs);
		out.rx = out.naxes[0] = dst_rx;
		out.ry = out.naxes[1] = dst_ry;
		out.naxes[2] = driz->is_bayer ? 3 : 1;
		size_t chansize = out.rx * out.ry * sizeof(float);
		out.fdata = calloc(out.naxes[2] * chansize, 1);
		out.fpdata[0] = out.fdata;
		out.fpdata[1] = out.naxes[2] == 1 ? out.fdata : out.fdata + chansize;
		out.fpdata[2] = out.naxes[2] == 1 ? out.fdata : out.fdata + 2 * chansize;
		p->output_data = &out;

		// Set up the output_counts fits to store pixel hit counts
		fits *output_counts = calloc(1, sizeof(fits));
		copyfits(&out, output_counts, CP_FORMAT, -1);
		output_counts->fdata = calloc(output_counts->rx * output_counts->ry * output_counts->naxes[2], sizeof(float));
		p->output_counts = output_counts;

		p->weights = driz->flat;

		if (dobox(p)) { // Do the drizzle
			siril_log_color_message("s\n", p->error->last_message);
			return 1;
		}
		clearfits(fit);
		// copy the DATE_OBS
		copyfits(&out, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		fit->keywords.date_obs = g_date_time_ref(out.keywords.date_obs);
		clearfits(&out);
		if (args->seq->type == SEQ_SER || com.pref.force_16bit) {
			fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, fit->rx * fit->ry * fit->naxes[2]), DATA_USHORT);
		}
		if (driz->is_bayer) {
			/* we need to do something special here because it's a 1-channel sequence and
			* image that will become 3-channel after this call, so stats caching will
			* not be correct. Preprocessing does not compute new stats so we don't need
			* to save them into cache now.
			* Destroying the stats from the fit will also prevent saving them in the
			* sequence, which may be done automatically by the caller.
			*/
			full_stats_invalidation_from_fit(fit);
			fit->keywords.lo = 0;
			clear_Bayer_information(fit); // we also reset the bayerpattern
		}

		free(p->pixmap->xmap);
		free(p->pixmap);

		// Get rid of the output_counts image, no longer required
		clearfits(output_counts);
		free(output_counts);
		output_counts = NULL;
	}
	else {
		if (regargs->interpolation <= OPENCV_LANCZOS4) {
			if (cvTransformImage(fit, dst_rx, dst_ry, H, scale, regargs->interpolation, regargs->clamp, disto)) {
				return 1;
			}
		} else {
			if (shift_fit_from_reg(fit, H)) {
				return 1;
			}
		}
	}
	if (in_index == regargs->reference_image)
		new_ref_index = out_index; // keeping track of the new ref index in output sequence

	regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
	regargs->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	regargs->imgparam[out_index].rx = dst_rx;
	regargs->imgparam[out_index].ry = dst_ry;
	regargs->regparam[out_index].fwhm = sadata->current_regdata[in_index].fwhm;
	regargs->regparam[out_index].weighted_fwhm = sadata->current_regdata[in_index].weighted_fwhm;
	regargs->regparam[out_index].roundness = sadata->current_regdata[in_index].roundness;
	regargs->regparam[out_index].background_lvl = sadata->current_regdata[in_index].background_lvl;
	regargs->regparam[out_index].number_of_stars = sadata->current_regdata[in_index].number_of_stars;
	regargs->regparam[out_index].H = Hs;

	if (scale != 1.f) {
		fit->keywords.pixel_size_x /= scale;
		fit->keywords.pixel_size_y /= scale;
		regargs->regparam[out_index].fwhm *= scale;
		regargs->regparam[out_index].weighted_fwhm *= scale;
	}
	if (has_wcs(fit))
		update_wcsdata_from_wcs(fit); // we need to do it when image dimensions have been changed

	sadata->success[out_index] = 1;
	return 0;
}

int apply_reg_finalize_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int failed = 0;

	// images may have been excluded but selnum wasn't updated
	fix_selnum(args->seq, FALSE);

	if (!args->retval) {
		for (int i = 0; i < args->nb_filtered_images; i++)
			if (!sadata->success[i])
				failed++;
		regargs->new_total = args->nb_filtered_images - failed;
		if (failed) {
			// regargs->imgparam and regargs->regparam may have holes caused by images
			// that failed to be registered - compact them
			for (int i = 0, j = 0; i < regargs->new_total; i++, j++) {
				while (!sadata->success[j] && j < args->nb_filtered_images) j++;
				g_assert(sadata->success[j]);
				if (i != j) {
					regargs->imgparam[i] = regargs->imgparam[j];
					regargs->regparam[i] = regargs->regparam[j];
				}
			}
		}
		seq_finalize_hook(args);
	} else {
		regargs->new_total = 0;

		// same as seq_finalize_hook but with file deletion
		if ((args->force_ser_output || args->seq->type == SEQ_SER) && args->new_ser) {
			ser_close_and_delete_file(args->new_ser);
			free(args->new_ser);
		} else if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
			fitseq_close_and_delete_file(args->new_fitseq);
			free(args->new_fitseq);
		} else if (args->seq->type == SEQ_REGULAR) {
			remove_prefixed_sequence_files(regargs->seq, regargs->prefix);
		}
	}

	if (sadata->success)
		free(sadata->success);
	free(sadata);
	args->user = NULL;

	if (!args->retval) {
		siril_log_message(_("Applying registration completed.\n"));
		gchar *str = ngettext("%d image processed.\n", "%d images processed.\n", args->nb_filtered_images);
		str = g_strdup_printf(str, args->nb_filtered_images);
		siril_log_color_message(str, "green");
		siril_log_color_message(_("Total: %d failed, %d exported.\n"), "green", failed, regargs->new_total);

		g_free(str);
		if (!(args->seq->type == SEQ_INTERNAL)) {
			// explicit sequence creation to copy imgparam and regparam
			create_output_sequence_for_registration(regargs, new_ref_index);
			// will be loaded in the idle function if (load_new_sequence)
			regargs->load_new_sequence = TRUE; // only case where a new sequence must be loaded
		}
	}
	else {
		siril_log_message(_("Transformation aborted.\n"));
	}
	return regargs->new_total == 0;
}

/* this function is called by both:
- apply_reg_compute_mem_limits
- star_align_compute_mem_limits
*/
int apply_reg_compute_mem_consumption(struct generic_seq_args *args, unsigned int *total_required_per_image_MB, unsigned int *required_per_dst_image_MB, unsigned int *max_mem_MB) {
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail, MB_per_mono_float_orig_image, MB_per_mono_float_scaled_image;
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;

	/* The apply reg memory consumption is:

		* the original image
		* the transformed image, including force_float, scale and x 3 for bayer drizzle if required
			Note: If drizzle, we should always assume force_float as the output of drizzle is always float (realloc'ed to 16b if set in pref)

		*** Maps ***
		* For interpolation, if undistortion is included:
			- two 32b mono images per image for its re-mapping
		* For drizzle, in all cases:
			- two 32b mono images per image for its re-mapping
		Note: the size of the maps depends of the transformation method:
			- destination size for Interpolation (so mono float upscaled)
			- source size for Drizzle (so mono float original)

		*** Drizzle specifics ***
		* the output counts image, the same size as the scaled image

		*** Interpolation specifics ***
		* If color, the consumption for openCV is trickier:
		* we have at some point:
			- 2 * orig (occurs during fits_to_bgrbgr_xxx)
			- 2 * dest - orig (occurs during Mat_to_image)
		* We should compare to orig + dest
			- if orig > dest (scale < 1.), we have at worst 2 * orig instead of orig + dest
			- if dest > orig (scale > sqrt(2)), we have at worst 2 * dest - orig instead of orig + dest
		* If clamping is enabled, one scaled image (the guide) and a 8b copy (the mask)

		*** Others not to be accounted for:
		- if Drizzle and master flat is specified, because drizzle init loads it
			before we compute available memory
		- If Interpolation and undistorion with pre-computed maps (DISTO_MAP_D2S or DISTO_MAP_S2D),
			those two maps are init before this hook

		TODO: there are a few approximations:
		- using seq->rx/seq->ry for original image size while for interpolation, the seq can be variable (not for drizzle though, drizzle is not allowed)
		- assuming the output scaled image is rx*ry*scale^2. This does not work for min or max framing. min is safe, max is not...

	*/

	// using drizzle, we force memory consumption as float
	gboolean force_float = args->force_float || regargs->driz;

	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;

	MB_per_mono_float_orig_image = max(1, (unsigned int)(args->seq->rx * args->seq->ry * sizeof(float) / BYTES_IN_A_MB));
	MB_per_mono_float_scaled_image = max(1, (unsigned int)(args->seq->rx * args->seq->ry * args->upscale_ratio * args->upscale_ratio * sizeof(float) / BYTES_IN_A_MB));

	if (regargs->driz && regargs->driz->is_bayer)
		MB_per_scaled_image *= 3; // we have a mono in and a color out

	// image in and out
	unsigned int required = MB_per_orig_image + MB_per_scaled_image;

	// maps
	if (regargs->driz)
		required += 2 * MB_per_mono_float_orig_image; // maps
	else if (regargs->undistort)
		required += 2 * MB_per_mono_float_scaled_image; // maps if undistortion, otherwise we directly use openCV warpPerspective() which computes maps iteratively on 32x32 chunks

	// drizzle specifics
	if (regargs->driz)
		required += MB_per_scaled_image;

	// interpolation specifics
	if (!regargs->driz && regargs->seq->nb_layers == 3) {
		if (regargs->output_scale < 1.f)
			required += MB_per_orig_image - MB_per_scaled_image; // we had orig + dest, now we have 2 * orig
		else if (regargs->output_scale > (float)SQRT2)
			required += MB_per_scaled_image - 2 * MB_per_orig_image; // we had orig + dest, now we have 2 * dest - orig
	}
	if (!regargs->driz && regargs->clamp && (regargs->interpolation == OPENCV_CUBIC ||
			regargs->interpolation == OPENCV_LANCZOS4)) {
		float factor = (is_float) ? 0.25 : 0.5;
		required += (1 + factor) * MB_per_scaled_image; // we need one scaled image (the guide) and a 8b copy (the mask)
	}

	// storing the returned values
	*total_required_per_image_MB = required;
	*max_mem_MB = MB_avail;
	*required_per_dst_image_MB = MB_per_scaled_image;

	limit = (int)(MB_avail / required);
	return limit;
}

int apply_reg_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int required = 0, MB_avail = 0, MB_per_scaled_image = 0;
	int limit = apply_reg_compute_mem_consumption(args, &required, &MB_per_scaled_image, &MB_avail);

	if (limit > 0) {

		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_scaled_image;
		} else limit = thread_limit;
	}

	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per thread, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_scaled_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

// Used by both apply_reg and global methods: definition is in registration.h
int initialize_drizzle_params(struct generic_seq_args *args, struct registration_args *regargs) {
	struct driz_args_t *driz = regargs->driz;
	set_progress_bar_data(_("Initializing drizzle data..."), PROGRESS_PULSATE);
	driz->scale = regargs->output_scale;
	driz_param_dump(driz); // Print some info to the log
	/* preparing reference data from reference fit and making sanity checks*/
	fits fit = { 0 };

	/* fit will now hold the reference frame */
	if (seq_read_frame_metadata(regargs->seq, regargs->reference_image, &fit)) {
		siril_log_message(_("Could not load reference image\n"));
	}
	sensor_pattern pattern;
	if (args->seq->type == SEQ_SER && args->seq->ser_file ) {
		driz->is_bayer = TRUE;
		switch (args->seq->ser_file->color_id) {
			case SER_MONO:
				pattern = get_cfa_pattern_index_from_string("");
				driz->is_bayer = FALSE;
				break;
			case SER_BAYER_RGGB:
				pattern = get_cfa_pattern_index_from_string("RGGB");
				break;
			case SER_BAYER_GRBG:
				pattern = get_cfa_pattern_index_from_string("GRBG");
				break;
			case SER_BAYER_GBRG:
				pattern = get_cfa_pattern_index_from_string("GBRG");
				break;
			case SER_BAYER_BGGR:
				pattern = get_cfa_pattern_index_from_string("BGGR");
				break;
			default:
				siril_log_message(_("Unsupported SER CFA pattern detected. Treating as mono.\n"));
				driz->is_bayer = FALSE;
		}
	} else {
		const gchar bayertest = fit.keywords.bayer_pattern[0];
		driz->is_bayer = (bayertest == 'R' || bayertest == 'G' || bayertest == 'B'); // If there is a CFA pattern we need to CFA drizzle
		pattern = com.pref.debayer.use_bayer_header ? get_cfa_pattern_index_from_string(fit.keywords.bayer_pattern) : com.pref.debayer.bayer_pattern;
	}
	if (driz->is_bayer) {
		adjust_Bayer_pattern(&fit, &pattern);
		driz->cfa = get_cfa_from_pattern(pattern);
		if (!driz->cfa) // if fit.bayer_pattern exists and get_cfa_from_pattern returns NULL then there is a problem.
			return 1;
	}
	return 0;
}

int register_apply_reg(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	args->force_float = !com.pref.force_16bit && regargs->seq->type != SEQ_SER;
	control_window_switch_to_tab(OUTPUT_LOGS);
	int retval = 0;
	if (!check_before_applyreg(regargs)) { // checks for input arguments wrong combinations
		retval = -1;
		goto END;
	}

	// If this is an astrometric aligned sequence,
	// we need to collect the WCS structures and
	// recompute the homographies based on selected frames
	// We will also unselected unsolved images before recomputing the filters
	regargs->undistort = (layer_has_distortion(regargs->seq, regargs->layer)) ? regargs->seq->distoparam[regargs->layer].index : DISTO_UNDEF;
	if (regargs->undistort == DISTO_FILES) {
		regargs->WCSDATA = calloc(regargs->seq->number, sizeof(struct wcsprm));
		if (collect_sequence_astrometry(regargs)) {
			retval = -1;
			goto END;
		}
		Homography Href = { 0 };
		if (compute_Hs_from_astrometry(regargs->seq, regargs->WCSDATA, regargs->framing, regargs->layer, &Href, &regargs->wcsref)) {
			retval = -1;
			goto END;
		}
		regargs->framingd.Htransf = Href;
	} else {
		regargs->framingd.Htransf = regargs->seq->regparam[regargs->layer][regargs->reference_image].H;
	}

	// we need to update filtering before check_applyreg_output that will check required space
	if (!regargs->filtering_criterion &&
				convert_parsed_filter_to_filter(&regargs->filters,
				regargs->seq, &regargs->filtering_criterion,
				&regargs->filtering_parameter)) {
		retval = -1;
		goto END;

	}

	// Prepare sequence filtering
	int nb_frames = compute_nb_filtered_images(regargs->seq,
			regargs->filtering_criterion, regargs->filtering_parameter);
	regargs->new_total = nb_frames;	// to avoid recomputing it later
	gchar *str = describe_filter(regargs->seq, regargs->filtering_criterion,
			regargs->filtering_parameter);
	siril_log_message(str);
	g_free(str);

	// We can now compute the framing and check the output size
	if (!check_applyreg_output(regargs)) {
		retval = -1;
		goto END;

	}

	if (regargs->no_output) {
		retval = 0;
		goto END;

	}

	if (regargs->driz && initialize_drizzle_params(args, regargs)) {
		retval = -1;
		goto END;

	}

	args->upscale_ratio = regargs->output_scale;
	args->prepare_hook = apply_reg_prepare_hook;
	args->compute_mem_limits_hook = apply_reg_compute_mem_limits;
	args->finalize_hook = apply_reg_finalize_hook;
	args->filtering_criterion = regargs->filtering_criterion;
	args->filtering_parameter = regargs->filtering_parameter;
	args->nb_filtered_images = regargs->new_total;
	args->image_hook = apply_reg_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Apply registration");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = regargs->prefix ? strdup(regargs->prefix) : NULL;
	args->load_new_sequence = TRUE;
	args->already_in_a_thread = TRUE;

	struct star_align_data *sadata = calloc(1, sizeof(struct star_align_data));
	if (!sadata) {
		retval = -1;
		goto END;

	}
	sadata->regargs = regargs;
	args->user = sadata;

	disto_source index = DISTO_UNDEF;

	if (regargs->undistort) {
		regargs->distoparam = regargs->seq->distoparam[regargs->layer];
		int status = 1;
		index = regargs->distoparam.index; // to keep track if DISTO_FILES as if no distorsion is effectively present, it will be set to UNDEF
		regargs->disto = init_disto_data(&regargs->distoparam, regargs->seq, regargs->WCSDATA, regargs->driz != NULL, &status);
		free(regargs->WCSDATA); // init_disto_data has freed each individual wcs, we can now free the array
		if (status) {
			free_generic_seq_args(args, FALSE);
			siril_log_color_message(_("Could not initialize distortion data, aborting\n"), "red");
			free(sadata);
			args->user = NULL;
			retval = 1;
			goto END;
		}
		if (!regargs->disto) {
			regargs->undistort = DISTO_UNDEF;
		}
	}
	if (!regargs->wcsref)
		regargs->wcsref = get_wcs_ref(regargs->seq);
	if (regargs->wcsref) {
		if (regargs->undistort && regargs->wcsref->lin.dispre)
			remove_dis_from_wcs(regargs->wcsref); // we remove distortions as the output is undistorted
		if (!regargs->undistort && regargs->wcsref->lin.dispre)
			siril_log_color_message(_("Distortion was found in reference image astrometry but you did not include distortion correction when registering the images\n"), "salmon");
		if (index == DISTO_FILES && image_is_flipped_from_wcs(regargs->wcsref)) { // we are in astrometric reg, we will need to flip the solution if required
			Homography H = { 0 };
			cvGetEye(&H);
			H.h11 = -1.;
			H.h12 = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
			reframe_wcs(regargs->wcsref, &H);
		}
		// we can't apply same for max framing as the individual image size is different from regargs->framingd.roi_out.h / scale
		if (regargs->framing == FRAMING_MIN || regargs->framing == FRAMING_COG) {
			Homography H = regargs->framingd.Hshift;
			cvInvertH(&H);
			int orig_ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
			cvApplyFlips(&H, orig_ry, regargs->framingd.roi_out.h / regargs->output_scale);
			reframe_wcs(regargs->wcsref, &H);
		}
	}
	if (regargs->seq->distoparam[regargs->layer].index == DISTO_FILE_COMET) {
		// we fetch the reference date to compute back comet shifts
		fits ref = { 0 };
		if (seq_read_frame_metadata(args->seq, regargs->reference_image, &ref)) {
			siril_log_message(_("Could not load reference image\n"));
			free(sadata);
			args->user = NULL;
			retval = -1;
			goto END;
		}
		regargs->reference_date = g_date_time_ref(ref.keywords.date_obs);
		clearfits(&ref);
	}

	generic_sequence_worker(args);
	regargs->retval = args->retval;
	retval = regargs->retval;

END:
	free_generic_seq_args(args, args->seq->type != SEQ_INTERNAL);

	return retval;
}

transformation_type guess_transform_from_H(Homography H) {
	if (fabs(H.h00 + H.h01 + H.h02 + H.h10 + H.h11 + H.h12 + H.h20 + H.h21 + H.h22) < __DBL_EPSILON__)
		return NULL_TRANSFORMATION;
	if (fabs(H.h20) > __DBL_EPSILON__ || fabs(H.h21) > __DBL_EPSILON__)
		return HOMOGRAPHY_TRANSFORMATION;
	if (fabs(H.h00 - 1.) < __DBL_EPSILON__ && fabs(H.h11 - 1.) < __DBL_EPSILON__ &&
			fabs(H.h10) < __DBL_EPSILON__ && fabs(H.h01) < __DBL_EPSILON__) {
		if (fabs(H.h02) > __DBL_EPSILON__ || fabs(H.h12) > __DBL_EPSILON__)
			return SHIFT_TRANSFORMATION;
		return IDENTITY_TRANSFORMATION;
	}
	if (fabs(H.h10 - H.h00  + H.h01 + H.h11) < __DBL_EPSILON__)
		return SIMILARITY_TRANSFORMATION;
	return AFFINE_TRANSFORMATION;
}

void guess_transform_from_seq(sequence *seq, int layer,
		transformation_type *min, transformation_type *max, gboolean excludenull) {
	*min = HOMOGRAPHY_TRANSFORMATION; // highest value
	*max = UNDEFINED_TRANSFORMATION;  // lowest value
	gboolean needs_sel_update = FALSE;

	if (!layer_has_registration(seq, layer)) {
		siril_debug_print("No registration data found in sequence or layer\n");
		return;
	}
	for (int i = 0; i < seq->number; i++){
		transformation_type val = guess_transform_from_H(seq->regparam[layer][i].H);
		//siril_debug_print("Image #%d - transf = %d\n", i+1, val);
		if (*max < val) *max = val;
		if (*min > val) *min = val;
		if (val == NULL_TRANSFORMATION && excludenull) {
			seq->imgparam[i].incl = FALSE;
			needs_sel_update = TRUE;
		}
	}
	if (excludenull && needs_sel_update)
		fix_selnum(seq, FALSE);
}

static int confirm_exceed_cairomaxdim(gpointer user_data) {
	struct registration_args *regargs = (struct registration_args*)user_data;
	gchar *msg = g_strdup(_("Images will be larger than what Siril can display for now. "
		"Tune scale ratio to get max dimension smaller than 32767 pixels "
		"or proceed and post process your images with an external program."));
	if (regargs->no_output) { // Estimate button was pressed, we just warn
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Output too large"), msg);
	} else if (!siril_confirm_dialog(_("Output too large"),  // Register button was pressed, we ask for confirmation
				msg, _("Proceed"))) {
		regargs->retval = 1;
	}
	g_free(msg);
	return 0;
}

static gboolean check_applyreg_output(struct registration_args *regargs) {
	/* compute_framing uses the filtered list of images, so we compute the filter here */

	// determines the reference homography (including framing shift) and output size
	if (!compute_framing(regargs)) {
		siril_log_color_message(_("Unselect the images generating the error or change framing method to max\n"), "red");
		return FALSE;
	}

	// TODO: temp check for final image size
	// if larger than cairo image buffer, pop a warning that image will not display at all
	int max_dim = max(regargs->framingd.roi_out.w, regargs->framingd.roi_out.h);
	if (max_dim > 32767) {
		if (!(com.script || com.python_script)) { // trhough GUI, we warn with GTK objects
			execute_idle_and_wait_for_it(confirm_exceed_cairomaxdim, regargs);
			if (regargs->retval)
				return FALSE;
		} else {
			siril_log_color_message(_("Images will be larger than what Siril can display for now."), "salmon");
			siril_log_color_message(_("Tune scale ratio to get max dimension smaller than 32767 pixels"), "salmon");
		}
	}

	// make sure we apply registration only if the output sequence has a meaningful size
	int rx0 = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].rx : regargs->seq->rx;
	int ry0 = (regargs->seq->is_variable) ? regargs->seq->imgparam[regargs->reference_image].ry : regargs->seq->ry;
	if (regargs->framingd.roi_out.w < rx0 * MIN_RATIO || regargs->framingd.roi_out.h < ry0 * MIN_RATIO) {
		siril_log_color_message(_("The output sequence is too small compared to reference image (too much rotation or little overlap?)\n"), "red");
		siril_log_color_message(_("You should change framing method, aborting\n"), "red");
		return FALSE;
	}

	int nb_frames = regargs->seq->selnum;
	// cannot use seq_compute_size as rx_out/ry_out are not necessarily consistent with seq->rx/ry
	int64_t size = (int64_t) regargs->framingd.roi_out.w * regargs->framingd.roi_out.h * regargs->seq->nb_layers;
	if (regargs->seq->type == SEQ_SER) {
		size *= regargs->seq->ser_file->byte_pixel_depth;
		size *= nb_frames;
		size += SER_HEADER_LEN;
	} else {
		size *= (get_data_type(regargs->seq->bitpix) == DATA_USHORT) ? sizeof(WORD) : sizeof(float);
		size += FITS_DOUBLE_BLOC_SIZE; // FITS double HDU size
		size *= nb_frames;
	}
	gchar* size_msg = g_format_size_full(size, G_FORMAT_SIZE_IEC_UNITS);
	siril_debug_print("Apply Registration: sequence out size: %s\n", size_msg);
	g_free(size_msg);
	if (test_available_space(size)) {
		siril_log_color_message(_("Not enough space to save the output images\n"), "red");
		return FALSE;
	}
	return TRUE;
}

gboolean check_before_applyreg(struct registration_args *regargs) {
	struct driz_args_t *driz = regargs->driz;
		// check the reference image matrix is not null
	transformation_type checkH = guess_transform_from_H(regargs->seq->regparam[regargs->layer][regargs->seq->reference_image].H);
	if (checkH == NULL_TRANSFORMATION) {
		siril_log_color_message(_("The reference image has a null matrix and was not previously aligned, choose another one, aborting\n"), "red");
		return FALSE;
	}
	// check the number of dof if -interp=none
	transformation_type min, max;
	guess_transform_from_seq(regargs->seq, regargs->layer, &min, &max, TRUE);
	if (!driz && max > SHIFT_TRANSFORMATION && regargs->interpolation == OPENCV_NONE) {
		siril_log_color_message(_("Applying registration computed with higher degree of freedom (%d) than shift is not allowed when interpolation is set to none, aborting\n"), "red", ((int)max + 1) * 2);
		return FALSE;
	}

	// check the consistency of output images size if -interp=none
	if (!driz && regargs->interpolation == OPENCV_NONE && (regargs->output_scale != 1.f || regargs->framing == FRAMING_MAX || regargs->framing == FRAMING_MIN)) {
		siril_log_color_message(_("Applying registration with changes in output image size is not allowed when interpolation is set to none , aborting\n"), "red");
		return FALSE;
	}

	// check the consistency of images size if -interp=none
	if (!driz && regargs->interpolation == OPENCV_NONE && regargs->seq->is_variable) {
		siril_log_color_message(_("Applying registration on images with different sizes when interpolation is set to none is not allowed, aborting\n"), "red");
		return FALSE;
	}

	if (!driz && regargs->interpolation == OPENCV_NONE && regargs->undistort) {
		siril_log_color_message(_("Applying registration on images with distortions when interpolation is set to none is not allowed, aborting\n"), "red");
		return FALSE;
	}

	// check that we are not trying to apply identity transform to all the images
	if (max == IDENTITY_TRANSFORMATION) {
		siril_log_color_message(_("Existing registration data is a set of identity matrices, no transformation would be applied, aborting\n"), "red");
		return FALSE;
	}

	// check that we are not trying to apply null transform to all the images
	if (max == NULL_TRANSFORMATION || (regargs->seq->selnum <= 1) ) {
		siril_log_color_message(_("Existing registration data is a set of null matrices, no transformation would be applied, aborting\n"), "red");
		return FALSE;
	}

	// force -selected if some matrices were null
	if (min == NULL_TRANSFORMATION) {
		siril_log_color_message(_("Some images were not registered, excluding them\n"), "salmon");
		regargs->filters.filter_included = TRUE;
	}

	// cog framing method requires all images to be of same size
	if (regargs->framing == FRAMING_COG && regargs->seq->is_variable) {
		siril_log_color_message(_("Framing method \"cog\" requires all images to be of same size, aborting\n"), "red");
		return FALSE;
	}
	return TRUE;
}

