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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_math.h>


#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_world_cs.h"
#include "core/siril_log.h"
#include "algos/siril_wcs.h"
#include "drizzle/cdrizzlebox.h"
#include "drizzle/cdrizzleutil.h"
#include "io/sequence.h"
#include "io/siril_catalogues.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "registration/matching/degtorad.h"


#define DEBUG_ASTROREG

// computing center cog dealing with range jumps
// https://stackoverflow.com/questions/6671183/calculate-the-center-point-of-multiple-latitude-longitude-coordinate-pairs
// ra is lon, dec is lat
static void compute_center_cog(double *ra, double *dec, int n, gboolean *incl, double *RA, double *DEC) {
	double x = 0.;
	double y = 0.;
	double z = 0.;
	int nb = 0;
	for (int i = 0; i < n; i++) {
		if (!incl[i])
			continue;
		double ra_rad = ra[i] * DEGTORAD;
		double dec_rad = dec[i] * DEGTORAD;
		x += cos(dec_rad) * cos(ra_rad);
		y += cos(dec_rad) * sin(ra_rad);
		z += sin(dec_rad);
		nb++;
	}
	if (nb == 0)
		return;
	x /= (double)nb;
	y /= (double)nb;
	z /= (double)nb;
	double norm = sqrt( x * x + y * y + z * z);
	x /= norm;
	y /= norm;
	z /= norm;
	*RA = atan2(y, x) * RADTODEG;
	*DEC = atan2(z, sqrt(x * x + y * y)) * RADTODEG;
}


static gboolean get_scales_and_framing(struct wcsprm *WCSDATA, Homography *K, double *framing) {
	K->h00 = -1;
	K->h11 = -1;
	*framing = -1;
	double pc[2][2];
	double cdelt[2];
	double *pcij = WCSDATA->pc;
	if (!pcij)
		return FALSE;
	for (int i = 0; i < 2; i++) {
		cdelt[i] = WCSDATA->cdelt[i];
		for (int j = 0; j <2; j++) {
			pc[i][j] = *(pcij++);
		}
	}
	int flipped = (image_is_flipped_from_wcs(WCSDATA)) ? -1 : 1;

	K->h00 = -RADTODEG / cdelt[0];
	K->h11 = flipped * RADTODEG / cdelt[1];
	double rotation1 = atan2(flipped * pc[0][1], pc[0][0]);
	double rotation2 = atan2(-pc[1][0], flipped * pc[1][1]);
	*framing = -0.5 * (rotation1 + rotation2) * RADTODEG;
	return TRUE;
}

static int calcH_from_corners(int rx, int ry, int ry_ref, wcsprm_t *prm_img, wcsprm_t *prm_ref, Homography *H) {
	// corners in siril coords
	double corner_x[4] = { 0., (double)rx, (double)rx, 0. };
	double corner_y[4] = { 0., 0., (double)ry, (double)ry };
	// other inits
	double corner_ra[4] = { 0. };
	double corner_dec[4] = { 0. };
	double corner_ref_x[4] = { 0. };
	double corner_ref_y[4] = { 0. };
	// we compute the world coordinates of the corners of the image with its own wcs solution
	for (int i = 0; i < 4; i++)
		// uses x/y in siril coords
		pix2wcs2(prm_img, corner_x[i], corner_y[i], &corner_ra[i], &corner_dec[i]);

	// then we compute the corners coordinates from their world position in the reference framing
	for (int i = 0; i < 4; i++)
		// returns x/y in siril coords
		wcs2pix2(prm_ref, corner_ra[i], corner_dec[i], &corner_ref_x[i], &corner_ref_y[i]);

	// we convert all the corner coordinates from siril to opencv coordinates
	for (int i = 0; i < 4; i++) {
		corner_x[i] = corner_x[i] - 0.5;
		corner_y[i] = ry - corner_y[i] - 0.5;
		corner_ref_x[i] = corner_ref_x[i] - 0.5;
		corner_ref_y[i] = ry_ref - corner_ref_y[i] - 0.5;
	}
	int ret = cvCalcH_from_corners(corner_x, corner_y, corner_ref_x, corner_ref_y, H);
	if (ret < 0) {
		siril_log_color_message(_("Failed to compute homography from corners\n"), "red");
		return -1;
	}
	return 0;
}

typedef struct {
	double ra;
	double dec;
	sequence *seq;
	gboolean *incl;
	Homography *Kref;
	Homography *Ks;
	Homography *Rstmp;
} optim_params_t;

static void get_maxframing(double framing, void *params, int *width,  int *height, double *area) {
	optim_params_t *p = (optim_params_t *)params;
	int n = p->seq->number;
	// computing relative rotations wrt to proj center + central image framing
	Homography Rref = { 0 };
	double angles[3] = { 90. - p->ra,  90. - p->dec, framing};
	rotation_type rottypes[3] = { ROTZ, ROTX, ROTZ};
	cvRotMat3(angles, rottypes, TRUE, &Rref);
	int xmin = INT_MAX, xmax = -INT_MAX;
	int ymin = INT_MAX, ymax = -INT_MAX;
	for (int i = 0;  i < n; i++) {
		if (!p->incl[i])
			continue;
		Homography H = { 0 }, R = { 0 };
		cvRelRot(&Rref, p->Rstmp + i, &R);
		cvcalcH_fromKKR(p->Kref, p->Ks + i, &R, &H);
		framing_roi roi = { 0 };
		int rx = ((p->seq->is_variable) ? p->seq->imgparam[i].rx : p->seq->rx);
		int ry = ((p->seq->is_variable) ? p->seq->imgparam[i].ry : p->seq->ry);
		compute_roi(&H, rx, ry, &roi);
		xmin = min(xmin, roi.x);
		ymin = min(ymin, roi.y);
		xmax = max(xmax, roi.x + roi.w);
		ymax = max(ymax, roi.y + roi.h);
	}
	if (width)
		*width = (xmax - xmin);
	if (height)
		*height = (ymax - ymin);
	if (area)
		*area = (double)(xmax - xmin) * (ymax - ymin);
}

static double get_maxframing_area(double framing, void *params) {
	double area = 0.;
	get_maxframing(framing,params, NULL,  NULL, &area);
	return area;
}

// this function minimizes the area by searching the optimal framing
static int optimize_max_framing(double *framingref, optim_params_t *params) {
	int status;
	int iter = 0, max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m = (*framingref < 0) ? *framingref + 180. : *framingref;
	double a = 0.0, b = 180.;
	gsl_function F;
	// optim_params_t params = { ra0, dec0, seq, incl, Kref, Ks, Rstmp};
	F.function = &get_maxframing_area;
	F.params = params;
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, &F, m, a, b);

	siril_debug_print("%5s [%9s, %9s] %9s %10s %9s %10s\n",
		"iter", "lower", "upper", "min",
		"err", "err(est)", "val");

	do {
		iter++;
		status = gsl_min_fminimizer_iterate(s);
		m = gsl_min_fminimizer_x_minimum(s);
		a = gsl_min_fminimizer_x_lower(s);
		b = gsl_min_fminimizer_x_upper(s);
		double v = GSL_FN_EVAL(&F, m);
		status = gsl_min_test_interval(a, b, 0.001, 0.0);
		siril_debug_print("%5d [%.7f, %.7f] "
				"%.7f %.7f %.0f\n",
				iter, a, b,
				m, b - a, v);
	} while (status == GSL_CONTINUE && iter < max_iter);
	
	if (status == GSL_SUCCESS) {
		siril_debug_print("Converged\n");
		double optimval = gsl_min_fminimizer_x_minimum(s);
		// TODO: should we correct to make sure the image is always horizontal?
		// int width = 0, height = 0;
		// get_maxframing(optimval, params, &width, &height, NULL);
		// if (width < height) // we try to have an image roughly horizontal
		// 	optimval += 90.;
		if (*framingref < 0)
			*framingref = optimval - 180.;
		else
			*framingref= optimval;
	}
	gsl_min_fminimizer_free(s);
	return status;
}

int compute_Hs_from_astrometry(sequence *seq, int *included, int ref_index, struct wcsprm *WCSDATA, framing_type framing, int layer, Homography *Hout, struct wcsprm **prmout) {
	int retval = 0;
	double ra0 = 0., dec0 = 0.;
	int n = seq->number;
	double *RA = calloc(n, sizeof(double)); // calloc to prevent possibility of uninit variable highlighted by coverity
	double *DEC = calloc(n, sizeof(double)); // as above
	double *dist = calloc(n, sizeof(double));
	Homography *Ks = NULL, *Rstmp = NULL;
	Rstmp = calloc(n, sizeof(Homography)); // camera tmp rotation matrices
	Ks = calloc(n, sizeof(Homography)); // camera intrinsic matrices
	gboolean *incl = calloc(n, sizeof(gboolean));
	gboolean included_passed = included != NULL;
	if (!RA || !DEC || !dist || !WCSDATA || !Rstmp || !Ks || !incl) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto free_all;
	}
	// collect all the centers
	for (int i = 0; i < n; i++) {
		if (!included_passed && !seq->imgparam[i].incl) // if included array is not passed, we use incl from seq
			continue;
		else if (included_passed && !included[i]) // if included array is passed, we use it instead of seq
			continue;
		int width = (seq->is_variable) ? seq->imgparam[i].rx : seq->rx;
		int height = (seq->is_variable) ? seq->imgparam[i].ry : seq->ry;
		center2wcs2(WCSDATA + i, width, height, RA + i, DEC + i);
		incl[i] = TRUE;
	}

	// find sequence cog and closest image to use as reference
	compute_center_cog(RA, DEC, n, incl, &ra0, &dec0);
	ra0 = (ra0 < 0) ? ra0 + 360. : ra0;
	siril_log_message(_("Sequence COG - RA:%7.3f - DEC:%+7.3f\n"), ra0, dec0);

	int refindex = (ref_index < 0) ? seq->reference_image : ref_index;

	// Obtaining Camera extrinsic and instrinsic matrices (resp. R and K)
	// ##################################################################
	// We will define the rotations made by the camera wrt World:
	// - We start with the camera pointing towards North Pole, X (sensor width) aligned with RA = 0deg
	// - We first make a rotation around cam Z axis of 90 + RA degrees
	// - We then make a rotation around cam X axis of 90 - DEC degrees (0 to 180 from N to S pole)
	// - We finally rotate the camera around its Z axis to get the correct framing (called "framing" angle below)
	// We repeat this calc for each image
	// We can then compute the relative rotation matrix between the images and refimage
	// This is made by computing Rref.t()*Rimg, this is the extrinsic matrix for a pure rotation wrt. mosaic ref
	// The intrisic camera matrix is obtained thanks to the CDij or PCij+CDELTi WCS values:
	// The focal length is 180 / pi / average of abs(CDELTi) - Note: CDELTi * 3600 = image sampling in "/px
	// The "framing" angle is obtained by making the PC matrix a true rotation matrix (averaging of CROTAi)

	rotation_type rottypes[3] = { ROTZ, ROTX, ROTZ};
	double framingref = 0.;
	int anglecount = 0;

	for (int i = 0; i < n; i++) {
		if (!incl[i])
			continue;
		double angles[3];
		angles[0] = 90 - RA[i];
		angles[1] = 90. - DEC[i];
		// initializing K with center point (-0.5 is for opencv convention)
		cvGetEye(Ks + i); // initializing to unity
		Ks[i].h02 = 0.5 * ((seq->is_variable) ? seq->imgparam[i].rx : seq->rx) - 0.5;
		Ks[i].h12 = 0.5 * ((seq->is_variable) ? seq->imgparam[i].ry : seq->ry) - 0.5;
		if (!get_scales_and_framing(WCSDATA + i, Ks + i, angles + 2)) {
			siril_log_message(_("Could not compute camera parameters of image %d from sequence %s\n"),
			i + 1, seq->seqname);
			retval = 1;
			goto free_all;
		}
		cvRotMat3(angles, rottypes, TRUE, Rstmp + i);
		siril_debug_print("Image #%d - rot:%.3f\n", i + 1, angles[2]);
		if (framing == FRAMING_CURRENT && i == refindex) {
			framingref = angles[2];
			anglecount++;
		}
		if (framing != FRAMING_CURRENT) {
			framingref += (angles[2] < -90.) ? angles[2] + 180. : (angles[2] > 90.) ? angles[2] - 180 : angles[2];
			anglecount++;
		}
		siril_log_message(_("Image #%5d - RA:%7.3f - DEC:%+7.3f - Rotation:%+6.1f\n"), i + 1, RA[i], DEC[i], angles[2]);
	}
	if (!anglecount) {
		siril_log_color_message(_("No image selected after computing transformations, aborting\n"), "red");
		goto free_all;
	}
	framingref /= anglecount;

	if (framing == FRAMING_MAX) {
		// We compute the H matrices wrt to ref as Kref * Rrel * Kimg^-1
		// The K matrices have the focals on the diag and (rx/2, ry/2) for the translation terms
		// We use that to optimize the framing angle of the resulting mosaic
		// For the final mosaics, we will use instead computing H from the 4 corners of the image
		// and the reference image, as it is more accurate if there's a bit of skew in the solution.
		float fscale = 0.5 * (fabs(Ks[refindex].h00) + fabs(Ks[refindex].h11));
		Homography Kref = { 0 };
		cvGetEye(&Kref);
		Kref.h00 = fscale;
		Kref.h11 = fscale;
		Kref.h02 = Ks[refindex].h02;
		Kref.h12 = Ks[refindex].h12;
#ifdef DEBUG_ASTROREG
		siril_debug_print("Scale: %.3f\n", fscale);
		print_H(&Kref);
#endif
		optim_params_t params = { ra0, dec0, seq, incl, &Kref, Ks, Rstmp};
		int status = optimize_max_framing(&framingref, &params);
		if (!status)
			siril_log_message(_("Sequence optimal framing: %.1f\n"), framingref);
		else
			siril_log_message(_("Sequence framing: %.1f\n"), framingref);
	} else {
		siril_log_message(_("Sequence framing: %.1f\n"), framingref);
	}

	// computing relative rotations wrt to proj center + central image framing
	wcsprm_t *wcsref = NULL;
	if (framing == FRAMING_CURRENT) {
		wcsref = wcs_copy_linear(WCSDATA + refindex);
		if (image_is_flipped_from_wcs(wcsref)) {
			Homography H = { 0 };
			cvGetEye(&H);
			int ry = (seq->is_variable) ? seq->imgparam[refindex].ry : seq->ry;
			H.h11 = -1.;
			H.h12 = (double)ry;
			reframe_wcs(wcsref, &H);
		}
	} else {
		wcsref = calloc(1, sizeof(wcsprm_t));
		double scale = 0.5 * (fabs((WCSDATA + refindex)->cdelt[0])+fabs((WCSDATA + refindex)->cdelt[1]));
		create_wcs(ra0, dec0, scale, framingref, 
			(seq->is_variable) ? seq->imgparam[refindex].rx : seq->rx, 
			(seq->is_variable) ? seq->imgparam[refindex].ry : seq->ry, 
			wcsref);
	}

	if (!wcsref) {
		siril_log_color_message(_("Could not create WCS structure for framing, aborting\n"), "red");
		retval = 1;
		goto free_all;
	}

	regdata *current_regdata = seq->regparam[layer];
	int xmin = INT_MAX, xmax = -INT_MAX;
	int ymin = INT_MAX, ymax = -INT_MAX;
	int rx_ref = ((seq->is_variable) ? seq->imgparam[ref_index].rx : seq->rx);
	int ry_ref = ((seq->is_variable) ? seq->imgparam[ref_index].ry : seq->ry);
	for (int i = 0;  i < n; i++) {
		if (!incl[i])
			continue;
		Homography H = { 0 };
		int rx = ((seq->is_variable) ? seq->imgparam[i].rx : seq->rx);
		int ry = ((seq->is_variable) ? seq->imgparam[i].ry : seq->ry);
		wcsprm_t *wcscurr = wcs_copy_linear(WCSDATA + i);
		calcH_from_corners(rx, ry, ry_ref, wcscurr, wcsref, &H);
		wcsfree(wcscurr);
		framing_roi roi = { 0 };
		compute_roi(&H, rx, ry, &roi);
		current_regdata[i].H = H;
		siril_debug_print("Image %d\n", i + 1);
#ifdef DEBUG_ASTROREG
		print_H(&H);
		siril_debug_print("%d,%d,%d,%d,%d\n", i + 1, roi.x, roi.y, roi.w, roi.h);
#endif
		xmin = min(xmin, roi.x);
		ymin = min(ymin, roi.y);
		xmax = max(xmax, roi.x + roi.w);
		ymax = max(ymax, roi.y + roi.h);
	}
	siril_debug_print("Framing: %d,%d,%d,%d\n", xmin, ymin, xmax - xmin, ymax - ymin);
	seq->reference_image = refindex;

	if (prmout)
		*prmout = wcsref;
	else {
		wcsfree(wcsref);
		wcsref = NULL;
	}

free_all:
	free(RA);
	free(DEC);
	free(dist);
	free(incl);
	free(Ks);
	free(Rstmp);
	siril_log_message(_("Astrometric registration computed.\n"));
	if (Hout) {
		cvGetEye(Hout);
	}
	// we don't free WCSDATA as it will be further used to initialize distortion data
	return retval;
}

int collect_sequence_astrometry(struct registration_args *regargs, int *included) {
	int n = regargs->seq->number;
	int retval = 0;
	fits fit = { 0 };
	for (int i = 0; i < n; i++) {
		if (regargs->filtering_criterion && !regargs->filtering_criterion(regargs->seq, i, regargs->filtering_parameter))
			continue;
		if (seq_read_frame_metadata(regargs->seq, i, &fit)) {
			siril_log_message(_("Could not load image %d from sequence %s\n"),
			i + 1, regargs->seq->seqname);
			retval = 1;
			break;
		}
		if (!has_wcs(&fit)) {
			siril_log_message(_("Image %d has not been plate-solved, unselecting\n"), i + 1);
			regargs->seq->imgparam[i].incl = FALSE;
			regargs->filters.filter_included = TRUE;
			convert_parsed_filter_to_filter(&regargs->filters,
				regargs->seq, &regargs->filtering_criterion,
				&regargs->filtering_parameter);
			continue;
		}
		regargs->WCSDATA[i].flag = -1;
		wcssub(1, fit.keywords.wcslib, NULL, NULL, regargs->WCSDATA + i); // copying wcsprm structure for each fit to avoid reopening
		clearfits(&fit);
		included[i] = TRUE;
	}
	return retval;
}

