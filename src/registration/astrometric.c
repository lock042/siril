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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_world_cs.h"
#include "core/siril_log.h"
#include "algos/astrometry_solver.h"
#include "algos/siril_wcs.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "drizzle/cdrizzlebox.h"
#include "drizzle/cdrizzlemap.h"
#include "drizzle/cdrizzleutil.h"
#include "io/sequence.h"
#include "io/siril_catalogues.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "registration/matching/degtorad.h"


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

int compute_Hs_from_astrometry(sequence *seq, struct wcsprm *WCSDATA, framing_type framing, int layer, Homography *Hout, struct wcsprm **prmout) {
	int retval = 0;
	double ra0 = 0., dec0 = 0.;
	int n = seq->number;
	double *RA = calloc(n, sizeof(double)); // calloc to prevent possibility of uninit variable highlighted by coverity
	double *DEC = calloc(n, sizeof(double)); // as above
	double *dist = calloc(n, sizeof(double));
	Homography *Rs = NULL, *Ks = NULL;
	Rs = calloc(n, sizeof(Homography)); // camera rotation matrices
	Ks = calloc(n, sizeof(Homography)); // camera intrinsic matrices
	gboolean *incl = calloc(n, sizeof(gboolean));
	if (!RA || !DEC || !dist || !WCSDATA || !Rs || !Ks || !incl) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto free_all;
	}
	// collect all the centers
	for (int i = 0; i < n; i++) {
		if (!seq->imgparam[i].incl)
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
	int refindex = seq->reference_image;

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
		cvRotMat3(angles, rottypes, TRUE, Rs + i);
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
		siril_log_color_message(_("No image selected after computing transformations, aborting%s\n"), "red");
		goto free_all;
	}
	framingref /= anglecount;
	siril_log_message("Sequence framing: %.3f\n", framingref);

	// computing relative rotations wrt to proj center + central image framing
	Homography Rref = { 0 };
	double angles[3] = { 0 };
	switch (framing) {
		case FRAMING_CURRENT:
		default:
			angles[0] = 90. - RA[seq->reference_image];
			angles[1] = 90. - DEC[seq->reference_image];
			break;
		case FRAMING_COG:
		case FRAMING_MAX:
		case FRAMING_MIN:
			angles[0] = 90. - ra0;
			angles[1] = 90. - dec0;
			break;
	}
	angles[2] = framingref;
	cvRotMat3(angles, rottypes, TRUE, &Rref);
	for (int i = 0; i < n; i++) {
		if (!incl[i])
			continue;
		cvRelRot(&Rref, Rs + i);
	}
	regdata *current_regdata = seq->regparam[layer];

	// We compute the H matrices wrt to ref as Kref * Rrel * Kimg^-1
	// The K matrices have the focals on the diag and (rx/2, ry/2) for the translation terms
	// We add to the seq in order to display the alignment
	float fscale = 0.5 * (fabs(Ks[refindex].h00) + fabs(Ks[refindex].h11));
	Homography Kref = { 0 };
	cvGetEye(&Kref);
	Kref.h00 = fscale;
	Kref.h11 = fscale;
	Kref.h02 = Ks[refindex].h02;
	Kref.h12 = Ks[refindex].h12;
	siril_debug_print("Scale: %.3f\n", fscale);
	print_H(&Kref);

	// We compute the H matrices wrt to ref as Kref * Rrel * Kimg^-1
	// The K matrices have the focals on the diag and (rx/2, ry/2) for the translation terms
	// We add to the seq in order to display a first alignment
	for (int i = 0;  i < n; i++) {
		if (!incl[i])
			continue;
		Homography H = { 0 };
		cvcalcH_fromKKR(Kref, Ks[i], Rs[i], &H);
		framing_roi roi = { 0 };
		int rx = ((seq->is_variable) ? seq->imgparam[i].rx : seq->rx);
		int ry = ((seq->is_variable) ? seq->imgparam[i].ry : seq->ry);
		compute_roi(&H, rx, ry, &roi);
		current_regdata[i].H = H;
		siril_debug_print("Image %d\n", i + 1);
		print_H(Ks + i);
		print_H(Rs + i);
		print_H(&H);
		siril_debug_print("%d,%d,%d,%d,%d\n", i + 1, roi.x, roi.y, roi.w, roi.h);
	}
	seq->reference_image = refindex;

	if (prmout) {
		if (framing != FRAMING_CURRENT) {
			*prmout = calloc(1, sizeof(wcsprm_t));
			double scale = 0.5 * (fabs((WCSDATA + refindex)->cdelt[0])+fabs((WCSDATA + refindex)->cdelt[1]));
			create_wcs(ra0, dec0, scale, framingref, (seq->is_variable) ? seq->imgparam[refindex].rx : seq->rx, (seq->is_variable) ? seq->imgparam[refindex].ry : seq->ry, *prmout);
		} else {
			*prmout = wcs_deepcopy(WCSDATA + refindex, NULL);
		}
	}

free_all:
	free(RA);
	free(DEC);
	free(dist);
	free(incl);
	free(Ks);
	free(Rs);
	siril_log_message(_("Astrometric registration computed.\n"));
	if (Hout) {
		cvGetEye(Hout);
	}
	// we don't free WCSDATA as it will be further used to initialize distortion data
	return retval;
}

int collect_sequence_astrometry(struct registration_args *regargs) {
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
	}
	return retval;
}

