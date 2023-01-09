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

#define HAVE_WCSLIB 1

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
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"

#ifdef HAVE_WCSLIB
#include <wcslib.h>
#include <wcsfix.h>
#endif


#include "registration/registration.h"
#include "registration/matching/degtorad.h"

// computing center cog dealing with range jumps
// https://stackoverflow.com/questions/6671183/calculate-the-center-point-of-multiple-latitude-longitude-coordinate-pairs
// ra is lon, dec is lat
static void compute_center_cog(double *ra, double *dec, int n, double *RA, double *DEC) {
	double x = 0.;
	double y = 0.;
	double z = 0.;
	for (int i = 0; i < n; i++) {
		double ra_rad = ra[i] * DEGTORAD;
		double dec_rad = dec[i] * DEGTORAD;
		x += cos(dec_rad) * cos(ra_rad);
		y += cos(dec_rad) * sin(ra_rad);
		z += sin(dec_rad);
	}
	x /= (double)n;
	y /= (double)n;
	z /= (double)n;
	*RA = atan2(y, x) * RADTODEG;
	*DEC = atan2(z, sqrt(x * x + y * y)) * RADTODEG;
}

// Haversine formula on unit sphere
// https://en.wikipedia.org/wiki/Haversine_formula
// dec is phi, ra is lambda
static double compute_haversine_distance(double ra1, double dec1, double ra2, double dec2) {
	double dra_2 = 0.5 * (ra2 - ra1) * DEGTORAD;
	double ddec_2 = 0.5 * (dec2 - dec1) * DEGTORAD;
	double h = pow(sin(ddec_2), 2.) + cos(dec1 * DEGTORAD) * cos(dec2 * DEGTORAD) * pow(sin(dra_2), 2.);
	if (h > 1.) h = 1.;
	return 2 * asin(pow(h, 0.5));
}

static gboolean get_focs_and_framing(struct wcsprm *WCSDATA, double *focx, double *focy, double *framing) {
#ifndef HAVE_WCSLIB
	siril_log_color_message(_("Mosaic registration requires to have wcslib dependency installed, aborting\n"), "red");
	*focx = -1;
	*focy = -1;
	*framing = -1;
	return FALSE;
#endif
	*focx = -1;
	*focy = -1;
	*framing = -1;
	double cd[2][2], pc[2][2];
	double cdelt[2];
	double *pcij = WCSDATA->pc;
	double *cdij = WCSDATA->cd;
	for (int i = 0; i < 2; i++) {
		cdelt[i] = WCSDATA->cdelt[i];
		for (int j = 0; j <2; j++) {
			pc[i][j] = *(pcij++);
			cd[i][j] = *(cdij++);
		}
	}
	if (cdelt[0] == 1. && cdelt[1] == 1.)// CD formalism
		wcs_cd_to_pc(cd, pc, cdelt);

	*focx = RADTODEG / cdelt[0];
	*focy = RADTODEG / cdelt[1];
	double rotation1 = atan2(pc[0][1], pc[0][0]);
	double rotation2 = atan2(-pc[1][0], pc[1][1]);
	*framing = 0.5 * (rotation1 + rotation2) * RADTODEG;
	return TRUE;
}

int register_mosaic(struct registration_args *regargs) {
#ifndef HAVE_WCSLIB
	siril_log_color_message(_("Mosaic registration requires to have wcslib dependency installed, aborting\n"), "red");
	return 1;
#else

	int retval = 0;
	struct timeval t_start, t_end;

	regdata *current_regdata = star_align_get_current_regdata(regargs); // clean the structure if it exists, allocates otherwise
	if (!current_regdata) return -1;
	regargs->type = HOMOGRAPHY_TRANSFORMATION; // forcing homography calc

	gettimeofday(&t_start, NULL);
	fits fit = { 0 };
	double ra0 = 0., dec0 = 0.;
	int n = regargs->seq->number;
	double *RA = malloc(n * sizeof(double));
	double *DEC = malloc(n * sizeof(double));
	double *dist = malloc(n * sizeof(double));
	struct wcsprm *WCSDATA = malloc(n * sizeof(struct wcsprm));
	int failed = 0, nb_aligned = 0;

	if (!RA || !DEC || !dist || !WCSDATA) {
		PRINT_ALLOC_ERR;
		retval = 1;
		goto free_all;
	}

	for (int i = 0; i < n; i++) {
		if (seq_read_frame_metadata(regargs->seq, i, &fit)) {
			siril_log_message(_("Could not load image %d from sequence %s\n"),
			i + 1, regargs->seq->seqname);
			retval = 1;
			goto free_all;
		}
		if (!has_wcs(&fit)) {
			siril_log_message(_("Image %d has not been plate-solved, aborting\n"), i + 1);
			retval = 1;
			goto free_all;
		}
		center2wcs(&fit, RA + i, DEC + i);
		WCSDATA[i].flag = -1;
		wcssub(1, fit.wcslib, NULL, NULL, WCSDATA + i); // copying wcsprm structure for each fit to avoid reopening
		clearfits(&fit);
		siril_debug_print("Image #%d - RA:%.3f - DEC:%.3f\n", i, RA[i], DEC[i]);
	}

	// find sequence cog and closes image to use as projection point
	compute_center_cog(RA, DEC, n, &ra0, &dec0);
	siril_debug_print("Sequence COG - RA:%.3f - DEC:%.3f\n", ra0, dec0);
	int refindex = -1;
	double mindist = DBL_MAX;
	for (int i = 0; i < n; i++) {
		dist[i] = compute_haversine_distance(RA[i], DEC[i], ra0, dec0);
		if (dist[i] < mindist) {
			mindist = dist[i];
			refindex = i;
		}
	}
	siril_debug_print("Closest image  is #%d - RA:%.3f - DEC:%.3f\n", refindex + 1, RA[refindex], DEC[refindex]);

	// TODO: should check allocation success
	Homography *Rs = calloc(n, sizeof(Homography)); // camera rotation matrices
	double *focx = calloc(n, sizeof(double)); // camera focal length on x
	double *focy = calloc(n, sizeof(double)); // camera focal length on y

	// Obtaining Camera extrinsic and instrinsic matrices (resp. R and K)
	// ##################################################################
	// We will define the rotations made by the camera wrt World:
	// - We start with the camera pointing towards North Pole, X (sensor width) aligned with RA = 0deg
	// - We first make a rotation around cam Z axis of 90 + RA degrees
	// - We then make a rotation around cam X axis of 90 - DEC degrees (0 to 180 from N to S pole)
	// - We finally rotate the camera around its Z axis to get the correct framing (called "framing" angle below)
	// We repeat this calc for each image and for the cog of the sequence
	// We can then compute the relative rotation matrix between the image and the cog
	// This is made by computing Rref.t()*Rimg, this is the extrinsic matrix for a pure rotation wrt. mosaic center
	// The intrisic camera matrix is obtained thanks to the CDij or PCij+CDELTi WCS values:
	// The focal length is 180 / pi / average of abs(CDELTi) - Note: CDELTi * 3600 = image sampling in "/px
	// The "framing" angle is obtained by making the PC matrix a true rotation matrix (averaging of CROTAi)

	rotation_type rottypes[3] = { ROTZ, ROTX, ROTZ};


	for (int i = 0; i < n; i++) {
		double angles[3];
		angles[0] = 90. + RA[i]; 
		angles[1] = 90. - DEC[i];
		if (!get_focs_and_framing(WCSDATA + i, focx + i, focy + i, angles + 2)) {
			siril_log_message(_("Could not compute camera parameters of image %d from sequence %s\n"),
			i + 1, regargs->seq->seqname);
			retval = 1;
			goto free_all;
		}
		cvRotMat3(angles, rottypes, TRUE, Rs + i);
	}

	// computing relative rotations wrt to ref image
	Homography Rref = { 0 };
	memcpy(&Rref, Rs + refindex, sizeof(Homography));
	for (int i = 0; i < n; i++) {
		cvRelRot(&Rref, Rs + i);
	}

	// now we store the relative rotations matrices and focals
	//writing outputs to *.smf file, a.k.a Siril Mosaic File
	char *filename;
	FILE *mscfile;
	if (!regargs->seq->seqname || regargs->seq->seqname[0] == '\0') {
		retval = 1;
		goto free_all;
	}
	filename = malloc(strlen(regargs->seq->seqname)+5);
	sprintf(filename, "%s.smf", regargs->seq->seqname);
	mscfile = g_fopen(filename, "w+t");
	for (int i = 0; i < n; i++) {
		fprintf(mscfile, "%2d %12.1f %12.1f R %g %g %g %g %g %g %g %g %g\n",
		i + 1,
		focx[i],
		focy[i],
		Rs[i].h00,
		Rs[i].h01,
		Rs[i].h02,
		Rs[i].h10,
		Rs[i].h11,
		Rs[i].h12,
		Rs[i].h20,
		Rs[i].h21,
		Rs[i].h22
		);
	}
	fclose(mscfile);

	// and adding to the seq in order to display some kind of alignment
	for (int i = 0; i < n; i++) {
		Homography H = { 0 };
		if (i != refindex) {
			cvcalcH_fromKR(Rs[i], focx, focy, refindex, i, &H);
		} else {
			cvGetEye(&H);
		}
		nb_aligned++;
		current_regdata[i].roundness = 1.;
		current_regdata[i].fwhm = 0.;
		current_regdata[i].weighted_fwhm = 0.;
		current_regdata[i].background_lvl = 0.;
		current_regdata[i].number_of_stars = 1;
		current_regdata[i].H = H;
		regargs->seq->imgparam[i].incl = TRUE;
	}

	// images may have been excluded but selnum wasn't updated
	regargs->seq->reference_image = refindex;
	fix_selnum(regargs->seq, FALSE);
	siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, nb_aligned);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);


free_all:

	free(RA);
	free(DEC);
	free(dist);
	// TODO: properly free WCS array
	// TODO: properly free Hs array
	// TODO: properly free Rs array
	return retval;
#endif
}

