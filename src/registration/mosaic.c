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
#include "algos/siril_wcs.h"
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

int register_mosaic(struct registration_args *regargs) {
#ifndef HAVE_WCSLIB
	siril_log_color_message(_("Mosaic registration requires to have wcslib dependency installed, aborting\n"), "red");
	return 1;
#else

	int retval = 0;
	struct timeval t_start, t_end;

	regdata *current_regdata = star_align_get_current_regdata(regargs); // clean the structure if it exists, allocates otherwise
	if (!current_regdata) return -1;

	gettimeofday(&t_start, NULL);
	fits fit = { 0 };
	double ra0 = 0., dec0 = 0.;
	int n = regargs->seq->number;
	double *RA = malloc(n * sizeof(double));
	double *DEC = malloc(n * sizeof(double));
	double *dist = malloc(n * sizeof(double));
	struct wcsprm *WCSDATA = malloc(n * sizeof(struct wcsprm));
	double *X = malloc(4 * n * sizeof(double));
	double *Y = malloc(4 * n * sizeof(double));
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
		RA[i] = fit.wcsdata.ra;
		DEC[i] = fit.wcsdata.dec;
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

	// Now computing each image corners ra/dec then calculating homography on refimage
	for (int i = 0; i < n; i++) {
		double rx = (regargs->seq->is_variable) ? regargs->seq->imgparam[i].rx : regargs->seq->rx;
		double ry = (regargs->seq->is_variable) ? regargs->seq->imgparam[i].ry : regargs->seq->ry;
		// double x[4] = {0., rx, rx, 0.};
		// double y[4] = {0., 0., ry, ry};
		double x[4] = {0., rx, rx, 0.};
		double y[4] = {0., 0., ry, ry};
		double phi, theta;
		for (int j = 0; j < 4; j++) {
			int status, stat[NWCSFIX];
			double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];
			pixcrd[0] = x[j];
			pixcrd[1] = y[j];
			// converting corners to ra/dec
			status = wcsp2s(WCSDATA + i, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
			if (status) {
				siril_log_color_message(_("Mosaic: Error when projecting point #%d of image #%d, aborting\n"), "red", j, i + 1);
				retval = 1;
				goto free_all;
			}
			// converting ra/dec to pixels in refimage projection
			status = wcss2p(WCSDATA + refindex, 1, 0, world, &phi, &theta, imgcrd, pixcrd, stat);
			if (status) {
				siril_log_color_message(_("Mosaic: Error when re-projecting point #%d of image #%d, aborting\n"), "red", j, i + 1);
				retval = 1;
				goto free_all;
			}
			// storing pixels coords to compute H matrices
			X[i * 4 + j] = pixcrd[0];
			Y[i * 4 + j] = ry - pixcrd[1];
			SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(world[0], world[1]);
			gchar *rap = siril_world_cs_alpha_format(world_cs, "%02dh%02dm%02ds");
			gchar *decp = siril_world_cs_delta_format(world_cs, "%c%02dd%02d\'%02d\"");
			siril_debug_print("x0=%8.1f px, y0=%8.1f px, x1=%8.1f px, y1=%8.1f px (%s , %s)\n", x[j], y[j], X[i * 4 + j], Y[i * 4 + j], rap, decp);

			// siril_world_cs_unref(world_cs);
			// //ra2hms(double ra, int *h, int *m, double *s)
			// //dec2dms(double dec, int *sign, int *d, int *m, double *s)
			// int rah, ram, decsign, decd, decm;
			// double ras, decs;
			// ra2hms(world[0], &rah, &ram, &ras);
			// dec2dms(world[1], &decsign, &decd, &decm, &decs);
			// printf("%8.1f, %8.1f,", X[i * 4 + j], Y[i * 4 + j]);
			// printf("RA: %d.%d.%.0f DEC:%s%d.%d.%.0f\n", rah, ram, ras, (decsign==1 ? "+":"-"), decd, decm, decs);
		}
	}

	for (int i = 0; i < n; i++) {
		Homography H = { 0 };
		if (i != refindex) {
			if (cvCalculH_exact(X, Y, refindex, i, &H)) {
				failed++;
				regargs->seq->imgparam[i].incl = FALSE;
				continue;
			}
		} else {
			cvGetEye(&H);
		}
		nb_aligned++;
		current_regdata[i].roundness = 1.;
		current_regdata[i].fwhm = 0.;
		current_regdata[i].weighted_fwhm = 0.;
		current_regdata[i].background_lvl = 0.;
		current_regdata[i].number_of_stars = 4;
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
// 	for (int i = 0; i < regargs->seq->number; i++)
// 		free_fitted_stars(mosaic_args->stars[i]);
// 	free(mosaic_args->stars);
// 	free(mosaic_args->nb_stars);
// 	free(mosaic_args);
// 	free(fwhm);
// 	free(roundness);
// 	free(B);
// 	free(A);
// 	free(Acut);
// 	free(included);
// 	free(tmp_included);
// 	free(meaningful);
// 	free(scores);
	free(RA);
	free(DEC);
  	free(dist);
	free(X);
	free(Y);
	// TODO: properly free WCS array
	return retval;
#endif
}

