/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <gio/gwin32inputstream.h>
#else
#include <gio/gunixinputstream.h>
#endif

#include "astrometry_solver.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "io/annotation_catalogues.h"
#include "algos/photometry.h"
#include "algos/photometric_cc.h"
#include "algos/spcc.h"
#include "algos/siril_wcs.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/siril_catalogues.h"
#include "io/local_catalogues.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "registration/matching/match.h"
#include "registration/matching/apply_match.h"
#include "registration/matching/atpmatch.h"
#include "gui/message_dialog.h"


#define DOWNSAMPLE_FACTOR 0.25
#define CONV_TOLERANCE 1E-3 // convergence tolerance in arcsec from the projection center
#define NB_GRID_POINTS 6 // the number of points in one direction to crete the X,Y meshgrid for inverse polynomial fiiting

#define CHECK_FOR_CANCELLATION_RET if (!get_thread_run()) { args->message = g_strdup(_("Cancelled")); args->ret = 1; goto clearup;}
#define CHECK_FOR_CANCELLATION if (!get_thread_run()) { ret = SOLVE_CANCELLED; goto clearup; }

#undef DEBUG		/* get some of diagnostic output */
#define ASTROMETRY_DEBUG 0

static gchar *asnet_version = NULL;

typedef struct {
	struct wcsprm *wcslib;
} solve_results;

static void debug_print_catalog_files(s_star *star_list_A, s_star *star_list_B) {
#if ASTROMETRY_DEBUG
	GFile *file = g_file_new_for_path("ABtars.csv");
	g_autoptr(GError) error = NULL;
	GOutputStream* output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE, G_FILE_CREATE_NONE, NULL, &error);
	if (error) {
		siril_debug_print("%s\n", error->message);
		return;
	}
	const gchar *header;
	header = "xA,yA,magA,indA,xB,yB,magB,indB\n";
	g_output_stream_write_all(output_stream, header, strlen(header), NULL, NULL,NULL);
	static gchar buffer[256] = { 0 };
	s_star *sA, *sB;
	for (sA = star_list_A, sB = star_list_B; ; sA = sA->next, sB = sB->next) {
		if (!sA) break;
		if (!sB) break;
		g_sprintf(buffer, "%.8f,%.8f,%g,%d,%.8f,%.8f,%g,%d\n", sA->x, sA->y, sA->mag, sA->id, sB->x, sB->y, sB->mag, sB->id);
		g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, NULL);
		memset(buffer, 0, 256);
	}
	g_object_unref(output_stream);
#endif
}

static struct astrometry_data *copy_astrometry_args(struct astrometry_data *args) {
	struct astrometry_data *ret = malloc(sizeof(struct astrometry_data));
	if (!ret) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	memcpy(ret, args, sizeof(struct astrometry_data));
	if (args->cat_center)
		ret->cat_center = siril_world_cs_copy(args->cat_center);
	if (args->ref_stars) {
		ret->ref_stars = calloc(1, sizeof(siril_catalogue));
		siril_catalogue_copy(args->ref_stars, ret->ref_stars);
	}
	ret->fit = NULL;
	ret->filename = NULL;
	return ret;
}

static void fov_in_DHMS(double var, gchar *fov) {
	int deg, decM;
	double decS;

	if (var < 0) {
		fprintf(stdout, "fov_in_DHMS: negative value, should not happen\n");
		return;
	}
	deg = (int) var;
	decM = abs((int) ((var - deg) * 60));
	decS = (fabs((var - deg) * 60) - decM) * 60;
	if (deg > 0)
		g_snprintf(fov, 256, "%02dd %02dm %.2lfs", deg, decM, decS);
	else if (decM > 0)
		g_snprintf(fov, 256, "%02d\' %.2lf\"", decM, decS);
	else if (decS > 0.0)
		g_snprintf(fov, 256, "%.2lf\"", decS);
}

/* get resolution in arcsec per pixel */
double get_resolution(double focal, double pixel) {
	if (focal <= 0.0 || pixel <= 0.0)
		return 0.0;
	return RADCONV / focal * pixel;
}

/* get diagonal field of view in arcmin, resolution in arcsec/px */
double get_fov_arcmin(double resolution, int rx, int ry) {
	uint64_t sqr_radius = (uint64_t) rx * (uint64_t) rx + (uint64_t) ry * (uint64_t) ry;
	double radius = resolution * sqrt((double)sqr_radius);	// in arcsec
	return radius / 60.0;	// in arcminutes
}

/* get half field of view in degrees, or angle from image centre, resolution in arcsec/px */
double get_radius_deg(double resolution, int rx, int ry) {
	uint64_t sqr_radius = ((uint64_t) rx * (uint64_t) rx + (uint64_t) ry * (uint64_t) ry) / 4;
	double radius = resolution * sqrt((double)sqr_radius);	// in arcsec
	return radius / 3600.0;	// in degrees
}

double compute_mag_limit_from_fov(double fov_degrees) {
	// Empiric formula for 1000 stars at 20 deg of galactic latitude
	double autoLimitMagnitudeFactor = 14.5;
	double m = autoLimitMagnitudeFactor * pow(fov_degrees, -0.179);
	// for astrometry, it can be useful to go down to mag 20, for
	// photometry the catalog's limit is 17 for APASS and 18 for NOMAD
	return round(100.0 * min(20.0, max(7.0, m))) / 100.;
}

static void compute_limit_mag(struct astrometry_data *args) {
	if (args->mag_mode == LIMIT_MAG_ABSOLUTE)
		args->ref_stars->limitmag = args->magnitude_arg;
	else {
		args->ref_stars->limitmag = compute_mag_limit_from_fov(args->used_fov / 60.0);
		if (args->mag_mode == LIMIT_MAG_AUTO_WITH_OFFSET)
			args->ref_stars->limitmag += args->magnitude_arg;
	}
	siril_debug_print("using limit magnitude %f\n", args->ref_stars->limitmag);
}

static void solve_is_blind(struct astrometry_data *args) {
	args-> blind = args->solver == SOLVER_LOCALASNET || (args->ref_stars && args->ref_stars->cat_index == CAT_LOCAL);
}

SirilWorldCS *get_eqs_from_header(fits *fit) {
	if (fit->wcsdata.objctra[0] != '\0' && fit->wcsdata.objctdec[0] != '\0')
		return siril_world_cs_new_from_objct_ra_dec(fit->wcsdata.objctra, fit->wcsdata.objctdec);

	else if (has_wcs(fit) && (fit->wcslib->crval[0] != 0.0 || fit->wcslib->crval[1] != 0.0))
		return siril_world_cs_new_from_a_d(fit->wcslib->crval[0], fit->wcslib->crval[1]);

	else if (fit->wcsdata.ra != 0.0 || fit->wcsdata.dec != 0.0)
		return siril_world_cs_new_from_a_d(fit->wcsdata.ra, fit->wcsdata.dec);
	return NULL;
}

static void update_wcsdata_from_wcs(struct astrometry_data *args) {
	if (has_wcsdata(args->fit))
		reset_wcsdata(args->fit);
	double ra0, dec0;
	if (args->solver == SOLVER_LOCALASNET) {
		center2wcs(args->fit, &ra0, &dec0); // asnet sometimes does not return solution at center despite the center flag
		snprintf(args->fit->wcsdata.pltsolvd_comment, FLEN_COMMENT, "Solved by Astrometry.net (%s)", asnet_version);
	} else {
		ra0 = args->fit->wcslib->crval[0];
		dec0 = args->fit->wcslib->crval[1];
		g_snprintf(args->fit->wcsdata.pltsolvd_comment, FLEN_COMMENT, "Siril internal solver");
	}
	args->fit->wcsdata.ra = ra0;
	args->fit->wcsdata.dec = dec0;
	args->fit->wcsdata.pltsolvd = TRUE;
	gchar *ra = siril_world_cs_alpha_format_from_double(ra0, "%02d %02d %.3lf");
	gchar *dec = siril_world_cs_delta_format_from_double(dec0, "%c%02d %02d %.3lf");
	g_sprintf(args->fit->wcsdata.objctra, "%s", ra);
	g_sprintf(args->fit->wcsdata.objctdec, "%s", dec);
	g_free(ra);
	g_free(dec);
	
}

static gboolean check_affine_TRANS_sanity(TRANS *trans) {
	double var1 = 0., var2 = 0.;
	var1 = fabs(trans->x10) - fabs(trans->y01);
	var2 = fabs(trans->y10) - fabs(trans->x01);
	siril_debug_print("abs(diff_cos)=%f et abs(diff_sin)=%f\n", var1, var2);
	return (fabs(var1) < 0.01 && fabs(var2) < 0.01);
}

static double get_det_from_trans(TRANS *trans) {
	return (trans->x10 * trans->y01 - trans->y10 * trans->x01);
}

static void get_cd_from_trans(TRANS *trans, double cd[][2]) {
	cd[0][0] = trans->x10;
	cd[0][1] = trans->x01;
	cd[1][0] = trans->y10;
	cd[1][1] = trans->y01;
}

static double get_center_offset_from_trans(TRANS *trans) {
	// this measures the remaining offset, expressed in pseudo-arcsec
	return sqrt(trans->x00 * trans->x00 + trans->y00 * trans->y00);
}
// 2x2 matrix vector multiplication
static void Mv(double matrix[2][2], double vector[2], double result[2]) {
	for (int i = 0; i < 2; i++) {
		result[i] = 0;
		for (int j = 0; j < 2; j++) {
			result[i] += matrix[i][j] * vector[j];
		}
	}
}
// and a helper function to simplify calling the M*v multiplication and decompose the result to its 2 values
static void Mvdecomp(double matrix[2][2], double v0, double v1, double *r0, double *r1) {
	double vector[2] = { v0, v1 };
	double result[2];
	Mv(matrix, vector, result);
	*r0 = result[0];
	*r1 = result[1];
}

static int add_disto_to_wcslib(struct wcsprm *wcslib, TRANS *trans, int rx, int ry) {
	if (!wcslib || !trans)
		return 1;

	// we start by zero-ing the constant terms because we have set CRVAL at the center of the solution
	trans->x00 = 0.;
	trans->y00 = 0.;

	// Definitions from https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf

	// We will create the inverse transform:
	// This is done by generating a grid of pixels coords (u,v)
	// and using the forward transform coeffs, generating the intermediate world coordinates grid (x,y)
	int nbpoints = NB_GRID_POINTS * NB_GRID_POINTS;
	struct s_star *uvgrid = create_grid_list(rx, ry, NB_GRID_POINTS);
	struct s_star *xygrid = create_grid_list(rx, ry, NB_GRID_POINTS); // these coords are then converted using the pixel-to-sky transform
	atApplyTrans(nbpoints, xygrid, trans);
	// we then extract the CD matrix and invert it
	double cd[2][2] = {{ 0. }}, cd_inv[2][2] = {{ 0. }};
	get_cd_from_trans(trans, cd);
	double det = get_det_from_trans(trans);
	double invdet = 1. / det;
	// CD^-1 is simply 1/det(CD)*[[d -b][-c a]] if CD is [[a b][c d]]
	cd_inv[0][0] =  invdet * cd[1][1];
	cd_inv[0][1] = -invdet * cd[0][1];
	cd_inv[1][0] = -invdet * cd[1][0];
	cd_inv[1][1] =  invdet * cd[0][0];

	// we form a linear trans structure which sends the iwc to corrected pixels coordinates U,V using CD^-1
	// see the definition in eq (4)
	TRANS transUV = { 0 };
	transUV.order = 1;
	transUV.x10 = cd_inv[0][0];
	transUV.x01 = cd_inv[0][1];
	transUV.y10 = cd_inv[1][0];
	transUV.y01 = cd_inv[1][1];
	atApplyTrans(nbpoints, xygrid, &transUV);
	// xygrid now holds the UV grid

	// we can then find the polynomials that send U,V to u,v the original pixel coordinates of the grid
	TRANS revtrans = { 0 };
	revtrans.order = trans->order;
	int status = atRecalcTrans(nbpoints, xygrid, nbpoints, uvgrid, AT_MATCH_MAXITER, AT_MATCH_HALTSIGMA, &revtrans);
	if (status) {
		siril_log_color_message(_("Could not invert the SIP distortion coefficients, try using a lower order or a linear solution\n"), "red");
		free_stars(&uvgrid);
		free_stars(&xygrid);
		return 1;
	}

	// We will apply CD^-1 to each pair in the forward trans structure to obtain the Aij/Bij SIP coeffs (using defs of eq (1), (2) and (3))
	// The inverse APij/BPij coeffs are simply the values of revtrans
	// For the inverse _10 and _01 terms, we need to substract 1., see definitions in eq (5) and (6))
	double A[5][5] = {{ 0. }}, B[5][5] = {{ 0. }}, AP[5][5]  = {{ 0. }}, BP[5][5]  = {{ 0. }}; // we deal with images up to order 4
	int N = trans->order;

	Mvdecomp(cd_inv, trans->x20, trans->y20, &A[2][0], &B[2][0]);
	Mvdecomp(cd_inv, trans->x11, trans->y11, &A[1][1], &B[1][1]);
	Mvdecomp(cd_inv, trans->x02, trans->y02, &A[0][2], &B[0][2]);

	AP[0][0] = revtrans.x00;
	AP[1][0] = revtrans.x10 - 1.;
	AP[0][1] = revtrans.x01;
	AP[2][0] = revtrans.x20;
	AP[1][1] = revtrans.x11;
	AP[0][2] = revtrans.x02;

	BP[0][0] = revtrans.y00;
	BP[1][0] = revtrans.y10;
	BP[0][1] = revtrans.y01 - 1.;
	BP[2][0] = revtrans.y20;
	BP[1][1] = revtrans.y11;
	BP[0][2] = revtrans.y02;

	if (trans->order >= AT_TRANS_CUBIC) {
		Mvdecomp(cd_inv, trans->x30, trans->y30, &A[3][0], &B[3][0]);
		Mvdecomp(cd_inv, trans->x21, trans->y21, &A[2][1], &B[2][1]);
		Mvdecomp(cd_inv, trans->x12, trans->y12, &A[1][2], &B[1][2]);
		Mvdecomp(cd_inv, trans->x03, trans->y03, &A[0][3], &B[0][3]);

		AP[3][0] = revtrans.x30;
		AP[2][1] = revtrans.x21;
		AP[1][2] = revtrans.x12;
		AP[0][3] = revtrans.x03;

		BP[3][0] = revtrans.y30;
		BP[2][1] = revtrans.y21;
		BP[1][2] = revtrans.y12;
		BP[0][3] = revtrans.y03;
	}

	if (trans->order >= AT_TRANS_QUARTIC) {
		Mvdecomp(cd_inv, trans->x40, trans->y40, &A[4][0], &B[4][0]);
		Mvdecomp(cd_inv, trans->x31, trans->y31, &A[3][1], &B[3][1]);
		Mvdecomp(cd_inv, trans->x22, trans->y22, &A[2][2], &B[2][2]);
		Mvdecomp(cd_inv, trans->x13, trans->y13, &A[1][3], &B[1][3]);
		Mvdecomp(cd_inv, trans->x04, trans->y04, &A[0][4], &B[0][4]);

		AP[4][0] = revtrans.x40;
		AP[3][1] = revtrans.x31;
		AP[2][2] = revtrans.x22;
		AP[1][3] = revtrans.x13;
		AP[0][4] = revtrans.x04;

		BP[4][0] = revtrans.y40;
		BP[3][1] = revtrans.y31;
		BP[2][2] = revtrans.y22;
		BP[1][3] = revtrans.y13;
		BP[0][4] = revtrans.y04;
	}

	// We can now fill the disprm structure and assign it to wcslib->lin
	struct disprm *dis = calloc(1, sizeof(struct disprm));
	int ipx = 0;
	int dpmax = 10 + 2 * (N + 1) * (N + 2);
	dis->flag = -1;
	// thread-safe version of disini as per WCSLIB documentation.
	// Anyway, dpmax should not change if we solve multiple images in parralel
	// as the trans order should be the same for all images of the sequence
	disinit(1, 2, dis, dpmax);
	char keyword[4];
	char field[30];
	for (int i = 0; i < 2; i++) {
		strcpy(dis->dtype[i], "SIP");
		snprintf(keyword, 4, "DP%d", i + 1);
		dpfill(dis->dp + ipx, keyword, "NAXES", i + 1, 0, 2, 0.0);
		ipx ++;
		dpfill(dis->dp + ipx, keyword, "AXIS.1", i + 1, 0, 1, 0.0);
		ipx ++;
		dpfill(dis->dp + ipx, keyword, "AXIS.2", i + 1, 0, 2, 0.0);
		ipx++;
		dpfill(dis->dp + ipx, keyword, "OFFSET.1", i + 1, 1, 0, wcslib->crpix[0]);
		ipx ++;
		dpfill(dis->dp + ipx, keyword, "OFFSET.2", i + 1, 1, 0, wcslib->crpix[1]);
		ipx++;
		for (int sipflag = 1; sipflag <= 2; sipflag++) { // fwd or rev
			for (int p = 0; p <= N; p++) {
				for (int q = 0; q <= N - p; q++) {
					snprintf(field, 30, "SIP.%s.%d_%d", (sipflag == 1) ? "FWD" : "REV", p, q);
					if (i == 0)
						dpfill(dis->dp + ipx, keyword, field, i + 1, 1, 0, (sipflag == 1) ? A[p][q] : AP[p][q]);
					else
						dpfill(dis->dp + ipx, keyword, field, i + 1, 1, 0, (sipflag == 1) ? B[p][q] : BP[p][q]);
					ipx++;
				}
			}
		}
	}
	dis->ndp = dpmax;
	wcslib->lin.dispre = dis;
	dis->flag = 0;
	wcslib->lin.flag = 0;
	disset(dis);
	linset(&wcslib->lin);
	free_stars(&uvgrid);
	free_stars(&xygrid);
	return 0;
}

static gboolean image_is_flipped_from_wcs(struct wcsprm *wcslib) {
	double cd[2][2];
	wcs_cd2mat(wcslib, cd);
	double det = (cd[0][0] * cd[1][1] - cd[1][0] * cd[0][1]); // determinant of rotation matrix (ad - bc)
	return det > 0; // convention is that angles are positive clockwise when image is not flipped
}

static void flip_bottom_up_astrometry_data(fits *fit) {
	Homography H = { 0 };
	cvGetEye(&H);
	H.h11 = -1.;
	H.h12 = (double)fit->ry;
	reframe_astrometry_data(fit, H);
}

static void print_updated_wcs(struct wcsprm *wcslib) {
	if (!wcslib)
		return;
	/* debug output */
	siril_debug_print("****Current WCS data*************\n");
	siril_debug_print("crpix1 = %*.12e\n", 20, wcslib->crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, wcslib->crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, wcslib->crval[0]);
	siril_debug_print("crval2 = %*.12e\n", 20, wcslib->crval[1]);
	siril_debug_print("cdelt1 = %*.12e\n", 20, wcslib->cdelt[0]);
	siril_debug_print("cdelt2 = %*.12e\n", 20, wcslib->cdelt[1]);
	siril_debug_print("pc1_1  = %*.12e\n", 20, wcslib->pc[0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, wcslib->pc[1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, wcslib->pc[2]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, wcslib->pc[3]);
	if (wcslib->lin.dispre != NULL)
		siril_debug_print("+ SIP terms\n");
	siril_debug_print("******************************************\n");
}

static void print_platesolving_results_from_wcs(struct astrometry_data *args) {
	char field_x[256] = "";
	char field_y[256] = "";
	gboolean report_flip = FALSE;

	if (90. - fabs(args->fit->wcsdata.dec) < 2.78e-3) // center is less than 10" off from a pole
		siril_log_message(_("Up position wrt. N is undetermined (too close to a Pole)\n"));
	else {
		// We move 10" to the North and we'll figure out the angle from there....
		// For some unknown reason, asnet may return a solution with the reference point not at center
		// We need to handle that case by passing ra and dec from args->fit->wcsdata which has been updated to take it into account
		double xN, yN, rotation;
		int status = wcs2pix(args->fit, args->fit->wcsdata.ra, args->fit->wcsdata.dec + 2.78e-3, &xN, &yN);
		xN -= args->fit->rx * 0.5;
		yN -= args->fit->ry * 0.5;
		if (!status) {
			rotation = -atan2(xN, yN) * RADTODEG; // we measure clockwise wrt. +y axis
			if (image_is_flipped_from_wcs(args->fit->wcslib)) {
				if (args->flip_image) {
					rotation = 180.0 - rotation;
				} else {
					report_flip = TRUE; // we only report a flip if the image is not flipped afterwards
				}
			}
			if (rotation < -180.0)
				rotation += 360.0;
			if (rotation > 180.0)
				rotation -= 360.0;
			siril_log_message(_("Up is %+.2lf deg ClockWise wrt. N%s\n"), rotation, report_flip ? _(" (flipped)") : "");
		}
	}
	/* Plate Solving */
	double resolution = get_wcs_image_resolution(args->fit) * 3600.0;
	siril_log_message(_("Resolution:%*.3lf arcsec/px\n"), 11, resolution);
	double focal_length = RADCONV * args->pixel_size / resolution;
	siril_log_message(_("Focal length:%*.2lf mm\n"), 8, focal_length);
	siril_log_message(_("Pixel size:%*.2lf µm\n"), 10, args->pixel_size);
	fov_in_DHMS(resolution * (double)args->fit->rx / 3600.0, field_x);
	fov_in_DHMS(resolution * (double)args->fit->ry / 3600.0, field_y);
	siril_log_message(_("Field of view:    %s x %s\n"), field_x, field_y);
	siril_log_message(_("Image center: alpha: %s, delta: %s\n"), args->fit->wcsdata.objctra, args->fit->wcsdata.objctdec);
	double dist = compute_coords_distance(args->fit->wcsdata.ra, args->fit->wcsdata.dec,
										siril_world_cs_get_alpha(args->cat_center), siril_world_cs_get_delta(args->cat_center));
	siril_log_message(_("Was %.2f arcmin from previous value\n"), dist * 60.);
}

// binomial coeffs up to n=6
// hard-coded for faster execution
static int Cnk(int n, int k)
{
	if(k == 0) return 1;
	if(n < k)  return 0;
	if (n > 6) return 0;
	switch (n) {
		case 0:
			return 1;
		case 1:;
			int v1[2] = { 1, 1 };
			return v1[k];
		case 2:;
			int v2[3] = { 1, 2, 1 };
			return v2[k];
		case 3:;
			int v3[4] = { 1, 3, 3, 1 };
			return v3[k];
		case 4:;
			int v4[5] = { 1, 4, 6, 4, 1 };
			return v4[k];
		case 5:;
			int v5[6] = { 1, 5, 10, 10, 5, 1 };
			return v5[k];
		case 6:;
			int v6[7] = { 1, 6, 15, 20, 15, 6, 1 };
			return v6[k];
		default:
			return 0;
	}
	// Full implementation for any n (https://stackoverflow.com/questions/35121401/binomial-coefficient-in-c-program-explanation)
	// In case we need it someday
	// 	if(k == 0) return 1;
	// 	if(n < k)  return 0;
	// 	return (n * Cnk(n - 1, k - 1)) / k;
}

static void transform_disto_coeff(struct disprm *dis, Homography *H) {
	double A[MAX_SIP_SIZE][MAX_SIP_SIZE] = {{ 0. }};
	double B[MAX_SIP_SIZE][MAX_SIP_SIZE] = {{ 0. }};
	double AP[MAX_SIP_SIZE][MAX_SIP_SIZE] = {{ 0. }};
	double BP[MAX_SIP_SIZE][MAX_SIP_SIZE] = {{ 0. }};
	int N = extract_SIP_order_and_matrices(dis, A, B, AP, BP);
	if (!N)
		return;
	double a = H->h00;
	double b = H->h01;
	double c = H->h10;
	double d = H->h11;
	double det = a * d - b * c;
	double detinv = 1. / det;
	double CD[2][2] = { {detinv * d, -detinv * b}, {-detinv * c, detinv * a}};
	double CDinv[2][2] = { {a, b}, {c, d}};
	a = CD[0][0];
	b = CD[0][1];
	c = CD[1][0];
	d = CD[1][1];
	// pre-computing powers up to N of the a, b, c, d terms
	double powa[MAX_SIP_SIZE], powb[MAX_SIP_SIZE], powc[MAX_SIP_SIZE], powd[MAX_SIP_SIZE];
	for (int i = 0; i <= N; i++) {
		powa[i] = pow(a, i);
		powb[i] = pow(b, i);
		powc[i] = pow(c, i);
		powd[i] = pow(d, i);
	}

	// allocations
	int size = (N + 1) * (N + 2) / 2;
	double *va = calloc(size, sizeof(double));
	double *vb = calloc(size, sizeof(double));
	double *vap = calloc(size, sizeof(double));
	double *vbp = calloc(size, sizeof(double));
	double *tva = calloc(size, sizeof(double));
	double *tvb = calloc(size, sizeof(double));
	double *tvap = calloc(size, sizeof(double));
	double *tvbp = calloc(size, sizeof(double));
	double **t = malloc(sizeof(double*) * size);
	for (int i = 0; i < size; ++i) {
		t[i] = calloc(size, sizeof(double));
	}
	// we will form the matrix T to transform the Aij/Bij/APij/BPij
	// terms to the new coordinates system
	int r = 0;
	for (int i = 0; i <= N; i++) {
		for (int j = i; j >= 0; j--) {
			int k = i - j;
			va[r]  = A[j][k];
			vb[r]  = B[j][k];
			vap[r] = AP[j][k];
			vbp[r] = BP[j][k];
			r++;
		}
	}
	r = 0;
	for (int i = 0; i <= N; i++) {
		// printf("i:%d\n", i);
		r += i;
		for (int j = 0; j <= i; j++) {
			int n = i - j;
			// printf("j(cd):%d , n(ab):%d\n", j, n);
			for (int l = 0; l <= j; l++) {
				double m = Cnk(j, l) * powc[j - l] * powd[l];
				// printf("%d: %d c^%d d^%d\n", l, Cnk(j, l), j - l, l);
				for (int k = 0; k <= i - j; k++) {
					t[r + l + k][r + j] += m * Cnk(n, k) * powa[n - k] * powb[k];
					// printf("%d,%d: %d a^%d b^%d\n", r + l + k, r + j, Cnk(n, k), n - k, k);
				}
			}
		}
	}
	// we can now make the products of matrix with the 4 vectors
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			// printf("%g ", t[i][j]);
			tva[i]  += t[i][j] * va[j];
			tvb[i]  += t[i][j] * vb[j];
			tvap[i] += t[i][j] * vap[j];
			tvbp[i] += t[i][j] * vbp[j];
		}
		// printf("\n");
	}
	// and reapply the inverse linear transform (we reuse the initial vectors)
	for (int i = 0; i < size; i++) {
		double v1[2] = { tva[i], tvb[i] };
		double v2[2] = { tvap[i], tvbp[i] };
		double v1r[2] = { 0., 0. };
		double v2r[2] = { 0., 0. };
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				v1r[j] += CDinv[j][k] * v1[k];
				v2r[j] += CDinv[j][k] * v2[k];
			}
		}
		va[i]  = v1r[0];
		vb[i]  = v1r[1];
		vap[i] = v2r[0];
		vbp[i] = v2r[1];
	}

	r = 0;
	// we redispatch to the matrices
	for (int i = 0; i <= N; i++) {
		for (int j = i; j >= 0; j--) {
			int k = i - j;
			// printf("A%d%d:%g, %g\n", j, k, A[j][k], va[r]);
			// printf("B%d%d:%g, %g\n", j, k, B[j][k], vb[r]);
			A[j][k]  = va[r];
			B[j][k]  = vb[r];
			AP[j][k] = vap[r];
			BP[j][k] = vbp[r];
			r++;
		}
	}
	// and update the dis structure
	update_SIP_keys(dis, A, B, AP, BP);
	dis->flag = 0;
	disset(dis);
	// and free
	free(va);
	free(vb);
	free(vap);
	free(vbp);
	free(tva);
	free(tvb);
	free(tvap);
	free(tvbp);
	for (int i = 0; i < size; ++i) {
		free(t[i]);
	}
	free(t);
}

/******
 *
 * Public functions
 */


void reframe_astrometry_data(fits *fit, Homography H) {
	double pc1_1, pc1_2, pc2_1, pc2_2;
	point refpointout;
	pc1_1 = H.h00 * fit->wcslib->pc[0] + H.h01 * fit->wcslib->pc[1];
	pc1_2 = H.h10 * fit->wcslib->pc[0] + H.h11 * fit->wcslib->pc[1];
	pc2_1 = H.h00 * fit->wcslib->pc[2] + H.h01 * fit->wcslib->pc[3];
	pc2_2 = H.h10 * fit->wcslib->pc[2] + H.h11 * fit->wcslib->pc[3];
	// we go back to cd formulation just to separate back again cdelt and pc
	double cd[2][2], pc[2][2];
	pc[0][0] = pc1_1;
	pc[0][1] = pc1_2;
	pc[1][0] = pc2_1;
	pc[1][1] = pc2_2;
	// we recombine pc and cdelt, and let wcslib decompose from cd
	// pc and cd are set to the new cd
	// cdelt is set to unity
	// we then call wcspcx() to do the decomposition
	wcs_pc_to_cd(pc, fit->wcslib->cdelt, cd);
	wcs_mat2pc(fit->wcslib, cd);
	wcs_mat2cd(fit->wcslib, cd);
	fit->wcslib->altlin = 2;
	wcs_cdelt2unity(fit->wcslib);
	wcspcx(fit->wcslib, 0, 0, NULL);

	// we fetch the refpoint in siril convention
	point refpointin = {fit->wcslib->crpix[0] - 0.5, fit->wcslib->crpix[1] - 0.5};
	cvTransformImageRefPoint(H, refpointin, &refpointout);
	// and convert it back to FITS/WCS convention
	fit->wcslib->crpix[0] = refpointout.x + 0.5;
	fit->wcslib->crpix[1] = refpointout.y + 0.5;

	// and we update all the wcslib structures
	if (fit->wcslib->lin.dispre) {
		transform_disto_coeff(fit->wcslib->lin.dispre, &H);
		struct disprm *dis = fit->wcslib->lin.dispre;
		for (int n = 0; n < dis->ndp; n++) { // update the OFFSET keywords to new CRPIX values
			if (g_str_has_prefix(dis->dp[n].field + 4, "OFFSET.1"))
				dis->dp[n].value.f = fit->wcslib->crpix[0];
			else if (g_str_has_prefix(dis->dp[n].field + 4, "OFFSET.2"))
				dis->dp[n].value.f = fit->wcslib->crpix[1];
		}
		dis->flag = 0; // to update the structure
		disset(dis);
	}
	fit->wcslib->lin.flag = 0; // to update the structure
	linset(&fit->wcslib->lin);
	fit->wcslib->flag = 0; // to update the structure
	wcsset(fit->wcslib);
	wcs_print(fit->wcslib);
	print_updated_wcs(fit->wcslib);

	// Update the center position in fit->wcsdata //
	double rac, decc;
	center2wcs(fit, &rac, &decc);
	if (rac != -1) {
		fit->wcsdata.ra = rac;
		fit->wcsdata.dec = decc;
		gchar *ra = siril_world_cs_alpha_format_from_double(rac, "%02d %02d %.3lf");
		gchar *dec = siril_world_cs_delta_format_from_double(decc, "%c%02d %02d %.3lf");
		g_sprintf(fit->wcsdata.objctra, "%s", ra);
		g_sprintf(fit->wcsdata.objctdec, "%s", dec);
		g_free(ra);
		g_free(dec);
	}
}

void wcs_pc_to_cd(double pc[][2], const double cdelt[2], double cd[][2]) {
	cd[0][0] = pc[0][0] * cdelt[0];
	cd[0][1] = pc[0][1] * cdelt[0];
	cd[1][0] = pc[1][0] * cdelt[1];
	cd[1][1] = pc[1][1] * cdelt[1];
}

static int siril_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution);
static int match_catalog(psf_star **stars, int nb_stars, siril_catalogue *siril_cat, double scale, int order, TRANS *trans_out, double *ra_out, double *dec_out);
static int local_asnet_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution);

static int get_catalog_stars(siril_catalogue *siril_cat) {
	if (siril_catalog_conesearch(siril_cat) <= 0)
		return 1;
	sort_cat_items_by_mag(siril_cat);
	siril_log_message(_("Fetched %d stars from %s catalogue\n"), siril_cat->nbitems, catalog_to_str(siril_cat->cat_index));
	return 0;
}

static psf_star **project_catalog_stars(siril_catalogue *siril_cat, double ra0, double dec0) {
	if (siril_cat->projected > CAT_PROJ_NONE) {
		siril_catalog_reset_projection(siril_cat);
	}
	psf_star **cstars = NULL;
	if (!siril_catalog_project_gnomonic(siril_cat, ra0, dec0, FALSE, NULL)) {
		cstars = convert_siril_cat_to_psf_stars(siril_cat);
	}
	return cstars;
}

/* entry point for plate solving */
gpointer plate_solver(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *) p;
	psf_star **stars = NULL;	// image stars
	int nb_stars = 0;	// number of image and catalogue stars

	args->ret = 0;
	args->message = NULL;
	solve_results solution = { 0 }; // used in the clean-up, init at the beginning

	if (args->verbose) {
		if (args->solver == SOLVER_LOCALASNET) {
			siril_log_message(_("Plate solving image with astrometry.net for a field of view of %.2f degrees\n"), args->used_fov / 60.0);
		} else if (args->ref_stars->cat_index == CAT_LOCAL) {
			siril_log_message(_("Plate solving image from local catalogues for a field of view of %.2f"
						" degrees%s, using a limit magnitude of %.2f\n"),
					args->used_fov / 60.0,
					args->uncentered ? _(" (uncentered)") : "", args->ref_stars->limitmag);
		} else {
			siril_log_message(_("Plate solving image from an online catalogue for a field of view of %.2f"
						" degrees%s, using a limit magnitude of %.2f\n"),
					args->used_fov / 60.0,
					args->uncentered ? _(" (uncentered)") : "", args->ref_stars->limitmag);
		}
	}

	/* 1. Get image stars */
	// store the size of the image being solved for later use in case of downscale
	args->rx_solver = args->fit->rx;
	args->ry_solver = args->fit->ry;
	args->scalefactor = 1.;
	if (!args->manual) {
		int detection_layer = args->fit->naxes[2] == 1 ? 0 : 1;
		fits fit_backup = { 0 };	// original image in case of downscale
		if (args->downsample) {
			int retval;
			fits tmp = { 0 };
			siril_log_message(_("Down-sampling image for faster star detection by a factor %.2f\n"),
					DOWNSAMPLE_FACTOR);
			retval = extract_fits(args->fit, &tmp, detection_layer, FALSE);
			if (!retval) {
				//copy_fits_metadata(args->fit, &tmp);
				args->rx_solver = round_to_int(DOWNSAMPLE_FACTOR * args->fit->rx);
				args->ry_solver = round_to_int(DOWNSAMPLE_FACTOR * args->fit->ry);
				retval = cvResizeGaussian(&tmp, args->rx_solver, args->ry_solver,
						OPENCV_AREA, FALSE);
			}
			if (retval) {
				clearfits(&tmp);
				siril_log_color_message(_("Failed to downsample image, aborting\n"), "red");
				args->message = g_strdup(_("Not enough memory"));
				args->ret = ERROR_PLATESOLVE;
				goto clearup;
			}
			memcpy(&fit_backup, args->fit, sizeof(fits));
			memcpy(args->fit, &tmp, sizeof(fits));

			// TODO: should we average x and y or even better separate scales on x and y?
			args->scalefactor = (double)fit_backup.rx / (double)args->fit->rx;
			detection_layer = 0;
		}

		image im = { .fit = args->fit, .from_seq = NULL, .index_in_seq = -1 };

#ifdef _WIN32
		// on Windows, asnet is not run in parallel neither on single image nor sequence, we can use all threads
		int nthreads = (!args->for_sequence || args->solver == SOLVER_LOCALASNET) ? com.max_thread : 1;
#else
		// on UNIX, asnet is in parallel for sequences, we need to restrain to one per worker
		int nthreads = (!args->for_sequence) ? com.max_thread : 1;
#endif

		stars = peaker(&im, detection_layer, &com.pref.starfinder_conf, &nb_stars,
				&(args->solvearea), FALSE, TRUE,
				BRIGHTEST_STARS, com.pref.starfinder_conf.profile, nthreads);

		if (args->downsample) {
			clearfits(args->fit);
			memcpy(args->fit, &fit_backup, sizeof(fits));
			// we go back to original scale by multiplying stars x/y pos by scalefactor
			if (stars) {
				for (int i = 0; i < nb_stars; i++) {
					stars[i]->xpos *= args->scalefactor;
					stars[i]->ypos *= args->scalefactor;
				}
			}
			args->rx_solver = args->fit->rx;
			args->ry_solver = args->fit->ry;
			args->scalefactor = 1.0;
		}
	} else {
		if (args->stars) {
			stars = args->stars;
			if (stars)
				while (stars[nb_stars])
					nb_stars++;
		} else { // we need to make a copy of com.stars as we will alter the coordinates
			stars = com.stars;
			if (stars) {
				while (com.stars[nb_stars])
					nb_stars++;
				stars = new_fitted_stars(nb_stars);
				for (int s = 0; s < nb_stars; s++) {
					stars[s] = new_psf_star();
					stars[s]->xpos = com.stars[s]->xpos;
					stars[s]->ypos = com.stars[s]->ypos;
				}
			}
		}
	}
	CHECK_FOR_CANCELLATION_RET;

	if (!stars || nb_stars < AT_MATCH_STARTN_LINEAR) {
		args->message = g_strdup_printf(_("There are not enough stars picked in the image. "
				"At least %d are needed."), AT_MATCH_STARTN_LINEAR);
		args->ret = ERROR_PLATESOLVE;
		goto clearup;
	}
	if (args->verbose)
		siril_log_message(_("Using %d detected stars from image.\n"), nb_stars);

	/* 2. Plate solving 
		No modifications done to args or args->fit
		We can read from them but cannot write
	*/
	if (args->solver == SOLVER_LOCALASNET) {
		if (!args->for_sequence) {
			com.child_is_running = EXT_ASNET;
			g_unlink("stop"); // make sure the flag file for cancel is not already in the folder
		}
		if (local_asnet_platesolve(stars, nb_stars, args, &solution)) {
			args->ret = ERROR_PLATESOLVE;
		}
	} else {
		double x0 = args->fit->rx * 0.5;
		double y0 = args->fit->ry * 0.5;
		for (int s = 0; s < nb_stars; s++) {
			stars[s]->xpos -= x0;
			stars[s]->ypos = y0 - stars[s]->ypos;
		}
		if (siril_platesolve(stars, nb_stars, args, &solution)) {
			args->ret = ERROR_PLATESOLVE;
		}
	}
	if (args->ret)
		goto clearup;

	/* 4. Print and store some results in args->fit*/
	// saving state for undo before modifying fit structure
	if (!com.script) {
		const char *undo_str = _("Plate Solve");
		undo_save_state(args->fit, undo_str);
	}
	// we copy the solution and update wcsdata
	if (has_wcs(args->fit))
		wcsfree(args->fit->wcslib);
	args->fit->wcslib = solution.wcslib;
	print_updated_wcs(args->fit->wcslib);
	update_wcsdata_from_wcs(args);
	if (args->verbose)
		print_platesolving_results_from_wcs(args);
	double resolution = get_wcs_image_resolution(args->fit) * 3600.0;
	double focal_length = RADCONV * args->pixel_size / resolution;
	args->fit->focal_length = focal_length;
	args->fit->pixel_size_x = args->pixel_size;
	args->fit->pixel_size_y = args->pixel_size;
	if (!args->for_sequence && com.pref.astrometry.update_default_scale) {
		com.pref.starfinder_conf.focal_length = focal_length;
		com.pref.starfinder_conf.pixel_size_x = args->pixel_size;
		siril_log_message(_("Saved focal length %.2f and pixel size %.2f as default values\n"), focal_length, args->pixel_size);
	}

	/* 5. Flip image if needed */
	if (args->flip_image && image_is_flipped_from_wcs(args->fit->wcslib)) {
		if (args->verbose)
			siril_log_color_message(_("Flipping image and updating astrometry data.\n"), "salmon");
		fits_flip_top_to_bottom(args->fit);
		flip_bottom_up_astrometry_data(args->fit);
		args->image_flipped = TRUE;
	}

clearup:
	if (stars) {
		for (int i = 0; i < nb_stars; i++)
			free_psf(stars[i]);
		free(stars);
	}
	if (args->cat_center)
		siril_world_cs_unref(args->cat_center);
	if (args->ref_stars)
		siril_catalog_free(args->ref_stars);
	g_free(args->filename);
	int retval = args->ret;
	if (com.script && retval) {
		siril_log_message(_("Plate solving failed: %s\n"), args->message);
		g_free(args->message);
	}
	if (!args->for_sequence) {
		if (com.child_is_running == EXT_ASNET && g_unlink("stop"))
			siril_debug_print("g_unlink() failed\n");
		com.child_is_running = EXT_NONE;
		siril_add_idle(end_plate_solver, args);
	}
	else free(args);
	return GINT_TO_POINTER(retval);
}

typedef enum {
	SOLVE_LASTSOLVE = -1,
	SOLVE_OK,
	SOLVE_NO_MATCH, //solution->message = g_strdup(_("Could not match stars from the catalogue"));
	SOLVE_CANCELLED, //solution->message = g_strdup(_("Cancelled"));
	SOLVE_INVALID_TRANS, 		//solution->message = g_strdup(_("Transformation matrix is invalid, solve failed"));
	SOLVE_PROJ, //solution->message = g_strdup(_("Reprojecting catalog failed."));
	SOLVE_TRANS_UPDATE, //g_strdup(_("Updating trans failed."));
} siril_solve_errcode;

// returns a list of center points to be searched
static point *get_centers(double fov_deg, double search_radius_deg, double ra0, double dec0, int *nb) {
	point *centers;
	int N;
	double radius = fov_deg * 0.5 * DEGTORAD;
	int n = (int)ceil(2. * search_radius_deg / fov_deg); // number of dec rings
	double *d = malloc(n * sizeof(double));
	int *nl = malloc(n * sizeof(int));
	N = 0;
	for (int i = 0; i < n; i++) {
		d[i] = M_PI_2 - i * radius;
		nl[i] = (int)ceil(2. * M_PI * sin(i * radius) / radius); // number of points of ith ring
		N += nl[i];
	}
	nl[0] = 1;
	N++;
	point *centers0 = malloc(N * sizeof(point)); // points in native coordinates
	centers = malloc(N * sizeof(point)); // points in celestial coordinates
	int s = 0;
	double init = 0.;
	// setting the points in native coordinates
	for (int i = 0; i < n; i++) {
		double pace = 2. * M_PI / (double)nl[i];
		for (int j = 0; j < nl[i]; j++) {
			centers0[s + j].x = init + j * pace;
			centers0[s + j].y = d[i];
		}
		if (s > 0) {
			init = 0.5 * (centers0[s].x + centers0[s + 1].x);
		}
		s += nl[i];
	}
	//converting to celestial coordinates
	ra0 *= DEGTORAD;
	dec0 *= DEGTORAD;
	double cos_d0 = cos(dec0);
	double sin_d0 = sin(dec0);
	for (int i = 0; i < N; i++) {
		double cos_p = cos(centers0[i].x);
		double sin_p = sin(centers0[i].x);
		double cos_t = cos(centers0[i].y);
		double sin_t = sin(centers0[i].y);
		// Using eq (2) of paper II from Calabretta and Greisen
		// with phi_p = 180
		centers[i].x = (ra0 + atan2(cos_t * sin_p, sin_t * cos_d0 + cos_t * sin_d0 * cos_p)) * RADTODEG;
		centers[i].y = asin(sin_t * sin_d0 - cos_t * cos_d0 * cos_p) * RADTODEG;
	}
	free(d);
	free(nl);
	free(centers0);
	*nb = N;
	siril_debug_print("%d center points to test\n", N);
	return centers;
}

/* entry point for siril's plate solver based on catalogue matching
   args->fit is not modified by this function
   Only solution is filled
*/
static int siril_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution) {
	point *centers;
	int N;
	double ra0 = siril_world_cs_get_alpha(args->cat_center);
	double dec0 = siril_world_cs_get_delta(args->cat_center);
	if (args->blind) {
		if (!args->blindradius)
			args->blindradius = 3.;
		centers = get_centers(args->used_fov / 60., args->blindradius, ra0, dec0, &N);
	} else {
		centers = malloc(sizeof(point));
		centers[0].x = ra0;
		centers[0].y = dec0;
		N = 1;
	}
	gboolean solved = FALSE;
	siril_catalogue *siril_cat = NULL;
	int n = 0;
	TRANS t = { 0 };
	int ret = 1;
	double ra = -1., dec = -1.;
	while (!solved && n < N) {
		if (args->for_sequence && !args->blind) {
			siril_cat = args->ref_stars; // the catalogue has already been fetched in the prepare_hook and copied in the image_hook
		} else {
			if (n == 0) {
				siril_cat = calloc(1, sizeof(siril_catalogue));
				siril_catalogue_copy(args->ref_stars, siril_cat); // this only copies the query parameters, the catalogue has not yet been fetched
			}
			siril_catalog_free_items(siril_cat);
			siril_cat->center_ra = centers[n].x;
			siril_cat->center_dec = centers[n].y;
			get_catalog_stars(siril_cat);
		}
		ret = match_catalog(stars, nb_stars, siril_cat, args->scale, args->trans_order, &t, &ra, &dec);
		if (!ret) { // we update the solution - but if blind, we do a last solve with a new catalogue fetched at the center
			if (args->blind) {
				siril_catalog_free_items(siril_cat);
				siril_cat->center_ra = ra;
				siril_cat->center_dec = dec;
				get_catalog_stars(siril_cat);
				int ret2;
				TRANS t2 = { 0 };
				double ra2 = -1., dec2 = -1.;
				ret2 = match_catalog(stars, nb_stars, siril_cat, args->scale, args->trans_order, &t2, &ra2, &dec2);
				if (!ret2) {
					t = t2;
					ra = ra2;
					dec = dec2;
				} else {
					ret = SOLVE_LASTSOLVE;
				}
			}
			double cd[2][2];
			get_cd_from_trans(&t, cd);
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					cd[i][j] *= ASECTODEG;
				}
			}
			/**** Fill solution wcslib structure ***/
			wcsprm_t *prm = calloc(1, sizeof(wcsprm_t));
			prm->flag = -1;
			wcsinit(1, 2, prm, 0, 0, 0);
			prm->equinox = 2000.0;
			prm->crpix[0] = (double)args->rx_solver * 0.5;
			prm->crpix[1] = (double)args->ry_solver * 0.5;
			// we now go back to FITS convention (from siril)
			prm->crpix[0] += 0.5;
			prm->crpix[1] += 0.5;
			prm->crval[0] = ra;
			prm->crval[1] = dec;
			prm->lonpole = 180.;
			if (t.order == AT_TRANS_LINEAR) {
					const char CTYPE[2][9] = { "RA---TAN", "DEC--TAN" };
					const char CUNIT[2][4] = { "deg", "deg" };
					for (int i = 0; i < NAXIS; i++) {
						strncpy(prm->cunit[i], &CUNIT[i][0], 71); // 72 char fixed buffer, keep 1 for the NULL
					}
					for (int i = 0; i < NAXIS; i++) {
						strncpy(prm->ctype[i], &CTYPE[i][0], 71); // 72 byte buffer, leave 1 byte for the NULL
					}
			} else {
					const char CTYPE[2][12] = { "RA---TAN-SIP", "DEC--TAN-SIP" };
					const char CUNIT[2][4] = { "deg", "deg" };
					for (int i = 0; i < NAXIS; i++) {
						strncpy(prm->cunit[i], &CUNIT[i][0], 71); // 72 char fixed buffer, keep 1 for the NULL
					}
					for (int i = 0; i < NAXIS; i++) {
						strncpy(prm->ctype[i], &CTYPE[i][0], 71); // 72 byte buffer, leave 1 byte for the NULL
					}
				}
			// we pass cd and let wcslib decompose it
			wcs_mat2cd(prm, cd);
			prm->altlin = 2;
			wcspcx(prm, 0, 0, NULL);
			if (t.order > AT_TRANS_LINEAR)
				add_disto_to_wcslib(prm, &t, args->rx_solver, args->ry_solver);
			wcs_print(prm);
			solution->wcslib = prm;
			if (args->verbose)
				siril_log_color_message(_("Siril internal solver succeeded.\n"), "green");
			solved = TRUE;
		} else if (ret == SOLVE_CANCELLED) {
			solved = TRUE; // just to break the while
		}
		n++;
	}
	if (ret) {
		args->message = g_strdup("failed\n");
	}
	free(centers);
	if (siril_cat != args->ref_stars)
		siril_catalog_free(siril_cat);
	return (ret > 0);
}

static int match_catalog(psf_star **stars, int nb_stars, siril_catalogue *siril_cat, double scale, int order, TRANS *trans_out, double *ra_out, double *dec_out) {
	TRANS trans = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;
	int max_trials = 0;
	s_star *star_list_A = NULL, *star_list_B = NULL;
	psf_star **cstars = NULL;
	int n_cat = siril_cat->nbitems;
	double ra0 = siril_cat->center_ra;
	double dec0 = siril_cat->center_dec;
	int ret = 1;

	max_trials = 20; //always try to converge to the center

	/* make sure that arrays are not too small
	 * make sure that the max of stars is BRIGHTEST_STARS */
	int n = min(min(nb_stars, n_cat), BRIGHTEST_STARS);

	double a = 1.0 + (com.pref.astrometry.percent_scale_range / 100.0);
	double b = 1.0 - (com.pref.astrometry.percent_scale_range / 100.0);
	double scale_min = 1.0 / (scale * a);
	double scale_max = 1.0 / (scale * b);
	int attempt = 1;

	cstars = project_catalog_stars(siril_cat, ra0, dec0);
	if (!cstars) {
		return SOLVE_PROJ;
	}
	while (ret && attempt <= 1) {
		free_stars(&star_list_A);
		free_stars(&star_list_B);
		ret = new_star_match(stars, cstars, n, nobj,
				scale_min, scale_max, NULL, &trans, TRUE,
				UNDEFINED_TRANSFORMATION, AT_TRANS_LINEAR, &star_list_A, &star_list_B);
		if (attempt == 2) {
			scale_min = -1.0;
			scale_max = -1.0;
		} else {
			nobj += 50;
		}
		attempt++;
		CHECK_FOR_CANCELLATION;
	}
	if (ret) {
		goto clearup;
	}
	if (!check_affine_TRANS_sanity(&trans)) {
		ret = SOLVE_INVALID_TRANS;
		goto clearup;
	}

	int num_matched = trans.nm;
	int trial = 0;
	double conv = DBL_MAX;

	/* try to get a better solution */
	conv = get_center_offset_from_trans(&trans);
	siril_debug_print("iteration %d - offset: %.3f, number of matches: %d\n", trial, conv, trans.nr);
	while (conv > CONV_TOLERANCE && trial < max_trials){
		// we get the new projection center
		apply_match(ra0, dec0, 0., 0., &trans, &ra0, &dec0);

		// we will reproject the catalog at the new image center
		siril_debug_print("Reprojecting to: alpha: %s, delta: %s\n",
				siril_world_cs_alpha_format_from_double(ra0, "%02d %02d %.3lf"),
				siril_world_cs_delta_format_from_double(dec0, "%c%02d %02d %.3lf"));

		// Uses the indexes in star_list_B to update the stars positions according to the new projection
		free_fitted_stars(cstars);
		cstars = project_catalog_stars(siril_cat, ra0, dec0);
		if (!cstars) {
			ret = SOLVE_PROJ;
			break;
		}
		update_stars_positions(&star_list_B, num_matched, cstars);
		// and recompute the trans structure
		if (atRecalcTrans(num_matched, star_list_A, num_matched, star_list_B, AT_MATCH_MAXITER, AT_MATCH_HALTSIGMA, &trans)) {
			ret = SOLVE_TRANS_UPDATE;
			break;
		}
		num_matched = trans.nm;
		conv = get_center_offset_from_trans(&trans);
		trial++;
		siril_debug_print("iteration %d - offset: %.3f, number of matches: %d\n", trial, conv, trans.nr);
	}
	if (ret)	// after the break
		goto clearup;
	debug_print_catalog_files(star_list_A, star_list_B);

	// updating solution to higher order if required
	if (order > AT_TRANS_LINEAR) {
		double ra1 = ra0;
		double dec1 = dec0;
		siril_debug_print("starting non linear match at order %d\n", order);
		TRANS newtrans = { 0 }; // we will fall back on the linear solution if not suceessfull at higher order
		newtrans.order = order;
		int ret2 = atRecalcTrans(num_matched, star_list_A, num_matched, star_list_B, AT_MATCH_MAXITER, AT_MATCH_HALTSIGMA, &newtrans);
		if (!ret2) {
			conv = get_center_offset_from_trans(&newtrans);
			trial = 0;
			siril_debug_print("iteration %d - offset: %.3f, number of matches: %d\n", trial, conv, newtrans.nr);
			while (conv > CONV_TOLERANCE && trial < max_trials && !ret2) {
				// we get the new projection center
				apply_match(ra1, dec1, 0., 0., &newtrans, &ra1, &dec1);
				// we will reproject the catalog at the new image center
				siril_debug_print("Reprojecting to: alpha: %s, delta: %s\n",
					siril_world_cs_alpha_format_from_double(ra1, "%02d %02d %.3lf"),
					siril_world_cs_delta_format_from_double(dec1, "%c%02d %02d %.3lf"));

				// Uses the indexes in star_list_B to update the stars positions according to the new projection
				free_fitted_stars(cstars);
				cstars = project_catalog_stars(siril_cat, ra1, dec1);
				if (!cstars) {
					ret2 = SOLVE_PROJ;
					break;
				}
				update_stars_positions(&star_list_B, num_matched, cstars);
				// and recompute the trans structure
				if (atRecalcTrans(num_matched, star_list_A, num_matched, star_list_B, AT_MATCH_MAXITER, AT_MATCH_HALTSIGMA, &newtrans)) {
					ret2 = SOLVE_TRANS_UPDATE;
					break;
				}
				conv = get_center_offset_from_trans(&newtrans);
				trial++;
				print_trans(&newtrans);
				siril_debug_print("iteration %d - offset: %.3f, number of matches: %d\n", trial, conv, newtrans.nr);
			}
		}
		if (!ret2) { // higher-order solve was successful, we update
			trans = newtrans;
			ra0 = ra1;
			dec0 = dec1;
		} else {
			siril_log_color_message(_("%s could not find distortion polynomials for the order specified (%d) and returned a linear solution, try with a lower order\n"), "red", "Siril", order);
		}
	}

	if (trial == max_trials) {
		siril_debug_print("No convergence found: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
	} else {
		siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f at iteration #%d\n", ra0, dec0, trial);
	}

	*trans_out = trans;
	*ra_out = ra0;
	*dec_out = dec0;

clearup:
	free_fitted_stars(cstars);
	free_stars(&star_list_A);
	free_stars(&star_list_B);
	return ret;
}

/*********************** finding asnet bash first **********************/

// Retrieves and caches asnet_version
static void get_asnet_version(gchar *path) {
	gchar* bin[3];
	gchar* child_stdout = NULL;
	bin[0] = path;
	bin[1] = "--version";
	bin[2] = NULL;
	g_autoptr(GError) error = NULL;
	g_spawn_sync(NULL, bin, NULL, G_SPAWN_SEARCH_PATH | G_SPAWN_STDERR_TO_DEV_NULL, NULL, NULL, &child_stdout, NULL, NULL, &error);
	if (error) {
		// This will happen on Windows, for users using ansvr but having set the path to cygwin_ansvr in the prefs (instead of leaving blank)
		// We'll not make this over complicated and fallback assuming this is the only case
		siril_debug_print("%s\n",error->message);
		asnet_version = g_strdup("ansvr");
	} else {
		gchar** chunks = g_strsplit(child_stdout, "\n", 2);
		if (strlen(chunks[1]) == 0) { // if the second chunk is void, it means there was only one line which contains the version
			asnet_version = g_strdup(chunks[0]);
		} else {
			asnet_version = g_strdup("<0.88"); // version switch was introduced from 0.88 (https://github.com/dstndstn/astrometry.net/commit/25f0b829d80c57984d404de50f9c645cb3da2223)
		}
		g_strfreev(chunks);
	}
	g_free(child_stdout);
	siril_debug_print("Running asnet version %s\n", asnet_version);
}

void reset_asnet_version() {
	if (asnet_version) {
		g_free(asnet_version);
		asnet_version = NULL;
	}
}

#ifdef _WIN32
static gchar *siril_get_asnet_bash() {
	// searching user-defined path if any
	if (com.pref.asnet_dir && com.pref.asnet_dir[0] != '\0') {
		gchar *testdir = g_build_filename(com.pref.asnet_dir, "bin", NULL);
		// only testing for dir existence, which will catch most path defintion errors
		// this is lighter than testing for existence of bash.exe with G_FILE_TEST_IS_EXECUTABLE flag
		if (!g_file_test(testdir, G_FILE_TEST_IS_DIR)) {
			siril_log_color_message(_("cygwin/bin was not found at %s - ignoring\n"), "red", testdir);
			g_free(testdir);
		} else {
			siril_debug_print("cygwin/bin found at %s\n", testdir);
			g_free(testdir);
			gchar *versionpath = g_build_filename(com.pref.asnet_dir, "bin", "solve-field", NULL);
			if (!asnet_version)
				get_asnet_version(versionpath);
			g_free(versionpath);
			return g_build_filename(com.pref.asnet_dir, NULL);
		}
	}
	// searching default location %localappdata%/cygwin_ansvr
	const gchar *localappdata = g_get_user_data_dir();
	gchar *testdir = g_build_filename(localappdata, "cygwin_ansvr", "bin", NULL);
	if (g_file_test(testdir, G_FILE_TEST_IS_DIR)) {
		siril_debug_print("cygwin/bin found at %s\n", testdir);
		g_free(testdir);
		if (!asnet_version) {
			asnet_version = g_strdup("ansvr");
			siril_debug_print("Running asnet version %s\n", asnet_version);
		}
		return g_build_filename(localappdata, "cygwin_ansvr", NULL);
	}
	siril_log_color_message(_("cygwin/bin was not found at %s - ignoring\n"), "red", testdir);
	g_free(testdir);
	return NULL;
}

gboolean asnet_is_available() {
	gchar *path = siril_get_asnet_bash();
	gboolean retval = path != NULL;
	g_free(path);
	return retval;
}

#else
static gboolean solvefield_is_in_path = FALSE;
static gchar *siril_get_asnet_bin() {
	if (solvefield_is_in_path)
		return g_strdup("solve-field");
	if (!com.pref.asnet_dir || com.pref.asnet_dir[0] == '\0')
		return NULL;
	return g_build_filename(com.pref.asnet_dir, "solve-field", NULL);
}

/* returns true if the command solve-field is available */
gboolean asnet_is_available() {
	const char *str = "solve-field -h > /dev/null 2>&1";
	int retval = system(str);
	if (WIFEXITED(retval) && (0 == WEXITSTATUS(retval))) {
		solvefield_is_in_path = TRUE;
		siril_debug_print("solve-field found in PATH\n");
		if (!asnet_version)
			get_asnet_version("solve-field");
		return TRUE;
	}
	siril_debug_print("solve-field not found in PATH\n");
	gchar *bin = siril_get_asnet_bin();
	if (!bin) return FALSE;
	gboolean is_available = g_file_test(bin, G_FILE_TEST_EXISTS);
	if (!asnet_version)
		get_asnet_version(bin);
	g_free(bin);
	return is_available;
}
#endif

static int local_asnet_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution) {
#ifdef _WIN32
	gchar *asnet_shell = siril_get_asnet_bash();
	if (!asnet_shell) {
		return 1;
	}
#else
	if (!args->asnet_checked) {
		if (!asnet_is_available()) {
			siril_log_color_message(_("solve-field was not found, set its path in the preferences\n"), "red");
			return 1;
		}
	}
#endif

	gchar *table_filename = replace_ext(args->filename, ".xyls");
#ifdef _WIN32
	gchar *stopfile = g_build_filename(com.wd, "stop", NULL);
	if (!g_path_is_absolute(table_filename)) {
		gchar *tmp = g_build_filename(com.wd, table_filename, NULL);
		g_free(table_filename);
		table_filename = tmp;
	}
#else
	gchar *stopfile = g_strdup("stop");
#endif
	if (save_list_as_FITS_table(table_filename, stars, nb_stars, args->rx_solver, args->ry_solver)) {
		siril_log_message(_("Failed to create the input data for solve-field\n"));
		g_free(table_filename);
		return 1;
	}

	char low_scale[16], high_scale[16], time_limit[16];
	double a = 1.0 + (com.pref.astrometry.percent_scale_range / 100.0);
	double b = 1.0 - (com.pref.astrometry.percent_scale_range / 100.0);
	sprintf(low_scale, "%.3f", args->scale * b);
	sprintf(high_scale, "%.3f", args->scale * a);
	sprintf(time_limit, "%d", com.pref.astrometry.max_seconds_run);
#ifndef _WIN32
	gchar *asnet_path = siril_get_asnet_bin();
	g_assert(asnet_path);
#endif

	char *sfargs[50] = {
#ifdef _WIN32
		"solve-field", "-C", "\"$c\"",
		// the stop file must be passed in asnet.sh to be properly quoted and called
		// in case there are spaces in cwd
#else
		asnet_path, "-C", stopfile,
		"--temp-axy",	// not available in the old version of ansvr
#endif
		"-p", "-O", "-N", "none", "-R", "none", "-M", "none", "-B", "none",
		"-U", "none", "-S", "none", "--crpix-center", "-l", time_limit,
		"-u", "arcsecperpix", "-L", low_scale, "-H", high_scale, "-s", "FLUX", NULL };

	char order[12];	// referenced in sfargs, needs the same scope
	if (com.pref.astrometry.sip_correction_order > 1) {
		sprintf(order, "%d", com.pref.astrometry.sip_correction_order);
		char *tweak_args[] = { "-t", order, NULL };
		append_elements_to_array(sfargs, tweak_args);
	} else {
		char *tweak_args[] = { "-T", NULL };
		append_elements_to_array(sfargs, tweak_args);
	}

	char start_ra[16], start_dec[16], radius[16];
	if (args->cat_center) {
		sprintf(start_ra, "%f", siril_world_cs_get_alpha(args->cat_center));
		sprintf(start_dec, "%f", siril_world_cs_get_delta(args->cat_center));
		sprintf(radius, "%.1f", com.pref.astrometry.radius_degrees);
		char *additional_args[] = { "--ra", start_ra, "--dec", start_dec,
			"--radius", radius, NULL};
		append_elements_to_array(sfargs, additional_args);
		siril_log_message(_("Astrometry.net solving with a search field at RA: %s, Dec: %s,"
					" within a %s degrees radius for scales [%s, %s]\n"),
				start_ra, start_dec, radius, low_scale, high_scale);
	} else {
		siril_log_message(_("Astrometry.net solving blindly for scales [%s, %s]\n"),
				low_scale, high_scale);
	}
#ifdef _WIN32
	char *file_args[] = { "\"$p\"", NULL };
#else
	char *file_args[] = { (char*)table_filename, NULL };
#endif
	append_elements_to_array(sfargs, file_args);

	gchar *command = build_string_from_words(sfargs);
	siril_debug_print("Calling solve-field:\n%s\n", command);

#ifdef _WIN32
	// in order to be compatible with different asnet cygwin builds
	// we need to send the command through a bash script
	// the script is written to the /tmp folder (in cygwin env)
	// and called with: /path/to/cygwin/bin/bash -l -c /tmp/asnet.sh
	gchar *asnetscript = g_build_filename(asnet_shell, "tmp", "asnet.sh", NULL);
	g_unlink(asnetscript);
	FILE* tmpfd = g_fopen(asnetscript, "wb+");
	if (tmpfd == NULL) {
		fprintf(stderr,"cannot create temporary file: exiting solve-field");
		g_free(asnetscript);
		g_free(command);
		return 1;
	}
	/* Write data to this file  */
	fprintf(tmpfd, "p=\"%s\"\n", (char*)table_filename);
	fprintf(tmpfd, "c=\"%s\"\n", (char*)stopfile);
	fprintf(tmpfd, "%s\n", command);
	fclose(tmpfd);
	g_free(asnetscript);
	gchar *asnet_bash = g_build_filename(asnet_shell, "bin", "bash", NULL);
	memset(sfargs, '\0', sizeof(sfargs));
	char *newargs[] = {asnet_bash, "-l", "-c", "/tmp/asnet.sh", NULL};
	append_elements_to_array(sfargs, newargs);
#endif
	g_free(command);

	/* call solve-field */
	gint child_stdout;
	g_autoptr(GError) error = NULL;

	g_spawn_async_with_pipes(NULL, sfargs, NULL,
			G_SPAWN_LEAVE_DESCRIPTORS_OPEN | G_SPAWN_SEARCH_PATH,
			NULL, NULL, NULL, NULL, &child_stdout, NULL, &error);
	if (error != NULL) {
		siril_log_color_message("Spawning solve-field failed: %s\n", "red", error->message);
		if (!com.pref.astrometry.keep_xyls_files)
			if (g_unlink(table_filename))
				siril_debug_print("Error unlinking table_filename\n");
		g_free(table_filename);
		g_free(stopfile);
#ifdef _WIN32
		g_free(asnet_shell);
#else
		g_free(asnet_path);
#endif
		return 1;
	}

	GInputStream *stream = NULL;
#ifdef _WIN32
	stream = g_win32_input_stream_new((HANDLE)_get_osfhandle(child_stdout), FALSE);
#else
	stream = g_unix_input_stream_new(child_stdout, FALSE);
#endif
	gboolean success = FALSE;
	gchar *buffer;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
					NULL, NULL))) {
		if (com.pref.astrometry.show_asnet_output)
			siril_log_message("solve-field: %s\n", buffer);
		else siril_debug_print("solver: %s\n", buffer);
		if (g_str_has_prefix(buffer, "Did not solve")) {
			siril_log_color_message(_("No astrometric solution found\n"), "red");
			g_free(buffer);
			break;
		}
		if (g_str_has_prefix(buffer, "Field center: (RA,Dec)")) {
			siril_debug_print("Found a solution, waiting for EOF and exit\n");
			success = TRUE;
		}
		g_free(buffer);
	}
	g_object_unref(data_input);
	g_object_unref(stream);
	if (!g_close(child_stdout, &error))
		siril_debug_print("%s\n", error->message);
	if (!com.pref.astrometry.keep_xyls_files)
		if (g_unlink(table_filename)) {
			siril_debug_print("Error unlinking table_filename\n");
		}
	g_free(table_filename);
	g_free(stopfile);
#ifdef _WIN32
	g_free(asnet_shell);
#else
	g_free(asnet_path);
#endif
	if (!success)
		return 1;

	/* get the results from the .wcs file */
	gchar *wcs_filename = replace_ext(args->filename, ".wcs");
	fits result = { 0 };
	if (read_fits_metadata_from_path_first_HDU(wcs_filename, &result)) {
		siril_log_color_message(_("Could not read the solution from solve-field (expected in file %s)\n"), "red", wcs_filename);
		return 1;
	}

	solution->wcslib = result.wcslib;
	result.wcslib = NULL;
	clearfits(&result);
	if (!com.pref.astrometry.keep_wcs_files)
		g_unlink(wcs_filename);
	g_free(wcs_filename);

	// In some cases asnet returns a dis struct with all coeffs null
	if (solution->wcslib->lin.dispre) { // some distortions were calculated, checked that the terms are not all null
		int N = extract_SIP_order_and_matrices(solution->wcslib->lin.dispre, NULL, NULL, NULL, NULL);
		if (!N) { // the computation of the distortions has failed for the order specified, we remove it and warn the user
			disfree(solution->wcslib->lin.dispre);
			solution->wcslib->lin.dispre = NULL;
			solution->wcslib->flag = 0;
			wcsset(solution->wcslib);
			siril_log_color_message(_("%s could not find distortion polynomials for the order specified (%d) and returned a linear solution, try with a lower order\n"), "red", "astrometry.net", com.pref.astrometry.sip_correction_order);
		}
	}
	// In other cases, the dis struct is empyty, we still need to warn the user
	if (com.pref.astrometry.sip_correction_order > 1 && !solution->wcslib->lin.dispre) {
		siril_log_color_message(_("%s could not find distortion polynomials for the order specified (%d) and returned a linear solution, try with a lower order\n"), "red", "astrometry.net", com.pref.astrometry.sip_correction_order);
	}
	if (args->verbose)
		siril_log_color_message(_("Local astrometry.net solve succeeded.\n"), "green");

	args->ret = 0;
	return 0;
}

// inputs: focal length, pixel size, manual, fit, autocrop, downsample, mag_mode and mag_arg
// outputs: scale, used_fov, uncentered, solvearea, limit_mag
void process_plate_solver_input(struct astrometry_data *args) {
	args->scale = get_resolution(args->focal_length, args->pixel_size);

	rectangle croparea = { 0 };
	if (!args->manual) {
		// first checking if there is a selection or if the full field is to be used
		if (com.selection.w != 0 && com.selection.h != 0) {
			memcpy(&croparea, &com.selection, sizeof(rectangle));
			siril_log_color_message(_("Warning: using the current selection to detect stars\n"), "salmon");
		} else {
			croparea.x = 0;
			croparea.y = 0;
			croparea.w = args->fit->rx;
			croparea.h = args->fit->ry;
		}
		double fov_arcmin = get_fov_arcmin(args->scale, croparea.w, croparea.h);
		siril_debug_print("image fov for given sampling: %f arcmin\n", fov_arcmin);

		// then apply or not autocropping to 5deg (300 arcmin)
		args->used_fov = args->autocrop ? min(fov_arcmin, 300.) : fov_arcmin;
		double cropfactor = (args->used_fov < fov_arcmin) ? args->used_fov / fov_arcmin : 1.0;
		if (cropfactor != 1.0) {
			croparea.x += (int) ((croparea.w - croparea.w * cropfactor) / 2);
			croparea.y += (int) ((croparea.h - croparea.h * cropfactor) / 2);
			croparea.w = (int) (cropfactor * croparea.w);
			croparea.h = (int) (cropfactor * croparea.h);
			siril_debug_print("Auto-crop factor: %.2f\n", cropfactor);
		}

		if (com.selection.w != 0 && com.selection.h != 0) {
			// detect if the selection is not centered enough that it matters
			double thr = max(args->fit->rx, args->fit->ry) / 10.0;
			args->uncentered =
				fabs(croparea.x + 0.5 * croparea.w - 0.5 * args->fit->rx) > thr ||
				fabs(croparea.y + 0.5 * croparea.h - 0.5 * args->fit->ry) > thr;
			if (args->uncentered)
				siril_debug_print("detected uncentered selection\n");
			else siril_debug_print("selection considered centered\n");
		} else {
			args->uncentered = FALSE;
		}

		if (args->downsample) {
			croparea.w *= DOWNSAMPLE_FACTOR;
			croparea.h *= DOWNSAMPLE_FACTOR;
			croparea.x *= DOWNSAMPLE_FACTOR;
			croparea.y *= DOWNSAMPLE_FACTOR;
		}
	} else { //stars manual selection - use full field centered
		args->used_fov = get_fov_arcmin(args->scale, args->fit->rx, args->fit->ry);
		args->uncentered = FALSE;
		if (com.selection.w != 0 && com.selection.h != 0)
			siril_log_message(_("Selection is not used in manual star selection mode\n"));
		// TODO: we could actually check if stars are in the selection
	}
	args->ref_stars->radius = args->used_fov * 0.5;

	if (croparea.w == args->fit->rx && croparea.h == args->fit->ry)
		memset(&croparea, 0, sizeof(rectangle));
	else siril_debug_print("reduced area for the solve: %d, %d, %d x %d%s\n",
			croparea.x, croparea.y, croparea.w, croparea.h,
			args->downsample ? " (down-sampled)" : "");
	memcpy(&(args->solvearea), &croparea, sizeof(rectangle));

	compute_limit_mag(args); // to call after having set args->used_fov
	if (args->ref_stars->cat_index == CAT_AUTO) {
		if (args->ref_stars->limitmag <= 12.5)
			args->ref_stars->cat_index = CAT_TYCHO2;
		else if (args->ref_stars->limitmag <= 17.0)
			args->ref_stars->cat_index = CAT_NOMAD;
		else args->ref_stars->cat_index = CAT_GAIADR3;
	} else if (args->ref_stars->cat_index == CAT_LOCAL && args->ref_stars->limitmag > 17.0) { // TODO: not sure about this one, but otherwise, we can't solve long focals with local cats installed
		args->ref_stars->cat_index = CAT_GAIADR3;
	}
	solve_is_blind(args); // sets flag to check if the catalogue needs to be fetched at the start or if it will be fetched on demand
}

static int astrometry_prepare_hook(struct generic_seq_args *arg) {
	struct astrometry_data *args = (struct astrometry_data *)arg->user;
	fits fit = { 0 };
	// load ref metadata in fit
	if (seq_read_frame_metadata(arg->seq, sequence_find_refimage(arg->seq), &fit))
		return 1;
	seq_prepare_hook(arg);
	args->fit = &fit;
	process_plate_solver_input(args); // compute required data to get the catalog
	clearfits(&fit);
	args->fit = NULL;
	if (args->solver == SOLVER_LOCALASNET) {
		com.child_is_running = EXT_ASNET;
		g_unlink("stop"); // make sure the flag file for cancel is not already in the folder
	}
	if (!args->blind)
		return get_catalog_stars(args->ref_stars);
	return 0;
}

static int astrometry_image_hook(struct generic_seq_args *arg, int o, int i, fits *fit, rectangle *area, int threads) {
	struct astrometry_data *aargs = (struct astrometry_data *)arg->user;
	aargs = copy_astrometry_args(aargs);
	if (!aargs)
		return 1;
	aargs->fit = fit;

	char root[256];
	if (!fit_sequence_get_image_filename(arg->seq, i, root, FALSE)) {
		free(aargs);
		return 1;
	}
	aargs->filename = g_strdup(root);	// for localasnet
	// process_plate_solver_input(aargs);	// depends on aargs->fit
	int retval = GPOINTER_TO_INT(plate_solver(aargs));

	if (retval)
		siril_log_color_message(_("Image %s did not solve\n"), "salmon", root);
	return retval;
}

static int astrometry_finalize_hook(struct generic_seq_args *arg) {
	struct astrometry_data *aargs = (struct astrometry_data *)arg->user;
	seq_finalize_hook(arg);
	if (aargs->cat_center)
		siril_world_cs_unref(aargs->cat_center);
	if (aargs->ref_stars)
		siril_catalog_free(aargs->ref_stars);
	free (aargs);
	if (com.child_is_running == EXT_ASNET && g_unlink("stop"))
		siril_debug_print("g_unlink() failed\n");
	com.child_is_running = EXT_NONE;
	return 0;
}

void start_sequence_astrometry(sequence *seq, struct astrometry_data *args) {
	struct generic_seq_args *seqargs = create_default_seqargs(seq);
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = seq->selnum;
	seqargs->stop_on_error = FALSE;
#ifdef _WIN32
	seqargs->parallel = args->solver != SOLVER_LOCALASNET;		// for now crashes on Cancel if parallel is enabled for asnet on windows
#else
	seqargs->parallel = TRUE;
#endif
	seqargs->prepare_hook = astrometry_prepare_hook;
	seqargs->image_hook = astrometry_image_hook;
	seqargs->finalize_hook = astrometry_finalize_hook;
	seqargs->has_output = TRUE;
	seqargs->output_type = get_data_type(seq->bitpix);
	seqargs->new_seq_prefix = strdup("ps_");
	seqargs->load_new_sequence = TRUE;
	seqargs->description = "plate solving";
	if (seq->type == SEQ_SER)
		seqargs->force_fitseq_output = TRUE;
	seqargs->user = args;

	siril_log_message(_("Running sequence plate solving using the %s catalogue\n"),
			catalog_to_str(args->ref_stars->cat_index));
	start_in_new_thread(generic_sequence_worker, seqargs);
}

