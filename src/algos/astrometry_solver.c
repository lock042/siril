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
#include "core/siril_spawn.h"
#include "core/undo.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "io/annotation_catalogues.h"
#include "algos/photometry.h"
#include "algos/photometric_cc.h"
#include "algos/spcc.h"
#include "algos/siril_wcs.h"
#include "io/image_format_fits.h"
#include "io/fits_keywords.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/siril_catalogues.h"
#include "io/local_catalogues.h"
#include "io/path_parse.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "registration/matching/match.h"
#include "registration/matching/apply_match.h"
#include "registration/matching/atpmatch.h"
#include "gui/message_dialog.h"


#define DOWNSAMPLE_FACTOR 0.25
#define CONV_TOLERANCE 1E-2 // convergence tolerance in arcsec from the projection center
#define TRANS_SANITY_CHECK 0.1 // TRANS sanity check to validate the first TRANS structure
#define NB_GRID_POINTS 7 // the number of points in one direction to create the X,Y meshgrid for inverse polynomial fiiting

#define CHECK_FOR_CANCELLATION_RET if (!get_thread_run()) { args->ret = SOLVE_CANCELLED; goto clearup;}
#define CHECK_FOR_CANCELLATION if (!get_thread_run()) { ret = SOLVE_CANCELLED; goto clearup; }

#undef DEBUG		/* get some of diagnostic output */
#define ASTROMETRY_DEBUG 0

static gchar *asnet_version = NULL;

typedef struct {
	struct wcsprm *wcslib;
} solve_results;

static void debug_print_catalog_files(TRANS *trans, s_star *star_list_A, s_star *star_list_B) {
#if ASTROMETRY_DEBUG
	GFile *file = g_file_new_for_path("ABtars.csv");
	g_autoptr(GError) error = NULL;
	GOutputStream* output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE, G_FILE_CREATE_NONE, NULL, &error);
	if (error) {
		siril_debug_print("%s\n", error->message);
		return;
	}
	gchar bufferx[1024] = { 0 }, buffery[1024] = { 0 };
	// x00, x10, x01, x20, x11, x02, x30, x21, x12, x03, x40, x31, x22, x13, x04, x50, x41, x32, x23, x14, x05
	g_sprintf(bufferx, "%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n", 
	trans->x00, trans->x10, trans->x01, trans->x20, trans->x11, trans->x02, 
	trans->x30, trans->x21, trans->x12, trans->x03, trans->x40, trans->x31, trans->x22, trans->x13,
	trans->x04, trans->x50, trans->x41, trans->x32, trans->x23, trans->x14, trans->x05);
	g_sprintf(buffery, "%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n",
	trans->y00, trans->y10, trans->y01, trans->y20, trans->y11, trans->y02,
	trans->y30, trans->y21, trans->y12, trans->y03, trans->y40, trans->y31, trans->y22, trans->y13,
	trans->y04, trans->y50, trans->y41, trans->y32, trans->y23, trans->y14, trans->y05);
	g_output_stream_write_all(output_stream, bufferx, strlen(bufferx), NULL, NULL,NULL);
	g_output_stream_write_all(output_stream, buffery, strlen(buffery), NULL, NULL,NULL);
	const gchar *header = "xA,yA,magA,indA,xB,yB,magB,indB\n";
	g_output_stream_write_all(output_stream, header, strlen(header), NULL, NULL,NULL);
	gchar buffer[256] = { 0 };
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
	struct astrometry_data *ret = calloc(1, sizeof(struct astrometry_data));
	if (!ret) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	memcpy(ret, args, sizeof(struct astrometry_data));
	if (args->cat_center)
		ret->cat_center = siril_world_cs_copy(args->cat_center);
	if (args->ref_stars) {
		ret->ref_stars = calloc(1, sizeof(siril_catalogue));
		siril_catalogue_copy(args->ref_stars, ret->ref_stars, args->nocache); // if nocache, we only copy metadata
	}
	ret->fit = NULL;
	ret->filename = NULL;
	ret->distofilename = NULL;
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

/* computes the limiting magnitude to have approx Nstars within fov
   this accounts for the position in the sky (converted to galactic coordinates)
   as star density is much scarcer when close to the galactic pole
*/
double compute_mag_limit_from_position_and_fov(double ra, double dec, double fov_degrees, int Nstars) {
	// Formulas as per https://siril-contrib-doc.readthedocs.io/en/latest/UsefulStuff.html#limit-magnitude
	// computing galactic coordinates
	double l0 = 122.9320 * DEGTORAD;
	double a0 = 192.8595 * DEGTORAD;
	double d0 =  27.1284 * DEGTORAD;
	ra *= DEGTORAD;
	dec *= DEGTORAD;
	double ml = (l0 - atan2(cos(dec) * sin(ra - a0), sin(dec) * cos(d0) - cos(dec) * sin(d0) * cos(ra - a0))) * RADTODEG;
	double mb = asin(sin(dec) * sin(d0) + cos(dec) * cos(d0) * cos(ra - a0)) * RADTODEG;
	if (ml > 180.)
		ml -= 360;
	// fov area in deg^2
	double S = 2 * (1 - cos(0.5 * fov_degrees * DEGTORAD)) * 180. * 180. / M_PI;
	// mag intercept
	double m0 = 11.68 + 2.66 * sin(fabs(mb) * DEGTORAD);
	// mag slope
	double a = 2.36 + (fabs(ml) - 90) * 0.0073 * (fabs(ml) < 90.);
	double b = 0.88 - (fabs(ml) - 90) * 0.0065 * (fabs(ml) < 90.);
	double s = a + b * sin(fabs(mb) * DEGTORAD);
	// limit mag
	double limit = m0 + s * (log10f((float)Nstars / S) - 2.);
	return max(limit, 7.);
}

static void compute_limit_mag(struct astrometry_data *args) {
	g_assert(args->ref_stars != NULL);
	if (args->mag_mode == LIMIT_MAG_ABSOLUTE)
		args->ref_stars->limitmag = args->magnitude_arg;
	else {
		// compute limit mag to have approx BRIGHTEST_STARS (i.e. 500) stars in the fov
		// same as the number of stars detected in an image
		args->ref_stars->limitmag = compute_mag_limit_from_position_and_fov(
		siril_world_cs_get_alpha(args->cat_center),
		siril_world_cs_get_delta(args->cat_center),
		args->used_fov / 60.0, BRIGHTEST_STARS);
		if (args->mag_mode == LIMIT_MAG_AUTO_WITH_OFFSET)
			args->ref_stars->limitmag += args->magnitude_arg;
	}
	siril_debug_print("using limit magnitude %f\n", args->ref_stars->limitmag);
}

static gboolean solve_is_near(struct astrometry_data *args) {
	return args->solver == SOLVER_LOCALASNET || (args->ref_stars && args->ref_stars->cat_index == CAT_LOCAL && args->searchradius > 0);
}

SirilWorldCS *get_eqs_from_header(fits *fit) {
	if (fit->keywords.wcsdata.objctra[0] != '\0' && fit->keywords.wcsdata.objctdec[0] != '\0')
		return siril_world_cs_new_from_objct_ra_dec(fit->keywords.wcsdata.objctra, fit->keywords.wcsdata.objctdec);

	else if (has_wcs(fit) && (fit->keywords.wcslib->crval[0] != 0.0 || fit->keywords.wcslib->crval[1] != 0.0))
		return siril_world_cs_new_from_a_d(fit->keywords.wcslib->crval[0], fit->keywords.wcslib->crval[1]);

	else if (fit->keywords.wcsdata.ra > DEFAULT_DOUBLE_VALUE || fit->keywords.wcsdata.dec > DEFAULT_DOUBLE_VALUE)
		return siril_world_cs_new_from_a_d(fit->keywords.wcsdata.ra, fit->keywords.wcsdata.dec);
	return NULL;
}

static void update_wcsdata_after_ps(struct astrometry_data *args) {
	if (has_wcsdata(args->fit))
		reset_wcsdata(args->fit);
	args->fit->keywords.wcsdata.pltsolvd = TRUE;
	if (args->solver == SOLVER_LOCALASNET) {
		snprintf(args->fit->keywords.wcsdata.pltsolvd_comment, FLEN_COMMENT, "Solved by Astrometry.net (%s)", asnet_version);
	} else {
		g_snprintf(args->fit->keywords.wcsdata.pltsolvd_comment, FLEN_COMMENT, "Siril internal solver");
	}
	update_wcsdata_from_wcs(args->fit);
}

static gboolean check_affine_TRANS_sanity(TRANS *trans) {
	double var1 = 0., var2 = 0.;
	var1 = fabs(trans->x10) - fabs(trans->y01);
	var2 = fabs(trans->y10) - fabs(trans->x01);
	siril_debug_print("abs(diff_cos)=%f et abs(diff_sin)=%f\n", var1, var2);
	return (fabs(var1) < TRANS_SANITY_CHECK && fabs(var2) < TRANS_SANITY_CHECK);
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
	double A[6][6] = {{ 0. }}, B[6][6] = {{ 0. }}, AP[6][6]  = {{ 0. }}, BP[6][6]  = {{ 0. }}; // we deal with images up to order 4
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

	if (trans->order >= AT_TRANS_QUINTIC) {
		Mvdecomp(cd_inv, trans->x50, trans->y50, &A[5][0], &B[5][0]);
		Mvdecomp(cd_inv, trans->x41, trans->y41, &A[4][1], &B[4][1]);
		Mvdecomp(cd_inv, trans->x32, trans->y32, &A[3][2], &B[3][2]);
		Mvdecomp(cd_inv, trans->x23, trans->y23, &A[2][3], &B[2][3]);
		Mvdecomp(cd_inv, trans->x14, trans->y14, &A[1][4], &B[1][4]);
		Mvdecomp(cd_inv, trans->x05, trans->y05, &A[0][5], &B[0][5]);

		AP[5][0] = revtrans.x50;
		AP[4][1] = revtrans.x41;
		AP[3][2] = revtrans.x32;
		AP[2][3] = revtrans.x23;
		AP[1][4] = revtrans.x14;
		AP[0][5] = revtrans.x05;

		BP[5][0] = revtrans.y50;
		BP[4][1] = revtrans.y41;
		BP[3][2] = revtrans.y32;
		BP[2][3] = revtrans.y23;
		BP[1][4] = revtrans.y14;
		BP[0][5] = revtrans.y05;
	}

	// We can now fill the disprm structure and assign it to wcslib->lin
	struct disprm *dis = calloc(1, sizeof(struct disprm));
	int ipx = 0;
	int dpmax = 10 + 2 * (N + 1) * (N + 2);
	dis->flag = -1;
	// thread-safe version of lindis as per WCSLIB documentation.
	// Anyway, dpmax should not change if we solve multiple images in parralel
	// as the trans order should be the same for all images of the sequence
	lindist(1, &wcslib->lin, dis, dpmax);
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
	dis->flag = 0;
	disset(dis);
	wcslib->lin.flag = 0;
	linset(&wcslib->lin);
	free_stars(&uvgrid);
	free_stars(&xygrid);
	return 0;
}

void flip_bottom_up_astrometry_data(fits *fit) {
	Homography H = { 0 };
	cvGetEye(&H);
	H.h11 = -1.;
	H.h12 = (double)fit->ry;
	reframe_astrometry_data(fit, &H);
}

gchar *platesolve_msg(struct astrometry_data *args) {
	const gchar *solvername = (args->solver == SOLVER_SIRIL) ? "Siril" : "Local Astrometry.net";
	switch (args->ret) {
	// warnings
		case SOLVE_LINONLY:
			return g_strdup_printf(_("%s could not find distortion polynomials for the order specified (%d) and returned a linear solution, try with a lower order or increase magnitude\n"),
			solvername, args->trans_order);
		case SOLVE_LASTSOLVE:
			return g_strdup(_("Siril solver could not find a final solution at the updated center, solution may be inaccurate\n"));
	// generic success
		case SOLVE_OK:
			return g_strdup_printf(_("%s solve succeeded.\n"), solvername);
	// generic platesolve common errors
		case SOLVE_NO_MATCH:
			return g_strdup(_("The image could not be aligned with the reference stars\n"));
		case SOLVE_CANCELLED:
			return g_strdup(_("Cancelled\n"));
		case SOLVE_DOWNSAMPLE:
			return g_strdup(_("Not enough memory\n"));
		case SOLVE_NOTENOUGHSTARS:
			return g_strdup_printf(_("There are not enough stars picked in the image. least %d are needed\n"), AT_MATCH_STARTN_LINEAR);
	// siril solver
		case SOLVE_INVALID_TRANS:
			return g_strdup(_("Transformation matrix is invalid, solve failed\n"));
		case  SOLVE_PROJ:
			return g_strdup(_("Reprojecting catalog failed\n"));
		case SOLVE_TRANS_UPDATE:
			return g_strdup(_("Updating trans failed\n"));
		case SOLVE_NEAR_NO_MATCH:
			return g_strdup_printf(_("Siril near solver could not find a solution over the search radius (%.1f deg)\n"), args->searchradius);
		case SOLVE_NEAR_TIMEOUT:
			return g_strdup_printf(_("Siril near solver time-out after %ds\n"), com.pref.astrometry.max_seconds_run);
		case SOLVE_NEAR_THREADS:
			return g_strdup(_("Siril near solver encountered a problem when allocating threads\n"));
		case SOLVE_ASNET_PROC:
			return g_strdup(_("Local Astrometry.net encountered a problem while launching its process\n"));
		default:
			return NULL;
	}
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

	if (90. - fabs(args->fit->keywords.wcsdata.dec) < 2.78e-3) // center is less than 10" off from a pole
		siril_log_message(_("Up position wrt. N is undetermined (too close to a Pole)\n"));
	else {
		// We move 10" to the North and we'll figure out the angle from there....
		// For some unknown reason, asnet may return a solution with the reference point not at center
		// We need to handle that case by passing ra and dec from args->fit->wcsdata which has been updated to take it into account
		double xN, yN, rotation;
		int status = wcs2pix(args->fit, args->fit->keywords.wcsdata.ra, args->fit->keywords.wcsdata.dec + 2.78e-3, &xN, &yN);
		xN -= args->fit->rx * 0.5;
		yN -= args->fit->ry * 0.5;
		if (!status) {
			rotation = atan2(xN, yN) * RADTODEG; // we measure clockwise wrt. +y axis
			if (image_is_flipped_from_wcs(args->fit->keywords.wcslib)) {
				if (args->flip_image) {
					rotation = 180.0 - rotation;
				} else {
					report_flip = TRUE; // we only report a flip if the image is not flipped afterwards
				}
			}
			if (rotation < 0.0)
				rotation += 360.0;
			siril_log_message(_("Up is %+.2lf deg CounterclockWise wrt. N%s\n"), rotation, report_flip ? _(" (flipped)") : "");
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
	siril_log_message(_("Image center: alpha: %s, delta: %s\n"), args->fit->keywords.wcsdata.objctra, args->fit->keywords.wcsdata.objctdec);
	if (args->cat_center) { // not true for asnet blind solve
		double dist = compute_coords_distance(args->fit->keywords.wcsdata.ra, args->fit->keywords.wcsdata.dec,
											siril_world_cs_get_alpha(args->cat_center), siril_world_cs_get_delta(args->cat_center));
		siril_log_message(_("Was %.2f arcmin from initial value\n"), dist * 60.);
	}
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

void reframe_wcs(struct wcsprm *wcslib, Homography *H) {
	if (!wcslib)
		return;
	double pc1_1, pc1_2, pc2_1, pc2_2;
	point refpointout;
	pc1_1 = H->h00 * wcslib->pc[0] + H->h01 * wcslib->pc[1];
	pc1_2 = H->h10 * wcslib->pc[0] + H->h11 * wcslib->pc[1];
	pc2_1 = H->h00 * wcslib->pc[2] + H->h01 * wcslib->pc[3];
	pc2_2 = H->h10 * wcslib->pc[2] + H->h11 * wcslib->pc[3];
	// we go back to cd formulation just to separate back again cdelt and pc
	double cd[2][2], pc[2][2];
	double scale_2 = fabs(H->h00 * H->h11 - H->h10 * H->h01);
	pc[0][0] = pc1_1;
	pc[0][1] = pc1_2;
	pc[1][0] = pc2_1;
	pc[1][1] = pc2_2;
	wcslib->cdelt[0] /= scale_2;
	wcslib->cdelt[1] /= scale_2;
	// we recombine pc and cdelt, and decompose it
	wcs_pc_to_cd(pc, wcslib->cdelt, cd);
	wcs_decompose_cd(wcslib, cd);

	// we fetch the refpoint in siril convention
	point refpointin = {wcslib->crpix[0] - 0.5, wcslib->crpix[1] - 0.5};
	cvTransformImageRefPoint(*H, refpointin, &refpointout);
	// and convert it back to FITS/WCS convention
	wcslib->crpix[0] = refpointout.x + 0.5;
	wcslib->crpix[1] = refpointout.y + 0.5;

	// and we update all the wcslib structures
	if (wcslib->lin.dispre) {
		transform_disto_coeff(wcslib->lin.dispre, H);
		struct disprm *dis = wcslib->lin.dispre;
		for (int n = 0; n < dis->ndp; n++) { // update the OFFSET keywords to new CRPIX values
			if (g_str_has_prefix(dis->dp[n].field + 4, "OFFSET.1"))
				dis->dp[n].value.f = wcslib->crpix[0];
			else if (g_str_has_prefix(dis->dp[n].field + 4, "OFFSET.2"))
				dis->dp[n].value.f = wcslib->crpix[1];
		}
		dis->flag = 0; // to update the structure
		disset(dis);
	}
	wcslib->lin.flag = 0; // to update the structure
	linset(&wcslib->lin);
	wcslib->flag = 0; // to update the structure
	wcsset(wcslib);
	wcs_print(wcslib);
	print_updated_wcs(wcslib);
}

void update_wcsdata_from_wcs(fits *fit) {
	double rac, decc;
	center2wcs(fit, &rac, &decc);
	if (rac != -1) {
		fit->keywords.wcsdata.ra = rac;
		fit->keywords.wcsdata.dec = decc;
		gchar *ra = siril_world_cs_alpha_format_from_double(rac, "%02d %02d %.3lf");
		gchar *dec = siril_world_cs_delta_format_from_double(decc, "%c%02d %02d %.3lf");
		g_sprintf(fit->keywords.wcsdata.objctra, "%s", ra);
		g_sprintf(fit->keywords.wcsdata.objctdec, "%s", dec);
		g_free(ra);
		g_free(dec);
	}
}

void reframe_astrometry_data(fits *fit, Homography *H) {
	if (!fit)
		return;
	reframe_wcs(fit->keywords.wcslib, H);
	// Update the center position in fit->wcsdata
	update_wcsdata_from_wcs(fit);
}

void wcs_pc_to_cd(double pc[][2], const double cdelt[2], double cd[][2]) {
	cd[0][0] = pc[0][0] * cdelt[0];
	cd[0][1] = pc[0][1] * cdelt[0];
	cd[1][0] = pc[1][0] * cdelt[1];
	cd[1][1] = pc[1][1] * cdelt[1];
}

static int siril_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution);
static int siril_near_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution);
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
	gboolean asnet_running = FALSE;

	args->ret = SOLVE_OK;
	solve_results solution = { 0 }; // used in the clean-up, init at the beginning

	if (args->verbose) {
		if (args->solver == SOLVER_LOCALASNET) {
			siril_log_message(_("Plate solving image with astrometry.net for a field of view of %.2f degrees\n"), args->used_fov / 60.0);
		} else {
			const char *catstr = args->ref_stars->cat_index == CAT_LOCAL ? _("local catalogues") :
				(args->ref_stars->cat_index == CAT_LOCAL_GAIA_ASTRO || args->ref_stars->cat_index == CAT_LOCAL_GAIA_XPSAMP) ? _("local Gaia catalogue") :
					_("an online catalogue");
			siril_log_message(_("Plate solving image from %s for a field of view of %.2f"
						" degrees%s, using a limit magnitude of %.2f\n"),
					catstr,
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
				args->ret = SOLVE_DOWNSAMPLE;
				goto clearup;
			}
			memcpy(&fit_backup, args->fit, sizeof(fits));
			memcpy(args->fit, &tmp, sizeof(fits));

			// TODO: should we average x and y or even better separate scales on x and y?
			args->scalefactor = (double)fit_backup.rx / (double)args->fit->rx;
			detection_layer = 0;
		}
		fits *green_fit = NULL;
		image im =  (image){ .fit = NULL, .from_seq = NULL, .index_in_seq = -1 };
		if (args->fit->keywords.bayer_pattern[0] != '\0') {
			green_fit = calloc(1, sizeof(fits));
			copyfits(args->fit, green_fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
			interpolate_nongreen(green_fit);
			im.fit = green_fit;
		} else {
			im.fit = args->fit;
		}

		stars = peaker(&im, detection_layer, &com.pref.starfinder_conf, &nb_stars,
				&(args->solvearea), FALSE, TRUE,
				BRIGHTEST_STARS, com.pref.starfinder_conf.profile, args->numthreads);

		clearfits(green_fit);
		free(green_fit);
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
		args->ret = SOLVE_NOTENOUGHSTARS;
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
			asnet_running = TRUE;
			if (g_unlink("stop")) // make sure the flag file for cancel is not already in the folder
				siril_debug_print("g_unlink() failed\n");
		}
		args->ret = local_asnet_platesolve(stars, nb_stars, args, &solution);
	} else {
		double x0 = args->fit->rx * 0.5;
		double y0 = args->fit->ry * 0.5;
		for (int s = 0; s < nb_stars; s++) {
			stars[s]->xpos -= x0;
			stars[s]->ypos = y0 - stars[s]->ypos;
		}
		if ((args->ret = siril_platesolve(stars, nb_stars, args, &solution)) > 0) {
			if (args->near_solve) {
				if (args->verbose) {
					siril_log_color_message(_("Initial solve failed\n"), "salmon", "salmon");
					siril_log_color_message(_("Attempting a near solve with a radius of %3.1f degrees\n"), "salmon", args->searchradius);
				}
				args->ret = siril_near_platesolve(stars, nb_stars, args, &solution);
			}
		}
	}
	if (args->ret > 0) // positive error codes indicate failed process
		goto clearup;

	/* 4. Print and store some results in args->fit*/
	// saving state for undo before modifying fit structure
	if (!com.script) {
		const char *undo_str = _("Plate Solve");
		undo_save_state(args->fit, undo_str);
	}
	// we copy the solution and update wcsdata
	if (has_wcs(args->fit))
		wcsfree(args->fit->keywords.wcslib);
	args->fit->keywords.wcslib = solution.wcslib;
	print_updated_wcs(args->fit->keywords.wcslib);
	update_wcsdata_after_ps(args);
	if (args->verbose)
		print_platesolving_results_from_wcs(args);
	double resolution = get_wcs_image_resolution(args->fit) * 3600.0;
	double focal_length = RADCONV * args->pixel_size / resolution;
	args->fit->keywords.focal_length = focal_length;
	args->fit->keywords.pixel_size_x = args->pixel_size;
	args->fit->keywords.pixel_size_y = args->pixel_size;
	if (com.pref.binning_update && args->fit->keywords.binning_x > 1) {
		args->fit->keywords.pixel_size_x /= args->fit->keywords.binning_x;
		args->fit->keywords.pixel_size_y /= args->fit->keywords.binning_x;
	}
	args->fit->pixelkey = TRUE;
	args->fit->focalkey = TRUE;
	if (!args->for_sequence && com.pref.astrometry.update_default_scale) {
		com.pref.starfinder_conf.focal_length = focal_length;
		com.pref.starfinder_conf.pixel_size_x = args->pixel_size;
		siril_log_message(_("Saved focal length %.2f and pixel size %.2f as default values\n"), focal_length, args->pixel_size);
	}

	/* 5.a Save distortion master
		We need to do it before flipping
	*/
	if (args->distofilename) {
		if (args->fit->keywords.wcslib->lin.dispre) {
			int mstatus = 0;
			gchar *mastername = path_parse(args->fit, args->distofilename, PATHPARSE_MODE_WRITE, &mstatus);
			if (mstatus) {
				siril_log_color_message(_("Could not save distortion master file, skipping\n"), "salmon");
			}
			if (!mstatus && save_wcs_fits(args->fit, mastername)) {
				siril_log_color_message(_("Could not save distortion master file, skipping\n"), "salmon");
			}
			g_free(mastername);
		} else {
			siril_log_color_message(_("Solution has no distortion\n"), "salmon");
			siril_log_color_message(_("Could not save distortion master file, skipping\n"), "salmon");
		}
		g_free(args->distofilename);
	}

	/* 5. Flip image if needed */
	if (args->flip_image && image_is_flipped_from_wcs(args->fit->keywords.wcslib)) {
		if (args->verbose)
			siril_log_color_message(_("Flipping image and updating astrometry data.\n"), "salmon");
		fits_flip_top_to_bottom(args->fit);
		flip_bottom_up_astrometry_data(args->fit);
		update_wcsdata_after_ps(args);
		args->image_flipped = TRUE;
	}

	/* 6. updating header */
	if (!args->for_sequence) {
		update_fits_header(args->fit);
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
	if (args->verbose && com.script) {
		gchar *msg = platesolve_msg(args);
		if (args->ret < 0)
			siril_log_color_message(_("Plate solving warning: %s"), "salmon", msg);
		else if (args->ret > 0)
			siril_log_color_message(_("Plate solving failed: %s"), "red", msg);
		else
			siril_log_color_message("%s", "green", msg);
		g_free(msg);
	}
	int ret = args->ret;
	gboolean is_verbose = args->verbose;
	if (!args->for_sequence) {
		if (asnet_running && g_unlink("stop"))
			siril_debug_print("g_unlink() failed\n");
		siril_add_idle(end_plate_solver, args);
	}
	else {
		free(args);
	}
	if (is_verbose)
		set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	return GINT_TO_POINTER(ret > 0);
}

static void nearsolve_pool_worker(gpointer data, gpointer user_data) {
	near_solve_data *nswdata = (near_solve_data *)data;
	int n = g_atomic_int_add(&nswdata->progress, 1) + 1; // g_atomic_int_add returns the atomic before the add
	if (!get_thread_run() || g_atomic_int_get(&nswdata->solved) || n == nswdata->N) {
		return;
	}
	siril_catalogue *siril_cat = calloc(1, sizeof(siril_catalogue));
	siril_catalogue_copy(nswdata->search_cat, siril_cat, TRUE); // this only copies the query parameters
	siril_cat->center_ra = nswdata->centers[n].x;
	siril_cat->center_dec = nswdata->centers[n].y;
	siril_cat->radius = nswdata->radius;
	int nbfound = siril_catalog_inner_conesearch(nswdata->search_cat, siril_cat);
	int ret;
	TRANS t = { 0 };
	double ra = -1, dec = -1;
	if (nbfound) {
		ret = match_catalog(nswdata->stars, nswdata->nb_stars, siril_cat, nswdata->scale, AT_TRANS_LINEAR, &t, &ra, &dec);
		if (ret == SOLVE_OK) {
			g_atomic_int_inc(&nswdata->solved);
			// we assign the center as a point so that if more than one worker
			// succeed, we are sure we have a consistent (ra,dec) pair
			point *center = calloc(1, sizeof(point));
			*center = (point){ra, dec};
			g_atomic_pointer_set(&(nswdata->center), center);
		}
	}
	siril_catalog_free(siril_cat);
	if (nswdata->verbose) {
		double percent = (double)n / (double)nswdata->N;
		set_progress_bar_data(NULL, percent);
	}
}

// returns a list of center points to be searched
static point *get_centers(double fov_deg, double search_radius_deg, double ra0, double dec0, int *nb) {
	point *centers;
	int N;
	double radius = fov_deg * 0.5 * DEGTORAD;
	int n = (int)ceil(2. * search_radius_deg / fov_deg) - 1; // number of dec rings
	double *d = malloc(n * sizeof(double));
	int *nl = malloc(n * sizeof(int));
	N = 0;
	for (int i = 1; i <= n; i++) {  // we don't search the initial point which has already failed
		d[i - 1] = M_PI_2 - i * radius;
		nl[i - 1] = (int)ceil(2. * M_PI * sin(i * radius) / radius); // number of points of ith ring
		N += nl[i - 1];
	}
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

static int siril_near_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution) {
	if (args->verbose) {
		set_progress_bar_data(_("Near solver started"), PROGRESS_RESET);
	}
	point *centers;
	int N, n = 0;
	int ret = SOLVE_NEAR_NO_MATCH;
	double ra = -1., dec = -1.;
	double ra0 = siril_world_cs_get_alpha(args->cat_center);
	double dec0 = siril_world_cs_get_delta(args->cat_center);
	// We prepare the larger catalogue containing the search cone
	// we limit the magnitude as we don't want as many stars (50 instead of 500 per fov)
	// we just want to have a good enough linear solution to resend a normal solve afterwards
	siril_catalogue *siril_search_cat = calloc(1, sizeof(siril_catalogue));
	siril_catalogue_copy(args->ref_stars, siril_search_cat, TRUE);
	siril_search_cat->radius = args->searchradius * 60.;
	if (args->mag_mode == LIMIT_MAG_ABSOLUTE) {
		double mag1 = compute_mag_limit_from_position_and_fov(ra0, dec0, args->used_fov / 60., BRIGHTEST_STARS / 10);
		double mag2 = compute_mag_limit_from_position_and_fov(ra0, dec0, args->used_fov / 60., BRIGHTEST_STARS);
		siril_search_cat->limitmag = args->magnitude_arg + mag1 - mag2;
		
	} else {
		siril_search_cat->limitmag = compute_mag_limit_from_position_and_fov(ra0, dec0, args->used_fov / 60., BRIGHTEST_STARS / 10);	
	}
	// we fetch all the stars within the search cone
	get_catalog_stars(siril_search_cat);

	// We prepare the list of centers to be searched
	centers = get_centers(args->used_fov / 60., args->searchradius, ra0, dec0, &N);
	// and we loop with a pool of workers...
	near_solve_data nsdata = { 0 };
	nsdata.search_cat = siril_search_cat;
	nsdata.centers = centers;
	nsdata.N = N;
	nsdata.progress = -1;
	nsdata.stars = stars;
	nsdata.nb_stars = nb_stars;
	nsdata.radius = args->ref_stars->radius;
	nsdata.scale = args->scale;
	nsdata.verbose = args->verbose;
	nsdata.center = NULL;
	gboolean immediate = FALSE;
	GThreadPool *pool = g_thread_pool_new(nearsolve_pool_worker, &nsdata, args->numthreads, FALSE, NULL);
	do {
		if (!g_thread_pool_push(pool, &nsdata, NULL)) {
			immediate = TRUE;
			break;
		}
		n++;
	} while (n < N);
	// we will wait till:
	// - we exhaust the list or
	// - one of the workers is successful or
	// - we hit timeout
	// except if we have a problem pushing new jobs, in which case we stop immediately
	if (immediate) {
		g_thread_pool_free(pool, TRUE, TRUE);
		ret = SOLVE_NEAR_THREADS;
	} else if (com.pref.astrometry.max_seconds_run > 0) {// we have a time-out specified
		guint64 timer = 0;
		guint64 timeout = com.pref.astrometry.max_seconds_run * G_TIME_SPAN_SECOND;
		while (get_thread_run() && !nsdata.solved && g_thread_pool_unprocessed(pool) > 0 && timer < timeout) {
			g_usleep(G_TIME_SPAN_SECOND);
			timer += G_TIME_SPAN_SECOND;
		}
		if (timer > timeout) {
			ret = SOLVE_NEAR_TIMEOUT;
		}
		if (!com.run_thread) {
			ret = SOLVE_CANCELLED;
		}
		g_thread_pool_free(pool, TRUE, TRUE);
	} else {
		g_thread_pool_free(pool, FALSE, TRUE);
	}

	if (nsdata.solved) {
		ret = SOLVE_OK;
		ra = nsdata.center->x;
		dec = nsdata.center->y;
	}
	if (args->verbose) {
		set_progress_bar_data(_("Near solver done"), PROGRESS_DONE);
	}
	free(centers);
	free(nsdata.center);
	siril_catalog_free(siril_search_cat);
	if (ret != SOLVE_OK)
		return ret; // the near search has failed, we stop there
	// Otherwise, we do a normal platesolve with the updated ra and dec
	if (args->verbose)
		siril_log_color_message(_("Solving again at updated center\n"), "green");
	siril_catalog_free_items(args->ref_stars);
	args->ref_stars->center_ra = ra;
	args->ref_stars->center_dec = dec;
	return siril_platesolve(stars, nb_stars, args, solution);
}

/* entry point for siril's plate solver based on catalogue matching
   args->fit is not modified by this function
   Only solution is filled
*/
static int siril_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution) {
	if (!args->ref_stars->cat_items)
		get_catalog_stars(args->ref_stars);
	TRANS t = { 0 };
	int ret = SOLVE_NO_MATCH;
	double ra = -1., dec = -1.;
	double ra0 = args->ref_stars->center_ra;
	double dec0 = args->ref_stars->center_dec;
	ret = match_catalog(stars, nb_stars, args->ref_stars, args->scale, args->trans_order, &t, &ra, &dec);
	if (ret <= 0) { // we update the solution - but if near_solve, we do a last solve with a new catalogue fetched at the center if it was too far away
		double dist = compute_coords_distance(ra0, dec0, ra, dec) * 60.; // distance from last fetched catalogue and solution center in arcmin
		if (!ret && args->near_solve && dist > 0.15 * args->used_fov) {
			siril_catalog_free_items(args->ref_stars);
			args->ref_stars->center_ra = ra;
			args->ref_stars->center_dec = dec;
			get_catalog_stars(args->ref_stars);
			TRANS t2 = { 0 };
			double ra2 = -1., dec2 = -1.;
			int ret2 = match_catalog(stars, nb_stars, args->ref_stars, args->scale, args->trans_order, &t2, &ra2, &dec2);
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
		siril_debug_print("order: %d, nmatched: %d, sigx:%g, sigy:%g\n", t.order, t.nr, t.sx, t.sy);
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
				const char CTYPE[2][13] = { "RA---TAN-SIP", "DEC--TAN-SIP" };
				const char CUNIT[2][4] = { "deg", "deg" };
				for (int i = 0; i < NAXIS; i++) {
					strncpy(prm->cunit[i], &CUNIT[i][0], 71); // 72 char fixed buffer, keep 1 for the NULL
				}
				for (int i = 0; i < NAXIS; i++) {
					strncpy(prm->ctype[i], &CTYPE[i][0], 71); // 72 byte buffer, leave 1 byte for the NULL
				}
			}
		// we pass cd and let wcslib decompose it
		wcs_decompose_cd(prm, cd);
		if (args->trans_order > AT_TRANS_LINEAR && t.order > AT_TRANS_LINEAR && add_disto_to_wcslib(prm, &t, args->rx_solver, args->ry_solver))
			ret = SOLVE_LINONLY;
		wcs_print(prm);
		solution->wcslib = prm;
	}
	return ret;
}

static int match_catalog(psf_star **stars, int nb_stars, siril_catalogue *siril_cat, double scale, int order, TRANS *trans_out, double *ra_out, double *dec_out) {
	TRANS trans = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;
	int max_trials = 5;
	s_star *star_list_A = NULL, *star_list_B = NULL;
	psf_star **cstars = NULL;
	double ra0 = siril_cat->center_ra;
	double dec0 = siril_cat->center_dec;
	int ret = SOLVE_NO_MATCH;

	double a = 1.0 + (com.pref.astrometry.percent_scale_range / 100.0);
	double b = 1.0 - (com.pref.astrometry.percent_scale_range / 100.0);
	double scale_min = 1.0 / (scale * a);
	double scale_max = 1.0 / (scale * b);

	cstars = project_catalog_stars(siril_cat, ra0, dec0);
	if (!cstars) {
		return SOLVE_PROJ;
	}
	ret = new_star_match(stars, cstars, nb_stars, siril_cat->nbitems, nobj,
			scale_min, scale_max, NULL, &trans, TRUE,
			UNDEFINED_TRANSFORMATION, AT_TRANS_LINEAR, &star_list_A, &star_list_B);

	CHECK_FOR_CANCELLATION;

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
		siril_debug_print("Reprojecting to: alpha: %.4f, delta: %+.4f\n", ra0, dec0);

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
	debug_print_catalog_files(&trans, star_list_A, star_list_B);

	// updating solution to higher order if required
	if (order > AT_TRANS_LINEAR) {
		double ra1 = ra0;
		double dec1 = dec0;
		siril_debug_print("starting non linear match at order %d\n", order);
		TRANS newtrans = { 0 }; // we will fall back on the linear solution if not successfull at higher order
		memcpy(&newtrans, &trans, sizeof(TRANS));
		newtrans.order = order;
		int ret2 = re_star_match(stars, cstars, nb_stars, siril_cat->nbitems, &newtrans, NULL, NULL);
		if (!ret2) {
			conv = get_center_offset_from_trans(&newtrans);
			trial = 0;
			siril_debug_print("iteration %d - offset: %.3f, number of matches: %d\n", trial, conv, newtrans.nr);
			while (conv > CONV_TOLERANCE && trial < max_trials && !ret2) {
				// we get the new projection center
				apply_match(ra1, dec1, 0., 0., &newtrans, &ra1, &dec1);
				// we will reproject the catalog at the new image center
				siril_debug_print("Reprojecting to: alpha: %.4f, delta: %+.4f\n", ra1, dec1);

				free_fitted_stars(cstars);
				cstars = project_catalog_stars(siril_cat, ra1, dec1);
				if (!cstars) {
					ret2 = SOLVE_PROJ;
					break;
				}
				// we have reprojected at the center, the shift terms should very close to null
				newtrans.x00 = 0.;
				newtrans.y00 = 0.;
				// and recompute the trans structure with all the stars
				if (re_star_match(stars, cstars, nb_stars, siril_cat->nbitems, &newtrans, NULL, NULL)) {
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
			ret = SOLVE_LINONLY;
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

// Retrieves and caches asnet_version. Returns true if asnet works
static gboolean get_asnet_version(gchar* exe_name) {
	g_autofree gchar* child_stdout = NULL;
	g_autofree gchar* child_stderr = NULL;
	g_autoptr(GError) error = NULL;

	gint child_exit = 127;
#ifdef _WIN32
	gchar* argv[] = {
		exe_name,
		"-l",
		"-c",
		"solve-field --version",
		NULL
	};
#else
	gchar* argv[] = {
		exe_name,
		"--version",
		NULL
	};
#endif

	siril_spawn_host_sync(
		NULL,
		argv,
		NULL,
		G_SPAWN_SEARCH_PATH,
		NULL,
		NULL,
		&child_stdout,
		&child_stderr,
		&child_exit,
		&error
	);

	if (error) {
		// Can't start asnet
		printf("%s\n",error->message);
		return FALSE;
	} else if (!siril_spawn_check_wait_status(child_exit, &error)) {
		// Abnormal termination
		if (error && error->domain == G_SPAWN_EXIT_ERROR && error->code == child_exit && child_stderr) {
			// this case happens when version is unknown, like for ansvr (returns an error code).
			// Posterior versions (>ansvr but <0.88) do not necessarily throw an error code, hence why <0.88 is defined here and below
			// Or they do, but we can't test every single version between 0.24 and 0.87...
			printf("%s\n", error->message);
			gchar *test_unknown_option = g_strstr_len(child_stderr, -1, "unknown option");
			if (test_unknown_option) {
				asnet_version = g_strdup("<0.88");
			} else
				return FALSE;
		} else
			return FALSE;
	} else {
		gchar** chunks = g_strsplit(child_stdout, "\n", 2);
		if (strlen(chunks[1]) == 0) { // if the second chunk is void, it means there was only one line which contains the version
			asnet_version = g_strdup(chunks[0]);
		} else {
			asnet_version = g_strdup("<0.88"); // version switch was introduced from 0.88 (https://github.com/dstndstn/astrometry.net/commit/25f0b829d80c57984d404de50f9c645cb3da2223)
		}
		g_strfreev(chunks);
	}
	siril_log_message(_("Running asnet version %s\n"), asnet_version);
	return TRUE;
}

void reset_asnet_version() {
	if (asnet_version) {
		g_free(asnet_version);
		asnet_version = NULL;
	}
}

#ifdef _WIN32
static gchar *siril_get_asnet_bash(gboolean *path_in_pref) {
	*path_in_pref = FALSE;
	// searching user-defined path if any
	if (com.pref.asnet_dir && com.pref.asnet_dir[0] != '\0') {
		*path_in_pref = TRUE;
		return g_build_filename(com.pref.asnet_dir, "bin", "bash", NULL);
	}
	const gchar *localappdata = g_get_user_data_dir();
	return  g_build_filename(localappdata, "cygwin_ansvr", "bin", "bash", NULL);
}

static gchar *siril_get_asnet_shell() {
	// searching user-defined path if any
	if (com.pref.asnet_dir && com.pref.asnet_dir[0] != '\0') {
		return g_strdup(com.pref.asnet_dir);
	}
	const gchar *localappdata = g_get_user_data_dir();
	return  g_build_filename(localappdata, "cygwin_ansvr", NULL);
}

/* returns true if the command solve-field is available */
gboolean asnet_is_available() {
	gboolean path_in_pref = FALSE;
	g_autofree gchar* bash_path = siril_get_asnet_bash(&path_in_pref);
	gboolean success = get_asnet_version(bash_path);
	if (!success && path_in_pref) {
		siril_log_color_message(_("Astrometry.net path (%s) incorrectly set in preferences\n"), "red", com.pref.asnet_dir);
	}
	return success;
}

#else

static gchar *siril_get_asnet_bin() {
	if (com.pref.asnet_dir && com.pref.asnet_dir[0] != '\0') {
		return g_build_filename(com.pref.asnet_dir, "solve-field", NULL);
	}

	return g_strdup("solve-field");
}

/* returns true if the command solve-field is available */
gboolean asnet_is_available() {
	g_autofree gchar* exe_path = siril_get_asnet_bin();
	return get_asnet_version(exe_path);
}
#endif

static void child_watch_cb(GPid pid, gint status, gpointer user_data) {
	siril_debug_print("asnet exited with status %d\n", status);
	g_spawn_close_pid(pid);
	// GraXpert has exited, reset the stored pid
	remove_child_from_children(pid);
}

static int local_asnet_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution) {
	if (!args->asnet_checked) {
		if (!asnet_is_available()) {
			siril_log_color_message(_("solve-field was not found, set its path in the preferences\n"), "red");
			return SOLVE_ASNET_PROC;
		}
	}

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
		g_free(stopfile);
		return SOLVE_ASNET_PROC;
	}

#ifndef _WIN32
	gchar *asnet_path = siril_get_asnet_bin();
	g_assert(asnet_path);
#else
	gchar *asnet_shell = siril_get_asnet_shell();
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
		"-U", "none", "-S", "none", "--crpix-center",
		"-s", "FLUX", NULL };

	char low_scale[16], high_scale[16];
	if (!args->asnet_blind_res) {
		double a = 1.0 + (com.pref.astrometry.percent_scale_range / 100.0);
		double b = 1.0 - (com.pref.astrometry.percent_scale_range / 100.0);
		sprintf(low_scale, "%.3f", args->scale * b);
		sprintf(high_scale, "%.3f", args->scale * a);
		char *scale_args[] = { "-u", "arcsecperpix", "-L", low_scale, "-H", high_scale, NULL };
		append_elements_to_array(sfargs, scale_args);
	}

	if (com.pref.astrometry.max_seconds_run > 0) {
		char time_limit[16];
		sprintf(time_limit, "%d", com.pref.astrometry.max_seconds_run);
		char *timeout_args[] = { "-l", time_limit, NULL };
		append_elements_to_array(sfargs, timeout_args);
	}
	char order[12];	// referenced in sfargs, needs the same scope
	if (args->trans_order > 1) {
		sprintf(order, "%d", args->trans_order);
		char *tweak_args[] = { "-t", order, NULL };
		append_elements_to_array(sfargs, tweak_args);
	} else {
		char *tweak_args[] = { "-T", NULL };
		append_elements_to_array(sfargs, tweak_args);
	}

	char start_ra[16], start_dec[16], radius[16];
	if (!args->asnet_blind_pos) {
		sprintf(start_ra, "%f", siril_world_cs_get_alpha(args->cat_center));
		sprintf(start_dec, "%f", siril_world_cs_get_delta(args->cat_center));
		sprintf(radius, "%.1f", args->searchradius);
		char *additional_args[] = { "--ra", start_ra, "--dec", start_dec,
			"--radius", radius, NULL};
		append_elements_to_array(sfargs, additional_args);
	}

	if (!args->asnet_blind_pos && !args->asnet_blind_res)
		siril_log_message(_("Astrometry.net solving with a search field at RA: %s, Dec: %s,"
			" within a %s degrees radius for scales [%s, %s]\n"),
		start_ra, start_dec, radius, low_scale, high_scale);
	else if (!args->asnet_blind_pos && args->asnet_blind_res)
		siril_log_message(_("Astrometry.net solving with a search field at RA: %s, Dec: %s,"
			" within a %s degrees radius, blindly for scale\n"),
			start_ra, start_dec, radius);
	else if (args->asnet_blind_pos && !args->asnet_blind_res)
		siril_log_message(_("Astrometry.net solving blindly in position for scales [%s, %s]\n"),
				low_scale, high_scale);
	else
		siril_log_message(_("Astrometry.net solving blindly\n"),
				low_scale, high_scale);

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
		return SOLVE_ASNET_PROC;
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
	GPid child_pid;
	child_info *child = g_malloc(sizeof(child_info));
	siril_spawn_host_async_with_pipes(NULL,
				sfargs,
				NULL,
				G_SPAWN_LEAVE_DESCRIPTORS_OPEN | G_SPAWN_SEARCH_PATH,
				NULL,
				NULL,
				&child_pid,
				NULL,
				&child_stdout,
				NULL,
				&error
	);
	// At this point, remove the processing thread from the list of children and replace it
	// with the asnet process. This avoids tracking two children for the same task.
	// Note: the pid isn't strictly needed here as it isn't used to kill the process
	// should we need to, but it makes a handy index to search for in com.children, so
	// we record it anyway.
	remove_child_from_children((GPid)-2);
	child->childpid = child_pid;
	child->program = EXT_ASNET;
	child->name = g_strdup(_("Astrometry.net local solver"));
	child->datetime = g_date_time_new_now_local();
	com.children = g_slist_prepend(com.children, child);

	// Required in order to remove the child from com.children on exit
	g_child_watch_add(child_pid, child_watch_cb, NULL);

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
		return SOLVE_ASNET_PROC;
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
	if (!success) {
		return SOLVE_NO_MATCH;
	}

	/* get the results from the .wcs file */
	gchar *wcs_filename = replace_ext(args->filename, ".wcs");
	fits result = { 0 };
	if (read_fits_metadata_from_path_first_HDU(wcs_filename, &result)) {
		siril_log_color_message(_("Could not read the solution from solve-field (expected in file %s)\n"), "red", wcs_filename);
		return SOLVE_NO_MATCH;
	}

	solution->wcslib = result.keywords.wcslib;
	wcsset(solution->wcslib);
	result.keywords.wcslib = NULL;
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
			return SOLVE_LINONLY;
		}
	}
	// In other cases, the dis struct is empyty, we still need to warn the user
	if (args->trans_order > 1 && !solution->wcslib->lin.dispre) {
		return SOLVE_LINONLY;
	}
	return SOLVE_OK;
}

// inputs: focal length, pixel size, manual, fit, autocrop, downsample, mag_mode and mag_arg
// outputs: scale, used_fov, uncentered, solvearea, limit_mag, near_solve, catalogue
void process_plate_solver_input(struct astrometry_data *args) {
	args->scale = get_resolution(args->focal_length, args->pixel_size);
	gboolean selected = FALSE;

	rectangle croparea = { 0 };
	if (!args->manual) {
		// first checking if there is a selection or if the full field is to be used (siril only)
		if (args->solver == SOLVER_SIRIL && com.selection.w != 0 && com.selection.h != 0) {
			memcpy(&croparea, &com.selection, sizeof(rectangle));
			siril_log_color_message(_("Warning: using the current selection to detect stars\n"), "salmon");
			selected = TRUE;
		} else {
			croparea.x = 0;
			croparea.y = 0;
			croparea.w = args->fit->rx;
			croparea.h = args->fit->ry;
		}
		double fov_arcmin = get_fov_arcmin(args->scale, croparea.w, croparea.h);
		siril_debug_print("image fov for given sampling: %f arcmin\n", fov_arcmin);

		// then apply or not autocropping to 5deg (300 arcmin)
		args->used_fov = args->autocrop ? min(fov_arcmin, 300.) : fov_arcmin; // autocrop is false for asnet or siril with local cat
		double cropfactor = (args->used_fov < fov_arcmin) ? args->used_fov / fov_arcmin : 1.0;
		if (cropfactor != 1.0) {
			croparea.x += (int) ((croparea.w - croparea.w * cropfactor) / 2);
			croparea.y += (int) ((croparea.h - croparea.h * cropfactor) / 2);
			croparea.w = (int) (cropfactor * croparea.w);
			croparea.h = (int) (cropfactor * croparea.h);
			siril_debug_print("Auto-crop factor: %.2f\n", cropfactor);
		}

		if (args->solver == SOLVER_SIRIL && com.selection.w != 0 && com.selection.h != 0) {
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

	if (args->solver == SOLVER_LOCALASNET) {
		args->near_solve = TRUE;
		return;
	}
	if (args->ref_stars) {
		args->ref_stars->radius = args->used_fov * 0.5;
	} else {
		siril_debug_print("Warning: args->ref_stars is NULL\n");
		return;
	}
	if (croparea.w == args->fit->rx && croparea.h == args->fit->ry)
		memset(&croparea, 0, sizeof(rectangle));
	else siril_debug_print("reduced area for the solve: %d, %d, %d x %d%s\n",
			croparea.x, croparea.y, croparea.w, croparea.h,
			args->downsample ? " (down-sampled)" : "");
	memcpy(&(args->solvearea), &croparea, sizeof(rectangle));

	compute_limit_mag(args); // to call after having set args->used_fov

	if (args->ref_stars->cat_index == CAT_AUTO) {
		if (args->ref_stars->limitmag <= 6.5) {
			args->ref_stars->cat_index = CAT_BSC;
		} else if (args->ref_stars->limitmag <= 12.5)
			args->ref_stars->cat_index = CAT_TYCHO2;
		else if (args->ref_stars->limitmag <= 17.0 || args->used_fov > 180.) // we should probably limit the use of GAIA to smaller fov, <60 as it is very long to query
			args->ref_stars->cat_index = CAT_NOMAD;
		else args->ref_stars->cat_index = CAT_GAIADR3;
	}
	args->ref_stars->columns = siril_catalog_columns(args->ref_stars->cat_index);
	if (selected)
		args->near_solve = FALSE;
	else
		args->near_solve = solve_is_near(args); // sets flag to check if the catalogue needs to be fetched at the start or if it will be fetched on demand
}

static int astrometry_prepare_hook(struct generic_seq_args *arg) {
	struct astrometry_data *args = (struct astrometry_data *)arg->user;
	fits fit = { 0 };
	// load ref metadata in fit
	if (seq_read_frame_metadata(arg->seq, sequence_find_refimage(arg->seq), &fit))
		return 1;
	args->fit = &fit;
	process_plate_solver_input(args); // compute required data to get the catalog
	if (args->solver == SOLVER_SIRIL && !args->ref_stars) {
		siril_log_color_message(_("Error: no reference stars available\n"), "red");
		return 1;
	}
	args->layer = fit.naxes[2] == 1 ? 0 : 1;
	args->sfargs->layer = fit.naxes[2] == 1 ? 0 : 1;
	args->fit = NULL;
	clearfits(&fit);
	// prepare for astrometric registration
	if (args->update_reg) {
		if (!arg->seq->regparam[args->layer]) {
			regdata *current_regdata = calloc(arg->seq->number, sizeof(regdata));
			if (current_regdata == NULL) {
				PRINT_ALLOC_ERR;
				return 1;
			}
			arg->seq->regparam[args->layer] = current_regdata;
		} else {
			siril_log_message(_("Recomputing already existing registration for this layer\n"));
			/* we reset all values as we may register different images */
			memset(arg->seq->regparam[args->layer], 0, arg->seq->number * sizeof(regdata));
		}
		args->WCSDATA = calloc(arg->seq->number, sizeof(struct wcsprm));
		arg->seq->distoparam[args->layer].index = DISTO_FILES;
	}
	if (arg->has_output && seq_prepare_hook(arg))
		return 1; // bail if seq_prepare_hook fails
	if (args->solver == SOLVER_LOCALASNET) {
		g_unlink("stop"); // make sure the flag file for cancel is not already in the folder
	}
	if (!args->nocache)
		return get_catalog_stars(args->ref_stars);
	args->seqprogress = 0; // initialize success counter
	args->seqskipped = 0; // initialize skipped counter
	return 0;
}

static int astrometry_image_hook(struct generic_seq_args *arg, int o, int i, fits *fit, rectangle *area, int threads) {
	struct astrometry_data *aargs_master = (struct astrometry_data *)arg->user;
	if (!aargs_master->force && has_wcs(fit) && !aargs_master->update_reg) { // if we need to store regdata, we still need to detect stars
		g_atomic_int_inc(&aargs_master->seqskipped);
		g_atomic_int_inc(&aargs_master->seqprogress);
		siril_log_color_message(_("Image %d already platesolved, skipping\n"), "salmon", i + 1);
		return 0;
	}

	// we detect stars first (using aargs_master)
	int nb_stars = 0;
	psf_star **stars = NULL;

	struct starfinder_data *sf_data = findstar_image_worker(aargs_master->sfargs, -1, i, fit, area, threads);
	if (sf_data) {
		stars = *sf_data->stars;
		nb_stars = *sf_data->nb_stars;
		free(sf_data);
	}

	if (aargs_master->update_reg && nb_stars) {
		float FWHMx, FWHMy, B;
		char *units;
		FWHM_stats(stars, nb_stars, arg->seq->bitpix, &FWHMx, &FWHMy, &units, &B, NULL, 0.);
		regdata *current_regdata = arg->seq->regparam[aargs_master->layer];
		current_regdata[i].roundness = FWHMy/FWHMx;
		current_regdata[i].fwhm = FWHMx;
		current_regdata[i].weighted_fwhm = FWHMx;
		current_regdata[i].background_lvl = B;
		current_regdata[i].number_of_stars = nb_stars;
		current_regdata[i].H = (Homography){ 0 }; // we will update it in the finalize_hook
	}

	if (!nb_stars) {
		siril_log_color_message(_("Image %d: no stars found\n"), "red", i + 1);
		arg->seq->imgparam[o].incl = FALSE;
		return 1;
	}

	if (!aargs_master->force && has_wcs(fit)) { // we have filled the regparam and allow to skip, we skip now
		int status = 0;
		struct wcsprm *wcs = wcs_deepcopy(fit->keywords.wcslib, &status);
		if (status) {
			siril_log_color_message(_("Could not copy WCS data, skipping image %d\n"), "salmon", i + 1);
		} else {
			memcpy(aargs_master->WCSDATA + i, wcs, sizeof(*wcs));
			g_atomic_int_inc(&aargs_master->seqskipped);
			siril_log_color_message(_("Image %d already platesolved, skipping\n"), "salmon", i + 1);
		}
		free_fitted_stars(stars);
		g_atomic_int_inc(&aargs_master->seqprogress);
		return status;
	}

	// We prepare to platesolve and collect the catalog inputs
	struct astrometry_data *aargs = copy_astrometry_args(aargs_master);
	if (!aargs)
		return 1;
	aargs->fit = fit;
	aargs->manual = TRUE;
	aargs->stars = stars;

	char root[256];
	if (!fit_sequence_get_image_filename(arg->seq, i, root, FALSE)) {
		siril_catalog_free(aargs->ref_stars);
		if (aargs->cat_center) // not filled if asnet blind solve
			siril_world_cs_unref(aargs->cat_center);
		free(aargs);
		return 1;
	}

	if (aargs->nocache) {
		if (!aargs->forced_metadata[FORCED_CENTER] && aargs->cat_center) { // center coordinates where not forced, we need to read new ones from header or skip if blind for asnet
			siril_world_cs_unref(aargs->cat_center);
			SirilWorldCS *target_coords = get_eqs_from_header(fit);
			aargs->cat_center = target_coords;
			if (aargs->ref_stars && target_coords) {
				aargs->ref_stars->center_ra = siril_world_cs_get_alpha(target_coords);
				aargs->ref_stars->center_dec = siril_world_cs_get_delta(target_coords);
			}
		}
		if (aargs->solver == SOLVER_SIRIL && !aargs->cat_center) { // we need a cat_center for Siril
			siril_debug_print("Could not set cat_center, skipping\n");
			siril_catalog_free(aargs->ref_stars);
			free(aargs);
			return 1;
		}
		if (!aargs->forced_metadata[FORCED_PIXEL]) { // pixel size was not forced, we need to read new one from header/settings
			aargs->pixel_size = max(fit->keywords.pixel_size_x, fit->keywords.pixel_size_y);
			if (aargs->pixel_size <= 0.) {
				aargs->pixel_size = com.pref.starfinder_conf.pixel_size_x;
				if (aargs->pixel_size <= 0.0) {
					siril_log_color_message(_("Could not retrieve pixel size from image %s metadata or settings, skipping\n"), "red", root);
					siril_catalog_free(aargs->ref_stars);
					if (aargs->cat_center) // not filled if asnet blind solve
						siril_world_cs_unref(aargs->cat_center);
					free(aargs);
					return 1;
				}
			}
		}
		if (!aargs->forced_metadata[FORCED_FOCAL]) { // focal was not forced, we need to read new one from header/settings
			aargs->focal_length = fit->keywords.focal_length;
			if (aargs->focal_length <= 0.0) {
				aargs->focal_length = com.pref.starfinder_conf.focal_length;
				if (aargs->focal_length <= 0.0) {
					siril_log_color_message(_("Could not retrieve focal length from image %s metadata or settings, skipping\n"), "red", root);
					siril_catalog_free(aargs->ref_stars);
					if (aargs->cat_center) // not filled if asnet blind solve
						siril_world_cs_unref(aargs->cat_center);
					free(aargs);
					return 1;
				}
			}
		}
	}
	if (aargs->solver == SOLVER_LOCALASNET)
		aargs->filename = g_strdup(root);	// for localasnet

	if (i == arg->seq->reference_image) { // we want to save master distortion only once for ref image
		aargs->distofilename = g_strdup(aargs_master->distofilename);
	}
	process_plate_solver_input(aargs);

	int retval = GPOINTER_TO_INT(plate_solver(aargs));

	if (retval) {
		siril_log_color_message(_("Image %s did not solve\n"), "red", root);
		arg->seq->imgparam[o].incl = FALSE;
	}

	if (!retval && !arg->has_output) { // SEQ_REGULAR
		fit_sequence_get_image_filename(arg->seq, i, root, TRUE);
		int status = 0;
		// we don't want to overwrite original files, so we test for symlinks
		if (is_symlink_file(root))
			siril_debug_print("Image %s was a symlink, creating a new file to keep original untouched\n", root);
		status = savefits(root, fit);
		if (!status) {
			siril_log_color_message(_("Image %s platesolved and updated\n"), "salmon", root);
		} else {
			siril_log_color_message(_("Image %s platesolved but could not be saved\n"), "red", root);
			free(aargs);
			g_atomic_int_inc(&aargs_master->seqprogress);
			return 1;
		}
	}
	if (!retval && aargs_master->update_reg) {
		int status = 0;
		struct wcsprm *wcs = wcs_deepcopy(fit->keywords.wcslib, &status);
		if (status) {
			siril_log_color_message(_("Could not copy WCS data, skipping image %d\n"), "salmon", i + 1);
		} else {
			memcpy(aargs_master->WCSDATA + i, wcs, sizeof(*wcs));
		}
	}
	g_atomic_int_inc(&aargs_master->seqprogress);
	return retval;
}

static int astrometry_finalize_hook(struct generic_seq_args *arg) {
	struct astrometry_data *aargs = (struct astrometry_data *)arg->user;
	int retval = 0;
	siril_log_color_message(_("%d images successfully platesolved out of %d included\n"), "green", aargs->seqprogress, arg->nb_filtered_images);
	if (aargs->seqskipped > 0)
		siril_log_color_message(_("(%d were already solved and skipped)\n"), "green", aargs->seqskipped);
	if (arg->has_output)
		seq_finalize_hook(arg);
	if (aargs->update_reg && !arg->retval) {
		siril_log_color_message(_("Computing astrometric registration...\n"), "green");
		arg->retval = compute_Hs_from_astrometry(arg->seq, aargs->WCSDATA, FRAMING_CURRENT, aargs->layer, NULL, NULL);
	}
	if (!arg->retval)
		writeseqfile(arg->seq);
	if (aargs->cat_center)
		siril_world_cs_unref(aargs->cat_center);
	if (aargs->ref_stars)
		siril_catalog_free(aargs->ref_stars);
	if (aargs->distofilename)
		g_free(aargs->distofilename);
	if (aargs->solver == SOLVER_LOCALASNET && g_unlink("stop"))
		siril_debug_print("g_unlink() failed\n");
	free(aargs);
	return retval;
}

void free_astrometry_data(struct astrometry_data *args) {
	if (!args)
		return;
	if (args->cat_center)
		siril_world_cs_unref(args->cat_center);
	if (args->stars)
		free_fitted_stars(args->stars);
	if (args->ref_stars)
		siril_catalog_free(args->ref_stars);
	if (args->filename)
		g_free(args->filename);
	if (args->distofilename)
		g_free(args->distofilename);
	free(args);
}


void start_sequence_astrometry(sequence *seq, struct astrometry_data *args) {
	struct generic_seq_args *seqargs = create_default_seqargs(seq);
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = seq->selnum;
	seqargs->stop_on_error = FALSE;
#ifdef _WIN32
	seqargs->parallel = args->solver != SOLVER_LOCALASNET && seq->type != SEQ_FITSEQ;		// for now crashes on Cancel if parallel is enabled for asnet on windows
#else
	seqargs->parallel = seq->type != SEQ_FITSEQ;
#endif
	args->numthreads = (seqargs->parallel) ? 1 : com.max_thread;
	seqargs->prepare_hook = astrometry_prepare_hook;
	seqargs->image_hook = astrometry_image_hook;
	seqargs->finalize_hook = astrometry_finalize_hook;
	seqargs->idle_function = end_platesolve_sequence;
	seqargs->has_output = seq->type == SEQ_FITSEQ || seq->type == SEQ_SER; // we don't save a new sequence for sequence of fits files
	seqargs->output_type = get_data_type(seq->bitpix);
	seqargs->description = "plate solving";
	if (seqargs->has_output) {
		seqargs->force_fitseq_output = TRUE;
		seqargs->new_seq_prefix = strdup("ps_");
		seqargs->load_new_sequence = TRUE;
		args->force = TRUE;
	}
	seqargs->user = args;
	if (args->solver == SOLVER_SIRIL)
		siril_log_message(_("Running sequence plate solving using the %s catalogue\n"),
				catalog_to_str(args->ref_stars->cat_index));
	start_in_new_thread(generic_sequence_worker, seqargs);
}

