/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#include "registration/matching/project_coords.h"
#include "gui/message_dialog.h"


#define DOWNSAMPLE_FACTOR 0.25
#define CONV_TOLERANCE 1E-1 // convergence tolerance in arcsec from the projection center
#define PLATESOLVE_STEP 100. // step made on CRPIX axes to compute the CD matrix

#undef DEBUG		/* get some of diagnostic output */
#define ASTROMETRY_DEBUG 1

static gchar *asnet_version = NULL;

typedef struct {
	point size;
	SirilWorldCS *px_cat_center;	// the original target first, but can get refined
	SirilWorldCS *image_center;
	double crpix[2];
	double pixel_size;		// in µm
	double focal_length;		// in mm
	Homography H;			// for matching results printing
	gboolean image_is_flipped;
} solve_results;

static void debug_print_catalog_files(s_star *star_list_A, s_star *star_list_B) {
#ifdef ASTROMETRY_DEBUG
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
		ret->cat_center = siril_world_cs_ref(args->cat_center);
	ret->fit = NULL;
	ret->filename = NULL;
	/* assuming catalog stays the same */
	return ret;
}

static void fov_in_DHMS(double var, gchar *fov) {
	int deg, decM;
	double decS;

	if (var < 0) {
		fprintf(stdout, "fov_in_DHMS: negative value, should not happened\n");
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

gboolean has_any_keywords() {
	return (gfit.focal_length > 0.0 ||
			gfit.pixel_size_x > 0.f ||
			gfit.pixel_size_y > 0.f ||
			(has_wcs(&gfit) && gfit.wcslib->crval[0] != 0.0 && gfit.wcslib->crval[0] != 0.0) ||
			(gfit.wcsdata.objctra[0] != '\0' && gfit.wcsdata.objctdec[0] != '\0') ||
			(gfit.wcsdata.ra != 0.0 && gfit.wcsdata.dec != 0.0));
}

SirilWorldCS *get_eqs_from_header(fits *fit) {
	if (fit->wcsdata.ra != 0.0 || fit->wcsdata.dec != 0.0)
		return siril_world_cs_new_from_a_d(fit->wcsdata.ra, fit->wcsdata.dec);

	else if (fit->wcsdata.objctra[0] != '\0' && fit->wcsdata.objctdec[0] != '\0')
		return siril_world_cs_new_from_objct_ra_dec(fit->wcsdata.objctra, fit->wcsdata.objctdec);

	else if (has_wcs(fit) && (fit->wcslib->crval[0] != 0.0 || fit->wcslib->crval[1] != 0.0))
		return siril_world_cs_new_from_a_d(fit->wcslib->crval[0], fit->wcslib->crval[1]);
	return NULL;
}

static void print_platesolving_results_from_wcs(struct astrometry_data *args) {
	char field_x[256] = "";
	char field_y[256] = "";
	gboolean report_flip = FALSE;

	double cd[2][2];
	wcs_cd2mat(args->fit->wcslib, cd);
	double det = (cd[0][0] * cd[1][1] - cd[1][0] * cd[0][1]); // determinant of rotation matrix (ad - bc)

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
			if (det > 0) {
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
}

static void print_image_center(solve_results *image) {
	gchar *alpha = siril_world_cs_alpha_format(image->image_center, "%02dh%02dm%02ds");
	gchar *delta = siril_world_cs_delta_format(image->image_center, "%c%02d°%02d\'%02d\"");
	siril_log_message(_("Image center: alpha: %s, delta: %s\n"), alpha, delta);
	g_free(alpha);
	g_free(delta);
}

static gboolean check_affine_TRANS_sanity(TRANS *trans) {
	double var1, var2;
	switch (trans->order) {
		case AT_TRANS_LINEAR:
			var1 = fabs(trans->b) - fabs(trans->f);
			var2 = fabs(trans->c) - fabs(trans->e);
			break;
		case AT_TRANS_QUADRATIC:
			var1 = fabs(trans->b) - fabs(trans->i);
			var2 = fabs(trans->c) - fabs(trans->h);
			break;
		case AT_TRANS_CUBIC:
			// var1 = fabs(trans->b) - fabs(trans->k);
			// var2 = fabs(trans->c) - fabs(trans->j);
			var1 = fabs(trans->b) - fabs(trans->m);
			var2 = fabs(trans->c) - fabs(trans->l);
			break;
		default:
			break;
	}
	siril_debug_print("abs(diff_cos)=%f et abs(diff_sin)=%f\n", var1, var2);
	return ((fabs(var1) < 0.01) && fabs(var2) < 0.01);
}

static gboolean image_is_flipped_from_trans(TRANS *trans) {
	double det = 1.;
	switch (trans->order) {
		case AT_TRANS_LINEAR:
			det = (trans->b * trans->f - trans->c * trans->e);
			break;
		case AT_TRANS_QUADRATIC:
			det = (trans->b * trans->i - trans->c * trans->h);
			break;
		case AT_TRANS_CUBIC:
			// det = (trans->b * trans->k - trans->c * trans->j);
			det = (trans->b * trans->m - trans->c * trans->l);
			break;
		default:
			break;
	}
	return det > 0;
}

static double get_resolution_from_trans(TRANS *trans) {
	switch (trans->order) {
		case AT_TRANS_LINEAR:
			return sqrt(fabs(trans->b * trans->f - trans->c * trans->e));
		case AT_TRANS_QUADRATIC:
			return sqrt(fabs(trans->b * trans->i - trans->c * trans->h));
		case AT_TRANS_CUBIC:
			// return sqrt(fabs(trans->b * trans->k - trans->c * trans->j));
			return sqrt(fabs(trans->b * trans->m - trans->c * trans->l));
		default:
			return 0.;
	}
}

static void get_cdelt_from_trans(TRANS *trans, double flip, double *cdelt1, double *cdelt2) {
	*cdelt1 = sqrt(trans->b * trans->b + trans->c * trans->c) / 3600 * -flip; // cdelt1 is < 0 when the image is not flipped
	switch (trans->order) {
		case AT_TRANS_LINEAR:
			*cdelt2 = sqrt(trans->e * trans->e + trans->f * trans->f) / 3600;
			break;
		case AT_TRANS_QUADRATIC:
			*cdelt2 = sqrt(trans->h * trans->h + trans->i * trans->i) / 3600;
			break;
		case AT_TRANS_CUBIC:
			// *cdelt2 = sqrt(trans->j * trans->j + trans->k * trans->k) / 3600;
			*cdelt2 = sqrt(trans->l * trans->l + trans->m * trans->m) / 3600;
			break;
		default:
			break;
	}
}

static double get_center_offset_from_trans(TRANS *trans) {
	// this measures the remaining offset, expressed in pseudo-arcsec
	switch (trans->order) {
		case AT_TRANS_LINEAR:
			return sqrt(trans->a * trans->a + trans->d * trans->d);
		case AT_TRANS_QUADRATIC:
			return sqrt(trans->a * trans->a + trans->g * trans->g);
		case AT_TRANS_CUBIC:
			// return sqrt(trans->a * trans->a + trans->i * trans->i);
			return sqrt(trans->a * trans->a + trans->k * trans->k);
		default:
			return 0.;
	}
}

static gboolean image_is_flipped_from_wcs(fits *fit) {
	double cd[2][2];
	wcs_cd2mat(fit->wcslib, cd);
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

void print_updated_wcs_data(fits *fit) {
	/* debug output */
	siril_debug_print("****Updated WCS data*************\n");
	siril_debug_print("crpix1 = %*.12e\n", 20, fit->wcslib->crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, fit->wcslib->crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, fit->wcslib->crval[0]);
	siril_debug_print("crval2 = %*.12e\n", 20, fit->wcslib->crval[1]);
	siril_debug_print("cdelt1 = %*.12e\n", 20, fit->wcslib->cdelt[0]);
	siril_debug_print("cdelt2 = %*.12e\n", 20, fit->wcslib->cdelt[1]);
	siril_debug_print("pc1_1  = %*.12e\n", 20, fit->wcslib->pc[0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, fit->wcslib->pc[1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, fit->wcslib->pc[2]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, fit->wcslib->pc[3]);
	siril_debug_print("******************************************\n");
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
	fit->wcslib->flag = 0; // to update the structure
	wcsset(fit->wcslib);

	print_updated_wcs_data(fit);
}


void wcs_pc_to_cd(double pc[][2], const double cdelt[2], double cd[][2]) {
	cd[0][0] = pc[0][0] * cdelt[0];
	cd[0][1] = pc[0][1] * cdelt[0];
	cd[1][0] = pc[1][0] * cdelt[1];
	cd[1][1] = pc[1][1] * cdelt[1];
}

static int match_catalog(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution);
static int local_asnet_platesolve(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution);

#define CHECK_FOR_CANCELLATION_RET if (!get_thread_run()) { args->message = g_strdup(_("Cancelled")); args->ret = 1; return 1; }
static int get_catalog_stars(struct astrometry_data *args, gboolean do_fetch) {
	if (args->ref_stars->cat_index == CAT_ASNET)
		return 0;
	if (do_fetch) {
		if (siril_catalog_conesearch(args->ref_stars) <= 0)
			return 1;
		siril_log_message(_("Fetched %d stars from %s catalogue\n"), args->ref_stars->nbitems, catalog_to_str(args->ref_stars->cat_index));
	} else { // the stars were already fetched, we just reset the projection
		siril_catalog_reset_projection(args->ref_stars);
	}

	CHECK_FOR_CANCELLATION_RET;
	double ra0 = siril_world_cs_get_alpha(args->cat_center);
	double dec0 = siril_world_cs_get_delta(args->cat_center);
	GDateTime *dateobs = NULL;
	if (args->fit && args->fit->date_obs)
		dateobs = args->fit->date_obs;
	if (args->cstars) {
		free_fitted_stars(args->cstars);
		args->cstars = NULL; // next step may fail and we may try to use it again
	}
	if (!siril_catalog_project_at_center(args->ref_stars, ra0, dec0, TRUE, dateobs)) {
		args->cstars = convert_siril_cat_to_psf_stars(args->ref_stars, &args->n_cat);
		args->n_cat = args->ref_stars->nbitems;
		return 0;
	}
	siril_debug_print("Could not convert catalog to a psf_star list\n");
	return 1;
}

#define CHECK_FOR_CANCELLATION if (!get_thread_run()) { args->message = g_strdup(_("Cancelled")); args->ret = 1; goto clearup; }

/* entry point for plate solving */
gpointer plate_solver(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *) p;
	psf_star **stars = NULL;	// image stars
	int nb_stars = 0;	// number of image and catalogue stars

	args->ret = ERROR_PLATESOLVE;
	args->message = NULL;
	solve_results solution = { 0 }; // used in the clean-up, init at the beginning

	if (args->verbose) {
		if (args->ref_stars->cat_index == CAT_ASNET) {
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

	/* 1. Get catalogue stars for the field of view (for sequences, see the prepare hook) */
	if (!args->for_sequence && get_catalog_stars(args, TRUE)) {
		goto clearup;
	}
	CHECK_FOR_CANCELLATION;

	/* 2. Get image stars */
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
		// capping the detection to max usable number of stars
		if (args->n_cat == 0)
				args->n_cat = BRIGHTEST_STARS;
		int max_stars = args->for_photometry_cc ? args->n_cat : min(args->n_cat, BRIGHTEST_STARS);

#ifdef _WIN32
		// on Windows, asnet is not run in parallel neither on single image nor sequence, we can use all threads
		int nthreads = (!args->for_sequence || args->ref_stars->cat_index == CAT_ASNET) ? com.max_thread : 1;
#else
		// on UNIX, asnet is in parallel for sequences, we need to restrain to one per worker
		int nthreads = (!args->for_sequence) ? com.max_thread : 1;
#endif

		stars = peaker(&im, detection_layer, &com.pref.starfinder_conf, &nb_stars,
				&(args->solvearea), FALSE, TRUE,
				max_stars, com.pref.starfinder_conf.profile, nthreads);

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
		stars = args->stars ? args->stars : com.stars;
		if (stars)
			while (stars[nb_stars])
				nb_stars++;

	}
	CHECK_FOR_CANCELLATION;

	if (!stars || nb_stars < AT_MATCH_STARTN_LINEAR) {
		args->message = g_strdup_printf(_("There are not enough stars picked in the image. "
				"At least %d are needed."), AT_MATCH_STARTN_LINEAR);
		args->ret = ERROR_PLATESOLVE;
		goto clearup;
	}
	if (args->verbose)
		siril_log_message(_("Using %d detected stars from image.\n"), nb_stars);

	/* 3. Plate solving */
	solution.size.x = args->fit->rx;
	solution.size.y = args->fit->ry;
	solution.pixel_size = args->pixel_size;

	if (args->ref_stars->cat_index == CAT_ASNET) {
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
		if (match_catalog(stars, nb_stars, args, &solution)) {
			args->ret = ERROR_PLATESOLVE;
		}
	}
	if (args->ret)
		goto clearup;

	/* 4. Print and store some results */
	args->fit->focal_length = solution.focal_length;
	args->fit->pixel_size_x = args->fit->pixel_size_y = solution.pixel_size;
	if (!args->for_sequence && com.pref.astrometry.update_default_scale) {
		com.pref.starfinder_conf.focal_length = solution.focal_length;
		com.pref.starfinder_conf.pixel_size_x = solution.pixel_size;
		siril_log_message(_("Saved focal length %.2f and pixel size %.2f as default values\n"), solution.focal_length, solution.pixel_size);
	}
	print_image_center(&solution);

	/* 5. Run photometric color correction, if enabled */
	if (args->for_photometry_cc) {
		pcc_star *pcc_stars = NULL;
		int nb_pcc_stars;
		// We relaunch the conesearch with phot flag set to TRUE
		// For local catalogue, we fetch a whole new set at updated center and res
		// For online catalogue, we re-read the same catalogue from cache with the new phot flag
		// In both cases, we free the cat_items members before fetching
		// and we project with the platesolve wcs
		args->ref_stars->phot = TRUE;
		if (args->ref_stars->cat_index == CAT_LOCAL) {
			// we update the ref_stars structure to query again the local catalogues
			args->ref_stars->center_ra = siril_world_cs_get_alpha(solution.image_center);
			args->ref_stars->center_dec = siril_world_cs_get_delta(solution.image_center);
			double res = get_resolution(solution.focal_length, args->pixel_size);
			args->ref_stars->radius = get_radius_deg(res, args->fit->rx, args->fit->ry) * 60.;
			// for photometry, we can use fainter stars, 1.5 seems ok above instead of 2.0
			if (args->verbose)
				siril_log_message(_("Getting stars from local catalogues for PCC, limit magnitude %.2f\n"), args->ref_stars->limitmag);
		}
		siril_catalog_free_items(args->ref_stars);
		siril_catalog_conesearch(args->ref_stars);
		siril_catalog_project_with_WCS(args->ref_stars, args->fit, TRUE, FALSE);
		pcc_stars = convert_siril_cat_to_pcc_stars(args->ref_stars, &nb_pcc_stars);
		args->ret = nb_pcc_stars == 0;

		if (args->ret) {
			if (pcc_stars)
				free(pcc_stars);
			args->message = g_strdup(_("Using plate solving to identify catalogue stars in the image failed, is plate solving wrong?\n"));
			args->ret = ERROR_PHOTOMETRY;
			goto clearup;
		}
		args->pcc->stars = pcc_stars;
		args->pcc->nb_stars = nb_pcc_stars;
		args->pcc->fwhm = filtered_FWHM_average(stars, nb_stars);
		if (args->downsample)
			args->pcc->fwhm /= DOWNSAMPLE_FACTOR;

		if (photometric_cc(args->pcc)) {
			args->ret = ERROR_PHOTOMETRY;
		}

		args->pcc = NULL; // freed in PCC code
		free(pcc_stars);
		pcc_stars = NULL;
		if (args->ret) {
			args->message = g_strdup_printf(_("An astrometric solution was found but photometry analysis of the %d stars failed. This generally happens if they are saturated in the image or if they are too faint to have B-V index information (mag > 18)\n"), nb_pcc_stars);
			//goto clearup; // still flip
		} else {
			if (!args->for_sequence) {
				set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
				siril_log_color_message(_("Photometric Color Calibration succeeded.\n"), "green");
			}
		}
	}

	/* 6. Flip image if needed */
	if (args->flip_image && solution.image_is_flipped) {
		if (args->verbose)
			siril_log_color_message(_("Flipping image and updating astrometry data.\n"), "salmon");
		fits_flip_top_to_bottom(args->fit);
		flip_bottom_up_astrometry_data(args->fit);
		args->image_flipped = TRUE;
	}

	/* 7. Clean-up */
	args->new_center = solution.image_center;

clearup:
	if (stars && !args->manual) {
		for (int i = 0; i < nb_stars; i++)
			free_psf(stars[i]);
		free(stars);
	}
	if (solution.px_cat_center)
		siril_world_cs_unref(solution.px_cat_center);
	if (args->cat_center)
		siril_world_cs_unref(args->cat_center);
	if (!args->for_sequence) {
		if (args->cstars)
			free_fitted_stars(args->cstars);
	}
	g_free(args->filename);

	int retval = args->ret;
	if (com.script && retval) {
		if (retval == ERROR_PHOTOMETRY) {
			siril_log_message(_("Photometry failed: %s\n"), args->message);
		} else {
			siril_log_message(_("Plate solving failed: %s\n"), args->message);
		}
		g_free(args->message);
	}
	if (!args->for_sequence) {
		com.child_is_running = EXT_NONE;
		if (g_unlink("stop"))
			siril_debug_print("g_unlink() failed\n");
		siril_add_idle(end_plate_solver, args);
	}
	else free(args);
	return GINT_TO_POINTER(retval);
}

/* entry point for siril's plate solver based on catalogue matching */
static int match_catalog(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution) {
	TRANS trans = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;
	int max_trials = 0;
	s_star *star_list_A = NULL, *star_list_B = NULL;
	args->trans_order = AT_TRANS_CUBIC; // TODO: for debugging purposes only - to be removed

	max_trials = 20; //retry to converge if solve is done at an offset from the center

	/* make sure that arrays are not too small
	 * make sure that the max of stars is BRIGHTEST_STARS */
	int n = min(min(nb_stars, args->n_cat), BRIGHTEST_STARS);

	double a = 1.0 + (com.pref.astrometry.percent_scale_range / 100.0);
	double b = 1.0 - (com.pref.astrometry.percent_scale_range / 100.0);
	double scale_min = 1.0 / (args->scale * a);
	double scale_max = 1.0 / (args->scale * b);
	int attempt = 1;
	while (args->ret && attempt <= 3) {
		free_stars(&star_list_A);
		free_stars(&star_list_B);
		args->ret = new_star_match(stars, args->cstars, n, nobj,
				scale_min, scale_max, NULL, &trans, TRUE,
				UNDEFINED_TRANSFORMATION, args->trans_order, &star_list_A, &star_list_B);
		if (attempt == 2) {
			scale_min = -1.0;
			scale_max = -1.0;
		} else {
			nobj += 50;
		}
		attempt++;
		CHECK_FOR_CANCELLATION;
	}
	if (args->ret) {
		args->message = g_strdup(_("Could not match stars from the catalogue"));
		goto clearup;
	}

	double conv = DBL_MAX;
	solution->px_cat_center = siril_world_cs_ref(args->cat_center);

	if (!check_affine_TRANS_sanity(&trans)) {
		args->message = g_strdup(_("Transformation matrix is invalid, solve failed"));
		args->ret = 1;
		goto clearup;
	}

	double ra0, dec0;
	// star coordinates were set with the origin at the grid center and y upwards
	double center[2] = {0., 0.};

	apply_match(solution->px_cat_center, center, &trans, &ra0, &dec0);
	int num_matched = trans.nm;
	int trial = 0;

	/* try to get a better solution in case of uncentered solution */
	conv = get_center_offset_from_trans(&trans);
	siril_debug_print("iteration %d - offset: %.3f, number of matches: %d\n", trial, conv, num_matched);
	while (conv > CONV_TOLERANCE && trial < max_trials){
		double resolution = get_resolution_from_trans(&trans);
		double focal = RADCONV * solution->pixel_size / resolution;
		siril_debug_print("Current focal: %0.2fmm\n", focal);

		// we will reproject the catalog at the new image center
		siril_world_cs_unref(args->cat_center);
		siril_world_cs_unref(solution->px_cat_center);
		args->cat_center = siril_world_cs_new_from_a_d(ra0, dec0);
		solution->px_cat_center = siril_world_cs_new_from_a_d(ra0, dec0);
		siril_debug_print("Reprojecting to: alpha: %s, delta: %s\n",
				siril_world_cs_alpha_format(args->cat_center, "%02d %02d %.3lf"),
				siril_world_cs_delta_format(args->cat_center, "%c%02d %02d %.3lf"));
		// this simply reprojects the catalog to the new center (no fetch)
		if (get_catalog_stars(args, FALSE)) {
			args->message = g_strdup(_("Reprojecting catalog failed."));
			args->ret = 1;
			break;
		}
		// Uses the indexes in star_list_B to update the stars positions according to the new projection
		update_stars_positions(&star_list_B, num_matched, args->cstars);
		// and recompute the trans structure
		if (atRecalcTrans(num_matched, star_list_A, num_matched, star_list_B, AT_MATCH_MAXITER, AT_MATCH_HALTSIGMA, &trans)) {
			args->message = g_strdup(_("Updating trans failed."));
			args->ret = 1;
			break;
		}
		print_trans(&trans);
		num_matched = trans.nm;
		apply_match(solution->px_cat_center, center, &trans, &ra0, &dec0);
		conv = get_center_offset_from_trans(&trans);
		trial++;
		siril_debug_print("iteration %d - offset: %.3f, number of matches: %d\n", trial, conv, num_matched);
		CHECK_FOR_CANCELLATION;
	}
	if (args->ret)	// after the break
		goto clearup;
	debug_print_catalog_files(star_list_A, star_list_B);

	double resolution = get_resolution_from_trans(&trans);
	solution->focal_length = RADCONV * solution->pixel_size / resolution;
	solution->image_center = siril_world_cs_new_from_a_d(ra0, dec0);
	if (max_trials == 0) {
		siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
	} else if (trial == max_trials) {
		siril_debug_print("No convergence found: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
	} else {
		siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f at iteration #%d\n", ra0, dec0, trial);
	}

	if (args->downsample)
		solution->focal_length *= args->scalefactor;

	solution->image_is_flipped = image_is_flipped_from_trans(&trans);
	double flip = (solution->image_is_flipped) ? -1. : 1.; // when the image is not flipped, x positive in native coords is to the left (East)

	/* compute cd matrix */
	double ra1, dec1, ra2, dec2, cd[2][2];

	/* first, convert center coordinates from deg to rad: */
	dec0 *= DEGTORAD;
	ra0 *= DEGTORAD;

	// We will use 2 points at -flip*dx and dy from the computed center away from the image center
	// to define the left(point1)  and the up (point2) directions
	// (the -flip factor is used to make the step towards increasing ra)
	// We use Appendix B.1 from WCS paper II (http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf):
	// - we compute l', m', n' of each point (the 3D coordinates of the point on the unit sphere)
	// - we then compute its l and m (the coordinates of the point on the unit sphere shifted to the native pole)
	// l is the N coordinate (increasing dec, pointing to NCP), m is the E coordinate (increasing ra)
	// we can then compute the 2 rotations wrt. N

	// Note: the reason for all this is that close to Celestial poles, ra varies largely for small angular distances
	// it's even undetermined (see https://www.atnf.csiro.au/people/mcalabre/WCS/notes_20040211.pdf) when exactly at poles
	// So that if dec0 ≈ ±90deg, ra0 returned can take any value and this will go into the header as CRVAL1,
	// which will in turn determine the first Euler angle of the spherical rotation to the native pole (that is ra0+90)
	// We therefore need to have a CD matrix consistant with this ra0 value
	// hence why we express the left and up vectors at the native pole though their (l,m) coordinates

	/* make a step in direction crpix1 accounting for flip or not*/
	double crpix1[] = {-flip * PLATESOLVE_STEP, 0. };
	apply_match(solution->px_cat_center, crpix1, &trans, &ra1, &dec1);
	siril_debug_print("alpha1: %0.8f, delta1: %0.8f\n", ra1, dec1);
	dec1 *= DEGTORAD;
	ra1 *= DEGTORAD;
	double l1 = sin(dec1) * cos(dec0) - cos(dec1) * sin(dec0) * cos(ra1 - ra0);
	double m1 = cos(dec1) * sin(ra1 - ra0);
	siril_debug_print("l1: %0.8f, m1: %0.8f\n", l1, m1);
	double crota1 = -atan2(l1, m1);

	/* make a step in direction crpix2*/
	double crpix2[] = { 0., PLATESOLVE_STEP};
	apply_match(solution->px_cat_center, crpix2, &trans, &ra2, &dec2);
	siril_debug_print("alpha2: %0.8f, delta2: %0.8f\n", ra2, dec2);
	dec2 *= DEGTORAD;
	ra2 *= DEGTORAD;
	double l2 = sin(dec2) * cos(dec0) - cos(dec2) * sin(dec0) * cos(ra2 - ra0);
	double m2 = cos(dec2) * sin(ra2 - ra0);
	siril_debug_print("l2: %0.8f, m2: %0.8f\n", l2, m2);
	double crota2 = atan2(m2, l2);
	siril_debug_print("crota1: %0.2f, crota2: %0.2f\n", crota1 * RADTODEG, crota2 * RADTODEG);

	double cdelt1 = 0., cdelt2 = 0.;
	get_cdelt_from_trans(&trans, flip, &cdelt1, &cdelt2);

	cd[0][0] =  cdelt1*cos(crota1);
	cd[0][1] = -cdelt1*sin(crota1) * flip;
	cd[1][0] =  cdelt2*sin(crota2) * flip;
	cd[1][1] =  cdelt2*cos(crota2);

	CHECK_FOR_CANCELLATION;

	// saving state for undo before modifying fit structure
	if (!com.script) {
		const char *undo_str = args->for_photometry_cc ? _("Photometric CC") : _("Plate Solve");
		undo_save_state(args->fit, undo_str);
	}

	solution->crpix[0] = args->rx_solver * 0.5;
	solution->crpix[1] = args->ry_solver * 0.5;
	solution->crpix[0] *= args->scalefactor;
	solution->crpix[1] *= args->scalefactor;
	// we now go back to FITS convention (from siril)
	solution->crpix[0] += 0.5;
	solution->crpix[1] += 0.5;

	/**** Fill wcsdata fit structure ***/
	args->fit->wcsdata.ra = siril_world_cs_get_alpha(solution->image_center);
	args->fit->wcsdata.dec = siril_world_cs_get_delta(solution->image_center);
	args->fit->wcsdata.pltsolvd = TRUE;
	g_snprintf(args->fit->wcsdata.pltsolvd_comment, FLEN_COMMENT, "Siril internal solver");
	gchar *ra = siril_world_cs_alpha_format(solution->image_center, "%02d %02d %.3lf");
	gchar *dec = siril_world_cs_delta_format(solution->image_center, "%c%02d %02d %.3lf");
	g_sprintf(args->fit->wcsdata.objctra, "%s", ra);
	g_sprintf(args->fit->wcsdata.objctdec, "%s", dec);
	g_free(ra);
	g_free(dec);


	/**** Fill wcslib fit structure ***/
	wcsprm_t *prm = calloc(1, sizeof(wcsprm_t));
	prm->flag = -1;
	wcsinit(1, 2, prm, 0, 0, 0);
	prm->equinox = 2000.0;
	prm->crpix[0] = solution->crpix[0];
	prm->crpix[1] = solution->crpix[1];
	prm->crval[0] = args->fit->wcsdata.ra;
	prm->crval[1] = args->fit->wcsdata.dec;
	// if (args->trans_order == AT_TRANS_LINEAR) {
		const char CTYPE[2][9] = { "RA---TAN", "DEC--TAN" };
		const char CUNIT[2][4] = { "deg", "deg" };
		for (int i = 0; i < NAXIS; i++) {
			strncpy(prm->cunit[i], &CUNIT[i][0], 71); // 72 char fixed buffer, keep 1 for the NULL
		}
		for (int i = 0; i < NAXIS; i++) {
			strncpy(prm->ctype[i], &CTYPE[i][0], 71); // 72 byte buffer, leave 1 byte for the NULL
		}
	// } else {
	// 	const char CTYPE[2][12] = { "RA---TAN-SIP", "DEC--TAN-SIP" };
	// 	const char CUNIT[2][4] = { "deg", "deg" };
	// 	for (int i = 0; i < NAXIS; i++) {
	// 		strncpy(prm->cunit[i], &CUNIT[i][0], 71); // 72 char fixed buffer, keep 1 for the NULL
	// 	}
	// 	for (int i = 0; i < NAXIS; i++) {
	// 		strncpy(prm->ctype[i], &CTYPE[i][0], 71); // 72 byte buffer, leave 1 byte for the NULL
	// 	}	
	}
	/* PC + CDELT seems to be the preferred approach
	 * according to Calabretta private discussion
	 *
	 *    |cd11 cd12|  = |cdelt1      0| * |pc11 pc12|
	 *    |cd21 cd22|    |0      cdelt2|   |pc21 pc22|
	 */
	// we pass cd and let wcslib decompose it
	wcs_mat2cd(prm, cd);
	prm->altlin = 2;
	wcspcx(prm, 0, 0, NULL);
	if (args->trans_order > 1) {
		disinit(TRUE, 2, prm->lin.dispre, 34);

	}
	wcs_print(prm);
	if (args->fit->wcslib)
		wcsfree(args->fit->wcslib);
	args->fit->wcslib = prm;

	siril_debug_print("****Solution found: WCS data*************\n");
	siril_debug_print("crpix1 = %*.12e\n", 20, solution->crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, solution->crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, args->fit->wcsdata.ra);
	siril_debug_print("crval2 = %*.12e\n", 20, args->fit->wcsdata.dec);
	siril_debug_print("cdelt1 = %*.12e\n", 20, args->fit->wcslib->cdelt[0]);
	siril_debug_print("cdelt2 = %*.12e\n", 20, args->fit->wcslib->cdelt[1]);
	siril_debug_print("pc1_1  = %*.12e\n", 20, args->fit->wcslib->pc[0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, args->fit->wcslib->pc[1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, args->fit->wcslib->pc[2]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, args->fit->wcslib->pc[3]);
	siril_debug_print("******************************************\n");

	CHECK_FOR_CANCELLATION;

	if (args->verbose)
		print_platesolving_results_from_wcs(args);
clearup:
	free_stars(&star_list_A);
	free_stars(&star_list_B);
	return args->ret;
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
		siril_debug_print(error->message);
		asnet_version = g_strdup("ansvr");
	} else {
		gchar** chunks = g_strsplit(child_stdout, "\n", 2);
		if (strlen(chunks[1]) == 0) { // if the second chunk is void, it means there was only one line which contains the version
			asnet_version = g_strdup(chunks[0]);
		} else {
			asnet_version = g_strdup("<0.88"); // version switch was introduced from 0.88 (https://github.com/dstndstn/astrometry.net/commit/25f0b829d80c57984d404de50f9c645cb3da2223)
		}
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
	g_free(bin);
	if (!asnet_version)
		get_asnet_version(bin);
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

	// saving state for undo before modifying fit structure
	if (!com.script) {
		undo_save_state(args->fit, _("Plate Solve"));
	}

	memcpy(&args->fit->wcsdata, &result.wcsdata, sizeof(wcs_info));
	memset(&result.wcsdata, 0, sizeof(wcs_info));
	args->fit->wcslib = result.wcslib;
	result.wcslib = NULL;
	clearfits(&result);
	if (!com.pref.astrometry.keep_wcs_files)
		g_unlink(wcs_filename);
	g_free(wcs_filename);

	solution->image_is_flipped = image_is_flipped_from_wcs(args->fit);

	// For some reason, asnet may not return a solution with the ref point at the center
	// We need to account for that
	double ra0, dec0;
	center2wcs(args->fit, &ra0, &dec0);
	args->fit->wcsdata.ra  = ra0;
	args->fit->wcsdata.dec = dec0;

	double resolution = get_wcs_image_resolution(args->fit) * 3600.0;
	solution->focal_length = RADCONV * args->pixel_size / resolution;

	args->fit->wcsdata.pltsolvd = TRUE;
	if (args->fit->wcsdata.pltsolvd_comment[0] != '\0')
		memset(args->fit->wcsdata.pltsolvd_comment, 0, FLEN_COMMENT);
	snprintf(args->fit->wcsdata.pltsolvd_comment, FLEN_COMMENT, "Solved by Astrometry.net (%s)", asnet_version);
	if (args->verbose)
		siril_log_color_message(_("Local astrometry.net solve succeeded.\n"), "green");

	// asnet puts more info in the HISTORY and the console log in COMMENT fields
	solution->image_center = siril_world_cs_new_from_a_d(
			ra0,
			dec0);
	/* print results from WCS data */
	print_updated_wcs_data(args->fit);

	if (args->verbose)
		print_platesolving_results_from_wcs(args);
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
	}
}

static int astrometry_prepare_hook(struct generic_seq_args *arg) {
	struct astrometry_data *args = (struct astrometry_data *)arg->user;
	fits fit = { 0 };
	// load ref metadata in fit
	if (seq_read_frame_metadata(arg->seq, sequence_find_refimage(arg->seq), &fit))
		return 1;
	if (!args->cat_center) {
		args->cat_center = get_eqs_from_header(&fit);
		if (args->cat_center) {
			args->ref_stars->center_ra = siril_world_cs_get_alpha(args->cat_center);
			args->ref_stars->center_dec = siril_world_cs_get_delta(args->cat_center);
		}
	}
	if (args->ref_stars->cat_index != CAT_ASNET && !args->cat_center) {
		siril_log_color_message(_("Cannot plate solve, no target coordinates passed and image header doesn't contain any either\n"), "red");
		return 1;
	}
	if (args->pixel_size <= 0.0) {
		args->pixel_size = max(fit.pixel_size_x, fit.pixel_size_y);
		if (args->pixel_size <= 0.0) {
			args->pixel_size = com.pref.starfinder_conf.pixel_size_x;
			if (args->pixel_size <= 0.0) {
				siril_log_color_message(_("Pixel size not found in image or in settings, cannot proceed\n"), "red");
				return 1;
			}
		}
	}
	if (args->focal_length <= 0.0) {
		args->focal_length = fit.focal_length;
		if (args->focal_length <= 0.0) {
			args->focal_length = com.pref.starfinder_conf.focal_length;
			if (args->focal_length <= 0.0) {
				siril_log_color_message(_("Focal length not found in image or in settings, cannot proceed\n"), "red");
				return 1;
			}
		}
	}

	seq_prepare_hook(arg);
	args->fit = &fit;
	process_plate_solver_input(args); // compute required data to get the catalog
	clearfits(&fit);
	if (args->ref_stars->cat_index == CAT_ASNET) {
		com.child_is_running = EXT_ASNET;
		g_unlink("stop"); // make sure the flag file for cancel is not already in the folder
	}
	return get_catalog_stars(args, TRUE);
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
	process_plate_solver_input(aargs);	// depends on aargs->fit
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
	if (aargs->cstars)
		free_fitted_stars(aargs->cstars);
	free (aargs);
	com.child_is_running = EXT_NONE;
	if (g_unlink("stop"))
		siril_debug_print("g_unlink() failed\n");
	return 0;
}

void start_sequence_astrometry(sequence *seq, struct astrometry_data *args) {
	struct generic_seq_args *seqargs = create_default_seqargs(seq);
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = seq->selnum;
	seqargs->stop_on_error = FALSE;
#ifdef _WIN32
	seqargs->parallel = args->ref_stars->cat_index != CAT_ASNET;		// for now crashes on Cancel if parallel is enabled for asnet on windows
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

