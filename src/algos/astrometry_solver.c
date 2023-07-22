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
#include "algos/annotate.h"
#include "algos/photometry.h"
#include "algos/photometric_cc.h"
#include "algos/siril_wcs.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/local_catalogues.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "registration/matching/match.h"
#include "registration/matching/apply_match.h"
#include "registration/matching/atpmatch.h"
#include "registration/matching/project_coords.h"
#include "gui/message_dialog.h"


#define DOWNSAMPLE_FACTOR 0.25
#define CONV_TOLERANCE 1E-8

#undef DEBUG		/* get some of diagnostic output */

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
	return round(100.0 * min(20.0, max(7.0, m))) / 100;
}

static void compute_limit_mag(struct astrometry_data *args) {
	if (args->mag_mode == LIMIT_MAG_ABSOLUTE)
		args->limit_mag = args->magnitude_arg;
	else {
		args->limit_mag = compute_mag_limit_from_fov(args->used_fov / 60.0);
		if (args->mag_mode == LIMIT_MAG_AUTO_WITH_OFFSET)
			args->limit_mag += args->magnitude_arg;
	}
	siril_debug_print("using limit magnitude %f\n", args->limit_mag);
}

static gchar *project_catalog(GFile *catalogue_name, SirilWorldCS *catalog_center) {
	GError *error = NULL;
	gchar *foutput = NULL;
	/* --------- Project coords of Vizier catalog and save it into catalog.proj ------- */

	GFile *fproj = g_file_new_build_filename(g_get_tmp_dir(), "catalog.proj", NULL);

	/* We want to remove the file if already exisit */
	if (!g_file_delete(fproj, NULL, &error)
			&& !g_error_matches(error, G_IO_ERROR, G_IO_ERROR_NOT_FOUND)) {
		// deletion failed for some reason other than the file not existing:
		// so report the error
		g_warning("Failed to delete %s: %s", g_file_peek_path(fproj),
				error->message);
	}

	convert_catalog_coords(catalogue_name, catalog_center, fproj);
	foutput = g_file_get_path(fproj);
	g_object_unref(fproj);
	return foutput;
}

gboolean has_any_keywords() {
	return (gfit.focal_length > 0.0 ||
			gfit.pixel_size_x > 0.f ||
			gfit.pixel_size_y > 0.f ||
			(gfit.wcsdata.crval[0] > 0.0 && gfit.wcsdata.crval[1] != 0.0) ||
			(gfit.wcsdata.objctra[0] != '\0' && gfit.wcsdata.objctdec[0] != '\0') ||
			(gfit.wcsdata.ra != 0.0 && gfit.wcsdata.dec != 0.0));
}

SirilWorldCS *get_eqs_from_header(fits *fit) {
	if (fit->wcsdata.ra != 0.0 || fit->wcsdata.dec != 0.0)
		return siril_world_cs_new_from_a_d(fit->wcsdata.ra, fit->wcsdata.dec);

	else if (fit->wcsdata.objctra[0] != '\0' && fit->wcsdata.objctdec[0] != '\0')
		return siril_world_cs_new_from_objct_ra_dec(fit->wcsdata.objctra, fit->wcsdata.objctdec);

	else if (fit->wcsdata.crval[0] != 0.0 || fit->wcsdata.crval[1] != 0.0)
		return siril_world_cs_new_from_a_d(fit->wcsdata.crval[0], fit->wcsdata.crval[1]);
	return NULL;
}

/* Extract CDELT from CD matrix.*/
static void extract_cdelt_from_cd(double cd1_1, double cd1_2, double cd2_1,
		double cd2_2, double *cdelt1, double *cdelt2) {
	int sign;
	if ((cd1_1 * cd2_2 - cd1_2 * cd2_1) >= 0)
		sign = +1;
	else
		sign = -1;

	*cdelt1 = sqrt((cd1_1 * cd1_1) + (cd2_1 * cd2_1)) * sign;
	*cdelt2 = sqrt((cd1_2 * cd1_2) + (cd2_2 * cd2_2));
}

static void print_platesolving_results_from_wcs(struct astrometry_data *args) {
	double rotationa, rotationb, rotation;
	char field_x[256] = "";
	char field_y[256] = "";

	double cd[2][2];
	wcs_pc_to_cd(args->fit->wcsdata.pc, args->fit->wcsdata.cdelt, cd);
	rotationa = atan2(-args->fit->wcsdata.pc[1][0], args->fit->wcsdata.pc[0][0]);
	rotationb = atan2(args->fit->wcsdata.pc[0][1], args->fit->wcsdata.pc[1][1]);
	rotation = 0.5 * (rotationa + rotationb) * RADTODEG;

	double det = (cd[0][0] * cd[1][1] - cd[1][0] * cd[0][1]); // determinant of rotation matrix (ad - bc)
	/* If the determinant of the top-left 2x2 rotation matrix is < 0
	 * the transformation is orientation-preserving. */

	if (det > 0 && args->flip_image)
		rotation = 180.0 - rotation;
	if (rotation < -180.0)
		rotation += 360.0;
	if (rotation > 180.0)
		rotation -= 360.0;
	siril_log_message(_("Up is %+.2lf deg ClockWise wrt. N%s\n"), rotation, det > 0.0 ? _(" (flipped)") : "");

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

static TRANS H_to_linear_TRANS(Homography H) {
	TRANS trans = { 0 };

	trans.order = AT_TRANS_LINEAR;

	trans.a = H.h02;
	trans.b = H.h00;
	trans.c = H.h01;
	trans.d = H.h12;
	trans.e = H.h10;
	trans.f = H.h11;

	return trans;
}

static gboolean check_affine_TRANS_sanity(TRANS *trans) {
	double var1 = fabs(trans->b) - fabs(trans->f);
	double var2 = fabs(trans->c) - fabs(trans->e);
	siril_debug_print("abs(b+f)=%f et abs(c+e)=%f\n", var1, var2);

	return ((fabs(var1) < 0.3) && fabs(var2) < 0.3);
}

static gboolean image_is_flipped(Homography H) {
	double det = (H.h00 * H.h11 - H.h01 * H.h10); // determinant of rotation matrix (ad - bc)
	return det < 0;
}

static gboolean image_is_flipped_from_wcs(fits *fit) {
	double cd[2][2];
	wcs_pc_to_cd(fit->wcsdata.pc, fit->wcsdata.cdelt, cd);
	double det = (cd[0][0] * cd[1][1] - cd[1][0] * cd[0][1]); // determinant of rotation matrix (ad - bc)
	return det > 0; // convention is that angles are positive clockwise when image is not flipped
}

// From projected starlist and center (ra,dec), go back to original ra and dec
// All formulas from AIPS memo #27 III.A.ii
// https://library.nrao.edu/public/memos/aips/memos/AIPSM_027.pdf

static void deproject_starlist(int num_stars, s_star *star_list, double ra0, double dec0, int doASEC) {
	ra0 *= DEGTORAD;
	dec0 *= DEGTORAD;
	s_star *currstar;
	currstar = star_list;
	for (int i = 0; i < num_stars; i++) {
		double xi = currstar->x;
		double eta = currstar->y;
		if (doASEC > 0) {
			xi /= RADtoASEC;
			eta /= RADtoASEC;
		}
		double delta_ra = atan(xi / (cos(dec0) - eta * sin(dec0)));
		double ra = ra0 + delta_ra;
		double dec = atan(cos(delta_ra) * (eta * cos(dec0) + sin(dec0)) / (cos(dec0) - eta * sin(dec0)));
		currstar->x = ra / DEGTORAD;
		currstar->y = dec / DEGTORAD;
		currstar = currstar->next;
	}
}

// From starlist in (ra,dec) and center (ra,dec), project in "pixels" (in arcsec)
// All formulas from AIPS memo #27 III.A.i
// https://library.nrao.edu/public/memos/aips/memos/AIPSM_027.pdf

static void project_starlist(int num_stars, s_star *star_list, double ra0, double dec0, int doASEC) {
	double delta_ra;
	dec0 *= DEGTORAD;
	s_star *currstar;
	currstar = star_list;
	for (int i = 0; i < num_stars; i++) {
		double ra = currstar->x;
		double dec = currstar->y;
		if ((ra < 10) && (ra0 > 350)) {
			delta_ra = (ra + 360) - ra0;
		} else if ((ra > 350) && (ra0 < 10)) {
			delta_ra = (ra - 360) - ra0;
		} else {
			delta_ra = ra - ra0;
		}
		delta_ra *= DEGTORAD;
		dec *= DEGTORAD;

		/*
		 * let's transform from (delta_RA, delta_Dec) to (xi, eta),
		 */
		double xx = cos(dec) * sin(delta_ra);
		double yy = sin(dec0) * sin(dec) + cos(dec0) * cos(dec) * cos(delta_ra);
		double xi = (xx / yy);
		xx = cos(dec0) * sin(dec) - sin(dec0) * cos(dec) * cos(delta_ra);
		double eta = (xx / yy);

		if (doASEC > 0) {
			xi *= RADtoASEC;
			eta *= RADtoASEC;
		}
		currstar->x = xi;
		currstar->y = eta;
		currstar = currstar->next;
	}
}

void print_updated_wcs_data(fits *fit) {
	/* debug output */
	siril_debug_print("****Updated WCS data*************\n");
	siril_debug_print("crpix1 = %*.12e\n", 20, fit->wcsdata.crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, fit->wcsdata.crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, fit->wcsdata.crval[0]);
	siril_debug_print("crval2 = %*.12e\n", 20, fit->wcsdata.crval[1]);
	siril_debug_print("cdelt1 = %*.12e\n", 20, fit->wcsdata.cdelt[0]);
	siril_debug_print("cdelt2 = %*.12e\n", 20, fit->wcsdata.cdelt[1]);
	siril_debug_print("pc1_1  = %*.12e\n", 20, fit->wcsdata.pc[0][0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, fit->wcsdata.pc[0][1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, fit->wcsdata.pc[1][0]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, fit->wcsdata.pc[1][1]);
	siril_debug_print("******************************************\n");
}

/******
 *
 * Public functions
 */

void flip_bottom_up_astrometry_data(fits *fit) {
	/* flip pc matrix */
	fit->wcsdata.pc[0][1] = -fit->wcsdata.pc[0][1];
	fit->wcsdata.pc[1][1] = -fit->wcsdata.pc[1][1];
	double cd[2][2];
	wcs_pc_to_cd(fit->wcsdata.pc, fit->wcsdata.cdelt, cd);
	wcs_cd_to_pc(cd, fit->wcsdata.pc, fit->wcsdata.cdelt);

	/* update crpix */
	fit->wcsdata.crpix[1] = fit->ry - fit->wcsdata.crpix[1];

	print_updated_wcs_data(fit);
}

void reframe_astrometry_data(fits *fit, Homography H) {
	double pc1_1, pc1_2, pc2_1, pc2_2;
	point refpointout;

	pc1_1 = H.h00 * fit->wcsdata.pc[0][0] + H.h01 * fit->wcsdata.pc[0][1];
	pc1_2 = H.h10 * fit->wcsdata.pc[0][0] + H.h11 * fit->wcsdata.pc[0][1];
	pc2_1 = H.h00 * fit->wcsdata.pc[1][0] + H.h01 * fit->wcsdata.pc[1][1];
	pc2_2 = H.h10 * fit->wcsdata.pc[1][0] + H.h11 * fit->wcsdata.pc[1][1];
	// we go back to cd formulation just to separate back again cdelt and pc
	double cd[2][2];
	fit->wcsdata.pc[0][0] = pc1_1;
	fit->wcsdata.pc[0][1] = pc1_2;
	fit->wcsdata.pc[1][0] = pc2_1;
	fit->wcsdata.pc[1][1] = pc2_2;
	wcs_pc_to_cd(fit->wcsdata.pc, fit->wcsdata.cdelt, cd);
	wcs_cd_to_pc(cd, fit->wcsdata.pc, fit->wcsdata.cdelt);

	point refpointin = {fit->wcsdata.crpix[0], fit->wcsdata.crpix[1]};
	cvTransformImageRefPoint(H, refpointin, &refpointout);

	fit->wcsdata.crpix[0] = refpointout.x;
	fit->wcsdata.crpix[1] = refpointout.y;

	print_updated_wcs_data(fit);
}

void wcs_cd_to_pc(double cd[][2], double pc[][2], double cdelt[2]) {
	extract_cdelt_from_cd(cd[0][0], cd[0][1], cd[1][0], cd[1][1], &cdelt[0], &cdelt[1]);

	pc[0][0] = cd[0][0] / cdelt[0];
	pc[0][1] = cd[0][1] / cdelt[0];
	pc[1][0] = cd[1][0] / cdelt[1];
	pc[1][1] = cd[1][1] / cdelt[1];
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
static int get_catalog_stars(struct astrometry_data *args) {
	if (args->onlineCatalog == CAT_ASNET)
		return 0;

	/* obtaining a star catalogue */
	if (!args->catalog_file && !args->use_local_cat) {
		args->catalog_file = download_catalog(args->onlineCatalog, args->cat_center,
				args->used_fov * 0.5, args->limit_mag);
		if (!args->catalog_file) {
			args->message = g_strdup(_("Could not download the online star catalogue."));
			return 1;
		}
	}
	CHECK_FOR_CANCELLATION_RET;

	args->cstars = new_fitted_stars(MAX_STARS);
	if (!args->cstars) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	/* project and open the file */
	gchar *catalogStars;	// file name of the projected catalog
	if (args->use_local_cat) {
		catalogStars = get_and_project_local_catalog(args->cat_center,
				args->used_fov / 120.0, args->limit_mag, FALSE);
	} else {
		catalogStars = project_catalog(args->catalog_file, args->cat_center);
		if (!catalogStars) {
			args->message = g_strdup(_("Cannot project the star catalog."));
			return 1;
		}
	}
	CHECK_FOR_CANCELLATION_RET;
	GFile *catalog = g_file_new_for_path(catalogStars);
	GError *error = NULL;
	GInputStream *input_stream = (GInputStream*) g_file_read(catalog, NULL, &error);
	if (!input_stream) {
		if (error != NULL) {
			args->message = g_strdup_printf(_("Could not load the star catalog (%s)."), error->message);
			g_clear_error(&error);
		}
		args->message = g_strdup_printf(_("Could not load the star catalog (%s)."), "generic error");
		return 1;
	}

	args->n_cat = read_projected_catalog(input_stream, args->cstars, args->onlineCatalog);
	g_object_unref(input_stream);
	g_object_unref(catalog);
	g_free(catalogStars);
	return 0;
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
		if (args->onlineCatalog == CAT_ASNET) {
			siril_log_message(_("Plate solving image from local catalogues for a field of view of %.2f degrees\n"), args->used_fov / 60.0);
		} else if (args->use_local_cat) {
			siril_log_message(_("Plate solving image from local catalogues for a field of view of %.2f"
						" degrees%s, using a limit magnitude of %.2f\n"),
					args->used_fov / 60.0,
					args->uncentered ? _(" (uncentered)") : "", args->limit_mag);
		} else {
			siril_log_message(_("Plate solving image from an online catalogue for a field of view of %.2f"
						" degrees%s, using a limit magnitude of %.2f\n"),
					args->used_fov / 60.0,
					args->uncentered ? _(" (uncentered)") : "", args->limit_mag);
		}
	}

	/* 1. Get catalogue stars for the field of view (for sequences, see the prepare hook) */
	if (!args->for_sequence && get_catalog_stars(args)) {
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
		int nthreads = (!args->for_sequence || args->onlineCatalog == CAT_ASNET) ? com.max_thread : 1;
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

	if (args->onlineCatalog == CAT_ASNET) {
		if (!args->for_sequence) {
			com.child_is_running = EXT_ASNET;
			g_unlink("stop"); // make sure the flag file for cancel is not already in the folder
		}
		if (local_asnet_platesolve(stars, nb_stars, args, &solution)) {
			args->ret = ERROR_PLATESOLVE;
		}
	} else
		if (match_catalog(stars, nb_stars, args, &solution)) {
			args->ret = ERROR_PLATESOLVE;
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
#ifndef HAVE_WCSLIB
		siril_log_color_message(_("This operation (PCC) relies on the missing WCSLIB software, cannot continue.\n"), "red");
		args->ret = ERROR_PLATESOLVE;
		goto clearup;
#endif
		if (args->use_local_cat) {
			double tra = siril_world_cs_get_alpha(solution.image_center);
			double tdec = siril_world_cs_get_delta(solution.image_center);
			double res = get_resolution(solution.focal_length, args->pixel_size);
			double radius = get_radius_deg(res, args->fit->rx, args->fit->ry);
			// for photometry, we can use fainter stars, 1.5 seems ok above instead of 2.0
			if (args->verbose)
				siril_log_message(_("Getting stars from local catalogues for PCC, limit magnitude %.2f\n"), args->limit_mag);
			if (get_photo_stars_from_local_catalogues(tra, tdec, radius, args->fit, args->limit_mag, &pcc_stars, &nb_pcc_stars)) {
				siril_log_color_message(_("Failed to get data from the local catalogue, is it installed?\n"), "red");
				args->ret = ERROR_PHOTOMETRY;
			}
		} else {
			args->ret = project_catalog_with_WCS(args->catalog_file, args->fit, TRUE,
					&pcc_stars, &nb_pcc_stars);
		}
		if (args->ret) {
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
		load_WCS_from_memory(args->fit);
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
		if (args->catalog_file)
			g_object_unref(args->catalog_file);
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
			siril_debug_print("g_unlink() failed");
		siril_add_idle(end_plate_solver, args);
	}
	else free(args);
	return GINT_TO_POINTER(retval);
}

/* entry point for siril's plate solver based on catalogue matching */
static int match_catalog(psf_star **stars, int nb_stars, struct astrometry_data *args, solve_results *solution) {
	Homography H = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;
	int max_trials = 0;
	s_star *star_list_A = NULL, *star_list_B = NULL;

	if (args->uncentered)
		max_trials = 20; //retry to converge if solve is done at an offset from the center

	/* make sure that arrays are not too small
	 * make sure that the max of stars is BRIGHTEST_STARS */
	int n = min(min(nb_stars, args->n_cat), BRIGHTEST_STARS);

	double scale_min = args->scale - 0.2;
	double scale_max = args->scale + 0.2;
	int attempt = 1;
	while (args->ret && attempt <= 3) {
		free_stars(&star_list_A);
		free_stars(&star_list_B);
		args->ret = new_star_match(stars, args->cstars, n, nobj,
				scale_min, scale_max, &H, TRUE,
				AFFINE_TRANSFORMATION, &star_list_A, &star_list_B);
		if (attempt == 1) {
			scale_min = -1.0;
			scale_max = -1.0;
			nobj += 10;
		} else {
			nobj += 30;
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
	/* we only want to compare with linear function
	 * Maybe one day we will apply match with homography matrix
	 */
	TRANS trans = H_to_linear_TRANS(H);
	if (!check_affine_TRANS_sanity(&trans)) {
		args->message = g_strdup(_("Transformation matrix is invalid, solve failed"));
		args->ret = 1;
		goto clearup;
	}

	double ra0, dec0;
	// using siril convention as were set by the peaker
	solution->crpix[0] = args->rx_solver * 0.5;
	solution->crpix[1] = args->ry_solver * 0.5;

	apply_match(solution->px_cat_center, solution->crpix, &trans, &ra0, &dec0);
	int num_matched = H.pair_matched;
	int trial = 0;

	/* try to get a better solution in case of uncentered selection */
	while (conv > CONV_TOLERANCE && trial < max_trials){
		double rainit = siril_world_cs_get_alpha(args->cat_center);
		double decinit = siril_world_cs_get_delta(args->cat_center);
		double orig_ra0 = ra0;
		double orig_dec0 = dec0;

		deproject_starlist(num_matched, star_list_B, rainit, decinit, 1);
		siril_debug_print("Deprojecting from: alpha: %s, delta: %s\n",
				siril_world_cs_alpha_format(args->cat_center, "%02d %02d %.3lf"),
				siril_world_cs_delta_format(args->cat_center, "%c%02d %02d %.3lf"));
		args->cat_center = siril_world_cs_new_from_a_d(ra0, dec0);
		siril_world_cs_unref(solution->px_cat_center);
		solution->px_cat_center = siril_world_cs_new_from_a_d(ra0, dec0);

		project_starlist(num_matched, star_list_B, ra0, dec0, 1);
		siril_debug_print("Reprojecting to: alpha: %s, delta: %s\n",
				siril_world_cs_alpha_format(args->cat_center, "%02d %02d %.3lf"),
				siril_world_cs_delta_format(args->cat_center, "%c%02d %02d %.3lf"));

		double scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
		double scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
		double resolution = (scaleX + scaleY) * 0.5; // we assume square pixels

		double focal = RADCONV * solution->pixel_size / resolution;
		siril_debug_print("Current focal: %0.2fmm\n", focal);

		if (atPrepareHomography(num_matched, star_list_A, num_matched, star_list_B, &H, TRUE, AFFINE_TRANSFORMATION)){
			args->message = g_strdup(_("Updating homography failed."));
			args->ret = 1;
			break;
		}
		trans = H_to_linear_TRANS(H);
		apply_match(solution->px_cat_center, solution->crpix, &trans, &ra0, &dec0);

		conv = fabs((dec0 - orig_dec0) / args->used_fov / 60.) + fabs((ra0 - orig_ra0) / args->used_fov / 60.);

		trial++;
		CHECK_FOR_CANCELLATION;
	}
	if (args->ret)	// after the break
		goto clearup;

	memcpy(&solution->H, &H, sizeof(Homography));
	double scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	double scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	double resolution = (scaleX + scaleY) * 0.5; // we assume square pixels
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

	solution->image_is_flipped = image_is_flipped(H);

	/* compute cd matrix */
	double ra7, dec7, delta_ra;

	/* first, convert center coordinates from deg to rad: */
	dec0 *= DEGTORAD;
	ra0 *= DEGTORAD;

	/* make 1 step in direction crpix1 */
	double crpix1[] = { solution->crpix[0] + 1.0 / args->scalefactor, solution->crpix[1] };
	apply_match(solution->px_cat_center, crpix1, &trans, &ra7, &dec7);

	dec7 *= DEGTORAD;
	ra7 *= DEGTORAD;

	delta_ra = ra7 - ra0;
	if (delta_ra > +M_PI)
		delta_ra = 2.0 * M_PI - delta_ra;
	if (delta_ra < -M_PI)
		delta_ra = delta_ra - 2.0 * M_PI;
	double cd1_1 = (delta_ra) * cos(dec0) * RADTODEG;
	double cd2_1 = (dec7 - dec0) * RADTODEG;

	/* make 1 step in direction crpix2
	 * WARNING: we use -1 because of the Y axis reversing */
	double crpix2[] = { solution->crpix[0], solution->crpix[1] - 1.0 / args->scalefactor };
	apply_match(solution->px_cat_center, crpix2, &trans, &ra7, &dec7);

	dec7 *= DEGTORAD;
	ra7 *= DEGTORAD;

	delta_ra = ra7 - ra0;
	if (delta_ra > +M_PI)
		delta_ra = 2.0 * M_PI - delta_ra;
	if (delta_ra < -M_PI)
		delta_ra = delta_ra - 2.0 * M_PI;
	double cd1_2 = (delta_ra) * cos(dec0) * RADTODEG;
	double cd2_2 = (dec7 - dec0) * RADTODEG;

	CHECK_FOR_CANCELLATION;

	// saving state for undo before modifying fit structure
	if (!com.script) {
		const char *undo_str = args->for_photometry_cc ? _("Photometric CC") : _("Plate Solve");
		undo_save_state(args->fit, undo_str);
	}

	/**** Fill wcsdata fit structure ***/
	args->fit->wcsdata.equinox = 2000.0;

	solution->crpix[0] = args->rx_solver * 0.5;
	solution->crpix[1] = args->ry_solver * 0.5;

	solution->crpix[0] *= args->scalefactor;
	solution->crpix[1] *= args->scalefactor;

	args->fit->wcsdata.ra = siril_world_cs_get_alpha(solution->image_center);
	args->fit->wcsdata.dec = siril_world_cs_get_delta(solution->image_center);

	args->fit->wcsdata.crpix[0] = solution->crpix[0];
	args->fit->wcsdata.crpix[1] = solution->crpix[1];
	args->fit->wcsdata.crval[0] = args->fit->wcsdata.ra;
	args->fit->wcsdata.crval[1] = args->fit->wcsdata.dec;

	args->fit->wcsdata.pltsolvd = TRUE;
	g_snprintf(args->fit->wcsdata.pltsolvd_comment, FLEN_COMMENT, "Siril internal solver");

	gchar *ra = siril_world_cs_alpha_format(solution->image_center, "%02d %02d %.3lf");
	gchar *dec = siril_world_cs_delta_format(solution->image_center, "%c%02d %02d %.3lf");

	g_sprintf(args->fit->wcsdata.objctra, "%s", ra);
	g_sprintf(args->fit->wcsdata.objctdec, "%s", dec);

	g_free(ra);
	g_free(dec);

	CHECK_FOR_CANCELLATION;
	double cdelt1, cdelt2;

	extract_cdelt_from_cd(cd1_1, cd1_2, cd2_1, cd2_2, &cdelt1, &cdelt2);

	args->fit->wcsdata.cdelt[0] = cdelt1;
	args->fit->wcsdata.cdelt[1] = cdelt2;

	/* PC + CDELT seems to be the preferred approach
	 * according to Calabretta private discussion
	 *
	 *    |cd11 cd12|  = |cdelt1      0| * |pc11 pc12|
	 *    |cd21 cd22|    |0      cdelt2|   |pc21 pc22|
	 */

	args->fit->wcsdata.pc[0][0] = cd1_1 / cdelt1;
	args->fit->wcsdata.pc[0][1] = cd1_2 / cdelt1;
	args->fit->wcsdata.pc[1][0] = cd2_1 / cdelt2;
	args->fit->wcsdata.pc[1][1] = cd2_2 / cdelt2;

	siril_debug_print("****Solution found: WCS data*************\n");
	siril_debug_print("crpix1 = %*.12e\n", 20, solution->crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, solution->crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, args->fit->wcsdata.ra);
	siril_debug_print("crval2 = %*.12e\n", 20, args->fit->wcsdata.dec);
	siril_debug_print("cdelt1 = %*.12e\n", 20, args->fit->wcsdata.cdelt[0]);
	siril_debug_print("cdelt2 = %*.12e\n", 20, args->fit->wcsdata.cdelt[1]);
	siril_debug_print("pc1_1  = %*.12e\n", 20, args->fit->wcsdata.pc[0][0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, args->fit->wcsdata.pc[0][1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, args->fit->wcsdata.pc[1][0]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, args->fit->wcsdata.pc[1][1]);
	siril_debug_print("******************************************\n");

	load_WCS_from_memory(args->fit);

	if (args->verbose)
		// print_platesolving_results(solution, args->downsample);
		print_platesolving_results_from_wcs(args);
clearup:
	free_stars(&star_list_A);
	free_stars(&star_list_B);
	return args->ret;
}

/*********************** finding asnet bash first **********************/
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
			return g_build_filename(com.pref.asnet_dir, NULL);
		}
	}
	// searching default location %localappdata%/cygwin_ansvr
	const gchar *localappdata = g_get_user_data_dir();
	gchar *testdir = g_build_filename(localappdata, "cygwin_ansvr", "bin", NULL);
	if (g_file_test(testdir, G_FILE_TEST_IS_DIR)) {
		siril_debug_print("cygwin/bin found at %s\n", testdir);
		g_free(testdir);
		return g_build_filename(localappdata, "cygwin_ansvr", NULL);
	}
	siril_log_color_message(_("cygwin/bin was not found at %s - ignoring\n"), "red", testdir);
	g_free(testdir);
	return NULL;
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
		return TRUE;
	}
	siril_debug_print("solve-field not found in PATH\n");
	gchar *bin = siril_get_asnet_bin();
	if (!bin) return FALSE;
	gboolean is_available = g_file_test(bin, G_FILE_TEST_EXISTS);
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
	if (!asnet_is_available()) {
		siril_log_color_message(_("solve-field was not found, set its path in the preferences\n"), "red");
		return 1;
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
	sprintf(low_scale, "%.3f", args->scale / a);
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
		"-u", "arcsecperpix", "-L", low_scale, "-H", high_scale, NULL };

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

	memcpy(&args->fit->wcsdata, &result.wcsdata, sizeof(wcs_info));
	memset(&result.wcsdata, 0, sizeof(wcs_info));
#ifdef HAVE_WCSLIB
	args->fit->wcslib = result.wcslib;
	result.wcslib = NULL;
#endif
	clearfits(&result);
	if (!com.pref.astrometry.keep_wcs_files)
		g_unlink(wcs_filename);
	g_free(wcs_filename);

	solution->image_is_flipped = image_is_flipped_from_wcs(args->fit);

	// we go back to siril convention
	args->fit->wcsdata.crpix[0] = (double)args->rx_solver * 0.5;
	args->fit->wcsdata.crpix[1] = (double)args->ry_solver * 0.5;
	args->fit->wcsdata.ra = args->fit->wcsdata.crval[0];
	args->fit->wcsdata.dec = args->fit->wcsdata.crval[1];

	double resolution = get_wcs_image_resolution(args->fit) * 3600.0;
	solution->focal_length = RADCONV * args->pixel_size / resolution;

	if (args->downsample) {
		solution->focal_length *= args->scalefactor;

		Homography S;
		cvGetMatrixResize(args->fit->wcsdata.crpix[0], args->fit->wcsdata.crpix[1],
				(double)args->fit->rx * 0.5, (double)args->fit->ry * 0.5, args->scalefactor, &S);
		reframe_astrometry_data(args->fit, S);
	}
	// we need to reload here to make sure everything in fit->wcslib is updated
	// TODO: this is where we loose the SIP info, will need to be smarter than this
	load_WCS_from_memory(args->fit);

	args->fit->wcsdata.pltsolvd = TRUE;
	strcpy(args->fit->wcsdata.pltsolvd_comment, "This WCS header was created by Astrometry.net.");
	if (args->verbose)
		siril_log_color_message(_("Local astrometry.net solve succeeded.\n"), "green");

	// asnet puts more info in the HISTORY and the console log in COMMENT fields
	solution->image_center = siril_world_cs_new_from_a_d(
			args->fit->wcsdata.crval[0],
			args->fit->wcsdata.crval[1]);
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

	if (croparea.w == args->fit->rx && croparea.h == args->fit->ry)
		memset(&croparea, 0, sizeof(rectangle));
	else siril_debug_print("reduced area for the solve: %d, %d, %d x %d%s\n",
			croparea.x, croparea.y, croparea.w, croparea.h,
			args->downsample ? " (down-sampled)" : "");
	memcpy(&(args->solvearea), &croparea, sizeof(rectangle));

	compute_limit_mag(args); // to call after having set args->used_fov
}

static int astrometry_prepare_hook(struct generic_seq_args *arg) {
	struct astrometry_data *args = (struct astrometry_data *)arg->user;
	fits fit = { 0 };
	// load ref metadata in fit
	if (seq_read_frame_metadata(arg->seq, sequence_find_refimage(arg->seq), &fit))
		return 1;
	if (!args->cat_center)
		args->cat_center = get_eqs_from_header(&fit);
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

	args->fit = &fit;
	remove_prefixed_sequence_files(arg->seq, arg->new_seq_prefix);
	remove_prefixed_star_files(arg->seq, arg->new_seq_prefix);
	process_plate_solver_input(args); // compute required data to get the catalog
	clearfits(&fit);
	if (args->onlineCatalog == CAT_ASNET) {
		com.child_is_running = EXT_ASNET;
		g_unlink("stop"); // make sure the flag file for cancel is not already in the folder
	}
	return get_catalog_stars(args);
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
	aargs->filename = g_strdup(root);
	process_plate_solver_input(aargs); // depends on args->fit
	int retval = GPOINTER_TO_INT(plate_solver(aargs));

	if (retval)
		siril_log_color_message(_("Image %s did not solve\n"), "salmon", root);
	return retval;
}

static int astrometry_finalize_hook(struct generic_seq_args *arg) {
	struct astrometry_data *aargs = (struct astrometry_data *)arg->user;
	if (aargs->cstars)
		free_fitted_stars(aargs->cstars);
	if (aargs->catalog_file)
		g_object_unref(aargs->catalog_file);
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
	seqargs->parallel = args->onlineCatalog != CAT_ASNET;		// for now crashes on Cancel if parallel is enabled for asnet on windows
#else
	seqargs->parallel = TRUE;
#endif
	seqargs->prepare_hook = astrometry_prepare_hook;
	seqargs->image_hook = astrometry_image_hook;
	seqargs->finalize_hook = astrometry_finalize_hook;
	seqargs->has_output = TRUE;
	seqargs->new_seq_prefix = strdup("ps_");
	seqargs->description = "plate solving";
	seqargs->user = args;

	siril_log_message(_("Running sequence plate solving using the %s catalogue\n"),
			catalog_to_str(args->onlineCatalog));
	start_in_new_thread(generic_sequence_worker, seqargs);
}

