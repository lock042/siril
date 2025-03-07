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

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/colors.h"
#include "algos/fitting.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "algos/statistics_float.h"
#include "algos/photometry.h"
#include "algos/spcc.h"
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
#include "algos/star_finder.h"
#include "algos/siril_wcs.h"
#include "io/single_image.h"
#include "io/image_format_fits.h" // For the datalink FITS functions
#include "io/local_catalogues.h"
#include "io/remote_catalogues.h"
#include "gui/siril_plot.h"
#include "gui/progress_and_log.h"
#include "gui/photometric_cc.h"
#include "registration/matching/misc.h" // for catalogue parsing helpers
#include "photometric_cc.h"

static const cmsCIEXYZ D65 = {0.95045471, 1.0, 1.08905029};
static const cmsCIEXYZ D50 = {0.964199999, 1.000000000, 0.824899998};

// Returns a valid xyY for 1667K <= t <= 25000K, otherwise xyY = { 0.0 }
// Uses Kim et al's cubic spline Planckian locus (https://en.wikipedia.org/wiki/Planckian_locus)
static void temp_to_xyY(cmsCIExyY *xyY, cmsFloat64Number t) {
	// Calculate x
	if (t < 1667.0)
		xyY->x = 0.0;
	else if (t < 4000.0)
		xyY->x = (-0.2661239e9 / (t * t * t)) - (0.2343589e6 / (t * t)) +( 0.8776956e3 / t) + 0.179910;
	else if (t < 25000.0)
		xyY->x = (-3.0258469e9 / (t * t * t)) + (2.1070379e6 / (t * t)) + (0.2226347e3 / t) + 0.240390;
	else
		xyY->x = 0.0;

	cmsFloat64Number x = xyY->x;
	// Calculate y
	if (t < 1667)
		xyY->y = 0.0;
	else if (t < 2222.0)
		xyY->y = (-1.1063814 * x * x * x) - (1.34811020 * x * x) + (2.18555832 * x) - 0.20219683;
	else if (t < 4000.0)
		xyY->y = (-0.9549476 * x * x * x) - (1.37418593 * x * x) + (2.09137015 * x) - 0.16748867;
	else if (t < 25000.0)
		xyY->y = (3.0817580 * x * x * x) - (5.87338670 * x * x) + (3.75112997 * x) - 0.37001483;
	else
		xyY->y = 0.0;

	if (!(xyY->x == 0.0 && xyY->y == 0.0))
		xyY->Y = 1.0;
	else
		xyY->Y = 0.0;
}

// Mitchell Charity's table of temperature to CIE (x,y) extends down to 1000K
// We don't use the whole table but only for temperatures < 1600K where the Kim
// splines don't work

// Returns a valid xyY for 1000K and up, otherwise xyY = { 0.0 }
// Uses Mitchell Charity's tabulation of black body xy values from
// http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html
// Cubic spline interpolation is used between each value; we only use this below
// 1650K approaching the region where the Kim splines don't work.
static const double tK[] = { 1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0 };
static const double x_1931_2deg_jv[] = { 0.6499,0.6361,0.6226,0.6095,0.5966,0.5841,0.572,0.5601 };
static const double y_1931_2deg_jv[] = { 0.3474,0.3594,0.3703,0.3801,0.3887,0.3962,0.4025,0.4076 };
#define NUM_POINTS_CHARITY 8

static void charity_temp_to_xyY(cmsCIExyY *xyY, cmsFloat64Number t) {
	if (t < 1000.0) {
		memset(xyY, 0.0, sizeof(cmsCIExyY));
		return;
	} else if (t > 1650.0) {
		siril_debug_print("Error, excessive temperature value %f passed to charity_temp_to_xyY\n", t);
		t = 1650.0;
	}
	gsl_interp *interp = gsl_interp_alloc(gsl_interp_cspline, NUM_POINTS_CHARITY);
	gsl_interp_init(interp, tK, x_1931_2deg_jv, NUM_POINTS_CHARITY);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	xyY->x = gsl_interp_eval(interp, tK, x_1931_2deg_jv, t, acc);
	gsl_interp_init(interp, tK, y_1931_2deg_jv, NUM_POINTS_CHARITY);
	gsl_interp_accel_reset(acc);
	xyY->y = gsl_interp_eval(interp, tK, y_1931_2deg_jv, t, acc);
	gsl_interp_free(interp);
	gsl_interp_accel_free(acc);
	xyY->Y = 1.0;
}

// Makes use of lcms2 to get the RGB values correct
// transform is calculated in get_white_balance_coeff below
// It provides the transform from XYZ to the required image colorspace
static void TempK2rgb(float *r, float *g, float *b, float TempK, cmsHTRANSFORM transform) { // RGB <0,1> <- BV <-0.4,+2.0> [-]
	cmsCIExyY WhitePoint;
	cmsCIEXYZ XYZ, XYZ_adapted;
	float xyz[3], rgb[3] = { 0.f };
	if (TempK > 1650.f)
		temp_to_xyY(&WhitePoint, TempK);
	else
		charity_temp_to_xyY(&WhitePoint, TempK);
	cmsxyY2XYZ(&XYZ, &WhitePoint);
	// Adapt the source (D65 Planckian locus) to the destination (lcms2 xyz profile with D50 whitepoint)
	cmsAdaptToIlluminant(&XYZ_adapted, &D65, &D50, &XYZ);
	xyz[0] = (float) XYZ_adapted.X;
	xyz[1] = (float) XYZ_adapted.Y;
	xyz[2] = (float) XYZ_adapted.Z;
	cmsDoTransform(transform, &xyz, &rgb, 1);
	cmsFloat64Number maxval = max(max(rgb[0], rgb[1]), rgb[2]);
	*r = rgb[0] / maxval;
	*g = rgb[1] / maxval;
	*b = rgb[2] / maxval;
}

int make_selection_around_a_star(cat_item *star, rectangle *area, fits *fit) {
	/* make a selection around the star, coordinates are in display reference frame */
	double fx = star->x, fy = star->y;
	double dx, dy;
	siril_to_display(fx, fy, &dx, &dy, fit->ry);

	double outer = com.pref.phot_set.outer;
	area->x = round_to_int(dx - outer);
	area->y = round_to_int(dy - outer);
	area->w = area->h = (int)ceil(outer * 2.0);

	/* Don't want stars too close to the edge */
	if (area->x <= 0 || area->x + area->w >= fit->rx - 1)
		return 1;
	if (area->y <= 0 || area->y + area->h >= fit->ry - 1)
		return 1;

	return 0;
}

static gchar *generate_title(const gchar *type, double arg, double br, double sig, gchar *wr, int nb_stars, int nb_excl, float *kw) {
	return g_strdup_printf(_("White Balance summary\n"
			"<span size=\"small\">"
			"%s SPCC Linear Fit: y = %f + %f·x, &#x03C3; = %f\n"
			"White reference: %s\n"
			"Number of stars: %d (removed %d outliers)\n"
			"White balance factors: %.3f %.3f %.3f"
			"</span>"),
			type, arg, br, sig, wr, nb_stars, nb_excl, kw[RLAYER], kw[GLAYER], kw[BLAYER]);
}

int filterArrays(double *x, double *y, int n) {
	int newSize = 0;

	for (int i = 0; i < n; i++) {
		if (x[i] != DBL_MAX && y[i] != DBL_MAX) {
			x[newSize] = x[i];
			y[newSize] = y[i];
			newSize++;
		}
	}
	return newSize;
}

int filtermaskArrays(double *x, double *y, gboolean *mask, int n) {
	int newSize = 0;

	for (int i = 0; i < n; i++) {
		if (mask[i]) {
			x[newSize] = x[i];
			y[newSize] = y[i];
			newSize++;
		}
	}
	return newSize;
}

static int get_spcc_white_balance_coeffs(struct photometric_cc_data *args, float* kw) {
	int nb_stars = args->nb_stars;
	fits *fit = args->fit;
	cat_item *stars = args->ref_stars->cat_items;
	double *irg = malloc(sizeof(double) * nb_stars);
	double *ibg = malloc(sizeof(double) * nb_stars);
	double *crg = malloc(sizeof(double) * nb_stars);
	double *cbg = malloc(sizeof(double) * nb_stars);
	double wrg = 0.f, wbg = 0.f;
	gboolean *maskrg = NULL, *maskbg = NULL;
	xpsampled response[3] = { init_xpsampled(), init_xpsampled(), init_xpsampled() };

	for (int k = 0 ; k < nb_stars; k++) {
		irg[k] = DBL_MAX;
		ibg[k] = DBL_MAX;
		crg[k] = DBL_MAX;
		cbg[k] = DBL_MAX;
	}
	gchar *str = ngettext("Applying aperture photometry to %d star.\n", "Applying aperture photometry to %d stars.\n", nb_stars);
	str = g_strdup_printf(str, nb_stars);
	siril_log_message(str);
	g_free(str);

	struct phot_config *ps = phot_set_adjusted_for_image(fit);
	siril_debug_print("aperture: %2.1f%s\tinner: %2.1f\touter: %2.1f\n", ps->aperture, ps->force_radius?"":" (auto)", ps->inner, ps->outer);
	gint ngood = 0, progress = 0;
	gint errors[PSF_ERR_MAX_VALUE] = { 0 };
	double minwl[3], maxwl[3];
	for (int chan = 0 ; chan < 3 ; chan++) {
		get_spectrum_from_args(args, &response[chan], chan);
		/* The idea here is that in narrowband mode we integrate the interpolated response (with no filtering
		 * included) over a very precise wavelength range, so as to get an accurate value. In broadband mode
		 * we include the effect of the filter in the response and we integrate over the full xp_sampled
		 * wavelength range. This principle is used in the flux and WB calcs too. */
		minwl[chan] = args->nb_mode ? args->nb_center[chan] - (args->nb_bandwidth[chan]/2) : XPSAMPLED_MIN_WL;
		maxwl[chan] = args->nb_mode ? args->nb_center[chan] + (args->nb_bandwidth[chan]/2) : XPSAMPLED_MAX_WL;
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) shared(progress, ngood)
#endif
	for (int i = 0; i < nb_stars; i++) {
		if (!get_thread_run())
			continue;
		rectangle area = { 0 };
		double flux[3] = { 0.0, 0.0, 0.0 };

		// Update the progress bar
		if (!(g_atomic_int_get(&progress) % 16))	// every 16 iterations
			set_progress_bar_data(NULL, (double) progress / (double) nb_stars);
		g_atomic_int_inc(&progress);

		// Make a selection 'area' around the ith star in the list of cat_items passed to the function
		if (make_selection_around_a_star(&stars[i], &area, fit)) {
			siril_debug_print("star %d is outside image or too close to border\n", i);
			g_atomic_int_inc(errors+PSF_ERR_OUT_OF_WINDOW);
			continue;
		}

		// Do photometry on each channel of the ith star; we compute the flux
		gboolean no_phot = FALSE;
		psf_error error = PSF_NO_ERR;
		for (int chan = 0; chan < 3 && !no_phot; chan ++) {
			// Photometry
			psf_star *photometry = psf_get_minimisation(fit, chan, &area, TRUE, ps, FALSE, com.pref.starfinder_conf.profile, &error);
			if (!photometry || !photometry->phot_is_valid || error != PSF_NO_ERR) {
				no_phot = TRUE;
			} else {
				// Flux calculation
				flux[chan] = pow(10., -0.4 * photometry->mag);
			}
			if (photometry)
				free_psf(photometry);
		}
		if (no_phot) {
			g_atomic_int_inc(errors+error);
			siril_debug_print("photometry failed for star %d, error %d\n", i, error);
			continue;
		}

		// Compute the image red/green rations from the ratios of r/b flux and b/g flux
		irg[i] = flux[RLAYER]/flux[GLAYER];
		ibg[i] = flux[BLAYER]/flux[GLAYER];

		// Get the Gaia DR3 xp_sampled spectrum from the downloaded datalink product
		double ref_flux[3];
		xpsampled star_spectrum = init_xpsampled();
//		switch (args->catalog) {
//			case CAT_GAIADR3_DIRECT:;
//			// Get the xp_sampled data from the datalink path
//				get_xpsampled(&star_spectrum, args->datalink_path, stars[i].index);
//				break;
//			case CAT_LOCAL_GAIA_XPSAMP:
//			default:;
//				memcpy(&star_spectrum.y, stars[i].flux, XPSAMPLED_LEN * sizeof(double));
		//				break;
//		}
		memcpy(&star_spectrum.y, stars[i].xp_sampled, XPSAMPLED_LEN * sizeof(double));

		// Convert flux to relative photon count normalized at 550nm
		flux_to_relcount(&star_spectrum);

		// Compute the expected response (product of catalogue flux and instrument response)
		xpsampled flux_expected = init_xpsampled();
		for (int chan = 0 ; chan < 3 ; chan++) {
			multiply_xpsampled(&flux_expected, &response[chan], &star_spectrum);
			ref_flux[chan] = integrate_xpsampled(&flux_expected, minwl[chan], maxwl[chan]);
		}

		// Compute the catalogue r/g and b/g ratios
		crg[i] = ref_flux[RLAYER]/ref_flux[GLAYER];
		cbg[i] = ref_flux[BLAYER]/ref_flux[GLAYER];

		// Remove any results with NaNs
		if (xisnanf(irg[i]) || xisnanf(ibg[i]) || xisnanf(crg[i]) || xisnanf(cbg[i])) {
			siril_debug_print("flux ratio NAN for star %d\n", i);
			irg[i] = DBL_MAX;
			ibg[i] = DBL_MAX;
			crg[i] = DBL_MAX;
			cbg[i] = DBL_MAX;
			g_atomic_int_inc(errors + PSF_ERR_FLUX_RATIO);
			continue;
		}
		g_atomic_int_inc(errors + PSF_NO_ERR);
		g_atomic_int_inc(&ngood);
	}
	free(ps);
	// Calculate white reference ratios
	double white_flux[3];
	xpsampled white_spectrum = init_xpsampled();
	GList *selected_white = g_list_nth(com.spcc_data.wb_ref, args->selected_white_ref);
	spcc_object *white = (spcc_object *) selected_white->data;
	load_spcc_object_arrays(white);
	init_xpsampled_from_library(&white_spectrum, white);
	spcc_object_free_arrays(white);
	xpsampled white_expected[3] = { init_xpsampled(), init_xpsampled(), init_xpsampled() };
	for (int chan = 0 ; chan < 3 ; chan++) {
		multiply_xpsampled(&white_expected[chan], &response[chan], &white_spectrum);
		white_flux[chan] = integrate_xpsampled(&white_expected[chan], minwl[chan], maxwl[chan]);
	}
	wrg = white_flux[RLAYER]/white_flux[GLAYER];
	wbg = white_flux[BLAYER]/white_flux[GLAYER];

	// Calculate effective primaries for later ** TODO: this doesn't work. To be revisited in a follow-on MR
	float white_flux_sum = white_flux[RLAYER] + white_flux[GLAYER] + white_flux[BLAYER];
	for (int chan = 0 ; chan < 3 ; chan++) {
		multiply_xpsampled_scalar(&white_expected[chan], 1.f / white_flux_sum);
	}
	args->primaries.Red = xpsampled_to_xyY(&response[RLAYER], com.pref.icc.cmf, minwl[RLAYER], maxwl[RLAYER]);
	args->primaries.Green = xpsampled_to_xyY(&response[GLAYER], com.pref.icc.cmf, minwl[GLAYER], maxwl[GLAYER]);
	args->primaries.Blue = xpsampled_to_xyY(&response[BLAYER], com.pref.icc.cmf, minwl[BLAYER], maxwl[BLAYER]);

	// Robust estimation of linear best fit
	double arg, brg, abg, bbg, deviation[2] = { 0.0 };
	// Remove any stars that failed photometry
	int n_rg = filterArrays(crg, irg, nb_stars);
	int n_bg = filterArrays(cbg, ibg, nb_stars);
	if (n_rg != n_bg) {
		free(irg);
		free(ibg);
		free(crg);
		free(cbg);
		siril_log_message(_("Array mismatch after discarding photometrically invalid stars\n"));
		return 1;
	}
	ngood = n_rg;
	int excl = nb_stars - ngood;
	str = ngettext("%d star excluded from the calculation\n", "%d stars excluded from the calculation\n", excl);
	str = g_strdup_printf(str, excl);
	siril_log_message(str);
	g_free(str);
	if (excl > 0)
		print_psf_error_summary(errors);
	if (ngood < 3) {
		free(irg);
		free(ibg);
		free(crg);
		free(cbg);
		siril_log_message(_("Error: insufficient photometrically valid stars\n"));
		return 1;
	}
	maskrg = calloc(ngood, sizeof(gboolean));
	if (robust_linear_fit(crg, irg, ngood, &arg, &brg, &deviation[0], maskrg)) {
		free(irg);
		free(ibg);
		free(crg);
		free(cbg);
		free(maskrg);
		siril_log_color_message(_("Error: unable to compute a fit to the data. "
				"Check your sensor and filter selections are correct.\n"), "red");
		return 1;
	}
	maskbg = calloc(ngood, sizeof(gboolean));
	if (robust_linear_fit(cbg, ibg, ngood, &abg, &bbg, &deviation[1], maskbg)) {
		free(irg);
		free(ibg);
		free(crg);
		free(cbg);
		free(maskrg);
		free(maskbg);
		siril_log_color_message(_("Error: unable to compute a fit to the data. "
				"Check your sensor and filter selections are correct.\n"), "red");
		return 1;
	}
	siril_log_color_message(_("SPCC Linear Fits\n"), "green");
	siril_log_message(_("Image R/G = %f + %f * Catalog R/G (sigma: %f)\n"), arg, brg, deviation[0]);
	siril_log_message(_("Image B/G = %f + %f * Catalog B/G (sigma: %f)\n"), abg, bbg, deviation[1]);
	double kr = 1.f / (arg + brg * wrg);
	double kg = 1.f;
	double kb = 1.f / (abg + bbg * wbg);
	double maxk = max(max(kr, kg), kb);
	kw[RLAYER] = kr / maxk;
	kw[GLAYER] = kg / maxk;
	kw[BLAYER] = kb / maxk;

	if (args->do_plot) {
		double stat_min, stat_max;
		int ngoodrg = filtermaskArrays(crg, irg, maskrg, ngood);
		gsl_stats_minmax(&stat_min, &stat_max, crg, 1, ngoodrg);
		double best_fit_rgx[2] = {stat_min, stat_max};
		double best_fit_rgy[2] = {arg + brg * best_fit_rgx[0], arg + brg * best_fit_rgx[1]};
		spcc_object *object = (spcc_object*) selected_white->data;

		siril_plot_data *spl_datarg = init_siril_plot_data();
		if (spl_datarg) {
			siril_plot_set_xlabel(spl_datarg, _("Catalog R/G (flux)"));
			siril_plot_set_savename(spl_datarg, "SPCC_RG_fit");
			gchar *title1 = generate_title("R/G", arg, brg, deviation[0], object->name, ngoodrg, ngood - ngoodrg, kw);
			siril_plot_set_title(spl_datarg, title1);
			g_free(title1);
			siril_plot_set_ylabel(spl_datarg, _("Image R/G (flux)"));
			siril_plot_add_xydata(spl_datarg, _("R/G"), ngoodrg, crg, irg, NULL, NULL);
			siril_plot_add_xydata(spl_datarg, _("Best fit"), 2, best_fit_rgx, best_fit_rgy, NULL, NULL);
			siril_plot_set_nth_plot_type(spl_datarg, 1, KPLOT_POINTS);
			siril_plot_set_nth_plot_type(spl_datarg, 2, KPLOT_LINES);
			siril_plot_set_yfmt(spl_datarg, "%.1lf");
			spl_datarg->cfgdata.point.radius = 1;
			spl_datarg->cfgdata.point.sz = 2;
			spl_datarg->cfgdata.line.sz = 2;
			siril_add_idle(create_new_siril_plot_window, spl_datarg);
		}

		int ngoodbg = filtermaskArrays(cbg, ibg, maskbg, ngood);
		gsl_stats_minmax(&stat_min, &stat_max, cbg, 1, ngoodbg);
		double best_fit_bgx[2] = {stat_min, stat_max};
		double best_fit_bgy[2] = {abg + bbg * best_fit_bgx[0], abg + bbg * best_fit_bgx[1]};
		siril_plot_data *spl_databg = init_siril_plot_data();
		if (spl_databg) {
			siril_plot_set_xlabel(spl_databg, _("Catalog B/G (flux)"));
			siril_plot_set_savename(spl_databg, "SPCC_BG_fit");
			gchar *title2 = generate_title("B/G", abg, bbg, deviation[1], object->name, ngoodbg, ngood - ngoodbg, kw);
			siril_plot_set_title(spl_databg, title2);
			g_free(title2);
			siril_plot_set_ylabel(spl_databg, _("Image B/G (flux)"));
			gchar *spl_legendbg = _("B/G");
			siril_plot_add_xydata(spl_databg, spl_legendbg, ngoodbg, cbg, ibg, NULL, NULL);
			siril_plot_add_xydata(spl_databg, _("Best fit"), 2, best_fit_bgx, best_fit_bgy, NULL, NULL);
			siril_plot_set_nth_plot_type(spl_databg, 1, KPLOT_POINTS);
			siril_plot_set_nth_plot_type(spl_databg, 2, KPLOT_LINES);
			siril_plot_set_yfmt(spl_databg, "%.1lf");
			spl_databg->cfgdata.point.radius = 1;
			spl_databg->cfgdata.point.sz = 2;
			spl_databg->cfgdata.line.sz = 2;
			siril_add_idle(create_new_siril_plot_window, spl_databg);
		}

		siril_add_idle(end_generic, NULL);
	}
	free(irg);
	free(ibg);
	free(crg);
	free(cbg);
	free(maskrg);
	free(maskbg);
	if (kw[RLAYER] < 0.f || kw[GLAYER] < 0.f || kw[BLAYER] < 0.f) {
		siril_log_color_message(_("Error calculating white balance coefficients: kw contains negative values.\n"),"red");
		return 1;
	}
	siril_log_message(_("Found a solution for color calibration using %d stars. Factors:\n"), ngood);
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3lf\n", chan, kw[chan]);
	}

	if (ngood < 20)
		siril_log_color_message(_("The photometric color calibration has found a solution which may not be perfect because it did not rely on many stars\n"), ngood < 5 ? "red" : "salmon");
	else if (deviation[0] > 0.1 || deviation[1] > 0.1)
		siril_log_message(_("The photometric color calibration seems to have found an imprecise solution, consider correcting the image gradient first\n"));

	return 0;
}

static int get_pcc_white_balance_coeffs(struct photometric_cc_data *args, float *kw) {
	int nb_stars = args->nb_stars;
	fits *fit = args->fit;
	cat_item *stars = args->ref_stars->cat_items;
	float *data[3];
	data[RLAYER] = malloc(sizeof(float) * nb_stars);
	data[GLAYER] = malloc(sizeof(float) * nb_stars);
	data[BLAYER] = malloc(sizeof(float) * nb_stars);
	if (!data[RLAYER] || !data[GLAYER] || !data[BLAYER]) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	for (int k = 0; k < nb_stars; k++) {
		data[RLAYER][k] = FLT_MAX;
		data[GLAYER][k] = FLT_MAX;
		data[BLAYER][k] = FLT_MAX;
	}
	gchar *str = ngettext("Applying aperture photometry to %d star.\n", "Applying aperture photometry to %d stars.\n", nb_stars);
	str = g_strdup_printf(str, nb_stars);
	siril_log_message(str);
	g_free(str);

	struct phot_config *ps = phot_set_adjusted_for_image(fit);
	siril_debug_print("aperture: %2.1f%s\tinner: %2.1f\touter: %2.1f\n", ps->aperture, ps->force_radius?"":" (auto)", ps->inner, ps->outer);
	gint ngood = 0, progress = 0;
	gint errors[PSF_ERR_MAX_VALUE] = { 0 };

	cmsHPROFILE xyzprofile = NULL;
	cmsHPROFILE profile; // This is initialised in either the if or else branch
	cmsHTRANSFORM transform = NULL;

	// This transform is for normal PCC to transform the reference
	//star XYZ to RGB: the SPCC source->working transform is dealt
	//with later.
	xyzprofile = cmsCreateXYZProfile();
	if (fit->icc_profile) {
		profile = copyICCProfile(fit->icc_profile);
		if (!fit_icc_is_linear(fit)) {
			siril_log_color_message(_("Image color space is nonlinear. It is recommended to "
					"apply photometric color calibration to linear images.\n"), "salmon");
		}
	} else {
		profile = siril_color_profile_linear_from_color_profile(com.icc.working_standard);
		fit->icc_profile = copyICCProfile(profile);
		color_manage(fit, (fit->icc_profile != NULL));
	}
	if (xyzprofile && profile) {
		transform = cmsCreateTransformTHR(com.icc.context_single,
					xyzprofile,
					TYPE_XYZ_FLT,
					profile,
					TYPE_RGB_FLT,
					INTENT_RELATIVE_COLORIMETRIC,
					cmsFLAGS_NONEGATIVES);
		cmsCloseProfile(profile);
		cmsCloseProfile(xyzprofile);
	}
	if (!transform) {
		siril_log_color_message(_("Error: failed to set up colorspace transform.\n"), "red");
		free(ps);
		return 1;
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) shared(progress, ngood)
#endif
	for (int i = 0; i < nb_stars; i++) {
		if (!get_thread_run())
			continue;
		rectangle area = { 0 };
		float flux[3] = { 0.f, 0.f, 0.f };
		float r, g, b, bv;
		if (!(g_atomic_int_get(&progress) % 16))	// every 16 iterations
			set_progress_bar_data(NULL, (double) progress / (double) nb_stars);
		g_atomic_int_inc(&progress);

		if (make_selection_around_a_star(&stars[i], &area, fit)) {
			siril_debug_print("star %d is outside image or too close to border\n", i);
			g_atomic_int_inc(errors+PSF_ERR_OUT_OF_WINDOW);
			continue;
		}

		gboolean no_phot = FALSE;
		psf_error error = PSF_NO_ERR;
		for (int chan = 0; chan < 3 && !no_phot; chan ++) {
			psf_star *photometry = psf_get_minimisation(fit, chan, &area, TRUE, ps, FALSE, com.pref.starfinder_conf.profile, &error);
			if (!photometry || !photometry->phot_is_valid || error != PSF_NO_ERR)
				no_phot = TRUE;
			else flux[chan] = powf(10.f, -0.4f * (float) photometry->mag);
			if (photometry)
				free_psf(photometry);
		}
		if (no_phot) {
			g_atomic_int_inc(errors + error);
			siril_debug_print("photometry failed for star %d, error %d\n", i, error);
			continue;
		}
		// get r g b coefficient
		// If the Gaia Teff field is populated (CAT_GAIADR3, CAT_LOCAL_GAIA_ASTRO and
		// CAT_GAIADR3_DIRECT), use that as it should be more accurate.
		// Otherwise, we convert from Johnson B-V
		// Note, with CAT_LOCAL_GAIA_ASTRO we have already excluded stars without Teff
		// in get_raw_stars_from_local_gaia_astro_catalogue()
		if (stars[i].teff < 0.5f) {
			bv = min(max(stars[i].BV, -0.4f), 2.f);
			cmsFloat64Number TempK = BV_to_T(bv);
			stars[i].teff = (float) TempK;
		}
		// TempK2rgb converts the temperature to RGB values in the working colorspace
		TempK2rgb(&r, &g, &b, stars[i].teff, transform);
		/* get Color calibration factors for current star */

		data[RLAYER][i] = (1.f / flux[RLAYER]) * r;
		data[GLAYER][i] = (1.f / flux[GLAYER]) * g;
		data[BLAYER][i] = (1.f / flux[BLAYER]) * b;

		if (xisnanf(data[RLAYER][i]) || xisnanf(data[GLAYER][i]) || xisnanf(data[BLAYER][i])) {
			data[RLAYER][i] = FLT_MAX;
			data[GLAYER][i] = FLT_MAX;
			data[BLAYER][i] = FLT_MAX;
			siril_debug_print("flux ratio NAN for star %d\n", i);
			g_atomic_int_inc(errors + PSF_ERR_FLUX_RATIO);
			continue;
		}
		g_atomic_int_inc(errors + PSF_NO_ERR);
		g_atomic_int_inc(&ngood);
	}
	if (transform)
		cmsDeleteTransform(transform);
	free(ps);
	if (!get_thread_run()) {
		free(data[RLAYER]);
		free(data[GLAYER]);
		free(data[BLAYER]);
		return 1;
	}
	int excl = nb_stars - ngood;
	str = ngettext("%d star excluded from the calculation\n", "%d stars excluded from the calculation\n", excl);
	str = g_strdup_printf(str, excl);
	siril_log_message(str);
	g_free(str);
	if (excl > 0)
		print_psf_error_summary(errors);

	if (ngood == 0) {
		siril_log_message(_("No valid stars found.\n"));
		free(data[RLAYER]);
		free(data[GLAYER]);
		free(data[BLAYER]);
		return 1;
	}
	/* sort in ascending order before using siril_stats_mean_from_linearFit
	 * Hence, DBL_MAX are at the end of the tab */
	quicksort_f(data[RLAYER], nb_stars);
	quicksort_f(data[GLAYER], nb_stars);
	quicksort_f(data[BLAYER], nb_stars);

	double deviation[3];
	/* we do not take into account FLT_MAX values */
	kw[RLAYER] = siril_stats_robust_mean(data[RLAYER], 1, ngood, &(deviation[RLAYER]));
	kw[GLAYER] = siril_stats_robust_mean(data[GLAYER], 1, ngood, &(deviation[GLAYER]));
	kw[BLAYER] = siril_stats_robust_mean(data[BLAYER], 1, ngood, &(deviation[BLAYER]));
	if (kw[RLAYER] < 0.f || kw[GLAYER] < 0.f || kw[BLAYER] < 0.f) {
		free(data[RLAYER]);
		free(data[GLAYER]);
		free(data[BLAYER]);
		return 1;
	}
	/* normalize factors */
	double maxk = max(max(kw[RLAYER], kw[GLAYER]), kw[BLAYER]);
	kw[RLAYER] /= maxk;
	kw[GLAYER] /= maxk;
	kw[BLAYER] /= maxk;
	siril_log_message(_("Found a solution for color calibration using %d stars. Factors:\n"), ngood);
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3lf\t(deviation: %.3f)\n", chan, kw[chan], deviation[chan]);
	}
	if (ngood < 20)
		siril_log_color_message(_("The photometric color calibration has found a solution which may not be perfect because it did not rely on many stars\n"), ngood < 5 ? "red" : "salmon");
	else if (deviation[RLAYER] > 0.1 || deviation[GLAYER] > 0.1 || deviation[BLAYER] > 0.1)
		siril_log_message(_("The photometric color calibration seems to have found an imprecise solution, consider correcting the image gradient first\n"));
	free(data[RLAYER]);
	free(data[GLAYER]);
	free(data[BLAYER]);
	return 0;
}

/*
Gets bg, min and max values per channel and sets the chennel with middle bg value
*/
int get_stats_coefficients(fits *fit, rectangle *area, float *bg, float t0, float t1) {
	// we cannot use compute_all_channels_statistics_single_image because of the area
	siril_log_message(_("Computing background reference with tolerance +%.2fσ / -%.2fσ.\n"), t1, -t0);
	for (int chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, area, STATS_BASIC | STATS_MAD, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		// Compute the robust median as the median of all points within t{0,1} * MAD of the channel median
		double median_c = stat->median;
		double sigma_c = MAD_NORM * stat->mad;
		double lower_c = median_c + t0 * sigma_c;
		double upper_c = median_c + t1 * sigma_c;
		double robust_median;
		// The robust median is calculated only using pixel values between lower_c and upper_c
		if (fit->type == DATA_USHORT) {
			robust_median = robust_median_w(fit, area, chan, (float) lower_c, (float) upper_c);
		} else {
			robust_median = robust_median_f(fit, area, chan, (float) lower_c, (float) upper_c);
		}
		bg[chan] = robust_median;

		free_stats(stat);
	}

	return 0;
}

int apply_photometric_color_correction(fits *fit, const float *kw, const float *bg) {
	float offset[3];
	float bg_mean = (bg[RLAYER] + bg[GLAYER] + bg[BLAYER]) / 3.f;

	for (int chan = 0; chan < 3; chan++) {
		if (isnan(kw[chan])) { // Check for NaN... If they are NaN the result image is junk
			siril_log_color_message(_("Error computing coefficients: aborting...\n"), "red");
			return 1;
		}
		offset[chan] = (-bg[chan] * kw[chan] + bg_mean);
	}
	siril_log_message("After renormalization, the following coefficients are applied\n");
	siril_log_color_message(_("White balance factors:\n"), "green");
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3f\n", chan, kw[chan]);
	}
	siril_log_color_message(_("Background reference:\n"), "green");
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("B%d: %+.5e\n", chan, offset[chan]);
	}

	for (int chan = 0; chan < 3; chan++) {
		size_t n = fit->naxes[0] * fit->naxes[1];
		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[chan];
			for (size_t i = 0; i < n; ++i) {
				buf[i] = roundf_to_WORD((float)buf[i] * kw[chan] + offset[chan]);
			}
		}
		else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[chan];
			for (size_t i = 0; i < n; ++i) {
				buf[i] = buf[i] * kw[chan] + offset[chan];
			}
		}
		else return 1;
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/* run the PCC using the existing star list of the image from the provided file */
int photometric_cc(struct photometric_cc_data *args) {
	float kw[3];
	float bg[3];

	// Moved here from gui/photometric_cc.c so it always applies. Once the command always calls photometric_cc_standalone()
	// it can move there, rather than after the catalog download
	if (args->fit->keywords.wcslib->lin.dispre == NULL) {
		siril_log_color_message(_("Found linear plate solve data. For better result you should redo platesolving\n"), "salmon");
	}

	if (!isrgb(args->fit)) {
		siril_log_message(_("Photometric color correction will do nothing for monochrome images\n"));
		return 0;
	}

	rectangle *bkg_sel = NULL;
	if (!args->bg_auto)
		bkg_sel = &(args->bg_area);

	/* we use the median of each channel to sort them by level and select
	 * the reference channel expressed in terms of order of middle median value */
	if (get_stats_coefficients(args->fit, bkg_sel, bg, args->t0, args->t1)) {
		siril_log_message(_("failed to compute statistics on image, aborting\n"));
		free(args);
		return 1;
	}

	/* set photometry parameters to values adapted to the image */
	struct phot_config backup = com.pref.phot_set;
	com.pref.phot_set.force_radius = FALSE;
	com.pref.phot_set.inner = max(7.0, 3.0 * args->fwhm);
	com.pref.phot_set.outer = com.pref.phot_set.inner + 10;
	siril_log_message(_("Photometry radii set to %.1f for inner and %.1f for outer\n"),
			com.pref.phot_set.inner, com.pref.phot_set.outer);

	set_progress_bar_data(_("Photometric color calibration in progress..."), PROGRESS_RESET);
	int ret;
	if (args->spcc)
		ret = get_spcc_white_balance_coeffs(args, kw);
	else
		ret = get_pcc_white_balance_coeffs(args, kw);

	if (!ret) {
		ret = apply_photometric_color_correction(args->fit, kw, bg);
		if (!ret) {
/*
 *	This is temporarily removed pending fixing the source profile calc in a separate MR
			if (args->spcc) {
				ret = spcc_set_source_profile(args);
			}
*/
			invalidate_stats_from_fit(args->fit);
		}
	} else {
		set_progress_bar_data(_("Photometric Color Calibration failed"), PROGRESS_DONE);
	}

	com.pref.phot_set = backup;
	free(args);
	return ret;
}

/* photometric_cc is the entry point for the PCC following the call to the
 * plate solver which gives the star list. We can also run the PCC on a
 * plate-solved image without running plate solving again, this is what this
 * function does.
 */
gpointer photometric_cc_standalone(gpointer p) {
	struct photometric_cc_data *args = (struct photometric_cc_data *)p;

	if (!has_wcs(args->fit)) {
		siril_log_color_message(_("Cannot run the standalone photometric color calibration on this image because it has no WCS data or it is not supported\n"), "red");
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	/* run peaker to measure FWHM of the image to adjust photometry settings */
	args->fwhm = measure_image_FWHM(args->fit, -1);
	if (args->fwhm <= 0.0f) {
		siril_log_message(_("Error computing FWHM for photometry settings adjustment\n"));
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	/* get stars from a photometric catalog */
	double ra, dec;
	center2wcs(args->fit, &ra, &dec);
	double resolution = get_wcs_image_resolution(args->fit);
	if ((ra == -1.0 && dec == -1.0) || resolution <= 0.0) {
		siril_log_color_message(_("Cannot run the standalone photometric color calibration on this image because it has no WCS data or it is not supported\n"), "red");
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	int nb_stars = 0;
	gboolean image_is_gfit = args->fit == &gfit;

	uint64_t sqr_radius = ((uint64_t) gfit.rx * (uint64_t) gfit.rx + (uint64_t) gfit.ry * (uint64_t) gfit.ry) / 4;
	double radius = resolution * sqrt((double)sqr_radius);	// in degrees
	double mag = args->mag_mode == LIMIT_MAG_ABSOLUTE ?
		args->magnitude_arg : compute_mag_limit_from_position_and_fov(ra, dec, radius * 2.0, BRIGHTEST_STARS * 2); // factor 2 on brightest_stars to account for the fact not all stars will be within image
	if (args->mag_mode == LIMIT_MAG_AUTO_WITH_OFFSET)
		mag += args->magnitude_arg;

	int retval = 0;
	if (args->catalog == CAT_LOCAL_KSTARS || args->catalog == CAT_LOCAL_GAIA_ASTRO || args->catalog == CAT_LOCAL_GAIA_XPSAMP) {
		siril_log_message(_("Getting stars from local catalogue %s for %s, with a radius of %.2f degrees and limit magnitude %.2f\n"), catalog_to_str(args->catalog), args->spcc ? _("SPCC") : _("PCC"), radius * 2.0,  mag);
	} else {
		switch (args->catalog) {
			case CAT_GAIADR3:
				mag = min(mag, 18.0);
				break;
			case CAT_GAIADR3_DIRECT:
				mag = min(mag, 17.6);	// most Gaia XP_SAMPLED spectra are for mag < 17.6
				break;
			case CAT_APASS:
				mag = min(mag, 17.0);	// in APASS, B is available for V < 17
				break;
			case CAT_NOMAD:
				mag = min(mag, 18.0);	// in NOMAD, B is available for V < 18
				break;
			default:
				siril_log_color_message(_("No valid catalog found.\n"), "red");
				return GINT_TO_POINTER(1);
		}
		siril_log_message(_("Getting stars from online catalogue %s for PCC, with a radius of %.2f degrees and limit magnitude %.2f\n"), catalog_to_str(args->catalog),radius * 2.0,  mag);
	}

	// preparing the catalogue query
	siril_catalogue *siril_cat = siril_catalog_fill_from_fit(args->fit, args->catalog, mag);
	// don't set phot if we are using GAIADR3, we use this catalog in a different way
	siril_cat->phot = !(siril_cat->cat_index == CAT_GAIADR3_DIRECT);

	/* Fetching the catalog*/
	if (args->spcc && siril_cat->cat_index == CAT_GAIADR3_DIRECT) {
		retval = siril_gaiadr3_datalink_query(siril_cat, XP_SAMPLED, &args->datalink_path, 5000);
		for (int i = 0 ; i < siril_cat->nbitems ; i++) {
			// Read the xp_sampled data from the RAW-structured FITS returned from Gaia datalink
			siril_cat->cat_items[i].xp_sampled = malloc(XPSAMPLED_LEN * sizeof(double));
			get_xpsampled(siril_cat->cat_items[i].xp_sampled, args->datalink_path, i);
		}
	} else if (siril_catalog_conesearch(siril_cat) <= 0) {
		retval = 1;
	}
	// At this point siril_cat contains an array of cat_items (with with xp_sampled populated if doing SPCC)
	nb_stars = siril_cat->nbitems;

	/* project using WCS */
	siril_catalog_project_with_WCS(siril_cat, args->fit, TRUE, FALSE);
	retval |= nb_stars == 0;

	gboolean spcc = args->spcc; // Needed for the success message after args has been freed in photometric_cc()
	if (!retval) {
		if (!com.script) {
			const gchar *algo = args->spcc ? _("SPCC") : _("PCC");
			undo_save_state(args->fit, _("Photometric CC (algorithm: %s)"), algo);
		}
		args->ref_stars = siril_cat;
		args->nb_stars = nb_stars;
		retval = photometric_cc(args);	// args is freed from here
	} else {
		free(args);
		siril_log_color_message(_("Catalog error, no stars identified!\n"), "red");
	}
	args = NULL;
	siril_catalog_free(siril_cat);

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);

	if (!retval && image_is_gfit) {
		if (spcc)
			siril_log_color_message(_("Spectrophotometric Color Calibration succeeded.\n"), "green");
		else
			siril_log_color_message(_("Photometric Color Calibration succeeded.\n"), "green");
		notify_gfit_modified();
	}
	else siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}
