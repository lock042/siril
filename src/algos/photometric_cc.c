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

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_interp.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/colors.h"
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
#include "gui/progress_and_log.h"
#include "gui/photometric_cc.h"
#include "registration/matching/misc.h" // for catalogue parsing helpers
#include "photometric_cc.h"

static gboolean spcc_filters_initialized = FALSE;
const double xpsampled_wl[163] = {378, 380, 382, 384, 386, 388, 390, 392, 394, 396, 398, 400, 402, 404, 406, 408, 410, 412, 414, 416, 418, 420, 422, 424, 426, 428, 430, 432, 434, 436, 438, 440, 442, 444, 446, 448, 450, 452, 454, 456, 458, 460,
  462, 464, 466, 468, 470, 472, 474, 476, 478, 480, 482, 484, 486, 488, 490, 492, 494, 496, 498, 500, 502, 504, 506, 508, 510, 512, 514, 516, 518, 520, 522, 524, 526, 528, 530, 532, 534, 536, 538, 540, 542, 544, 546, 548, 550, 552,
  554, 556, 558, 560, 562, 564, 566, 568, 570, 572, 574, 576, 578, 580, 582, 584, 586, 588, 590, 592, 594, 596, 598, 600, 602, 604, 606, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 628, 630, 632, 634, 636, 638, 640, 642, 644,
  646, 648, 650, 652, 654, 656, 658, 660, 662, 664, 666, 668, 670, 672, 674, 676, 678, 680, 682, 684, 686, 688, 690, 692, 694, 696, 698, 700, 702};

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

// Alternative temp_to_xyY converters, for comparison with Kim et al above
// Only one should be kept once this MR is complete

// Returns a valid xyY for 1000K and up, otherwise xyY = { 0.0 }
// Uses Mitchell Charity's tabulation of black body xy values from
// http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html
// Linear interpolation is used between each value.
// Commented out until / unless needed
/*
static void charity_temp_to_xyY(cmsCIExyY *xyY, cmsFloat64Number t) {
	int i = 0;
	if (t < 1000.0) {
		memset(xyY, 0.0, sizeof(cmsCIExyY));
		return;
	} else if (t > 40000.0)
		t = 40000.0;
	while (t > tK[i+1]) {
		i++;
	}
	float t1 = tK[i];
	float t2 = tK[i+1];
	float x1 = x_1931_2deg_jv[i];
	float x2 = x_1931_2deg_jv[i+1];
	float y1 = y_1931_2deg_jv[i];
	float y2 = y_1931_2deg_jv[i+1];
	xyY->x = (cmsFloat64Number) x1 + ((t - t1) / (t2 - t1)) * (x2 - x1);
	xyY->y = (cmsFloat64Number) y1 + ((t - t1) / (t2 - t1)) * (y2 - y1);
	xyY->Y = 1.0;
}
*/

// Returns a valid xyY for 2000K <= t, otherwise xyY = { 0.0 }
// Uses BQ Octantis's 6th order best fit for D65 values down to 2000K
// (https://www.cloudynights.com/topic/849382-generating-a-planckian-ccm/)
// Commented out until / unless needed
/*
static void BQ_temp_to_xyY(cmsCIExyY *xyY, cmsFloat64Number t) {
	if (t < 2000.0) {
		xyY->x = xyY->y = xyY->Y = 0.0;
	} else {
		xyY->x = -1.3737e-25 * pow(t,6.0) + 3.0985e-22 * pow(t, 5.0) + 1.5842e-16 * pow(t, 4.0)
				- 4.0138e-12 * pow(t, 3.0) + 4.3777e-08 * pow(t, 2.0) - 2.4363e-04 * t + 8.7309e-01;
		cmsFloat64Number x = xyY->x;
		xyY->y =  2.5110e+01 * pow(x, 5.0) - 5.9883e+01 * pow(x, 4.0) + 5.5545e+01 * pow(x, 3.0)
				- 2.7667e+01 * pow(x, 2.0) + 8.1167e+00 * x - 7.0462e-01;
		xyY->Y = 1.0;
	}
}
*/

// Makes use of lcms2 to get the RGB values correct
// transform is calculated in get_white_balance_coeff below
// It provides the transform from XYZ to the required image colorspace
static void TempK2rgb(float *r, float *g, float *b, float TempK, cmsHTRANSFORM transform) { // RGB <0,1> <- BV <-0.4,+2.0> [-]
	cmsCIExyY WhitePoint;
	cmsCIEXYZ XYZ, XYZ_adapted;
	float xyz[3], rgb[3] = { 0.f };
	temp_to_xyY(&WhitePoint, TempK);
//	charity_temp_to_xyY(&WhitePoint, TempK);
//	BQ_temp_to_xyY(&WhitePoint, TempK);
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

static int make_selection_around_a_star(pcc_star star, rectangle *area, fits *fit) {
	/* make a selection around the star, coordinates are in display reference frame */
	double fx = star.x, fy = star.y;
	double dx, dy;
	fits_to_display(fx, fy, &dx, &dy, fit->ry);

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

static int get_white_balance_coeff(pcc_star *stars, int nb_stars, fits *fit, float *kw, int norm_channel, struct photometric_cc_data *args) {
	float *data[3];
	data[RED] = malloc(sizeof(float) * nb_stars);
	data[GREEN] = malloc(sizeof(float) * nb_stars);
	data[BLUE] = malloc(sizeof(float) * nb_stars);
	if (!data[RED] || !data[GREEN] || !data[BLUE]) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	for (int k = 0; k < nb_stars; k++) {
		data[RED][k] = FLT_MAX;
		data[GREEN][k] = FLT_MAX;
		data[BLUE][k] = FLT_MAX;
	}
	gchar *str = ngettext("Applying aperture photometry to %d star.\n", "Applying aperture photometry to %d stars.\n", nb_stars);
	str = g_strdup_printf(str, nb_stars);
	siril_log_message(str);
	g_free(str);

	struct phot_config *ps = phot_set_adjusted_for_image(fit);
	siril_debug_print("aperture: %2.1f%s\tinner: %2.1f\touter: %2.1f\n", ps->aperture, ps->force_radius?"":" (auto)", ps->inner, ps->outer);
	gint ngood = 0, progress = 0;
	gint errors[PSF_ERR_MAX_VALUE] = { 0 };

	xpsampled response[3] = { { xpsampled_wl, { 0.0 } }, { xpsampled_wl, { 0.0 } }, { xpsampled_wl, { 0.0 } } };

	cmsHPROFILE xyzprofile = NULL;
	cmsHPROFILE profile = NULL;
	cmsHTRANSFORM transform = NULL;

	if (args->spcc) {
		for (int chan = 0 ; chan < 3 ; chan++) {
			get_spectrum_from_args(args, &response[chan], chan);
		}
		// Calculate effective primaries for later
		args->primaries.Red = xpsampled_to_xyY(&response[RED], CMF_1931);
		args->primaries.Green = xpsampled_to_xyY(&response[GREEN], CMF_1931);
		args->primaries.Blue = xpsampled_to_xyY(&response[BLUE], CMF_1931);
	} else {
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
			return 1;
		}
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

		if (make_selection_around_a_star(stars[i], &area, fit)) {
			siril_debug_print("star %d is outside image or too close to border\n", i);
			g_atomic_int_inc(errors+PSF_ERR_OUT_OF_WINDOW);
			continue;
		}

		gboolean no_phot = FALSE;
		psf_error error = PSF_NO_ERR;
		for (int chan = 0; chan < 3 && !no_phot; chan ++) {
			psf_star *photometry = psf_get_minimisation(fit, chan, &area, TRUE, ps, FALSE, com.pref.starfinder_conf.profile, &error);
			g_atomic_int_inc(errors+error);
			if (!photometry || !photometry->phot_is_valid || error != PSF_NO_ERR)
				no_phot = TRUE;
			else flux[chan] = powf(10.f, -0.4f * (float) photometry->mag);
			if (photometry)
				free_psf(photometry);
		}
		if (no_phot) {
			siril_debug_print("photometry failed for star %d, error %d\n", i, error);
			continue;
		}
		if (!args->spcc) {
			// get r g b coefficient
			// If the Gaia Teff field is populated (CAT_GAIADR3 and
			// CAT_GAIADR3_DIRECT), use that as it should be more accurate.
			// Otherwise, we convert from Johnson B-V
			if (stars[i].teff == 0.f) {
				bv = min(max(stars[i].BV, -0.4f), 2.f);
				cmsFloat64Number TempK = BV_to_T(bv);
				stars[i].teff = (float) TempK;
			}
			// TempK2rgb converts the temperature to RGB values in the working colorspace
			TempK2rgb(&r, &g, &b, stars[i].teff, transform);
			/* get Color calibration factors for current star */
		} else {
			float ref_flux[3];
			// Get the XP_sampled spectrum for this star
			xpsampled xp_sampled = { xpsampled_wl, { 0.0 } };
			get_xpsampled(&xp_sampled, args->datalink_path, stars[i].index);

			// Multiply the stellar spectrum by the channel response and integrate
			xpsampled flux_expected = { xpsampled_wl, { 0.0 } };
			for (int chan = 0 ; chan < 3 ; chan++) {
				multiply_xpsampled(&flux_expected, &response[chan], &xp_sampled);
				ref_flux[chan] = integrate_xpsampled(&flux_expected);
			}
			r = ref_flux[RED]; g = ref_flux[GREEN]; b = ref_flux[BLUE];
		}
		data[RED][i] = (flux[norm_channel] / flux[RED]) * r;
		data[GREEN][i] = (flux[norm_channel] / flux[GREEN]) * g;
		data[BLUE][i] = (flux[norm_channel] / flux[BLUE]) * b;

		if (xisnanf(data[RED][i]) || xisnanf(data[GREEN][i]) || xisnanf(data[BLUE][i])) {
			data[RED][i] = FLT_MAX;
			data[GREEN][i] = FLT_MAX;
			data[BLUE][i] = FLT_MAX;
			continue;
		}
		g_atomic_int_inc(&ngood);
	}
	if (transform)
		cmsDeleteTransform(transform);
	free(ps);
	if (!get_thread_run()) {
		free(data[RED]);
		free(data[GREEN]);
		free(data[BLUE]);
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
		free(data[RED]);
		free(data[GREEN]);
		free(data[BLUE]);
		return 1;
	}
	/* sort in ascending order before using siril_stats_mean_from_linearFit
	 * Hence, DBL_MAX are at the end of the tab */
	quicksort_f(data[RED], nb_stars);
	quicksort_f(data[GREEN], nb_stars);
	quicksort_f(data[BLUE], nb_stars);

	double deviation[3];
	/* we do not take into account FLT_MAX values */
	kw[RED] = siril_stats_robust_mean(data[RED], 1, ngood, &(deviation[RED]));
	kw[GREEN] = siril_stats_robust_mean(data[GREEN], 1, ngood, &(deviation[GREEN]));
	kw[BLUE] = siril_stats_robust_mean(data[BLUE], 1, ngood, &(deviation[BLUE]));
	if (kw[RED] < 0.f || kw[GREEN] < 0.f || kw[BLUE] < 0.f) {
		free(data[RED]);
		free(data[GREEN]);
		free(data[BLUE]);
		return 1;
	}

	/* normalize factors */
	float norm = kw[norm_channel];
	kw[RED] /= norm;
	kw[GREEN] /= norm;
	kw[BLUE] /= norm;
	siril_log_message(_("Found a solution for color calibration using %d stars. Factors:\n"), ngood);
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3lf\t(deviation: %.3f)\n", chan, kw[chan], deviation[chan]);
	}

	if (ngood < 20)
		siril_log_color_message(_("The photometric color calibration has found a solution which may not be perfect because it did not rely on many stars\n"), ngood < 5 ? "red" : "salmon");
	else if (deviation[RED] > 0.1 || deviation[GREEN] > 0.1 || deviation[BLUE] > 0.1)
		siril_log_message(_("The photometric color calibration seems to have found an imprecise solution, consider correcting the image gradient first\n"));
	free(data[RED]);
	free(data[GREEN]);
	free(data[BLUE]);
	return 0;
}

static int cmp_coeff(const void *a, const void *b) {
	const coeff *a1 = (const coeff *) a;
	const coeff *a2 = (const coeff *) b;
	if (a1->value > a2->value)
		return 1;
	if (a1->value < a2->value)
		return -1;
	return 0;
}

/*
Gets bg, min and max values per channel and sets the chennel with middle bg value
*/
int get_stats_coefficients(fits *fit, rectangle *area, coeff *bg, float *mins, float *maxs, int *norm_channel) {
	// we cannot use compute_all_channels_statistics_single_image because of the area
	for (int chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, area, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg[chan].value = stat->median;
		bg[chan].channel = chan;
		if (!area) {
			mins[chan] = stat->min;
			maxs[chan] = stat->max;
		}
		free_stats(stat);
	}
	coeff tmp[3];
	memcpy(tmp, bg, 3 * sizeof(coeff));
	/* ascending order */
	qsort(tmp, 3, sizeof(tmp[0]), cmp_coeff);
	//selecting middle channel for norm
	*norm_channel = tmp[1].channel;

	// if no selection for background we have all the stats required
	if (!area) return 0;

	// otherwise, we compute image min/max
	for (int chan = 0; chan < 3; chan++) {
		const imstats *stat = statistics(NULL, -1, fit, chan, NULL, STATS_MINMAX, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		mins[chan] = stat->min;
		maxs[chan] = stat->max;
	}
	return 0;
}

int apply_photometric_color_correction(fits *fit, const float *kw, const coeff *bg, const float *mins, const float *maxs, int norm_channel) {
	float maximum = -FLT_MAX;
	float minimum = FLT_MAX;
	float scale[3];
	float offset[3];
	float invrange;

	for (int chan = 0; chan < 3; chan++) {
		maximum = max(maximum, kw[chan] * (maxs[chan] - bg[chan].value) + bg[norm_channel].value);
		minimum = min(minimum, kw[chan] * (mins[chan] - bg[chan].value) + bg[norm_channel].value);
	}
	invrange = ((fit->type == DATA_USHORT) ? USHRT_MAX_SINGLE : 1.f) / (maximum - minimum);

	for (int chan = 0; chan < 3; chan++) {
		scale[chan] = kw[chan] * invrange;
		if (scale[chan] != scale[chan]) { // Check for NaN... If they are NaN the result image is junk
			siril_log_color_message(_("Error computing coefficients: aborting...\n"), "red");
			return 1;
		}
		offset[chan] = (-bg[chan].value * kw[chan] + bg[norm_channel].value  - minimum) * invrange;
	}
	siril_log_message("After renormalization, the following coefficients are applied\n");
	siril_log_color_message(_("White balance factors:\n"), "green");
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3f\n", chan, scale[chan]);
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
				buf[i] = roundf_to_WORD((float)buf[i] * scale[chan] + offset[chan]);
			}
		}
		else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[chan];
			for (size_t i = 0; i < n; ++i) {
				buf[i] = buf[i] * scale[chan] + offset[chan];
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
	coeff bg[3];
	float mins[3];
	float maxs[3];
	int norm_channel;

	if (!isrgb(args->fit)) {
		siril_log_message(_("Photometric color correction will do nothing for monochrome images\n"));
		return 0;
	}

	rectangle *bkg_sel = NULL;
	if (!args->bg_auto)
		bkg_sel = &(args->bg_area);

	/* we use the median of each channel to sort them by level and select
	 * the reference channel expressed in terms of order of middle median value */
	if (get_stats_coefficients(args->fit, bkg_sel, bg, mins, maxs, &norm_channel)) {
		siril_log_message(_("failed to compute statistics on image, aborting\n"));
		free(args);
		return 1;
	}
	siril_log_message(_("Normalizing on %s channel.\n"), (norm_channel == 0) ? _("red") : ((norm_channel == 1) ? _("green") : _("blue")));

	/* set photometry parameters to values adapted to the image */
	struct phot_config backup = com.pref.phot_set;
	com.pref.phot_set.force_radius = FALSE;
	com.pref.phot_set.inner = max(7.0, 3.0 * args->fwhm);
	com.pref.phot_set.outer = com.pref.phot_set.inner + 10;
	siril_log_message(_("Photometry radii set to %.1f for inner and %.1f for outer\n"),
			com.pref.phot_set.inner, com.pref.phot_set.outer);

	set_progress_bar_data(_("Photometric color calibration in progress..."), PROGRESS_RESET);
	int ret = get_white_balance_coeff(args->stars, args->nb_stars, args->fit, kw, norm_channel, args);

	if (!ret) {
		ret = apply_photometric_color_correction(args->fit, kw, bg, mins, maxs, norm_channel);
		if (args->spcc) {
			ret = spcc_colorspace_transform(args);
		}
		if (args->spcc && !ret) {
			invalidate_stats_from_fit(args->fit);
			if (!ret) {
				if (args->spcc) {
					args->fit->spcc_applied = TRUE;
				}
			}
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

		// Initialize filters if required
	if (args->spcc && !spcc_filters_initialized) {
		init_spcc_filters();
		spcc_filters_initialized = TRUE;
	}

	if (!has_wcs(args->fit)) {
		siril_log_color_message(_("Cannot run the standalone photometric color calibration on this image because it has no WCS data or it is not supported\n"), "red");
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	if (check_prior_spcc(args->fit))
		return GINT_TO_POINTER(1);

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

	pcc_star *stars = NULL;
	int nb_stars = 0;
	gboolean image_is_gfit = args->fit == &gfit;

	uint64_t sqr_radius = ((uint64_t) gfit.rx * (uint64_t) gfit.rx + (uint64_t) gfit.ry * (uint64_t) gfit.ry) / 4;
	double radius = resolution * sqrt((double)sqr_radius);	// in degrees
	double mag = args->mag_mode == LIMIT_MAG_ABSOLUTE ?
		args->magnitude_arg : compute_mag_limit_from_fov(radius * 2.0);
	if (args->mag_mode == LIMIT_MAG_AUTO_WITH_OFFSET)
		mag += args->magnitude_arg;

	int retval = 0;
	if (args->catalog == CAT_LOCAL) {
		siril_log_message(_("Getting stars from local catalogues for PCC, with a radius of %.2f degrees and limit magnitude %.2f\n"), radius * 2.0,  mag);
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
	if (args->spcc) {
		retval = siril_gaiadr3_datalink_query(siril_cat, XP_SAMPLED, &args->datalink_path);
	} else if (siril_catalog_conesearch(siril_cat) <= 0) {
		retval = 1;
	}
	/* project using WCS */
	siril_catalog_project_with_WCS(siril_cat, args->fit, TRUE, FALSE);
	stars = convert_siril_cat_to_pcc_stars(siril_cat, &nb_stars);
	retval = nb_stars == 0;

	siril_catalog_free(siril_cat);

	if (!retval) {
		if (!com.script) {
			// WARNING: Do not make this "algo" string translatable: it is used to
			// check whether SPCC has previously been applied
			const gchar *algo = args->spcc ? "SPCC" : "PCC";
			undo_save_state(args->fit, _("Photometric CC (algorithm: %s)"), algo);
		}
		args->stars = stars;
		args->nb_stars = nb_stars;
		retval = photometric_cc(args);	// args is freed from here
	}
	free(stars);
	args = NULL;

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);

	if (!retval && image_is_gfit) {
		siril_log_color_message(_("Photometric Color Calibration succeeded.\n"), "green");
		notify_gfit_modified();
	}
	else siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

// TODO: remove pcc_star?
// This interface enables for now to use new catalogues and pcc_star where required
pcc_star *convert_siril_cat_to_pcc_stars(siril_catalogue *siril_cat, int *nbstars) {
	*nbstars = 0;

	if (!siril_cat || !siril_cat->nbincluded)
		return NULL;
	if (siril_cat->projected == CAT_PROJ_NONE) {
		siril_debug_print("Catalog has not been projected\n");
	}
	if (!has_field(siril_cat, RA) || !has_field(siril_cat, DEC) || !has_field(siril_cat, MAG) || !has_field(siril_cat, BMAG))
		return NULL;
	pcc_star *results = malloc(siril_cat->nbincluded * sizeof(pcc_star));

	int n = 0;
	for (int i = 0; i < siril_cat->nbitems; i++) {
		if (n >= siril_cat->nbincluded) {
			siril_debug_print("problem when converting siril_cat to pcc_stars, more than allocated");
			break;
		}
		if (siril_cat->cat_items[i].included) {
			results[n].x = siril_cat->cat_items[i].x;
			results[n].y = siril_cat->cat_items[i].y;
			results[n].mag = siril_cat->cat_items[i].mag;
			results[n].BV = siril_cat->cat_items[i].bmag - siril_cat->cat_items[i].mag; // check for valid values was done at catalog readout
			results[n].teff = siril_cat->cat_items[i].teff; // Gaia Teff / K, computed from the sampled spectrum (better than from B-V)
			results[n].index = i; // For matching the right HDU in the Datalink query, in case of excluded stars
			n++;
		}
	}
	if (n != siril_cat->nbincluded) {
		siril_debug_print("problem when converting siril_cat to pcc_stars, number differs from catalogue info");
		free(results);
		return NULL;
	}
	*nbstars = n;
	return results;
}
