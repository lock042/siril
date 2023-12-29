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

#include <gsl/gsl_statistics.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "algos/statistics_float.h"
#include "algos/photometry.h"
#include "algos/spcc_filters.h"
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
#include "algos/star_finder.h"
#include "algos/siril_wcs.h"
#include "io/single_image.h"
#include "io/image_format_fits.h" // For the datalink FITS functions
#include "io/local_catalogues.h"
#include "io/remote_catalogues.h"
#include "gui/progress_and_log.h"
#include "registration/matching/misc.h" // for catalogue parsing helpers
#include "photometric_cc.h"

//#define DEBUG_PCC

static gboolean spcc_filters_initialized = FALSE;

static const cmsCIEXYZ D65 = {0.95045471, 1.0, 1.08905029};
static const cmsCIEXYZ D50 = {0.964199999, 1.000000000, 0.824899998};

enum {
	RED, GREEN, BLUE
};

static void init_spcc_filters() {
	Optolong_Blue.wl = Optolong_Blue_wl;
	Optolong_Blue.si = Optolong_Blue_sr;
	Optolong_Blue.n = 71;
	Optolong_Green.wl = Optolong_Green_wl;
	Optolong_Green.si = Optolong_Green_sr;
	Optolong_Green.n = 50;
	Optolong_Red.wl = Optolong_Red_wl;
	Optolong_Red.si = Optolong_Red_sr;
	Optolong_Red.n = 65;
	Sony_IMX571M.wl = Sony_IMX571_wl;
	Sony_IMX571M.si = Sony_IMX571_qe;
	Sony_IMX571M.n = 29;
}

// CIE XYZ Color Matching Functions
// Ref: Wylie / Sloan / Shirley, Simple Analytic Approximations to the CIE XYZ
// Color Matching Functions, Journal of Computer Graphics Techniques, Vol. 2,
// No. 2, 2013.

// CIE 1931 2-degree CMF using multi-lobe piecewise Gaussian
float x1931(float w) {
	float t1 = (w-442.0f)*((w<442.0f)?0.0624f:0.0374f);
	float t2 = (w-599.8f)*((w<599.8f)?0.0264f:0.0323f);
	float t3 = (w-501.1f)*((w<501.1f)?0.0490f:0.0382f);
	return 0.362f*expf(-0.5f*t1*t1) + 1.056f*expf(-0.5f*t2*t2)
	- 0.065f*expf(-0.5f*t3*t3);
}

float y1931( float w) {
	float t1 = (w-568.8f)*((w<568.8f)?0.0213f:0.0247f);
	float t2 = (w-530.9f)*((w<530.9f)?0.0613f:0.0322f);
	return 0.821f*expf(-0.5f*t1*t1) + 0.286f*expf(-0.5f*t2*t2);
}

float z1931(float w) {
	float t1 = (w-437.0f)*((w<437.0f)?0.0845f:0.0278f);
	float t2 = (w-459.0f)*((w<459.0f)?0.0385f:0.0725f);
	return 1.217f*expf(-0.5f*t1*t1) + 0.681f*expf(-0.5f*t2*t2);
}

// CIE 1964 10-degree CMF using single-lobe fit
float x1964(float w) {
	float i1 = 0.4f*expf(-1250.f*powf(logf((w+570.f)/1014.0f),2.f));
	float i2 = 1.13f*expf(-234.0f*powf(logf((1338.0f-w)/743.5f),2.f));
	return i1+i2;
}

float y1964(float w) {
	return 1.011f*expf(-0.5f*pow((w-556.1f)/46.14f, 2.f));
}

float z1964(float w) {
	return 2.06f*expf( -32.0f*powf(logf((w-266.0f)/180.4f),2.f));
}

void si_free(spectral_intensity *foo) {
	free(foo->wl);
	free(foo->si);
	free(foo);
	return;
}

/* Returns the spectral response of a filter (or QE of a sensor) at wavelength wl.
 * wavelengths and filter are float* arrays defined in spcc_filters.h
 * They MUST have the same length and all arrays must be terminated by FLT_MAX */

float sr(const float wl, const float *wavelengths, const float *filter) {
	if (wl < wavelengths[0]) // Provided wl is lower than the min wavelengths range
		return 0.f;
	int i = 0;
	while (wl > wavelengths[i+1]) // Advance to the required values
		i++;
	if (wavelengths[i+1] == FLT_MAX) // Provided wl is greater than the max wavelengths range
		return 0.f;

	float wl1 = wavelengths[i], wl2 = wavelengths[i+1], t1, t2, sr;
	t1 = filter[i];
	t2 = filter[i+1];
	sr = t1 + ((wl - wl1) / (wl2 - wl1)) * (t2 - t1);
	return (sr);
}

// Returns the emittance of a Planckian black body spectrum at wavelength
// wl and temperature bbTemp
// from "Colour Rendering of Spectra", John Walker, Fourmilab. Public domain code, last updated March 9 2003
float bb_spectrum(float wl, float bbTemp)
{
    float wlm = wl * 1e-9;   /* Wavelength in meters */

    return (3.74183e-16f * pow(wlm, -5.f)) /
           (expf(1.4388e-2f / (wlm * bbTemp)) - 1.f);
}

// Convert a spectrum to CIE xy
/* spec_intens is a float buffer containing spectral intensities, terminated by
 * -FLT_MAX
 * wavelength is a float buffer  containing wavelengths such that spec_intens[i]
 * is the intensity at wavelength[i]
 * wavelengths should be evenly spaced.
 */
enum {
	CMF_1931,
	CMF_1964
};

void spectrum_to_XYZ(const pointf *data, uint32_t n,
					 cmsCIEXYZ *xyz, const int cmf)
{
    double X = 0., Y = 0., Z = 0.;//, XYZ;

	if (cmf == CMF_1931) {
		for (uint32_t i = 0 ; i < n ; i++) {
			X += data[i].y * x1931(data[i].x);
			Y += data[i].y * y1931(data[i].x);
			Z += data[i].y * z1931(data[i].x);
		}
	} else {
		for (uint32_t i = 0 ; i < n ; i++) {
			X += data[i].y * x1964(data[i].x);
			Y += data[i].y * y1964(data[i].x);
			Z += data[i].y * z1964(data[i].x);
		}
	}
//    XYZ = (X + Y + Z);
    xyz->X = X;
	xyz->Y = Y;
    xyz->Z = Z;
}

// Covers 400nm to 700nm at 1nm resolution
// Takes a NULL-terminated list of spectral_intensity structs
cmsCIExyY pipeline_to_xyY(spectral_pipeline* pipeline, const int cmf) {
	cmsCIEXYZ XYZ;
	cmsCIExyY xyY;
	pointf* data = malloc(301 * sizeof(pointf));
	for (int i = 0 ; i < 301 ; i++) {
		float wl = 400.f + (float) i;
		pointf p = {wl, 1.f};
		for (int filter = 0 ; filter < pipeline->n ; filter++) {
			p.y *= sr(wl, pipeline->si[filter].wl, pipeline->si[filter].si);
			filter++;
		}
		memcpy(&data[i], &p, sizeof(pointf));
	}
	spectrum_to_XYZ(data, 301, &XYZ, cmf);
	cmsXYZ2xyY(&xyY, &XYZ);
	xyY.Y = 1.f;
	return xyY;
}

// Reference: https://en.wikipedia.org/wiki/Color_index and https://arxiv.org/abs/1201.1809 (Ballesteros, F. J., 2012)
// Uses Ballesteros' formula based on considering stars as black bodies
cmsFloat64Number bvToT(float bv) {
	cmsFloat64Number t = 4600. * ((1. / ((0.92 * bv) + 1.7)) + (1. / ((0.92 * bv) + 0.62)));
	return t;
}

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
static void TempK2rgb(float *r, float *g, float *b, float TempK, cmsHTRANSFORM transform) { //, cmsCIEXYZ *profile_whitepoint) { // RGB <0,1> <- BV <-0.4,+2.0> [-]
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
#ifdef DEBUG_PCC
	fprintf(stderr,"%f %f %f %f %f %f %f\n", bv, TempK, WhitePoint.x, WhitePoint.y, *r, *g, *b);
#endif
}

static int make_selection_around_a_star(double fx, double fy, rectangle *area, fits *fit) {
	/* make a selection around the star, coordinates are in display reference frame */
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

static int get_white_balance_coeff(pcc_star *stars, int nb_stars, fits *fit, float *kw, int norm_channel) {
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

	cmsHPROFILE xyzprofile = cmsCreateXYZProfile();
	cmsHPROFILE profile = NULL;
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
	cmsHTRANSFORM transform = NULL;
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
		siril_log_color_message(_("Error: failed to set up colorspace transform. This is a bug, please report it!\n"), "red");
		return 1;
	}
#ifdef DEBUG_PCC
	fprintf(stderr,"BV tempK x y r g b\n");
#endif
	int gaia_xpsamp = 0, gaia_teff = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) shared(progress, ngood)
#endif
	for (int i = 0; i < nb_stars; i++) {
		if (!get_thread_run())
			continue;
		if (stars[i].teff)
			gaia_teff++;
		rectangle area = { 0 };
		float flux[3] = { 0.f, 0.f, 0.f };
		float r, g, b, bv;
		if (!(g_atomic_int_get(&progress) % 16))	// every 16 iterations
			set_progress_bar_data(NULL, (double) progress / (double) nb_stars);
		g_atomic_int_inc(&progress);

		if (make_selection_around_a_star(stars[i].x, stars[i].y, &area, fit)) {
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
		// get r g b coefficient
		// If the Gaia Teff field is populated, use that as it is more accurate.
		// Otherwise, we convert from Johnson B-V
		if (stars[i].teff == 0.f) {
			bv = min(max(stars[i].BV, -0.4f), 2.f);
			cmsFloat64Number TempK = bvToT(bv);
			stars[i].teff = (float) TempK;
		}
		// TempK2rgb converts the temperature to RGB values in the working colorspace
		TempK2rgb(&r, &g, &b, stars[i].teff, transform);
		/* get Color calibration factors for current star */
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
	printf("nb good stars %d, nb with Gaia TempK %d, nb with Gaia sampled spectrum %d\n", ngood, gaia_teff, gaia_xpsamp);
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
	kw[RED] /= (kw[norm_channel]);
	kw[GREEN] /= (kw[norm_channel]);
	kw[BLUE] /= (kw[norm_channel]);
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
static int get_stats_coefficients(fits *fit, rectangle *area, coeff *bg, float *mins, float *maxs, int *norm_channel) {
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

static int apply_photometric_color_correction(fits *fit, const float *kw, const coeff *bg, const float *mins, const float *maxs, int norm_channel) {
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
	int ret = get_white_balance_coeff(args->stars, args->nb_stars, args->fit, kw, norm_channel);

	if (!ret) {
		ret = apply_photometric_color_correction(args->fit, kw, bg, mins, maxs, norm_channel);
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
			case CAT_GAIADR3_DIRECT:
				mag = min(mag, 18.0);
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
	siril_cat->phot = TRUE;

	/* Fetching the catalog*/
	if (siril_catalog_conesearch(siril_cat) <= 0) {
		retval = 1;
	} else {
		/* project using WCS */
		siril_catalog_project_with_WCS(siril_cat, args->fit, TRUE, FALSE);
		stars = convert_siril_cat_to_pcc_stars(siril_cat, &nb_stars);
		retval = nb_stars == 0;
	}
	siril_catalog_free(siril_cat);

	if (!retval) {
		if (!com.script) {
			undo_save_state(args->fit, _("Photometric CC"));
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

/* run the SPCC using the existing star list of the image from the provided file */
int spectrophotometric_cc(struct photometric_cc_data *args) {
	float kw[3];
	coeff bg[3];
	float mins[3];
	float maxs[3];
	int norm_channel;

	// Initialize filters if required
	if (!spcc_filters_initialized) {
		init_spcc_filters();
		spcc_filters_initialized = TRUE;
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

	set_progress_bar_data(_("Spectrophotometric color calibration in progress..."), PROGRESS_RESET);

	//	Calculate filter responses at 2nm spacing from 380nm to 700nm
	//	Either look up the responses (for OSC) or calculate as a product of filter_resp * QE
	spectral_pipeline *pipeline[3] = { NULL };
	for (int chan = 0 ; chan < 3 ; chan++) {
		pipeline[chan] = get_pipeline_from_ui(chan, 380, 700);
	}

	// Calculate effective primaries (and get the desired white point while we're
	// at it) for later
	cmsCIExyYTRIPLE primaries;
	primaries.Red = pipeline_to_xyY(pipeline[RED], CMF_1931);
	primaries.Green = pipeline_to_xyY(pipeline[GREEN], CMF_1931);
	primaries.Blue = pipeline_to_xyY(pipeline[BLUE], CMF_1931);
	cmsCIExyY spcc_whitepoint = get_spcc_whitepoint_from_UI();

	//	Do a Gaia DR3 cone search and retrieve the reference stars CSV and the
	//	xp_sampled FITS
	gchar *datalink_path = NULL;
	siril_gaiadr3_datalink_query(args->ref_stars, XP_SAMPLED, &datalink_path);

	struct phot_config *ps = phot_set_adjusted_for_image(args->fit);
	siril_debug_print("aperture: %2.1f%s\tinner: %2.1f\touter: %2.1f\n", ps->aperture, ps->force_radius?"":" (auto)", ps->inner, ps->outer);
	gint ngood = 0, progress = 0;
	gint errors[PSF_ERR_MAX_VALUE] = { 0 };

	//	Obtain spectrophotometric white balance
	for (int i = 0 ; i < args->ref_stars->nbitems ; i++) {
		// Obtain sampled spectrum between 380-700nm from the xp_sampled FITS
		// Each call to get_xpsampled() opens a new fptr so it is safe to call
		// from multiple threads
		spectral_intensity *xp_sampled = get_xpsampled(datalink_path, i, 380.f, 700.f);
		float ref_flux[3] = { 0.f };
		for (int chan = 0 ; chan < 3 ; chan++) {
			pipeline[chan]->si = realloc(pipeline[chan]->si, pipeline[chan]->n+1);
			pipeline[chan]->si[pipeline[chan]->n] = *xp_sampled;
			pipeline[chan]->n++;
			// Multiply by the filter responses
			ref_flux[chan] = integrate_pipeline(pipeline[chan]);
			si_free(&pipeline[chan]->si[pipeline[chan]->n--]);
		}
		// Obtain the expected ratios of R, G and B
		float max_flux = max(max(ref_flux[0], ref_flux[1]), ref_flux[2]);
		for (int chan = 0 ; chan < 3 ; chan++)
			ref_flux[chan] /= max_flux;

		// Measure actual flux by running peaker around ref star coordinates
		rectangle area = { 0 };
		float flux[3] = { 0.f };
		if (make_selection_around_a_star(args->ref_stars->cat_items[i].x, args->ref_stars->cat_items[i].y, &area, args->fit)) {
			siril_debug_print("star %d is outside image or too close to border\n", i);
			g_atomic_int_inc(errors+PSF_ERR_OUT_OF_WINDOW);
			continue;
		}
		gboolean no_phot = FALSE;
		psf_error error = PSF_NO_ERR;
		for (int chan = 0; chan < 3 && !no_phot; chan ++) {
			psf_star *photometry = psf_get_minimisation(args->fit, chan, &area, TRUE, ps, FALSE, com.pref.starfinder_conf.profile, &error);
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

		// Compare ref_flux with actual flux
		// (rest is as per PCC)
	}

	// Apply white balance as per PCC
	int ret = 0;
	if (!ret) {
		ret = apply_photometric_color_correction(args->fit, kw, bg, mins, maxs, norm_channel);
	} else {
		set_progress_bar_data(_("Spectrophotometric Color Calibration failed"), PROGRESS_DONE);
	}

	// Transform image from source profile to working profile using
	// INTENT_ABSOLUTE_COLORIMETRIC
	cmsToneCurve *curve[3], *tonecurve;
	tonecurve = cmsBuildGamma (NULL, 1.00);
	curve[0] = curve[1] = curve[2] = tonecurve;
	cmsHPROFILE source_profile = cmsCreateRGBProfile(&spcc_whitepoint, &primaries, curve);
	cmsFreeToneCurve(tonecurve);

	cmsHPROFILE profile = NULL;
	if (args->fit->icc_profile) {
		profile = copyICCProfile(args->fit->icc_profile);
		if (!fit_icc_is_linear(args->fit)) {
			siril_log_color_message(_("Image color space is nonlinear. It is recommended to "
					"apply photometric color calibration to linear images.\n"), "salmon");
		}
	} else {
		profile = siril_color_profile_linear_from_color_profile(com.icc.working_standard);
		args->fit->icc_profile = copyICCProfile(profile);
		color_manage(args->fit, (args->fit->icc_profile != NULL));
	}

	// Create transform from source profile to image-based profile, clean up the profiles
	// TODO...

	// Apply transform and clean up
	// TODO...

	com.pref.phot_set = backup;
	free(args);
	return ret;
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
			results[n].gaiasourceid = siril_cat->cat_items[i].gaiasourceid; // For building a Gaia DR3 Datalink query
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
