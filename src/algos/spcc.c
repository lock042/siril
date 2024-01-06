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
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
#include "algos/star_finder.h"
#include "algos/siril_wcs.h"
#include "io/single_image.h"
#include "io/image_format_fits.h" // For the datalink FITS functions
#include "io/local_catalogues.h"
#include "io/remote_catalogues.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/photometric_cc.h"
#include "registration/matching/misc.h" // for catalogue parsing helpers
#include "algos/photometric_cc.h"
#include "algos/spcc.h"
#include "algos/spcc_filters.h"


// SPCC White Points
const cmsCIExyY Whitepoint_D58 = {0.32598, 0.33532, 1.0}; // Sun as white point, modelled as D58 black body
// TODO: could model this better based on a real solar spectrum

const cmsCIExyY Whitepoint_average_galaxy = {0.345702915, 0.358538597, 1.0}; // D50
// TODO: for testing, D50 is used. Needs replacing with a value computed from real galactic spectra.

void init_spcc_filters() {
	Optolong_Blue.x = Optolong_Blue_wl;
	Optolong_Blue.y = Optolong_Blue_sr;
	Optolong_Blue.n = 72;
	Optolong_Green.x = Optolong_Green_wl;
	Optolong_Green.y = Optolong_Green_sr;
	Optolong_Green.n = 52;
	Optolong_Red.x = Optolong_Red_wl;
	Optolong_Red.y = Optolong_Red_sr;
	Optolong_Red.n = 66;
	Astrodon_RE.x = full_2nm_wl;
	Astrodon_RE.y = Astrodon_Red_E_sr;
	Astrodon_RE.n = 163;
	Astrodon_RI.x = full_2nm_wl;
	Astrodon_RI.y = Astrodon_Red_I_sr;
	Astrodon_RI.n = 163;
	Astrodon_GE.x = full_2nm_wl;
	Astrodon_GE.y = Astrodon_Green_E_sr;
	Astrodon_GE.n = 163;
	Astrodon_GI.x = full_2nm_wl;
	Astrodon_GI.y = Astrodon_Green_I_sr;
	Astrodon_GI.n = 163;
	Astrodon_B.x = full_2nm_wl;
	Astrodon_B.y = Astrodon_Blue_sr;
	Astrodon_B.n = 163;
	Chroma_Red.x = full_2nm_wl;
	Chroma_Red.y = Chroma_Red_sr;
	Chroma_Red.n = 163;
	Chroma_Green.x = full_2nm_wl;
	Chroma_Green.y = Chroma_Green_sr;
	Chroma_Green.n = 163;
	Chroma_Blue.x = full_2nm_wl;
	Chroma_Blue.y = Chroma_Blue_sr;
	Chroma_Blue.n = 163;
	Sony_IMX571M.x = Sony_IMX571_wl;
	Sony_IMX571M.y = Sony_IMX571_qe;
	Sony_IMX571M.n = 32;
	ZWO_1600M.x = ZWO_ASI1600_wl;
	ZWO_1600M.y = ZWO_ASI1600_qe;
	ZWO_1600M.n = 64;
	KAF3200.x = KAF_3200_wl;
	KAF3200.y = KAF_3200_qe;
	KAF3200.n = 27;
	KAF8300.x = KAF_8300_wl;
	KAF8300.y = KAF_8300_qe;
	KAF8300.n = 22;
	KAF1603ME.x = KAF_1603ME_wl;
	KAF1603ME.y = KAF_1603ME_qe;
	KAF1603ME.n = 25;
	Sony_ICX694.x = Sony_ICX694_wl;
	Sony_ICX694.y = Sony_ICX694_qe;
	Sony_ICX694.n = 22;
}

cmsCIExyY xpsampled_to_xyY(xpsampled* xps, const int cmf) {
	cmsCIEXYZ XYZ;
	cmsCIExyY xyY;
	double	*dbl_si_x = malloc(XPSAMPLED_LEN * sizeof(double)),
			*dbl_si_y = malloc(XPSAMPLED_LEN * sizeof(double)),
			*dbl_si_z = malloc(XPSAMPLED_LEN * sizeof(double));

	if (cmf == CMF_1931) {
		for (int i = 0 ; i < XPSAMPLED_LEN ; i++) {
			dbl_si_x[i] = xps->y[i] * x1931(xps->x[i]);
			dbl_si_y[i] = xps->y[i] * y1931(xps->x[i]);
			dbl_si_z[i] = xps->y[i] * z1931(xps->x[i]);
		}
	} else {
		for (int i = 0 ; i < XPSAMPLED_LEN ; i++) {
			dbl_si_x[i] = xps->y[i] * x1964(xps->x[i]);
			dbl_si_y[i] = xps->y[i] * y1964(xps->x[i]);
			dbl_si_z[i] = xps->y[i] * z1964(xps->x[i]);
		}
	}
	gsl_interp *interp = gsl_interp_alloc(gsl_interp_akima, XPSAMPLED_LEN);
	gsl_interp_init(interp, xps->x, dbl_si_x, XPSAMPLED_LEN);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	XYZ.X = gsl_interp_eval_integ(interp, xps->x, dbl_si_x, 380.0, 700.0, acc);
	gsl_interp_accel_reset(acc);
	free(dbl_si_x);

	gsl_interp_init(interp, xps->x, dbl_si_y, XPSAMPLED_LEN);
	XYZ.Y = gsl_interp_eval_integ(interp, xps->x, dbl_si_y, 380.0, 700.0, acc);
	gsl_interp_accel_reset(acc);
	free(dbl_si_y);

	gsl_interp_init(interp, xps->x, dbl_si_z, XPSAMPLED_LEN);
	XYZ.Z = gsl_interp_eval_integ(interp, xps->x, dbl_si_z, 380.0, 700.0, acc);
	free(dbl_si_z);
	gsl_interp_free(interp);
	gsl_interp_accel_free(acc);
	cmsXYZ2xyY(&xyY, &XYZ);
	xyY.Y = 1.f;
	return xyY;
}

void si_free(spectral_intensity *foo, gboolean free_struct) {
	free(foo->x);
	free(foo->y);
	if (free_struct)
		free(foo);
	return;
}

/* Fills a destination spectral_intensity at evenly spaced wavelength intervals
 * from a source spectral_intensity. This allows library sensor / filter
 * spectral_intensities to be stored as unevenly spaced data points that suit the
 * data, but interpolated to the same spacings as the Gaia DR3 data. */

void init_xpsampled_from_library(xpsampled *out, spectral_intensity *in) {
	const int n = in->n;
	double *dbl_x = malloc(n * sizeof(double));
	double *dbl_y = malloc(n * sizeof(double));
	for (int i = 0 ; i < n ; i++) {
		dbl_x[i] = in->x[i];
	}
	for (int i = 0 ; i < n ; i++) {
		dbl_y[i] = in->y[i];
	}
	gsl_interp *interp = gsl_interp_alloc(gsl_interp_akima, (size_t) n);
	gsl_interp_init(interp, dbl_x, dbl_y, n);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	for (int i = 0 ; i < XPSAMPLED_LEN ; i++) {
		if (out->x[i] < in->x[0] || out->x[i] > in->x[in->n-1])
			out->y[i] = 0.0;
		else
			out->y[i] = max(0.0, gsl_interp_eval(interp, dbl_x, dbl_y, out->x[i], acc));
	}
	free(dbl_x);
	free(dbl_y);
	gsl_interp_free(interp);
	gsl_interp_accel_free(acc);
	return;
}

// Takes one spectral_intensity and multiplies each value by the (interpolated)
// value of a second spectral intensity at each of the wavelength values of the
// first one. Result is returned as a new spectral_intensity*
// The two spectral_intensities must be compatible (i.e. a->n = b->n,
// a->x[i] = b->x[i] for all i.

void multiply_xpsampled(xpsampled *result, const xpsampled *a, const xpsampled *b) {
	for (int i = 0 ; i < XPSAMPLED_LEN ; i++) {
		result->y[i] = a->y[i] * b->y[i];
	}
	return;
}

// Uses the gsl interp_integ routine to evaluate the integral
// of the product of a spectral_intensity and a given CIE color
// matching function. This is returned as a cmsCIExyY with the
// Y component set to 1.0

double integrate_xpsampled(const xpsampled *xps) {
	gsl_interp *interp = gsl_interp_alloc(gsl_interp_akima, (size_t) XPSAMPLED_LEN);
	gsl_interp_init(interp, xps->x, xps->y, XPSAMPLED_LEN);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	double result = gsl_interp_eval_integ(interp, xps->x, xps->y, 380.0, 700.0, acc);
	gsl_interp_free(interp);
	gsl_interp_accel_free(acc);
	return result;
}

void get_spectrum_from_args(struct photometric_cc_data *args, xpsampled* spectrum, int chan) {
	xpsampled spectrum2 = { spectrum->x, { 0.0 } };
	mono_sensor_t selected_sensor_m = (mono_sensor_t) args->selected_sensor_m;
	rgb_sensor_t selected_sensor_rgb = (rgb_sensor_t) args->selected_sensor_rgb;
	filter_t selected_filters = (filter_t) args->selected_filters;

	if (args->spcc_mono_sensor) {
		switch (selected_sensor_m) {
			case IMX571M:
				init_xpsampled_from_library(spectrum, &Sony_IMX571M);
				break;
			case ZWO1600M:
				init_xpsampled_from_library(spectrum, &ZWO_1600M);
				break;
			case KAF_1603ME:
				init_xpsampled_from_library(spectrum, &KAF1603ME);
				break;
			case KAF_3200:
				init_xpsampled_from_library(spectrum, &KAF3200);
				break;
			case KAF_8300:
				init_xpsampled_from_library(spectrum, &KAF8300);
				break;
			case ICX_694:
				init_xpsampled_from_library(spectrum, &Sony_ICX694);
				break;
			// Add other mono sensors here
		}
		switch (selected_filters) {
			case FILTER_NONE:
				break;
			// TODO: Add filter data
			case FILTER_L_ENHANCE:
			case FILTER_DUAL:
			case FILTER_QUAD:
			case ANTLIA:
			case ASTRODON_E:
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Astrodon_RE : chan == 1 ? &Astrodon_GE : &Astrodon_B);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case ASTRODON_I:
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Astrodon_RI : chan == 1 ? &Astrodon_GI : &Astrodon_B);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case ASTRONOMIK:
			case BAADER:
			case CHROMA:
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Chroma_Red : chan == 1 ? &Chroma_Green : &Chroma_Blue);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case OPTOLONG:
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Optolong_Red : chan == 1 ? &Optolong_Green : &Optolong_Blue);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case ZWO:
				// TODO: populate with the correct data once OSC camera response curves are scanned in
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Optolong_Red : chan == 1 ? &Optolong_Green : &Optolong_Blue);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
		}
	} else {
		switch (selected_sensor_rgb) {
			case CANONT3I:
				// TODO: populate with the correct data once OSC camera response curves are scanned in
				init_xpsampled_from_library(spectrum, chan == 0 ? &Optolong_Red : chan == 1 ? &Optolong_Green : &Optolong_Blue);
				break;
			// Add other RGB sensors here
		}
		switch (selected_filters) {
			// TODO: Currently all these fall through to Optolong RGB, need to address this once more filter data is available
			case FILTER_L_ENHANCE:
			case FILTER_DUAL:
			case FILTER_QUAD:
				// TODO: populate with the correct data once I have obtained it
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Optolong_Red : chan == 1 ? &Optolong_Green : &Optolong_Blue);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case ANTLIA:
			case ASTRODON_E:
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Astrodon_RE : chan == 1 ? &Astrodon_GE : &Astrodon_B);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case ASTRODON_I:
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Astrodon_RI : chan == 1 ? &Astrodon_GI : &Astrodon_B);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case ASTRONOMIK:
			case BAADER:
			case CHROMA:
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Chroma_Red : chan == 1 ? &Chroma_Green : &Chroma_Blue);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case OPTOLONG:
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Optolong_Red : chan == 1 ? &Optolong_Green : &Optolong_Blue);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			case ZWO:
				// TODO: populate with the correct data once OSC camera response curves are scanned in
				init_xpsampled_from_library(&spectrum2, chan == 0 ? &Optolong_Red : chan == 1 ? &Optolong_Green : &Optolong_Blue);
				multiply_xpsampled(spectrum, spectrum, &spectrum2);
				break;
			default:
				// Do nothing for FILTER_NONE
				break;
		}
	}
}

// SPCC color space transform from a source profile constructed from the
// computed primaries of the sensor / filter setup and the chosen white point
// to a linear version of the working colorspace
int spcc_colorspace_transform(struct photometric_cc_data *args) {
	cmsToneCurve *curve[3], *tonecurve;
	tonecurve = cmsBuildGamma(NULL, 1.00);
	curve[0] = curve[1] = curve[2] = tonecurve;
	cmsHPROFILE source_profile = cmsCreateRGBProfile(&args->whitepoint, &args->primaries, curve);
	if (!source_profile)
		return 1;
	cmsHPROFILE dest_profile = siril_color_profile_linear_from_color_profile(com.icc.working_standard);
	if (!dest_profile) {
		cmsCloseProfile(source_profile);
		return 1;
	}

	uint32_t npixels = args->fit->rx * args->fit->ry;
	gboolean threaded = !get_thread_run();
	cmsUInt32Number fit_colorspace = cmsSigRgbData;
	cmsUInt32Number type = get_planar_formatter_type(fit_colorspace, args->fit->type, FALSE);

	cmsHTRANSFORM transform = cmsCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), source_profile, type, dest_profile, type, INTENT_ABSOLUTE_COLORIMETRIC, com.icc.rendering_flags);
	cmsCloseProfile(dest_profile);
	if (!transform) {
		cmsCloseProfile(source_profile);
		return 1;
	}
	if (args->fit->icc_profile)
		cmsCloseProfile(args->fit->icc_profile);
	args->fit->icc_profile = copyICCProfile(source_profile);
	cmsCloseProfile(source_profile);

	void *data = (args->fit->type == DATA_FLOAT) ? (void *) args->fit->fdata : (void *) args->fit->data;
	cmsUInt32Number datasize = args->fit->type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
	cmsUInt32Number bytesperline = args->fit->rx * datasize;
	cmsUInt32Number bytesperplane = npixels * datasize;
	cmsDoTransformLineStride(transform, data, data, args->fit->rx, args->fit->ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
	cmsDeleteTransform(transform);
	refresh_icc_transforms();
	color_manage(args->fit, TRUE);
	return 0;
}

int check_prior_spcc(fits *fit) {
	// Check SPCC hasn't been applied already
	GSList* entry = NULL;
	if (fit->spcc_applied)
		return 1;
	if (fit->history) {
		entry = fit->history;
		while (entry) {
			if (strstr(entry->data, "SPCC")) {
				gchar *msg = g_strdup("Spectrophotometric Color Correction "
							"has already been applied to this image. Re-applying it will "
							"result in incorrect results!");
				if (!com.script) {
					if (!siril_confirm_dialog(_("Warning!"), _(msg), _("Continue"))) {
						g_free(msg);
						return 1;
					}
					break;
				} else {
					siril_log_color_message(_("Warning! %s\n"), "red", msg);
					g_free(msg);
				}
			}
			entry = entry->next;
		}
	}
	return 0;
}
