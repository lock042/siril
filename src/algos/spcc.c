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

/* Fills a destination xpsampled at evenly spaced wavelength intervals
 * from a source spcc_object. This allows library sensor / filter
 * spectral_intensities to be stored as unevenly spaced data points that suit the
 * data, but interpolated to the same spacings as the Gaia DR3 data. */

void init_xpsampled_from_library(xpsampled *out, spcc_object *in) {
	const int n = in->n;
	gsl_interp *interp = gsl_interp_alloc(gsl_interp_akima, (size_t) n);
	gsl_interp_init(interp, in->x, in->y, n);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	for (int i = 0 ; i < XPSAMPLED_LEN ; i++) {
		if (out->x[i] < in->x[0] || out->x[i] > in->x[n-1])
			out->y[i] = 0.0;
		else
			out->y[i] = max(0.0, gsl_interp_eval(interp, in->x, in->y, out->x[i], acc));
	}
	gsl_interp_free(interp);
	gsl_interp_accel_free(acc);
	return;
}

// Takes one xpsampled and multiplies each value by the (interpolated)
// value of a second spectral intensity at each of spcc_objectthe wavelength values of the
// first one. Result is passed back in *result.
// The two spectral_intensities must be compatible (i.e. a->n = b->n,
// a->x[i] = b->x[i] for all i.

void multiply_xpsampled(xpsampled *result, const xpsampled *a, const xpsampled *b) {
	for (int i = 0 ; i < XPSAMPLED_LEN ; i++) {
		result->y[i] = a->y[i] * b->y[i];
	}
	return;
}

// Uses the gsl interp_integ routine to evaluate the integral
// of the product of an xpsampled and a given CIE color
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

	if (args->spcc_mono_sensor) {
		GList *sensor = g_list_nth(com.spcc_data.mono_sensors, args->selected_sensor_m);
		load_spcc_object_arrays( (spcc_object*) sensor->data);
		init_xpsampled_from_library(spectrum, (spcc_object*) sensor->data);
		int selected_filter = chan == 0 ? args->selected_filter_r : chan == 1 ? args->selected_filter_g : args->selected_filter_b;
		GList *filter = g_list_nth(com.spcc_data.mono_filters[chan], selected_filter);
		load_spcc_object_arrays( (spcc_object*) filter->data);
		init_xpsampled_from_library(&spectrum2, (spcc_object*) filter->data);
		spcc_object_free_arrays( (spcc_object*) filter->data);
		multiply_xpsampled(spectrum, spectrum, &spectrum2);
	} else {
		// The 3 channels of an OSC sensor are stored in RGB order in the JSON file and will be in order in the GList.
		GList *object = g_list_nth(com.spcc_data.osc_sensors, args->selected_sensor_osc);
		osc_sensor *osc = (osc_sensor*) object->data;
		spcc_object *sensor = &osc->channel[chan];
		load_spcc_object_arrays( (spcc_object*) sensor);
		init_xpsampled_from_library(spectrum, (spcc_object*) sensor);
		spcc_object_free_arrays( (spcc_object*) sensor);
		if (args->use_osc_filter) {
			GList *filter = g_list_nth(com.spcc_data.osc_sensors, args->selected_filter_osc);
			load_spcc_object_arrays( (spcc_object*) filter->data);
			init_xpsampled_from_library(&spectrum2, (spcc_object*) filter->data);
			spcc_object_free_arrays( (spcc_object*) filter->data);
			multiply_xpsampled(spectrum, spectrum, &spectrum2);
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
	cmsDoTransformLineStride(transform, data, data, args->fit->rx, args->fit->ry,
							 bytesperline, bytesperline, bytesperplane, bytesperplane);
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
