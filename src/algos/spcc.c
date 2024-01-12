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

const double xpsampled_wl[XPSAMPLED_LEN] = {336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,402,404,406,408,410,412,414,416,418,420,422,424,426,428,430,432,434,436,438,440,442,444,446,448,450,452,454,456,458,460,462,464,466,468,470,472,474,476,478,480,482,484,486,488,490,492,494,496,498,500,502,504,506,508,510,512,514,516,518,520,522,524,526,528,530,532,534,536,538,540,542,544,546,548,550,552,554,556,558,560,562,564,566,568,570,572,574,576,578,580,582,584,586,588,590,592,594,596,598,600,602,604,606,608,610,612,614,616,618,620,622,624,626,628,630,632,634,636,638,640,642,644,646,648,650,652,654,656,658,660,662,664,666,668,670,672,674,676,678,680,682,684,686,688,690,692,694,696,698,700,702,704,706,708,710,712,714,716,718,720,722,724,726,728,730,732,734,736,738,740,742,744,746,748,750,752,754,756,758,760,762,764,766,768,770,772,774,776,778,780,782,784,786,788,790,792,794,796,798,800,802,804,806,808,810,812,814,816,818,820,822,824,826,828,830,832,834,836,838,840,842,844,846,848,850,852,854,856,858,860,862,864,866,868,870,872,874,876,878,880,882,884,886,888,890,892,894,896,898,900,902,904,906,908,910,912,914,916,918,920,922,924,926,928,930,932,934,936,938,940,942,944,946,948,950,952,954,956,958,960,962,964,966,968,970,972,974,976,978,980,982,984,986,988,990,992,994,996,998,1000,1002,1004,1006,1008,1010,1012,1014,1016,1018,1020};

xpsampled init_xpsampled() {
	xpsampled xps = { xpsampled_wl, { 0.0 } };
	return xps;
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
//	xyY.Y = 1.f;
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
		// TODO: Also add a LPF filter if the OSC sensor is a DSLR
		GList *filter = g_list_nth(com.spcc_data.osc_sensors, args->selected_filter_osc);
		load_spcc_object_arrays( (spcc_object*) filter->data);
		init_xpsampled_from_library(&spectrum2, (spcc_object*) filter->data);
		spcc_object_free_arrays( (spcc_object*) filter->data);
		multiply_xpsampled(spectrum, spectrum, &spectrum2);
	}
}

// Set a source profile constructed from the computed primaries of the sensor /
// filter setup. A D50 white point is nominally assigned but this would only be
// used with INTENT_ABSOLUTE_COLORIMETRIC. Since white point correction is
// already carried out, that intent should generally not be used when converting
// from this profile.
int spcc_set_source_profile(struct photometric_cc_data *args) {
	cmsCIExyY d50_illuminant_specs = {0.345702915, 0.358538597, 1.0};
	cmsCIEXYZ d50_illuminant_specs_media_whitepoint = {0.964199999, 1.000000000, 0.824899998};
	cmsMLU *copyright = cmsMLUalloc(NULL, 1);
	cmsMLUsetASCII(copyright, "en", "US", "Copyright 2024, Team free-astro (https://free-astro.org/index.php/Siril), CC-BY-SA 3.0 Unported (https://creativecommons.org/licenses/by-sa/3.0/");
	cmsToneCurve *curve[3], *tonecurve;
	tonecurve = cmsBuildGamma(NULL, 1.00);
	curve[0] = curve[1] = curve[2] = tonecurve;
	cmsHPROFILE profile = cmsCreateRGBProfile(&d50_illuminant_specs, &args->primaries, curve);
	if (!profile)
		return 1;
	cmsWriteTag(profile, cmsSigCopyrightTag, copyright);
	cmsWriteTag (profile, cmsSigMediaWhitePointTag, &d50_illuminant_specs_media_whitepoint);
	cmsMLU *description = cmsMLUalloc(NULL, 1);
	GList *sensorlistitem, *filterlistitem;
	if (args->spcc_mono_sensor) {
		sensorlistitem = g_list_nth(com.spcc_data.mono_sensors, args->selected_sensor_m);
		filterlistitem = g_list_nth(com.spcc_data.mono_filters[0], args->selected_filter_r);
	} else {
		sensorlistitem = g_list_nth(com.spcc_data.osc_sensors, args->selected_sensor_osc);
		filterlistitem = g_list_nth(com.spcc_data.osc_sensors, args->selected_filter_osc);
	}
	spcc_object *sensor = (spcc_object*) sensorlistitem->data;
	spcc_object *filter = (spcc_object*) filterlistitem->data;
	gchar *description_text = g_strdup_printf("Siril SPCC sensor source profile (linear). Sensor: %s, filters: %s", sensor->name, filter->model);
	cmsMLUsetASCII(description, "en", "US", description_text);
	cmsWriteTag(profile, cmsSigProfileDescriptionTag, description);
	cmsMLUfree(description);
	cmsMLUfree(copyright);
	g_free(description_text);
	if (args->fit->icc_profile) {
		cmsCloseProfile(args->fit->icc_profile);
		args->fit->icc_profile = NULL;
	}
	// As the existing profile is NULL, we are just assigning here.
	siril_colorspace_transform(args->fit, profile);
	return 0;
}
