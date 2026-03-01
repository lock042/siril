/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#include "filters.h"
#include "core/siril.h"
#include "core/icc_profile.h"
#include <math.h>

/* A common narrow-band filter list. */
narrow_filter narrow_band_filters[] = {
	/* narrow band filters */
	{ "H-alpha", 656.28 },
	{ "H-beta", 486.1 },
	{ "O III", 500.7 },
	{ "S II", 671.7 },
	{ "N II", 658.35 },
	{ "Ca II", 393.37 }
};

int get_nb_narrow_filters() {
	return G_N_ELEMENTS(narrow_band_filters);
}

/*
 * TODO: Add common broad band filters with their number
 * example list: http://www.myastroshop.com.au/guides/filters.asp
 */
broad_filter broad_band_filters[] = {
	/* broad band filters */
	{ "#1 (red)", "#ff0000" },
	{ "#2 (blue)", "#0000ff" },
};

// TODO: this function should be replaced by the CIE cmfs from algos/photometric_cc.c for consistency
// Still need to decide on CIE 1931 2 degree or 1964 10 degree
void wavelength_to_XYZ(float wavelength, float *X, float *Y, float *Z) {
    float Xt1 = (wavelength - 442.0f) *
                ((wavelength < 442.0f) ? 0.0624f : 0.0374f);
    float Xt2 = (wavelength - 599.8f) *
                ((wavelength < 599.8f) ? 0.0264f : 0.0323f);
    float Xt3 = (wavelength - 501.1f) *
                ((wavelength < 501.1f) ? 0.0490f : 0.0382f);
    *X =   0.362f * expf(-0.5f*Xt1*Xt1) +
                1.056f * expf(-0.5f*Xt2*Xt2) -
                0.065f * expf(-0.5f*Xt3*Xt3);

    float Yt1 = (wavelength - 568.8f) *
                ((wavelength < 568.8f) ? 0.0213f : 0.0247f);
    float Yt2 = (wavelength - 530.9f) *
                ((wavelength < 530.9f) ? 0.0613f : 0.0322f);
    *Y = 0.821f * expf(-0.5f*Yt1*Yt1) +
                0.286f * expf(-0.5f*Yt2*Yt2);

    float Zt1 = (wavelength - 437.0f) *
                ((wavelength < 437.0f) ? 0.0845f : 0.0278f);
    float Zt2 = (wavelength - 459.0f) *
                ((wavelength - 459.0f) ? 0.0385f : 0.0725f);
    *Z =   1.217f * expf(-0.5f*Zt1*Zt1) +
                0.681f * expf(-0.5f*Zt2*Zt2);
}

void wavelength_to_display_RGB(double wavelength, GdkRGBA *rgb) {
	float XYZ[3] = { 0.f };
	float RGB[3] = { 0.f };
	wavelength_to_XYZ((float) wavelength, &XYZ[0], &XYZ[1], &XYZ[2]);
	if (gui.icc.monitor) {
		cmsHPROFILE profile_xyz = cmsCreateXYZProfile();
		// This transform is unbounded: Gdk crops the negative values returned
		// for the displayed colors, but retaining the unbounded values allows
		// them to be used to convert to the image color space later.
		cmsHTRANSFORM transform = cmsCreateTransformTHR(com.icc.context_single,
														profile_xyz,
														TYPE_XYZ_FLT_PLANAR,
														gui.icc.monitor,
														TYPE_RGB_FLT_PLANAR,
														INTENT_RELATIVE_COLORIMETRIC,
														0);
		cmsDoTransform(transform, &XYZ, &RGB, 1);
		cmsDeleteTransform(transform);
		cmsCloseProfile(profile_xyz);
	} else {
		siril_debug_print("Error: no monitor profile exists. This is a bug!\n");
	}
	rgb->red = (double) RGB[0];
	rgb->green = (double) RGB[1];
	rgb->blue = (double) RGB[2];
	rgb->alpha = 1.0;
}
