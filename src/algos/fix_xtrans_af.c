/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"

#include "fix_xtrans_af.h"

supported_xtrans_list supported_xtrans[] =
		{
		// Camera Name      AF Pixels x,y,w,h        Sample x,y,w,h
		{ "Fujifilm X-T2", { 1510, 505, 3010, 3017 }, { 1992, 990, 2048, 2048 } },
		{ "Fujifilm X-T20", { 1510, 505, 3010, 3017 }, { 1992, 990, 2048, 2048 } },
		{ "Fujifilm X-Pro2", { 1500, 505, 3010, 3017 }, { 1992, 990, 2048, 2048 } }
};

static int get_nb_xtrans_supported() {
	return G_N_ELEMENTS(supported_xtrans);
}

static int get_model(const char *model) {
	int nb_xtrans = get_nb_xtrans_supported();
	for (int i = 0; i < nb_xtrans; i++) {
		if (!g_ascii_strcasecmp(model, supported_xtrans[i].model))
			return i;
	}
	return -1;
}

static void set_af_matrix(gchar *pattern, af_pixel_matrix af_matrix) {
	// af_pixel_matrix is [12][6].
        // Lowercase are AF pixels.  Uppercase are regular.

	if (!g_ascii_strcasecmp("GGRGGBGGBGGRBRGRBGGGBGGRGGRGGBRBGBRG", pattern)) {
		memcpy(af_matrix[0], "GGRGGB", 6);
		memcpy(af_matrix[1], "GGBGGR", 6);
		memcpy(af_matrix[2], "BRGRBG", 6);
		memcpy(af_matrix[3], "GgBGgR", 6);
		memcpy(af_matrix[4], "GGRGGB", 6);
		memcpy(af_matrix[5], "RBGBRG", 6);
		memcpy(af_matrix[6], "GGRGGB", 6);
		memcpy(af_matrix[7], "GgBGgR", 6);
		memcpy(af_matrix[8], "BRGRBG", 6);
		memcpy(af_matrix[9], "GGBGGR", 6);
		memcpy(af_matrix[10], "GGRGGB", 6);
		memcpy(af_matrix[11], "RBGBRG", 6);
	} else if (!g_ascii_strcasecmp("RBGBRGGGRGGBGGBGGRBRGRBGGGBGGRGGRGGB", pattern)) {
		memcpy(af_matrix[0], "RBGBRG", 6);
		memcpy(af_matrix[1], "GGRGGB", 6);
		memcpy(af_matrix[2], "GGBGGR", 6);
		memcpy(af_matrix[3], "BRGRBG", 6);
		memcpy(af_matrix[4], "GgBGgR", 6);
		memcpy(af_matrix[5], "GGRGGB", 6);
		memcpy(af_matrix[6], "RBGBRG", 6);
		memcpy(af_matrix[7], "GGRGGB", 6);
		memcpy(af_matrix[8], "GgBGgR", 6);
		memcpy(af_matrix[9], "BRGRBG", 6);
		memcpy(af_matrix[10], "GGBGGR", 6);
		memcpy(af_matrix[11], "GGRGGB", 6);
	}
}

// This returns the pixel type based on our AF matrix if we are within the AF rectangle.
// It returns an X if we are outside of the AF rectangle.
char get_pixel_type(rectangle af, int x, int y, af_pixel_matrix *af_matrix) {

	if (x >= af.x && x <= (af.x + af.w) && y >= af.y && y <= (af.y + af.h)) {
		// We are within the AF rectangle.
		// This is written assuming we don't know the size of the matrix.
		int matrix_cols = sizeof((*af_matrix)[0]);
		int matrix_rows = sizeof((*af_matrix)) / sizeof((*af_matrix)[0]);

		// This will return the corresponding pixel type.
	//	printf("matrix_cols=%d\tmatrix_rows=%d\n", matrix_cols, matrix_rows);
		return (*af_matrix)[y % matrix_rows][x % matrix_cols];
	} else {
		// We are outside of the AF rectangle.
		return 'X';
	}
}

static int subtract_fudge(fits *fit, rectangle af, float fudge, af_pixel_matrix *af_matrix ) {
	int width = fit->rx;
	int height = fit->ry;

	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->pdata[RLAYER];

		for (unsigned int y = 0; y < height; y++) {
			for (unsigned int x = 0; x < width; x++) {
				if (get_pixel_type(af, x, y, af_matrix) == 'g') {
					// This is an auto focus pixel.  Subtract the fudge.
					buf[x + y * width] -= roundf_to_WORD(fudge);
				}
			}
		}
	} else if (fit->type == DATA_FLOAT) {
		float *buf = fit->fpdata[RLAYER];

		for (unsigned int y = 0; y < height; y++) {
			for (unsigned int x = 0; x < width; x++) {
				if (get_pixel_type(af, x, y, af_matrix) == 'g') {
					// This is an auto focus pixel.  Subtract the fudge.
					buf[x + y * width] -= fudge;
				}
			}
		}
	}
	// the caller should call invalidate_stats_from_fit(fit);
	return 0;
}

int fix_xtrans_ac(fits *fit) {
	rectangle af, sam;

	int model = get_model(fit->instrume);
	if (model < 0) {
		siril_log_color_message(_("Fix X-Trans: Unknown camera %s, trying to read information from preferences.\n"), "red", fit->instrume);
		if (com.pref.xtrans_af.w != 0 && com.pref.xtrans_af.h != 0) {
			if (com.pref.xtrans_sample.w > fit->rx || com.pref.xtrans_sample.h > fit->ry) {
				siril_log_color_message(_("Sample box cannot be bigger than the image.\n"), "red");
				return 1;
			}
			if (com.pref.xtrans_af.w > fit->rx || com.pref.xtrans_af.h > fit->ry) {
				siril_log_color_message(_("AF box cannot be bigger than the image.\n"), "red");
				return 1;
			}
			af = com.pref.xtrans_af;
			if (com.pref.xtrans_sample.w != 0 && com.pref.xtrans_sample.h != 0) {
				sam = com.pref.xtrans_sample;
			} else {
				sam.x = 0;
				sam.y = 0;
				sam.w = fit->rx;
				sam.h = fit->ry;
			}
		} else {
			siril_log_color_message(_("No information available in preferences.\n"), "red");
			return 1;
		}
	} else {
		af = supported_xtrans[model].af;
		sam = supported_xtrans[model].sample;
	}


	// non-focus pixels
	double nfsum = 0.0;
	float nfmean = 0.f;
	long nfcount = 0L;

	// auto focus pixels
	double afsum = 0.0;
	float afmean = 0.f;
	long afcount = 0L;

	// The fudge amount to apply to auto focus pixels. (computed)
	float fudge;

	// af_matrix is an RGB pattern where lowercase letters represent AF pixels and their color.
	af_pixel_matrix af_matrix = { 0 };
	set_af_matrix(fit->bayer_pattern, af_matrix);

	WORD *buf = fit->pdata[RLAYER];
	float *fbuf = fit->fpdata[RLAYER];

	if (af_matrix[0][0] == 0) {
		siril_log_color_message(_("This CFA pattern cannot be handled by fix_xtrans_ac.\n"), "red");
		return 1;
	}

	// Loop through sample rectangle and count/sum AF and non-AF pixels.
	for (unsigned int y = sam.y; y <= (sam.y + sam.h); y++) {
		for (unsigned int x = sam.x; x <= (sam.x + sam.w); x++) {
			float pixel = fit->type == DATA_FLOAT ?
							fbuf[x + y * fit->rx] :
							(float) buf[x + y * fit->rx];

			switch (get_pixel_type(af, x, y, &af_matrix)) {
			case 'G': // This is a Green (non-AF) pixel.
				nfcount++;
				nfsum += (double) pixel;
				break;
			case 'g': // This is a Green AF pixel.
				afcount++;
				afsum += (double) pixel;
				break;
			default: // We don't care about other colors... yet.
				break;
			}

		}
	}

	// Make sure we have a valid sample.
	if (nfcount == 0 || afcount == 0) {
		siril_log_message(_("Failed to sample enough pixels.\n"));
		return -1.f;
	}

	// Compute averages and fudge amount.
	nfmean = nfsum / nfcount;
	afmean = afsum / afcount;
	fudge = afmean - nfmean;

	// Debug statements.
	siril_debug_print("XTRANS non-AF Mean... %.10f (%ld pixels)\n", nfmean, nfcount);
	siril_debug_print("XTRANS AF Mean....... %.10f (%ld pixels)\n", afmean, afcount);
	siril_debug_print("XTRANS AF Adjust..... %.10f\n", fudge);

	// Stay FIT, Subtract the fudge!
	subtract_fudge(fit, af, fudge, &af_matrix);

	invalidate_stats_from_fit(fit);

	return 0;
}
