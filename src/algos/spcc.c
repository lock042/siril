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
#include "gui/photometric_cc.h"
#include "registration/matching/misc.h" // for catalogue parsing helpers
#include "photometric_cc.h"
#include "spcc.h"

static gboolean spcc_filters_initialized = FALSE;

// SPCC White Points
const cmsCIExyY Whitepoint_D58 = {0.32598, 0.33532, 1.0}; // Sun as white point, modelled as D58 black body
// TODO: could model this better based on a real solar spectrum

const cmsCIExyY Whitepoint_average_galaxy = {0.345702915, 0.358538597, 1.0}; // D50
// TODO: for testing, D50 is used. Needs replacing with a value computed from real galactic spectra.
/*
const double xpsampled_wl[163] = {378, 380, 382, 384, 386, 388, 390, 392, 394, 396, 398, 400, 402, 404, 406, 408, 410, 412, 414, 416, 418, 420, 422, 424, 426, 428, 430, 432, 434, 436, 438, 440, 442, 444, 446, 448, 450, 452, 454, 456, 458, 460,
  462, 464, 466, 468, 470, 472, 474, 476, 478, 480, 482, 484, 486, 488, 490, 492, 494, 496, 498, 500, 502, 504, 506, 508, 510, 512, 514, 516, 518, 520, 522, 524, 526, 528, 530, 532, 534, 536, 538, 540, 542, 544, 546, 548, 550, 552,
  554, 556, 558, 560, 562, 564, 566, 568, 570, 572, 574, 576, 578, 580, 582, 584, 586, 588, 590, 592, 594, 596, 598, 600, 602, 604, 606, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 628, 630, 632, 634, 636, 638, 640, 642, 644,
  646, 648, 650, 652, 654, 656, 658, 660, 662, 664, 666, 668, 670, 672, 674, 676, 678, 680, 682, 684, 686, 688, 690, 692, 694, 696, 698, 700, 702};
*/
enum {
	RED, GREEN, BLUE
};

enum {
	CMF_1931,
	CMF_1964
};
/*
static void init_spcc_filters() {
	Optolong_Blue.x = Optolong_Blue_wl;
	Optolong_Blue.y = Optolong_Blue_sr;
	Optolong_Blue.n = 72;
	Optolong_Green.x = Optolong_Green_wl;
	Optolong_Green.y = Optolong_Green_sr;
	Optolong_Green.n = 52;
	Optolong_Red.x = Optolong_Red_wl;
	Optolong_Red.y = Optolong_Red_sr;
	Optolong_Red.n = 66;
	Sony_IMX571M.x = Sony_IMX571_wl;
	Sony_IMX571M.y = Sony_IMX571_qe;
	Sony_IMX571M.n = 32;
}
*/
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



