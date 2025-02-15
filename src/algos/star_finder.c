/*
 * This file is part of Siril, an astronomy image processor.
 *
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

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <locale.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "algos/siril_wcs.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/PSF_list.h"
#include "gui/utils.h"
#include "registration/registration.h"
#include "opencv/opencv.h"
#include <wcslib.h>
#include <wcsfix.h>

#define _SQRT_EXP1 1.6487212707
#define KERNEL_SIZE 2.  // sigma of the gaussian smoothing kernel
#define DENSITY_THRESHOLD 0.001 // energy density threshold at which the PSF theoretical radius is cut. Value of 0.001 means we set the box radius to enclose 99.9% of the volume below the psf
#define SAT_THRESHOLD 0.7 // fraction of the dynamic range (frame max - bg) above which pixels are expected to saturate or be close to saturation
#define SAT_DETECTION_RANGE 0.1 // fraction of the dynamic range (frame max - bg) below local max value above which the 8 adjacent pixels must remain to consider we have a saturation plateau
#define MAX_BOX_RADIUS 200 // max allowable value for R (the radius of the box that is passed to PSF fitting)
#define MAX_RADIUS_RATIO_DUP 0.2 // The fraction of the box radius to classify as a duplicate

// Use this flag to print canditates rejection output (0 or 1, only works if SIRIL_OUTPUT_DEBUG is on)
#define DEBUG_STAR_DETECTION 0

static float compute_threshold(image *image, double ksigma, int layer, rectangle *area, float *norm, double *bg, double *bgnoise, double *max, int threads) {
	float threshold;
	imstats *stat;

	assert(layer <= 3);

	stat = statistics(image->from_seq, image->index_in_seq, image->fit,
			layer, area, STATS_BASIC, threads);
	if (!stat) {
		siril_log_message(_("Error: statistics computation failed.\n"));
		*norm = 0;
		*bg = 0.0;
		return 0;
	}
	threshold = (float)(stat->median + ksigma * stat->bgnoise);
	*norm = stat->normValue;
	*bg = stat->median;
	*bgnoise = stat->bgnoise;
	*max = stat->max;
	free_stats(stat);

	return threshold;
}

static sf_errors reject_star(psf_star *result, const star_finder_params *sf, starc *se, double dynrange, double minA, double maxA, gchar *errmsg) {
	if (isnan(result->fwhmx) || isnan(result->fwhmy))
		return SF_NO_FWHM; //crit 11
	if (isnan(result->x0) || isnan(result->y0))
		return SF_NO_POS; //crit 12
	if (isnan(result->mag))
		return SF_NO_MAG; //crit 13
	if (se->iscolor && (result->fwhmx <= 1.0 || result->fwhmy <= 1.0)) {
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "fwhmx: %3.1f, fwhmy: %3.1f\n", result->fwhmx, result->fwhmy);
		return SF_FWHM_TOO_SMALL; //crit 14
	}
	if (result->fwhmx <= 0.0 || result->fwhmy <= 0.0) {
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "fwhmx: %3.1f, fwhmy: %3.1f\n", result->fwhmx, result->fwhmy);
		return SF_FWHM_NEG; //crit 15
	}
	double r = result->fwhmy / result->fwhmx;
	if (r < sf->roundness || r > sf->max_r) {
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "fwhmx: %3.1f, fwhmy: %3.1f\n", result->fwhmx, result->fwhmy);
		return SF_ROUNDNESS_BELOW_CRIT; //crit 16
	}
	// if ((fabs(result->x0 - (double)se->R) >= se->sx) || (fabs(result->y0 - (double)se->R) >= se->sy)) {
	// 	// if star center off from original candidate detection by more than sigma radius
	// 	if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "x0: %3.1f, y0: %3.1f, R:%d\n", result->x0, result->y0, se->R);
	// 	return SF_CENTER_OFF; //crit 10
	// }
	if (result->fwhmx > max(se->sx, se->sy) * _2_SQRT_2_LOG2 * (1 + 0.5 * log(max(se->sx, se->sy) / KERNEL_SIZE))) {
		// criteria gets looser as guessed fwhm gets larger than kernel
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "fwhm: %3.1f, s: %3.1f, m: %3.1f, R: %3d\n", result->fwhmx, max(se->sx, se->sy), _2_SQRT_2_LOG2 * (1 + 0.5 * log(max(se->sx, se->sy) / KERNEL_SIZE)), se->R);
		return SF_FWHM_TOO_LARGE; //crit 2
	}
	if (((result->rmse * sf->sigma / result->A) > 0.2) && (!(result->A > dynrange))) {
		//  do not apply for saturated stars to keep them for alignement purposes
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "RMSE: %4.3e, A: %4.3e, B: %4.3e\n", result->rmse, result->A, result->B);
		return SF_RMSE_TOO_LARGE; //crit 3
	}
	if ((minA > 0.0 || maxA > 0.0) && (result->A < minA || result->A > maxA)) {
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "A: %4.3e, allowed [%4.3f, %4.3f]\n", result->A, minA, maxA);
		return SF_AMPLITUDE_OUTSIDE_RANGE;
	}
	if ((result->profile != PSF_GAUSSIAN) && (result->beta <  sf->min_beta)) {
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "Beta: %4.3e\n", result->beta);
		return SF_MOFFAT_BETA_TOO_SMALL; // crit 17
	}
	return SF_OK;
}

static int star_cmp_by_mag_est(const void *a, const void *b) {
	starc *a1 = (starc *)a;
	starc *a2 = (starc *)b;
	if ((*a1).mag_est > (*a2).mag_est)
		return -1;
	else if ((*a1).mag_est < (*a2).mag_est)
		return 1;
	return 0;
}

// Finds if the point xx,yy is already listed in the candidates within matchradius
// The list is crawled in reverse order to return TRUE faster
// When the next candidate is further away (along y axis) than the current box radius, the search breaks and returns FALSE
static gboolean candidate_is_duplicate(int xx, int yy, int boxradius, starc *candidates, int nbstars) {
	int matchradius = (int)((double)boxradius * MAX_RADIUS_RATIO_DUP);
	matchradius = max(matchradius, 1); // just in case we find zero
					   // we go though the list in reverse order as we may break out faster
	for (int i = nbstars - 1; i >= 0; i--) {
		if (abs(xx - candidates[i].x) + abs(yy - candidates[i].y) <= matchradius)
			return TRUE;
		if (yy - candidates[i].y > boxradius) // if yy is further away than the box radius, we won't find a duplicate anymore
			break;
	}
	return FALSE;
}

/* peaker is the function that searches for stars in an image.
 * It is based on two old implementations that have been refined here over time:
 * 1. Copyleft (L) 1998 Kenneth J. Mighell (Kitt Peak National Observatory)
 * 2. DAOFIND by Peter Stetson, 1987.
 * This is a peak detector on a smoother image (now using Gaussian blur) which
 * identifies any pixel greater than its eight neighbors. The candidates are
 * then checked for consistency before being fitted to Gaussian or Moffat star
 * profiles (PSF).
 */

static int minimize_candidates(fits *image, star_finder_params *sf, starc *candidates, int nb_candidates, int layer, double dynrange, psf_star ***retval, gboolean limit_nbstars, int maxstars, starprofile profile, int threads);

psf_star **peaker(image *image, int layer, star_finder_params *sf, int *nb_stars, rectangle *area, gboolean showtime, gboolean limit_nbstars, int maxstars, starprofile profile, int threads) {
	int nx = image->fit->rx;
	int ny = image->fit->ry;
	int areaX0 = 0;
	int areaY0 = 0;
	int areaX1 = nx;
	int areaY1 = ny;
	int nbstars = 0;
	double bg, bgnoise, maxi;
	float threshold, norm;
	float **smooth_image;
	fits smooth_fit = { 0 };
	starc *candidates;
	struct timeval t_start, t_end;
	WORD **image_ushort = NULL;
	float **image_float = NULL;

	assert(nx > 0 && ny > 0);
	gboolean ismono = (image->fit->naxes[2] == 1);

	if (showtime) {
		siril_log_color_message(_("Findstar: processing for channel %d...\n"), "green", layer);
		gettimeofday(&t_start, NULL);
	}
	else siril_log_message(_("Findstar: processing for channel %d...\n"), layer);

	/* running statistics on the input image is best as it caches them */
	threshold = compute_threshold(image, sf->sigma * 5.0, layer, area, &norm, &bg, &bgnoise, &maxi, threads);
	if (norm == 0.0f)
		return NULL;

	siril_debug_print("Threshold: %f (background level: %f, noise: %f, norm: %f)\n", threshold, bg, bgnoise, norm);

	/* Applying a Gaussian filter to select candidates */
	if (extract_fits(image->fit, &smooth_fit, layer, TRUE)) {
		siril_log_color_message(_("Failed to copy the image for processing\n"), "red");
		return NULL;
	}

	//if (cvUnsharpFilter(&smooth_fit, 3, 0)) {
	if (gaussian_blur_RT(&smooth_fit, KERNEL_SIZE, threads)) {
		siril_log_color_message(_("Could not apply Gaussian filter, aborting\n"), "red");
		clearfits(&smooth_fit);
		return NULL;
	}

	/* Build 2D representation of smoothed image upside-down */
	smooth_image = malloc(ny * sizeof(float *));
	if (!smooth_image) {
		PRINT_ALLOC_ERR;
		clearfits(&smooth_fit);
		return NULL;
	}
	for (int k = 0; k < ny; k++) {
		smooth_image[ny - k - 1] = smooth_fit.fdata + k * nx;
	}

	if (area && area->w != 0 && area->h != 0) {
		areaX0 = area->x;
		areaY0 = area->y;
		areaX1 = area->w + areaX0;
		areaY1 = area->h + areaY0;

		if (areaX1 > nx || areaY1 > ny) {
			siril_log_color_message(_("Selection is larger than image\n"), "red");
			clearfits(&smooth_fit);
			free(smooth_image);
			return NULL;
		}
	}

	/* Build 2D representation of input image upside-down */
	data_type itype = image->fit->type;
	if (itype == DATA_USHORT) {
		image_ushort = malloc(ny * sizeof(WORD *));
		for (int k = 0; k < ny; k++)
			image_ushort[ny - k - 1] = image->fit->pdata[layer] + k * nx;	}
	else if (itype == DATA_FLOAT) {
		image_float = malloc(ny * sizeof(float *));
		for (int k = 0; k < ny; k++)
			image_float[ny - k - 1] = image->fit->fpdata[layer] + k * nx;
	}
	else return 0;

	candidates = malloc(MAX_STARS * sizeof(starc));
	if (!candidates) {
		free(image_ushort);
		free(image_float);
		clearfits(&smooth_fit);
		free(smooth_image);
		PRINT_ALLOC_ERR;
		return NULL;
	}

	int r = sf->radius;
	int boxsize = (2 * r + 1)*(2 * r + 1);
	double locthreshold = sf->sigma * 5.0 * bgnoise;
	double sat = 0.;
	double dynrange = (min(maxi, norm) - bg);
	double minsatlevel = dynrange * SAT_THRESHOLD; // the level above background at which pixels are considered to have saturated or be close to saturation
	double satrange = dynrange * SAT_DETECTION_RANGE; // the max variation level that defines the plateau of saturation
	double s_factor = sqrt(-2. * log(DENSITY_THRESHOLD));
	siril_debug_print("Min saturation level: %3.1f\n", minsatlevel);
	/* Search for candidate stars in the filtered image */
	for (int y = r + areaY0; y < areaY1 - r; y++) {
		for (int x = r + areaX0; x < areaX1 - r; x++) {
			float pixel = smooth_image[y][x];
			float pixel0 = 0.f;
			if (pixel > threshold) {
				gboolean bingo = TRUE;
				float neighbor;
				double meanhigh = 0., minhigh = DBL_MAX;
				int count = 0;
				// making sure the central pixel is a local max compared to all neighbors in the search box
				for (int yy = y - r; yy <= y + r; yy++) {
					for (int xx = x - r; xx <= x + r; xx++) {
						if (xx == x && yy == y)
							continue;
						neighbor = smooth_image[yy][xx];
						if (neighbor > pixel) {
							bingo = FALSE;
							break;
						} else if (neighbor == pixel) {
							if ((xx <= x && yy <= y) || (xx > x && yy < y)) {
								bingo = FALSE;
								break;
							}
						}
						if (!bingo) break;
						count ++;
					}
				}
				if (count < boxsize - 1) {
					bingo = FALSE;
					continue;
				}
				x += r; //no other local max can be found in the next r pixels band

				count = 0;
				for (int yy = y - 1; yy <= y + 1; yy++) {
					for (int xx = x - r - 1; xx <= x - r + 1; xx++) { // x corrected by the anticipated shift
						if (xx == x - r && yy == y) {
							pixel0 = (itype == DATA_USHORT) ? (float)image_ushort[yy][xx] : image_float[yy][xx];
							continue;
						}
						neighbor = (itype == DATA_USHORT) ? (float)image_ushort[yy][xx] : image_float[yy][xx];
						if (neighbor >= threshold) {
							if (neighbor < minhigh) minhigh = neighbor;
							meanhigh += neighbor;
							count++;
						}
					}
				}
				// if the image is mono, there is no leakage possible between channels, so we consider a block of 4 a valid candidate
				// if the image is debayered from color, we need to find a block of 9
				if (count == 0 || (ismono && count < 3) || (!ismono && count < 8)) {
					bingo = FALSE;
					continue;
				}
				meanhigh /= (double)count;

				// first derivatives. r c is for row and column. l, r, u, d for left, right, up and down
				int xx = x - r, yy = y;
				int xr = 0, xl = 0, yu = 0, yd = 0;
				// siril_debug_print("%d: %d - %d\n", nbstars, xx, yy);
				float d1rl, d1rr, d1cu, d1cd, r0, c0;
				int i = 0, j = 0;
				gboolean has_saturated = FALSE;

				// detect if we've come too close to the edges
				if (xx - 2 < areaX0 + 1 || xx + 2 > areaX1 - 1 || yy - 2 < areaY0 + 1 || yy + 2 > areaY1 - 1) continue;

				if (meanhigh - bg < minsatlevel || pixel0 - minhigh > satrange) {
					d1rl = pixel - smooth_image[yy][xx - 1];
					d1rr = smooth_image[yy][xx + 1] - pixel;
					d1cu = pixel - smooth_image[yy - 1][xx];
					d1cd = smooth_image[yy + 1][xx] - pixel;
					// computing zero-crossing to guess correctly psf widths later on
					r0 = -0.5 - d1rl/(d1rr - d1rl);
					c0 = -0.5 - d1cu/(d1cd - d1cu);
					// siril_debug_print("%.2f, %.2f\n\n", r0, c0);
				} else {
					// we need to refine the center position when pixels have saturated
					// to make sure the sanity checks do not fail
					// and to correctly center the box for psf fitting
					has_saturated = TRUE;
					sat = min(pixel0, norm) - satrange;
					while ((xx + i < areaX1 - 1) && smooth_image[yy][xx + i] > sat) i++; // we move right to find the edge
					// and now we'll do some edge-walking
					// we move SW or S
					while ((xx + i < areaX1 - 1) && (yy + j < areaY1 - 1)
						&& (smooth_image[yy + j + 1][xx + i + 1] > sat
						||  smooth_image[yy + j + 1][xx + i    ] > sat)) {
						if (smooth_image[yy + j + 1][xx + i + 1] > sat) i++;
						j++;
					}
					xr = i; // we've reached the Western border
					// we move SE or E
					while ((xx + i > areaX0 + 1) && (yy + j < areaY1 - 1)
						&& (smooth_image[yy + j + 1][xx + i - 1] > sat
						||  smooth_image[yy + j    ][xx + i - 1] > sat)) {
						if (smooth_image[yy + j + 1][xx + i - 1] > sat) j++;
						i--;
					}
					yd = j; // we've reached the Southern border
					// we move NE or N
					while ((xx + i > areaX0 + 1) && (yy + j > areaY0 + 1)
						&& (smooth_image[yy + j - 1][xx + i - 1] > sat
						||  smooth_image[yy + j - 1][xx + i    ] > sat)) {
						if (smooth_image[yy + j - 1][xx + i - 1] > sat) i--;
						j--;
					}
					xl = i; // we've reached the Eastern border
					// we move NW or W
					while ((xx + i < areaX1 - 1) && (yy + j > areaY0 + 1)
						&& (smooth_image[yy + j - 1][xx + i + 1] > sat
						||  smooth_image[yy + j    ][xx + i + 1] > sat)) {
						if (smooth_image[yy + j - 1][xx + i + 1] > sat) j--;
						i++;
					}
					yu = j; // we've reached the Northern border

					// we move SW or S again to be sure we close the loop
					while ((xx + i < areaX1 - 1) && (yy + j < areaY1 - 1)
						&& (smooth_image[yy + j + 1][xx + i + 1] > sat
						||  smooth_image[yy + j + 1][xx + i    ] > sat)) {
						if (smooth_image[yy + j + 1][xx + i + 1] > sat) i++;
						j++;
					}
					if (i > xr) xr = i;
					xx += (xr + xl) / 2;
					yy += (yu + yd) / 2;
					r0 = -0.5;
					c0 = -0.5;
					x += xr;
				}

				float d2rr, d2rrr, d2rl, d2rll, d2cu, d2cuu, d2cd, d2cdd; // second dev
				float srr, srl, scd, scu;
				float Arr, Arl, Acd, Acu;

				// Computing zero-upcrossing of the 2nd deriv row-wise moving right
				i = 0;
				if (has_saturated) while (xx + i < nx && smooth_image[yy][xx + i] > sat) i++;
				if (xx + i >= nx - 2) { // largely saturated star close to border - no chances the fit will be meaningful - discarding
					bingo = FALSE;
					continue;
				}
				d2rr  = smooth_image[yy][xx + i + 1] + smooth_image[yy][xx + i - 1] - 2 * smooth_image[yy][xx + i    ];
				d2rrr = smooth_image[yy][xx + i + 2] + smooth_image[yy][xx + i    ] - 2 * smooth_image[yy][xx + i + 1];
				while ((d2rrr < 0) && ((xx + i + 2) < areaX1 - 1)) {
					i++;
					d2rr = d2rrr;
					d2rrr = smooth_image[yy][xx + i + 2] + smooth_image[yy][xx + i] - 2 * smooth_image[yy][xx + i + 1];
				}
				srr = i  -  d2rr / (d2rrr - d2rr) - r0; // right width estimate
				d1rr = (smooth_image[yy][xx + i] - smooth_image[yy][xx + i - 1]); // first deriv at 2nd deriv upcrossing
				Arr =  -d1rr * srr * _SQRT_EXP1; // right Amplitude estimate

				// Computing zero-upcrossing of the 2nd deriv row-wise moving left
				i = 0;
				if (has_saturated) while (xx - i > 0 && smooth_image[yy][xx - i] > sat) i++;
				if (xx - i <= 2) { // largely saturated star close to border - no chances the fit will be meaningful - discarding
					bingo = FALSE;
					continue;
				}
				d2rl  = smooth_image[yy][xx - i - 1] + smooth_image[yy][xx - i + 1] - 2 * smooth_image[yy][xx - i    ];
				d2rll = smooth_image[yy][xx - i - 2] + smooth_image[yy][xx - i    ] - 2 * smooth_image[yy][xx - i - 1];
				while ((d2rll < 0) && ((xx - i - 2) > areaX0 + 1)) {
					i++;
					d2rl = d2rll;
					d2rll = smooth_image[yy][xx - i - 2] + smooth_image[yy][xx - i] - 2 * smooth_image[yy][xx - i - 1];
				}
				srl = -(i - d2rl / (d2rll - d2rl)) - r0; // left width estimate
				d1rl = (smooth_image[yy][xx - i] - smooth_image[yy][xx - i + 1]); // first deriv at 2nd deriv upcrossing
				Arl =  d1rl * srl * _SQRT_EXP1; // left Amplitude estimate

				// Computing zero-upcrossing of the 2nd deriv column-wise moving down
				i = 0;
				if (has_saturated) while (yy + i < ny && smooth_image[yy + i][xx] > sat) i++;
				if (yy + i >= ny - 2) { // largely saturated star close to border - no chances the fit will be meaningful - discarding
					bingo = FALSE;
					continue;
				}
				d2cd  = smooth_image[yy + i + 1][xx] + smooth_image[yy + i - 1][xx] - 2 * smooth_image[yy + i    ][xx];
				d2cdd = smooth_image[yy + i + 2][xx] + smooth_image[yy + i    ][xx] - 2 * smooth_image[yy + i + 1][xx];
				while ((d2cdd < 0) && ((yy + i + 2) < areaY1 - 1)) {
					i++;
					d2cd = d2cdd;
					d2cdd = smooth_image[yy + i + 2][xx] + smooth_image[yy + i][xx] - 2 * smooth_image[yy + i + 1][xx];
				}
				scd = i  -  d2cd / (d2cdd - d2cd) - c0; // down width estimate
				d1cd = (smooth_image[yy + i][xx] - smooth_image[yy + i - 1][xx]); // first deriv at 2nd deriv upcrossing
				Acd =  -d1cd * scd * _SQRT_EXP1; // down Amplitude estimate

				// Computing zero-upcrossing of the 2nd deriv column-wise moving up
				i = 0;
				if (has_saturated) while (yy - i > 0 && smooth_image[yy - i][xx] > sat) i++;
				if (yy - i <= 2) { // largely saturated star close to border - no chances the fit will be meaningful - discarding
					bingo = FALSE;
					continue;
				}
				d2cu =  smooth_image[yy - i - 1][xx] + smooth_image[yy - i + 1][xx] - 2 * smooth_image[yy - i    ][xx];
				d2cuu = smooth_image[yy - i - 2][xx] + smooth_image[yy - i    ][xx] - 2 * smooth_image[yy - i - 1][xx];
				while ((d2cuu < 0) && ((yy - i - 2) > areaY0 + 1)) {
					i++;
					d2cu = d2cuu;
					d2cuu = smooth_image[yy - i - 2][xx] + smooth_image[yy - i][xx] - 2 * smooth_image[yy - i - 1][xx];
				}
				scu = -(i - d2cu / (d2cuu - d2cu)) - c0; // up width estimate
				d1cu = (smooth_image[yy - i][xx] - smooth_image[yy - i + 1][xx]); // first deriv at 2nd deriv upcrossing
				Acu =  d1cu * scu * _SQRT_EXP1; // up Amplitude estimate

				// computing smoothed psf estimators
				float Sr, Ar, Sc, Ac;
				Sr = 0.5 * (-srl + srr);
				Ar = 0.5 * (Arl + Arr);
				Sc = 0.5 * (-scu + scd);
				Ac = 0.5 * (Acu + Acd);

				// restimate new box size to enclose enough background
				int Rr = (int) ceil(s_factor * Sr);
				int Rc = (int) ceil(s_factor * Sc);
				int Rm = max(Rr, Rc);
				Rm = min(Rm, MAX_BOX_RADIUS);
				int R = max(Rm, r);

				// avoid enlarging outside frame width
				if (xx - R < 0)
					R = xx;
				if (xx + R >= nx)
					R = nx - xx - 1;
				// avoid enlarging outside frame height
				if (yy - R < 0)
					R = yy;
				if (yy + R >= ny)
					R = ny - yy - 1;

				// Quality checks
				float dA = max(Ar,Ac)/min(Ar,Ac);
				float dSr = max(-srl,srr)/min(-srl,srr);
				float dSc = max(-scu,scd)/min(-scu,scd);
				if (!com.pref.starfinder_conf.relax_checks)
					if (dA > 2. || dSr > 2. || dSc > 2. || max(Ar,Ac) < locthreshold)
						bingo = FALSE;

				if (bingo && nbstars < MAX_STARS) {
					if (nbstars > 0 && candidate_is_duplicate(xx, yy, R, candidates, nbstars)) {
						if (DEBUG_STAR_DETECTION)
							siril_debug_print("candidate is a duplicate\n");
						continue; // avoid duplicates for large saturated stars
					}
					candidates[nbstars].x = xx;
					candidates[nbstars].y = yy;
					candidates[nbstars].mag_est = meanhigh;
					candidates[nbstars].R = R; // revised box size
					candidates[nbstars].sx = Sr;
					candidates[nbstars].sy = Sc;
					candidates[nbstars].sat = (has_saturated) ? sat : norm;
					candidates[nbstars].has_saturated = (has_saturated);
					candidates[nbstars].iscolor = !ismono;
					nbstars++;
					if (has_saturated && DEBUG_STAR_DETECTION)
						siril_debug_print("%d: %d - %d is saturated with R = %d\n",
								nbstars, xx, yy, R);
					if (nbstars == MAX_STARS) break;
				}
			}
		}
		if (nbstars == MAX_STARS) break;
	}
	free(smooth_image);
	clearfits(&smooth_fit);
	siril_debug_print("Candidates for stars: %d\n", nbstars);
	/* Check if candidates are stars by minimizing a PSF on each */
	psf_star **results;
	nbstars = minimize_candidates(image->fit, sf, candidates, nbstars, layer, dynrange, &results, limit_nbstars, maxstars, profile, threads);
	if (nbstars == 0)
		results = NULL;
	sort_stars_by_mag(results, nbstars);
	free(candidates);

	if (showtime) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}

	if (image_ushort) free(image_ushort);
	if (image_float) free(image_float);

	if (nb_stars)
		*nb_stars = nbstars;
	return results;
}

/* returns number of stars found, result is in parameters */
static int minimize_candidates(fits *image, star_finder_params *sf, starc *candidates, int nb_candidates, int layer, double dynrange, psf_star ***retval, gboolean limit_nbstars, int maxstars, starprofile profile, int threads) {
	int nx = image->rx;
	int ny = image->ry;
	WORD **image_ushort = NULL;
	float **image_float = NULL;
	gint nbstars = 0;
	sf_errors accepted_level = (com.pref.starfinder_conf.relax_checks) ? SF_RMSE_TOO_LARGE : SF_OK;
	double bg = background(image, layer, NULL, SINGLE_THREADED);
	int psf_failure = 0;
	double minA = 0.0, maxA = 0.0;
	if (sf->min_A > 0.0 || sf->max_A > 0.0) {
		if (image->type == DATA_USHORT) {
			double mult = (image->bitpix == BYTE_IMG) ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
			minA = sf->min_A * mult;
			maxA = sf->max_A * mult;
		}
		else if (image->type == DATA_FLOAT) {
			minA = sf->min_A;
			maxA = sf->max_A;
		}
		else return 0;
		siril_debug_print("using the amplitude range [%f, %f]\n", minA, maxA);
	}

	if (image->type == DATA_USHORT) {
		image_ushort = malloc(ny * sizeof(WORD *));
		for (int k = 0; k < ny; k++)
			image_ushort[ny - k - 1] = image->pdata[layer] + k * nx;
	}
	else {
		image_float = malloc(ny * sizeof(float *));
		for (int k = 0; k < ny; k++)
			image_float[ny - k - 1] = image->fpdata[layer] + k * nx;
	}

	psf_star **results = new_fitted_stars(nb_candidates);
	if (!results) {
		PRINT_ALLOC_ERR;
		if (image_float)
			free(image_float);
		if (image_ushort)
			free(image_ushort);
		return 0;
	}

	//sorting candidates by starc.mag_est values
	qsort(candidates, nb_candidates, sizeof(starc), star_cmp_by_mag_est);

	int round = 0;
	if (limit_nbstars)
		siril_debug_print("limiting stars to %d for %d candidates\n", maxstars, nb_candidates);
	else siril_debug_print("not limiting stars for %d candidates\n", nb_candidates);
	int number_per_round = limit_nbstars ? maxstars + maxstars / 4 : nb_candidates;
	int lower_limit_for_this_round, upper_limit;
	do {
		lower_limit_for_this_round = round * number_per_round;
		upper_limit = (round + 1) * number_per_round;
		if (upper_limit > nb_candidates)
			upper_limit = nb_candidates;
		//siril_debug_print("round %d from %d to %d candidates\n", round, lower_limit_for_this_round, upper_limit);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1) shared(psf_failure)
#endif
		for (int candidate = lower_limit_for_this_round; candidate < upper_limit; candidate++) {
			if (limit_nbstars && nbstars >= maxstars)
				continue;
			int x = candidates[candidate].x, y = candidates[candidate].y;
			int R = candidates[candidate].R;
			int ii, jj, i, j;
			gsl_matrix *z = gsl_matrix_alloc(R * 2 + 1, R * 2 + 1);
			/* FILL z */
			if (image->type == DATA_USHORT) {
				for (jj = 0, j = y - R; j <= y + R; j++, jj++) {
					for (ii = 0, i = x - R; i <= x + R; i++, ii++) {
						gsl_matrix_set(z, jj, ii, (double)image_ushort[j][i]);
					}
				}
			} else {
				for (jj = 0, j = y - R; j <= y + R; j++, jj++) {
					for (ii = 0, i = x - R; i <= x + R; i++, ii++) {
						gsl_matrix_set(z, jj, ii, (double)image_float[j][i]);
					}
				}
			}
			psf_error error;
			psf_star *cur_star = psf_global_minimisation(z, bg, candidates[candidate].sat, com.pref.starfinder_conf.convergence, TRUE, FALSE, NULL, FALSE, profile, &error);
			gsl_matrix_free(z);
			if (cur_star) {
				gchar errmsg[SF_ERRMSG_LEN] = "";
				sf_errors star_invalidated = reject_star(cur_star, sf, &candidates[candidate], dynrange, minA, maxA, (DEBUG_STAR_DETECTION) ? errmsg : NULL);
				if (error != PSF_ERR_DIVERGED && star_invalidated <= accepted_level) { // we don't return NULL on convergence errors so we need to catch that PSF has not diverged
					cur_star->layer = layer;
					cur_star->xpos = (x - R) + cur_star->x0;
					cur_star->ypos = (y - R) + cur_star->y0;
					cur_star->sat = candidates[candidate].sat;
					cur_star->R = candidates[candidate].R;
					cur_star->has_saturated = (cur_star->A > dynrange);
#if DEBUG_STAR_DETECTION
					if (star_invalidated > SF_OK)
						siril_debug_print("Candidate #%5d: X: %5d, Y: %5d - criterion #%2d failed (but star kept)\n%s", candidate, x, y, star_invalidated, errmsg);
#endif
					if (threads > 1)
						results[candidate] = cur_star;
					else results[nbstars++] = cur_star;
				} else {
#if DEBUG_STAR_DETECTION
					siril_debug_print("Candidate #%5d: X: %5d, Y: %5d - criterion #%2d failed\n%s", candidate, x, y, star_invalidated, errmsg);
#endif
					free_psf(cur_star);
					if (threads > 1)
						results[candidate] = NULL;
				}
			} else {
				if (threads > 1) results[candidate] = NULL;
#if DEBUG_STAR_DETECTION
				siril_debug_print("Candidate #%5d: X: %5d, Y: %5d - PSF fit failed with error %d\n", candidate, x, y, error);
#endif
				g_atomic_int_inc(&psf_failure);
			}
		}
		if (threads > 1) {
			// we kept the candidates at the same indices to keep the list ordered, now we compact it
			for (int candidate = lower_limit_for_this_round; candidate < upper_limit; candidate++) {
				if (limit_nbstars && nbstars >= maxstars) {
					if (results[candidate])
						free_psf(results[candidate]);
					continue;
				}
				if (results[candidate] && candidate >= nbstars)
					results[nbstars++] = results[candidate];
			}
		}
		results[nbstars] = NULL;
		//siril_debug_print("after round %d, found %d stars\n", round, nbstars);
		round++;
	} while (limit_nbstars && nbstars < maxstars && upper_limit < nb_candidates);
	float psf_failure_rate = (float)psf_failure / (float)upper_limit;
	if (psf_failure_rate > 0.5)
		siril_log_color_message(_("More than half of PSF fits have failed - try increasing the convergence criterion\n"), "red");
	if (retval)
		*retval = results;
	if (image_ushort) free(image_ushort);
	if (image_float) free(image_float);
	return nbstars;
}

int compare_stars_by_mag(const void* star1, const void* star2) {
	const psf_star *s1 = *(psf_star**) star1;
	const psf_star *s2 = *(psf_star**) star2;
	if (s1->mag < s2->mag)
		return -1;
	if (s1->mag > s2->mag)
		return 1;
	return 0;
}

void sort_stars_by_mag(psf_star **stars, int total) {
	if (stars)
		qsort(stars, total, sizeof(psf_star*), compare_stars_by_mag);
}

/* allocates a new psf_star structure with a size of n + 1. First element is initialized to NULL */
psf_star **new_fitted_stars(size_t n) {
	psf_star **stars = malloc((n + 1) * sizeof(psf_star *));
	if (stars) stars[0] = NULL;

	return stars;
}

void free_fitted_stars(psf_star **stars) {
	int i = 0;
	while (stars && stars[i])
		free_psf(stars[i++]);
	free(stars);
}

psf_star **filter_stars_by_amplitude(psf_star **stars, float threshold, int *nbfilteredstars) {
	int i = 0;
	int nb = 0;
	while (stars && stars[i])
		if (stars[i++]->A >= threshold) nb++;
	*nbfilteredstars = nb;
	if (nb == 0) {
		free_fitted_stars(stars);
		return NULL;
	}
	psf_star **filtered_stars = new_fitted_stars(nb);
	i = 0;
	nb = 0;
	while (stars && stars[i]) {
		if (stars[i]->A >= threshold) {
			filtered_stars[nb] = new_psf_star();
			memcpy(filtered_stars[nb], stars[i], sizeof(psf_star));
			nb++;
		}
		i++;
	}
	filtered_stars[nb] = NULL;
	free_fitted_stars(stars);
	return filtered_stars;
}

void FWHM_stats(psf_star **stars, int nb, int bitpix, float *FWHMx, float *FWHMy, char **units, float *B, float *Acut, double Acutp) {
	*FWHMx = 0.0f;
	*FWHMy = 0.0f;
	*B = 0.0f;
	int n = 0;
	if (stars && stars[0]) {
		double fwhmx = 0.0, fwhmy = 0.0, b = 0.0;
		*units = stars[0]->units;
		for (int i = 0; i < nb; i++) {
			if (!stars[i]->has_saturated) { // removing saturated stars
				fwhmx += stars[i]->fwhmx;
				fwhmy += stars[i]->fwhmy;
				b += stars[i]->B;
				n++;
			}
		}
		n = (n == 0) ? 1 : n;
		*FWHMx = (float)(fwhmx / (double)n);
		*FWHMy = (float)(fwhmy / (double)n);
		*B = (float)(b / (double)n);

		if (Acut) {
			float *A = malloc(nb * sizeof(float));
			if (!A) {
				PRINT_ALLOC_ERR;
				return;
			}
			for (int i = 0; i < nb; i++)
				A[i] = stars[i]->A;
			quicksort_f(A, nb);
			*Acut = (float)gsl_stats_float_quantile_from_sorted_data(A, 1, nb, Acutp);
			free(A);
		}
	}
}

float filtered_FWHM_average(psf_star **stars, int nb) {
	if (!stars || !stars[0])
		return 0.0f;
	float *fwhms = malloc(nb * sizeof(float));
	if (!fwhms) {
		PRINT_ALLOC_ERR;
		return 0.0f;
	}
	for (int i = 0; i < nb; i++)
		fwhms[i] = stars[i]->fwhmx;
	quicksort_f(fwhms, nb);
	float retval = siril_stats_trmean_from_sorted_data(0.15, fwhms, 1, nb);
	free(fwhms);
	return retval;
}

/* looks for the list file, tries to load it, returns TRUE on success */
static gboolean check_star_list(gchar *filename, struct starfinder_data *sfargs) {
	FILE *fd = g_fopen(filename, "r");
	if (!fd)
		return FALSE;
	siril_debug_print("star list file %s found, checking...\n", filename);
	char buffer[300];
	gboolean params_ok = FALSE, read_failure = FALSE, discard_file = FALSE;;
	const star_finder_params *sf = &com.pref.starfinder_conf;
	int star = 0, nb_stars = -1;
	while (fgets(buffer, 300, fd)) {
		if (buffer[0] != '#' && !params_ok) {
			read_failure = TRUE;
			break;
		}
		if (!strncmp(buffer, "# ", 2)) {
			if (sscanf(buffer, "# %d stars found ", &nb_stars) == 1) { // nbstars tainted as based on data from file, needs extra checking
				siril_debug_print("nb stars: %d\n", nb_stars);
				if (nb_stars > 0 && nb_stars < MAX_STARS && sfargs->stars) { // check nbstars against lower and upper bounds
					*sfargs->stars = malloc((nb_stars + 1) * sizeof(struct psf_star *));
					if (!(*sfargs->stars)) {
						PRINT_ALLOC_ERR;
						read_failure = TRUE;
						break;
					}
				}
			}
		}
		if (!strncmp(buffer, "# sigma=", 8)) {
			star_finder_params fparams = { 0 };
			int fmax_stars;
			int prof, layer;
			if (sscanf(buffer, "# sigma=%lf roundness=%lf radius=%d relax=%d profile=%d minbeta=%lf max_stars=%d layer=%d minA=%lf maxA=%lf maxR=%lf",
						&fparams.sigma, &fparams.roundness, &fparams.radius,
						&fparams.relax_checks, &prof, &fparams.min_beta, &fmax_stars,
						&layer, &fparams.min_A, &fparams.max_A, &fparams.max_r) < 8) {
				// minA, maxA and maxR are newer and not mandatory
				read_failure = TRUE;
				break;
			}
			fparams.profile = prof;
			params_ok = fparams.sigma == sf->sigma && fparams.roundness == sf->roundness &&
				fparams.radius == sf->radius &&	fparams.relax_checks == sf->relax_checks &&
				fparams.profile == sf->profile && fparams.min_beta == sf->min_beta &&
				(fmax_stars >= sfargs->max_stars_fitted) && layer == sfargs->layer &&
				fparams.min_A == sf->min_A && fparams.max_A == sf->max_A && fparams.max_r == sf->max_r;
			if (fmax_stars > sfargs->max_stars_fitted) sfargs->max_stars_fitted = fmax_stars;
			siril_debug_print("params check: %d\n", params_ok);
			if (!params_ok) {
				read_failure = TRUE;
				break;
			}
			if (!sfargs->stars) {// if cannot store them, no need to read them
				star = nb_stars; // for message display
				break;
			}
		}
		if (buffer[0] == '#' || buffer[0] == '\0')
			continue;

		if (star >= nb_stars) {
			siril_debug_print("more star in the list than reported\n");
			read_failure = TRUE;
			break;
		}

		psf_star *s = new_psf_star();
		memset(s, 0, sizeof(psf_star));
		int fi, tokens;
		char dump[256];
		tokens = sscanf(buffer,
				"%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%le\t%lf\t%d\t%s",
				&fi, &s->layer, &s->B, &s->A, &s->beta, &s->xpos, &s->ypos,
				&s->fwhmx, &s->fwhmy, &s->fwhmx_arcsec, &s->fwhmy_arcsec, &s->angle, &s->rmse, &s->mag, &s->has_saturated, dump);
		if (tokens != 16) {
			siril_debug_print("malformed line: %s", buffer);
			read_failure = TRUE;
			discard_file = TRUE;
			free(s);
			s = NULL;
			break;
		}
		if (fi != star + 1)
			siril_debug_print("star index mismatch\n");
		s->units = (s->fwhmx_arcsec > 0. && s->fwhmy_arcsec > 0.) ? "\"" : "px";
		(*sfargs->stars)[star++] = s;
		(*sfargs->stars)[star] = NULL;
	}
	if (!read_failure) {
		siril_log_message(_("Found %d stars with same settings in %s, skipping detection\n"), star, filename);
		if (sfargs->nb_stars) *sfargs->nb_stars = star;
	} else {
		if (sfargs->stars) {
			free(*sfargs->stars);
			(*sfargs->stars) = NULL;
		}
	}
	fclose(fd);
	if (discard_file && g_unlink(filename)) {
		siril_debug_print("g_unlink failed\n");
	}
	return !read_failure;
}

#define HANDLE_WRITE_ERR \
	g_warning("%s\n", error->message); \
	g_clear_error(&error); \
	g_object_unref(output_stream); \
	g_object_unref(file); \
	return 1

int save_list(gchar *filename, int max_stars_fitted, psf_star **stars, int nbstars, const star_finder_params *sf, int layer, gboolean verbose) {
	int i = 0;
	GError *error = NULL;
	gchar *dirname = g_path_get_dirname(filename);
	GDir *dir = g_dir_open(dirname, 0, NULL);
	if (!dir && g_mkdir_with_parents(dirname, 0755) < 0) {
		siril_log_color_message(_("Cannot create output folder: %s\n"), "red", dirname);
		g_free(dirname);
		return 1;
	}
	if (dir)
		g_dir_close(dir);
	g_free(dirname);

	GFile *file = g_file_new_for_path(filename);
	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (!output_stream) {
		if (error) {
			siril_log_message(_("Cannot save star list %s: %s\n"), filename, error->message);
			g_clear_error(&error);
		}
		g_object_unref(file);
		return 1;
	}
	if (nbstars <= 0) {
		// unknown by caller
		nbstars = 0;
		if (stars)
			while (stars[nbstars++]);
	}

	char buffer[320];
	char gausstr[] = "Gaussian";
	char moffstr[] = "Moffat";
	char *starprof;
	gdouble beta;
	int len = snprintf(buffer, 320, "# %d stars found using the following parameters:%s", nbstars, SIRIL_EOL);
	if (!g_output_stream_write_all(output_stream, buffer, len, NULL, NULL, &error)) {
		HANDLE_WRITE_ERR;
	}
	len = snprintf(buffer, 320, "# sigma=%3.2f roundness=%3.2f radius=%d relax=%d profile=%d minbeta=%3.1f max_stars=%d layer=%d minA=%3.2f maxA=%3.2f maxR=%3.2f%s",
			sf->sigma, sf->roundness, sf->radius, sf->relax_checks,sf->profile, sf->min_beta, max_stars_fitted, layer, sf->min_A, sf->max_A, sf->max_r, SIRIL_EOL);
	if (!g_output_stream_write_all(output_stream, buffer, len, NULL, NULL, &error)) {
		HANDLE_WRITE_ERR;
	}
	len = snprintf(buffer, 320,
			"# star#\tlayer\tB\tA\tbeta\tX\tY\tFWHMx [px]\tFWHMy [px]\tFWHMx [\"]\tFWHMy [\"]\tangle\tRMSE\tmag\tSat\tProfile\tRA\tDec%s",
			SIRIL_EOL);
	if (!g_output_stream_write_all(output_stream, buffer, len, NULL, NULL, &error)) {
		HANDLE_WRITE_ERR;
	}
	if (stars) {
		while (stars[i]) {
			if (stars[i]->profile == PSF_GAUSSIAN) {
				beta = -1;
				starprof = N_(gausstr);
			} else {
				beta = stars[i]->beta;
				starprof = N_(moffstr);
			}
			len = snprintf(buffer, 320,
					"%d\t%d\t%10.6f\t%10.6f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%3.2f\t%10.3e\t%10.2f\t%d\t%s\t%f\t%f%s",
					i + 1, stars[i]->layer, stars[i]->B, stars[i]->A, beta,
					stars[i]->xpos, stars[i]->ypos, stars[i]->fwhmx,
					stars[i]->fwhmy, stars[i]->fwhmx_arcsec ,stars[i]->fwhmy_arcsec,
					stars[i]->angle, stars[i]->rmse, stars[i]->mag + com.magOffset,
					stars[i]->has_saturated, starprof, stars[i]->ra, stars[i]->dec, SIRIL_EOL);
			if (!g_output_stream_write_all(output_stream, buffer, len, NULL, NULL, &error)) {
				HANDLE_WRITE_ERR;
			}
			i++;
		}
	}
	if (verbose) siril_log_message(_("The file %s has been created.\n"), filename);
	g_object_unref(output_stream);
	g_object_unref(file);

	return 0;
}

/* saving a list of sources to the input format of solve-field, the astrometry.net
 * plate solver, as described here:
 * http://astrometry.net/doc/readme.html#source-lists-xylists
 * Convention is to have the columns for pixel coordinates named IMAGEX and IMAGEY and
 * having in the header the keywords IMAGEW and IMAGEH giving rx and ry.
 * Correction by 0.5 pixel to conform to asnet convention (compared to siril) is also added
 */
int save_list_as_FITS_table(const char *filename, psf_star **stars, int nbstars, int rx, int ry) {
	if (g_unlink(filename)) {
		siril_debug_print("g_unlink failure\n");
	} /* Delete old file if it already exists */

	fitsfile *fptr = NULL;
	int status = 0;
	if (siril_fits_create_diskfile(&fptr, filename, &status)) { /* create new FITS file */
		report_fits_error(status);
		return 1;
	}

	char* rownames[] = { "X", "Y", "FLUX", "BACKGROUND"};
	char* rowtypes[] = { "1E", "1E", "1E", "1E"};	// rE is float, rD is double
	if (fits_create_tbl(fptr, BINARY_TBL, nbstars, 4, rownames, rowtypes, NULL, NULL, &status)) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fptr, &status);
		return 1;
	}

	fits_update_key(fptr, TINT, "IMAGEW", &rx, "Image width in pixels", &status);
	fits_update_key(fptr, TINT, "IMAGEH", &ry, "Image height in pixels", &status);
	if (status)
		siril_debug_print("Failed to write the IMAGEW and IMAGEH headers to the FITS table\n");

	/* reorganize data per row for saving */
	float *data = malloc(nbstars * sizeof(float));
	if (!data) {
		PRINT_ALLOC_ERR;
		fits_close_file(fptr, &status);
		return 1;
	}

	status = 0;
	for (int i = 0; i < nbstars; i++)
		data[i] = stars[i]->xpos + 0.5; // asnet convention
	if (fits_write_col(fptr, TFLOAT, 1, 1, 1, nbstars, data, &status)) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fptr, &status);
		free(data);
		return 1;
	}

	for (int i = 0; i < nbstars; i++)
		data[i] = ry - stars[i]->ypos + 0.5; // asnet convention
	if (fits_write_col(fptr, TFLOAT, 2, 1, 1, nbstars, data, &status)) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fptr, &status);
		free(data);
		return 1;
	}

	for (int i = 0; i < nbstars; i++)
		data[i] = stars[i]->A;
	if (fits_write_col(fptr, TFLOAT, 3, 1, 1, nbstars, data, &status)) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fptr, &status);
		free(data);
		return 1;
	}

	for (int i = 0; i < nbstars; i++)
		data[i] = stars[i]->B;
	if (fits_write_col(fptr, TFLOAT, 4, 1, 1, nbstars, data, &status)) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fptr, &status);
		free(data);
		return 1;
	}

	int retval = status;
	status = 0;
	fits_close_file(fptr, &status);
	free(data);
	return retval;
}

static int findstar_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_image, MB_avail, required = 0;
	// using output scale of 0. as this does not produce output image
	int limit = compute_nb_images_fit_memory(args->seq, 0.0, FALSE, &MB_per_image, NULL, &MB_avail);

	if (limit > 0) {
		int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;
		int float_multiplier = (is_float) ? 1 : 2;
		int chan_multiplier = (args->seq->nb_layers == 3) ? 3 : 1;
		int MB_per_float_image = MB_per_image * float_multiplier;
		int MB_per_float_chan = MB_per_float_image / chan_multiplier;

		// Allocations:
		// ------------
		// * extract_fits() allocates 1 * float_channel in all cases.
		// * As the 1-channel fits populated by extract_fits() is DATA_FLOAT no
		//   memory is required for the Gaussian blur as the in-place RT algorithm is
		//   used;
		// * If the sequence is CFA, we also need to make a copy of the input image
		// to interpolate non green pixels
		// * Also allow MAX_STARS * sizeof(psf_star) + 1MB margin for gslsolver data;
		//   and the indexing arrays e.g. smooth_array (ry * sizeof(float*)).
		int stars_and_overhead = 1 + (int) ceilf(((float) MAX_STARS * (float) sizeof(psf_star)) / (float) BYTES_IN_A_MB);
		required = MB_per_image + MB_per_float_chan + stars_and_overhead;

		// we need to read the bayer pattern of ref image...
		fits fit = { 0 };
		// load ref metadata in fit
		if (seq_read_frame_metadata(args->seq, sequence_find_refimage(args->seq), &fit))
			return 0;
		if (fit.keywords.bayer_pattern[0] != '\0')
			required += MB_per_image;
		clearfits(&fit);

		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
				thread_limit = com.max_thread;
		limit = thread_limit;
	}
	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_image, limit, for_writer ? "images" : "threads");
#else
		limit = 1;
#endif
	}
	return limit;
}

/* return FALSE to avoid reading image */
static gboolean findstar_image_read_hook(struct generic_seq_args *args, int index) {
	struct starfinder_data *findstar_args = (struct starfinder_data *)args->user;

	struct starfinder_data *curr_findstar_args = calloc(1, sizeof(struct starfinder_data));
	memcpy(curr_findstar_args, findstar_args, sizeof(struct starfinder_data));
	curr_findstar_args->im.index_in_seq = index;
	curr_findstar_args->im.fit = NULL;
	if (findstar_args->stars && findstar_args->nb_stars) {
		curr_findstar_args->stars = findstar_args->stars + index;
		curr_findstar_args->nb_stars = findstar_args->nb_stars + index;
	}
	curr_findstar_args->threading = SINGLE_THREADED;

	gchar *star_filename = get_sequence_cache_filename(args->seq, index, "lst", NULL);
	if (!star_filename) {
		free(curr_findstar_args);
		return TRUE;
	}

	if (findstar_args->save_to_file)
		curr_findstar_args->starfile = star_filename;

	gboolean status = check_star_list(star_filename, curr_findstar_args);
	free(curr_findstar_args);
	g_free(star_filename);
	return !status; // check_star_list returns TRUE on success
}

// contrarily to findstar_worker, this function first checks if a lst file exists:
// - if yes and star finder params are consistent, the list is used to populate stars
// - else, it calls peaker
// It is standalone to be:
// - wrapped by findstar_image_hook
// - called by registration
struct starfinder_data *findstar_image_worker(const struct starfinder_data *findstar_args, int o, int i, fits *fit, rectangle *_, int threads) {
	struct starfinder_data *curr_findstar_args = calloc(1, sizeof(struct starfinder_data));
	memcpy(curr_findstar_args, findstar_args, sizeof(struct starfinder_data));
	curr_findstar_args->im.index_in_seq = i;
	curr_findstar_args->im.fit = fit;
	sequence *seq = findstar_args->im.from_seq;
	if (findstar_args->stars && findstar_args->nb_stars) { // used by 2pass reg which needs to store all star lists
		curr_findstar_args->stars = findstar_args->stars + i;
		curr_findstar_args->nb_stars = findstar_args->nb_stars + i;
		curr_findstar_args->onepass = FALSE;
	} else if (findstar_args->keep_stars) { // used by global reg which needs to store the current star list
		curr_findstar_args->stars = calloc(1, sizeof(psf_star **));
		curr_findstar_args->nb_stars = calloc(1, sizeof(int));
		curr_findstar_args->onepass = TRUE;
	}
	curr_findstar_args->threading = threads;
	gboolean can_use_cache = !(com.selection.w != 0 && com.selection.h != 0); //TODO: not ideal, would be better to be passed as args to starfinder_data (to be used in findstar_worker)

	int retval = 0;
	gchar *star_filename = NULL;

	if (can_use_cache) {// otherwise, we don't try to read the lst nor save it
		// build the star list file name in all cases to try reading it
		star_filename = get_sequence_cache_filename(seq, i, "lst", NULL);
		if (!star_filename) {
			if (curr_findstar_args->onepass == TRUE) {
				free(curr_findstar_args->nb_stars);
				free(curr_findstar_args->stars);
			}
			free(curr_findstar_args);
			curr_findstar_args = NULL;
			return curr_findstar_args;
		}

		if (seq->type == SEQ_INTERNAL || !check_cachefile_date(seq, i, star_filename) ||
				!check_star_list(star_filename, curr_findstar_args))
				can_use_cache = FALSE;
		if (findstar_args->save_to_file)
			curr_findstar_args->starfile = star_filename;
		else
			g_free(star_filename);
	}

	if (!can_use_cache) {
		fits *green_fit = NULL;
		if (fit->keywords.bayer_pattern[0] != '\0') {
			green_fit = calloc(1, sizeof(fits));
			copyfits(fit, green_fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
			interpolate_nongreen(green_fit);
			curr_findstar_args->im.fit = green_fit;
			// uncomment these lines to save the green fit for each image
			// const gchar *green_filename = get_sequence_cache_filename(seq, i, "fit", "green_");
			// savefits(green_filename, green_fit);
		}
		retval = GPOINTER_TO_INT(findstar_worker(curr_findstar_args));
		clearfits(green_fit);
		free(green_fit);
		if (retval) {
			if (curr_findstar_args->onepass == TRUE) {
				free(curr_findstar_args->nb_stars);
				free(curr_findstar_args->stars);
			}
			free(curr_findstar_args);
			curr_findstar_args = NULL;
		}
	}
	return curr_findstar_args;
}

static int findstar_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	const struct starfinder_data *findstar_args = (struct starfinder_data *)args->user;
	struct starfinder_data *curr_findstar_args = findstar_image_worker(findstar_args, o, i, fit, _, threads);
	gboolean retval = !curr_findstar_args;
	free(curr_findstar_args);
	return retval;
}

gboolean end_findstar_sequence(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;
	struct starfinder_data *findstar_args = (struct starfinder_data*) args->user;
	if (args->has_output && args->load_new_sequence &&
			args->new_seq_prefix && !args->retval) {
		gchar *basename = g_path_get_basename(args->seq->seqname);
		gchar *seqname = g_strdup_printf("%s%s.seq", args->new_seq_prefix, basename);
		check_seq();
		update_sequences_list(seqname);
		g_free(seqname);
		g_free(basename);
	}
	if (!check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);

	free(findstar_args);
	free(p);
	return end_generic(NULL);
}

int findstar_finalize_hook(struct generic_seq_args *args) {
	struct starfinder_data *data = (struct starfinder_data *) args->user;
	if (data->ref_wcs) {
		if (!wcsfree(data->ref_wcs))
			free(data->ref_wcs);
	}
	if (data->startable)
		g_free(data->startable);
	if (data->starfile)
		g_free(data->starfile);
	printf("findstar_args have been freed\n");
	return 0;
}

int apply_findstar_to_sequence(struct starfinder_data *findstar_args) {
	struct generic_seq_args *args = create_default_seqargs(findstar_args->im.from_seq);
	if (!findstar_args->process_all_images) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = args->seq->selnum;
	}
	args->image_read_hook = findstar_image_read_hook;
	args->image_hook = findstar_image_hook;
	args->finalize_hook = findstar_finalize_hook;
	args->compute_mem_limits_hook = findstar_compute_mem_limits;
	args->stop_on_error = FALSE;
	args->description = _("FindStar");
	args->has_output = FALSE;
	args->load_new_sequence = FALSE;
	args->user = findstar_args;
	args->idle_function = end_findstar_sequence;
	args->already_in_a_thread = findstar_args->already_in_thread;

	if (findstar_args->save_to_file && findstar_args->save_eqcoords) {
		/* saving equatorial coordinates is possible if all images are
		 * plate solved or if the reference is plate solved and
		 * transformation matrices H from registration are available */
		fits ref = { 0 };
		int refidx = sequence_find_refimage(args->seq);
		if (seq_read_frame_metadata(args->seq, refidx, &ref)) {
			siril_log_message(_("Could not load reference image\n"));
			free(args);
			return 1;
		}
		if (!has_wcs(&ref))
			findstar_args->save_eqcoords = FALSE;
		else {
			findstar_args->reference_image = refidx;
			findstar_args->ref_wcs = ref.keywords.wcslib;
			if (args->seq->regparam[findstar_args->layer] &&
					guess_transform_from_H(args->seq->regparam[findstar_args->layer][refidx].H) != NULL_TRANSFORMATION) {
				findstar_args->reference_H = args->seq->regparam[findstar_args->layer][refidx].H;
			}
		}
		ref.keywords.wcslib = NULL;	// don't free it
		clearfits(&ref);
	}
	if (findstar_args->already_in_thread) {
		int retval = GPOINTER_TO_INT(generic_sequence_worker(args));
		free(args);
		return retval;
	}
	if (!start_in_new_thread(generic_sequence_worker, args)) {
		free(args->user);
		free_generic_seq_args(args);
	}
	return 0;
}

gboolean end_findstar(gpointer p);	// in the GUI file

// for a single image
gpointer findstar_worker(gpointer p) {
	struct starfinder_data *args = (struct starfinder_data *)p;
	int retval = 0;
	int nbstars = 0;
	rectangle *selection = NULL;
	if (com.selection.w != 0 && com.selection.h != 0) //TODO: not ideal, would be better to be passed as args to starfinder_data
		selection = &com.selection;
	gboolean limit_stars = (args->max_stars_fitted > 0);
	int threads = check_threading(&args->threading);
	fits *green_fit = NULL;
	if (args->im.fit->keywords.bayer_pattern[0] != '\0') {
		green_fit = calloc(1, sizeof(fits));
		copyfits(args->im.fit, green_fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		interpolate_nongreen(green_fit);
		args->im.fit = green_fit;
		siril_log_color_message(_("Undebayered CFA image. Detection is done on green pixels only, using interpolation\n"), "salmon");
	}
	psf_star **stars = peaker(&args->im, args->layer, &com.pref.starfinder_conf, &nbstars,
			selection, args->update_GUI, limit_stars, args->max_stars_fitted, com.pref.starfinder_conf.profile, threads);
	if (green_fit) {
		clearfits(green_fit);
		free(green_fit);
	}

	double fwhm = 0.0;
	if (stars) {
		int i = 0;
		double sum = 0.0;
		while (stars[i]) {
			sum += stars[i]->fwhmx;
			fwhm_to_arcsec_if_needed(args->im.fit, stars[i]);

			if (args->starfile && args->save_eqcoords) {
				sequence *seq = args->im.from_seq;
				double dx = stars[i]->xpos, dy = stars[i]->ypos;
				if (has_wcs(args->im.fit)) {
					double ra = 0.0, dec = 0.0;
					// coordinates of the star in FITS/WCS coordinates
					double fx, fy;
					display_to_siril(dx, dy, &fx, &fy, args->im.fit->ry);
					pix2wcs2(args->ref_wcs, fx, fy, &ra, &dec);
					// ra and dec = -1 is the error code
					stars[i]->ra = ra;
					stars[i]->dec = dec;
				}
				else if (!seq->regparam[args->layer] ||
						guess_transform_from_H(seq->regparam[args->layer][args->im.index_in_seq].H) == NULL_TRANSFORMATION) {
					// image was not registered, ignore
				}
				else {
					cvTransfPoint(&dx, &dy, seq->regparam[args->layer][args->im.index_in_seq].H, args->reference_H, 1.);
					double ra = 0.0, dec = 0.0;
					// coordinates of the star in FITS/WCS coordinates
					double fx, fy;
					display_to_siril(dx, dy, &fx, &fy, args->im.fit->ry);
					pix2wcs2(args->ref_wcs, fx, fy, &ra, &dec);
					// ra and dec = -1 is the error code
					stars[i]->ra = ra;
					stars[i]->dec = dec;
				}
			}
			i++;
		}
		if (i > 0)
			fwhm = sum / i;
		else fwhm = 0.0;
	} else {
		goto END;
	}

	if (args->update_GUI)
		update_star_list(stars, TRUE, TRUE);

	siril_log_message(_("Found %d %s profile stars in %s, channel #%d (FWHM %f)\n"), nbstars,
			com.pref.starfinder_conf.profile == PSF_GAUSSIAN ? _("Gaussian") : _("Moffat"),
			selection ? _("selection") : _("image"), args->layer, fwhm);
	if (args->starfile &&
			save_list(args->starfile, args->max_stars_fitted, stars, nbstars,
				&com.pref.starfinder_conf, args->layer, args->update_GUI)) {
		retval = 1;
	}

	if (args->startable &&
			save_list_as_FITS_table(args->startable, stars, nbstars,
				args->im.fit->rx, args->im.fit->ry)) {
		siril_debug_print("saving list as FITS table failed\n");
		retval = 1;
	}

	if (args->stars && args->nb_stars) {
		*args->stars = stars;
		*args->nb_stars = nbstars;
	}
	else if (!args->update_GUI)
		free_fitted_stars(stars);
END:
	if (args->update_GUI)
		siril_add_idle(end_findstar, args);
	/*gettimeofday(&t_end, NULL);
	gchar *msg = g_strdup_printf("findstar for image %d", args->im.index_in_seq);
	show_time_msg(t_start, t_end, msg);
	g_free(msg);*/

	return GINT_TO_POINTER(retval);
}

// if channel < 0, returns the max of each channel's robust mean
// if channel, returns the robust mean for this channel
// returns 0 for errors
float measure_image_FWHM(fits *fit, int channel) {
	float fwhm[3];
	image im = { .fit = fit, .from_seq = NULL, .index_in_seq = -1 };
	gboolean failed = FALSE;
	int nb_chan = channel < 0 ? (int)fit->naxes[2] : 1;
	g_assert(nb_chan == 1 || nb_chan == 3);
#ifdef _OPENMP
	int *threads = compute_thread_distribution(nb_chan, com.max_thread);
#pragma omp parallel for num_threads(com.max_thread) if(nb_chan > 1)
#endif
	for (int chan = 0; chan < nb_chan; chan++) {
		int nb_stars;
		int nb_subthreads;
#ifdef _OPENMP
		nb_subthreads = threads[chan];
#else
		nb_subthreads = com.max_thread;
#endif
		int real_chan = channel < 0 ? chan : channel;
		psf_star **stars = peaker(&im, real_chan, &com.pref.starfinder_conf, &nb_stars,
				NULL, FALSE, TRUE, 200, com.pref.starfinder_conf.profile, nb_subthreads);
		if (stars) {
			fwhm[chan] = filtered_FWHM_average(stars, nb_stars);
			siril_debug_print("FWHM for channel %d: %.3f\n", real_chan, fwhm[chan]);

			for (int i = 0; i < nb_stars; i++)
				free_psf(stars[i]);
			free(stars);
		}
		else failed = TRUE;
	}
#ifdef _OPENMP
	free(threads);
#endif
	if (failed)
		return 0.0f;
	if (nb_chan == 1)
		return fwhm[0];
	return max(fwhm[0], max(fwhm[1], fwhm[2]));
}
