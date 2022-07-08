/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <locale.h>
#include <gsl/gsl_matrix.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/Def_Wavelet.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "gui/image_interactions.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "core/OS_utils.h"
#include "opencv/opencv.h"

#define WAVELET_SCALE 3
#define _SQRT_EXP1 1.6487212707
#define KERNEL_SIZE 3.

// Use this flag to print canditates rejection output (0 or 1, only works if SIRIL_OUTPUT_DEBUG is on)
#define DEBUG_STAR_DETECTION 0

static double guess_resolution(fits *fit) {
	double focal = fit->focal_length;
	double size = fit->pixel_size_x;
	double bin;

	/* if we have no way to guess, we return the
		* flag -1
		*/
	if ((focal <= 0.0) || (size <= 0.0)) {  // try to read values that were filled only if coming from astrometry solver
		focal = com.pref.starfinder_conf.focal_length;
		size = com.pref.starfinder_conf.pixel_size_x;
	}

	if ((focal <= 0.0) || (size <= 0.0))
		return -1.0;

	bin = ((fit->binning_x + fit->binning_y) / 2.0);
	if (bin <= 0) bin = 1.0;

	double res = RADCONV / focal * size * bin;

	/* test for high value. In this case we increase
	 * the number of detected star in reject_star function
	 */
	/* if res > 1.0 we use default radius value */
	if (res > 1.0) return 1.0;
	/* if res is too small, we bound the value to 2.0 to avoid expanding the first search box (and decrease perf) */
	if (res < 0.5) return 0.5;
	return res;
}

static float compute_threshold(image *image, double ksigma, int layer, rectangle *area, float *norm, double *bg, double *bgnoise, int threads) {
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
	free_stats(stat);

	return threshold;
}

static sf_errors reject_star(psf_star *result, star_finder_params *sf, starc *se, gchar *errmsg) {
	if (isnan(result->fwhmx) || isnan(result->fwhmy))
		return SF_NO_FWHM;
	if (isnan(result->x0) || isnan(result->y0))
		return SF_NO_POS;
	if (isnan(result->mag))
		return SF_NO_MAG;
	if (result->fwhmx <= 1.0 || result->fwhmy <= 1.0) {
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "fwhmx: %3.1f, fwhmy: %3.1f\n", result->fwhmx, result->fwhmy);
		return SF_FWHM_TOO_SMALL;
	}
	if (result->fwhmx <= 0.0 || result->fwhmy <= 0.0) {
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "fwhmx: %3.1f, fwhmy: %3.1f\n", result->fwhmx, result->fwhmy);
		return SF_FWHM_NEG;
	}
	if ((result->fwhmy / result->fwhmx) < sf->roundness) {
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "fwhmx: %3.1f, fwhmy: %3.1f\n", result->fwhmx, result->fwhmy);
		return SF_ROUNDNESS_BELOW_CRIT;
	}
	if ((fabs(result->x0 - (double)se->R) >= se->sx) || (fabs(result->y0 - (double)se->R) >= se->sy)) { // if star center off from original candidate detection by more than sigma radius
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "x0: %3.1f, y0: %3.1f, R:%d\n", result->x0, result->y0, se->R);
		return SF_CENTER_OFF;
	}
	if (result->fwhmx > max(se->sx, se->sy) * _2_SQRT_2_LOG2 * (1 + 0.5 * log(max(se->sx, se->sy) / KERNEL_SIZE))) {// criteria gets looser as guessed fwhm gets larger than kernel
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "fwhm: %3.1f, s: %3.1f, m: %3.1f, R: %3d\n", result->fwhmx, max(se->sx, se->sy), _2_SQRT_2_LOG2 * (1 + 0.5 * log(max(se->sx, se->sy) / KERNEL_SIZE)), se->R);
		return SF_FWHM_TOO_LARGE;
	}
	if (((result->rmse * sf->sigma / result->A) > 0.1) && (result->A < ((result->B < 1.0) ? 1. : USHRT_MAX_DOUBLE) * 0.5)) {
	//  do not apply for bright stars (above 50% of bitdepth range) to avoid removing bright saturated stars
		if (errmsg) g_snprintf(errmsg, SF_ERRMSG_LEN, "RMSE: %4.3e, A: %4.3e, B: %4.3e\n", result->rmse, result->A, result->B);
		return SF_RMSE_TOO_LARGE;
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

void on_toggle_radius_adjust_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.pref.starfinder_conf.adjust = gtk_toggle_button_get_active(togglebutton);
}

void on_toggle_relax_checks_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.pref.starfinder_conf.relax_checks = gtk_toggle_button_get_active(togglebutton);
}

void on_spin_sf_radius_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.pref.starfinder_conf.radius = (int)gtk_spin_button_get_value(spinbutton);
}

void on_spin_sf_threshold_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.pref.starfinder_conf.sigma = gtk_spin_button_get_value(spinbutton);
}

void on_spin_sf_roundness_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.pref.starfinder_conf.roundness = gtk_spin_button_get_value(spinbutton);
}

void update_peaker_GUI() {
	static GtkSpinButton *spin_radius = NULL, *spin_sigma = NULL,
			*spin_roundness = NULL;
	static GtkToggleButton *toggle_adjust = NULL, *toggle_checks = NULL;

	if (spin_radius == NULL) {
		spin_radius = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_radius"));
		spin_sigma = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_threshold"));
		spin_roundness = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round"));
		toggle_adjust = GTK_TOGGLE_BUTTON(lookup_widget("toggle_radius_adjust"));
		toggle_checks = GTK_TOGGLE_BUTTON(lookup_widget("toggle_relax_checks"));
	}
	gtk_spin_button_set_value(spin_radius, (double) com.pref.starfinder_conf.radius);
	gtk_toggle_button_set_active(toggle_adjust, com.pref.starfinder_conf.adjust);
	gtk_spin_button_set_value(spin_sigma, com.pref.starfinder_conf.sigma);
	gtk_spin_button_set_value(spin_roundness, com.pref.starfinder_conf.roundness);
	gtk_toggle_button_set_active(toggle_checks, com.pref.starfinder_conf.relax_checks);
}

void confirm_peaker_GUI() {
	static GtkSpinButton *spin_radius = NULL, *spin_sigma = NULL,
			*spin_roundness = NULL;

	if (spin_radius == NULL) {
		spin_radius = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_radius"));
		spin_sigma = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_threshold"));
		spin_roundness = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round"));
	}
	gtk_spin_button_update(spin_radius);
	gtk_spin_button_update(spin_sigma);
	gtk_spin_button_update(spin_roundness);
}


/*
 This is an implementation of a simple peak detector algorithm which
 identifies any pixel that is greater than any of its eight neighbors.

 Original algorithm come from:
 Copyleft (L) 1998 Kenneth J. Mighell (Kitt Peak National Observatory)
 */

static int minimize_candidates(fits *image, star_finder_params *sf, starc *candidates, int nb_candidates, int layer, psf_star ***retval, gboolean limit_nbstars, int maxstars, int threads);

psf_star **peaker(image *image, int layer, star_finder_params *sf, int *nb_stars, rectangle *area, gboolean showtime, gboolean limit_nbstars, int maxstars, int threads) {
	int nx = image->fit->rx;
	int ny = image->fit->ry;
	int areaX0 = 0;
	int areaY0 = 0;
	int areaX1 = nx;
	int areaY1 = ny;
	int nbstars = 0;
	double bg, bgnoise;
	float threshold, norm;
	float **smooth_image;
	fits smooth_fit = { 0 };
	starc *candidates;
	struct timeval t_start, t_end;

	assert(nx > 0 && ny > 0);

	if (showtime) {
		siril_log_color_message(_("Findstar: processing for channel %d...\n"), "green", layer);
		gettimeofday(&t_start, NULL);
	}
	else siril_log_message(_("Findstar: processing for channel %d...\n"), layer);

	/* running statistics on the input image is best as it caches them */
	threshold = compute_threshold(image, sf->sigma * 5.0, layer, area, &norm, &bg, &bgnoise, threads);
	if (norm == 0.0f)
		return NULL;

	siril_debug_print("Threshold: %f (background: %f, norm: %f)\n", threshold, bg, norm);

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

	candidates = malloc(MAX_STARS * sizeof(starc));
	if (!candidates) {
		clearfits(&smooth_fit);
		free(smooth_image);
		PRINT_ALLOC_ERR;
		return NULL;
	}

	double res = guess_resolution(image->fit);

	if (res < 0) {
		res = 1.0;
	}

	sf->adj_radius = sf->adjust ? sf->radius / res : sf->radius;
	siril_debug_print("Adjusted radius: %d\n", sf->adj_radius);
	int r = sf->adj_radius;
	int boxsize = (2 * r + 1)*(2 * r + 1);
	double locthreshold = sf->sigma * 5.0 * bgnoise;

	/* Search for candidate stars in the filtered image */
	for (int y = r + areaY0; y < areaY1 - r; y++) {
		for (int x = r + areaX0; x < areaX1 - r; x++) {
			float pixel = smooth_image[y][x];
			if (pixel > threshold) {
				gboolean bingo = TRUE;
				float neighbor;
				double mean = 0., meanhigh = 0.;
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
						mean += neighbor;
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
							continue;
						}
						neighbor = smooth_image[yy][xx]; 
						if (neighbor <= threshold) {
							bingo = FALSE;
							break;
						}
						if (!bingo) break;
						meanhigh += neighbor;
						count++;
					}
				}
				if (count < 8) {
					bingo = FALSE;
					continue;
				}
				mean = (mean - meanhigh) / (boxsize - 8); // (boxsize - 1) pix in mean and 8 pix in meanhigh
				/* trying to remove false positives in nebs
				mean of 9 central pixels must be above mean of the whole search box excluding them
				i.e, an approximation of the local background, by a significant amount
				*/
				// meanhigh = (meanhigh + pixel) / 9;
				// if (meanhigh - mean <= locthreshold) { 
				// 	bingo = FALSE;
				// 	continue;
				// }
				/*This check has been moved later on and done with the amplitude 
				so that the initial search box size does not influence the outcome
				Could miss weak stars when to much of the "mean" variable was polluted
				by large stars for small sampling images*/

				// first derivatives. r c is for row and column. l, r, u, d for left, right, up and down
				int xx = x - r, yy = y;
				// siril_debug_print("%d: %d - %d\n", nbstars, xx, yy);
				float d1rl, d1rr, d1cu, d1cd, r0, c0;
				d1rl = pixel - smooth_image[yy][xx - 1];
				d1rr = smooth_image[yy][xx + 1] - pixel;
				d1cu = pixel - smooth_image[yy - 1][xx];
				d1cd = smooth_image[yy + 1][xx] - pixel;
				// computing zero-crossing to guess correctly psf widths later on
				r0 = -0.5 - d1rl/(d1rr - d1rl);
				c0 = -0.5 - d1cu/(d1cd - d1cu);
				// siril_debug_print("%.2f, %.2f\n\n", r0, c0);

				float d2rr, d2rrr, d2rl, d2rll, d2cu, d2cuu, d2cd, d2cdd; // second dev
				float srr, srl, scd, scu;
				float Arr, Arl, Acd, Acu;
				int i;
				
				// Computing zero-upcrossing of the 2nd deriv row-wise moving right
				i = 0;
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
				d2rl  = smooth_image[yy][xx - i - 1] + smooth_image[yy][xx - i + 1] - 2 * smooth_image[yy][xx + i    ];
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
				d2cu =  smooth_image[yy - i - 1][xx] + smooth_image[yy - i + 1][xx] - 2 * smooth_image[yy + i    ][xx];
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
				float Sr, Ar, Br, Sc, Ac, Bc, B;
				Sr = 0.5 * (-srl + srr);
				Ar = 0.5 * (Arl + Arr);
				Br = pixel - Ar;
				Sc = 0.5 * (-scu + scd);
				Ac = 0.5 * (Acu + Acd);
				Bc = pixel - Ac;
				B = 0.5 * (Br + Bc);

				// restimate new box size to enclose enough background
				// term S / 3 increases the box radius when the guessed fwhm is larger than the smoothing kernel size
				int Rr = (int) ceil(2 * Sr * Sr / KERNEL_SIZE);
				int Rc = (int) ceil(2 * Sc * Sc / KERNEL_SIZE);
				int Rm = max(Rr, Rc);
				if (Rm > r) {
				// avoid enlarging outside frame width
					if (xx - Rm < 0)
						Rm = xx;
					if (xx + Rm >= nx)
						Rm = nx - xx - 1;
				// avoid enlarging outside frame height
					if (yy - Rm < 0)
						Rm = yy;
					if (yy + Rm >= ny)
						Rm = ny - yy - 1;
				}
				int R = max(Rm, r);

				// Quality checks
				float dA = max(Ar,Ac)/min(Ar,Ac);
				float dSr = max(-srl,srr)/min(-srl,srr);
				float dSc = max(-scu,scd)/min(-scu,scd);
				if (!com.pref.starfinder_conf.relax_checks)
					if ((dA > 2.) || (dSr > 2.) || ( dSc > 2.) || (max(Ar,Ac) < locthreshold))  bingo = FALSE;

				if (bingo && nbstars < MAX_STARS) {
					candidates[nbstars].x = xx;
					candidates[nbstars].y = yy;
					candidates[nbstars].mag_est = meanhigh;
					candidates[nbstars].bg = mean; //using local background
					candidates[nbstars].B = B;
					candidates[nbstars].R = R; // revised box size
					candidates[nbstars].sx = Sr;
					candidates[nbstars].sy = Sc;
					nbstars++;
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
	nbstars = minimize_candidates(image->fit, sf, candidates, nbstars, layer, &results, limit_nbstars, maxstars, threads);
	if (nbstars == 0)
		results = NULL;
	sort_stars_by_mag(results, nbstars);
	free(candidates);

	if (showtime) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
	if (nb_stars)
		*nb_stars = nbstars;
	return results;
}

/* returns number of stars found, result is in parameters */
static int minimize_candidates(fits *image, star_finder_params *sf, starc *candidates, int nb_candidates, int layer, psf_star ***retval, gboolean limit_nbstars, int maxstars, int threads) {
	int nx = image->rx;
	int ny = image->ry;
	WORD **image_ushort = NULL;
	float **image_float = NULL;
	gint nbstars = 0;
	sf_errors accepted_level = (com.pref.starfinder_conf.relax_checks) ? SF_RMSE_TOO_LARGE : SF_OK;

	if (image->type == DATA_USHORT) {
		image_ushort = malloc(ny * sizeof(WORD *));
		for (int k = 0; k < ny; k++)
			image_ushort[ny - k - 1] = image->pdata[layer] + k * nx;
	}
	else if (image->type == DATA_FLOAT) {
		image_float = malloc(ny * sizeof(float *));
		for (int k = 0; k < ny; k++)
			image_float[ny - k - 1] = image->fpdata[layer] + k * nx;
	}
	else return 0;

	psf_star **results = new_fitted_stars(nb_candidates);
	if (!results) {
		PRINT_ALLOC_ERR;
		return 0;
	}

	//sorting candidates by starc.mag_est values
	qsort(candidates, nb_candidates, sizeof(starc), star_cmp_by_mag_est);

	int round = 0;
	siril_debug_print("limiting stars (%d) to %d for %d candidates\n", limit_nbstars, maxstars, nb_candidates);
	int number_per_round = limit_nbstars ? maxstars + maxstars / 4 : nb_candidates;
	int lower_limit_for_this_round, upper_limit;
	do {
		lower_limit_for_this_round = round * number_per_round;
		upper_limit = (round + 1) * number_per_round;
		if (upper_limit > nb_candidates)
			upper_limit = nb_candidates;
		//siril_debug_print("round %d from %d to %d candidates\n", round, lower_limit_for_this_round, upper_limit);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
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
					for (ii = 0, i = x - R; i <= x + R;
							i++, ii++) {
						gsl_matrix_set(z, ii, jj, (double)image_ushort[j][i]);
					}
				}
			} else {
				for (jj = 0, j = y - R; j <= y + R; j++, jj++) {
					for (ii = 0, i = x - R; i <= x + R;
							i++, ii++) {
						gsl_matrix_set(z, ii, jj, (double)image_float[j][i]);
					}
				}
			}

			psf_star *cur_star = psf_global_minimisation(z, candidates[candidate].B, FALSE, FALSE, 1.0, FALSE, FALSE, NULL);
			gsl_matrix_free(z);
			if (cur_star) {
				gchar errmsg[SF_ERRMSG_LEN] = "";
				sf_errors star_invalidated = reject_star(cur_star, sf, &candidates[candidate], (DEBUG_STAR_DETECTION) ? errmsg : NULL);
				if (star_invalidated <= accepted_level) {
					//fwhm_to_arcsec_if_needed(image, cur_star);	// should we do this here?
					cur_star->layer = layer;
					cur_star->xpos = (x - R) + cur_star->x0 - 1.0;
					cur_star->ypos = (y - R) + cur_star->y0 - 1.0;
#if DEBUG_STAR_DETECTION
					if (star_invalidated > SF_OK)
						siril_debug_print("Candidate #%5d: X: %4d, Y: %4d - criterion #%2d failed (but star kept)\n%s", candidate, x, y, star_invalidated, errmsg);
#endif
					if (threads > 1)
						results[candidate] = cur_star;
					else results[nbstars++] = cur_star;
				} else {
#if DEBUG_STAR_DETECTION
					siril_debug_print("Candidate #%5d: X: %4d, Y: %4d - criterion #%2d failed\n%s", candidate, x, y, star_invalidated, errmsg);
#endif
					free_psf(cur_star);
					if (threads > 1)
						results[candidate] = NULL;
				}
			}
			else if (threads > 1)
				results[candidate] = NULL;
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

	if (retval)
		*retval = results;
	if (image_ushort) free(image_ushort);
	if (image_float) free(image_float);
	return nbstars;
}

/* Function to add star one by one, from the selection rectangle, the
 * minimization is run and the star is detected and added to the list of stars.
 *
 * IF A STAR IS FOUND and not already present in com.stars, the return value is
 * the new star and index is set to the index of the new star in com.stars.
 * IF NO NEW STAR WAS FOUND, either because it was already in the list, or a
 * star failed to be detected in the selection, or any other error, the return
 * value is NULL and index is set to -1.
 */
psf_star *add_star(fits *fit, int layer, int *index) {
	int i = 0;
	gboolean already_found = FALSE;

	*index = -1;
	psf_star *result = psf_get_minimisation(&gfit, layer, &com.selection, FALSE, FALSE, TRUE, FALSE);
	if (!result)
		return NULL;
	/* We do not check if it's matching with the "reject_star()" criteria.
	 * Indeed, in this case the user can add manually stars missed by star_finder */

	if (com.stars && !com.star_is_seqdata) {
		// check if the star was already detected/peaked
		while (com.stars[i]) {
			if (fabs(result->x0 + com.selection.x - com.stars[i]->xpos) < 0.9
					&& fabs(com.selection.y + com.selection.h - result->y0
									- com.stars[i]->ypos) < 0.9)
				already_found = TRUE;
			i++;
		}
	} else {
		if (com.star_is_seqdata) {
			/* com.stars was allocated with a size of 2, we need to free it before reallocating */
			clear_stars_list(TRUE);
		}
		com.stars = new_fitted_stars(MAX_STARS);
		if (!com.stars) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		com.star_is_seqdata = FALSE;
	}

	if (already_found) {
		free_psf(result);
		result = NULL;
		char *msg = siril_log_message(_("This star has already been picked !\n"));
		siril_message_dialog( GTK_MESSAGE_INFO, _("Peaker"), msg);
	} else {
		if (i < MAX_STARS) {
			result->xpos = result->x0 + com.selection.x - 0.5;
			result->ypos = com.selection.y + com.selection.h - result->y0 - 0.5;
			psf_star **newstars = realloc(com.stars, (i + 2) * sizeof(psf_star *));
			if (!newstars)
				PRINT_ALLOC_ERR;
			else {
				com.stars = newstars;
				com.stars[i] = result;
				com.stars[i + 1] = NULL;
				*index = i;
			}
		} else {
			free_psf(result);
			result = NULL;
		}
	}
	return result;
}

int get_size_star_tab() {
	int i = 0;
	while (com.stars[i])
		i++;
	return i;
}

/* Remove a star from com.stars, at index index. The star is freed. */
int remove_star(int index) {
	if (index < 0 || !com.stars || !com.stars[index])
		return 1;

	int N = get_size_star_tab() + 1;

	free_psf(com.stars[index]);
	memmove(&com.stars[index], &com.stars[index + 1],
			(N - index - 1) * sizeof(*com.stars));
	redraw(REDRAW_OVERLAY);
	return 0;
}

int compare_stars_by_mag(const void* star1, const void* star2) {
	psf_star *s1 = *(psf_star**) star1;
	psf_star *s2 = *(psf_star**) star2;
	if (s1->mag < s2->mag)
		return -1;
	if (s1->mag > s2->mag)
		return 1;
	return 0;
}

void sort_stars_by_mag(psf_star **stars, int total) {
	if (*(&stars))
		qsort(*(&stars), total, sizeof(psf_star*), compare_stars_by_mag);
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

void FWHM_average(psf_star **stars, int nb, float *FWHMx, float *FWHMy, char **units, float *B) {
	*FWHMx = 0.0f;
	*FWHMy = 0.0f;
	*B = 0.0f;
	if (stars && stars[0]) {
		double fwhmx = 0.0, fwhmy = 0.0, b = 0.0;
		*units = stars[0]->units;
		for (int i = 0; i < nb; i++) {
			fwhmx += stars[i]->fwhmx;
			fwhmy += stars[i]->fwhmy;
			b += stars[i]->B;
		}
		*FWHMx = (float)(fwhmx / (double)nb);
		*FWHMy = (float)(fwhmy / (double)nb);
		*B = (float)(b / (double)nb);
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

	float retval = siril_stats_trmean_from_sorted_data(0.15, fwhms, 1, nb);
	free(fwhms);
	return retval;
}

static gboolean end_findstar(gpointer p) {
	struct starfinder_data *args = (struct starfinder_data *) p;
	stop_processing_thread();
	if (com.stars)
		refresh_star_list(com.stars);

	set_cursor_waiting(FALSE);

	free(args);
	return FALSE;
}

gpointer findstar(gpointer p) {
	struct starfinder_data *args = (struct starfinder_data *)p;

	int nbstars = 0;
	rectangle *selection = NULL;
	if (com.selection.w != 0 && com.selection.h != 0)
		selection = &com.selection;
	psf_star **stars = peaker(&args->im, args->layer, &com.pref.starfinder_conf, &nbstars, selection, TRUE, FALSE, MAX_STARS_FITTED, com.max_thread);
	if (stars) {
		clear_stars_list(FALSE);
		com.stars = stars;
	}
	siril_log_message(_("Found %d stars in %s, channel #%d\n"), nbstars,
			selection ? _("selection") : _("image"), args->layer);

	siril_add_idle(end_findstar, args);

	return GINT_TO_POINTER(0);
}

void on_process_starfinder_button_clicked(GtkButton *button, gpointer user_data) {
	int layer = gui.cvport == RGB_VPORT ? GLAYER : gui.cvport;
	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		siril_log_color_message(_("Load an image first, aborted.\n"), "red");
		return;
	}
	confirm_peaker_GUI(); //making sure the spin buttons values are read even without confirmation

	struct starfinder_data *args = malloc(sizeof(struct starfinder_data));
	args->im.fit = &gfit;
	if (sequence_is_loaded() && com.seq.current >= 0) {
		args->im.from_seq = &com.seq;
		args->im.index_in_seq = com.seq.current;
	} else {
		args->im.from_seq = NULL;
		args->im.index_in_seq = -1;
	}
	args->layer = layer;

	start_in_new_thread(findstar, args);
}
