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

#include <assert.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/colors.h"
#include "algos/star_finder.h"
#include "algos/PSF.h"
#include "algos/extraction.h"
#include "algos/geometry.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "filters/mtf.h"
#include "filters/synthstar.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "opencv/opencv.h"

float calculate_mean_box_f(float* buf, double minfwhm, double xpos, double ypos, int dimx, int dimy) {
	double tot = 0.0;
	int dim = (int) ((minfwhm+1) / 2);
	int num = 0;
	int rsq = dim * dim;
	for (int x=xpos-dim; x<xpos+dim; x++) {
		for (int y=ypos-dim; y<ypos+dim; y++) {
			if (x > 0 && x < dimx && y > 0 && y < dimy) {
				if (((y-ypos)*(y-ypos)) + ((x-xpos)*(x-xpos)) < rsq) {
					tot += buf[x + ((dimy-y) * dimx)] * buf[x+((dimy-y)*dimx)];
					num++;
				}
			}
		}
	}
	return sqrt(tot / num);
}

double calculate_mean_box_W(WORD* buf, double minfwhm, double xpos, double ypos, int dimx, int dimy, double invnorm) {
	double tot = 0.0;
	int dim = (int) ((minfwhm+1) / 2);
	int num = 0;
	int rsq = dim * dim;
	for (int x=xpos-dim; x<xpos+dim; x++) {
		for (int y=ypos-dim; y<ypos+dim; y++) {
			if (x > 0 && x < dimx && y > 0 && y < dimy) {
				if (((y-ypos)*(y-ypos)) + ((x-xpos)*(x-xpos)) < rsq) {
					tot += (double) buf[x + ((dimy-y) * dimx)]/invnorm * (double) buf[x+((dimy-y)*dimx)]/invnorm;
					num++;
				}
			}
		}
	}
	return sqrt(tot / num);
}

void makemoffat(double* psf, int size, double fwhm, double lum, double xoff, double yoff, int type) {
	double beta = 2.2;
	double alpha = 0.6667 * fwhm;
	int halfpsfdim = (size - 1) / 2;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) collapse(2) num_threads(com.max_thread) if(com.max_thread > 1)
#endif
	for (int x=-halfpsfdim; x<=halfpsfdim; x++) {
		for (int y=-halfpsfdim; y<=halfpsfdim; y++) {
			double xf = x-xoff;
			double yf = y-yoff;
			psf[(x+halfpsfdim)+((y+halfpsfdim)*size)] = lum * pow(1 + ((xf*xf + yf*yf)/(alpha*alpha)),-beta);
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void add_star_to_rgb_buffer(double psfH, double psfS, double *psfL, int size, double *Hsynth, double *Ssynth, double *Lsynth, int x, int y, int dimx, int dimy) {
	int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#define EPSILON 1e-30
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
#endif
	for (int psfx=0; psfx<size; psfx++) {
		for (int psfy=0; psfy<size; psfy++) {
			xx=x+psfx-halfpsfdim;
			yy=y+psfy-halfpsfdim;
			if (xx > 0 && xx < dimx && yy > 0 && yy < dimy) {
				double factor = psfL[psfx+(psfy*size)] / (max(EPSILON,(Lsynth[xx+((dimy-yy)*dimx)] + psfL[psfx+(psfy*size)])));
				Hsynth[xx+((dimy-yy)*dimx)] = ((1 - factor) * Hsynth[xx+((dimy-yy)*dimx)]) + (factor * psfH);
				Ssynth[xx+((dimy-yy)*dimx)] = ((1 - factor) * Ssynth[xx+((dimy-yy)*dimx)]) + (factor * psfS);
				Lsynth[xx+((dimy-yy)*dimx)] += psfL[psfx+(psfy*size)];
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void add_star_to_mono_buffer(double *psfL, int size, double *Lsynth, int x, int y, int dimx, int dimy) {
	int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
#endif
	for (int psfx=0; psfx<size; psfx++) {
		for (int psfy=0; psfy<size; psfy++) {
			xx=x+psfx-halfpsfdim;
			yy=y+psfy-halfpsfdim;
			if (xx > 0 && xx < dimx && yy>0 && yy < dimy) {
				Lsynth[xx+((dimy-yy)*dimx)] += psfL[psfx+(psfy*size)];
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif

	return;
}

int generate_synthstars(fits *fit, fits *starless) {
	siril_log_message(_("Star intensive care: processing...\n"));
	gboolean is_RGB = TRUE;
	gboolean is_32bit = TRUE;
	double norm, invnorm;
	if (fit->type == DATA_USHORT) {
		is_32bit = FALSE;
		norm = get_normalized_value(fit);
		invnorm = 1.0 / norm;
	}
	if (fit->naxes[2] != 3)
		is_RGB = FALSE;

	int dimx = fit->naxes[0];
	int dimy = fit->naxes[1];
	int count = dimx * dimy;
	// Detect and replace as many stars as possible
	int orig_starfind_sigma = com.pref.starfinder_conf.sigma;
	com.pref.starfinder_conf.sigma = 0.05;
	double *H, *S, *Hsynth, *Ssynth, *Lsynth, junk;
	Lsynth = (double *) calloc(count, sizeof(double));

	// For RGB images, convert pixel colour data from fit into H and S arrays. L is irrelevant as we will synthesize L.
	if (is_RGB) {
		H = (double *) calloc(count, sizeof(double));
		S = (double *) calloc(count, sizeof(double));
		Hsynth = (double *) calloc(count, sizeof(double));
		Ssynth = (double *) calloc(count, sizeof(double));
		for (int n = 0; n < count; n++) {
			if (is_32bit)
				rgb_to_hsl(fit->fpdata[0][n], fit->fpdata[1][n], fit->fpdata[2][n], &H[n], &S[n], &junk);
			else
				rgb_to_hsl(fit->pdata[0][n] * invnorm, fit->pdata[1][n] * invnorm, fit->pdata[2][n] * invnorm, &H[n], &S[n], &junk);
		}
	}

	// Detect stars in the original image to replace with synthesized ones
	image im = { .fit = fit, .from_seq = NULL, .index_in_seq = -1 };
	int nb_stars = 0;
	psf_star **s = peaker(&im, 0, &com.pref.starfinder_conf, &nb_stars, NULL, TRUE, FALSE, MAX_STARS, com.max_thread);

	// Synthesize a PSF for each star in the star array s, based on its measured parameters
	siril_log_message(_("Resynthesizing %d stars...\n"), nb_stars);
	for (int n = 0; n < nb_stars; n++) {
		double lum = s[n]->A;
		if (lum < 0.0)
			lum = 0.0;
		if (!is_32bit)
			lum *= invnorm;
		assert(lum > 0.0);
		double xoff = s[n]->xpos - (int) s[n]->xpos;
		double yoff = s[n]->ypos - (int) s[n]->ypos;
		int size=14 * max(s[n]->sx, s[n]->sy); // This is big enough that even under extreme stretching the synthesized psf tails off smoothly
		if (!(size %2))
			size++;
		double minfwhm = min(s[n]->fwhmx, s[n]->fwhmy);
		double psfH, psfS;

		// For RGB images, obtain colour information to ensure that synthesized star profiles are the same average hue as the stars they replace.
		if (is_RGB) {
			double psfR, psfG, psfB;
			if (is_32bit) {
				psfR = (double) calculate_mean_box_f(fit->fpdata[0], minfwhm, s[n]->xpos, s[n]->ypos, dimx, dimy);
				psfG = (double) calculate_mean_box_f(fit->fpdata[1], minfwhm, s[n]->xpos, s[n]->ypos, dimx, dimy);
				psfB = (double) calculate_mean_box_f(fit->fpdata[2], minfwhm, s[n]->xpos, s[n]->ypos, dimx, dimy);
			} else {
				psfR = calculate_mean_box_W(fit->pdata[0], minfwhm, s[n]->xpos, s[n]->ypos, dimx, dimy, invnorm);
				psfG = calculate_mean_box_W(fit->pdata[1], minfwhm, s[n]->xpos, s[n]->ypos, dimx, dimy, invnorm);
				psfB = calculate_mean_box_W(fit->pdata[2], minfwhm, s[n]->xpos, s[n]->ypos, dimx, dimy, invnorm);
			}
			rgb_to_hsl(psfR, psfG, psfB, &psfH, &psfS, &junk);
		}

		// Synthesize the Moffat luminance profile and add to the star mask in HSL colourspace
		double *psfL = (double *) calloc(size * size, sizeof(double));
		makemoffat(psfL, size, minfwhm, lum, xoff, yoff, SYNTHESIZE_MOFFAT);
		if (is_RGB)
			add_star_to_rgb_buffer(psfH, psfS, psfL, size, Hsynth, Ssynth, Lsynth, s[n]->xpos, s[n]->ypos, dimx, dimy);
		else
			add_star_to_mono_buffer(psfL, size, Lsynth, s[n]->xpos, s[n]->ypos, dimx, dimy);
		free(psfL);
	}

	// Construct the RGB from synthetic L (and for RGB images, also the H and S values from the orginal image thus giving our synthesized stars the correct colour)
	if (is_RGB) {
		double *R, *G, *B;
		R = (double *) calloc(count, sizeof(double));
		G = (double *) calloc(count, sizeof(double));
		B = (double *) calloc(count, sizeof(double));
		for (int n = 0; n < count; n++) {
			hsl_to_rgb(Hsynth[n], Ssynth[n], Lsynth[n], &R[n], &G[n], &B[n]);
			if (is_32bit) {
				fit->fpdata[0][n] = (float) R[n];
				fit->fpdata[1][n] = (float) G[n];
				fit->fpdata[2][n] = (float) B[n];
			} else {
				fit->pdata[0][n] = roundf_to_WORD((float) R[n] * norm);
				fit->pdata[1][n] = roundf_to_WORD((float) G[n] * norm);
				fit->pdata[2][n] = roundf_to_WORD((float) B[n] * norm);
			}
		}
		// Free memory
		free(R);
		free(G);
		free(B);
		free(H);
		free(S);
		free(Hsynth);
		free(Ssynth);
	} else {
		for (int n = 0; n < count; n++) {
			if (is_32bit) {
				fit->fdata[n] = (float) Lsynth[n];
				if (com.pref.force_16bit) {
					const size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
					fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
				}
			}
			else
				fit->data[n] = roundf_to_WORD((float) Lsynth[n] * norm);
		}
	}
	free(Lsynth);
	com.pref.starfinder_conf.sigma = orig_starfind_sigma;
	return 0;
}




