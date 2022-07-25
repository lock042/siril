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
/*
void makegaussian(double* psf, int size, double fwhm, double lum, double xoffset, double yoffset, int type) {
	int halfpsfdim = (size - 1) / 2;
	double sigma = fwhm / _2_SQRT_2_LOG2;
	double tss = 2 * sigma * sigma;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2) num_threads(com.max_thread) if(com.max_thread > 1)
#endif
	for (int x=-halfpsfdim; x<=halfpsfdim; x++) {
		for (int y=-halfpsfdim; y<=halfpsfdim; y++) {
			double xf = x-xoffset;
			double yf = y-yoffset;
			psf[(x+halfpsfdim)+((y+halfpsfdim)*size)] = lum * exp(-(((xf*xf)/tss) + ((yf*yf)/tss)));
		}
	}
}
*/

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
}

void add_star_to_rgb_buffer(double psfH, double psfS, double *psfV, int size, double *Hsynth, double *Ssynth, double *Vsynth, int x, int y, int dimx, int dimy) {
	int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) collapse(2) num_threads(com.max_thread) if(com.max_thread > 1)
#endif
	for (int psfx=0; psfx<size; psfx++) {
		for (int psfy=0; psfy<size; psfy++) {
			xx=x+psfx-halfpsfdim;
			yy=y+psfy-halfpsfdim;
			if (xx > 0 && xx < dimx && yy>0 && yy < dimy) {
				double factor = psfV[psfx+(psfy*size)] / (Vsynth[xx+((dimy-yy)*dimx)] + psfV[psfx+(psfy*size)]);
				Hsynth[xx+((dimy-yy)*dimx)] = ((1 - factor) * Hsynth[xx+((dimy-yy)*dimx)]) + (factor * psfH);
				Ssynth[xx+((dimy-yy)*dimx)] = ((1 - factor) * Ssynth[xx+((dimy-yy)*dimx)]) + (factor * psfS);
				Vsynth[xx+((dimy-yy)*dimx)] += psfV[psfx+(psfy*size)];
			}
		}
	}
}
void add_star_to_mono_buffer(double *psfV, int size, double *Vsynth, int x, int y, int dimx, int dimy) {
	int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) collapse(2) num_threads(com.max_thread) if(com.max_thread > 1)
#endif
	for (int psfx=0; psfx<size; psfx++) {
		for (int psfy=0; psfy<size; psfy++) {
			xx=x+psfx-halfpsfdim;
			yy=y+psfy-halfpsfdim;
			if (xx > 0 && xx < dimx && yy>0 && yy < dimy) {
				Vsynth[xx+((dimy-yy)*dimx)] += psfV[psfx+(psfy*size)];
			}
		}
	}
}

int generate_synthstars(fits *fit, fits *starless) {
	struct timeval t_start, t_end;
	siril_log_message(_("Star intensive care: processing...\n"));
//	gettimeofday(&t_start, NULL);
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

	fprintf(stdout, "naxes %d\n", fit->naxes[2]);
	int dimx = fit->naxes[0];
	int dimy = fit->naxes[1];
	int count = dimx * dimy;
	// Detect and replace as many stars as possible
	int orig_starfind_sigma = com.pref.starfinder_conf.sigma;
	com.pref.starfinder_conf.sigma = 0.05;
	double *H, *S, *V, *Hsynth, *Ssynth, *Vsynth, junk;
	V = (double *) calloc(count, sizeof(double));
	Vsynth = (double *) calloc(count, sizeof(double));
	// Convert pixel data from fit into H and S arrays. V is irrelevant as we will synthesize V.
	if (is_RGB) {
		H = (double *) calloc(count, sizeof(double));
		S = (double *) calloc(count, sizeof(double));
		Hsynth = (double *) calloc(count, sizeof(double));
		Ssynth = (double *) calloc(count, sizeof(double));
		for (int n = 0; n < count; n++) {
			if (is_32bit)
				rgb_to_hsv(fit->fpdata[0][n], fit->fpdata[1][n], fit->fpdata[2][n], &H[n], &S[n], &V[n]);
			else
				rgb_to_hsv(fit->pdata[0][n] * invnorm, fit->pdata[1][n] * invnorm, fit->pdata[2][n] * invnorm, &H[n], &S[n], &V[n]);
		}
	} else {
		for (int n = 0; n < count; n++) {
			if (is_32bit)
				V[n] = (double) fit->fdata[n];
			else
				V[n] = (double) fit->data[n] * invnorm;
		}
	}

	// Detect stars in the original image to replace with synthesized ones
	image im = { .fit = fit, .from_seq = NULL, .index_in_seq = -1 };
	int nb_stars = 0;
	psf_star **s = peaker(&im, 0, &com.pref.starfinder_conf, &nb_stars, NULL, TRUE, FALSE, MAX_STARS, com.max_thread);

	// Synthesize a PSF for each star in the star array s, based on its measured minimum fwhm and
	// peak brightness, and add it to the synthetic V buffer in the correct place
	siril_log_message(_("Resynthesizing %d stars...\n"), nb_stars);
	for (int n = 0; n < nb_stars; n++) {
		double lum = s[n]->A;
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
			rgb_to_hsv(psfR, psfG, psfB, &psfH, &psfS, &junk);
		}
		double *psfV = (double *) calloc(size * size, sizeof(double));

		// Set parameter type for testing... Likely to settle on Moffat profiles only
		int type = SYNTHESIZE_MOFFAT;

		switch(type) {
			case SYNTHESIZE_GAUSSIAN:
				makegaussian(psfV, size, minfwhm, lum, xoff, yoff, SYNTHESIZE_GAUSSIAN);
				break;
			case SYNTHESIZE_MOFFAT:
				makemoffat(psfV, size, minfwhm, lum, xoff, yoff, SYNTHESIZE_MOFFAT);
				break;
		}
		if (is_RGB)
			add_star_to_rgb_buffer(psfH, psfS, psfV, size, Hsynth, Ssynth, Vsynth, s[n]->xpos, s[n]->ypos, dimx, dimy);
		else
			add_star_to_mono_buffer(psfV, size, Vsynth, s[n]->xpos, s[n]->ypos, dimx, dimy);
		free(psfV);
	}

	// Construct colour stars from synthetic V (and for RGB images, also the H and S values from the orginal image thus giving our synthesized stars the correct colour)
	if (is_RGB) {
		double *R, *G, *B;
		R = (double *) calloc(count, sizeof(double));
		G = (double *) calloc(count, sizeof(double));
		B = (double *) calloc(count, sizeof(double));
		for (int n = 0; n < count; n++) {
			hsl_to_rgb(Hsynth[n], Ssynth[n], Vsynth[n], &R[n], &G[n], &B[n]);
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
			if (is_32bit)
				fit->fdata[n] = (float) Vsynth[n];
			else
				fit->data[n] = roundf_to_WORD((float) Vsynth[n] * norm);
		}
	}

	free(V);
	free(Vsynth);

	com.pref.starfinder_conf.sigma = orig_starfind_sigma;

//	gettimeofday(&t_end, NULL);
//	show_time(t_start, t_end);

	return 0;
}




