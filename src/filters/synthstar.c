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
/*
float calculate_mean_box_f(float* buf, float minfwhm, float xpos, float ypos, int dimx, int dimy) {
	float tot = 0.0f;
	float dim = (float) minfwhm + 1;
	int num = 0;
	int rsq = dim * dim;
	float insq = (float) (minfwhm + 1) / 2.7;
	for (int x=xpos-dim; x<xpos+dim; x++) {
		for (int y=ypos-dim; y<ypos+dim; y++) {
			if (x > 0 && x < dimx && y > 0 && y < dimy) {
				int distsq = ((y-ypos)*(y-ypos)) + ((x-xpos)*(x-xpos));
				if (distsq >= insq && distsq <= rsq) {
					tot += buf[x + ((dimy-y) * dimx)] * buf[x+((dimy-y)*dimx)];
					num++;
				}
			}
		}
	}
	return sqrt(tot / num);
}

float calculate_mean_box_W(WORD* buf, float minfwhm, float xpos, float ypos, int dimx, int dimy, float invnorm) {
	float tot = 0.0f;
	int dim = (int) ((minfwhm+1) / 2);
	int num = 0;
	int rsq = dim * dim;
	for (int x=xpos-dim; x<xpos+dim; x++) {
		for (int y=ypos-dim; y<ypos+dim; y++) {
			if (x > 0 && x < dimx && y > 0 && y < dimy) {
				if (((y-ypos)*(y-ypos)) + ((x-xpos)*(x-xpos)) < rsq) {
					tot += (float) buf[x + ((dimy-y) * dimx)]/invnorm * (float) buf[x+((dimy-y)*dimx)]/invnorm;
					num++;
				}
			}
		}
	}
	return sqrt(tot / num);
}
*/
void makemoffat(float* psf, int size, float fwhm, float lum, float xoff, float yoff, int type) {
#define BETA 2.2f
	const float alpha = 0.6667f * fwhm;
	const int halfpsfdim = (size - 1) / 2;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,8) collapse(2) num_threads(com.max_thread) if(com.max_thread > 1)
#endif
	for (int x=-halfpsfdim; x<=halfpsfdim; x++) {
		for (int y=-halfpsfdim; y<=halfpsfdim; y++) {
			float xf = x-xoff;
			float yf = y-yoff;
			psf[(x+halfpsfdim)+((y+halfpsfdim)*size)] = lum * pow(1.0f + ((xf*xf + yf*yf)/(alpha*alpha)),-BETA);
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
#undef BETA
}

void add_star_to_rgb_buffer(float *H, float *S, float *psfL, int size, float *Hsynth, float *Ssynth, float *Lsynth, int x, int y, int dimx, int dimy) {
	int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#define EPSILON 1e-30
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,8) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
#endif
	for (int psfx=0; psfx<size; psfx++) {
		for (int psfy=0; psfy<size; psfy++) {
			xx=x+psfx-halfpsfdim;
			yy=y+psfy-halfpsfdim;
			if (xx > 0 && xx < dimx && yy > 0 && yy < dimy) {
//				float factor = psfL[psfx+(psfy*size)] / (max(EPSILON,(Lsynth[xx+((dimy-yy)*dimx)] + psfL[psfx+(psfy*size)])));
//				Hsynth[xx+((dimy-yy)*dimx)] = ((1.0f - factor) * Hsynth[xx+((dimy-yy)*dimx)]) + (factor * H[xx+((dimy-yy)*dimx)]);
//				Ssynth[xx+((dimy-yy)*dimx)] = ((1.0f - factor) * Ssynth[xx+((dimy-yy)*dimx)]) + (factor * S[xx+((dimy-yy)*dimx)]);
				Hsynth[xx+((dimy-yy)*dimx)] = H[xx+((dimy-yy)*dimx)];
				Ssynth[xx+((dimy-yy)*dimx)] = S[xx+((dimy-yy)*dimx)];
				Lsynth[xx+((dimy-yy)*dimx)] += psfL[psfx+(psfy*size)];
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void add_star_to_mono_buffer(float *psfL, int size, float *Lsynth, int x, int y, int dimx, int dimy) {
	const int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,8) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
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

int starcount(psf_star **stars) {
	int i=0;
	if (!(stars)) {
		return 0;
	} else {
		while (stars[i]) {
			i++;
		}
	}
	return i;
}

gpointer do_synthstar() {
	generate_synthstars(&gfit);
	return GINT_TO_POINTER(0);
}

int generate_synthstars(fits *fit) {
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	siril_log_color_message(_("Star intensive care: processing...\n"), "green");
	gboolean is_RGB = TRUE;
	gboolean is_32bit = TRUE;
	float norm = 1.0f, invnorm = 1.0f;
	int nb_stars = starcount(com.stars);
	if (nb_stars < 1) {
		siril_log_color_message(_("No stars found. Use the Dynamic PSF tool to detect stars in the image.\n"), "red");
		return -1;
	} else {
		siril_log_message(_("Resynthesizing %d stars...\n"), nb_stars);
	}
	if (fit->type == DATA_USHORT) {
		is_32bit = FALSE;
		norm = get_normalized_value(fit);
		invnorm = 1.0f / norm;
	}
	if (fit->naxes[2] != 3)
		is_RGB = FALSE;

	int dimx = fit->naxes[0];
	int dimy = fit->naxes[1];
	int count = dimx * dimy;
	float *H, *S, *Hsynth = NULL, *Ssynth = NULL, *Lsynth, junk;
	Lsynth = (float *) calloc(count, sizeof(float));

	// For RGB images, convert pixel colour data from fit into H and S arrays. L is irrelevant as we will synthesize L.
	if (is_RGB) {
		H = (float *) calloc(count, sizeof(float));
		S = (float *) calloc(count, sizeof(float));
		Hsynth = (float *) calloc(count, sizeof(float));
		Ssynth = (float *) calloc(count, sizeof(float));
		for (int n = 0; n < count; n++) {
			if (is_32bit)
				rgb_to_hsl_float_sat(fit->fpdata[0][n], fit->fpdata[1][n], fit->fpdata[2][n], 0.f, &H[n], &S[n], &junk);
			else
				rgb_to_hsl_float_sat(fit->pdata[0][n] * invnorm, fit->pdata[1][n] * invnorm, fit->pdata[2][n] * invnorm, 0.f, &H[n], &S[n], &junk);
		}
	}

	// Synthesize a PSF for each star in the star array s, based on its measured parameters
	for (int n = 0; n < nb_stars; n++) {
		float lum = com.stars[n]->A;
		if (lum < 0.0f)
			lum = 0.0f;
		if (!is_32bit)
			lum *= invnorm;
		assert(lum > 0.0f);
		float xoff = com.stars[n]->xpos - (int) com.stars[n]->xpos;
		float yoff = com.stars[n]->ypos - (int) com.stars[n]->ypos;
		int size=20 * max(com.stars[n]->sx, com.stars[n]->sy); // This is big enough that even under extreme stretching the synthesized psf tails off smoothly
		if (!(size %2))
			size++;
		float minfwhm = min(com.stars[n]->fwhmx, com.stars[n]->fwhmy);
//		float psfH, psfS;

		// For RGB images, obtain colour information to ensure that synthesized star profiles are the same average hue as the stars they replace.
//		if (is_RGB) {
//			psfH = calculate_mean_box_f(H, minfwhm, com.stars[n]->xpos, com.stars[n]->ypos, dimx, dimy);
//			psfS = calculate_mean_box_f(S, minfwhm, com.stars[n]->xpos, com.stars[n]->ypos, dimx, dimy);
//		}

		// Synthesize the Moffat luminance profile and add to the star mask in HSL colourspace
		float *psfL = (float *) calloc(size * size, sizeof(float));
		makemoffat(psfL, size, minfwhm, lum, xoff, yoff, SYNTHESIZE_MOFFAT);
		if (is_RGB)
			add_star_to_rgb_buffer(H, S, psfL, size, Hsynth, Ssynth, Lsynth, com.stars[n]->xpos, com.stars[n]->ypos, dimx, dimy);
		else
			add_star_to_mono_buffer(psfL, size, Lsynth, com.stars[n]->xpos, com.stars[n]->ypos, dimx, dimy);
		free(psfL);
	}

	// Construct the RGB from synthetic L (and for RGB images, also the H and S values from the orginal image thus giving our synthesized stars the correct colour)
	if (is_RGB) {
		float *R, *G, *B;
		R = (float *) calloc(count, sizeof(float));
		G = (float *) calloc(count, sizeof(float));
		B = (float *) calloc(count, sizeof(float));
#ifdef _OPENMP
#pragma omp parallel
{
#endif
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
		for (int n = 0; n < count; n++) {
			hsl_to_rgb_float_sat(Hsynth[n], Ssynth[n], Lsynth[n], &R[n], &G[n], &B[n]);
			// Trap NaNs and infinities
			if (isnan(R[n]))
				R[n] = 0.f;
			if (isinf(R[n]))
				R[n] = 0.f;
			if (R[n] < 0.f)
				R[n] = 0.f;
			if (isnan(G[n]))
				G[n] = 0.f;
			if (isinf(G[n]))
				G[n] = 0.f;
			if (G[n] < 0.f)
				G[n] = 0.f;
			if (isnan(B[n]))
				B[n] = 0.f;
			if (isinf(B[n]))
				B[n] = 0.f;
			if (B[n] < 0.f)
				B[n] = 0.f;
		}
#ifdef _OPENMP
#pragma omp barrier
#endif
		if (is_32bit) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (int n = 0; n < count ; n++) {
				fit->fpdata[0][n] = R[n];
				fit->fpdata[1][n] = G[n];
				fit->fpdata[2][n] = B[n];
			}
		} else {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (int n = 0; n < count ; n++) {
				fit->pdata[0][n] = roundf_to_WORD(R[n] * norm);
				fit->pdata[1][n] = roundf_to_WORD(G[n] * norm);
				fit->pdata[2][n] = roundf_to_WORD(B[n] * norm);
			}
		}
#ifdef _OPENMP
}
#endif

		// Free memory
		free(R);
		free(G);
		free(B);
		free(H);
		free(S);
		free(Hsynth);
		free(Ssynth);
	} else {
		// Mono image. Populate the L values into the fits WORD or float data array
#ifdef _OPENMP
#pragma omp parallel
{
#endif
		if (is_32bit) {
#ifdef _OPENMP
#pragma omp for simd schedule(static,8)
#endif
			for (int n = 0; n < count; n++) {
				fit->fdata[n] = (float) Lsynth[n];
			}
			if (com.pref.force_16bit) {
				const size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
				fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
			}
		} else {
#ifdef _OPENMP
#pragma omp for simd schedule(static,8)
#endif
			for (int n = 0; n < count; n++) {
				fit->data[n] = roundf_to_WORD(Lsynth[n] * norm);
			}
		}
	}
#ifdef _OPENMP
}
#endif
	free(Lsynth);
	if (fit == &gfit)
		notify_gfit_modified();
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, "Execution time");
	return 0;
}

void on_synthstar_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("synthstar_dialog");
}

void on_synthstar_dynpsf_clicked(GtkButton *button, gpointer user_data) {
	siril_open_dialog("stars_list_window");
}

void on_synthstar_apply_clicked(GtkButton *button, gpointer user_data) {
	start_in_new_thread(do_synthstar, NULL);
	siril_close_dialog("synthstar_dialog");
}

