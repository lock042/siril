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
#include "core/siril_log.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/colors.h"
#include "algos/median_fast.h"
#include "algos/star_finder.h"
#include "algos/PSF.h"
#include "algos/extraction.h"
#include "algos/geometry.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "algos/star_finder.h"
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

int generate_synthstars(fits *fit);
int reprofile_saturated_stars(fits *fit);

void makemoffat(float *psf, const int size, const float fwhm, const float lum, const float xoff, const float yoff, const float beta) {
	const float alpha = 0.6667f * fwhm;
	const int halfpsfdim = (size - 1) / 2;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,8) collapse(2) num_threads(com.max_thread) if(com.max_thread > 1)
#endif
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float xf = x - xoff + 0.5f;
			float yf = y - yoff - 0.5f;
			psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = lum
					* powf(1.0f + ((xf * xf + yf * yf) / (alpha * alpha)),
							-beta);
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
#undef BETA
}

void makegaussian(float *psf, int size, float fwhm, float lum, float xoffset, float yoffset) {
	int halfpsfdim = (size - 1) / 2;
	float sigma = fwhm / _2_SQRT_2_LOG2;
	float tss = 2 * sigma * sigma;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,8) collapse(2) num_threads(com.max_thread) if(com.max_thread > 1)
#endif
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float xf = x - xoffset + 0.5f;
			float yf = y - yoffset - 0.5f;
			psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = lum
					* expf(-(((xf * xf) / tss) + ((yf * yf) / tss)));
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void add_star_to_rgb_buffer(float *H, float *S, float *psfL, int size, float *Hsynth, float *Ssynth, float *Lsynth, int x, int y, int dimx, int dimy) {
	int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#define EPSILON 1e-30
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,8) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
#endif
	for (int psfx = 0; psfx < size; psfx++) {
		for (int psfy = 0; psfy < size; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
			if (xx > 0 && xx < dimx && yy > 0 && yy < dimy) {
				Hsynth[xx + ((dimy - yy) * dimx)] =
						H[xx + ((dimy - yy) * dimx)];
				Ssynth[xx + ((dimy - yy) * dimx)] =
						S[xx + ((dimy - yy) * dimx)];
				Lsynth[xx + ((dimy - yy) * dimx)] += psfL[psfx + (psfy * size)];
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
	for (int psfx = 0; psfx < size; psfx++) {
		for (int psfy = 0; psfy < size; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
			if (xx > 0 && xx < dimx && yy > 0 && yy < dimy) {
				Lsynth[xx + ((dimy - yy) * dimx)] += psfL[psfx + (psfy * size)];
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif

	return;
}

void replace_sat_star_in_buffer(float *psfL, int size, float *Lsynth, int x, int y, int dimx, int dimy, float sat, float bg, float noise) {
	const int halfpsfdim = (size - 1) / 2;
	int xx, yy;

#ifdef _OPENMP
#pragma omp parallel
{
omp_set_num_threads(com.max_thread);
#pragma omp for simd schedule(static) collapse(2) private(xx, yy)
#endif
	for (int psfx = 0; psfx < size; psfx++) {
		for (int psfy = 0; psfy < size; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
			if (xx >= 0 && xx < dimx && yy >= 0 && yy < dimy) {
				float orig = Lsynth[xx + ((dimy - yy) * dimx)];
				float synth = psfL[psfx + (psfy * size)];
				float synthfactor = (orig < sat) ? 1.f -
					((sat - orig) / sat) : 1.f;
				Lsynth[xx + ((dimy - yy) * dimx)] += max(synth * synthfactor
					+ orig * (1 - synthfactor) - sat, 0.f);
			}
		}
	}
#ifdef _OPENMP
#pragma omp for simd schedule(static,8) collapse(2) private(xx, yy)
#endif
// 3x3 median blur to smooth off the join between synthetic PSF and original data
	for (int psfx = 0; psfx < size; psfx++) {
		for (int psfy = 0; psfy < size; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
			if (xx > 0 && xx < dimx && yy > 0 && yy < dimy)
				if ((psfx - halfpsfdim) * (psfx - halfpsfdim) + (psfy - halfpsfdim) * (psfy - halfpsfdim) < (size / 6) * (size / 6))
					Lsynth[xx + ((dimy - yy) * dimx)] = median9f(
						Lsynth[xx - 1 + ((dimy - yy - 1) * dimx)],
						Lsynth[xx + ((dimy - yy - 1) * dimx)],
						Lsynth[xx + 1 + ((dimy - yy - 1) * dimx)],
						Lsynth[xx - 1 + ((dimy - yy) * dimx)],
						Lsynth[xx + ((dimy - yy) * dimx)],
						Lsynth[xx + 1 + ((dimy - yy) * dimx)],
						Lsynth[xx - 1 + ((dimy - yy + 1) * dimx)],
						Lsynth[xx + ((dimy - yy + 1) * dimx)],
						Lsynth[xx + 1 + ((dimy - yy + 1) * dimx)]);
		}
	}
#ifdef _OPENMP
#pragma omp barrier
}
#endif
	return;
}

int starcount(psf_star **stars) {
	int i = 0;
	if (!(stars)) {
		return 0;
	} else {
		while (stars[i]) {
			i++;
		}
	}
	return i;
}

gpointer fix_saturated_stars() {
	reprofile_saturated_stars(&gfit);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(0);
}

gpointer do_synthstar() {
	generate_synthstars(&gfit);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(0);
}

int generate_synthstars(fits *fit) {
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	siril_log_color_message(_("Star synthesis (full star mask creation): processing...\n"), "green");
	gboolean is_RGB = TRUE;
	gboolean is_32bit = TRUE;
	gboolean stars_needs_freeing = FALSE;
	float norm = 1.0f, invnorm = 1.0f;
	int nb_stars = starcount(com.stars);
	int channel = 1;
	psf_star **stars = NULL;
	if (nb_stars < 1) {
		image *input_image = NULL;
		input_image = calloc(1, sizeof(image));
		input_image->fit = fit;
		input_image->from_seq = NULL;
		input_image->index_in_seq = -1;
		if (fit->naxes[2] == 1)
			channel = 0;
		stars = peaker(input_image, channel, &com.pref.starfinder_conf, &nb_stars,
				NULL, FALSE, FALSE, MAX_STARS, MULTI_THREADED);
		free(input_image);
		stars_needs_freeing = TRUE;
	} else {
		stars = com.stars;
	}
	if (nb_stars < 1) {
		siril_log_color_message(_("No stars detected in the image.\n"), "red");
		return -1;
	} else {
		siril_log_message(_("Synthesizing %d stars...\n"), nb_stars);
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
	gboolean buf_needs_freeing = FALSE;
	// Regardless of 16/32bit store the data in a buffer, converting if needed
	float *buf[3];
	if (is_RGB) {
		if (is_32bit) {
			buf[0] = fit->fpdata[0];
			buf[1] = fit->fpdata[1];
			buf[2] = fit->fpdata[2];
		} else {
			buf[0] = (float*) calloc(count, sizeof(float));
			buf[1] = (float*) calloc(count, sizeof(float));
			buf[2] = (float*) calloc(count, sizeof(float));
			buf_needs_freeing = TRUE;
			for (size_t i = 0; i < count; i++) {
				buf[0][i] = (float) fit->pdata[0][i] * invnorm;
				buf[1][i] = (float) fit->pdata[1][i] * invnorm;
				buf[2][i] = (float) fit->pdata[2][i] * invnorm;
			}
		}
	} else { // mono
		if (is_32bit)
			buf[0] = fit->fdata;
		else {
			buf[0] = (float*) calloc(count, sizeof(float));
			buf_needs_freeing = TRUE;
			for (size_t i = 0; i < count; i++)
				buf[0][i] = (float) fit->data[i] * invnorm;
		}
	}

	// Normalize the buffer to avoid issues with colorspace conversion
	float bufmax = 0.f;
	for (size_t chan = 0; chan < fit->naxes[2]; chan++)
		for (size_t i = 0; i < count; i++)
			if (buf[chan][i] > bufmax)
				bufmax = buf[chan][i];
	for (size_t chan = 0; chan < fit->naxes[2]; chan++)
		for (size_t i = 0; i < count; i++)
			buf[chan][i] /= bufmax;

	float *H, *S, *Hsynth = NULL, *Ssynth = NULL, *Lsynth, junk;
	Lsynth = (float*) calloc(count, sizeof(float));

	// For RGB images, convert pixel colour data from fit into H and S arrays. L is irrelevant as we will synthesize L.
	if (is_RGB) {
		H = (float*) calloc(count, sizeof(float));
		S = (float*) calloc(count, sizeof(float));
		Hsynth = (float*) calloc(count, sizeof(float));
		Ssynth = (float*) calloc(count, sizeof(float));
		for (size_t i = 0; i < count; i++) {
			rgb_to_hsl_float_sat(buf[0][i], buf[1][i],
					buf[2][i], 0.f, &H[i], &S[i], &junk);
		}
	}

	// Synthesize a PSF for each star in the star array s, based on its measured parameters
	for (int n = 0; n < nb_stars; n++) {
		// Check if stop has been pressed
		if (!get_thread_run())
			return -1;
		float lum = (float) stars[n]->A;
		if (lum < 0.0f)
			lum = 0.0f;
		if (!is_32bit)
			lum *= invnorm;
		assert(lum >= 0.0f);
		float xoff = (float) stars[n]->xpos - (int) stars[n]->xpos;
		float yoff = (float) stars[n]->ypos - (int) stars[n]->ypos;
		int size = (int) 20 * max(stars[n]->fwhmx, stars[n]->fwhmy); // This is big enough that even under extreme stretching the synthesized psf tails off smoothly
		if (!(size % 2))
			size++;
		float minfwhm = min(stars[n]->fwhmx, stars[n]->fwhmy);

		// Synthesize the luminance profile and add to the star mask in HSL colourspace
		float *psfL = (float*) calloc(size * size, sizeof(float));
		float beta = 2.2f;
		if (stars[n]->has_saturated)
			makegaussian(psfL, size, minfwhm, lum, xoff, yoff);
		else
			makemoffat(psfL, size, minfwhm, lum, xoff, yoff, beta);
		if (is_RGB)
			add_star_to_rgb_buffer(H, S, psfL, size, Hsynth, Ssynth, Lsynth,
					(float) stars[n]->xpos, (float) stars[n]->ypos, dimx, dimy);
		else
			add_star_to_mono_buffer(psfL, size, Lsynth, (float) stars[n]->xpos,
					(float) stars[n]->ypos, dimx, dimy);
		free(psfL);
	}
	if (stars_needs_freeing)
		free(stars);
	// Construct the RGB from synthetic L (and for RGB images, also the H and S values from the orginal image thus giving our synthesized stars the correct colour)
	if (is_RGB) {
		float *R, *G, *B;
		R = (float*) calloc(count, sizeof(float));
		G = (float*) calloc(count, sizeof(float));
		B = (float*) calloc(count, sizeof(float));
#ifdef _OPENMP
#pragma omp parallel
		{
#endif
			float bufmax = 0.f;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (size_t chan = 0; chan < fit->naxes[2]; chan++)
				for (size_t i = 0; i < count; i++)
					if (buf[chan][i] > bufmax)
						bufmax = buf[chan][i];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (size_t chan = 0; chan < fit->naxes[2]; chan++)
				for (size_t i = 0; i < count; i++)
					buf[chan][i] /= bufmax;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (size_t n = 0; n < count; n++) {
				hsl_to_rgb_float_sat(Hsynth[n], Ssynth[n], Lsynth[n], &R[n],
						&G[n], &B[n]);
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
				for (size_t n = 0; n < count; n++) {
					fit->fpdata[0][n] = R[n];
					fit->fpdata[1][n] = G[n];
					fit->fpdata[2][n] = B[n];
				}
			} else {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
				for (size_t n = 0; n < count; n++) {
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
				for (size_t n = 0; n < count; n++) {
					fit->fdata[n] = (float) Lsynth[n];
				}
				if (com.pref.force_16bit) {
					const size_t ndata = fit->naxes[0] * fit->naxes[1]
							* fit->naxes[2];
					fit_replace_buffer(fit,
							float_buffer_to_ushort(fit->fdata, ndata),
							DATA_USHORT);
				}
			} else {
#ifdef _OPENMP
#pragma omp for simd schedule(static,8)
#endif
				for (size_t n = 0; n < count; n++) {
					fit->data[n] = roundf_to_WORD(Lsynth[n] * norm);
				}
			}
		}
#ifdef _OPENMP
	}
#endif
	free(Lsynth);
	if (buf_needs_freeing) {
		if (is_RGB) {
			for (size_t i = 0; i <3; i++)
				free(buf[i]);
		} else
			free(buf[0]);
	}
	if (fit == &gfit)
		notify_gfit_modified();
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, "Execution time");
	return 0;
}

// Fix up saturated stars only

int reprofile_saturated_stars(fits *fit) {
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	siril_log_color_message(_("Star synthesis (desaturating clipped star profiles): processing...\n"), "green");
	gboolean is_RGB = (fit->naxes[2] == 3) ? TRUE : FALSE;
	gboolean is_32bit = TRUE;
	gboolean buf_needs_freeing = FALSE;
	float norm = 1.0f, invnorm = 1.0f;
	if (fit->type == DATA_USHORT) {
		is_32bit = FALSE;
		norm = (float) get_normalized_value(fit);
		invnorm = 1.0f / norm;
	}
	siril_debug_print("norm %f, invnorm %f\n", (float) norm, (float) invnorm);
	int dimx = fit->naxes[0];
	int dimy = fit->naxes[1];
	int count = dimx * dimy;
	float *buf[3];

	if (is_RGB) {
		if (is_32bit) {
			buf[0] = fit->fpdata[0];
			buf[1] = fit->fpdata[1];
			buf[2] = fit->fpdata[2];
		} else {
			buf[0] = (float*) calloc(count, sizeof(float));
			buf[1] = (float*) calloc(count, sizeof(float));
			buf[2] = (float*) calloc(count, sizeof(float));
			buf_needs_freeing = TRUE;
			for (size_t i = 0; i < count; i++) {
				buf[0][i] = (float) fit->pdata[0][i] * invnorm;
				buf[1][i] = (float) fit->pdata[1][i] * invnorm;
				buf[2][i] = (float) fit->pdata[2][i] * invnorm;
			}
		}
	} else { // mono
		if (is_32bit)
			buf[0] = fit->fdata;
		else {
			buf[0] = (float*) calloc(count, sizeof(float));
			buf_needs_freeing = TRUE;
			for (size_t i = 0; i < count; i++)
				buf[0][i] = (float) fit->data[i] * invnorm;
		}
	}

	// Synthesize a PSF for each saturated star in the star array, based on its measured parameters. To fix saturated star profiles we have to do this for each color channel as we can't rely on the hue and saturation within the saturated area, whereas the profiles will be accurate.
	image *input_image = NULL;
	for (size_t chan = 0; chan < fit->naxes[2]; chan++) {
		input_image = calloc(1, sizeof(image));
		input_image->fit = fit;
		input_image->from_seq = NULL;
		input_image->index_in_seq = -1;
		int nb_stars;
		psf_star **stars = peaker(input_image, chan, &com.pref.starfinder_conf, &nb_stars, NULL, FALSE, FALSE, MAX_STARS, MULTI_THREADED);
		free(input_image);
		siril_log_message(_("Star synthesis: desaturating channel %u...\n"),
				chan);
		for (size_t n = 0; n < nb_stars; n++) {
			// Check if stop has been pressed
			if (!get_thread_run())
				return -1;
			if (stars[n]->has_saturated) {
				float lum = (float) stars[n]->A;
				float bg = (float) stars[n]->B;
				float sat = (float) stars[n]->sat;
				if (lum < 0.0f)
					lum = 0.0f;
				if (!is_32bit) {
					lum *= invnorm;
					bg *= invnorm;
					sat *= invnorm;
				}
				assert(lum >= 0.0f);
				float xoff = (float) stars[n]->xpos - (int) stars[n]->xpos;
				float yoff = (float) stars[n]->ypos - (int) stars[n]->ypos;
				int size = 5.f * max(stars[n]->fwhmx, stars[n]->fwhmy); // This is big enough that it should cover the saturated parts of the star
				if (!(size % 2))
					size++;
				float maxfwhm = (float) max(stars[n]->fwhmx, stars[n]->fwhmy);

				float *psfL = (float*) calloc(size * size, sizeof(float));
//				makegaussian(psfL, size, maxfwhm, (lum - bg), xoff, yoff);
				float beta = 1.5f;
				makemoffat(psfL, size, maxfwhm, (lum - bg), xoff, yoff, beta);

				// Replace the part of the profile above the sat threshold
				replace_sat_star_in_buffer(psfL, size, buf[chan],
						(float) stars[n]->xpos, (float) stars[n]->ypos, dimx,
						dimy, sat, bg, 0.f);
				free(psfL);
			}
		}
		free(stars);
	}

	// Desaturating stars will take their peak brightness over 1.f so we need to rescale the values of all pixels by a factor of (1 / maxbuf) where maxbuf is the maximum subpixel value across all channels
	siril_log_message(_("Remapping output to floating point range 0.0 to 1.0\n"));
	float bufmax = 0.f;
	for (size_t chan = 0; chan < fit->naxes[2]; chan++)
		for (size_t i = 0; i < count; i++)
			if (buf[chan][i] > bufmax)
				bufmax = buf[chan][i];
	for (size_t chan = 0; chan < fit->naxes[2]; chan++)
		for (size_t i = 0; i < count; i++)
			buf[chan][i] /= bufmax;

	if (!is_32bit) {
		for (size_t n = 0; n < count; n++) {
			if (is_RGB) {
				fit->pdata[0][n] = roundf_to_WORD(buf[0][n] * norm);
				fit->pdata[1][n] = roundf_to_WORD(buf[0][n] * norm);
				fit->pdata[2][n] = roundf_to_WORD(buf[0][n] * norm);
			}
			else
				fit->data[n] = roundf_to_WORD(buf[0][n] * norm);
		}
	}
	if (buf_needs_freeing) {
		if (is_RGB) {
			for (size_t i = 0; i <3; i++)
				free(buf[i]);
		} else
			free(buf[0]);
	}

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

void on_synthstar_desaturate_clicked(GtkButton *button, gpointer user_data) {
	undo_save_state(&gfit, "Synthetic stars: desaturate clipped stars");
	control_window_switch_to_tab(OUTPUT_LOGS);
	start_in_new_thread(fix_saturated_stars, NULL);
	siril_close_dialog("synthstar_dialog");
}

void on_synthstar_apply_clicked(GtkButton *button, gpointer user_data) {
	undo_save_state(&gfit, "Synthetic stars: full replacement");
	control_window_switch_to_tab(OUTPUT_LOGS);
	start_in_new_thread(do_synthstar, NULL);
	siril_close_dialog("synthstar_dialog");
}

