/*
 * This file is part of Siril, an astronomy image processor.
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

#include <assert.h>
#include <math.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/colors.h"
#include "algos/median_fast.h"
#include "algos/star_finder.h"
#include "algos/PSF.h"
#include "algos/extraction.h"
#include "algos/siril_random.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "filters/synthstar.h"
#include "gui/progress_and_log.h"
#include "opencv/opencv.h"

int generate_synthstars(fits *fit);
int reprofile_saturated_stars(fits *fit);

void makeairy(float *psf, const int size, const float lum, const float xoff, const float yoff, float wavelength, float aperture, float focal_length, float pixel_size, float obstruction) {
	// wavelength is given in nm; pixel size in microns; aperture and focal length are given in mm. Convert all to metres for consistency
	// obstruction is central obstruction ratio given between [0, 1[
	wavelength *= 1.e-9f;
	aperture *= 1.e-3f;
	focal_length *= 1.e-3f;
	pixel_size *= 1.e-6f;
	const int halfpsfdim = (size - 1) / 2;
	// obstruction = 0.5; // for testing purposes - can be removed
	float obscorr = (obstruction > 0.f) ? 1.f / pow(1 - obstruction * obstruction, 2.f) : 1.f;

	// Following the formulae at the Wikipedia "Airy disk" article
	const float constant = (2.f * M_PI * (aperture / 2.f) / wavelength) * (1.f / focal_length);
#ifdef _OPENMP
#pragma omp simd
#endif
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float xf = (x - xoff + 0.5f) * pixel_size;
			float yf = (y - yoff - 0.5f) * pixel_size;
			float q = constant * sqrtf(xf * xf + yf * yf);
			float bessel = j1(q);
			if (obstruction == 0.f) {
				psf[(x + halfpsfdim) + (y + halfpsfdim) * size] = (q != 0.f) ? lum * pow(2.f * bessel / q, 2.f) : lum;
			} else {
				psf[(x + halfpsfdim) + (y + halfpsfdim) * size] = (q != 0.f) ? lum * obscorr * pow(2.f / q * (bessel - obstruction * j1(obstruction * q)), 2.f) : lum;
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void makemoffat(float *psf, const int size, const float fwhm, const float lum, const float xoff,
				const float yoff, const float beta, const float ratio, const float angle) {
	float anglerad = angle * M_PI / 180.f;
	const float alpha = 0.6667f * fwhm;
	const float alphax = alpha;
	const float alphay = alpha / ratio;
	const int halfpsfdim = (size - 1) / 2;
	float a = powf(cosf(anglerad)/alphax, 2.f) + powf(sinf(anglerad)/alphay, 2.f);
	float b = powf(sinf(anglerad)/alphax, 2.f) + powf(cosf(anglerad)/alphay, 2.f);
	float c = 2.f * sinf(anglerad) * cosf(anglerad) * (1.f/(alphax * alphax) - 1.f/(alphay * alphay));
#ifdef _OPENMP
#pragma omp simd
#endif
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float xf = (x - xoff + 0.5f);
			float yf = (y - yoff - 0.5f);
			psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = lum
					* powf(1.0f + (a * xf * xf) + (b * yf * yf) + (c * xf * yf),
							-beta);
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void makegaussian(float *psf, int size, float fwhm, float lum, float xoffset, float yoffset, float ratio, float angle) {
	int halfpsfdim = (size - 1) / 2;
	float anglerad = angle * M_PI / 180.f;
	float sigmax = fwhm / _2_SQRT_2_LOG2;
	float sigmay = fwhm / (ratio * _2_SQRT_2_LOG2);
	float tssx = 2 * sigmax * sigmax;
	float tssy = 2 * sigmay * sigmay;
	float a = powf(cosf(anglerad), 2.f) / tssx + powf(sinf(anglerad), 2.f) / tssy;
	float b = sinf(2 * anglerad) / (2 * tssx) - sinf(2 * anglerad) / (2 * tssy);
	float c = powf(sinf(anglerad), 2.f) / tssx + powf(cosf(anglerad), 2.f) / tssy;
#ifdef _OPENMP
#pragma omp simd
#endif
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float xf = (x - xoffset + 0.5f);
			float yf = (y - yoffset - 0.5f);
			psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = lum
					* expf(-(a * xf * xf + 2 * b * xf * yf + c * yf * yf));
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void makedisc(float *psf, int size, float width, float lum, float xoffset, float yoffset) {
	int halfpsfdim = (size - 1) / 2;
	float radius = width / 2.f;
	// maxranditer big enough to get a good random survey of each pixel,
	// not enough to cause slowness with large kernels.
	const int maxranditer = 10000;
	float radiussq = radius * radius;
	float solidradsq = powf(radius - 0.5f, 2.f); // radius less than which pixel value is 1.f
	float zeroradsq =powf(radius + 0.5f, 2.f); // radius greater than which pixel value is 0.f
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float pixradsq = powf((x - xoffset + 0.5f), 2.f) + powf((y - yoffset - 0.5f), 2.f);
			if (pixradsq < solidradsq) {
				psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = 1.f;
			} else if (pixradsq > zeroradsq) {
				psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = 0.f;
			} else {
				int count = 0;
#ifdef _OPENMP
#pragma omp simd
#endif
				for (int randiter = 0 ; randiter < maxranditer; randiter++) {
					float xrandoff = siril_random_float();
					float yrandoff = siril_random_float();
					float xf = x - xoffset + xrandoff + 0.5f;
					float yf = y - yoffset + yrandoff - 0.5f;
					if ((xf * xf + yf * yf) < radiussq)
						count++;
				}
				psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = (lum * count) / maxranditer;
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

static void add_star_to_rgb_buffer(const float *H, const float *S, const float *psfL, int size, float *Hsynth, float *Ssynth, float *Lsynth, int x, int y, int dimx, int dimy) {
	int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#define EPSILON 1e-30
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
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

static void add_star_to_mono_buffer(const float *psfL, int size, float *Lsynth, int x, int y, int dimx, int dimy) {
	const int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
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

static void replace_sat_star_in_buffer(const float *psfL, int size, float *Lsynth, int x, int y, int dimx, int dimy, float sat, float bg, float noise) {
	float* buf = calloc(1, size * size * sizeof(float));
	float* resbuf = malloc(size * size * sizeof(float));
	const int halfpsfdim = (size - 1) / 2;
	int xx, yy;
// Blend synthetic data into Lsynth and make a copy in buf for filtering the join
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
{
#pragma omp for simd schedule(static) collapse(2) private(xx, yy)
#endif
	for (int psfx = 0; psfx < size; psfx++) {
		for (int psfy = 0; psfy < size; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
// Note the bounds, xx >= 0 and < dimx but yy > 0 and <= dimy
// This is correct, it is because Lsynth gets indexed by xx but by dimy - yy
// Same comment applies below in the copy back
			if (xx >= 0 && xx < dimx && yy > 0 && yy <= dimy) {
				float orig = Lsynth[xx + ((dimy - yy) * dimx)];
				float synth = psfL[psfx + (psfy * size)];
				float synthfactor = (orig < sat) ? 1.f -
					((sat - orig) / sat) : 1.f;
				Lsynth[xx + ((dimy - yy) * dimx)] += max(synth * synthfactor
					+ orig * (1 - synthfactor) - sat, 0.f);
				buf[psfx + psfy * size] = Lsynth[xx + ((dimy - yy) * dimx)];
			}
		}
	}
	memcpy(resbuf, buf, size * size * sizeof(float));

// Carry out median blur of middle part, storing the result in resbuf
// in order not to overwrite data in buf that is still needed as input
// for remaining pixel calculations
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
	for (int i = halfpsfdim * 2/3 ; i < halfpsfdim * 4/3 ; i++) {
		for (int j = halfpsfdim * 2/3 ; j < halfpsfdim * 4/3 ; j++) {
			int il = i - 1;
			int iu = i + 1;
			int jl = (j - 1) * size;
			int jm = j * size;
			int ju = (j + 1) * size;
			int offx = i - halfpsfdim;
			int offy = j - halfpsfdim;
			int rad = halfpsfdim * 2/3;
			// Only blur within a circle of radius rad
			if (offx * offx + offy * offy <= rad * rad)
				resbuf[i + jm] = median9f(
					buf[il + jl],
					buf[i + jl],
					buf[iu + jl],
					buf[il + jm],
					buf[i + jm],
					buf[iu + jm],
					buf[il + ju],
					buf[i + ju],
					buf[iu + ju]);
		}
	}

// Copy resbuf back into Lsynth
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
	for (int psfx = size * 2/3; psfx < size * 2/3; psfx++) {
		for (int psfy = size * 2/3; psfy < size * 4/3; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
			if (xx >= 0 && xx < dimx && yy > 0 && yy <= dimy) {
				Lsynth[xx + ((dimy - yy) * dimx)] = resbuf[psfx + psfy * size];
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
}
#endif
	free(buf);
	free(resbuf);
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

gpointer fix_saturated_stars(gpointer data) {
	// Remove unused argument warnings
	(void) data;
	reprofile_saturated_stars(&gfit);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(0);
}

gpointer do_synthstar(gpointer data) {
	// Remove unused argument warnings
	(void) data;
	generate_synthstars(&gfit);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(0);
}

int generate_synthstars(fits *fit) {
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	char *msg = siril_log_color_message(_("Star synthesis (full star mask creation): processing...\n"), "green");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	gboolean is_RGB = TRUE;
	gboolean is_32bit = TRUE;
	gboolean stars_needs_freeing = FALSE;
	float norm = 1.0f, invnorm = 1.0f;
	int nb_stars = 0;
	psf_star **stars = NULL;

	if (starcount(com.stars) < 1) {
		// Set up starfinder_data structure
		struct starfinder_data *sf_data = calloc(1, sizeof(struct starfinder_data));
		if (!sf_data) {
			siril_log_color_message(_("Memory allocation failed\n"), "red");
			set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
			return -1;
		}

		sf_data->im.fit = fit;
		sf_data->im.from_seq = NULL;
		sf_data->im.index_in_seq = -1;
		sf_data->layer = (fit->naxes[2] == 1) ? 0 : 1;
		sf_data->max_stars_fitted = MAX_STARS;
		sf_data->selection = (rectangle){0, 0, 0, 0}; // no selection
		sf_data->save_eqcoords = FALSE;
		sf_data->ref_wcs = NULL;
		sf_data->stars = &stars;
		sf_data->nb_stars = &nb_stars;
		sf_data->threading = MULTI_THREADED;
		sf_data->update_GUI = FALSE;
		sf_data->process_all_images = FALSE;
		sf_data->already_in_thread = TRUE;
		sf_data->keep_stars = FALSE;

		// Call the worker function
		int retval = GPOINTER_TO_INT(findstar_worker(sf_data));
		free(sf_data);

		if (retval != 0 || !stars) {
			siril_log_color_message(_("Star detection failed\n"), "red");
			set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
			if (stars)
				free_fitted_stars(stars);
			return -1;
		}
		stars_needs_freeing = TRUE;
	} else {
		stars = com.stars;
		nb_stars = starcount(com.stars);
	}

	if (nb_stars < 1 || !stars) {
		siril_log_color_message(_("No stars detected in the image.\n"), "red");
		if (stars_needs_freeing)
			free_fitted_stars(stars);
		set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
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
	int npixels = dimx * dimy;
	gboolean buf_needs_freeing = FALSE;
	// Regardless of 16/32bit store the data in a buffer, converting if needed
	float *buf[3];
	if (is_RGB) {
		if (is_32bit) {
			buf[RLAYER] = fit->fpdata[RLAYER];
			buf[GLAYER] = fit->fpdata[GLAYER];
			buf[BLAYER] = fit->fpdata[BLAYER];
		} else {
			buf[RLAYER] = (float*) calloc(npixels, sizeof(float));
			buf[GLAYER] = (float*) calloc(npixels, sizeof(float));
			buf[BLAYER] = (float*) calloc(npixels, sizeof(float));
			buf_needs_freeing = TRUE;
			for (size_t i = 0; i < npixels; i++) {
				buf[RLAYER][i] = (float) fit->pdata[RLAYER][i] * invnorm;
				buf[GLAYER][i] = (float) fit->pdata[GLAYER][i] * invnorm;
				buf[BLAYER][i] = (float) fit->pdata[BLAYER][i] * invnorm;
			}
		}
	} else { // mono
		if (is_32bit)
			buf[RLAYER] = fit->fdata;
		else {
			buf[RLAYER] = (float*) calloc(npixels, sizeof(float));
			buf_needs_freeing = TRUE;
			for (size_t i = 0; i < npixels; i++)
				buf[RLAYER][i] = (float) fit->data[i] * invnorm;
		}
	}

	// Normalize the buffer to avoid issues with colorspace conversion
	float bufmax = 1.f;
	for (size_t chan = 0; chan < fit->naxes[2]; chan++)
		for (size_t i = 0; i < npixels; i++)
			if (buf[chan][i] > bufmax)
				bufmax = buf[chan][i];
	for (size_t chan = 0; chan < fit->naxes[2]; chan++)
		for (size_t i = 0; i < npixels; i++)
			buf[chan][i] /= bufmax;

	float *H = NULL, *S = NULL, *Hsynth = NULL, *Ssynth = NULL, *Lsynth, junk;
	Lsynth = (float*) calloc(npixels, sizeof(float));

	// For RGB images, convert pixel colour data from fit into H and S arrays. L is irrelevant as we will synthesize L.
	if (is_RGB) {
		H = (float*) calloc(npixels, sizeof(float));
		S = (float*) calloc(npixels, sizeof(float));
		Hsynth = (float*) calloc(npixels, sizeof(float));
		Ssynth = (float*) calloc(npixels, sizeof(float));
		for (size_t i = 0; i < npixels; i++) {
			rgb_to_hsl_float_sat(buf[RLAYER][i], buf[GLAYER][i],
					buf[BLAYER][i], 0.f, &H[i], &S[i], &junk);
		}
	}
	// Calculate average Moffat beta
	size_t moffat_count = 0;
	double avg_moffat_beta = 0.;

	gboolean stopcalled = FALSE;
	// Synthesize a PSF for each star in the star array s, based on its measured parameters
	gboolean gaussian = TRUE;
	if (stars[0]->profile == PSF_MOFFAT_BFREE)
		gaussian = FALSE;
	if (!gaussian) {
		for (size_t n = 0 ; n < nb_stars ; n++) {
			moffat_count++;
			avg_moffat_beta += stars[n]->beta;
		}
		avg_moffat_beta /= moffat_count;
		siril_debug_print("# Moffat profile stars: %zd, average beta = %.3f\n", moffat_count, avg_moffat_beta);
	}
	for (int n = 0; n < nb_stars; n++) {
		// Check if stop has been pressed
		if (!get_thread_run())
			stopcalled = TRUE;
		set_progress_bar_data(NULL,	(double) n / (double) nb_stars);
		if (!stopcalled) {
			float lum = (float) stars[n]->A;
			if (lum < 0.0f)
				lum = 0.0f;
			if (!is_32bit)
				lum *= invnorm;
			assert(lum >= 0.0f);
			float xoff = (float) stars[n]->xpos - (int) stars[n]->xpos;
			float yoff = (float) stars[n]->ypos - (int) stars[n]->ypos;
			int size = (int) 5 * max(stars[n]->fwhmx, stars[n]->fwhmy); // This is big enough that even under extreme stretching the synthesized psf tails off smoothly
			if (!gaussian)
				size *= 10 / stars[n]->beta; // Increase the PSF size markedly for low beta stars
			if (!(size % 2))
				size++;
			if (size > 1024) // protect against excessive memory allocations due to bad parameters;
							 // 100px should be more than enough for the fwhm of even a very saturated star
				continue;
			float minfwhm = min(stars[n]->fwhmx, stars[n]->fwhmy);

			// Synthesize the luminance profile and add to the star mask in HSL colourspace
			float *psfL = (float*) calloc(size * size, sizeof(float));
			if (!psfL) // May happen if size is excessively large because of a bad fwhm value
				continue;
			float beta = 8.f;
			if (!gaussian) {
				if (stars[n]->beta > 0.0) {
					beta=stars[n]->beta;
					minfwhm = min(stars[n]->fwhmx, stars[n]->fwhmy);
				} else if (moffat_count > 0)
					beta = avg_moffat_beta;
			}
			if (stars[n]->has_saturated || gaussian)
				makegaussian(psfL, size, minfwhm, lum, xoff, yoff, 1.f, 0.f);
			else
				makemoffat(psfL, size, minfwhm, lum, xoff, yoff, beta, 1.f, 0.f);
			if (is_RGB)
				add_star_to_rgb_buffer(H, S, psfL, size, Hsynth, Ssynth, Lsynth,
						(float) stars[n]->xpos, (float) stars[n]->ypos, dimx, dimy);
			else
				add_star_to_mono_buffer(psfL, size, Lsynth, (float) stars[n]->xpos,
						(float) stars[n]->ypos, dimx, dimy);
			free(psfL);
		}
	}
	// Stars are only freed if they were *not* taken from com.stars: if the
	// user has made a specific selection of stars, we want to leave that
	// selection intact.
	if (stars_needs_freeing)
		free_fitted_stars(stars);

	// Construct the RGB from synthetic L (and for RGB images, also the H and S values from the orginal image thus giving our synthesized stars the correct colour)
	if (!stopcalled) {
		if (is_RGB) {
			float *R, *G, *B;
			R = (float*) calloc(npixels, sizeof(float));
			G = (float*) calloc(npixels, sizeof(float));
			B = (float*) calloc(npixels, sizeof(float));
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) if (com.max_thread > 1)
{
#endif
			float bufmaxx = 1.f;
			for (size_t i = 0; i < npixels; i++)
				if (Lsynth[i] > bufmaxx)
					bufmaxx = Lsynth[i];
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
			for (size_t i = 0; i < npixels; i++)
				Lsynth[i] /= bufmaxx;

#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
			for (size_t n = 0; n < npixels; n++) {
				hsl_to_rgb_float_sat(Hsynth[n], Ssynth[n], Lsynth[n], &R[n],
						&G[n], &B[n]);
				// Trap NaNs and infinities
				R[n] = (isnan(R[n]) || isinf(R[n]) || R[n] < 0.f) ? 0.f : R[n];
				G[n] = (isnan(G[n]) || isinf(G[n]) || G[n] < 0.f) ? 0.f : G[n];
				B[n] = (isnan(B[n]) || isinf(B[n]) || B[n] < 0.f) ? 0.f : B[n];
			}
			if (is_32bit) {
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
				for (size_t n = 0; n < npixels; n++) {
					fit->fpdata[RLAYER][n] = R[n];
					fit->fpdata[GLAYER][n] = G[n];
					fit->fpdata[BLAYER][n] = B[n];
				}
			} else {
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
				for (size_t n = 0; n < npixels; n++) {
					fit->pdata[RLAYER][n] = roundf_to_WORD(R[n] * norm);
					fit->pdata[GLAYER][n] = roundf_to_WORD(G[n] * norm);
					fit->pdata[BLAYER][n] = roundf_to_WORD(B[n] * norm);
				}
			}
#ifdef _OPENMP
}
#endif
			// Free memory
			free(R);
			free(G);
			free(B);
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
					for (size_t n = 0; n < npixels; n++) {
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
					for (size_t n = 0; n < npixels; n++) {
						fit->data[n] = roundf_to_WORD(Lsynth[n] * norm);
					}
				}
#ifdef _OPENMP
			}
#endif
		}
	}
	if (H != NULL)
		free(H);
	if (S != NULL)
		free(S);
	if (Hsynth != NULL)
		free(Hsynth);
	if (Ssynth != NULL)
		free(Ssynth);
	free(Lsynth);
	if (buf_needs_freeing) {
		if (is_RGB) {
			for (size_t i = 0; i < 3; i++)
				free(buf[i]);
		} else
			free(buf[RLAYER]);
	}
	update_filter_information(fit, "StarMask", TRUE);
	if (fit == &gfit && !stopcalled)
		notify_gfit_modified();
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, "Execution time");
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	return 0;
}

// Fix up saturated stars only

int reprofile_saturated_stars(fits *fit) {
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	char *msg = siril_log_color_message(_("Star synthesis (desaturating clipped star profiles): processing...\n"), "green");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	gboolean is_RGB = (fit->naxes[2] == 3) ? TRUE : FALSE;
	gboolean is_32bit = TRUE;
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

	buf[RLAYER] = malloc(count * sizeof(float));
	if (is_RGB) {
		buf[GLAYER] = malloc(count * sizeof(float));
		buf[BLAYER] = malloc(count * sizeof(float));
		if (is_32bit) {
			memcpy(buf[RLAYER], fit->fpdata[RLAYER], fit->rx * fit->ry * sizeof(float));
			memcpy(buf[GLAYER], fit->fpdata[GLAYER], fit->rx * fit->ry * sizeof(float));
			memcpy(buf[BLAYER], fit->fpdata[BLAYER], fit->rx * fit->ry * sizeof(float));
		} else {
			for (size_t i = 0; i < count; i++) {
				buf[RLAYER][i] = (float) fit->pdata[RLAYER][i] * invnorm;
				buf[GLAYER][i] = (float) fit->pdata[GLAYER][i] * invnorm;
				buf[BLAYER][i] = (float) fit->pdata[BLAYER][i] * invnorm;
			}
		}
	} else { // mono
		if (is_32bit)
			memcpy(buf[RLAYER], fit->fdata, fit->rx * fit->ry * sizeof(float));
		else {
			buf[RLAYER] = malloc(count * sizeof(float));
			for (size_t i = 0; i < count; i++)
				buf[RLAYER][i] = (float) fit->data[i] * invnorm;
		}
	}

	// Set up starfinder_data structure once, we will reuse it for each channel
	struct starfinder_data sf_data = { 0 };
	sf_data.im.fit = fit;
	sf_data.im.from_seq = NULL;
	sf_data.im.index_in_seq = -1;
	sf_data.max_stars_fitted = MAX_STARS;
	sf_data.selection = (rectangle){0, 0, 0, 0}; // no selection
	sf_data.save_eqcoords = FALSE;
	sf_data.ref_wcs = NULL;
	sf_data.threading = MULTI_THREADED;
	sf_data.update_GUI = FALSE;
	sf_data.process_all_images = FALSE;
	sf_data.already_in_thread = TRUE;
	sf_data.keep_stars = FALSE;

	// Synthesize a PSF for each saturated star in the star array, based on its measured parameters.
	// To fix saturated star profiles we have to do this for each color channel as we can't rely on
	// the hue and saturation within the saturated area, whereas the profiles will be accurate.
	gboolean stopcalled = FALSE;
	for (size_t chan = 0; chan < fit->naxes[2]; chan++) {
		if (stopcalled)
			break;

		psf_star **stars = NULL;
		int nb_stars = 0;

		// Update only the channel-specific fields
		sf_data.layer = chan;
		sf_data.stars = &stars;
		sf_data.nb_stars = &nb_stars;

		// Call the worker function
		int retval = GPOINTER_TO_INT(findstar_worker(&sf_data));

		if (retval != 0 || !stars) {
			siril_log_color_message(_("Star detection failed for channel %u\n"), "red", chan);
			if (stars)
				free_fitted_stars(stars);
			continue; // Skip this channel but continue with others
		}

		int sat_stars = 0;
		siril_log_message(_("Star synthesis: desaturating stars in channel %u...\n"), chan);
		double total = fit->naxes[2] * nb_stars;
		for (size_t n = 0; n < nb_stars; n++) {
			// Check if stop has been pressed
			if (!get_thread_run())
				stopcalled = TRUE;
			set_progress_bar_data(NULL, (double) (n * fit->naxes[2] + chan) / total);
			if (stars[n]->has_saturated && !stopcalled) {
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
				if (size > 1024)
					size = 1024; // Protection against bad star params
				float ratio = stars[n]->fwhmx / stars[n]->fwhmy;
				float angle = (float) stars[n]->angle;

				float *psfL = (float*) calloc(size * size, sizeof(float));
				if (!psfL)
					continue;
				makegaussian(psfL, size, stars[n]->fwhmx, (lum - bg), xoff, yoff, ratio, angle);

				// Replace the part of the profile above the sat threshold
				replace_sat_star_in_buffer(psfL, size, buf[chan],
						(float) stars[n]->xpos, (float) stars[n]->ypos, dimx,
						dimy, sat, bg, 0.f);
				free(psfL);
				sat_stars++;
			}
		}
		free_fitted_stars(stars);
		siril_log_message(_("Star synthesis: %d stars desaturated\n"), sat_stars);
	}

	// Desaturating stars will take their peak brightness over 1.f so we need to rescale the values
	// of all pixels by a factor of (1 / maxbuf) where maxbuf is the maximum subpixel value across all channels
	if (!stopcalled) {
		float bufmax = 1.f;
		for (size_t chan = 0; chan < fit->naxes[2]; chan++)
			for (size_t i = 0; i < count; i++)
				if (buf[chan][i] > bufmax)
					bufmax = buf[chan][i];
		if (bufmax > 1.f) {
			float invbufmax = 1.f / bufmax;
			siril_log_message(_("Remapping output to floating point range 0.0 to 1.0\n"));
			for (size_t chan = 0; chan < fit->naxes[2]; chan++)
				for (size_t i = 0; i < count; i++)
					buf[chan][i] *= invbufmax;
		}
		if (is_32bit) {
			if (is_RGB) {
				memcpy(fit->fpdata[RLAYER], buf[RLAYER], fit->rx * fit->ry * sizeof(float));
				memcpy(fit->fpdata[GLAYER], buf[GLAYER], fit->rx * fit->ry * sizeof(float));
				memcpy(fit->fpdata[BLAYER], buf[BLAYER], fit->rx * fit->ry * sizeof(float));
			}
			else {
				memcpy(fit->fdata, buf[RLAYER], fit->rx * fit->ry * sizeof(float));
			}
		} else {
			for (size_t n = 0; n < count; n++) {
				if (is_RGB) {
					fit->pdata[RLAYER][n] = roundf_to_WORD(buf[RLAYER][n] * norm);
					fit->pdata[GLAYER][n] = roundf_to_WORD(buf[GLAYER][n] * norm);
					fit->pdata[BLAYER][n] = roundf_to_WORD(buf[BLAYER][n] * norm);
				}
				else {
					fit->data[n] = roundf_to_WORD(buf[RLAYER][n] * norm);
				}
			}
		}
	}
	if (is_RGB) {
		for (size_t i = 0; i < 3; i++)
			free(buf[i]);
	} else
		free(buf[RLAYER]);

	if (fit == &gfit && !stopcalled)
		notify_gfit_modified();
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, "Execution time");
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	return 0;
}
