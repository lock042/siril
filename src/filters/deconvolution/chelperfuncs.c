/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define dontneedcppfftwmultithreaded
#define dontneedcppfftwtimelimit
#define dontneedcppfftwflags
#include "chelperfuncs.h"
#undef dontneedcppfftwmultithreaded
#undef dontneedcppfftwtimelimit
#undef dontneedcppfftwflags
#include "gui/progress_and_log.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "core/siril_log.h"

#ifndef RT_INCLUDE
#undef max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#undef min
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

// Helpful wrapper functions to call Siril functions without causing clashes with some of the img_expr functionality

void updateprogress(const char *text, double percent) {
	set_progress_bar_data(text, percent);
}

int is_thread_stopped() {
	return (get_thread_run() ? 0 : 1);
}

void sirillog(const char* text) {
	siril_log_message(_(text));
}

int updatenoise(float *array, int nx, int ny, int nchans, double *noise) {
	return sos_update_noise_float(array, (long) nx, (long) ny, (long) nchans, noise);
}

// Scaling functions

inline static float sample(const float *buf, int w, int h, int d, int i, int j, int c) {
    i = min(max(i, 0), w - 1);
    j = min(max(j, 0), h - 1);
    c = min(max(c, 0), d - 1);
	return buf[(i+j*w)*d + c];
}

static float cubic(float p[4][4], float x, float y)
{
	float v[4];
    v[0] = p[0][1] + 0.5f * y*(p[0][2] - p[0][0] + y*(2.f*p[0][0] - 5.f*p[0][1] + 4.f*p[0][2] - p[0][3] + y*(3.f*(p[0][1] - p[0][2]) + p[0][3] - p[0][0])));
    v[1] = p[1][1] + 0.5f * y*(p[1][2] - p[1][0] + y*(2.f*p[1][0] - 5.f*p[1][1] + 4.f*p[1][2] - p[1][3] + y*(3.f*(p[1][1] - p[1][2]) + p[1][3] - p[1][0])));
    v[2] = p[2][1] + 0.5f * y*(p[2][2] - p[2][0] + y*(2.f*p[2][0] - 5.f*p[2][1] + 4.f*p[2][2] - p[2][3] + y*(3.f*(p[2][1] - p[2][2]) + p[2][3] - p[2][0])));
    v[3] = p[3][1] + 0.5f * y*(p[3][2] - p[3][0] + y*(2.f*p[3][0] - 5.f*p[3][1] + 4.f*p[3][2] - p[3][3] + y*(3.f*(p[3][1] - p[3][2]) + p[3][3] - p[3][0])));
    return v[1] + 0.5f * x*(v[2] - v[0] + x*(2.f*v[0] - 5.f*v[1] + 4.f*v[2] - v[3] + x*(3.f*(v[1] - v[2]) + v[3] - v[0])));
}

void magnify(float *out, const float *in, int out_w, int out_h, int d, int in_w, int in_h, float factor) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2) num_threads(cppmaxthreads)
#endif
	for (int j = 0; j < out_h; j++) {
        for (int i = 0; i < out_w; i++) {
            float tmp[d];
            float x = i / factor - 1;
            float y = j / factor - 1;
            int ix = floor(x);
            int iy = floor(y);
            for (int c = 0; c < d; c++) {
                float block[4][4];
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        block[i][j] = sample(in, in_w, in_h, d, ix + i, iy + j, c);
					}
				}
                tmp[c] = cubic(block, x - ix, y - iy);
                if (!(i < 0 || i >= out_w || j < 0 || j >= out_h || c < 0 || c >= d)) {
					out[(i+j*out_w)*d + c] = tmp[c];
				}
            }
        }
    }
}

void gaussblur(float *y, float *x, int w, int h, float sigma) {
	fftwf_complex *fx = fftwf_malloc(w * h * sizeof(fftwf_complex));
	fftwf_complex *fk = fftwf_malloc(w * h * sizeof(fftwf_complex));
	fftwf_complex *a = fftwf_malloc(w * h * sizeof(fftwf_complex));
	fftwf_plan p = fftwf_plan_dft_2d(h, w, a, fx, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_plan q = fftwf_plan_dft_2d(h, w, a, fk, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_plan r = fftwf_plan_dft_2d(h, w, a, fk, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (int i = 0 ; i < (w * h); i++) {
		a[i] = x[i];
	}
	float *k = malloc(w * h * sizeof(float));
	float (*kk)[w] = (void *)k;
	float invss = 1/(sigma * sigma);
	float alpha = invss / (M_PI);
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,16) collapse(2) num_threads(cppmaxthreads)
#endif
	for (int j = 0 ; j < h ; j++) {
		for (int i = 0 ; i < w ; i++) {
			float r = hypotf(i < w / 2 ? i : i - w, j < h / 2 ? j : j - h);
			kk[j][i] = alpha * expf(-r*r*invss);
		}
	}
	float sum = 0;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,16) collapse(2) num_threads(cppmaxthreads)
#endif
	for (int j = 0 ; j < h ; j++) {
		for (int i = 0 ; i < w ; i++) {
			sum += kk[j][i];
		}
	}
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,16) collapse(2) num_threads(cppmaxthreads)
#endif
	for (int j = 0 ; j < h ; j++) {
		for (int i = 0 ; i < w ; i++) {
			kk[j][i] /= sum;
		}
	}
	for (int i = 0 ; i < (w * h); i++) {
		a[i] = k[i];
	}
	free(k);
	fftwf_execute(p);
	fftwf_destroy_plan(p);
	for (int i = 0 ; i < w * h ; i++) {
		fx[i] = fx[i] * fk[i];
	}
	fftwf_execute(q);
	fftwf_destroy_plan(q);
	for (int i = 0 ; i < w * h; i++) {
		a[i] = fx[i];
	}
	fftwf_free(fx);
	fftwf_execute(r); // Reuse fk here to hold the output
	fftwf_destroy_plan(r);
	float scale = 1.f / (w * h);
	for (int i = 0 ; i < w * h; i++) {
		fftwf_complex z = fk[i] * scale;
		y[i] = crealf(z);
	}
	fftwf_free(fk);
	fftwf_free(a);
}

static float bilinear(float *x, int w, int h, float i, float j) {
	int ii = (int) i;
	int jj = (int) j;
	float (*xx)[w] = (void*)x;
	float a = xx[min(max(jj, 0), h - 1)][min(max(ii, 0), w - 1)];
	float b = xx[min(max(jj + 1, 0), h - 1)][min(max(ii, 0), w - 1)];
	float c = xx[min(max(jj, 0), h - 1)][min(max(ii + 1, 0), w - 1)];
	float d = xx[min(max(jj + 1, 0), h - 1)][min(max(ii + 1, 0), w - 1)];
	float xoff = i-ii;
	float yoff = j-jj;
	return (a * (1 - xoff) * (1 - yoff)) + (b * (1 - xoff) * yoff) + (c * xoff * (1 - yoff)) + (d * xoff * yoff);
}

void shrink(float *out, float *in, int outw, int outh, int inw, int inh, float scale, float sigma) {
	if (scale == -2) {
		assert(2*outw == inw && 2 * outh == inh);
		float (*y)[outw] = (void*)out;
		float (*x)[inw] = (void*)in;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static, 16) collapse(2) num_threads(cppmaxthreads)
#endif
		for (int j = 0; j < outh; j++) {
			for (int i = 0; i < outw; i++) {
				float g = x[2 * j][2 * i] + x[2 * j][2 * i + 1] + x[ 2 * j + 1][2 * i] + x[2 * j + 1][ 2 * i + 1];
				y[j][i] = g/4;
			}
		}
		return;
	}
	float xf = inw/(float)outw;
	float yf = inh/(float)outh;
	float blur_size = sigma * sqrtf((xf * yf - 1) / 3);
	float *gaussian = malloc(inw * inh * sizeof(float));
	if (outw < inw || outh < inh) {
		gaussblur(gaussian, in, inw, inh, blur_size);
	} else {
		for (int i = 0; i < inw * inh; i++)
			gaussian[i] = in[i];
	}
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,16) collapse(2) num_threads(cppmaxthreads)
#endif
	for (int j = 0; j < outh; j++) {
		for (int i = 0; i < outw; i++) {
			float x = xf * i;
			float y = yf * j;
			out[outw * j + i] = bilinear(gaussian, inw, inh, x, y);
		}
	}
	free(gaussian);
}
