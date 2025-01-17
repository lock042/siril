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
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"

#include "arithm.h"

/*****************************************************************************
 *       S I R I L      A R I T H M E T I C      O P E R A T I O N S         *
 ****************************************************************************/

static int soper_ushort_to_ushort(fits *a, float scalar, image_operator oper) {
	if (!a) return 1;
	if (!a->data) return 1;
	WORD *data = NULL;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	data = a->data;
	float invnorm = (a->bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : 1.f;
	if (oper == OPER_DIV) {
		scalar = 1.0f / scalar;
		oper = OPER_MUL;
	}
	if (oper == OPER_MUL) scalar *= invnorm;

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				if (invnorm == 1.f) {
					data[i] = float_to_ushort_range(pixel + scalar);
				} else {
					data[i] = float_to_ushort_range(invnorm * (pixel + scalar));
				}
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				if (invnorm == 1.f) {
					data[i] = float_to_ushort_range(pixel - scalar);
				} else {
					data[i] = float_to_ushort_range(invnorm * (pixel - scalar));
				}
			}
			break;
		case OPER_MUL:
		default:
			for (i = 0; i < n; ++i) {
				data[i] = roundf_to_WORD((float)data[i] * scalar);
			}
			break;
	}
	invalidate_stats_from_fit(a);
	return 0;
}

static int soper_ushort_to_float(fits *a, float scalar, image_operator oper) {
	if (!a) return 1;
	if (!a->data) return 1;
	WORD *data = NULL;
	float *result= NULL;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	data = a->data;
	result = malloc(n * sizeof(float));
	if (!result) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	if (oper == OPER_DIV) {
		scalar = 1.0f / scalar;
		oper = OPER_MUL;
	}

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				result[i] = pixel + scalar;
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				result[i] = pixel - scalar;
			}
			break;
		case OPER_MUL:
		default:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				result[i] = pixel * scalar;
			}
			break;
	}
	fit_replace_buffer(a, result, DATA_FLOAT);
	return 0;
}

int soper_unscaled_div_ushort_to_float(fits *a, int scalar) {
	if (!a) return 1;
	if (!a->data) return 1;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	WORD *data = a->data;
	float *result = malloc(n * sizeof(float));
	if (!result) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	float operand = 1.0f / (float)scalar;
	for (i = 0; i < n; ++i) {
		result[i] = data[i] * operand;
	}
	fit_replace_buffer(a, result, DATA_FLOAT);
	return 0;
}

static int soper_float(fits *a, float scalar, image_operator oper) {
	if (!a) return 1;
	if (!a->fdata) return 1;
	float *data = a->fdata;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	if (oper == OPER_DIV) {
		scalar = 1.0f / scalar;
		oper = OPER_MUL;
	}

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				data[i] = data[i] + scalar;
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				data[i] = data[i] - scalar;
			}
			break;
		case OPER_MUL:
		default:
			for (i = 0; i < n; ++i) {
				data[i] = data[i] * scalar;
			}
			break;
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* equivalent to (map simple_operation a), with simple_operation being
 * (lambda (pixel) (oper pixel scalar))
 * scalar is applied on images with data in [0, 1]
 */
int soper(fits *a, float scalar, image_operator oper, gboolean conv_to_float) {
	if (oper == OPER_DIV && scalar == 0.f) {
		siril_log_message(_("Cannot divide by zero, aborting."));
		return 1;
	}
	if (a->type == DATA_USHORT) {
		if (conv_to_float)
			return soper_ushort_to_float(a, scalar, oper);
		return soper_ushort_to_ushort(a, scalar, oper);
	}
	if (a->type == DATA_FLOAT)
		return soper_float(a, scalar, oper);
	return 1;
}

static int imoper_to_ushort(fits *a, fits *b, image_operator oper, float factor) {
	if (!a) return 1;
	if (!a->data) return 1;
	if (!b) return 1;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		return 1;
	}

	if (b->type == DATA_USHORT) {
		if (!b->data) return 1;
		WORD *abuf = a->data, *bbuf = b->data;
		if (oper == OPER_DIV) {
			for (i = 0; i < n; ++i) {
				if (bbuf[i] == 0)
					abuf[i] = 0;
				else {
					float aval = (float) abuf[i];
					float bval = (float) bbuf[i];
					if (factor != 1.0f)
						abuf[i] = roundf_to_WORD(factor * (aval / bval));
					else
						abuf[i] = roundf_to_WORD(aval / bval);
				}
			}
		} else if (oper == OPER_MUL) {
			for (i = 0; i < n; ++i) {
				if (bbuf[i] == 0)
					abuf[i] = 0;
				else {
					float aval = (float) abuf[i];
					float bval = (float) bbuf[i];
					if (factor != 1.0f)
						abuf[i] = roundf_to_WORD(factor * (aval * bval));
					else
						abuf[i] = roundf_to_WORD(aval * bval);
				}
			}
		} else {
			for (i = 0; i < n; ++i) {
				int aval = (int) abuf[i];
				int bval = (int) bbuf[i];
				switch (oper) {
				case OPER_ADD:
					abuf[i] = truncate_to_WORD(aval + bval);
					break;
				case OPER_SUB:
					abuf[i] = truncate_to_WORD(aval - bval);
					break;
				default:	// OPER_MUL, OPER_DIV handled above
					break;
				}

				if (factor != 1.0f)
					abuf[i] = roundf_to_WORD(factor * (float) abuf[i]);
			}
		}
	} else if (b->type == DATA_FLOAT) {
		if (!b->fdata) return 1;
		WORD *abuf = a->data;
		float *bbuf = b->fdata;
		float norm = (a->bitpix == BYTE_IMG) ? UCHAR_MAX_SINGLE : USHRT_MAX_SINGLE;

		if (oper == OPER_DIV) {
			for (i = 0; i < n; ++i) {
				if (bbuf[i] == 0.f)
					abuf[i] = 0;
				else {
					float aval = (float) abuf[i];
					float bval = bbuf[i] * norm;
					if (factor != 1.0f)
						abuf[i] = roundf_to_WORD(factor * (aval / bval));
					else
						abuf[i] = roundf_to_WORD(aval / bval);
				}
			}
		} else if (oper == OPER_MUL) {
			for (i = 0; i < n; ++i) {
				if (bbuf[i] == 0.f)
					abuf[i] = 0;
				else {
					float aval = (float) abuf[i];
					float bval = bbuf[i] * norm;
					if (factor != 1.0f)
						abuf[i] = roundf_to_WORD(factor * (aval * bval));
					else
						abuf[i] = roundf_to_WORD(aval * bval);
				}
			}
		} else {
			for (i = 0; i < n; ++i) {
				int aval = (int) abuf[i];
				int bval = roundf_to_int(bbuf[i] * norm);
				switch (oper) {
				case OPER_ADD:
					abuf[i] = truncate_to_WORD(aval + bval);
					break;
				case OPER_SUB:
					abuf[i] = truncate_to_WORD(aval - bval);
					break;
				default:	// OPER_MUL and OPER_DIV handled above
					break;
				}
				if (factor != 1.0f)
					abuf[i] = roundf_to_WORD(factor * (float) abuf[i]);
			}
		}
	}
	invalidate_stats_from_fit(a);
	a->neg_ratio = 0.0f;
	return 0;
}

int imoper_to_float(fits *a, fits *b, image_operator oper, float factor) {
	if (!a) return 1;
	if (!b) return 1;
	size_t n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	float *result = NULL;

	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		return 1;
	}

	if (a->type == DATA_FLOAT) {
		if (!a->fdata) return 1;
		result = a->fdata;
	}
	else if (a->type == DATA_USHORT) {
		result = malloc(n * sizeof(float));
		if (!result) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}
	else return 1;

	if (b->type == DATA_USHORT && (!b->data)) { free(result); return 1; }
	if (b->type == DATA_FLOAT && (!b->fdata)) { free(result); return 1; }
	if (!(b->type == DATA_FLOAT || b->type == DATA_USHORT)) { free(result); return 1; }

	size_t nb_negative = 0;
	for (size_t i = 0; i < n; ++i) {
		float aval = a->type == DATA_USHORT ? ushort_to_float_bitpix(a, a->data[i]) : a->fdata[i];
		float bval = b->type == DATA_USHORT ? ushort_to_float_bitpix(b, b->data[i]) : b->fdata[i];
		switch (oper) {
			case OPER_ADD:
				result[i] = aval + bval;
				break;
			case OPER_SUB:
				result[i] = aval - bval;
				break;
			case OPER_MUL:
				result[i] = aval * bval;
				break;
			case OPER_DIV:
				if (bval == 0.0f)
					result[i] = 0.0f;
				else result[i] = aval / bval;
		}
		if (factor != 1.0f)
			result[i] *= factor;
		if (result[i] > 1.0f)	// should we truncate by default?
			result[i] = 1.0f;
		if (result[i] < -1.0f)	// Here we want to clip all garbages pixels.
			result[i] = 0.0f;
		if (result[i] < 0.0f)
			nb_negative++;
	}
	a->neg_ratio = (float)((double)nb_negative / n);
	if (a->type == DATA_USHORT)
		fit_replace_buffer(a, result, DATA_FLOAT);
	else invalidate_stats_from_fit(a);
	return 0;
}

/* applies operation of image a with image b, for all their layers:
 * a = factor * a oper b
 * returns 0 on success */
static int imoper_with_factor(fits *a, fits *b, image_operator oper, float factor, gboolean allow_32bits) {
	// ushort result can only be forced when both input images are ushort
	if (allow_32bits && !com.pref.force_16bit)
		return imoper_to_float(a, b, oper, factor);
	else {
		if (a->type == DATA_USHORT)
			return imoper_to_ushort(a, b, oper, factor);
		siril_log_color_message(_("Image operations can only be kept 16 bits if first input images are 16 bits. Aborting.\n"), "red");
	}
	return 1;
}

int imoper(fits *a, fits *b, image_operator oper, gboolean allow_32bits) {
	return imoper_with_factor(a, b, oper, 1.0f, allow_32bits);
}

/* a = coef * a / b
 * a is expected as USHORT, returned as FLOAT. b is expected as USHORT.
 * If overflow, siril_fdiv returns 1*/
int siril_fdiv(fits *a, fits *b, float coef, gboolean allow_32bits) {
	return imoper_with_factor(a, b, OPER_DIV, coef, allow_32bits);
}

// a = max(a, b)
int addmax(fits *a, fits *b) {
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];

	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		return 1;
	}
	if (a->type != b->type) {
		siril_log_color_message(_("Images must have same data type.\n"), "red");
		return 1;
	}
	g_assert(a->naxes[2] == 1 || a->naxes[2] == 3);

	if (a->type == DATA_USHORT) {
		WORD *abuf = a->data, *bbuf = b->data;
		for (i = 0; i < n; ++i) {
			if (bbuf[i] > abuf[i])
				abuf[i] = bbuf[i];
		}
	} else {
		float *abuf = a->fdata, *bbuf = b->fdata;
		for (i = 0; i < n; ++i) {
			if (bbuf[i] > abuf[i])
				abuf[i] = bbuf[i];
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

void rgbblend(blend_data *data, float* r, float* g, float* b, float m_CB) {
	float *out[3] = { r, g, b };

	// Compute the maximum values using fmaxf for better vectorization
	float sfmax = fmaxf(fmaxf(data->sf[0], data->sf[1]), data->sf[2]);
	float tfmax = fmaxf(fmaxf(data->tf[0], data->tf[1]), data->tf[2]);

	// Compute the difference
	float d = sfmax - tfmax;

	// Calculate k based on conditions
	float k = (tfmax + m_CB * d > 1.f) ? ((d != 0.f) ? fminf(m_CB, (1.f - tfmax) / d) : m_CB) : m_CB;

	// Loop over channels
	for (size_t chan = 0; chan < 3; chan++) {
		// Compute the output using the precomputed k
		*out[chan] = data->do_channel[chan] ? (1.f - k) * data->tf[chan] + k * data->sf[chan] : *out[chan];
	}
}

void clip(fits *fit) {
	if (fit->type == DATA_USHORT) return; // USHORT cannot be out of range
	size_t ndata = fit->rx * fit->ry * fit->naxes[2];
	for (size_t i = 0 ; i < ndata ; i++) {
		fit->fdata[i] = fit->fdata[i] < 0.f ? 0.f : fit->fdata[i] > 1.f ? 1.f : fit->fdata[i];
	}
}

void clipneg(fits *fit) {
	if (fit->type == DATA_USHORT) return; // USHORT cannot be out of range
	size_t ndata = fit->rx * fit->ry * fit->naxes[2];
	for (size_t i = 0 ; i < ndata ; i++) {
		fit->fdata[i] = fit->fdata[i] < 0.f ? 0.f : fit->fdata[i];
	}
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
	int available_threads = com.max_thread - omp_get_num_threads();
#pragma omp parallel for schedule(static) collapse(2) num_threads(available_threads) if (available_threads > 1)
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
	float *k = malloc(w * h * sizeof(float));
	float (*kk)[w] = (void *)k;
	float invss = 1/(sigma * sigma);
	float alpha = invss / (M_PI);
	float scale = 1.f / (w * h);
	float sum = 0;
#ifdef _OPENMP
	int available_threads = com.max_thread - omp_get_num_threads();
#pragma omp parallel num_threads(available_threads) if (available_threads > 1)
{
#pragma omp for simd schedule(static)
#endif
	for (int i = 0 ; i < (w * h); i++) {
		a[i] = x[i];
	}
#ifdef _OPENMP
#pragma omp for simd schedule(static,16) collapse(2)
#endif
	for (int j = 0 ; j < h ; j++) {
		for (int i = 0 ; i < w ; i++) {
			float r = hypotf(i < w / 2 ? i : i - w, j < h / 2 ? j : j - h);
			kk[j][i] = alpha * expf(-r*r*invss);
		}
	}
#ifdef _OPENMP
#pragma omp for simd schedule(static,16) collapse(2)
#endif
	for (int j = 0 ; j < h ; j++) {
		for (int i = 0 ; i < w ; i++) {
			sum += kk[j][i];
		}
	}
#ifdef _OPENMP
#pragma omp for simd schedule(static,16) collapse(2)
#endif
	for (int j = 0 ; j < h ; j++) {
		for (int i = 0 ; i < w ; i++) {
			kk[j][i] /= sum;
		}
	}
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
	for (int i = 0 ; i < (w * h); i++) {
		a[i] = k[i];
	}
#ifdef _OPENMP
#pragma omp single
{
#endif
	free(k);
	fftwf_execute(p);
	fftwf_destroy_plan(p);
#ifdef _OPENMP
}
#pragma omp for simd schedule(static)
#endif
	for (int i = 0 ; i < w * h ; i++) {
		fx[i] = fx[i] * fk[i];
	}
#ifdef _OPENMP
#pragma omp single
{
#endif
	fftwf_execute(q);
	fftwf_destroy_plan(q);
#ifdef _OPENMP
}
#pragma omp for simd schedule(static)
#endif
	for (int i = 0 ; i < w * h; i++) {
		a[i] = fx[i];
	}
#ifdef _OPENMP
#pragma omp single
{
#endif
	fftwf_free(fx);
	fftwf_execute(r); // Reuse fk here to hold the output
	fftwf_destroy_plan(r);
#ifdef _OPENMP
}
#pragma omp for simd schedule(static)
#endif
	for (int i = 0 ; i < w * h; i++) {
		fftwf_complex z = fk[i] * scale;
		y[i] = crealf(z);
	}
#ifdef _OPENMP
#pragma omp single
{
#endif
	fftwf_free(fk);
	fftwf_free(a);
#ifdef _OPENMP
}
}
#endif
}

// Returns a floating-point interpolated value between 4 pixels in a float* array.
float bilinear(float *x, int w, int h, float i, float j) {
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

// Returns a floating-point interpolated value between 4 pixels in a WORD* array.
// No scaling is done so the value is between 0.f and 65535.f
float bilinear_ushort(WORD *x, int w, int h, float i, float j) {
	int ii = (int) i;
	int jj = (int) j;
	WORD (*xx)[w] = (void*)x;
	float a = (float) xx[min(max(jj, 0), h - 1)][min(max(ii, 0), w - 1)];
	float b = (float) xx[min(max(jj + 1, 0), h - 1)][min(max(ii, 0), w - 1)];
	float c = (float) xx[min(max(jj, 0), h - 1)][min(max(ii + 1, 0), w - 1)];
	float d = (float) xx[min(max(jj + 1, 0), h - 1)][min(max(ii + 1, 0), w - 1)];
	float xoff = i-ii;
	float yoff = j-jj;
	return (a * (1 - xoff) * (1 - yoff)) + (b * (1 - xoff) * yoff) + (c * xoff * (1 - yoff)) + (d * xoff * yoff);
}

void shrink(float *out, float *in, int outw, int outh, int inw, int inh, float scale, float sigma) {
#ifdef _OPENMP
	int available_threads = com.max_thread - omp_get_num_threads();
#endif
	if (scale == -2) {
		g_assert(2*outw == inw && 2 * outh == inh);
		float (*y)[outw] = (void*)out;
		float (*x)[inw] = (void*)in;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static, 16) collapse(2) num_threads(available_threads) if (available_threads > 1)
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
		memcpy(gaussian, in, inw * inh * sizeof(float));
	}
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static,16) collapse(2) num_threads(available_threads) if (available_threads > 1)
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

#ifdef __has_builtin
#if __has_builtin(__builtin_clz)
#define HAS_BUILTIN_CLZ 1
#else
#define HAS_BUILTIN_CLZ 0
#endif
#else
#define HAS_BUILTIN_CLZ 0
#endif

float half_to_float(const uint16_t val) {
	// Extract the sign from the bits
	const uint32_t sign = (uint32_t)(val & 0x8000) << 16;
	// Extract the exponent from the bits
	const uint8_t exp16 = (val & 0x7c00) >> 10;
	// Extract the fraction from the bits
	uint16_t frac16 = val & 0x3ff;

	uint32_t exp32 = 0;
	if (exp16 == 0x1f) {
		exp32 = 0xff; // Infinity or NaN
	} else if (exp16 != 0) {
		exp32 = exp16 + 112; // Normal numbers
	} else {
		if (frac16 == 0) {
			// Zero
			return *(float*)&sign;
		}

		// Denormal number: normalize it
		uint8_t offset;
		#if HAS_BUILTIN_CLZ
		offset = __builtin_clz(frac16) - 21; // Offset correction for 32-bit __builtin_clz()
		#else
		// Portable alternative: Count leading zeros manually
		offset = 0;
		while ((frac16 & 0x400) == 0) {
			frac16 <<= 1;
			++offset;
		}
		#endif
		frac16 &= 0x3ff; // Mask off the implicit leading bit
		exp32 = 113 - offset;
	}

	uint32_t frac32 = (uint32_t)frac16 << 13;

	// Compose the final FP32 binary representation
	uint32_t bits = 0;
	bits |= sign;
	bits |= (exp32 << 23);
	bits |= frac32;

	return *(float *)&bits;
}
