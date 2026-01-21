/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#include <glib.h>
#include "curve_transform.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/image_format_fits.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*****************************************************************************
 *      MATH IMPLEMENTATIONS (Spline Fits)
 ****************************************************************************/

void linear_fit(GList *points, double *slopes) {
	if (!points || g_list_length(points) < 2) return;
	GList *current = points;
	for (int i = 0; i < g_list_length(points) - 1; i++) {
		point *point1 = (point *) current->data;
		point *point2 = (point *) current->next->data;
		double dx = point2->x - point1->x;
		slopes[i] = (dx != 0) ? (point2->y - point1->y) / dx : 0;
		current = current->next;
	}
}

float linear_interpolate(float x, GList *points, double *slopes) {
	if (!points) return x;
	if (x > ((point *) g_list_last(points)->data)->x)
		return ((point *) g_list_last(points)->data)->y;
	else if (x < ((point *) points->data)->x)
		return ((point *) points->data)->y;

	x = fmax(0, fmin(1, x));
	GList *point1 = points;
	int point_index = 0;
	while (point1->next && point1->next->next != NULL && ((point *) point1->next->data)->x < x) {
		point_index++;
		point1 = point1->next;
	}
	float x1 = ((point *) point1->data)->x;
	float y1 = ((point *) point1->data)->y;
	return y1 + slopes[point_index] * (x - x1);
}

void cubic_spline_fit(GList *points, cubic_spline_data *cspline_data) {
	if (!points || g_list_length(points) < 2) return;
	cspline_data->n = g_list_length(points);
	if (cspline_data->n > MAX_POINTS) cspline_data->n = MAX_POINTS;

	double h[MAX_POINTS], alpha[MAX_POINTS], l[MAX_POINTS], mu[MAX_POINTS], z[MAX_POINTS];

	GList *current = points;
	for (int i = 0; i < cspline_data->n; i++) {
		point *p = (point *) current->data;
		cspline_data->x_values[i] = p->x;
		cspline_data->y_values[i] = p->y;
		current = current->next;
	}

	l[0] = 1; mu[0] = 0; z[0] = 0;

	for (int i = 0; i < cspline_data->n - 1; i++)
		h[i] = cspline_data->x_values[i + 1] - cspline_data->x_values[i];

	for (int i = 1; i < cspline_data->n - 1; i++)
		alpha[i] = (3 / h[i]) * (cspline_data->y_values[i + 1] - cspline_data->y_values[i]) -
				   (3 / h[i - 1]) * (cspline_data->y_values[i] - cspline_data->y_values[i - 1]);

	for (int i = 1; i < cspline_data->n - 1; i++) {
		l[i] = 2 * (cspline_data->x_values[i + 1] - cspline_data->x_values[i - 1]) - h[i - 1] * mu[i - 1];
		mu[i] = h[i] / l[i];
		z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	l[cspline_data->n - 1] = 1; z[cspline_data->n - 1] = 0; cspline_data->c[cspline_data->n - 1] = 0;

	for (int j = cspline_data->n - 2; j >= 0; j--) {
		cspline_data->c[j] = z[j] - mu[j] * cspline_data->c[j + 1];
		cspline_data->b[j] = (cspline_data->y_values[j + 1] - cspline_data->y_values[j]) / h[j] -
							 h[j] * (cspline_data->c[j + 1] + 2 * cspline_data->c[j]) / 3.0;
		cspline_data->d[j] = (cspline_data->c[j + 1] - cspline_data->c[j]) / (3 * h[j]);
	}
}

float cubic_spline_interpolate(float x, cubic_spline_data *cspline_data) {
	if (x >= cspline_data->x_values[cspline_data->n - 1]) return cspline_data->y_values[cspline_data->n - 1];
	if (x <= cspline_data->x_values[0]) return cspline_data->y_values[0];

	int i = 0;
	while (i < cspline_data->n - 1 && x > cspline_data->x_values[i + 1]) i++;

	double diff = x - cspline_data->x_values[i];
	double val = cspline_data->y_values[i] + cspline_data->b[i] * diff +
	             cspline_data->c[i] * diff * diff + cspline_data->d[i] * diff * diff * diff;
	return (float)val;
}

void akima_spline_fit(GList *points, akima_spline_data *data) {
	int n = g_list_length(points);
	if (!points || n < 2) return;

	if (n < 5) {
		cubic_spline_data cb;
		cubic_spline_fit(points, &cb);
		memcpy(data->x_values, cb.x_values, sizeof(double)*MAX_POINTS);
		memcpy(data->y_values, cb.y_values, sizeof(double)*MAX_POINTS);
		memcpy(data->b, cb.b, sizeof(double)*MAX_POINTS);
		memcpy(data->c, cb.c, sizeof(double)*MAX_POINTS);
		memcpy(data->d, cb.d, sizeof(double)*MAX_POINTS);
		data->n = n;
		return;
	}

	if (n > MAX_POINTS) n = MAX_POINTS;
	data->n = n;

	GList *curr = points;
	for (int i = 0; i < n; i++) {
		point *p = (point *) curr->data;
		data->x_values[i] = p->x;
		data->y_values[i] = p->y;
		curr = curr->next;
	}

	double m[MAX_POINTS + 4];
	for (int i = 0; i < n - 1; i++) {
		m[i + 2] = (data->y_values[i + 1] - data->y_values[i]) / (data->x_values[i + 1] - data->x_values[i]);
	}

	m[1] = 2 * m[2] - m[3];
	m[0] = 2 * m[1] - m[2];
	m[n + 1] = 2 * m[n] - m[n - 1];
	m[n + 2] = 2 * m[n + 1] - m[n];

	double t[MAX_POINTS];
	for (int i = 0; i < n; i++) {
		double m_minus_2 = m[i];
		double m_minus_1 = m[i+1];
		double m_0       = m[i+2];
		double m_plus_1  = m[i+3];

		double w1 = fabs(m_plus_1 - m_0);
		double w2 = fabs(m_minus_1 - m_minus_2);

		if (w1 + w2 == 0.0) {
			t[i] = 0.5 * (m_minus_1 + m_0);
		} else {
			t[i] = (w1 * m_minus_1 + w2 * m_0) / (w1 + w2);
		}
	}

	for (int i = 0; i < n - 1; i++) {
		double dx = data->x_values[i+1] - data->x_values[i];
		if (dx == 0) dx = 1e-6;
		data->b[i] = t[i];
		data->c[i] = (3 * m[i+2] - 2 * t[i] - t[i+1]) / dx;
		data->d[i] = (t[i] + t[i+1] - 2 * m[i+2]) / (dx * dx);
	}
}

float akima_spline_interpolate(float x, akima_spline_data *data) {
	if (data->n == 0) return x;
	if (x >= data->x_values[data->n - 1]) return data->y_values[data->n - 1];
	if (x <= data->x_values[0]) return data->y_values[0];

	int i = 0;
	while (i < data->n - 1 && x > data->x_values[i + 1]) i++;

	double diff = x - data->x_values[i];
	double val = data->y_values[i] + data->b[i] * diff +
	             data->c[i] * diff * diff + data->d[i] * diff * diff * diff;
	return (float)val;
}

/*****************************************************************************
 *      COLOR SPACE HELPERS
 ****************************************************************************/

static inline float f_func(float t) {
	return (t > 0.008856f) ? powf(t, 1.0f/3.0f) : (7.787f * t + 16.0f/116.0f);
}

static inline float inv_f(float t) {
	return (t > 0.206893f) ? (t * t * t) : ((t - 16.0f/116.0f) / 7.787f);
}

static inline void rgb_to_hsv(float r, float g, float b, float *h, float *s, float *v) {
	float min_v = fminf(r, fminf(g, b));
	float max_v = fmaxf(r, fmaxf(g, b));
	float delta = max_v - min_v;

	*v = max_v;
	if (max_v != 0.0f)
		*s = delta / max_v;
	else {
		*s = 0.0f;
		*h = -1.0f;
		return;
	}

	if (r == max_v)
		*h = (g - b) / delta;
	else if (g == max_v)
		*h = 2.0f + (b - r) / delta;
	else
		*h = 4.0f + (r - g) / delta;

	*h *= 60.0f;
	if (*h < 0.0f) *h += 360.0f;
}

static inline void hsv_to_rgb(float h, float s, float v, float *r, float *g, float *b) {
	if (s == 0.0f) {
		*r = *g = *b = v;
		return;
	}
	if (h >= 360.0f) h = 0.0f;
	h /= 60.0f;
	int i = (int)h;
	float f = h - i;
	float p = v * (1.0f - s);
	float q = v * (1.0f - s * f);
	float t = v * (1.0f - s * (1.0f - f));

	switch (i) {
		case 0: *r = v; *g = t; *b = p; break;
		case 1: *r = q; *g = v; *b = p; break;
		case 2: *r = p; *g = v; *b = t; break;
		case 3: *r = p; *g = q; *b = v; break;
		case 4: *r = t; *g = p; *b = v; break;
		default: *r = v; *g = p; *b = q; break;
	}
}

static inline void rgb_to_lab(float r, float g, float b, float *L, float *A, float *B_val) {
	// D65 Standard Observer
	float X = 0.4124564f * r + 0.3575761f * g + 0.1804375f * b;
	float Y = 0.2126729f * r + 0.7151522f * g + 0.0721750f * b;
	float Z = 0.0193339f * r + 0.1191920f * g + 0.9503041f * b;

	X /= 0.95047f;
	Z /= 1.08883f;

	float fx = f_func(X);
	float fy = f_func(Y);
	float fz = f_func(Z);

	*L = (116.0f * fy) - 16.0f;
	*A = 500.0f * (fx - fy);
	*B_val = 200.0f * (fy - fz);
}

static inline void lab_to_rgb(float L, float A, float B_val, float *r, float *g, float *b) {
	float fy = (L + 16.0f) / 116.0f;
	float fx = A / 500.0f + fy;
	float fz = fy - B_val / 200.0f;

	float X = inv_f(fx) * 0.95047f;
	float Y = inv_f(fy);
	float Z = inv_f(fz) * 1.08883f;

	*r = 3.2404542f * X - 1.5371385f * Y - 0.4985314f * Z;
	*g = -0.9692660f * X + 1.8760108f * Y + 0.0415560f * Z;
	*b = 0.0556434f * X - 0.2040259f * Y + 1.0572252f * Z;
}

/*****************************************************************************
 *      MEMORY & HELPERS
 ****************************************************************************/

struct curve_params *new_curve_params() {
	struct curve_params *params = calloc(1, sizeof(struct curve_params));
	if (params) {
		params->destroy_fn = free_curve_params;
		params->target_channel = CHAN_RGB_K;
		for(int i=0; i<CHAN_COUNT; i++) {
			params->channels[i].points = NULL;
			params->channels[i].active = FALSE;
			params->channels[i].range_enabled = FALSE;
			params->channels[i].lum_min = 0.0f;
			params->channels[i].lum_max = 1.0f;
			params->channels[i].feather = 0.25f;
		}
	}
	return params;
}

void free_curve_params(void *ptr) {
	struct curve_params *params = (struct curve_params *)ptr;
	if (!params) return;
	for (int i=0; i<CHAN_COUNT; i++) {
		if (params->channels[i].points) {
			g_list_free_full(params->channels[i].points, g_free);
		}
	}
	free(ptr);
}

static inline float sigmoid(float x) {
	return 1.0f / (1.0f + expf(-x));
}

// ------------------------------------------------------------------------
// MASK GENERATION (Gaussian Blur)
// ------------------------------------------------------------------------

static void gaussian_blur_buffer(float *buffer, int w, int h, float sigma) {
	if (sigma <= 0.1f) return;

	int r = (int)ceilf(sigma * 3.0f);
	int ksize = 2 * r + 1;
	float *kernel = malloc(ksize * sizeof(float));
	float sum = 0.0f;
	for (int i = 0; i < ksize; i++) {
		float x = i - r;
		kernel[i] = expf(-(x * x) / (2.0f * sigma * sigma));
		sum += kernel[i];
	}
	for (int i = 0; i < ksize; i++) kernel[i] /= sum;

	float *temp = malloc(w * h * sizeof(float));
	if (!temp || !kernel) {
		free(kernel); free(temp);
		return;
	}

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			float val = 0.0f;
			for (int k = -r; k <= r; k++) {
				int px = x + k;
				px = (px < 0) ? -px : (px >= w) ? 2 * w - px - 2 : px;
				if (px < 0) px = 0;
				if (px >= w) px = w - 1;

				val += buffer[y * w + px] * kernel[k + r];
			}
			temp[y * w + x] = val;
		}
	}

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int x = 0; x < w; x++) {
		for (int y = 0; y < h; y++) {
			float val = 0.0f;
			for (int k = -r; k <= r; k++) {
				int py = y + k;
				py = (py < 0) ? -py : (py >= h) ? 2 * h - py - 2 : py;
				if (py < 0) py = 0;
				if (py >= h) py = h - 1;

				val += temp[py * w + x] * kernel[k + r];
			}
			buffer[y * w + x] = val;
		}
	}

	free(kernel);
	free(temp);
}

static float* generate_lum_mask(fits *fit, float l_min, float l_max, float feather) {
	size_t npixels = fit->rx * fit->ry;
	float *mask = malloc(npixels * sizeof(float));
	if (!mask) return NULL;

	float k_smooth = 2.5f;
	int channels = fit->naxes[2];
	float norm = 1.0f;
	if (fit->type == DATA_USHORT) norm = 1.0f / (float)get_normalized_value(fit);

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
	#endif
	for (size_t i = 0; i < npixels; i++) {
		float r, g, b;
		if (fit->type == DATA_USHORT) {
			r = fit->pdata[0][i] * norm;
			g = (channels > 1) ? fit->pdata[1][i] * norm : r;
			b = (channels > 2) ? fit->pdata[2][i] * norm : r;
		} else {
			r = fit->fpdata[0][i];
			g = (channels > 1) ? fit->fpdata[1][i] : r;
			b = (channels > 2) ? fit->fpdata[2][i] : r;
		}

		float lum = (r + g + b) / 3.0f;
		float val = 1.0f;

		if (feather < 1e-6f) {
			if (lum < l_min || lum > l_max) val = 0.0f;
		} else {
			if (l_min > 0.0f) {
				float dist = (lum - l_min) / feather;
				val = fminf(val, sigmoid(k_smooth * dist));
			}
			if (l_max < 1.0f) {
				float dist = (l_max - lum) / feather;
				val = fminf(val, sigmoid(k_smooth * dist));
			}
		}
		mask[i] = val;
	}

	if (feather > 1e-6f) {
		float blur_sigma = feather * 150.0f;
		if (blur_sigma > 100.0f) blur_sigma = 100.0f;
		if (blur_sigma < 1.0f) blur_sigma = 1.0f;
		gaussian_blur_buffer(mask, fit->rx, fit->ry, blur_sigma);
	}

	return mask;
}

// ------------------------------------------------------------------------
// LUT GENERATION
// ------------------------------------------------------------------------

static void generate_lut(GList *points, enum curve_algorithm algo, float *lut_out) {
	cubic_spline_data cubic;
	akima_spline_data akima;
	double slopes[MAX_POINTS];

	if (algo == CUBIC_SPLINE) cubic_spline_fit(points, &cubic);
	else if (algo == AKIMA_SPLINE) akima_spline_fit(points, &akima);
	else linear_fit(points, slopes);

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < LUT_SIZE; i++) {
		float x = (float)i / (float)(LUT_SIZE - 1);
		float y = x;

		if (algo == CUBIC_SPLINE) y = cubic_spline_interpolate(x, &cubic);
		else if (algo == AKIMA_SPLINE) y = akima_spline_interpolate(x, &akima);
		else y = linear_interpolate(x, points, slopes);

		if (y < 0.0f) y = 0.0f;
		if (y > 1.0f) y = 1.0f;
		lut_out[i] = y;
	}
}

// ------------------------------------------------------------------------
// CORE PIPELINE
// ------------------------------------------------------------------------

void apply_curve(fits *from, fits *to, struct curve_params *params, gboolean multithreaded) {
	g_assert(from->type == to->type);
	const size_t layersize = from->naxes[0] * from->naxes[1];
	int channels = from->naxes[2];

	// *** FIX: ALLOCATE CACHE ON HEAP (Stack Overflow Prevention) ***
	channel_runtime_cache *caches = calloc(CHAN_COUNT, sizeof(channel_runtime_cache));
	if (!caches) return; // Should log error

	for(int c=0; c<CHAN_COUNT; c++) {
		curve_channel_config *cfg = &params->channels[c];

		gboolean is_identity = TRUE;
		if (cfg->points && g_list_length(cfg->points) >= 2) {
			point *p1 = (point*)cfg->points->data;
			point *p2 = (point*)cfg->points->next->data;
			if (fabs(p1->x) > 1e-5 || fabs(p1->y) > 1e-5 ||
			    fabs(p2->x - 1.0) > 1e-5 || fabs(p2->y - 1.0) > 1e-5 ||
				g_list_length(cfg->points) > 2) {
				is_identity = FALSE;
			}
		}

		caches[c].active = !is_identity;

		if (caches[c].active) {
			generate_lut(cfg->points, params->algorithm, caches[c].lut);
			if (cfg->range_enabled) {
				caches[c].use_mask = TRUE;
				caches[c].mask_buffer = generate_lum_mask(from, cfg->lum_min, cfg->lum_max, cfg->feather);
			}
		}
	}

	#define APPLY_LUT(val, cache_ptr) \
		({ \
			float _v = (val); \
			int _idx = (int)(_v * (LUT_SIZE - 1)); \
			if (_idx < 0) _idx = 0; \
			if (_idx >= LUT_SIZE) _idx = LUT_SIZE - 1; \
			(cache_ptr)->lut[_idx]; \
		})

	#define MIX(orig, new_val, mask_val) ((orig) * (1.0f - (mask_val)) + (new_val) * (mask_val))

	float norm = 1.0f;
	float inv_norm = 1.0f;
	if (from->type == DATA_USHORT) {
		norm = (float) get_normalized_value(from);
		inv_norm = 1.f / norm;
	}

	long local_clipped_low = 0;
    long local_clipped_high = 0;

	#ifdef _OPENMP
    #pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded) reduction(+:local_clipped_low, local_clipped_high)
    #endif
    for (size_t j = 0; j < layersize; j++) {
        float r, g, b;

		if (from->type == DATA_USHORT) {
			r = from->pdata[0][j] * inv_norm;
			g = (channels > 1) ? from->pdata[1][j] * inv_norm : r;
			b = (channels > 2) ? from->pdata[2][j] * inv_norm : r;
		} else {
			r = from->fpdata[0][j];
			g = (channels > 1) ? from->fpdata[1][j] : r;
			b = (channels > 2) ? from->fpdata[2][j] : r;
		}

		// RGB Master
		if (caches[CHAN_RGB_K].active) {
			float nr = APPLY_LUT(r, &caches[CHAN_RGB_K]);
			float ng = (channels > 1) ? APPLY_LUT(g, &caches[CHAN_RGB_K]) : nr;
			float nb = (channels > 2) ? APPLY_LUT(b, &caches[CHAN_RGB_K]) : nr;

			if (caches[CHAN_RGB_K].use_mask && caches[CHAN_RGB_K].mask_buffer) {
				float m = caches[CHAN_RGB_K].mask_buffer[j];
				r = MIX(r, nr, m);
				if (channels > 1) { g = MIX(g, ng, m); b = MIX(b, nb, m); }
			} else {
				r = nr; if (channels > 1) { g = ng; b = nb; }
			}
		}

		// R, G, B
		if (caches[CHAN_R].active) {
			float nr = APPLY_LUT(r, &caches[CHAN_R]);
			if (caches[CHAN_R].use_mask && caches[CHAN_R].mask_buffer) r = MIX(r, nr, caches[CHAN_R].mask_buffer[j]);
			else r = nr;
		}
		if (channels > 1) {
			if (caches[CHAN_G].active) {
				float ng = APPLY_LUT(g, &caches[CHAN_G]);
				if (caches[CHAN_G].use_mask && caches[CHAN_G].mask_buffer) g = MIX(g, ng, caches[CHAN_G].mask_buffer[j]);
				else g = ng;
			}
			if (caches[CHAN_B].active) {
				float nb = APPLY_LUT(b, &caches[CHAN_B]);
				if (caches[CHAN_B].use_mask && caches[CHAN_B].mask_buffer) b = MIX(b, nb, caches[CHAN_B].mask_buffer[j]);
				else b = nb;
			}
		}

		// Lab
		if (channels > 1 && (caches[CHAN_L].active || caches[CHAN_C].active)) {
			float L, A, B_val;
			rgb_to_lab(r, g, b, &L, &A, &B_val);

			if (caches[CHAN_L].active) {
				float l_norm = L / 100.0f;
				float l_new = APPLY_LUT(l_norm, &caches[CHAN_L]);
				if (caches[CHAN_L].use_mask && caches[CHAN_L].mask_buffer) l_new = MIX(l_norm, l_new, caches[CHAN_L].mask_buffer[j]);
				L = l_new * 100.0f;
			}

			if (caches[CHAN_C].active) {
				float C = sqrtf(A*A + B_val*B_val);
				float c_norm = C / 128.0f;
				if (c_norm > 1.0f) c_norm = 1.0f;

				float c_new_norm = APPLY_LUT(c_norm, &caches[CHAN_C]);
				if (caches[CHAN_C].use_mask && caches[CHAN_C].mask_buffer) c_new_norm = MIX(c_norm, c_new_norm, caches[CHAN_C].mask_buffer[j]);

				float C_new = c_new_norm * 128.0f;
				float factor = (C > 1e-5f) ? (C_new / C) : 1.0f;
				A *= factor;
				B_val *= factor;
			}
			lab_to_rgb(L, A, B_val, &r, &g, &b);
		}

		// HSV
		if (channels > 1 && caches[CHAN_S].active) {
			float h, s, v;
			rgb_to_hsv(r, g, b, &h, &s, &v);

			float s_new = APPLY_LUT(s, &caches[CHAN_S]);
			if (caches[CHAN_S].use_mask && caches[CHAN_S].mask_buffer) s_new = MIX(s, s_new, caches[CHAN_S].mask_buffer[j]);

			s = s_new;
			hsv_to_rgb(h, s, v, &r, &g, &b);
		}

		r = fmaxf(0.0f, fminf(1.0f, r));
        g = fmaxf(0.0f, fminf(1.0f, g));
        b = fmaxf(0.0f, fminf(1.0f, b));

        if (r <= 0.0f) local_clipped_low++; else if (r >= 1.0f) local_clipped_high++;

        if (channels > 1) {
            if (g <= 0.0f) local_clipped_low++; else if (g >= 1.0f) local_clipped_high++;
            if (b <= 0.0f) local_clipped_low++; else if (b >= 1.0f) local_clipped_high++;
        }

		if (to->type == DATA_USHORT) {
			to->pdata[0][j] = roundf_to_WORD(r * norm);
			if (channels > 1) {
				to->pdata[1][j] = roundf_to_WORD(g * norm);
				to->pdata[2][j] = roundf_to_WORD(b * norm);
			}
		} else {
			to->fpdata[0][j] = r;
			if (channels > 1) {
				to->fpdata[1][j] = g;
				to->fpdata[2][j] = b;
			}
		}
	}

	if (params->clipped_count) {
        params->clipped_count[0] = local_clipped_low;
        params->clipped_count[1] = local_clipped_high;
    }

	// Cleanup Heap
	for(int c=0; c<CHAN_COUNT; c++) {
		if (caches[c].mask_buffer) free(caches[c].mask_buffer);
	}
	free(caches); // IMPORTANT
}

// Hooks...
int curve_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct curve_params *params = (struct curve_params *)args->user;
	if (!params) return 1;

	params->fit = fit;
	apply_curve(fit, fit, params, TRUE);
	return 0;
}

gchar *curves_log_hook(gpointer p, log_hook_detail detail) {
	struct curve_params *params = (struct curve_params *) p;
	gchar *algo_name = (params->algorithm == LINEAR) ? "linear" : ((params->algorithm == AKIMA_SPLINE) ? "Akima" : "cubic spline");

	if (detail == SUMMARY) {
		return g_strdup_printf(_("%s curve transformation (VeraLux Engine)"), algo_name);
	} else {
		return g_strdup_printf(_("Applying %s curve transformation with Range Masking..."), algo_name);
	}
}