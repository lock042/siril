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
	int n = g_list_length(points);
	g_assert(n >= 2);

	// Precalculate the slope between each pair of points
	GList *current = points;
	for (int i = 0; i < n - 1; i++) {
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

	/* Akima needs at least 3 points (the VeraLux reference uses linear
	 * interpolation for 2 points); a cubic spline through 2 points is
	 * exactly that straight line. */
	if (n < 3) {
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
		*h = 0.0f;
		return;
	}

	if (delta == 0.0f) {
		*h = 0.0f;
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
// LUMINANCE RANGE MASK
// ------------------------------------------------------------------------

/* Per-pixel luminance range mask, matching the VeraLux reference: a pure
 * sigmoid roll-off (k = 2.5) on the pixel luminance, with NO spatial blur
 * ("no spatial blur for maximum pixel-to-pixel fidelity").  The reference
 * recomputes the mask from the *current* intermediate image before each
 * channel is applied, which a pointwise formula reproduces exactly. */
static inline float lum_mask_value(float lum, const curve_channel_config *cfg) {
	if (cfg->feather < 1e-6f)
		return (lum < cfg->lum_min || lum > cfg->lum_max) ? 0.0f : 1.0f;

	float val = 1.0f;
	const float k_smooth = 2.5f;
	if (cfg->lum_min > 0.0f)
		val = fminf(val, sigmoid(k_smooth * (lum - cfg->lum_min) / cfg->feather));
	if (cfg->lum_max < 1.0f)
		val = fminf(val, sigmoid(k_smooth * (cfg->lum_max - lum) / cfg->feather));
	return val;
}

// Exported helpers so the GUI can build its L / C histograms and pipette
// readouts in the exact color spaces the transform operates in.
float curve_lab_L_norm(float r, float g, float b) {
	float L, A, B_val;
	rgb_to_lab(r, g, b, &L, &A, &B_val);
	L /= 100.0f;
	return fmaxf(0.0f, fminf(1.0f, L));
}

float curve_lab_chroma_norm(float r, float g, float b) {
	float L, A, B_val;
	rgb_to_lab(r, g, b, &L, &A, &B_val);
	float c = sqrtf(A * A + B_val * B_val) / 128.0f;
	return fmaxf(0.0f, fminf(1.0f, c));
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

// Runtime cache: one LUT per active channel
typedef struct {
	gboolean active;
	float lut[LUT_SIZE];
} channel_runtime_cache;

/* LUT lookup with linear interpolation between samples, equivalent to the
 * reference's np.interp() over the LUT domain. */
static inline float apply_lut(float v, const float *lut) {
	if (v <= 0.0f) return lut[0];
	if (v >= 1.0f) return lut[LUT_SIZE - 1];
	float pos = v * (float)(LUT_SIZE - 1);
	int idx = (int)pos;
	float frac = pos - (float)idx;
	return lut[idx] + (lut[idx + 1] - lut[idx]) * frac;
}

static inline float clamp01(float v) {
	return fmaxf(0.0f, fminf(1.0f, v));
}

/* Applies the full VeraLux curve pipeline.  Channel order, luminance-mask
 * source (the current intermediate pixel value, not the input image) and
 * the [0,1] clamp after every color space round trip all follow the
 * reference implementation: RGB/K -> R,G,B -> L (Lab) -> S (HSV) -> C (Lab). */
void apply_curve(fits *from, fits *to, struct curve_params *params, gboolean multithreaded) {
	g_assert(from->type == to->type);
	const size_t layersize = from->naxes[0] * from->naxes[1];
	int channels = from->naxes[2];

	// LUTs are 256 KB each: allocate the cache on the heap
	channel_runtime_cache *caches = calloc(CHAN_COUNT, sizeof(channel_runtime_cache));
	if (!caches) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (int c = 0; c < CHAN_COUNT; c++) {
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
		if (caches[c].active)
			generate_lut(cfg->points, params->algorithm, caches[c].lut);
	}

	#define MIX(orig, new_val, mask_val) ((orig) * (1.0f - (mask_val)) + (new_val) * (mask_val))
	// Mask of the channel being applied, from the current pixel state
	#define RANGE_MASK(chan) \
		(params->channels[chan].range_enabled \
			? lum_mask_value((r + g + b) / 3.0f, &params->channels[chan]) : 1.0f)

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

		r = clamp01(r);
		if (channels > 1) { g = clamp01(g); b = clamp01(b); }
		else { g = r; b = r; }

		// RGB Master
		if (caches[CHAN_RGB_K].active) {
			float m = RANGE_MASK(CHAN_RGB_K);
			float nr = apply_lut(r, caches[CHAN_RGB_K].lut);
			r = MIX(r, nr, m);
			if (channels > 1) {
				float ng = apply_lut(g, caches[CHAN_RGB_K].lut);
				float nb = apply_lut(b, caches[CHAN_RGB_K].lut);
				g = MIX(g, ng, m);
				b = MIX(b, nb, m);
			} else {
				g = r; b = r;
			}
		}

		// R, G, B — each mask sees the pixel as left by the previous channel
		if (caches[CHAN_R].active) {
			float m = RANGE_MASK(CHAN_R);
			r = MIX(r, apply_lut(r, caches[CHAN_R].lut), m);
			if (channels == 1) { g = r; b = r; }
		}
		if (channels > 1) {
			if (caches[CHAN_G].active) {
				float m = RANGE_MASK(CHAN_G);
				g = MIX(g, apply_lut(g, caches[CHAN_G].lut), m);
			}
			if (caches[CHAN_B].active) {
				float m = RANGE_MASK(CHAN_B);
				b = MIX(b, apply_lut(b, caches[CHAN_B].lut), m);
			}
		}

		// L: CIE L* adjusted in Lab space (safe for mono: r=g=b stays gray)
		if (caches[CHAN_L].active) {
			float m = RANGE_MASK(CHAN_L);
			float L, A, B_val;
			rgb_to_lab(r, g, b, &L, &A, &B_val);
			float l_norm = L / 100.0f;
			float l_new = MIX(l_norm, apply_lut(l_norm, caches[CHAN_L].lut), m);
			L = clamp01(l_new) * 100.0f;
			lab_to_rgb(L, A, B_val, &r, &g, &b);
			r = clamp01(r); g = clamp01(g); b = clamp01(b);
		}

		// S: HSV saturation
		if (channels > 1 && caches[CHAN_S].active) {
			float m = RANGE_MASK(CHAN_S);
			float h, s, v;
			rgb_to_hsv(r, g, b, &h, &s, &v);
			s = clamp01(MIX(s, apply_lut(s, caches[CHAN_S].lut), m));
			hsv_to_rgb(h, s, v, &r, &g, &b);
			r = clamp01(r); g = clamp01(g); b = clamp01(b);
		}

		// C: Lab chroma, scaled by C_new / C to preserve hue
		if (channels > 1 && caches[CHAN_C].active) {
			float m = RANGE_MASK(CHAN_C);
			float L, A, B_val;
			rgb_to_lab(r, g, b, &L, &A, &B_val);
			float C = sqrtf(A * A + B_val * B_val);
			float c_norm = fminf(C / 128.0f, 1.0f);
			float c_new_norm = MIX(c_norm, apply_lut(c_norm, caches[CHAN_C].lut), m);
			float factor = (C > 1e-5f) ? (c_new_norm * 128.0f / C) : 1.0f;
			A *= factor;
			B_val *= factor;
			lab_to_rgb(L, A, B_val, &r, &g, &b);
			r = clamp01(r); g = clamp01(g); b = clamp01(b);
		}

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

	#undef MIX
	#undef RANGE_MASK

	if (params->clipped_count) {
		params->clipped_count[0] = local_clipped_low;
		params->clipped_count[1] = local_clipped_high;
	}

	free(caches);
}

/* "Show Mask" preview: replace the image with the luminance range mask of
 * one channel, computed from the input pixels (white = full application,
 * black = excluded), as the VeraLux reference does. */
static void render_lum_mask(fits *from, fits *to, const curve_channel_config *cfg, gboolean multithreaded) {
	g_assert(from->type == to->type);
	const size_t layersize = from->naxes[0] * from->naxes[1];
	int channels = from->naxes[2];

	float norm = 1.0f;
	float inv_norm = 1.0f;
	if (from->type == DATA_USHORT) {
		norm = (float) get_normalized_value(from);
		inv_norm = 1.f / norm;
	}

	#ifdef _OPENMP
	#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithreaded)
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
		float lum = (clamp01(r) + clamp01(g) + clamp01(b)) / 3.0f;
		float m = lum_mask_value(lum, cfg);

		if (to->type == DATA_USHORT) {
			WORD w = roundf_to_WORD(m * norm);
			for (int c = 0; c < channels; c++)
				to->pdata[c][j] = w;
		} else {
			for (int c = 0; c < channels; c++)
				to->fpdata[c][j] = m;
		}
	}
}

// Hooks...
int curve_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct curve_params *params = (struct curve_params *)args->user;
	if (!params) return 1;

	params->fit = fit;
	if (params->show_mask && params->for_preview
			&& params->channels[params->target_channel].range_enabled)
		render_lum_mask(fit, fit, &params->channels[params->target_channel], TRUE);
	else
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