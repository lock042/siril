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
#ifndef SRC_CORE_ARITHM_H_
#define SRC_CORE_ARITHM_H_

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* operations on image data */
typedef enum {
	OPER_ADD,
	OPER_SUB,
	OPER_MUL,
	OPER_DIV
} image_operator;

typedef struct blend_data {
	float sf[3]; // Luminance stretched values
	float tf[3]; // Independently stretched values
	const gboolean *do_channel; // Whether or not to do each channel
} blend_data;

int soper(fits *a, float scalar, image_operator oper, gboolean conv_to_float);
int imoper(fits *a, fits *b, image_operator oper, gboolean allow_32bits);
int descreen(fits *a, fits *b, gboolean allow_32bits, int threads);
int screen(fits *a, fits *b, fits *result, gboolean allow_32bits, int threads);
int addmax(fits *a, fits *b);
int siril_fdiv(fits *a, fits *b, float scalar, gboolean allow_32bits);
int siril_ndiv(fits *a, fits *b);

int soper_unscaled_div_ushort_to_float(fits *a, int scalar);
void clip(fits *fit);
void clipneg(fits *fit);

// Scaling / blurring functions
void magnify(float *y, const float *x, int W, int H, int pd, int w, int h, float n);
void shrink(float *y, float *x, int outw, int outh, int inw, int inh, float scale, float sigma);
void gaussblur(float *y, float *x, int w, int h, float sigma);
float half_to_float(const uint16_t val);
int round_down_to_multiple(int value, int multiple);

// Static inline functions
// RGB clipping modes
static inline void rgbblend(blend_data *data, float* r, float* g, float* b, float m_CB) {
    // load channel components
    float sf0 = data->sf[0], sf1 = data->sf[1], sf2 = data->sf[2];
    float tf0 = data->tf[0], tf1 = data->tf[1], tf2 = data->tf[2];

    // maxima
#define FASTMAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
    float sfmax = FASTMAX(FASTMAX(sf0, sf1), sf2);
    float tfmax = FASTMAX(FASTMAX(tf0, tf1), tf2);
#undef FASTMAX

    // difference
    float d = sfmax - tfmax;

    // build masks as 0.0f / 1.0f without branching
    float mask_cond = (tfmax + m_CB * d > 1.0f) ? 1.0f : 0.0f;   // condition tfmax + m_CB*d > 1
    float mask_dnz  = (d != 0.0f)               ? 1.0f : 0.0f;   // d != 0

    // full_mask = mask_cond && mask_dnz  (still 0.0f or 1.0f)
    float full_mask = mask_cond * mask_dnz;

    // construct safe_d: if d!=0 then d else 1.0f  (done branchless)
    float safe_d = mask_dnz * d + (1.0f - mask_dnz) * 1.0f;

    // candidate and limited value
    float candidate = (1.0f - tfmax) / safe_d;       // safe even when original d==0 due to safe_d
    float limited   = fminf(m_CB, candidate);

    // final k = full_mask ? limited : m_CB  (blend)
    float k = full_mask * limited + (1.0f - full_mask) * m_CB;

    // precompute factor
    float one_minus_k = 1.0f - k;

    // channel masks (0 or 1) to decide whether to update each channel
    float mr = data->do_channel[0] ? 1.0f : 0.0f;
    float mg = data->do_channel[1] ? 1.0f : 0.0f;
    float mb = data->do_channel[2] ? 1.0f : 0.0f;

    // blend per-channel without branching:
    // if channel enabled: result = one_minus_k * tf + k * sf
    // else: keep old value (*r, *g, *b)
    *r = mr * (one_minus_k * tf0 + k * sf0) + (1.0f - mr) * (*r);
    *g = mg * (one_minus_k * tf1 + k * sf1) + (1.0f - mg) * (*g);
    *b = mb * (one_minus_k * tf2 + k * sf2) + (1.0f - mb) * (*b);
}
#ifdef __cplusplus
}
#endif

#endif /* SRC_CORE_ARITHM_H_ */
