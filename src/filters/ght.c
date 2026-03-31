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
#include "ght.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "algos/statistics.h"
#include "algos/colors.h"
#include "core/siril_log.h"

// For clarity when referring to HSL layers
enum {
	HLAYER, SLAYER, LLAYER
};

int GHTsetup(ght_compute_params *c, float B, float D, float LP, float SP, float HP, int stretchtype) {
	if (D == 0.0f || stretchtype == STRETCH_LINEAR) {
		c->qlp = c->q0 = c->qwp = c->q1 = c->q = c->b1 = c->a2 = c->b2 = c->c2 = c->d2 = c->a3 = c->b3 = c->c3 = c->d3 = c->a4 = c->b4 = 0.0f;
		return 0;
	}
	else if (stretchtype == STRETCH_PAYNE_NORMAL) {
		if (B == -1.0f) {
			//B = -B;
			c->qlp = -1.0f*log1p(D*(SP - LP));
			c->q0 = c->qlp - D * LP / (1.0f + D * (SP - LP));
			c->qwp = log1p(D * (HP - SP));
			c->q1 = c->qwp + D * (1.0f - HP) / (1.0f + D * (HP - SP));
			c->q = 1.0f / (c->q1 - c->q0);
			c->b1 = (1.0f + D * (SP - LP)) / (D * c->q);
			c->a2 = (-c->q0) * c->q;
			c->b2 = -c->q;
			c->c2 = 1.0f + D * SP;
			c->d2 = -D;
			c->a3 = (-c->q0) * c->q;
			c->b3 = c->q;
			c->c3 = 1.0f - D * SP;
			c->d3 = D;
			c->a4 = (c->qwp - c->q0 - D * HP / (1.0f + D * (HP - SP))) * c->q;
			c->b4 = c->q * D / (1.0f + D * (HP - SP));
		} else if (B < 0.0f) {
			B = -B;
			c->qlp = (1.0f - pow((1.0f + D * B * (SP - LP)), (B - 1.0f) / B)) / (B - 1.0f);
			c->q0 = c->qlp - D * LP * (pow((1.0f + D * B * (SP - LP)), -1.0f / B));
			c->qwp = (pow((1.0f + D * B * (HP - SP)), (B - 1.0f) / B) - 1.0f) / (B - 1.0f);
			c->q1 = c->qwp + D * (1.0f - HP) * (pow((1.0f + D * B * (HP - SP)), -1.0f / B));
			c->q = 1.0f / (c->q1 - c->q0);
			c->b1 = D * pow(1.0f + D * B * (SP - LP), -1.0f / B) *c->q;
			c->a2 = (1.0f/(B-1.0f)-c->q0) * c->q;
			c->b2 = -c->q/(B-1.0f);
			c->c2 = 1.0f + D * B * SP;
			c->d2 = -D * B;
			c->e2 = (B - 1.0f)/B;
			c->a3 = (-1.0f/(B-1.0f) - c->q0) *c->q;
			c->b3 = c->q/(B-1.0f);
			c->c3 = 1.0f - D * B * SP;
			c->d3 = D * B;
			c->e3 = (B - 1.0f) / B;
			c->a4 = (c->qwp - c->q0 - D * HP * pow((1.0f + D * B * (HP - SP)), -1.0f / B)) * c->q;
			c->b4 = D * pow((1.0f + D * B * (HP - SP)), -1.0f / B) * c->q;
		} else if (B == 0.0f) {
			c->qlp = exp(-D * (SP - LP));
			c->q0 = c->qlp - D * LP * exp(-D*(SP - LP));
			c->qwp = 2.0f - exp(-D * (HP -SP));
			c->q1 = c->qwp + D * (1.0f - HP) * exp (-D * (HP - SP));
			c->q = 1.0f / (c->q1 - c->q0);
			c->a1 = 0.0f;
			c->b1 = D * exp (-D * (SP - LP)) * c->q;
			c->a2 = -c->q0 * c->q;
			c->b2 = c->q;
			c->c2 = -D * SP;
			c->d2 = D;
			c->a3 = (2.0f - c->q0) * c->q;
			c->b3 = -c->q;
			c->c3 = D * SP;
			c->d3 = -D;
			c->a4 = (c->qwp - c->q0 - D * HP * exp(-D * (HP - SP))) * c->q;
			c->b4 = D * exp(-D * (HP - SP)) * c->q;
		} else if (B > 0.0f) {
			c->qlp = pow((1.0f + D * B * (SP - LP)), -1.0f/B);
			c->q0 = c->qlp - D * LP * pow((1 + D * B * (SP - LP)), -(1.0f + B) / B);
			c->qwp = 2.0f - pow(1.0f + D * B * (HP - SP), -1.0f / B);
			c->q1 = c->qwp + D * (1.0f - HP) * pow((1.0f + D * B * (HP - SP)), -(1.0f + B) / B);
			c->q = 1.0f / (c->q1 - c->q0);
			c->b1 = D * pow((1.0f + D * B * (SP - LP)), -(1.0f+B)/B) * c->q;
			c->a2 = -c->q0 * c->q;
			c->b2 = c->q;
			c->c2 = 1.0f + D * B * SP;
			c->d2 = -D * B;
			c->e2 = -1.0f / B;
			c->a3 = (2.0f - c->q0) * c->q;
			c->b3 = -c->q;
			c->c3 = 1.0f - D * B * SP;
			c->d3 = D * B;
			c->e3 = -1.0f / B;
			c->a4 = (c->qwp - c->q0 - D * HP * pow((1.0f + D * B * (HP - SP)), -(B + 1.0f) / B)) * c->q;
			c->b4 = (D * pow((1.0f + D * B * (HP - SP)), -(B + 1.0f) / B)) * c->q;
		}
	} else if (stretchtype == STRETCH_PAYNE_INVERSE) {
		if (B == -1.0f) {
			//B = -B;
			c->qlp = -1.0f * log1p(D * (SP - LP));
			c->q0 = c->qlp - D * LP / (1.0f + D * (SP - LP));
			c->qwp = log1p(D * (HP - SP));
			c->q1 = c->qwp + D * (1.0f - HP) / (1.0f + D * (HP - SP));
			c->q = 1.0f / (c->q1 - c->q0);
			c->LPT = (c->qlp-c->q0)*c->q;
			c->SPT = c->q0*c->q;
			c->HPT = (c->qwp-c->q0)*c->q;
			c->b1 = (1.0f + D * (SP - LP)) / (D * c->q);
			c->a2 = (1.0f + D * SP) / D;
			c->b2 = -1.0f / D;
			c->c2 = - c->q0;
			c->d2 = - 1.0f/c->q;
			c->a3 = - (1.0f - D * SP) / D;
			c->b3 = 1.0f / D;
			c->c3 = c->q0;
			c->d3 = 1.0f / c->q;
			c->a4 = HP + (c->q0-c->qwp) * (1+D*(HP-SP))/D;
			c->b4 = (1.0f + D * (HP - SP) )/(c->q * D) ;
		} else if (B < 0.0f) {
			B = -B;
			c->qlp = (1.0f - pow((1.0f + D * B * (SP - LP)), (B - 1.0f) / B)) / (B - 1.0f);
			c->q0 = c->qlp - D * LP * (pow((1.0f + D * B * (SP - LP)), -1.0f / B));
			c->qwp = (pow((1.0f + D * B * (HP - SP)), (B - 1.0f) / B) - 1.0f) / (B - 1.0f);
			c->q1 = c->qwp + D * (1.0f - HP) * (pow((1.0f + D * B * (HP - SP)), -1.0f / B));
			c->q = 1.0f / (c->q1 - c->q0);
			c->LPT = (c->qlp-c->q0)*c->q;
			c->SPT = -c->q0*c->q;
			c->HPT = (c->qwp-c->q0)*c->q;
			c->b1 = pow(1.0f + D * B * (SP - LP), 1.0f / B) / (c->q * D);
			c->a2 = (1.0f + D * B * SP) / (D * B);
			c->b2 = -1.0f / (D * B);
			c->c2 = -c->q0*(B-1.0f) + 1.0f;
			c->d2 = (1.0f - B) / c->q;
			c->e2 = B / (B - 1.0f);
			c->a3 = (D * B * SP - 1.0f) / (D * B);
			c->b3 = 1.0f / (D * B);
			c->c3 = 1.0f + c->q0 * (B - 1);
			c->d3 = (B - 1.0f) / c->q;
			c->e3 = B / (B - 1.0f);
			c->a4 = (c->q0-c->qwp)/(D * pow((1.0f + D * B * (HP - SP)), -1.0f / B))+HP;
			c->b4 = 1.0f / (D * pow((1.0f + D * B * (HP - SP)), -1.0f / B) * c->q) ;
		} else if (B == 0.0f) {
			c->qlp = exp(-D * (SP - LP));
			c->q0 = c->qlp - D * LP * exp(-D*(SP - LP));
			c->qwp = 2.0f - exp(-D * (HP -SP));
			c->q1 = c->qwp + D * (1.0f - HP) * exp (-D * (HP - SP));
			c->q = 1.0f / (c->q1 - c->q0);
			c->LPT = (c->qlp-c->q0)*c->q;
			c->SPT = (1.0f-c->q0)*c->q;
			c->HPT = (c->qwp-c->q0)*c->q;
			c->a1 = 0.0f;
			c->b1 = 1.0f/(D * exp(-D * (SP - LP)) * c->q);
			c->a2 = SP;
			c->b2 = 1.0f / D;
			c->c2 = c->q0;
			c->d2 = 1.0f/c->q;
			c->a3 = SP;
			c->b3 = -1.0f / D;
			c->c3 = (2.0f - c->q0);
			c->d3 = -1.0f / c->q;
			c->a4 = (c->q0 - c->qwp)/(D * exp(-D * (HP - SP))) + HP;
			c->b4 = 1.0f/(D * exp(-D * (HP - SP)) * c->q);
		} else if (B > 0.0f) {
			c->qlp = pow(( 1.0f + D * B * (SP - LP)), -1.0f/B);
			c->q0 = c->qlp - D * LP * pow((1.0f + D * B * (SP - LP)), -(1.0f + B) / B);
			c->qwp = 2.0f - pow(1.0f + D * B * (HP - SP), -1.0f / B);
			c->q1 = c->qwp + D * (1.0f - HP) * pow((1.0f + D * B * (HP - SP)), -(1.0f + B) / B);
			c->q = 1.0f / (c->q1 - c->q0);
			c->LPT = (c->qlp-c->q0)*c->q;
			c->SPT = (1.0f-c->q0)*c->q;
			c->HPT = (c->qwp-c->q0)*c->q;
			c->b1 = 1/(D * pow((1.0f + D * B * (SP - LP)), -(1.0f + B) / B) * c->q);
			c->a2 = 1.0f / (D * B) + SP;
			c->b2 = -1.0f/(D * B);
			c->c2 = c->q0;
			c->d2 = 1.0f / c->q;
			c->e2 = -B;
			c->a3 = -1.0f / (D * B) + SP;
			c->b3 = 1.0f / (D * B);
			c->c3 = (2.0f - c->q0);
			c->d3 = -1.0f / c->q;
			c->e3 = -B;
			c->a4 = (c->q0-c->qwp)/(D * pow((1.0f + D * B * (HP - SP)), -(B + 1.0f) / B))+HP;
			c->b4 = 1.0f/((D * pow((1.0f + D * B * (HP - SP)), -(B + 1.0f) / B)) * c->q);
		}
	} else if (stretchtype == STRETCH_ASINH) {
		c->qlp = -log(D * (SP - LP) + pow((D * D * (SP - LP) * (SP - LP) + 1.0f), 0.5f));
		c->q0 = c->qlp - LP * D * pow((D * D * (SP - LP) * (SP - LP) + 1.0f), -0.5f);
		c->qwp = log(D * (HP - SP) + pow((D * D * (HP - SP) * (HP - SP) + 1.0f), 0.5f));
		c->q1 = c->qwp + (1.0f - HP) * D * pow((D * D * (HP - SP) * (HP - SP) +1.0f), -0.5f);
		c->q = 1.0f / (c->q1 - c->q0);
		c->a1 = 0.0f;
		c->b1 = D * pow((D * D * (SP - LP) * (SP - LP) + 1.0f),-0.5f)*c->q;
		c->a2 = -c->q0 * c->q;
		c->b2 = -c->q;
		c->c2 = -D;
		c->d2 = D * D;
		c->e2 = SP;
		c->a3 = -c->q0 * c->q;
		c->b3 = c->q;
		c->c3 = D;
		c->d3 = D * D;
		c->e3 = SP;
		c->a4 = (c->qwp - HP * D * pow((D * D * (HP - SP) * (HP - SP) + 1.0f), -0.5f) - c->q0) * c->q;
		c->b4 = D * pow((D * D * (HP - SP) * (HP - SP) + 1.0f), -0.5f) * c->q;
	}
	else if (stretchtype == STRETCH_INVASINH) {
		c->qlp = -log(D * (SP - LP) + pow((D * D * (SP - LP) * (SP - LP) + 1.0f), 0.5f));
		c->q0 = c->qlp - LP * D * pow((D * D * (SP - LP) * (SP - LP) + 1.0f), -0.5f);
		c->qwp = log(D * (HP - SP) + pow((D * D * (HP - SP) * (HP - SP) + 1.0f), 0.5f));
		c->q1 = c->qwp + (1.0f - HP) * D * pow((D * D * (HP - SP) * (HP - SP) +1.0f), -0.5f);
		c->q = 1.0f / (c->q1 - c->q0);
		c->a1 = 0.0f;
		c->b1 = D * pow((D * D * (SP - LP) * (SP - LP) + 1.0f),-0.5f)*c->q;
		c->a2 = -c->q0 * c->q;
		c->b2 = -c->q;
		c->c2 = -D;
		c->d2 = D * D;
		c->e2 = SP;
		c->a3 = -c->q0 * c->q;
		c->b3 = c->q;
		c->c3 = D;
		c->d3 = D * D;
		c->e3 = SP;
		c->a4 = (c->qwp - HP * D * pow((D * D * (HP - SP) * (HP - SP) + 1.0f), -0.5f) - c->q0) * c->q;
		c->b4 = D * pow((D * D * (HP - SP) * (HP - SP) + 1.0f), -0.5f) * c->q;
		c->LPT = c->a1 + c->b1 * LP;
		c->SPT = c->a2 + c->b2 * logf(c->c2 * (SP - c->e2) + sqrt(c->d2 * (SP - c->e2) * (SP - c->e2) + 1.0f));
		c->HPT = c->a4 + c->b4 * HP;
	}
	return 0;
}

float GHT(float in, float B, float D, float LP, float SP, float HP, float BP, int stretchtype, const ght_compute_params *c) {
    float out;
    float ex, val, res1, res2;

    BP = (stretchtype == STRETCH_LINEAR) ? BP : 0.0f;

    if (stretchtype == STRETCH_PAYNE_NORMAL || stretchtype == STRETCH_PAYNE_INVERSE) {
        in = fmaxf(0.0f, (in - BP) / (1.0f - BP));

        if (D == 0.0f) {
            return in;
        }

        if (stretchtype == STRETCH_PAYNE_NORMAL) {
            if (B == -1.0f) {
                res1 = c->a2 + c->b2 * logf(c->c2 + c->d2 * in);
                res2 = c->a3 + c->b3 * logf(c->c3 + c->d3 * in);
            } else if (B < 0.0f || B > 0.0f) {
                res1 = c->a2 + c->b2 * powf(c->c2 + c->d2 * in, c->e2);
                res2 = c->a3 + c->b3 * powf(c->c3 + c->d3 * in, c->e3);
            } else {
                res1 = c->a2 + c->b2 * expf(c->c2 + c->d2 * in);
                res2 = c->a3 + c->b3 * expf(c->c3 + c->d3 * in);
            }
            out = (in < LP) ? c->b1 * in : (in < SP) ? res1 : (in < HP) ? res2 : c->a4 + c->b4 * in;
        } else {
            if (B == -1.0f) {
                res1 = c->a2 + c->b2 * expf(c->c2 + c->d2 * in);
                res2 = c->a3 + c->b3 * expf(c->c3 + c->d3 * in);
            } else if (B < 0.0f || B > 0.0f) {
                res1 = c->a2 + c->b2 * powf(c->c2 + c->d2 * in, c->e2);
                res2 = c->a3 + c->b3 * powf(c->c3 + c->d3 * in, c->e3);
            } else {
                res1 = c->a2 + c->b2 * logf(c->c2 + c->d2 * in);
                res2 = c->a3 + c->b3 * logf(c->c3 + c->d3 * in);
            }
            out = (in < c->LPT) ? c->b1 * in : (in < c->SPT) ? res1 : (in < c->HPT) ? res2 : c->a4 + c->b4 * in;
        }
    } else if (stretchtype == STRETCH_LINEAR) {
        out = fmaxf(0.0f, (in - BP) / (1.0f - BP));
    } else if (stretchtype == STRETCH_ASINH) {
        in = fmaxf(0.0f, (in - BP) / (1.0f - BP));

        if (D == 0.0f) {
            out = in;
        } else {
            val = c->c2 * (in - c->e2) + sqrtf(c->d2 * (in - c->e2) * (in - c->e2) + 1.0f);
            res1 = c->a2 + c->b2 * logf(val);
            val = c->c3 * (in - c->e3) + sqrtf(c->d3 * (in - c->e3) * (in - c->e3) + 1.0f);
            res2 = c->a3 + c->b3 * logf(val);

            out = (in < LP) ? c->a1 + c->b1 * in : (in < SP) ? res1 : (in < HP) ? res2 : c->a4 + c->b4 * in;
        }
    } else {
        if (D == 0.0f) {
            out = in;
        } else {
            ex = expf((c->a2 - in) / c->b2);
            res1 = c->e2 - (ex - (1.0f / ex)) / (2.0f * c->c2);
            ex = expf((c->a3 - in) / c->b3);
            res2 = c->e3 - (ex - (1.0f / ex)) / (2.0f * c->c3);

            out = (in < c->LPT) ? (in - c->a1) / c->b1 : (in < c->SPT) ? res1 : (in < c->HPT) ? res2 : (in - c->a4) / c->b4;
        }
    }

    return out;
}

float GHTp(float in, const ght_params *params, const ght_compute_params *compute_params) {
	return GHT(in, params->B, params->D, params->LP, params->SP, params->HP, params->BP, params->stretchtype, compute_params);
}

/* ── Variant-specific inline evaluators ────────────────────────────────────
 * stretchtype and B sign are resolved once at the call site so the compiler
 * sees a clean, branch-free inner loop for each pixel. All take an input x
 * that has already been clamped to [0,1] and BP-normalised.
 * ─────────────────────────────────────────────────────────────────────────*/

/* Returns 0 for B==-1 (log formula), 1 for B==0 (exp formula),
 * 2 for any other B (pow formula).  Matches the branch structure in GHT(). */
static int b_category(float B) {
	if (B == -1.0f) return 0;
	if (B ==  0.0f) return 1;
	return 2;
}

static inline float ght_payne_normal_log(float x, const ght_compute_params *c,
                                          float LP, float SP, float HP) {
	float res1 = c->a2 + c->b2 * logf(c->c2 + c->d2 * x);
	float res2 = c->a3 + c->b3 * logf(c->c3 + c->d3 * x);
	return (x < LP) ? c->b1 * x : (x < SP) ? res1 : (x < HP) ? res2 : c->a4 + c->b4 * x;
}

static inline float ght_payne_normal_exp(float x, const ght_compute_params *c,
                                          float LP, float SP, float HP) {
	float res1 = c->a2 + c->b2 * expf(c->c2 + c->d2 * x);
	float res2 = c->a3 + c->b3 * expf(c->c3 + c->d3 * x);
	return (x < LP) ? c->b1 * x : (x < SP) ? res1 : (x < HP) ? res2 : c->a4 + c->b4 * x;
}

static inline float ght_payne_normal_pow(float x, const ght_compute_params *c,
                                          float LP, float SP, float HP) {
	float res1 = c->a2 + c->b2 * powf(c->c2 + c->d2 * x, c->e2);
	float res2 = c->a3 + c->b3 * powf(c->c3 + c->d3 * x, c->e3);
	return (x < LP) ? c->b1 * x : (x < SP) ? res1 : (x < HP) ? res2 : c->a4 + c->b4 * x;
}

static inline float ght_payne_inverse_exp(float x, const ght_compute_params *c) {
	float res1 = c->a2 + c->b2 * expf(c->c2 + c->d2 * x);
	float res2 = c->a3 + c->b3 * expf(c->c3 + c->d3 * x);
	return (x < c->LPT) ? c->b1 * x : (x < c->SPT) ? res1 : (x < c->HPT) ? res2 : c->a4 + c->b4 * x;
}

static inline float ght_payne_inverse_log(float x, const ght_compute_params *c) {
	float res1 = c->a2 + c->b2 * logf(c->c2 + c->d2 * x);
	float res2 = c->a3 + c->b3 * logf(c->c3 + c->d3 * x);
	return (x < c->LPT) ? c->b1 * x : (x < c->SPT) ? res1 : (x < c->HPT) ? res2 : c->a4 + c->b4 * x;
}

static inline float ght_payne_inverse_pow(float x, const ght_compute_params *c) {
	float res1 = c->a2 + c->b2 * powf(c->c2 + c->d2 * x, c->e2);
	float res2 = c->a3 + c->b3 * powf(c->c3 + c->d3 * x, c->e3);
	return (x < c->LPT) ? c->b1 * x : (x < c->SPT) ? res1 : (x < c->HPT) ? res2 : c->a4 + c->b4 * x;
}

static inline float ght_asinh_eval(float x, const ght_compute_params *c,
                                    float LP, float SP, float HP) {
	float dx2 = x - c->e2;
	float res1 = c->a2 + c->b2 * logf(c->c2 * dx2 + sqrtf(c->d2 * dx2 * dx2 + 1.0f));
	float dx3 = x - c->e3;
	float res2 = c->a3 + c->b3 * logf(c->c3 * dx3 + sqrtf(c->d3 * dx3 * dx3 + 1.0f));
	return (x < LP) ? c->a1 + c->b1 * x : (x < SP) ? res1 : (x < HP) ? res2 : c->a4 + c->b4 * x;
}

static inline float ght_invasinh_eval(float x, const ght_compute_params *c) {
	float ex1 = expf((c->a2 - x) / c->b2);
	float res1 = c->e2 - (ex1 - 1.0f / ex1) / (2.0f * c->c2);
	float ex2 = expf((c->a3 - x) / c->b3);
	float res2 = c->e3 - (ex2 - 1.0f / ex2) / (2.0f * c->c3);
	return (x < c->LPT) ? (x - c->a1) / c->b1 : (x < c->SPT) ? res1 : (x < c->HPT) ? res2 : (x - c->a4) / c->b4;
}

/* Dispatch GHT evaluation using an already-resolved variant.
 * Avoids branching on stretchtype/B inside the per-pixel hot path. */
#define GHT_DISPATCH(x, c, stretchtype, bcat, LP, SP, HP) \
	((stretchtype) == STRETCH_PAYNE_NORMAL \
	    ? ((bcat) == 0 ? ght_payne_normal_log((x),(c),(LP),(SP),(HP)) \
	                   : (bcat) == 1 ? ght_payne_normal_exp((x),(c),(LP),(SP),(HP)) \
	                                 : ght_payne_normal_pow((x),(c),(LP),(SP),(HP))) \
	: (stretchtype) == STRETCH_PAYNE_INVERSE \
	    ? ((bcat) == 0 ? ght_payne_inverse_exp((x),(c)) \
	                   : (bcat) == 1 ? ght_payne_inverse_log((x),(c)) \
	                                 : ght_payne_inverse_pow((x),(c))) \
	: (stretchtype) == STRETCH_ASINH \
	    ? ght_asinh_eval((x),(c),(LP),(SP),(HP)) \
	    : ght_invasinh_eval((x),(c)))

void apply_linked_ght_to_fbuf_lum(float* fbuf, float* out, size_t layersize, size_t nchans, ght_params *params, gboolean multithreaded) {
	const clip_mode_t clip_mode = params->clip_mode;
	const gboolean do_channel[3] = { params->do_red, params->do_green, params->do_blue };
	int active_channels = 3;
	const float m_CB = 1.f;
	float* fpbuf[nchans];
	float* outp[nchans];
	for (size_t chan = 0; chan < nchans; chan++) {
		size_t ndata = chan * layersize;
		fpbuf[chan] = fbuf + ndata;
		outp[chan] = out + ndata;
		if (!do_channel[chan])
			active_channels--;
	}
	if (active_channels == 0) {
		siril_log_color_message(_("Error: no channels selected. Doing nothing.\n"), "red");
		return;
	}
	struct ght_compute_params compute_params = { 0 };
	// If only working with selected channels, human-weighted luminance no longer
	// makes sense so luminance-based stretches are forced to even weighting.
	const float invchans = 1.f / active_channels;
	if (!(do_channel[0] && do_channel[1] && do_channel[2]))
		if (params->payne_colourstretchmodel == COL_HUMANLUM)
			params->payne_colourstretchmodel = COL_EVENLUM;
	const float factor_red   = params->payne_colourstretchmodel == COL_EVENLUM ? invchans : 0.2126f;
	const float factor_green = params->payne_colourstretchmodel == COL_EVENLUM ? invchans : 0.7152f;
	const float factor_blue  = params->payne_colourstretchmodel == COL_EVENLUM ? invchans : 0.0722f;
	GHTsetup(&compute_params, params->B, params->D, params->LP, params->SP, params->HP, params->stretchtype);

	// Resolve the transform variant once; stretchtype and bcat never change
	// within a call, so the compiler can eliminate dead branches from each loop.
	const int stretchtype = params->stretchtype;
	const int bcat = b_category(params->B);
	const float LP = params->LP, SP = params->SP, HP = params->HP;
	float globalmax = -FLT_MAX;

	// Hoist clip_mode outside the pixel loop: each mode gets its own tight
	// parallel loop with no per-pixel switch overhead.
	switch (clip_mode) {
	case CLIP:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
		for (size_t i = 0; i < layersize; i++) {
			float f0 = fmaxf(0.f, fminf(1.f, fpbuf[0][i]));
			float f1 = fmaxf(0.f, fminf(1.f, fpbuf[1][i]));
			float f2 = fmaxf(0.f, fminf(1.f, fpbuf[2][i]));
			float fbar = do_channel[0] * factor_red * f0 + do_channel[1] * factor_green * f1 + do_channel[2] * factor_blue * f2;
			float sfbar = GHT_DISPATCH(fbar, &compute_params, stretchtype, bcat, LP, SP, HP);
			float stretch_factor = sfbar / fmaxf(fbar, FLT_MIN);
			outp[0][i] = do_channel[0] ? fmaxf(0.f, fminf(1.f, f0 * stretch_factor)) : f0;
			outp[1][i] = do_channel[1] ? fmaxf(0.f, fminf(1.f, f1 * stretch_factor)) : f1;
			outp[2][i] = do_channel[2] ? fmaxf(0.f, fminf(1.f, f2 * stretch_factor)) : f2;
		}
		break;

	case RGBBLEND:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
		for (size_t i = 0; i < layersize; i++) {
			float f0 = fmaxf(0.f, fminf(1.f, fpbuf[0][i]));
			float f1 = fmaxf(0.f, fminf(1.f, fpbuf[1][i]));
			float f2 = fmaxf(0.f, fminf(1.f, fpbuf[2][i]));
			float fbar = do_channel[0] * factor_red * f0 + do_channel[1] * factor_green * f1 + do_channel[2] * factor_blue * f2;
			float sfbar = GHT_DISPATCH(fbar, &compute_params, stretchtype, bcat, LP, SP, HP);
			float stretch_factor = sfbar / fmaxf(fbar, FLT_MIN);
			blend_data data = { .do_channel = do_channel };
			data.sf[0] = f0 * stretch_factor;
			data.sf[1] = f1 * stretch_factor;
			data.sf[2] = f2 * stretch_factor;
			data.tf[0] = do_channel[0] ? GHT_DISPATCH(f0, &compute_params, stretchtype, bcat, LP, SP, HP) : 0.f;
			data.tf[1] = do_channel[1] ? GHT_DISPATCH(f1, &compute_params, stretchtype, bcat, LP, SP, HP) : 0.f;
			data.tf[2] = do_channel[2] ? GHT_DISPATCH(f2, &compute_params, stretchtype, bcat, LP, SP, HP) : 0.f;
			rgbblend(&data, &data.sf[0], &data.sf[1], &data.sf[2], m_CB);
			outp[0][i] = do_channel[0] ? data.sf[0] : f0;
			outp[1][i] = do_channel[1] ? data.sf[1] : f1;
			outp[2][i] = do_channel[2] ? data.sf[2] : f2;
		}
		break;

	case RESCALE:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
		for (size_t i = 0; i < layersize; i++) {
			float f0 = fmaxf(0.f, fminf(1.f, fpbuf[0][i]));
			float f1 = fmaxf(0.f, fminf(1.f, fpbuf[1][i]));
			float f2 = fmaxf(0.f, fminf(1.f, fpbuf[2][i]));
			float fbar = do_channel[0] * factor_red * f0 + do_channel[1] * factor_green * f1 + do_channel[2] * factor_blue * f2;
			float sfbar = GHT_DISPATCH(fbar, &compute_params, stretchtype, bcat, LP, SP, HP);
			float stretch_factor = sfbar / fmaxf(fbar, FLT_MIN);
			float sf0 = f0 * stretch_factor;
			float sf1 = f1 * stretch_factor;
			float sf2 = f2 * stretch_factor;
			float maxval = fmaxf(sf0, fmaxf(sf1, sf2));
			if (maxval > 1.f) {
				float invmaxval = 1.f / maxval;
				sf0 *= invmaxval;
				sf1 *= invmaxval;
				sf2 *= invmaxval;
			}
			outp[0][i] = do_channel[0] ? sf0 : f0;
			outp[1][i] = do_channel[1] ? sf1 : f1;
			outp[2][i] = do_channel[2] ? sf2 : f2;
		}
		break;

	case RESCALEGLOBAL:
		// Pass 1: find global maximum.
		// reduction(max:globalmax) fixes the data race present in the original.
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(max:globalmax) if (multithreaded)
#endif
		for (size_t i = 0; i < layersize; i++) {
			float f0 = fmaxf(0.f, fminf(1.f, fpbuf[0][i]));
			float f1 = fmaxf(0.f, fminf(1.f, fpbuf[1][i]));
			float f2 = fmaxf(0.f, fminf(1.f, fpbuf[2][i]));
			float fbar = do_channel[0] * factor_red * f0 + do_channel[1] * factor_green * f1 + do_channel[2] * factor_blue * f2;
			float sfbar = GHT_DISPATCH(fbar, &compute_params, stretchtype, bcat, LP, SP, HP);
			float stretch_factor = sfbar / fmaxf(fbar, FLT_MIN);
			float maxval = fmaxf(f0 * stretch_factor, fmaxf(f1 * stretch_factor, f2 * stretch_factor));
			if (maxval > globalmax)
				globalmax = maxval;
		}
		// Pass 2: apply with global normalisation.
		{
			const float inv_globalmax = (globalmax > 0.f) ? 1.f / globalmax : 1.f;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
			for (size_t i = 0; i < layersize; i++) {
				float f0 = fmaxf(0.f, fpbuf[0][i]);
				float f1 = fmaxf(0.f, fpbuf[1][i]);
				float f2 = fmaxf(0.f, fpbuf[2][i]);
				float fbar = do_channel[0] * factor_red * f0 + do_channel[1] * factor_green * f1 + do_channel[2] * factor_blue * f2;
				float sfbar = GHT_DISPATCH(fbar, &compute_params, stretchtype, bcat, LP, SP, HP);
				float stretch_factor = sfbar / fmaxf(fbar, FLT_MIN) * inv_globalmax;
				outp[0][i] = do_channel[0] ? f0 * stretch_factor : f0;
				outp[1][i] = do_channel[1] ? f1 * stretch_factor : f1;
				outp[2][i] = do_channel[2] ? f2 * stretch_factor : f2;
			}
		}
		break;
	}
}

void apply_linked_ght_to_fbuf_indep(float* in, float* out, size_t layersize, size_t nchans, ght_params *params, gboolean multithreaded) {
	const gboolean do_channel[3] = { params->do_red, params->do_green, params->do_blue };
	int active_channels = 3;
	float *fpbuf[nchans], *fpout[nchans];
	for (size_t chan = 0; chan < nchans; chan++) {
		size_t ndata = chan * layersize;
		fpbuf[chan] = in + ndata;
		fpout[chan] = out + ndata;
		if (!do_channel[chan])
			active_channels--;
	}
	if (active_channels == 0) {
		siril_log_color_message(_("Error: no channels selected. Doing nothing.\n"), "red");
		return;
	}
	struct ght_compute_params compute_params;
	GHTsetup(&compute_params, params->B, params->D, params->LP, params->SP, params->HP, params->stretchtype);

	const int stretchtype = params->stretchtype;
	const int bcat = b_category(params->B);
	const float LP = params->LP, SP = params->SP, HP = params->HP;
	const float BP = params->BP;

	// D == 0 for non-linear stretches: identity (just clamp to [0,1]).
	if (params->D == 0.f && stretchtype != STRETCH_LINEAR) {
		for (size_t chan = 0; chan < nchans; chan++) {
			const float *src = fpbuf[chan];
			float *dst = fpout[chan];
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
			for (size_t i = 0; i < layersize; i++)
				dst[i] = fmaxf(0.f, fminf(1.f, src[i]));
		}
		return;
	}

	// do_channel[chan] is constant within the inner loop: hoist it to the
	// channel level so inactive channels skip the transform entirely.
	// The dispatch on stretchtype × bcat is also resolved once per channel,
	// giving the compiler a branch-free inner loop to vectorise.
	for (size_t chan = 0; chan < nchans; chan++) {
		const float *src = fpbuf[chan];
		float *dst = fpout[chan];

		if (!do_channel[chan]) {
			// Inactive channel: copy with clamping to match original behaviour.
			for (size_t i = 0; i < layersize; i++)
				dst[i] = fmaxf(0.f, fminf(1.f, src[i]));
			continue;
		}

		if (stretchtype == STRETCH_PAYNE_NORMAL) {
			if (bcat == 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
				for (size_t i = 0; i < layersize; i++) {
					float x = fmaxf(0.f, fminf(1.f, src[i]));
					dst[i] = ght_payne_normal_log(x, &compute_params, LP, SP, HP);
				}
			} else if (bcat == 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
				for (size_t i = 0; i < layersize; i++) {
					float x = fmaxf(0.f, fminf(1.f, src[i]));
					dst[i] = ght_payne_normal_exp(x, &compute_params, LP, SP, HP);
				}
			} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
				for (size_t i = 0; i < layersize; i++) {
					float x = fmaxf(0.f, fminf(1.f, src[i]));
					dst[i] = ght_payne_normal_pow(x, &compute_params, LP, SP, HP);
				}
			}
		} else if (stretchtype == STRETCH_PAYNE_INVERSE) {
			if (bcat == 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
				for (size_t i = 0; i < layersize; i++) {
					float x = fmaxf(0.f, fminf(1.f, src[i]));
					dst[i] = ght_payne_inverse_exp(x, &compute_params);
				}
			} else if (bcat == 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
				for (size_t i = 0; i < layersize; i++) {
					float x = fmaxf(0.f, fminf(1.f, src[i]));
					dst[i] = ght_payne_inverse_log(x, &compute_params);
				}
			} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
				for (size_t i = 0; i < layersize; i++) {
					float x = fmaxf(0.f, fminf(1.f, src[i]));
					dst[i] = ght_payne_inverse_pow(x, &compute_params);
				}
			}
		} else if (stretchtype == STRETCH_ASINH) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
			for (size_t i = 0; i < layersize; i++) {
				float x = fmaxf(0.f, fminf(1.f, src[i]));
				dst[i] = ght_asinh_eval(x, &compute_params, LP, SP, HP);
			}
		} else if (stretchtype == STRETCH_INVASINH) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
			for (size_t i = 0; i < layersize; i++) {
				float x = fmaxf(0.f, fminf(1.f, src[i]));
				dst[i] = ght_invasinh_eval(x, &compute_params);
			}
		} else { // STRETCH_LINEAR
			const float invscale = 1.0f / (1.0f - BP);
			const float offset   = BP * invscale;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
			for (size_t i = 0; i < layersize; i++)
				dst[i] = fmaxf(0.f, fminf(1.f, src[i]) * invscale - offset);
		}
	}
}

void apply_linked_ght_to_Wbuf_lum(WORD* buf, WORD* out, size_t layersize, size_t nchans, ght_params *params, gboolean multithreaded) {
	const clip_mode_t clip_mode = params->clip_mode;
	const gboolean do_channel[3] = { params->do_red, params->do_green, params->do_blue };
	const float invnorm = 1.0f / USHRT_MAX_SINGLE;
	int active_channels = 3;
	float m_CB = 1.f; // This could be set to 0.f to switch off RGBBlend
	WORD* Wpbuf[nchans];
	WORD* outp[nchans];
	for (size_t chan = 0; chan < nchans; chan++) {
		size_t ndata = chan * layersize;
		Wpbuf[chan] = buf + ndata;
		outp[chan] = out + ndata;
		if (!do_channel[chan])
			active_channels--;
	}
	if (active_channels == 0) {
		siril_log_color_message(_("Error: no channels selected. Doing nothing.\n"), "red");
		return;
	}
	struct ght_compute_params compute_params = { 0 };
	// Independent stretching of mono or RGB channels
	// If only working with selected channels, human-weighted luminance no longer makes sense so luminance
	// based stretches are forced to even weighting
	if (!(do_channel[0] && do_channel[1] && do_channel[2]))
		if (params->payne_colourstretchmodel == COL_HUMANLUM)
			params->payne_colourstretchmodel = COL_EVENLUM;
	float factor_red = params->payne_colourstretchmodel == COL_EVENLUM ? 1.f / active_channels : 0.2126f;
	float factor_green = params->payne_colourstretchmodel == COL_EVENLUM ? 1.f / active_channels : 0.7152f;
	float factor_blue = params->payne_colourstretchmodel == COL_EVENLUM ? 1.f / active_channels : 0.0722f;
	// Do calcs that can be done prior to the loop
	GHTsetup(&compute_params, params->B, params->D, params->LP, params->SP, params->HP, params->stretchtype);
	float globalmax = -FLT_MAX;

	// Set up a LUT
	WORD *lut = malloc((USHRT_MAX + 1) * sizeof(WORD));
#ifdef _OPENMP
	// This is only a small loop: 8 threads seems to be about as many as is worthwhile
	// because of the thread startup cost
	int threads = min(com.max_thread, 8);
#pragma omp parallel for num_threads(threads) schedule(static) if (multithreaded)
#endif
	for (int i = 0 ; i <= USHRT_MAX ; i++) {
		lut[i] = roundf_to_WORD(USHRT_MAX_SINGLE * GHTp(i * invnorm, params, &compute_params));
	}

	// Loop over pixels
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(max:globalmax) if (multithreaded)
#endif
	for (size_t i = 0 ; i < layersize ; i++) {
		float f[3];
		if (clip_mode == RGBBLEND) {
			for (size_t chan = 0; chan < 3 ; chan++)
				f[chan] = fmaxf(0.f, fminf(1.f, Wpbuf[chan][i] * invnorm));
		} else {
			for (size_t chan = 0; chan < 3 ; chan++)
				f[chan] = Wpbuf[chan][i] * invnorm;
		}
		WORD fbar = roundf_to_WORD(USHRT_MAX_SINGLE * (do_channel[0] * factor_red * f[0] + do_channel[1] * factor_green * f[1] + do_channel[2] * factor_blue * f[2]));
		WORD sfbar = lut[fbar];
		float stretch_factor = (fbar == 0) ? 0.f : (float) sfbar / (float) fbar;
		blend_data data = { .do_channel = do_channel };
		//Calculate the luminance and independent channel stretches for the pixel
		for (size_t chan = 0; chan < 3 ; chan++) {
			data.sf[chan] = f[chan] * stretch_factor;
		}
		if (clip_mode == RGBBLEND) {
			for (size_t chan = 0; chan < 3 ; chan++) {
				if (do_channel[chan]) {
					data.tf[chan] = invnorm * lut[roundf_to_WORD(USHRT_MAX_SINGLE * f[chan])];
				}
			}
		}
		float maxval;
		switch (clip_mode) {
			case CLIP:
				for (size_t chan = 0 ; chan < 3 ; chan++) {
					outp[chan][i] = (do_channel[chan]) ? roundf_to_WORD(USHRT_MAX_SINGLE * data.sf[chan]) : Wpbuf[chan][i];
				}
				break;
			case RGBBLEND:
				rgbblend(&data, &data.sf[0], &data.sf[1], &data.sf[2], m_CB);
				for (size_t chan = 0 ; chan < 3 ; chan++) {
					f[chan] = (do_channel[chan]) ? data.sf[chan] : f[chan];
					outp[chan][i] = roundf_to_WORD(USHRT_MAX_SINGLE * f[chan]);
				}
				break;
			case RESCALE:
				maxval = fmaxf(data.sf[0], fmaxf(data.sf[1], data.sf[2]));
				float invmaxval = 1.f / maxval;
				if (maxval > 1.f) {
					data.sf[0] *= invmaxval;
					data.sf[1] *= invmaxval;
					data.sf[2] *= invmaxval;
				}
				for (size_t chan = 0 ; chan < 3 ; chan++) {
					outp[chan][i] = (do_channel[chan]) ? roundf_to_WORD(USHRT_MAX_SINGLE * data.sf[chan]) : Wpbuf[chan][i];
				}
				break;
			case RESCALEGLOBAL:
				maxval = fmaxf(data.sf[0], fmaxf(data.sf[1], data.sf[2]));
				if (maxval > globalmax)
					globalmax = maxval;
				break;
		}
	}
	free(lut);
	if (clip_mode == RESCALEGLOBAL) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
		for (size_t i = 0 ; i < layersize ; i++) {
			float f[3];
			for (size_t chan = 0; chan < 3 ; chan++)
				f[chan] = Wpbuf[chan][i] * invnorm;
			float fbar = (int) do_channel[0] * factor_red * f[0] + (int) do_channel[1] * factor_green * f[1] + (int) do_channel[2] * factor_blue * f[2];
			float sfbar = GHTp(fbar, params, &compute_params);
			float stretch_factor = (fbar == 0) ? 0.f : sfbar / fbar;
			blend_data data = { .do_channel = do_channel };
			//Calculate the luminance and independent channel stretches for the pixel
			for (size_t chan = 0; chan < 3 ; chan++) {
				data.sf[chan] = f[chan] * stretch_factor;
			}
			for (size_t chan = 0 ; chan < 3 ; chan++) {
				data.sf[chan] /= globalmax;
				outp[chan][i] = (do_channel[chan]) ? roundf_to_WORD(USHRT_MAX_SINGLE * data.sf[chan]) : Wpbuf[chan][i];
			}
		}
	}
}

void apply_linked_ght_to_Wbuf_indep(WORD* in, WORD* out, size_t layersize, size_t nchans, ght_params *params, gboolean multithreaded) {
	const gboolean do_channel[3] = { params->do_red, params->do_green, params->do_blue };
	const float invnorm = 1.0f / USHRT_MAX_SINGLE;
	int active_channels = 3;
	WORD *Wpbuf[nchans], *Wpout[nchans];
	for (size_t chan = 0; chan < nchans; chan++) {
		size_t ndata = chan * layersize;
		Wpbuf[chan] = in + ndata;
		Wpout[chan] = out + ndata;
		if (!do_channel[chan])
			active_channels--;
	}
	if (active_channels == 0) {
		siril_log_color_message(_("Error: no channels selected. Doing nothing.\n"), "red");
		return;
	}
	struct ght_compute_params compute_params;

	// Do calcs that can be done prior to the loop
	GHTsetup(&compute_params, params->B, params->D, params->LP, params->SP, params->HP, params->stretchtype);

	// Set up a LUT
	WORD *lut = malloc((USHRT_MAX + 1) * sizeof(WORD));
#ifdef _OPENMP
	// This is only a small loop: 8 threads seems to be about as many as is worthwhile
	// because of the thread startup cost
	int threads = min(com.max_thread, 8);
#pragma omp parallel for num_threads(threads) schedule(static) if (multithreaded)
#endif
	for (int i = 0 ; i <= USHRT_MAX ; i++) {
		lut[i] = roundf_to_WORD(USHRT_MAX_SINGLE * GHTp(i * invnorm, params, &compute_params));
	}

	// Loop over pixels: independent stretching of mono or RGB channels
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) collapse(2) if (multithreaded)
#endif
	for (size_t chan=0; chan < nchans; chan++) {
		for (size_t i = 0; i < layersize; i++) {
			Wpout[chan][i] = do_channel[chan] ? lut[Wpbuf[chan][i]] : Wpbuf[chan][i];
		}
	}
	free(lut);
}

void apply_linked_ght_to_fits(fits *from, fits *to, ght_params *params, gboolean multithreaded) {
	g_assert(from);
	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	g_assert(from->type == to->type);
	size_t npixels = from->rx * from->ry;
	if (from->type == DATA_FLOAT) {
		if (from->naxes[2] == 3 && params->stretchtype != STRETCH_LINEAR && params->payne_colourstretchmodel != COL_INDEP) {
			apply_linked_ght_to_fbuf_lum(from->fdata, to->fdata, npixels, from->naxes[2], params, multithreaded);
		} else {
			apply_linked_ght_to_fbuf_indep(from->fdata, to->fdata, npixels, from->naxes[2], params, multithreaded);
		}
	} else if (from->type == DATA_USHORT) {
		if (from->naxes[2] == 3 && params->stretchtype != STRETCH_LINEAR && params->payne_colourstretchmodel != COL_INDEP) {
			apply_linked_ght_to_Wbuf_lum(from->data, to->data, npixels, from->naxes[2], params, multithreaded);
		} else {
			apply_linked_ght_to_Wbuf_indep(from->data, to->data, npixels, from->naxes[2], params, multithreaded);
		}
	}
	invalidate_stats_from_fit(to);
	return;
}

void apply_ght_to_fits_channel(fits *from, fits *to, int channel, ght_params *params, gboolean multithreaded) {
	g_assert(from);
	g_assert(from->naxes[2] <= 3 && from->naxes[2] > channel);
	g_assert(from->type == to->type);
	size_t npixels = from->rx * from->ry;
	if (from->type == DATA_FLOAT) {
		apply_linked_ght_to_fbuf_indep(from->fpdata[channel], to->fpdata[channel], npixels, 1, params, multithreaded);
	} else if (from->type == DATA_USHORT) {
		apply_linked_ght_to_Wbuf_indep(from->pdata[channel], to->pdata[channel], npixels, 1, params, multithreaded);
	}
	invalidate_stats_from_fit(to);
	return;
}

/* Can't avoid creation of a full-sized float buffer here owing to the need to convert to hsl and back */
void apply_sat_ght_to_fits(fits *fit, ght_params *params, gboolean multithreaded) {
	g_assert(fit);
	g_assert(fit->naxes[2] == 1 || fit->naxes[2] == 3);
	if (fit->naxes[2] == 1)
		return;
	size_t npixels = fit->rx * fit->ry;

	if (fit->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) if (multithreaded)
		{
#pragma omp for schedule(static)
#endif
			for (long i = 0; i < npixels; i++) {
				rgb_to_hslf(fit->fpdata[RLAYER][i], fit->fpdata[GLAYER][i], fit->fpdata[BLAYER][i], &fit->fpdata[HLAYER][i], &fit->fpdata[SLAYER][i], &fit->fpdata[LLAYER][i]);
			}

#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
			apply_linked_ght_to_fbuf_indep(fit->fpdata[SLAYER], fit->fpdata[SLAYER], npixels, 1, params, multithreaded);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (long i = 0; i < npixels; i++) {
				hsl_to_rgbf(fit->fpdata[HLAYER][i], fit->fpdata[SLAYER][i], fit->fpdata[LLAYER][i], &fit->fpdata[RLAYER][i], &fit->fpdata[GLAYER][i], &fit->fpdata[BLAYER][i]);
			}
#ifdef _OPENMP
		}
#endif
	} else if (fit->type == DATA_USHORT) {

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) if (multithreaded)
		{
#pragma omp for schedule(static)
#endif
			for (long i = 0; i < npixels; i++) {
				rgbw_to_hslw(fit->pdata[RLAYER][i], fit->pdata[GLAYER][i], fit->pdata[BLAYER][i], &fit->pdata[HLAYER][i], &fit->pdata[SLAYER][i], &fit->pdata[LLAYER][i]);
			}

#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
			apply_linked_ght_to_Wbuf_indep(fit->pdata[SLAYER], fit->pdata[SLAYER], npixels, 1, params, multithreaded);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (long i = 0; i < npixels; i++) {
				hslw_to_rgbw(fit->pdata[HLAYER][i], fit->pdata[SLAYER][i], fit->pdata[LLAYER][i], &fit->pdata[RLAYER][i], &fit->pdata[GLAYER][i], &fit->pdata[BLAYER][i]);
			}

#ifdef _OPENMP
		}
#endif
	}
	invalidate_stats_from_fit(fit);
	return;
}
