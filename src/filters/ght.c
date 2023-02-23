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

#include <glib.h>
#include "ght.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "core/siril_log.h"

float GHT(float in, float B, float D, float LP, float SP, float HP, float BP, int stretchtype, ght_compute_params *c) {
	float out;
	if (stretchtype != STRETCH_LINEAR) {
		BP = 0.0f;
	}
	if (stretchtype == STRETCH_PAYNE_NORMAL || stretchtype == STRETCH_PAYNE_INVERSE) {
		in = max (0.0f, (in - BP)/(1.0f - BP));
		if (D == 0.0f) {
			out = in;
		} else if (stretchtype == STRETCH_PAYNE_NORMAL) {
			if (B == -1.0f) {
				if (in < LP) {
					out = c->b1 * in;
				} else if (in < SP) {
					out = c->a2 + c->b2 * logf(c->c2 + c->d2 * in);
				} else if (in < HP) {
					out = c->a3 + c->b3 * logf(c->c3 + c->d3 * in);
				} else {
					out = c->a4 + c->b4 * in;
				}
			} else if (B < 0.0f) {
				if (in < LP) {
					out = c->b1 * in;
				} else if (in < SP) {
					out = c->a2 + c->b2 * pow((c->c2 + c->d2 * in), c->e2);
				} else if (in < HP) {
					out = c->a3 + c->b3 * pow((c->c3 + c->d3 * in), c->e3);
				} else {
					out = c->a4 + c->b4 * in;
				}
			} else if (B == 0.0f) {
				if (in < LP) {
					out = c->a1 + c->b1 * in;
				} else if (in < SP) {
					out = c->a2 + c->b2 * exp(c->c2 + c->d2 * in);
				} else if (in < HP) {
					out = c->a3 + c->b3 * exp(c->c3 + c->d3 * in);
				} else {
					out = c->a4 + c->b4 * in;
				}
			} else /*if (B > 0)*/ {
				if (in < LP) {
					out = c->b1 * in;
				} else if (in < SP) {
					out = c->a2 + c->b2 * pow((c->c2 + c->d2 * in), c->e2);
				} else if (in < HP) {
					out = c->a3 + c->b3 * pow((c->c3 + c->d3 * in), c->e3);
				} else {
					out = c->a4 + c->b4 * in;
				}
			}
		} else /*if (stretchtype == STRETCH_PAYNE_INVERSE)*/ {
			if (B == -1.0f) {
				if (in < c->LPT) {
					out = c->b1 * in;
				} else if (in < c->SPT) {
					out = c->a2 + c->b2 * exp(c->c2 + c->d2 * in);
				} else if (in < c->HPT) {
					out = c->a3 + c->b3 * exp(c->c3 + c->d3 * in);
				} else {
					out = c->a4 + c->b4 * in;
				}
			} else if (B < 0.0f) {
				if (in < c->LPT) {
					out = c->b1 * in;
				} else if (in < c->SPT) {
					out = c->a2 + c->b2 * pow((c->c2 + c->d2 * in), c->e2);
				} else if (in < c->HPT) {
					out = c->a3 + c->b3 * pow((c->c3 + c->d3 * in), c->e3);
				} else {
					out = c->a4 + c->b4 * in;
				}
			} else if (B == 0.0f) {
				if (in < c->LPT) {
					out = c->a1 + c->b1 * in;
				} else if (in < c->SPT) {
					out = c->a2 + c->b2 * logf(c->c2 + c->d2 * in);
				} else if (in < c->HPT) {
					out = c->a3 + c->b3 * logf(c->c3 + c->d3 * in);
				} else {
					out = c->a4 + c->b4 * in;
				}
			} else /* if (B > 0) */{
				if (in < c->LPT) {
					out = c->b1 * in;
				} else if (in < c->SPT) {
					out = c->a2 + c->b2 * pow((c->c2 + c->d2 * in), c->e2);
				} else if (in < c->HPT) {
					out = c->a3 + c->b3 * pow((c->c3 + c->d3 * in), c->e3);
				} else {
					out = c->a4 + c->b4 * in;
				}
			}
		}
	} else if (stretchtype == STRETCH_LINEAR) {
		out = max(0.0f, (in - BP) / (1.0f - BP));
	} else if (stretchtype == STRETCH_ASINH) {
		in = max(0.0f, (in - BP) / (1.0f - BP));
		if (D == 0.0f) {
			out = in;
		} else if (in < LP) {
			out = c->a1 + c->b1 * in;
		} else if (in < SP) {
			out = c->a2 + c->b2 * logf(c->c2 * (in - c->e2) + sqrt(c->d2 * (in - c->e2) * (in - c->e2) + 1.0f));
		} else if (in < HP) {
			out = c->a3 + c->b3 * logf(c->c3 * (in - c->e3) + sqrt(c->d3 * (in - c->e3) * (in - c->e3) + 1.0f));
		} else {
			out = c->a4 + c->b4 * in;
		}
	} else /* if (stretchtype == STRETCH_INVASINH)*/ {
		if (D == 0.0f) {
			out = in;
		} else if (in < c->LPT) {
			out = (in - c->a1) / c->b1;
		} else if (in < c->SPT) {
			float ex = exp((c->a2 - in) / c->b2);
			out = c->e2 - (ex - (1.0f / ex)) / (2.0f * c->c2);
		} else if (in < c->HPT) {
			float ex = exp((c->a3 - in) / c->b3);
			out = c->e3 - (ex - (1.0f / ex)) / (2.0f * c->c3);
		} else {
			out = (in - c->a4) / c->b4;
		}
	}
	return out;
}

float GHTp(float in, ght_params *params, ght_compute_params *compute_params) {
	return GHT(in, params->B, params->D, params->LP, params->SP, params->HP, params->BP, params->stretchtype, compute_params);
}

int GHTsetup(ght_compute_params *c, float B, float D, float LP, float SP, float HP, int stretchtype) {
	if (D == 0.0f || stretchtype == STRETCH_LINEAR) {
		c->qlp = c->q0 = c->qwp = c->q1 = c->q = c->b1 = c->a2 = c->b2 = c->c2 = c->d2 = c->a3 = c->b3 = c->c3 = c->d3 = c->a4 = c->b4 = 0.0f;
		return 0;
	}
	else if (stretchtype == STRETCH_PAYNE_NORMAL) {
		if (B == -1.0f) {
			//B = -B;
			c->qlp = -1.0f*log1pf(D*(SP - LP));
			c->q0 = c->qlp - D * LP / (1.0f + D * (SP - LP));
			c->qwp = log1pf(D * (HP - SP));
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
			c->qlp = -1.0f * log1pf(D * (SP - LP));
			c->q0 = c->qlp - D * LP / (1.0f + D * (SP - LP));
			c->qwp = log1pf(D * (HP - SP));
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
		c->qlp = -logf(D * (SP - LP) + pow((D * D * (SP - LP) * (SP - LP) + 1.0f), 0.5f));
		c->q0 = c->qlp - LP * D * pow((D * D * (SP - LP) * (SP - LP) + 1.0f), -0.5f);
		c->qwp = logf(D * (HP - SP) + pow((D * D * (HP - SP) * (HP - SP) + 1.0f), 0.5f));
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
		c->qlp = -logf(D * (SP - LP) + pow((D * D * (SP - LP) * (SP - LP) + 1.0f), 0.5f));
		c->q0 = c->qlp - LP * D * pow((D * D * (SP - LP) * (SP - LP) + 1.0f), -0.5f);
		c->qwp = logf(D * (HP - SP) + pow((D * D * (HP - SP) * (HP - SP) + 1.0f), 0.5f));
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

void apply_linked_ght_to_fits(fits *from, fits *to, ght_params *params, gboolean multithreaded) {
	const gboolean do_channel[3] = {params->do_red, params->do_green, params->do_blue};
	int active_channels = 3;
	for (size_t i=0;i<3;i++)
		if (!do_channel[i])
			active_channels--;
	if (active_channels == 0) {
		siril_log_color_message(_("Error: no channels selected. Doing nothing.\n"), "red");
		return;
	}
	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t layersize = from->naxes[0] * from->naxes[1];
	g_assert(from->type == to->type);
	struct ght_compute_params compute_params;
	float factor_red = 0.2126f;
	float factor_green = 0.7152f;
	float factor_blue = 0.0722f;
	// If only working with selected channels, colour mode is forced to independent
	if (!(do_channel[0] && do_channel[1] && do_channel[2]))
		if (params->payne_colourstretchmodel == COL_HUMANLUM)
			params->payne_colourstretchmodel = COL_EVENLUM;
	if (params->payne_colourstretchmodel == COL_EVENLUM) {
		factor_red = 1.0f/active_channels;
		factor_green = 1.0f/active_channels;
		factor_blue = 1.0f/active_channels;
	}
	// Do calcs that can be done prior to the loop
	GHTsetup(&compute_params, params->B, params->D, params->LP, params->SP, params->HP, params->stretchtype);
	if (from->type == DATA_USHORT) {
		float norm = get_normalized_value(from);
		float invnorm = 1.0f / norm;
		if (from->naxes[2] == 3 && params->payne_colourstretchmodel != COL_INDEP) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
			for (size_t i = 0 ; i < layersize ; i++) {
				float L[3] = { 0.0f };
				for (size_t chan = 0; chan < 3 ; chan++)
					L[chan] = (float) from->pdata[chan][i] * invnorm;
				float x = (int) do_channel[0] * factor_red * L[0] + (int) do_channel[1] * factor_green * L[1] + (int) do_channel[2] * factor_blue * L[2];
				float z = GHTp(x, params, &compute_params);
				for (size_t chan = 0; chan < 3 ; chan++) {
					if (do_channel[chan])
						to->pdata[chan][i] = (x == 0.0f) ? 0 : roundf_to_WORD(norm * min(1.0f, max(0.0f, L[chan] * (z / x))));
				}
			}
			for (size_t chan = 0 ; chan < 3; chan++) {
				if (!do_channel[chan])
					memcpy(to->pdata[chan], from->pdata[chan], layersize * sizeof(WORD));
			}
		} else {
			for (size_t chan = 0; chan < from->naxes[2]; chan++) {
				if (do_channel[chan]) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
					for (size_t i = 0; i < layersize; i++) {
						float x = from->pdata[chan][i] * invnorm;
						to->pdata[chan][i] = (x == 0.0f) ? 0 : roundf_to_WORD(norm * min(1.0f, max(0.0f, GHTp(x, params, &compute_params))));
					}
				} else
					memcpy(to->pdata[chan], from->pdata[chan], layersize * sizeof(WORD));
			}
		}
	} else if (from->type == DATA_FLOAT) {
		if (from->naxes[2] == 3 && params->payne_colourstretchmodel != COL_INDEP) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
			for (size_t i = 0 ; i < layersize ; i++) {
				float L[3] = { 0.0f };
				for (size_t chan = 0; chan < 3 ; chan++)
					L[chan] = (float)from->fpdata[chan][i];
				float x = (int) do_channel[0] * factor_red * L[0] + (int) do_channel[1] * factor_green * L[1] + (int) do_channel[2] * factor_blue * L[2];
				float z = GHTp(x, params, &compute_params);
				for (size_t chan = 0; chan < 3 ; chan++)
					if (do_channel[chan])
						to->fpdata[chan][i] = (x == 0.0f) ? 0.0f : min(1.0f, max(0.0f, L[chan] * (z / x)));
			}
			for (size_t chan = 0 ; chan < 3; chan++) {
				if (!do_channel[chan])
					memcpy(to->fpdata[chan], from->fpdata[chan], layersize * sizeof(float));
			}

		} else {
			for (size_t chan=0; chan < from->naxes[2]; chan++) {
				if (do_channel[chan]) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if (multithreaded)
#endif
					for (size_t i = 0; i < layersize; i++) {
						float x = (float)from->fpdata[chan][i];
						to->fpdata[chan][i] = (x == 0.0f) ? 0.0f : min(1.0f, max(0.0f, GHTp(x, params, &compute_params)));
					}
				} else
					memcpy(to->fpdata[chan], from->fpdata[chan], layersize * sizeof(float));
			}
		}
	}
	else return;
	invalidate_stats_from_fit(to);
}
