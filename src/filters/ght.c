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

#include <glib.h>
#include "ght.h"
#include "core/proto.h"
#include "algos/statistics.h"

double GHT(double in, double B, double D, double LP, double SP, double HP, double BP, int stretchtype, ght_compute_params c) {
	double out;
	if (stretchtype != STRETCH_LINEAR) {
		BP = 0.0;
	}
	if (stretchtype == STRETCH_PAYNE_NORMAL || stretchtype == STRETCH_PAYNE_INVERSE) {
		in = max (0, (in - BP)/(1 - BP));
		if (D == 0.0) {
			out = in;
		} else if (stretchtype == STRETCH_PAYNE_NORMAL) {
			if (B == -1.0) {
				if (in < LP) {
					out = c.b1 * in;
				} else if (in < SP) {
					out = c.a2 + c.b2 * log(c.c2 + c.d2 * in);
				} else if (in < HP) {
					out = c.a3 + c.b3 * log(c.c3 + c.d3 * in);
				} else {
					out = c.a4 + c.b4 * in;
				}
			} else if (B < 0) {
				if (in < LP) {
					out = c.b1 * in;
				} else if (in < SP) {
					out = c.a2 + c.b2 * pow((c.c2 + c.d2 * in), c.e2);
				} else if (in < HP) {
					out = c.a3 + c.b3 * pow((c.c3 + c.d3 * in), c.e3);
				} else {
					out = c.a4 + c.b4 * in;
				}
			} else if (B == 0) {
				if (in < LP) {
					out = c.a1 + c.b1 * in;
				} else if (in < SP) {
					out = c.a2 + c.b2 * exp(c.c2 + c.d2 * in);
				} else if (in < HP) {
					out = c.a3 + c.b3 * exp(c.c3 + c.d3 * in);
				} else {
					out = c.a4 + c.b4 * in;
				}
			} else /*if (B > 0)*/ {
				if (in < LP) {
					out = c.b1 * in;
				} else if (in < SP) {
					out = c.a2 + c.b2 * pow((c.c2 + c.d2 * in), c.e2);
				} else if (in < HP) {
					out = c.a3 + c.b3 * pow((c.c3 + c.d3 * in), c.e3);
				} else {
					out = c.a4 + c.b4 * in;
				}
			}
		} else /*if (stretchtype == STRETCH_PAYNE_INVERSE)*/ {
			if (B == -1.0) {
				if (in < c.LPT) {
					out = c.b1 * in;
				} else if (in < c.SPT) {
					out = c.a2 + c.b2 * exp(c.c2 + c.d2 * in);
				} else if (in < c.HPT) {
					out = c.a3 + c.b3 * exp(c.c3 + c.d3 * in);
				} else {
					out = c.a4 + c.b4 * in;
				}
			} else if (B < 0) {
				if (in < c.LPT) {
					out = c.b1 * in;
				} else if (in < c.SPT) {
					out = c.a2 + c.b2 * pow((c.c2 + c.d2 * in), c.e2);
				} else if (in < c.HPT) {
					out = c.a3 + c.b3 * pow((c.c3 + c.d3 * in), c.e3);
				} else {
					out = c.a4 + c.b4 * in;
				}
			} else if (B == 0) {
				if (in < c.LPT) {
					out = c.a1 + c.b1 * in;
				} else if (in < c.SPT) {
					out = c.a2 + c.b2 * log(c.c2 + c.d2 * in);
				} else if (in < c.HPT) {
					out = c.a3 + c.b3 * log(c.c3 + c.d3 * in);
				} else {
					out = c.a4 + c.b4 * in;
				}
			} else /* if (B > 0) */{
				if (in < c.LPT) {
					out = c.b1 * in;
				} else if (in < c.SPT) {
					out = c.a2 + c.b2 * pow((c.c2 + c.d2 * in), c.e2);
				} else if (in < c.HPT) {
					out = c.a3 + c.b3 * pow((c.c3 + c.d3 * in), c.e3);
				} else {
					out = c.a4 + c.b4 * in;
				}
			}
		}
	} else if (stretchtype == STRETCH_LINEAR) {
		out = max(0, (in - BP) / (1 - BP));
	} else if (stretchtype == STRETCH_ASINH) {
		in = max(0, (in - BP) / (1 - BP));
		if (D == 0.0) {
			out = in;
		} else if (in < LP) {
			out = c.b1 * in;
		} else if (in < SP) {
			out = c.a2 + c.b2 * log(c.c2 + (in - c.e2) + pow(c.d2 * (in - c.e2) * (in - c.e2) + 1, 0.5));
		} else if (in < HP) {
			out = c.a3 + c.b3 * log(c.c3 * (in - c.e3) + pow(c.d3 * (in - c.e3) * (in - c.e3) + 1, 0.5));
		} else {
			out = c.a4 + c.b4 * in;
		}
	} else /* if (stretchtype == STRETCH_INVASINH)*/ {
		if (D == 0.0) {
			out = in;
		} else if (in < c.LPT) {
			out = (in - c.a1) / c.b1;
		} else if (in < c.SPT) {
			double ex = exp((c.a2 - in) / c.b2);
			out = c.e2 - (ex - (1 / ex)) / (2 * c.c2);
		} else if (in < c.HPT) {
			double ex = exp((c.a3 - in) / c.b3);
			out = c.e3 - (ex - (1 / ex)) / (2 * c.c3);
		} else {
			out = (in - c.a4) / c.b4;
		}
	}
	return out;
}

double GHTp(double in, ght_params params, ght_compute_params compute_params) {
	return GHT(in, params.B, params.D, params.LP, params.SP, params.HP, params.BP, params.stretchtype, compute_params);
}

int GHTsetup(ght_compute_params *c, double B, double D, double LP, double SP, double HP, int stretchtype) {
	if (D == 0.0 || stretchtype == STRETCH_LINEAR) {
		c->qlp = c->q0 = c->qwp = c->q1 = c->q = c->b1 = c->a2 = c->b2 = c->c2 = c->d2 = c->a3 = c->b3 = c->c3 = c->d3 = c->a4 = c->b4 = 0.0;
		return 0;
	}
	else if (stretchtype == STRETCH_PAYNE_NORMAL) {
		if (B == -1.0) {
			//B = -B;
			c->qlp = -1.0*log1p(D*(SP - LP));
			c->q0 = c->qlp - D * LP / (1.0 + D * (SP - LP));
			c->qwp = log1p(D * (HP - SP));
			c->q1 = c->qwp + D * (1.0 - HP) / (1.0 + D * (HP - SP));
			c->q = 1.0 / (c->q1 - c->q0);
			c->b1 = (1.0 + D * (SP - LP)) / (D * c->q);
			c->a2 = (-c->q0) * c->q;
			c->b2 = -c->q;
			c->c2 = 1.0 + D * SP;
			c->d2 = -D;
			c->a3 = (-c->q0) * c->q;
			c->b3 = c->q;
			c->c3 = 1.0 - D * SP;
			c->d3 = D;
			c->a4 = (c->qwp - c->q0 - D * HP / (1.0 + D * (HP - SP))) * c->q;
			c->b4 = c->q * D / (1.0 + D * (HP - SP));
		} else if (B < 0.0) {
			B = -B;
			c->qlp = (1.0 - pow((1.0 + D * B * (SP - LP)), (B - 1.0) / B)) / (B - 1);
			c->q0 = c->qlp - D * LP * (pow((1.0 + D * B * (SP - LP)), -1.0 / B));
			c->qwp = (pow((1.0 + D * B * (HP - SP)), (B - 1.0) / B) - 1.0) / (B - 1);
			c->q1 = c->qwp + D * (1.0 - HP) * (pow((1.0 + D * B * (HP - SP)), -1.0 / B));
			c->q = 1.0 / (c->q1 - c->q0);
			c->b1 = D * pow(1.0 + D * B * (SP - LP), -1.0 / B) *c->q;
			c->a2 = (1/(B-1)-c->q0) * c->q;
			c->b2 = -c->q/(B-1);
			c->c2 = 1.0 + D * B * SP;
			c->d2 = -D * B;
			c->e2 = (B - 1.0)/B;
			c->a3 = (-1/(B-1) - c->q0) *c->q;
			c->b3 = c->q/(B-1);
			c->c3 = 1.0 - D * B * SP;
			c->d3 = D * B;
			c->e3 = (B - 1.0) / B;
			c->a4 = (c->qwp - c->q0 - D * HP * pow((1.0 + D * B * (HP - SP)), -1.0 / B)) * c->q;
			c->b4 = D * pow((1.0 + D * B * (HP - SP)), -1.0 / B) * c->q;
		} else if (B == 0.0) {
			c->qlp = exp(-D * (SP - LP));
			c->q0 = c->qlp - D * LP * exp(-D*(SP - LP));
			c->qwp = 2.0 - exp(-D * (HP -SP));
			c->q1 = c->qwp + D * (1.0 - HP) * exp (-D * (HP - SP));
			c->q = 1.0 / (c->q1 - c->q0);
			c->a1 = 0.0;
			c->b1 = D * exp (-D * (SP - LP)) * c->q;
			c->a2 = -c->q0 * c->q;
			c->b2 = c->q;
			c->c2 = -D * SP;
			c->d2 = D;
			c->a3 = (2.0 - c->q0) * c->q;
			c->b3 = -c->q;
			c->c3 = D * SP;
			c->d3 = -D;
			c->a4 = (c->qwp - c->q0 - D * HP * exp(-D * (HP - SP))) * c->q;
			c->b4 = D * exp(-D * (HP - SP)) * c->q;
		} else if (B > 0.0) {
			c->qlp = pow(( 1 + D * B * (SP - LP)), -1.0/B);
			c->q0 = c->qlp - D * LP * pow((1 + D * B * (SP - LP)), -(1.0 + B) / B);
			c->qwp = 2.0 - pow(1.0 + D * B * (HP - SP), -1.0 / B);
			c->q1 = c->qwp + D * (1.0 - HP) * pow((1.0 + D * B * (HP - SP)), -(1.0 + B) / B);
			c->q = 1.0 / (c->q1 - c->q0);
			c->b1 = D * pow((1 + D * B * (SP - LP)), -(1.0+B)/B) * c->q;
			c->a2 = -c->q0 * c->q;
			c->b2 = c->q;
			c->c2 = 1.0 + D * B * SP;
			c->d2 = -D * B;
			c->e2 = -1.0 / B;
			c->a3 = (2.0 - c->q0) * c->q;
			c->b3 = -c->q;
			c->c3 = 1.0 - D * B * SP;
			c->d3 = D * B;
			c->e3 = -1.0 / B;
			c->a4 = (c->qwp - c->q0 - D * HP * pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B)) * c->q;
			c->b4 = (D * pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B)) * c->q;
		}
	} else if (stretchtype == STRETCH_PAYNE_INVERSE) {
		if (B == -1.0) {
			//B = -B;
			c->qlp = -1.0 * log1p(D * (SP - LP));
			c->q0 = c->qlp - D * LP / (1.0 + D * (SP - LP));
			c->qwp = log1p(D * (HP - SP));
			c->q1 = c->qwp + D * (1.0 - HP) / (1.0 + D * (HP - SP));
			c->q = 1.0 / (c->q1 - c->q0);
			c->LPT = (c->qlp-c->q0)*c->q;
			c->SPT = c->q0*c->q;
			c->HPT = (c->qwp-c->q0)*c->q;
			c->b1 = (1.0 + D * (SP - LP)) / (D * c->q);
			c->a2 = (1 + D * SP) / D;
			c->b2 = -1 / D;
			c->c2 = - c->q0;
			c->d2 = - 1/c->q;
			c->a3 = - (1 - D * SP) / D;
			c->b3 = 1 / D;
			c->c3 = c->q0;
			c->d3 = 1 / c->q;
			c->a4 = HP + (c->q0-c->qwp) * (1+D*(HP-SP))/D;
			c->b4 = (1 + D * (HP - SP) )/(c->q * D) ;
		} else if (B < 0.0) {
			B = -B;
			c->qlp = (1.0 - pow((1.0 + D * B * (SP - LP)), (B - 1.0) / B)) / (B - 1);
			c->q0 = c->qlp - D * LP * (pow((1.0 + D * B * (SP - LP)), -1.0 / B));
			c->qwp = (pow((1.0 + D * B * (HP - SP)), (B - 1.0) / B) - 1.0) / (B - 1);
			c->q1 = c->qwp + D * (1.0 - HP) * (pow((1.0 + D * B * (HP - SP)), -1.0 / B));
			c->q = 1.0 / (c->q1 - c->q0);
			c->LPT = (c->qlp-c->q0)*c->q;
			c->SPT = -c->q0*c->q;
			c->HPT = (c->qwp-c->q0)*c->q;
			c->b1 = pow(1.0 + D * B * (SP - LP), 1.0 / B) / (c->q * D);
			c->a2 = (1.0 + D * B * SP) / (D * B);
			c->b2 = -1.0 / (D * B);
			c->c2 = -c->q0*(B-1) + 1;
			c->d2 = (1 - B) / c->q;
			c->e2 = B / (B - 1.0);
			c->a3 = (D * B * SP - 1) / (D * B);
			c->b3 = 1 / (D * B);
			c->c3 = 1.0 + c->q0 * (B - 1);
			c->d3 = (B - 1) / c->q;
			c->e3 = B / (B - 1.0);
			c->a4 = (c->q0-c->qwp)/(D * pow((1.0 + D * B * (HP - SP)), -1.0 / B))+HP;
			c->b4 = 1 / (D * pow((1.0 + D * B * (HP - SP)), -1.0 / B) * c->q) ;
		} else if (B == 0.0) {
			c->qlp = exp(-D * (SP - LP));
			c->q0 = c->qlp - D * LP * exp(-D*(SP - LP));
			c->qwp = 2.0 - exp(-D * (HP -SP));
			c->q1 = c->qwp + D * (1.0 - HP) * exp (-D * (HP - SP));
			c->q = 1.0 / (c->q1 - c->q0);
			c->LPT = (c->qlp-c->q0)*c->q;
			c->SPT = (1-c->q0)*c->q;
			c->HPT = (c->qwp-c->q0)*c->q;
			c->a1 = 0.0;
			c->b1 = 1/(D * exp(-D * (SP - LP)) * c->q);
			c->a2 = SP;
			c->b2 = 1 / D;
			c->c2 = c->q0;
			c->d2 = 1/c->q;
			c->a3 = SP;
			c->b3 = -1 / D;
			c->c3 = (2.0 - c->q0);
			c->d3 = -1 / c->q;
			c->a4 = (c->q0 - c->qwp)/(D * exp(-D * (HP - SP))) + HP;
			c->b4 = 1/(D * exp(-D * (HP - SP)) * c->q);
		} else if (B > 0.0) {
			c->qlp = pow(( 1 + D * B * (SP - LP)), -1.0/B);
			c->q0 = c->qlp - D * LP * pow((1 + D * B * (SP - LP)), -(1.0 + B) / B);
			c->qwp = 2.0 - pow(1.0 + D * B * (HP - SP), -1.0 / B);
			c->q1 = c->qwp + D * (1.0 - HP) * pow((1.0 + D * B * (HP - SP)), -(1.0 + B) / B);
			c->q = 1.0 / (c->q1 - c->q0);
			c->LPT = (c->qlp-c->q0)*c->q;
			c->SPT = (1-c->q0)*c->q;
			c->HPT = (c->qwp-c->q0)*c->q;
			c->b1 = 1/(D * pow((1 + D * B * (SP - LP)), -(1.0 + B) / B) * c->q);
			c->a2 = 1 / (D * B) + SP;
			c->b2 = -1/(D * B);
			c->c2 = c->q0;
			c->d2 = 1 / c->q;
			c->e2 = -B;
			c->a3 = -1 / (D * B) + SP;
			c->b3 = 1 / (D * B);
			c->c3 = (2.0 - c->q0);
			c->d3 = -1 / c->q;
			c->e3 = -B;
			c->a4 = (c->q0-c->qwp)/(D * pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B))+HP;
			c->b4 = 1/((D * pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B)) * c->q);
		}
	} else if (stretchtype == STRETCH_ASINH) {
		c->qlp = -log(D * (SP - LP) + pow((D * D * (SP - LP) * (SP - LP) + 1.0), 0.5));
		c->q0 = c->qlp - LP * D * pow((D * D * (SP - LP) * (SP - LP) + 1.0), -0.5);
		c->qwp = log(D * (HP - SP) + pow((D * D * (HP - SP) * (HP - SP) + 1.0), 0.5));
		c->q1 = c->qwp + (1.0 - HP) * D * pow((D * D * (HP - SP) * (HP - SP) +1.0), -0.5);
		c->q = 1.0 / (c->q1 - c->q0);
		c->a1 = 0.0;
		c->b1 = D * pow((D * D * (SP - LP) * (SP - LP) + 1),-0.5)*c->q;
		c->a2 = -c->q0 * c->q;
		c->b2 = c->q;
		c->c2 = -D;
		c->d2 = D * D;
		c->e2 = SP;
		c->a3 = -c->q0 * c->q;
		c->b3 = c->q;
		c->c3 = D;
		c->d3 = D * D;
		c->e3 = SP;
		c->a4 = (c->qwp - HP * D * pow((D * D * (HP - SP) * (HP - SP) + 1.0), -0.5) - c->q0) * c->q;
		c->b4 = D * pow((D * D * (HP - SP) * (HP - SP) + 1.0), -0.5) * c->q;
	}
	else if (stretchtype == STRETCH_INVASINH) {
		c->qlp = -log(D * (SP - LP) + pow((D * D * (SP - LP) * (SP - LP) + 1.0), 0.5));
		c->q0 = c->qlp - LP * D * pow((D * D * (SP - LP) * (SP - LP) + 1.0), -0.5);
		c->qwp = log(D * (HP - SP) + pow((D * D * (HP - SP) * (HP - SP) + 1.0), 0.5));
		c->q1 = c->qwp + (1.0 - HP) * D * pow((D * D * (HP - SP) * (HP - SP) +1.0), -0.5);
		c->q = 1.0 / (c->q1 - c->q0);
		c->a1 = 0.0;
		c->b1 = D * pow((D * D * (SP - LP) * (SP - LP) + 1),-0.5)*c->q;
		c->a2 = -c->q0 * c->q;
		c->b2 = c->q;
		c->c2 = -D;
		c->d2 = D * D;
		c->e2 = SP;
		c->a3 = -c->q0 * c->q;
		c->b3 = c->q;
		c->c3 = D;
		c->d3 = D * D;
		c->e3 = SP;
		c->a4 = (c->qwp - HP * D * pow((D * D * (HP - SP) * (HP - SP) + 1.0), -0.5) - c->q0) * c->q;
		c->b4 = D * pow((D * D * (HP - SP) * (HP - SP) + 1.0), -0.5) * c->q;
		c->LPT = c->a1 + c->b1 * LP;
		c->SPT = c->a2 + c->b2 * log(c->c2 * (SP - c->e2) + sqrt(c->d2 * (SP - c->e2) * (SP - c->e2) + 1));
		c->HPT = c->a4 + c->b4 * HP;
	}
	return 0;
}

void apply_linked_ght_to_fits(fits *from, fits *to, ght_params params, struct ght_compute_params compute_params) {
	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t ndata = from->naxes[0] * from->naxes[1] * from->naxes[2];
	const size_t layersize = from->naxes[0] * from->naxes[1];
	g_assert(from->type == to->type);
	double factor_red = 0.2126;
	double factor_green = 0.7152;
	double factor_blue = 0.0722;
	if (params.payne_colourstretchmodel == COL_EVENLUM) {
		factor_red = 1.0/3.0;
		factor_green = 1.0/3.0;
		factor_blue = 1.0/3.0;
	}

//	siril_log_message(_("Applying GHT with values %lf, %lf, %lf, %lf, %lf\n"),
//			params.B, params.D, params.LP, params.SP, params.HP);
	// Do calcs that can be done prior to the loop
	GHTsetup(&compute_params, params.B, params.D, params.LP, params.SP, params.HP, params.stretchtype);
	if (from->type == DATA_USHORT) {
		double norm = get_normalized_value(from);
		double invnorm = 1.0f / norm;
		if (from->naxes[2] == 3 && params.payne_colourstretchmodel != COL_INDEP) {
//			siril_log_message(_("Luminance stretch (16bit)\n"));
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0 ; i < layersize ; i++) {
				double r = (double)from->pdata[RLAYER][i] * invnorm;
				double g = (double)from->pdata[GLAYER][i] * invnorm;
				double b = (double)from->pdata[BLAYER][i] * invnorm;
				double x = factor_red * r + factor_green * g + factor_blue * b;
				double z = GHTp(x, params, compute_params);
				to->pdata[RLAYER][i] = (x == 0.0) ? 0 : round_to_WORD(norm * min(1.0, max(0.0, r * (z / x))));
				to->pdata[GLAYER][i] = (x == 0.0) ? 0 : round_to_WORD(norm * min(1.0, max(0.0, g * (z / x))));
				to->pdata[BLAYER][i] = (x == 0.0) ? 0 : round_to_WORD(norm * min(1.0, max(0.0, b * (z / x))));
			}
		} else {
//			siril_log_message(_("Independent stretch (16bit)\n"));
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0; i < ndata; i++) {
				double pxl = (double)from->data[i] * invnorm;
				to->data[i] = (pxl == 0.0) ? 0 : round_to_WORD(norm * min(1.0, max(0.0, GHTp(pxl, params, compute_params))));
			}
		}
	} else if (from->type == DATA_FLOAT) {
		if (from->naxes[2] == 3 && params.payne_colourstretchmodel != COL_INDEP) {
//			siril_log_message(_("Luminance stretch (32bit)\n"));
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0 ; i < layersize ; i++) {
				double r = (double)from->fpdata[RLAYER][i];
				double g = (double)from->fpdata[GLAYER][i];
				double b = (double)from->fpdata[BLAYER][i];
				double x = factor_red * r + factor_green * g + factor_blue * b;
				double z = GHTp(x, params, compute_params);
				to->fpdata[RLAYER][i] = (x == 0.0) ? 0.0 : (float)min(1.0, max(0.0, r * (z / x)));
				to->fpdata[GLAYER][i] = (x == 0.0) ? 0.0 : (float)min(1.0, max(0.0, g * (z / x)));
				to->fpdata[BLAYER][i] = (x == 0.0) ? 0.0 : (float)min(1.0, max(0.0, b * (z / x)));
			}
		} else {
//			siril_log_message(_("Independent stretch (32bit)\n"));
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
			for (size_t i = 0; i < ndata; i++) {
				to->fdata[i] = (from->fdata[i] == 0.0) ? 0.0 : min(1.0, max(0.0, GHTp(from->fdata[i], params, compute_params)));
			}
		}
	}
	else return;
	invalidate_stats_from_fit(to);
}
