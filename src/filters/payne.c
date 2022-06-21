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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"

#include "payne.h"

static int payne_colourstretchmodel = COL_HUMANLUM;
static double payne_b = 0.5, payne_D = 0.0, payne_SP = 0.00, payne_LP = 0.0, payne_HP = 1.00, payne_black_value = 0.0;
static gboolean payne_show_preview;
static int payne_stretchtype = STRETCH_PAYNE_NORMAL;
// Set up intermediate variables for stretch calculations
static double qlp = 0.0, q0 = 0.0, qwp = 0.0, q1 = 0.0, q = 0.0, b1 = 0.0, a1 = 0.0, a2 = 0.0, b2 = 0.0, c2 = 0.0, d2 = 0.0, e2 = 0.0, a3 = 0.0, b3 = 0.0, c3 = 0.0, d3 = 0.0, e3 = 0.0, a4 = 0.0, b4 = 0.0, LPT = 0.0, SPT = 0.0, HPT = 0.0;

static void payne_startup() {
	copy_gfit_to_backup();
}

static void payne_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				_("Generalised Hyperbolic Transformation"));
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int payne_update_preview() {
	copy_backup_to_gfit();
	paynelut(&gfit, payne_b, payne_D, payne_LP, payne_SP, payne_HP, payne_black_value, payne_colourstretchmodel, payne_stretchtype);
	return 0;
}

static int payne_precompute(double B, double D, double LP, double SP, double HP, int stretchtype) {
	if (stretchtype == STRETCH_PAYNE_NORMAL) {
		if (B == -1.0) {
			B = -B;
			qlp = -1.0*log(1.0+D*(SP - LP));
			q0 = qlp - D * LP / (1.0 + D * (SP - LP));
			qwp = log(1.0 + D * (HP - SP));
			q1 = qwp + D * (1.0 - HP) / (1.0 + D * (HP - SP));
			q = 1.0 / (q1 - q0);
			b1 = (1.0 + D * (SP - LP)) / (D * q);
			a2 = (-q0) * q;
			b2 = -q;
			c2 = 1.0 + D * SP;
			d2 = -D;
			a3 = (-q0) * q;
			b3 = q;
			c3 = 1.0 - D * SP;
			d3 = D;
			a4 = (qwp - q0 - D * HP / (1.0 + D * (HP - SP))) * q;
			b4 = q * D / (1.0 + D * (HP -SP));
		} else if ((B != -1.0) && (B < 0.0)) {
			B = -B;
			qlp = (1.0 - pow((1.0 + D * B * (SP - LP)), (B - 1.0) / B)) / (B - 1);
			q0 = qlp - D * LP * (pow((1.0 + D * B * (SP - LP)), -1.0 / B));
			qwp = (pow((1.0 + D * B * (HP - SP)), (B - 1.0) / B) - 1.0) / (B - 1);
			q1 = qwp + D * (1.0 - HP) * (pow((1.0 + D * B * (HP - SP)), -1.0 / B));
			q = 1.0 / (q1 - q0);
			b1 = D * pow(1.0 + D * B * (SP - LP), -1.0 / B) *q;
			a2 = (1/(B-1)-q0) * q;
			b2 = -q/(B-1);
			c2 = 1.0 + D * B * SP;
			d2 = -D * B;
			e2 = (B - 1.0)/B;
			a3 = (-1/(B-1) - q0) *q;
			b3 = q/(B-1);
			c3 = 1.0 - D * B * SP;
			d3 = D * B;
			e3 = (B - 1.0) / B;
			a4 = (qwp - q0 - D * HP * pow((1.0 + D * B * (HP - SP)), -1.0 / B)) * q;
			b4 = D * pow((1.0 + D * B * (HP - SP)), -1.0 / B) * q;
		} else if (B == 0.0) {
			qlp = exp(-D * (SP - LP));
			q0 = qlp - D * LP * exp(-D*(SP - LP));
			qwp = 2.0 - exp(-D * (HP -SP));
			q1 = qwp + D * (1.0 - HP) * exp (-D * (HP - SP));
			q = 1.0 / (q1 - q0);
			a1 = 0.0;
			b1 = D * exp (-D * (SP - LP)) * q;
			a2 = -q0 * q;
			b2 = q;
			c2 = -D * SP;
			d2 = D;
			a3 = (2.0 - q0) * q;
			b3 = -q;
			c3 = D * SP;
			d3 = -D;
			a4 = (qwp - q0 - D * HP * exp(-D * (HP - SP))) * q;
			b4 = D * exp(-D * (HP - SP)) * q;
		} else if (B > 0.0) {
			qlp = pow(( 1 + D * B * (SP - LP)), -1.0/B);
			q0 = qlp - D * LP * pow((1 + D * B * (SP - LP)), -(1.0 + B) / B);
			qwp = 2.0 - pow(1.0 + D * B * (HP - SP), -1.0 / B);
			q1 = qwp + D * (1.0 - HP) * pow((1.0 + D * B * (HP - SP)), -(1.0 + B) / B);
			q = 1.0 / (q1 - q0);
			b1 = D * pow((1 + D * B * (SP - LP)), -(1.0+B)/B) * q;
			a2 = -q0 * q;
			b2 = q;
			c2 = 1.0 + D * B * SP;
			d2 = -D * B;
			e2 = -1.0 / B;
			a3 = (2.0 - q0) * q;
			b3 = -q;
			c3 = 1.0 - D * B * SP;
			d3 = D * B;
			e3 = -1.0 / B;
			a4 = (qwp -q0 -D * HP * pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B)) * q;
			b4 = (D * pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B)) * q;
		}
	} else if (stretchtype == STRETCH_PAYNE_INVERSE) {
		if (B == -1.0) {
			B = -B;
			qlp = -1.0*log(1.0+D*(SP - LP));
			q0 = qlp - D * LP / (1.0 + D * (SP - LP));
			qwp = log(1.0 + D * (HP - SP));
			q1 = qwp + D * (1.0 - HP) / (1.0 + D * (HP - SP));
			q = 1.0 / (q1 - q0);
			LPT = (qlp-q0)*q;
			SPT = q0*q;
			HPT = (qwp-q0)*q;
			b1 = (1.0 + D * (SP - LP)) / (D * q);
			a2 = (1 + D * SP) / D;
			b2 = -1 / D;
			c2 = -q0;
			d2 = -1/q;
			a3 = - (1 - D * SP) / D;
			b3 = 1 / D;
			c3 = q0;
			d3 = 1 / q;
			a4 = HP + (q0-qwp) * (1+D*(HP-SP))/D;
			b4 = (1 + D * (HP - SP) )/(q * D) ;
		} else if ((B != -1.0) && (B < 0.0)) {
			B = -B;
			qlp = (1.0 - pow((1.0 + D * B * (SP - LP)), (B - 1.0) / B)) / (B - 1);
			q0 = qlp - D * LP * (pow((1.0 + D * B * (SP - LP)), -1.0 / B));
			qwp = (pow((1.0 + D * B * (HP - SP)), (B - 1.0) / B) - 1.0) / (B - 1);
			q1 = qwp + D * (1.0 - HP) * (pow((1.0 + D * B * (HP - SP)), -1.0 / B));
			q = 1.0 / (q1 - q0);
			LPT = (qlp-q0)*q;
			SPT = -q0*q;
			HPT = (qwp-q0)*q;
			b1 = pow(1.0 + D * B * (SP - LP), 1.0 / B) / (q * D);
			a2 = (1.0 + D * B * SP) / (D * B);
			b2 = -1.0 / (D * B);
			c2 = -q0*(B-1) + 1;
			d2 = (1 - B) / q;
			e2 = B / (B - 1.0);
			a3 = (D * B * SP - 1) / (D * B);
			b3 = 1 / (D * B);;
			c3 = 1.0 + q0 * (B - 1);
			d3 = (B - 1) / q;
			e3 = B / (B - 1.0);
			a4 = (q0-qwp)/(D * pow((1.0 + D * B * (HP - SP)), -1.0 / B))+HP;
			b4 = 1 / (D * pow((1.0 + D * B * (HP - SP)), -1.0 / B) * q) ;
		} else if (B == 0.0) {
			qlp = exp(-D * (SP - LP));
			q0 = qlp - D * LP * exp(-D*(SP - LP));
			qwp = 2.0 - exp(-D * (HP -SP));
			q1 = qwp + D * (1.0 - HP) * exp (-D * (HP - SP));
			q = 1.0 / (q1 - q0);
			LPT = (qlp-q0)*q;
			SPT = (1-q0)*q;
			HPT = (qwp-q0)*q;
			a1 = 0.0;
			b1 = 1/(D * exp(-D * (SP - LP)) * q);
			a2 = SP;
			b2 = 1 / D;
			c2 = q0;
			d2 = 1/q;
			a3 = SP;
			b3 = -1 / D;
			c3 = (2.0 - q0);
			d3 = -1 / q;
			a4 = (q0 - qwp)/(D * exp(-D * (HP - SP))) + HP;
			b4 = 1/(D * exp(-D * (HP - SP)) * q);
		} else if (B > 0.0) {
			qlp = pow(( 1 + D * B * (SP - LP)), -1.0/B);
			q0 = qlp - D * LP * pow((1 + D * B * (SP - LP)), -(1.0 + B) / B);
			qwp = 2.0 - pow(1.0 + D * B * (HP - SP), -1.0 / B);
			q1 = qwp + D * (1.0 - HP) * pow((1.0 + D * B * (HP - SP)), -(1.0 + B) / B);
			q = 1.0 / (q1 - q0);
			LPT = (qlp-q0)*q;
			SPT = (1-q0)*q;
			HPT = (qwp-q0)*q;
			b1 = 1/(D * pow((1 + D * B * (SP - LP)), -(1.0 + B) / B) * q);
			a2 = 1 / (D * B) + SP;
			b2 = -1/(D * B);
			c2 = q0;
			d2 = 1 / q;
			e2 = -B;
			a3 = -1 / (D * B) + SP;
			b3 = 1 / (D * B);
			c3 = (2.0 - q0);
			d3 = -1 / q;
			e3 = -B;
			a4 = (q0-qwp)/(D * pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B))+HP;
			b4 = 1/((D * pow((1.0 + D * B * (HP - SP)), -(B + 1.0) / B)) * q);
		}
	} else if (stretchtype == STRETCH_ASINH) {
		qlp = -log(D * (SP - LP) + pow((D * D * (SP - LP) * (SP - LP) + 1.0), 0.5));
		q0 = qlp - LP * D * pow((D * D * (SP - LP) * (SP - LP) + 1.0), -0.5);
		qwp = log(D * (HP - SP) + pow((D * D * (HP - SP) * (HP - SP) + 1.0), 0.5));
		q1 = qwp + (1.0 - HP) * D * pow((D * D * (HP - SP) * (HP - SP) +1.0), -0.5);
		q = 1.0 / (q1 - q0);
		a1 = 0.0;
		b1 = D * pow((D * D * (SP - LP) * (SP - LP) + 1),-0.5)*q;
		a2 = -q0 * q;
		b2 = q;
		c2 = -D;
		d2 = D * D;
		e2 = SP;
		a3 = -q0 * q;
		b3 = q;
		c3 = D;
		d3 = D * D;
		e3 = SP;
		a4 = (qwp - HP * D * pow((D * D * (HP - SP) * (HP - SP) + 1.0), -0.5) - q0) * q;
		b4 = D * pow((D * D * (HP - SP) * (HP - SP) + 1.0), -0.5) * q;
	}
	else if (stretchtype == STRETCH_INVASINH) {
		qlp = -log(D * (SP - LP) + pow((D * D * (SP - LP) * (SP - LP) + 1.0), 0.5));
		q0 = qlp - LP * D * pow((D * D * (SP - LP) * (SP - LP) + 1.0), -0.5);
		qwp = log(D * (HP - SP) + pow((D * D * (HP - SP) * (HP - SP) + 1.0), 0.5));
		q1 = qwp + (1.0 - HP) * D * pow((D * D * (HP - SP) * (HP - SP) +1.0), -0.5);
		q = 1.0 / (q1 - q0);
		a1 = 0.0;
		b1 = D * pow((D * D * (SP - LP) * (SP - LP) + 1),-0.5)*q;
		a2 = -q0 * q;
		b2 = q;
		c2 = -D;
		d2 = D * D;
		e2 = SP;
		a3 = -q0 * q;
		b3 = q;
		c3 = D;
		d3 = D * D;
		e3 = SP;
		a4 = (qwp - HP * D * pow((D * D * (HP - SP) * (HP - SP) + 1.0), -0.5) - q0) * q;
		b4 = D * pow((D * D * (HP - SP) * (HP - SP) + 1.0), -0.5) * q;
		LPT = a1 + b1 * LP;
		SPT = a2 + b2 * log(c2 * (SP - e2) + sqrt(d2 * (SP - e2) * (SP - e2) + 1));
		HPT = a4 + b4 * HP;
	}
	else if (stretchtype == STRETCH_LINEAR) {
		return 0;
	}
	return 0;
}

static double payne_compute(double in, double B, double D, double LP, double SP, double HP, double BP, int stretchtype) {
	double out;
	if (stretchtype == STRETCH_PAYNE_NORMAL || stretchtype == STRETCH_PAYNE_INVERSE) {
	in = max (0, (in - BP)/(1 - BP));
		if (D == 0.0) {
		out = in;
		} else if (stretchtype == STRETCH_PAYNE_NORMAL) {
			if (B == -1.0) {
				if (in < LP) {
					out = b1 * in;
				} else if (in < SP) {
					out = a2 + b2 * log(c2 + d2 * in);
				} else if (in < HP) {
					out = a3 + b3 * log(c3 + d3 * in);
				} else {
					out = a4 + b4 * in;
				}
			} else if ((B != -1.0) && (B < 0)) {
				if (in < LP) {
					out = b1 * in;
				} else if (in < SP) {
					out = a2 + b2 * pow((c2 + d2 * in), e2);
				} else if (in < HP) {
					out = a3 + b3 * pow((c3 + d3 * in), e3);
				} else {
					out = a4 + b4 * in;
				}
			} else if (B == 0) {
				if (in < LP) {
					out = a1 + b1 * in;
				} else if (in < SP) {
					out = a2 + b2 * exp(c2 + d2 * in);
				} else if (in < HP) {
					out = a3 + b3 * exp(c3 + d3 * in);
				} else {
					out = a4 + b4 * in;
				}
			} else if (B > 0) {
				if (in < LP) {
					out = b1 * in;
				} else if (in < SP) {
					out = a2 + b2 * pow((c2 + d2 * in), e2);
				} else if (in < HP) {
					out = a3 + b3 * pow((c3 + d3 * in), e3);
				} else {
					out = a4 + b4 * in;
				}
			}
		} else if (stretchtype == STRETCH_PAYNE_INVERSE) {
			if (B == -1.0) {
				if (in < LPT) {
					out = b1 * in;
				} else if (in < SPT) {
					out = a2 + b2 * exp(c2 + d2 * in);
				} else if (in < HPT) {
					out = a3 + b3 * exp(c3 + d3 * in);
				} else {
					out = a4 + b4 * in;
				}
			} else if ((B != -1.0) && (B < 0)) {
				if (in < LPT) {
					out = b1 * in;
				} else if (in < SPT) {
					out = a2 + b2 * pow((c2 + d2 * in), e2);
				} else if (in < HPT) {
					out = a3 + b3 * pow((c3 + d3 * in), e3);
				} else {
					out = a4 + b4 * in;
				}
			} else if (B == 0) {
				if (in < LPT) {
					out = a1 + b1 * in;
				} else if (in < SPT) {
					out = a2 + b2 * log(c2 + d2 * in);
				} else if (in < HPT) {
					out = a3 + b3 * log(c3 + d3 * in);
				} else {
					out = a4 + b4 * in;
				}
			} else if (B > 0) {
				if (in < LPT) {
					out = b1 * in;
				} else if (in < SPT) {
					out = a2 + b2 * pow((c2 + d2 * in), e2);
				} else if (in < HPT) {
					out = a3 + b3 * pow((c3 + d3 * in), e3);
				} else {
					out = a4 + b4 * in;
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
			out = b1 * in;
		} else if (in < SP) {
			out = a2 + b2 * log(c2 + (in - e2) + pow(d2 * (in - e2) * (in - e2) + 1, 0.5));
		} else if (in < HP) {
			out = a3 + b3 * log(c3 * (in - e3) + pow(d3 * (in - e3) * (in - e3) + 1, 0.5));
		} else {
			out = a4 + b4 * in;
		}
	} else if (stretchtype == STRETCH_INVASINH) {
		if (D == 0.0) {
			out = in;
		} else if (in < LPT) {
			out = (in - a1) / b1;
		} else if (in < SPT) {
			double ex = exp((a2 - in) / b2);
			out = e2 - (ex - (1 / ex)) / (2 * c2);
		} else if (in < HPT) {
			double ex = exp((a3 - in) / b3);
			out = e3 - (ex - (1 / ex)) / (2 * c3);
		} else {
			out = (in - a4) / b4;
		}
	}
	return out;
}

static int paynelut_ushort(fits *fit, double B, double D, double LP, double SP, double HP, double BP, gboolean human_luminance, int stretchtype) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	double expD = exp(D) - 1;
	double norm = get_normalized_value(fit);
	double invnorm = 1.0 / norm;
	double factor_red = 0.2126;
	double factor_green = 0.7152;
	double factor_blue = 0.0722;
	if (payne_colourstretchmodel == COL_EVENLUM) {
		factor_red = 0.3333;
		factor_green = 0.3333;
		factor_blue = 0.3333;
	}

	//Set up expressions
	payne_precompute(B, expD, LP, SP, HP, stretchtype);

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {

			double r = buf[RLAYER][i]*invnorm;
			double g = buf[GLAYER][i]*invnorm;
			double b = buf[BLAYER][i]*invnorm;
			if (payne_colourstretchmodel == COL_INDEP) {
				buf[RLAYER][i] = (r == 0.0) ? 0.0 : round_to_WORD(norm * min(1.0, max(0.0, payne_compute(r, B, expD, LP, SP, HP, BP, stretchtype))));
				buf[GLAYER][i] = (g == 0.0) ? 0.0 : round_to_WORD(norm * min(1.0, max(0.0, payne_compute(g, B, expD, LP, SP, HP, BP, stretchtype))));
				buf[BLAYER][i] = (b == 0.0) ? 0.0 : round_to_WORD(norm * min(1.0, max(0.0, payne_compute(b, B, expD, LP, SP, HP, BP, stretchtype))));
			} else {
				double x = factor_red * r + factor_green * g + factor_blue * b;
				double z = payne_compute(x, B, expD, LP, SP, HP, BP, stretchtype);
				buf[RLAYER][i] = (x == 0.0) ? 0.0 : round_to_WORD(norm * min(1.0, max(0.0, ((r - BP) / (1 - BP)) * (z / x))));
				buf[GLAYER][i] = (x == 0.0) ? 0.0 : round_to_WORD(norm * min(1.0, max(0.0, ((g - BP) / (1 - BP)) * (z / x))));
				buf[BLAYER][i] = (x == 0.0) ? 0.0 : round_to_WORD(norm * min(1.0, max(0.0, ((b - BP) / (1 - BP)) * (z / x))));
			}
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {

			double x = buf[RLAYER][i]*invnorm;
			buf[RLAYER][i] = round_to_WORD(norm * min(1.0, max(0.0, payne_compute(x, B, expD, LP, SP, HP, BP, stretchtype))));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int paynelut_float(fits *fit, double B, double D, double LP, double SP, double HP, double BP, gboolean human_luminance, int stretchtype) {
	float *buf[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };
	double expD = exp(D) - 1;

	double factor_red = 0.2126;
	double factor_green = 0.7152;
	double factor_blue = 0.0722;
	if (payne_colourstretchmodel == COL_EVENLUM) {
		factor_red = 0.3333;
		factor_green = 0.3333;
		factor_blue = 0.3333;
	}

	//Set up expressions
	payne_precompute(B, expD, LP, SP, HP, stretchtype);

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {
			double r = buf[RLAYER][i];
			double g = buf[GLAYER][i];
			double b = buf[BLAYER][i];
			if (payne_colourstretchmodel == COL_INDEP) {
				double rprime = (r - BP) / (1.0 - r);
				double gprime = (g - BP) / (1.0 - g);
				double bprime = (b - BP) / (1.0 - b);
				buf[RLAYER][i] = (r == 0.0) ? 0.0 : min(1.0, max(0.0, (rprime / r) * payne_compute(r, B, expD, LP, SP, HP, BP, stretchtype)));
				buf[GLAYER][i] = (g == 0.0) ? 0.0 : min(1.0, max(0.0, (gprime / g) * payne_compute(g, B, expD, LP, SP, HP, BP, stretchtype)));
				buf[BLAYER][i] = (b == 0.0) ? 0.0 : min(1.0, max(0.0, (bprime / b) * payne_compute(b, B, expD, LP, SP, HP, BP, stretchtype)));
			} else {
				double x = factor_red * r + factor_green * g + factor_blue * b;
				double z = payne_compute(x, B, expD, LP, SP, HP, BP, stretchtype);
				buf[RLAYER][i] = (x == 0.0) ? 0.0 : min(1.0, max(0.0, ((r - BP) / (1 - BP)) * (z / x)));
				buf[GLAYER][i] = (x == 0.0) ? 0.0 : min(1.0, max(0.0, ((g - BP) / (1 - BP)) * (z / x)));
				buf[BLAYER][i] = (x == 0.0) ? 0.0 : min(1.0, max(0.0, ((b - BP) / (1 - BP)) * (z / x)));
			}
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {
			double x = buf[RLAYER][i];
			buf[RLAYER][i] = payne_compute(x, B, expD, LP, SP, HP, BP, stretchtype);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int paynelut(fits *fit, double beta, double intensity, double lower, double shoulder, double headroom, double offset, gboolean human_luminance, int stretchtype) {
	if (fit->type == DATA_USHORT)
		return paynelut_ushort(fit, beta, intensity, lower, shoulder, headroom, offset, human_luminance, stretchtype);
	if (fit->type == DATA_FLOAT)
		return paynelut_float(fit, beta, intensity, lower, shoulder, headroom, offset, human_luminance, stretchtype);
	return 1;
}

static void apply_payne_changes() {
	gboolean status = (payne_b != 0.0) || (payne_D != 0.0) || (payne_SP != 0.0) || (payne_HP != 1.0) || (payne_LP != 0.0) || (payne_black_value != 0.0);
	payne_close(!status);
}

void apply_payne_cancel() {
	payne_close(TRUE);
	siril_close_dialog("payne_dialog");
}

/*** callbacks **/

void on_payne_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_D = GTK_SPIN_BUTTON(lookup_widget("spin_payneD"));
	GtkSpinButton *spin_b = GTK_SPIN_BUTTON(lookup_widget("spin_payneb"));
	GtkSpinButton *spin_SP = GTK_SPIN_BUTTON(lookup_widget("spin_payneSP"));
	GtkSpinButton *spin_HP = GTK_SPIN_BUTTON(lookup_widget("spin_payneHP"));
	GtkSpinButton *spin_black_p = GTK_SPIN_BUTTON(lookup_widget("black_point_spin_payne"));
	GtkComboBox *combo_payne_type = GTK_COMBO_BOX(lookup_widget("combo_payneType"));
	GtkComboBox *combo_payne_colourstretchmodel = GTK_COMBO_BOX(lookup_widget("combo_payne_colour_stretch_model"));

	payne_startup();
	payne_D = 0.0;
	payne_b = 0.0;
	payne_SP = 0.0;
	payne_HP = 1.0;
	payne_LP = 0.0;
	payne_black_value = 0.0;
	payne_colourstretchmodel = COL_HUMANLUM;

	set_notify_block(TRUE);
	if (com.uniq->fit->naxes[2] > 1) {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("combo_payne_colour_stretch_model")), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("payne_colourstretchlabel")), TRUE);
	} else {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("combo_payne_colour_stretch_model")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("payne_colourstretchlabel")), FALSE);
	}
	gtk_combo_box_set_active(combo_payne_colourstretchmodel, payne_colourstretchmodel);
	gtk_combo_box_set_active(combo_payne_type, STRETCH_PAYNE_NORMAL);
	gtk_spin_button_set_value(spin_D, payne_D);
	gtk_spin_button_set_value(spin_b, payne_b);
	gtk_spin_button_set_value(spin_SP, payne_SP);
	gtk_spin_button_set_value(spin_HP, payne_HP);
	gtk_spin_button_set_value(spin_black_p, payne_black_value);
	gtk_spin_button_set_increments(spin_black_p, 0.001, 0.01);
	set_notify_block(FALSE);

	payne_show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("payne_preview")));

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}

void on_payne_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_payne_cancel();
}

void on_payne_apply_clicked(GtkButton *button, gpointer user_data) {
	if (payne_show_preview == FALSE) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = payne_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}

	apply_payne_changes();
	siril_close_dialog("payne_dialog");
}

void on_payne_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_payne_changes();
}

void on_payne_undo_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_b = GTK_SPIN_BUTTON(lookup_widget("spin_payneb"));
	GtkSpinButton *spin_d = GTK_SPIN_BUTTON(lookup_widget("spin_payneD"));
	GtkSpinButton *spin_SP = GTK_SPIN_BUTTON(lookup_widget("spin_payneSP"));
	GtkSpinButton *spin_HP = GTK_SPIN_BUTTON(lookup_widget("spin_payneHP"));
	GtkSpinButton *spin_black_p = GTK_SPIN_BUTTON(lookup_widget("black_point_spin_payne"));
	GtkComboBox *combo_payne_colourstretchmodel = GTK_COMBO_BOX(lookup_widget("combo_payne_colour_stretch_model"));

	payne_D = 0.0;
	payne_b = 0.0;
	payne_black_value = 0.0;
	payne_SP = 0.0;
	payne_HP = 1.0;
	payne_colourstretchmodel = COL_HUMANLUM;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_d, payne_D);
	gtk_spin_button_set_value(spin_b, payne_b);
	gtk_spin_button_set_value(spin_SP, payne_SP);
	gtk_spin_button_set_value(spin_HP, payne_HP);
	gtk_spin_button_set_value(spin_black_p, payne_black_value);
	gtk_combo_box_set_active(combo_payne_colourstretchmodel,payne_colourstretchmodel);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}

/*** adjusters **/
void on_spin_payneD_value_changed(GtkSpinButton *button, gpointer user_data) {
	payne_D = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}

void on_spin_payneb_value_changed(GtkSpinButton *button, gpointer user_data) {
	payne_b = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}

void on_spin_payneSP_value_changed(GtkSpinButton *button, gpointer user_data) {
	GtkSpinButton *spin_HP = GTK_SPIN_BUTTON(lookup_widget("spin_payneHP"));
	GtkSpinButton *spin_LP = GTK_SPIN_BUTTON(lookup_widget("spin_payneLP"));
	payne_SP = gtk_spin_button_get_value(button);
// If SP goes higher than HP or lower than LP, drag HP or LP with it.
	if (payne_SP < payne_LP) {
		gtk_spin_button_set_value(spin_LP, payne_SP);
	}
	if (payne_SP > payne_HP) {
		gtk_spin_button_set_value(spin_HP, payne_SP);
	}
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}
void on_spin_payneLP_value_changed(GtkSpinButton *button, gpointer user_data) {
	payne_LP = gtk_spin_button_get_value(button);
	if (payne_LP > payne_SP) {
		gtk_spin_button_set_value(button,payne_SP);
		siril_log_message(_("Shadow preservation point cannot be set higher than stretch focal point\n"));
	}
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}

void on_spin_payneHP_value_changed(GtkSpinButton *button, gpointer user_data) {
	payne_HP = gtk_spin_button_get_value(button);
	if (payne_HP < payne_SP) {
		gtk_spin_button_set_value(button, payne_SP);
		siril_log_message(_("Highlight preservation point cannot be set lower than stretch focal point\n"));
	}
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}

void on_blackpoint_spin_payne_value_changed(GtkSpinButton *button, gpointer user_data) {
	payne_black_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}

void on_combo_payne_colour_stretch_model_changed(GtkComboBox *combo, gpointer user_data) {
	payne_colourstretchmodel = gtk_combo_box_get_active(combo);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}

void on_payne_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (payne_show_preview == TRUE) {
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = payne_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
	payne_show_preview = !payne_show_preview;
}


void on_combo_payneType_changed(GtkComboBox *combo, gpointer user_data) {
	payne_stretchtype = gtk_combo_box_get_active(combo);
	switch (payne_stretchtype) {
		case STRETCH_LINEAR:
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneD")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneD")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneD")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneB")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneb")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneb")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneLP")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneLP")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneLP")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneSP")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneSP")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneSP")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneHP")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneHP")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneHP")), FALSE);
			break;
		case STRETCH_PAYNE_NORMAL:
		case STRETCH_PAYNE_INVERSE:
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneD")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneD")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneD")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneB")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneb")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneb")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneLP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneLP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneLP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneSP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneSP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneSP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneHP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneHP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneHP")), TRUE);
			break;
		case STRETCH_ASINH:
		case STRETCH_INVASINH:
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneD")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneD")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneD")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneB")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneb")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneb")), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneLP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneLP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneLP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneSP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneSP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneSP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("label_payneHP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("spin_payneHP")), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(lookup_widget("scale_payneHP")), TRUE);
			break;
	}
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}
