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

static gboolean payne_rgb_space = FALSE;
static double payne_b = 0.5, payne_D = 0.0, payne_SP = 0.00, payne_LP = 0.0, payne_HP = 1.00, payne_black_value = 0.0;
static gboolean payne_show_preview;
static gboolean payne_inverse = FALSE;
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
	paynelut(&gfit, payne_b, payne_D, payne_LP, payne_SP, payne_HP, payne_black_value, payne_rgb_space, payne_inverse);
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
	}
	return 0;
}

static double payne_compute(double in, double x, double B, double D, double LP, double SP, double HP, int stretchtype) {
	double out;
	if (stretchtype == STRETCH_PAYNE_NORMAL || stretchtype == STRETCH_PAYNE_INVERSE) {
		if (D == 0.0) {
		out = in;
		} else if (stretchtype == STRETCH_PAYNE_NORMAL) {
			if (B == -1.0) {
				if (x < LP) {
					out = b1 * in;
				} else if (x < SP) {
					out = a2 + b2 * log(c2 + d2 * in);
				} else if (x < HP) {
					out = a3 + b3 * log(c3 + d3 * in);
				} else {
					out = a4 + b4 * in;
				}
			} else if ((B != -1.0) && (B < 0)) {
				if (x < LP) {
					out = b1 * in;
				} else if (x < SP) {
					out = a2 + b2 * pow((c2 + d2 * in), e2);
				} else if (x < HP) {
					out = a3 + b3 * pow((c3 + d3 * in), e3);
				} else {
					out = a4 + b4 * in;
				}
			} else if (B == 0) {
				if (x < LP) {
					out = a1 + b1 * in;
				} else if (x < SP) {
					out = a2 + b2 * exp(c2 + d2 * in);
				} else if (x < HP) {
					out = a3 + b3 * exp(c3 + d3 * in);
				} else {
					out = a4 + b4 * in;
				}
			} else if (B > 0) {
				if (x < LP) {
					out = b1 * in;
				} else if (x < SP) {
					out = a2 + b2 * pow((c2 + d2 * in), e2);
				} else if (x < HP) {
					out = a3 + b3 * pow((c3 + d3 * in), e3);
				} else {
					out = a4 + b4 * in;
				}
			}
		} else if (stretchtype == STRETCH_PAYNE_INVERSE) {
			if (B == -1.0) {
				if (x < LPT) {
					out = b1 * x;
				} else if (x < SPT) {
					out = a2 + b2 * exp(c2 + d2 * x);
				} else if (x < HPT) {
					out = a3 + b3 * exp(c3 + d3 * x);
				} else {
					out = a4 + b4 * x;
				}
			} else if ((B != -1.0) && (B < 0)) {
				if (x < LPT) {
					out = b1 * x;
				} else if (x < SPT) {
					out = a2 + b2 * pow((c2 + d2 * x), e2);
				} else if (x < HPT) {
					out = a3 + b3 * pow((c3 + d3 * x), e3);
				} else {
					out = a4 + b4 * x;
				}
			} else if (B == 0) {
				if (x < LPT) {
					out = a1 + b1 * x;
				} else if (x < SPT) {
					out = a2 + b2 * log(c2 + d2 * x);
				} else if (x < HPT) {
					out = a3 + b3 * log(c3 + d3 * x);
				} else {
					out = a4 + b4 * x;
				}
			} else if (B > 0) {
				if (x < LPT) {
					out = b1 * x;
				} else if (x < SPT) {
					out = a2 + b2 * pow((c2 + d2 * x), e2);
				} else if (x < HPT) {
					out = a3 + b3 * pow((c3 + d3 * x), e3);
				} else {
					out = a4 + b4 * x;
				}
			}
		}
	}
	return out;
}

static int paynelut_ushort(fits *fit, double B, double D, double LP, double SP, double HP, double BP, gboolean human_luminance, int stretchtype) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	double expD = exp(D) - 1;
	double norm = get_normalized_value(fit);
	double invnorm = 1.0 / norm;
	double factor_red = human_luminance ? 0.2126 : 0.3333;
	double factor_green = human_luminance ? 0.7152 : 0.3333;
	double factor_blue = human_luminance ? 0.0722 : 0.3333;

	//Set up expressions
	payne_precompute(B, expD, LP, SP, HP, stretchtype);

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {

			double r = buf[RLAYER][i]*invnorm - BP;
			double g = buf[GLAYER][i]*invnorm - BP;
			double b = buf[BLAYER][i]*invnorm - BP;
			double x = factor_red * r + factor_green * g + factor_blue * b;
			buf[RLAYER][i] = round_to_WORD(payne_compute(r, x, B, expD, LP, SP, HP, stretchtype) * norm);
			buf[GLAYER][i] = round_to_WORD(payne_compute(g, x, B, expD, LP, SP, HP, stretchtype) * norm);
			buf[BLAYER][i] = round_to_WORD(payne_compute(b, x, B, expD, LP, SP, HP, stretchtype) * norm);
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {

			double x = buf[RLAYER][i]*invnorm - BP;
			buf[RLAYER][i] = round_to_WORD(payne_compute(x, x, B, expD, LP, SP, HP, stretchtype) * norm);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int paynelut_float(fits *fit, double B, double D, double LP, double SP, double HP, double BP, gboolean human_luminance, int stretchtype) {
	float *buf[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };
	double expD = exp(D) - 1;

	double factor_red = human_luminance ? 0.2126 : 0.3333;
	double factor_green = human_luminance ? 0.7152 : 0.3333;
	double factor_blue = human_luminance ? 0.0722 : 0.3333;

	//Set up expressions
	payne_precompute(B, expD, LP, SP, HP, stretchtype);

	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {
			double r = buf[RLAYER][i] - BP;
			double g = buf[GLAYER][i] - BP;
			double b = buf[BLAYER][i] - BP;
			double x = factor_red * r + factor_green * g + factor_blue * b;
			buf[RLAYER][i] = payne_compute(r, x, B, expD, LP, SP, HP, stretchtype);
			buf[GLAYER][i] = payne_compute(g, x, B, expD, LP, SP, HP, stretchtype);
			buf[BLAYER][i] = payne_compute(b, x, B, expD, LP, SP, HP, stretchtype);
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->ry * 16)
#endif
		for (i = 0; i < n; i++) {
			double x = buf[RLAYER][i] - BP;
			buf[RLAYER][i] = payne_compute(x, x, B, expD, LP, SP, HP, stretchtype);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int paynelut(fits *fit, double beta, double intensity, double lower, double shoulder, double headroom, double offset, gboolean human_luminance, gboolean inverted) {
	int stretchtype;
	if (!inverted) {
		stretchtype = STRETCH_PAYNE_NORMAL;
	} else {
		stretchtype = STRETCH_PAYNE_INVERSE;
	}
	if (fit->type == DATA_USHORT)
		return paynelut_ushort(fit, beta, intensity, lower, shoulder, headroom, offset, human_luminance, stretchtype);
	if (fit->type == DATA_FLOAT)
		return paynelut_float(fit, beta, intensity, lower, shoulder, headroom, offset, human_luminance, stretchtype);
	return 1;
}

static void apply_payne_changes() {
	gboolean status = (payne_b != 0.5) || (payne_D != 0.0) || (payne_SP != 0.0) || (payne_HP != 1.0) || (payne_LP != 0.0) || !payne_rgb_space || payne_inverse;
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
	GtkToggleButton *toggle_rgb = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace"));

	payne_startup();
	payne_D = 0.0;
	payne_b = 0.5;
	payne_SP = 0.0;
	payne_HP = 1.0;
	payne_LP = 0.0;
	payne_black_value = 0.0;
	payne_rgb_space = TRUE;
//	payne_show_preview = TRUE;

	set_notify_block(TRUE);
	gtk_toggle_button_set_active(toggle_rgb, payne_rgb_space);
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
	GtkToggleButton *toggle_rgb = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_payneRGBspace"));
	payne_D = 0.0;
	payne_b = 0.5;
	payne_black_value = 0.0;
	payne_SP = 0.0;
	payne_HP = 1.0;
	payne_rgb_space = TRUE;

	set_notify_block(TRUE);
	gtk_toggle_button_set_active(toggle_rgb, payne_rgb_space);
	gtk_spin_button_set_value(spin_d, payne_D);
	gtk_spin_button_set_value(spin_b, payne_b);
	gtk_spin_button_set_value(spin_SP, payne_SP);
	gtk_spin_button_set_value(spin_HP, payne_HP);
	gtk_spin_button_set_value(spin_black_p, payne_black_value);
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

void on_payne_RGBspace_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	payne_rgb_space = gtk_toggle_button_get_active(togglebutton);
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

void on_payne_inverse_toggled(GtkToggleButton *button, gpointer user_data) {
	if (payne_inverse == TRUE) {
		payne_inverse = FALSE;
	} else payne_inverse = TRUE;
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = payne_update_preview;
	param->show_preview = payne_show_preview;
	notify_update((gpointer) param);
}
