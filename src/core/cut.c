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

/*
 * Non-GUI lifecycle helpers for cut_struct (intensity profile / spectrogram).
 * The GTK dialog callbacks remain in gui/cut.c.
 */

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/cut.h"

void initialize_cut_struct(cut_struct *arg) {
	arg->fit = gfit;
	arg->seq = NULL;
	arg->imgnumber = -1;
	arg->cut_start.x = -1;
	arg->cut_start.y = -1;
	arg->cut_end.x = -1;
	arg->cut_end.y = -1;
	arg->cut_wn1.x = -1;
	arg->cut_wn1.y = -1;
	arg->cut_wn2.x = -1;
	arg->cut_wn2.y = -1;
	arg->cut_measure = FALSE;
	arg->plot_as_wavenumber = FALSE;
	arg->wavenumber1 = -1;
	arg->wavenumber2 = -1;
	arg->plot_spectro_bg = FALSE;
	arg->bg_poly_order = 3;
	arg->tri = FALSE;
	arg->mode = CUT_MONO;
	arg->width = 1;
	arg->step = 1;
	arg->display_graph = TRUE;
	g_free(arg->filename);
	arg->filename = NULL;
	g_free(arg->title);
	arg->title = NULL;
	g_free(arg->user_title);
	arg->user_title = NULL;
	arg->title_has_sequence_numbers = FALSE;
	arg->save_dat = FALSE;
	arg->save_png_too = FALSE;
	arg->pref_as = gfit->keywords.wcsdata.pltsolvd;
	arg->vport = -1;
	/* gui_function(reset_cut_gui, NULL) is called by gui/cut.c when needed */
}

void free_cut_args(cut_struct *arg) {
	if (!arg) return;
	g_free(arg->filename);   arg->filename   = NULL;
	g_free(arg->title);      arg->title      = NULL;
	g_free(arg->user_title); arg->user_title = NULL;
	free(arg);
}

gboolean cut_struct_is_valid(cut_struct *arg) {
	int nb_layers = (arg->seq) ? arg->seq->nb_layers : arg->fit->naxes[2];
	if (arg->tri && arg->cfa) {
		siril_log_error(_("Error: CFA mode and tri-profile are mutually exclusive.\n"));
		return FALSE;
	}
	if (arg->tri && arg->mode == CUT_COLOR) {
		siril_log_error(_("Error: color plot and tri-profile are mutually exclusive.\n"));
		return FALSE;
	}
	if (nb_layers == 1 && arg->mode == CUT_COLOR) {
		siril_log_error(_("Error: mono image is loaded: color mode is unavailable.\n"));
		return FALSE;
	}
	if (nb_layers == 1 && arg->vport > 0) {
		siril_log_error(_("Error: mono image is loaded: colored layers cannot be specified.\n"));
		return FALSE;
	}
	if (arg->cfa && arg->mode == CUT_COLOR) {
		siril_log_error(_("Error: color plot and CFA mode are mutually exclusive.\n"));
		return FALSE;
	}
	if (arg->cut_start.x == arg->cut_end.x && arg->cut_start.y == arg->cut_end.y) {
		siril_log_error(_("Error: start and finish points are the same.\n"));
		return FALSE;
	}
	if (arg->seq && arg->seq->is_variable) {
		siril_log_error(_("Error: variable image size in sequence.\n"));
		return FALSE;
	}

	int rx = (arg->seq) ? arg->seq->rx : gfit->rx;
	int ry = (arg->seq) ? arg->seq->ry : gfit->ry;

	if (arg->cut_wn1.x > -1.0 && arg->cut_wn1.x == arg->cut_wn2.x && arg->cut_wn1.y == arg->cut_wn2.y) {
		siril_log_error(_("Error: wavenumber points are the same.\n"));
		return FALSE;
	}
	if (arg->cut_wn1.x > -1.0 && arg->wavenumber1 == -1.0) {
		siril_log_error(_("Error: wavenumber / wavelength for point 1 is not set.\n"));
		return FALSE;
	}
	if (arg->cut_wn2.x > -1.0 && arg->wavenumber2 == -1.0) {
		siril_log_error(_("Error: wavenumber / wavelength for point 2 is not set.\n"));
		return FALSE;
	}
	if (arg->wavenumber1 >= 0.0 && (arg->cut_wn1.x < 0.0 || arg->cut_wn1.y < 0.0)) {
		siril_log_error(_("Error: wavenumber / wavelength set for point 1 but corresponding location is not set.\n"));
		return FALSE;
	}
	if (arg->wavenumber2 >= 0.0 && (arg->cut_wn2.x < 0.0 || arg->cut_wn2.y < 0.0)) {
		siril_log_error(_("Error: wavenumber / wavelength set for point 2 but corresponding location is not set.\n"));
		return FALSE;
	}
	if ((arg->cut_wn1.x < 0.0 || arg->cut_wn2.y < 0.0 || arg->cut_wn2.x < 0.0 || arg->cut_wn2.y < 0.0) &&
		(!(arg->cut_wn1.x == -1.0 && arg->cut_wn1.y == -1.0 && arg->cut_wn2.x == -1.0 && arg->cut_wn2.y == -1.0))) {
		siril_log_error(_("Error: wavenumber / wavelength point outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->cut_wn1.x > rx - 1 || arg->cut_wn2.x > rx - 1 || arg->cut_wn1.y > ry - 1 || arg->cut_wn2.y > ry - 1) {
		siril_log_error(_("Error: wavenumber / wavelength point outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->cut_start.x > rx - 1 || arg->cut_start.y > ry - 1 || arg->cut_end.x > rx - 1 || arg->cut_end.y > ry - 1) {
		siril_log_error(_("Error: profile line endpoint outside image dimensions.\n"));
		return FALSE;
	}
	if (arg->step != 1 && !arg->tri) {
		siril_log_error(_("Error: spacing is set but tri-profile mode is not selected.\n"));
		return FALSE;
	}
	if (arg->vport < 0 || arg->vport > nb_layers) {
		siril_log_error(_("Error: layer out of range.\n"));
		return FALSE;
	}
	if (arg->bg_poly_order < 1 || arg->bg_poly_order > 6) {
		siril_log_error(_("Error: background removal polynomial degree should be >= 1 and <= 6.\n"));
		return FALSE;
	}
	return TRUE;
}
