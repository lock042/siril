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

#include "core/siril.h"
#include "core/command.h"
#include "algos/demosaicing.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/sequence_list.h"
#include "gui/PSF_list.h"
#include "gui/message_dialog.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "filters/deconvolution/deconvolution.h"

estk_data args = { 0 };

void reset_args() {
	args.fdata = NULL;
	args.rx = 0;
	args.ry = 0;
	args.ks = com.kernelsize;
	args.nchans = 1;
	args.lambda = 4e-3f;
	args.lambda_ratio = 1/1.1f;
	args.lambda_min = 1e-2f;
	args.gamma = 20.f;
	args.iterations = 2;
	args.multiscale = TRUE;
	args.scalefactor = 0.5f;
	args.kernel_threshold_max = 0.f;
	args.remove_isolated = FALSE;
	args.better_kernel = FALSE;
	args.upscaleblur = 0.f;
	args.downscaleblur = 1.6f;
	args.k_l1 = 0.5f;
}

void reset_kernel() {
	if (com.kernel != NULL) {
		free(com.kernel);
		com.kernel = NULL;
	}
	com.kernelsize = 0;
}

void reset_controls() {
//	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("merge_cfa_pattern")), 0);
//	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("merge_cfa_seqapply")), FALSE);
//	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("merge_cfa_seq_controls")), FALSE);
}

void reset_controls_and_args() {
	reset_controls();
	reset_args();
}

void on_bdeconv_close_clicked(GtkButton *button, gpointer user_data) {
	reset_controls_and_args();
	siril_close_dialog("bdeconv_dialog");
}

void on_bdeconv_reset_clicked(GtkButton *button, gpointer user_data) {
	reset_controls_and_args();
}

void on_bdeconv_show(GtkWidget *widget, gpointer user_data) {
	reset_controls_and_args();
}

void on_bdeconv_apply_clicked(GtkButton *button, gpointer user_data) {
	args.fdata = gfit.fdata;
	args.rx = gfit.rx;
	args.ry = gfit.ry;
	args.nchans = gfit.naxes[2];
	float *kernel = NULL;
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_isblind")))) {
		reset_kernel();
		estimate_kernel(&args, kernel);
		com.kernel = kernel;
		com.kernelsize = args.ks;
	}
	if (com.kernel == NULL)
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: no kernel defined"),
				_("Select blind deconvolution or define a kernel from selection or psf parameters"));
	else
		split_bregman(args.fdata, args.rx, args.ry, args.nchans, com.kernel, com.kernelsize, args.lambda);

}


