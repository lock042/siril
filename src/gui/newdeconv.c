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

#include <fftw3.h>
#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/command.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/sequence_list.h"
#include "gui/message_dialog.h"
#include "gui/registration_preview.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "filters/deconvolution/deconvolution.h"
#include "filters/synthstar.h"

estk_data args = { 0 };
int psftype = 0;
float *kernel;

void reset_conv_args() {
	siril_debug_print("Resetting deconvolution args\n");

	// Basic image and kernel parameters
	args.fdata = NULL;
	args.rx = 0;
	args.ry = 0;
	args.nchans = 1;
	args.ks = 15;

	// Process parameters
	args.blindtype = 0;

	// l0 Descent Kernel Estimation parameters
	args.lambda = 3000.f;
	args.lambda_ratio = 1/1.1f;
	args.lambda_min = 1e-2f; // (was 1e-3 but disagrees with header)
	args.gamma = 20.f;
	args.iterations = 2;
	args.multiscale = FALSE;
	args.scalefactor = 0.5f;
	args.kernel_threshold_max = 0.f;
	args.remove_isolated = FALSE;
	args.better_kernel = FALSE;
	args.upscaleblur = 0.f;
	args.downscaleblur = 1.6f;
	args.k_l1 = 0.5f;

	// Spectral irregularity kernel estimation parameters
	args.ninner = 300;
	args.ntries = 30;
	args.nouter = 3;
	args.compensationfactor = 2.1f;
	args.medianfilter = 1.f; // check
	args.intermediatedeconvolutionweight = 3000.f;
	args.finaldeconvolutionweight = 3000.f;

	// Synthetic kernel parameters
	args.psf_fwhm = 1.f;
	args.psf_beta = 4.5f;
	args.psf_angle = 0.f;
	args.psf_ratio = 1.f;

	// Non-blind deconvolution parameters
	args.alpha = 3000.f;
}

void reset_conv_kernel() {
	if (kernel != NULL) {
		free(kernel);
		kernel = NULL;
	}
}

void reset_conv_controls() {
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("bdeconv_profile")), 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_multiscale")), args.multiscale);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_lambda")), args.lambda);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_lambdaratio")), args.lambda_ratio);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_lambdamin")), args.lambda_min);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_gamma")), args.gamma);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_iters")), args.iterations);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_scalefactor")), args.scalefactor);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_kthresh")), args.kernel_threshold_max);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_removeisolated")), args.remove_isolated);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_betterkernel")), args.better_kernel);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_upsampleblur")), args.upscaleblur);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_downsampleblur")), args.downscaleblur);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_kl1")), args.k_l1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_alpha")), args.alpha);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfwhm")), args.psf_fwhm);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfbeta")), args.psf_beta);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfangle")), args.psf_angle);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfratio")), args.psf_ratio);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ninner")), args.ninner);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ntries")), args.ntries);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_nouter")), args.nouter);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ncomp")), args.compensationfactor);
}

void reset_conv_controls_and_args() {
	reset_conv_args();
	reset_conv_controls();
	reset_conv_args();
}

void on_bdeconv_psfblind_toggled(GtkToggleButton *button, gpointer user_data) {
	psftype = 0;
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_psfcontrols")), FALSE);
	if (args.blindtype == 0) {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_gfcontrols")), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_l0controls")), FALSE);
	} else if (args.blindtype == 1) {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_gfcontrols")), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_l0controls")), TRUE);
	}
}

void on_bdeconv_ks_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.ks = gtk_spin_button_get_value(button);
	if(!(args.ks % 2)) {
		args.ks++;
		gtk_spin_button_set_value(button, args.ks);
	}
}

void on_bdeconv_blindtype_changed(GtkComboBox *combo, gpointer user_data) {
	args.blindtype = gtk_combo_box_get_active(combo);
	printf("combo changed\n");
	if (psftype == 0) {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_psfcontrols")), FALSE);
		switch (args.blindtype) {
			case 0:
				gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_l0controls")), FALSE);
				gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_gfcontrols")), TRUE);
				break;
			case 1:
				gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_l0controls")), TRUE);
				gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_gfcontrols")), FALSE);
				break;
		}
	}
}

void on_bdeconv_expander_activate(GtkExpander *expander, gpointer user_data) {
		gtk_window_resize(GTK_WINDOW(lookup_widget("bdeconv_dialog")), 1, 1);
}

void on_bdeconv_psfselection_toggled(GtkToggleButton *button, gpointer user_data) {
	psftype = 1;
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_psfcontrols")), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_l0controls")), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_gfcontrols")), FALSE);
}

void on_bdeconv_psfstars_toggled(GtkToggleButton *button, gpointer user_data) {
	psftype = 2;
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_psfcontrols")), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_l0controls")), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_gfcontrols")), FALSE);
}

void on_bdeconv_psfmanual_toggled(GtkToggleButton *button, gpointer user_data) {
	psftype = 3;
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_psfcontrols")), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_l0controls")), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("bdeconv_gfcontrols")), FALSE);
}

void on_bdeconv_psfwhm_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.psf_fwhm = gtk_spin_button_get_value(button);
}

void on_bdeconv_psfratio_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.psf_ratio = gtk_spin_button_get_value(button);
}

void on_bdeconv_psfbeta_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.psf_beta = gtk_spin_button_get_value(button);
}

void on_bdeconv_psfangle_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.psf_angle = gtk_spin_button_get_value(button);
}

void on_bdeconv_gamma_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.gamma = gtk_spin_button_get_value(button);
}

void on_bdeconv_iters_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.iterations = gtk_spin_button_get_value(button);
}

void on_bdeconv_lambdaratio_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.lambda_ratio = gtk_spin_button_get_value(button);
}

void on_bdeconv_lambdamin_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.lambda_min = gtk_spin_button_get_value(button);
}

void on_bdeconv_scalefactor_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.scalefactor = gtk_spin_button_get_value(button);
}

void on_bdeconv_upsampleblur_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.upscaleblur = gtk_spin_button_get_value(button);
}

void on_bdeconv_downsampleblur_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.downscaleblur = gtk_spin_button_get_value(button);
}

void on_bdeconv_kthresh_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.kernel_threshold_max = gtk_spin_button_get_value(button);
}

void on_bdeconv_kl1_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.k_l1 = gtk_spin_button_get_value(button);
}

void on_bdeconv_ninner_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.ninner = gtk_spin_button_get_value(button);
}

void on_bdeconv_ntries_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.ntries = gtk_spin_button_get_value(button);
}

void on_bdeconv_nouter_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.nouter = gtk_spin_button_get_value(button);
}

void on_bdeconv_ncomp_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.compensationfactor = gtk_spin_button_get_value(button);
}

void on_bdeconv_gflambda_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.intermediatedeconvolutionweight = gtk_spin_button_get_value(button);
	args.finaldeconvolutionweight = gtk_spin_button_get_value(button);
	args.lambda = 1.f / gtk_spin_button_get_value(button);
}

void on_bdeconv_alpha_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.alpha = 1.f / gtk_spin_button_get_value(button);
}

void on_bdeconv_close_clicked(GtkButton *button, gpointer user_data) {
	reset_conv_controls_and_args();
	siril_close_dialog("bdeconv_dialog");
}

void on_bdeconv_reset_clicked(GtkButton *button, gpointer user_data) {
	reset_conv_controls_and_args();
}

void on_bdeconv_dialog_show(GtkWidget *widget, gpointer user_data) {
	reset_conv_controls_and_args();
}

void on_bdeconv_setlambdafromnoise_clicked(GtkButton *button, gpointer user_data) {

}

void on_bdeconv_apply_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	gchar *wisdom_file = g_build_filename(g_get_user_cache_dir(), "siril_fftw.wisdom", NULL);
    if (fftwf_import_wisdom_from_filename(wisdom_file) == 1) {
        printf("Siril wisdom imported successfully...\n");
	} else if (fftwf_import_system_wisdom() == 1) {
        printf("System wisdom imported successfully...\n");
	} else {
        printf("No wisdom found to import...\n");
	}
	undo_save_state(&gfit, _("Deconvolution"));
	gboolean out_16bit = (com.pref.force_16bit) ? TRUE : FALSE;
	unsigned ndata = gfit.rx * gfit.ry * gfit.naxes[2];

	args.fdata = malloc(ndata * sizeof(float));
	if (gfit.type == DATA_FLOAT)
		memcpy(args.fdata, gfit.fdata, ndata * sizeof(float));
	else {
		float invnorm = 1.f / USHRT_MAX_SINGLE;
		for (size_t i = 0 ; i < ndata ; i++) {
			args.fdata[i] = (float) gfit.data[i] * invnorm;
		}
	}
	args.rx = gfit.rx;
	args.ry = gfit.ry;
	args.nchans = gfit.naxes[2];
	float *kernel = NULL;
	switch (psftype) {
		case 0:
			reset_conv_kernel();
			args.better_kernel = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_betterkernel")));
			args.remove_isolated = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_removeisolated")));
			args.multiscale = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_multiscale")));

			siril_log_message(_("Starting kernel estimation...\n"));
			switch(args.blindtype) {
				case 0:
					kernel = gf_estimate_kernel(&args);
					break;
				case 1:
					kernel = estimate_kernel(&args);
					break;
			}
			siril_log_message(_("Kernel estimation complete.\n"));
			break;
		case 1:
			break;
		case 2:
			break;
		case 3:
			kernel = (float*) calloc(args.ks * args.ks, sizeof(float));
			// Should switch here and do gaussian, moffat, disc
			makegaussian(kernel, args.ks, args.psf_fwhm, 1.f, 0.f, 0.f);
			break;
	}
	if (kernel == NULL) {
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error: no kernel defined"),
				_("Select blind deconvolution or define a kernel from selection or psf parameters"));
		return;
	}
	else {
		siril_log_message(_("Starting non-blind deconvolution...\n"));
#ifdef SIRIL_OUTPUT_DEBUG
		for (int i = 0; i < args.ks; i++) {
			for (int j = 0 ; j < args.ks ; j++) {
				printf("%f ",kernel[i + j * args.ks]);
			}
			printf("\n");
		}
#endif
		split_bregman(args.fdata, args.rx, args.ry, args.nchans, kernel, args.ks, args.alpha);
//		richardson_lucy(args.fdata, args.rx,args.ry, args.nchans, kernel, args.ks, args.alpha, 8);
//		stochastic(args.fdata, args.rx, args.ry, args.nchans, kernel, args.ks, 0.0008f);
	}
	if (kernel) {
		free (kernel);
		kernel = NULL;
	}
	if (gfit.type == DATA_FLOAT)
		memcpy(gfit.fdata, args.fdata, ndata * sizeof(float));
	else {
		for (size_t i = 0 ; i < ndata ; i++) {
			gfit.data[i] = roundf_to_WORD(args.fdata[i] * USHRT_MAX_SINGLE);
		}
	}

	if (fftwf_export_wisdom_to_filename(wisdom_file) == 1) {
        printf("Siril wisdom updated successfully...\n");
	} else {
        printf("Siril wisdom update failed...\n");
	}
    g_free(wisdom_file);

	free(args.fdata);
	args.fdata = NULL;

	update_zoom_label();
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);

}


