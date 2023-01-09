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

#ifdef _WIN32
#include <windows.h> // for MAX_PATH
#endif
#include <fftw3.h>
#include <math.h>
#include <gdk/gdk.h>
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
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/newdeconv.h"
#include "gui/newdeconv_fit.h"
#include "filters/deconvolution/deconvolution.h"
#include "filters/deconvolution/rlstrings.h"
#include "filters/deconvolution/chelperfuncs.h"
#include "filters/synthstar.h"
#include "filters/starnet.h"
#include "algos/statistics.h"
#include "algos/PSF.h"
#include "io/sequence.h"
#include "io/ser.h"

gboolean aperture_warning_given = FALSE;

estk_data args = { 0 };
static GtkWidget *drawingPSF = NULL;
static void DrawPSF();
static cairo_surface_t *surface = NULL;

//set below flag to zero to avoid kernel printout to stdout
#define DEBUG_PSF 1

void reset_conv_args(estk_data* args) {
	siril_debug_print("Resetting deconvolution args\n");

	// Basic image and kernel parameters
	args->psftype = PSF_BLIND;
	the_fit = &gfit;
	args->fdata = NULL;
	args->rx = 0;
	args->ry = 0;
	args->nchans = 1;
	args->ks = 15;

	// Process parameters
	args->blindtype = BLIND_L0;

	// l0 Descent Kernel Estimation parameters
	args->lambda = 1.f / 3000.f;
	args->lambda_ratio = 1/1.1f;
	args->lambda_min = 1e-3f;
	args->gamma = 20.f;
	args->iterations = 3;
	args->multiscale = FALSE;
	args->scalefactor = 0.5f;
	args->kernel_threshold_max = 0.f;
	args->remove_isolated = FALSE;
	args->better_kernel = TRUE;
	args->upscaleblur = 0.f;
	args->downscaleblur = 1.6f;
	args->k_l1 = 0.5f;

	// Spectral irregularity kernel estimation parameters
	args->ninner = 300;
	args->ntries = 30;
	args->nouter = 3;
	args->compensationfactor = 2.1f;
	args->medianfilter = 1.f; // check
	args->intermediatedeconvolutionweight = 3000.f;
	args->finaldeconvolutionweight = 3000.f;

	// Synthetic kernel parameters
	args->profile = PROFILE_GAUSSIAN;
	args->psf_fwhm = 3.0f;
	args->psf_beta = 4.5f;
	args->psf_angle = 0.f;
	args->psf_ratio = 1.f;
	args->airy_wl = 525.f;
	args->airy_pixelsize = 0.1f;
	args->airy_fl = 1.f;
	args->airy_diameter = 1.f;
	args->airy_obstruction = 0.f;

	// Non-blind deconvolution parameters
	args->nonblindtype = DECONV_RL;
	args->finaliters = 10;
	args->alpha = 1.f / 3000.f;
	args->stopcriterion = 0.002f;
	args->rl_method = 1;
	args->stepsize = 0.0003f;
	args->regtype = REG_TV_GRAD;
}

void reset_conv_kernel() {
	if (com.kernel != NULL) {
		free(com.kernel);
		com.kernel = NULL;
	}
}

void reset_conv_controls() {
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("bdeconv_profile")), args.profile);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("bdeconv_blindtype")), args.blindtype);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("bdeconv_nonblindtype")), args.nonblindtype);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_multiscale")), args.multiscale);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_lambdaratio")), args.lambda_ratio);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_lambdamin")), args.lambda_min);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_gamma")), args.gamma);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_iters")), args.iterations);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_scalefactor")), args.scalefactor);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_kthresh")), args.kernel_threshold_max);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_betterkernel")), args.better_kernel);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_upsampleblur")), args.upscaleblur);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_downsampleblur")), args.downscaleblur);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_kl1")), args.k_l1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_alpha")), 1.f / args.alpha);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfwhm")), args.psf_fwhm);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfbeta")), args.psf_beta);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfratio")), args.psf_ratio);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ninner")), args.ninner);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ntries")), args.ntries);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_nouter")), args.nouter);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_finaliters")), args.finaliters);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ncomp")), args.compensationfactor);
}

void reset_conv_controls_and_args() {
	if (!get_thread_run())
		reset_conv_args(&args);
	reset_conv_controls();
}

void on_bdeconv_psfblind_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_BLIND;
	gtk_widget_set_visible(lookup_widget("bdeconv_blindcontrols"), TRUE);
	gtk_widget_set_visible(lookup_widget("bdeconv_psfcontrols"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_starpsf_details"), FALSE);
	if (args.blindtype == BLIND_SI) {
		gtk_widget_set_visible(lookup_widget("bdeconv_gfcontrols"), TRUE);
		gtk_widget_set_visible(lookup_widget("bdeconv_l0controls"), FALSE);
	} else if (args.blindtype == BLIND_L0) {
		gtk_widget_set_visible(lookup_widget("bdeconv_gfcontrols"), FALSE);
		gtk_widget_set_visible(lookup_widget("bdeconv_l0controls"), TRUE);
	}
}

void on_bdeconv_ks_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.ks = gtk_spin_button_get_value(button);
	// Prevent setting ks to an even value.
	if(!(args.ks % 2)) {
		args.ks++;
		gtk_spin_button_set_value(button, args.ks);
	}
}

void on_bdeconv_blindtype_changed(GtkComboBox *combo, gpointer user_data) {
	args.blindtype = gtk_combo_box_get_active(combo);
	if (args.psftype == PSF_BLIND) {
		gtk_widget_set_visible(lookup_widget("bdeconv_psfcontrols"), FALSE);
		switch (args.blindtype) {
			case 0:
				gtk_widget_set_visible(lookup_widget("bdeconv_l0controls"), FALSE);
				gtk_widget_set_visible(lookup_widget("bdeconv_gfcontrols"), TRUE);
				args.intermediatedeconvolutionweight = 15.f;
				args.finaldeconvolutionweight = 15.f;
				args.lambda = 1.f / 15.f;
				gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_gflambda")), 15.);
				break;
			case 1:
				gtk_widget_set_visible(lookup_widget("bdeconv_l0controls"), TRUE);
				gtk_widget_set_visible(lookup_widget("bdeconv_gfcontrols"), FALSE);
				args.lambda = 1.f / 3000.f;
				gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_gflambda")), 3000.);
				break;
		}
	}
}

void on_bdeconv_nonblindtype_changed(GtkComboBox *combo, gpointer user_data) {
	args.nonblindtype = gtk_combo_box_get_active(combo);
	switch (args.nonblindtype) {
		case 0:
			args.finaliters = 1; // Default niters for Split Bregman
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_finaliters")), (double) args.finaliters);
			args.alpha = 3000.f;
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_alpha")), (double) args.alpha);
			gtk_widget_set_visible(lookup_widget("regul_label"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_rl_regularization"), FALSE);
			gtk_widget_set_visible(lookup_widget("algo_method_label"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_rl_method"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stopping_toggle"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stopcriterion"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stepsize"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_steplabel"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_finaliters"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_iterlabel"), TRUE);
			break;
		case 1:
			args.finaliters = 10; // Default niters for RL. Reduce this if using the multiplicative
			// algorithm in order to avoid burning holes round your stars!
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_finaliters")), (double) args.finaliters);
			args.alpha = 3000.f;
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_alpha")), (double) args.alpha);
			gtk_widget_set_visible(lookup_widget("regul_label"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_rl_regularization"), TRUE);
			gtk_widget_set_visible(lookup_widget("algo_method_label"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_rl_method"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stopping_toggle"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stopcriterion"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stepsize"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_steplabel"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_finaliters"), TRUE);
			gtk_widget_set_visible(lookup_widget("bdeconv_iterlabel"), TRUE);
			break;
		case 2:
			args.alpha = 500.f;
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_alpha")), (double) args.alpha);
			gtk_widget_set_visible(lookup_widget("regul_label"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_rl_regularization"), FALSE);
			gtk_widget_set_visible(lookup_widget("algo_method_label"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_rl_method"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stopping_toggle"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stopcriterion"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_stepsize"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_steplabel"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_finaliters"), FALSE);
			gtk_widget_set_visible(lookup_widget("bdeconv_iterlabel"), FALSE);
			break;
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_finaliters")), args.finaliters);
}

void on_bdeconv_expander_activate(GtkExpander *expander, gpointer user_data) {
		gtk_window_resize(GTK_WINDOW(lookup_widget("bdeconv_dialog")), 1, 1);
// This callback is just to prevent excessive blank space after shrinking the expander
}

// This type of PSF generation doesn't exist yet. Maybe don't need it with PSF from stars working well.
/*void on_bdeconv_psfselection_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_SELECTION;
	gtk_widget_set_visible(lookup_widget("bdeconv_psfcontrols"), TRUE);
	gtk_widget_set_visible(lookup_widget("bdeconv_l0controls"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_gfcontrols"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_blindcontrols"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_starpsf_details"), FALSE);
}
*/

void on_airy_diameter_value_changed(GtkSpinButton *button, gpointer user_data) {
		GtkWidget *the_button = lookup_widget("airy_diameter");
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:@theme_color;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(the_button);

		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
}

void on_airy_fl_value_changed(GtkSpinButton *button, gpointer user_data) {
		GtkWidget *the_button = lookup_widget("airy_fl");
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:@theme_color;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(the_button);

		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
}

void on_airy_pixelsize_value_changed(GtkSpinButton *button, gpointer user_data) {
		GtkWidget *the_button = lookup_widget("airy_pixelsize");
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:@theme_color;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(the_button);

		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
}

static void initialize_airy_parameters() {
	// Get initial stab at parameters for Airy function from FITS header. Not essential they be correct as they are user-editable in the UI.
	args.airy_fl = the_fit->focal_length;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("airy_fl")), args.airy_fl);
	if (the_fit->focal_length < 10.) {
		GtkWidget *button = lookup_widget("airy_fl");
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:salmon;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(button);
		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
	}
	args.airy_diameter = the_fit->aperture;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("airy_diameter")), args.airy_diameter);
	if (the_fit->aperture < 10.) {
		GtkWidget *button = lookup_widget("airy_diameter");
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:salmon;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(button);
		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
	}
	args.airy_pixelsize = the_fit->pixel_size_x;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("airy_pixelsize")), args.airy_pixelsize);
	if (the_fit->pixel_size_x <=1.) {
		GtkWidget *button = lookup_widget("airy_pixelsize");
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:salmon;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(button);
		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
	}
	args.airy_wl = 525.f; // Approximation of the average wavelength of light accepted by an OSC camera covering 400-656nm. User editable to give precise values e.g. for narrowband.
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("airy_wl")), args.airy_wl);
}

void on_bdeconv_psfmanual_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_MANUAL;
	gtk_widget_set_visible(lookup_widget("bdeconv_psfcontrols"), TRUE);
	gtk_widget_set_visible(lookup_widget("bdeconv_l0controls"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_gfcontrols"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_blindcontrols"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_starpsf_details"), FALSE);
}

void on_bdeconv_profile_changed(GtkComboBox *combo, gpointer user_data) {
	int profile = gtk_combo_box_get_active(combo);
	if (profile != PROFILE_AIRY) {
		gtk_widget_set_visible(lookup_widget("bdeconv_manual_stars"), TRUE);
		gtk_widget_set_visible(lookup_widget("bdeconv_manual_airy"), FALSE);
		gtk_widget_set_sensitive(lookup_widget("bdeconv_psfangle"), (profile != PROFILE_MOFFAT && profile!= PROFILE_GAUSSIAN));
		gtk_widget_set_sensitive(lookup_widget("bdeconv_psfratio"), (profile != PROFILE_MOFFAT && profile != PROFILE_GAUSSIAN));
		gtk_widget_set_sensitive(lookup_widget("bdeconv_psfbeta"), profile == PROFILE_MOFFAT);
	} else {
		gtk_widget_set_visible(lookup_widget("bdeconv_manual_stars"), FALSE);
		gtk_widget_set_visible(lookup_widget("bdeconv_manual_airy"), TRUE);
		if (!aperture_warning_given) {
			aperture_warning_given = TRUE;
			if (args.airy_diameter < 10.f) {
				siril_log_color_message(_("Warning: telescope aperture obtained from FITS header data may be incorrect.\n"), "salmon");
			}
			if (args.airy_fl < 10.f) {
				siril_log_color_message(_("Warning: telescope focal length obtained from FITS header data may be incorrect.\n"), "salmon");
			}
			if (args.airy_pixelsize <=1.f) {
				siril_log_color_message(_("Warning: sensor pixel size obtained from FITS header data may be incorrect.\n"), "salmon");
			}
		}
	}
}

void on_bdeconv_rl_method_changed(GtkComboBox *combo, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("bdeconv_stepsize"), gtk_combo_box_get_active(combo) == 1);
}

void on_bdeconv_psfprevious_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_PREVIOUS;
}

void on_bdeconv_close_clicked(GtkButton *button, gpointer user_data) {
	if (sequence_is_running == 0)
		reset_conv_controls_and_args();
	siril_close_dialog("bdeconv_dialog");
}

void on_bdeconv_reset_clicked(GtkButton *button, gpointer user_data) {
	reset_conv_controls_and_args();
}

static void calculate_parameters() {
	if (com.stars && com.stars[0]) {
		int i = 0;
		double FWHMx = 0.0, FWHMy = 0.0, beta = 0.0, angle = 0.0;
		gboolean unit_is_arcsec;
		int n = 0, layer;
		starprofile profiletype;

		while (com.stars[i]) {
			double fwhmx, fwhmy;
			char *unit;
			gboolean is_as = get_fwhm_as_arcsec_if_possible(com.stars[i], &fwhmx, &fwhmy, &unit);
			if (i == 0) {
				unit_is_arcsec = is_as;
				layer = com.stars[i]->layer;
				profiletype = com.stars[i]->profile;
			}
			else if (is_as != unit_is_arcsec) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
						_("Stars FWHM must have the same units."));
				return;
			}
			else if (layer != com.stars[i]->layer ) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
						_("Stars properties must all be computed on the same layer"));
				return;
			}
			else if (profiletype != com.stars[i]->profile) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
						_("Stars must all be modeled with the same profile type"));
				return;
			}
			if (!com.stars[i]->has_saturated) {
				FWHMx += fwhmx;
				beta += com.stars[i]->beta;
				FWHMy += fwhmy;
				angle += com.stars[i]->angle;
				n++;
				}
			i++;
		}
		if (i <= 0 || n <= 0) return;
		/* compute average */
		args.psf_fwhm = (float) FWHMx / (float) n;
		FWHMy = (float) FWHMy / (float) n;
		args.psf_beta = (float) beta / (float) n;
		args.psf_ratio = args.symkern ? 1.f : args.psf_fwhm / FWHMy;
		if (unit_is_arcsec) {
			double bin_X = the_fit->unbinned ? (double) the_fit->binning_x : 1.0;
			double conversionfactor = (((3600.0 * 180.0) / G_PI) / 1.0E3 * (double)the_fit->pixel_size_x / the_fit->focal_length) * bin_X;
			args.psf_fwhm /= (float) conversionfactor;
		}
		args.psf_angle = (float) angle / (float) n;
		if (com.stars[0]->profile == PSF_GAUSSIAN)
			gtk_label_set_text(GTK_LABEL(lookup_widget("bdeconv_starprofile_text")), _("Gaussian"));
		else
			gtk_label_set_text(GTK_LABEL(lookup_widget("bdeconv_starprofile_text")), _("Moffat"));
		char temp[10];
		sprintf(temp, "%.2f", args.psf_fwhm);
		gtk_label_set_text(GTK_LABEL(lookup_widget("bdeconv_starfwhm_text")), temp);
		sprintf(temp, "%.2f", args.psf_ratio);
		gtk_label_set_text(GTK_LABEL(lookup_widget("bdeconv_starratio_text")), temp);
		sprintf(temp, "%.2f", args.psf_angle);
		gtk_label_set_text(GTK_LABEL(lookup_widget("bdeconv_starangle_text")), temp);
		sprintf(temp, "%.2f", args.psf_beta);
		gtk_label_set_text(GTK_LABEL(lookup_widget("bdeconv_starbeta_text")), temp);
		return;
	} else {
		gtk_label_set_text(GTK_LABEL(lookup_widget("bdeconv_starprofile_text")), _("No stars selected"));
	}
}

void on_bdeconv_psfstars_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_STARS;
	calculate_parameters();
	gtk_widget_set_visible(lookup_widget("bdeconv_psfcontrols"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_l0controls"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_gfcontrols"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_blindcontrols"), FALSE);
	gtk_widget_set_visible(lookup_widget("bdeconv_starpsf_details"), TRUE);
}

void on_bdeconv_dialog_show(GtkWidget *widget, gpointer user_data) {
	reset_conv_controls_and_args();
	calculate_parameters();
	initialize_airy_parameters();
	control_window_switch_to_tab(OUTPUT_LOGS);

}

int save_kernel(gchar* filename) {
	int retval = 0;
	fits *fit = NULL;
	if ((retval = new_fit_image_with_data(&fit, args.ks, args.ks, 1, DATA_FLOAT, com.kernel)))
		return retval;
#ifdef HAVE_LIBTIFF
	retval = savetif(filename, fit, 32, "Saved Siril deconvolution kernel", NULL, FALSE);
#else
	siril_log_color_message(_("This copy of Siril was compiled without libtiff support: saving kernel at reduced precision as a 16-bit greyscale PGM file.\n"), "red");
	retval = saveNetPBM(filename, fit);
#endif
	return retval;
}

int load_kernel(gchar* filename) {
	int retval = 0;
	int orig_size;
	fits fit = { 0 };
	if ((retval = read_single_image(filename, &fit, NULL, FALSE, NULL, FALSE, TRUE)))
		goto ENDSAVE;
	if (fit.rx != fit.ry){
		retval = 1;
		char *msg = siril_log_color_message(_("Error: kernel file does not contain a square kernel. Cannot load this file.\n"), "red");
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong kernel size"), msg);
		goto ENDSAVE;
	}
	if (!(fit.rx % 2)) {
		com.kernelsize = fit.rx - 1;
		orig_size = fit.rx;
		siril_log_color_message(_("Warning: kernel file is even (%d x %d). Kernels should always be odd. Cropping by 1 pixel in each direction. "
				"This may not produce optimum results.\n"), "salmon", fit.rx, fit.rx);
	} else {
		com.kernelsize = fit.rx;
		orig_size = com.kernelsize;
	}
	if (!com.headless)
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ks")), com.kernelsize);

	com.kernel = (float*) malloc(com.kernelsize * com.kernelsize * sizeof(float));
	if (fit.type == DATA_FLOAT) {
		for (int i = 0 ; i < com.kernelsize ; i++) {
			for (int j = 0 ; j < com.kernelsize ; j++) {
			com.kernel[i + j * com.kernelsize] = fit.fdata[i + j * orig_size];
			}
		}
	} else {
		for (int i = 0 ; i < com.kernelsize ; i++) {
			for (int j = 0 ; j < com.kernelsize ; j++) {
				com.kernel[i + j * com.kernelsize] = (float) fit.data[i + j * orig_size] / USHRT_MAX_SINGLE;
			}
		}
	}
	ENDSAVE:
	clearfits(&fit);
	if (filename)
		g_free(filename);
	return retval;
}

void on_bdeconv_savekernel_clicked(GtkButton *button, gpointer user_data) {
	// Only allocate as much space for filenames as required - we determine the max pathlength
#ifndef _WIN32
	long pathmax = get_pathmax();
#else
	long pathmax = MAX_PATH;	// On Windows use of MAX_PATH is fine as it is not
								// a configurable item
#endif
	gchar filename[pathmax];
	char *imagenoext;
	char *imagenoextorig;
	gchar kernelsuffix[10] = "_PSF";
	// Initialise the filename strings as empty strings
	memset(filename, 0, sizeof(filename));
	// Set up paths and filenames
	imagenoextorig = g_path_get_basename(com.uniq->filename);
	imagenoext = g_path_get_basename(com.uniq->filename);
	for (char *c = imagenoextorig, *q = imagenoext;  *c;  ++c, ++q)
        *q = *c == ' ' ? '_' : *c;
	if (g_strcmp0(imagenoext, imagenoextorig))
		siril_log_color_message(_("Deconvolution: spaces detected in filename. These have been replaced by underscores.\n"), "salmon");
	free(imagenoextorig);
	imagenoext = g_build_filename(com.wd, imagenoext, NULL);
	imagenoext = remove_ext_from_filename(imagenoext);
	strncat(filename, imagenoext, sizeof(filename) - strlen(imagenoext));
	strncat(filename, kernelsuffix, 10);
#ifdef HAVE_LIBTIFF
	strncat(filename, ".tif", 5);
#else
	strncat(filename, ".pgm", 5);
#endif
	save_kernel(filename);
}

void on_bdeconv_filechooser_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	gchar* filename = g_strdup(gtk_file_chooser_get_filename(filechooser));
	if (filename == NULL) {
		siril_log_color_message(_("No kernel file selected.\n"), "red");
	} else {
		load_kernel(filename);
	}
	args.psftype = PSF_PREVIOUS; // Set to use previous kernel
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_psfprevious")), TRUE);
	DrawPSF();
	return;
}

gboolean deconvolve_idle(gpointer arg) {
	set_progress_bar_data("Ready.", 0.);
	if (args.fdata) {
		free(args.fdata);
		args.fdata = NULL;
	}
	update_zoom_label();
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	siril_debug_print("Deconvolve idle stopping processing thread\n");
	stop_processing_thread();
	return FALSE;
}

gboolean estimate_idle(gpointer arg) {
	if (args.psftype == PSF_BLIND || args.psftype == PSF_STARS) {
		args.psftype = PSF_PREVIOUS; // If a blind estimate is done, switch to previous kernel in order to avoid recalculating it when Apply is clicked
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_psfprevious")), TRUE);
	}
	set_progress_bar_data("Ready.", 0.);
	set_cursor_waiting(FALSE);
	siril_debug_print("Estimate idle stopping processing thread\n");
	stop_processing_thread();
	free(args.fdata);
	args.fdata = NULL;
	DrawPSF();
	return FALSE;
}

void set_estimate_params() {
	if (com.headless)
		return;
	args.ks = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ks")));
	args.gamma = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_gamma")));
	args.iterations = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_iters")));
	args.lambda_ratio = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_lambdaratio")));
	args.lambda_min = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_lambdamin")));
	args.scalefactor = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_scalefactor")));
	args.upscaleblur = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_upsampleblur")));
	args.downscaleblur = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_downsampleblur")));
	args.kernel_threshold_max = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_kthresh")));
	args.k_l1 = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_kl1")));
	args.ninner = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ninner")));
	args.ntries = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ntries")));
	args.nouter = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ntries")));
	args.compensationfactor = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_ntries")));
	args.intermediatedeconvolutionweight = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_gflambda")));
	args.finaldeconvolutionweight = args.intermediatedeconvolutionweight;
	args.lambda = 1.f / args.intermediatedeconvolutionweight;
	args.profile = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("bdeconv_profile")));
	args.psf_fwhm = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfwhm")));
	args.psf_ratio = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfratio")));
	args.psf_beta = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfbeta")));
	args.psf_angle = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_psfangle")));
	args.airy_diameter = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("airy_diameter")));
	args.airy_fl = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("airy_fl")));
	args.airy_wl = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("airy_wl")));
	args.airy_pixelsize = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("airy_pixelsize")));
	args.airy_obstruction = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("airy_obstruction"))) / 100.f;
	args.better_kernel = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_betterkernel")));
	args.multiscale = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_multiscale")));
}

void set_deconvolve_params() {
	if (com.headless)
		return;
	args.finaliters = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_finaliters")));
	args.alpha = 1.f / gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_alpha")));
	args.stopcriterion = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_stopcriterion")));
	args.stopcriterion_active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_stopping_toggle"))) ? 1 : 0;
	args.rl_method = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("bdeconv_rl_method")));
	args.stepsize = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("bdeconv_stepsize")));
	args.regtype = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("bdeconv_rl_regularization")));
}

int get_kernel() {
	int retval = 0;
	args.rx = the_fit->rx;
	args.ry = the_fit->ry;
	args.nchans = the_fit->naxes[2];
	int fftw_max_thread = com.pref.fftw_conf.multithreaded ? com.max_thread : 1;
	switch (args.psftype) {
		case PSF_BLIND: // Blind deconvolution
			reset_conv_kernel();
			switch(args.blindtype) {
				case BLIND_SI:
					com.kernel = gf_estimate_kernel(&args, fftw_max_thread);
					break;
				case BLIND_L0:
					com.kernel = estimate_kernel(&args, fftw_max_thread);
					break;
			}
			if (get_thread_run())
				siril_log_message(_("Kernel estimation complete.\n"));
			break;
		case PSF_SELECTION: // Kernel from selection
			// Not implemented yet (maybe never)
			break;
		case PSF_STARS: // Kernel from com.stars
			if (com.kernel && sequence_is_running == 1) {
// We don't want to recompute stars on a running sequence, if we are using a kernel from stars
// detected in the first image then the only sensible approach for subsequent images is to reuse
// the kernel, as there is no way of automatically re-detecting the same stars in subsequent images.
				retval = 0;
				goto END;
			}
			if (!(com.stars && com.stars[0])) {
				siril_log_color_message(_("Error: no stars detected. Run findstar or select stars using Dynamic PSF dialog first.\n"), "red");
				retval = 1;
				goto END;
			}
			reset_conv_kernel();
			calculate_parameters();
			com.kernel = (float*) calloc(args.ks * args.ks, sizeof(float));
			if (com.stars[0]->profile == PSF_GAUSSIAN)
				makegaussian(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f,args.psf_ratio, args.psf_angle);
			else
				makemoffat(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f, args.psf_beta, args.psf_ratio, args.psf_angle);

			break;
		case PSF_MANUAL: // Kernel from provided parameters
			reset_conv_kernel();
			com.kernel = (float*) calloc(args.ks * args.ks, sizeof(float));
			switch (args.profile) {
				case PROFILE_GAUSSIAN:
					makegaussian(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f, args.psf_ratio, args.psf_angle);
					break;
				case PROFILE_MOFFAT:
					makemoffat(com.kernel, args.ks, args.psf_fwhm, 1.f, +0.5f, -0.5f, args.psf_beta, args.psf_ratio, args.psf_angle);
					break;
				case PROFILE_DISK:
					makedisc(com.kernel, args.ks, (args.psf_fwhm / 2.f), 1.f, +0.5f, -0.5f);
					break;
				case PROFILE_AIRY:
					if (args.airy_fl == 1.f)
						siril_log_message(_("Warning: focal length appears likely to be incorrect. Continuing anyway...\n"));
					if (args.airy_diameter == 1.f)
						siril_log_message(_("Warning: diameter appears likely to be incorrect. Continuing anyway...\n"));
					if (args.airy_fl == 0.1f)
						siril_log_message(_("Warning: sensor pixel size appears likely to be incorrect. Continuing anyway...\n"));

					makeairy(com.kernel, args.ks, 1.f, +0.5, -0.5, args.airy_wl, args.airy_diameter, args.airy_fl, args.airy_pixelsize, args.airy_obstruction);
					break;
			}
			break;
		case PSF_PREVIOUS:
			if (com.kernel == NULL) {
				siril_log_color_message(_("Error: no previous PSF found. A blind PSF estimator or calculation of PSF from stars or manual parameters must already have been done to use the Previous PSF option.\n"), "red");
				retval = 1;
				goto END;
			}
			break;
	}
	if (com.kernel == NULL) {
		com.kernelsize = 0;
		siril_log_color_message(_("Error: no kernel defined. Select blind deconvolution or define a kernel from selection or psf parameters.\n"), "red");
		retval = 1;
		goto END;
	}
#if DEBUG_PSF
	for (int i = 0; i < args.ks; i++) {
		for (int j = 0; j < args.ks; j++) {
			siril_debug_print("%0.2f\t", com.kernel[i * args.ks + j]);
		}
		siril_debug_print("\n");
	}
#endif
	com.kernelsize = (!com.kernel) ? 0 : args.ks;
END:
	return retval;
}

gpointer estimate_only(gpointer p) {
	if (p != NULL) {
		estk_data *command_data = (estk_data *) p;
		memcpy(&args, command_data, sizeof(estk_data));
		free(command_data);
	}
	int retval = 0;
	cppmaxthreads = com.max_thread;
	cppfftwflags = com.pref.fftw_conf.strategy;
	cppfftwtimelimit = com.pref.fftw_conf.timelimit;
	cppfftwmultithreaded = com.pref.fftw_conf.multithreaded;
	set_cursor_waiting(TRUE);
	if (args.psftype == PSF_BLIND) {
		if (fftwf_import_wisdom_from_filename(com.pref.fftw_conf.wisdom_file) == 1) {
			siril_log_message(_("Siril FFT wisdom imported successfully...\n"));
		} else if (fftwf_import_system_wisdom() == 1) {
			siril_log_message(_("System FFT wisdom imported successfully...\n"));
		} else {
			siril_log_message(_("No FFT wisdom found to import...\n"));
		}
	}
	args.ndata = the_fit->rx * the_fit->ry * the_fit->naxes[2];
	args.fdata = malloc(args.ndata * sizeof(float));
	if (the_fit->type == DATA_FLOAT)
		memcpy(args.fdata, the_fit->fdata, args.ndata * sizeof(float));
	else {
		float invnorm = 1.f / USHRT_MAX_SINGLE;
		for (size_t i = 0 ; i < args.ndata ; i++) {
			args.fdata[i] = (float) the_fit->data[i] * invnorm;
		}
	}

	set_progress_bar_data(_("Starting kernel computation..."), PROGRESS_PULSATE);
	get_kernel();

	if (args.psftype == PSF_BLIND) {
		if (fftwf_export_wisdom_to_filename(com.pref.fftw_conf.wisdom_file) == 1) {
			siril_log_message(_("Siril FFT wisdom updated successfully...\n"));
		} else {
			siril_log_message(_("Siril FFT wisdom update failed...\n"));
		}
	}
	siril_log_color_message(_("Deconvolution kernel generated.\n"), "green");
	siril_add_idle(estimate_idle, NULL);
	return GINT_TO_POINTER(retval);
}

gpointer deconvolve(gpointer p) {
	if (p != NULL) {
		estk_data *command_data = (estk_data *) p;
		memcpy(&args, command_data, sizeof(estk_data));
		free(command_data);
	}

	int retval = 0;
	int fftw_max_thread = com.pref.fftw_conf.multithreaded ? com.max_thread : 1;
	cppmaxthreads = fftw_max_thread;
	cppfftwflags = com.pref.fftw_conf.strategy;
	cppfftwtimelimit = com.pref.fftw_conf.timelimit;
	cppfftwmultithreaded = com.pref.fftw_conf.multithreaded;

	set_cursor_waiting(TRUE);
	if (fftwf_import_wisdom_from_filename(com.pref.fftw_conf.wisdom_file) == 1) {
		if (sequence_is_running == 0)
			siril_log_message(_("Siril FFT wisdom imported successfully...\n"));
	} else if (fftwf_import_system_wisdom() == 1) {
		if (sequence_is_running == 0)
			siril_log_message(_("System FFT wisdom imported successfully...\n"));
	} else {
		if (sequence_is_running == 0)
			siril_log_message(_("No FFT wisdom found to import...\n"));
	}
	if (the_fit == &gfit)
		undo_save_state(&gfit, _("Deconvolution"));
	args.ndata = the_fit->rx * the_fit->ry * the_fit->naxes[2];
	args.fdata = malloc(args.ndata * sizeof(float));
	if (the_fit->type == DATA_FLOAT)
		memcpy(args.fdata, the_fit->fdata, args.ndata * sizeof(float));
	else {
		float invnorm = 1.f / USHRT_MAX_SINGLE;
		for (size_t i = 0 ; i < args.ndata ; i++) {
			args.fdata[i] = (float) the_fit->data[i] * invnorm;
		}
	}

	// Get the kernel
	if (sequence_is_running == 0)
		set_progress_bar_data("Starting kernel estimation...", PROGRESS_PULSATE);
	siril_debug_print("Starting kernel estimation\n");
	get_kernel();
	if (!com.kernel) {
		siril_debug_print("Kernel missing!\n");
		retval = 1;
		goto ENDDECONV;
	}
	DrawPSF();

	if (get_thread_run() || sequence_is_running == 1) {
		if (sequence_is_running == 0)
			set_progress_bar_data("Starting non-blind deconvolution...", 0);
		siril_debug_print("Starting non-blind deconvolution\n");
		// Normalize the kernel
		float ksum = 0.f;
		for (unsigned i = 0 ; i < args.ks * args.ks ; i++)
			ksum += com.kernel[i];
		for (unsigned i = 0 ; i < args.ks * args.ks ; i++)
			com.kernel[i] /= ksum;

		// Non-blind deconvolution stage
		switch (args.nonblindtype) {
			case DECONV_SB:
				split_bregman(args.fdata, args.rx, args.ry, args.nchans, com.kernel, args.ks, args.alpha, args.finaliters, fftw_max_thread);
				break;
			case DECONV_RL:
				msg_earlystop = malloc(100*sizeof(BYTE));
				msg_rl = malloc(50*sizeof(BYTE));
				snprintf(msg_earlystop, 99, _("Richardson-Lucy halted early by the stopping criterion after iteration"));
				snprintf(msg_rl, 49, _("Richardson-Lucy deconvolution..."));
				if (args.rl_method == RL_MULT) {
					if (args.regtype == REG_TV_GRAD)
						args.regtype = REG_TV_MULT;
					else if (args.regtype == REG_FH_GRAD)
						args.regtype = REG_FH_MULT;
					else args.regtype = REG_NONE_MULT;
				}
				richardson_lucy(args.fdata, args.rx,args.ry, args.nchans, com.kernel, args.ks, args.alpha, args.finaliters, args.stopcriterion, fftw_max_thread, args.regtype, args.stepsize, args.stopcriterion_active);
				free(msg_rl);
				msg_rl = NULL;
				free(msg_earlystop);
				msg_earlystop = NULL;
				break;
			case DECONV_WIENER:
				wienerdec(args.fdata, args.rx, args.ry, args.nchans, com.kernel, args.ks, args.alpha, fftw_max_thread);
				break;
		}
	}

	// Update the_fit with the result
	if (get_thread_run()) {
		if (the_fit->type == DATA_FLOAT) {
			memcpy(the_fit->fdata, args.fdata, args.ndata * sizeof(float));
			if (com.pref.force_16bit)
				fit_replace_buffer(the_fit, float_buffer_to_ushort(the_fit->fdata, args.ndata), DATA_USHORT);
		} else {
			for (size_t i = 0 ; i < args.ndata ; i++) {
				the_fit->data[i] = roundf_to_WORD(args.fdata[i] * USHRT_MAX_SINGLE);
			}
		}
	}
ENDDECONV:
	// Do not free the kernel here as it is populated into com.kernel
	if (fftwf_export_wisdom_to_filename(com.pref.fftw_conf.wisdom_file) == 1) {
		if (sequence_is_running == 0)
			siril_log_message(_("Siril FFT wisdom updated successfully...\n"));
	} else {
		if (sequence_is_running == 0)
			siril_log_message(_("Siril FFT wisdom update failed...\n"));
	}

	if (sequence_is_running == 0)
		siril_add_idle(deconvolve_idle, NULL);
	else {
		if (args.fdata) {
			free(args.fdata);
			args.fdata = NULL;
		}
	}
	return GINT_TO_POINTER(retval);
}

void on_bdeconv_symkern_toggled(GtkToggleButton *button, gpointer user_data) {
	args.symkern = gtk_toggle_button_get_active(button);
	start_in_new_thread(estimate_only, NULL);

}

void on_bdeconv_apply_clicked(GtkButton *button, gpointer user_data) {
	sequence_is_running = 0;
	control_window_switch_to_tab(OUTPUT_LOGS);
	GtkToggleButton* seq = GTK_TOGGLE_BUTTON(lookup_widget("bdeconv_seqapply"));
	deconvolution_sequence_data* seqargs = NULL;
	GtkEntry* deconvolutionSeqEntry = GTK_ENTRY(lookup_widget("bdeconv_seq_prefix"));
	set_estimate_params(); // Do this before entering the thread as it contains GTK functions
	set_deconvolve_params();
	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		seqargs = malloc(sizeof(deconvolution_sequence_data));
		seqargs->seqEntry = gtk_entry_get_text(deconvolutionSeqEntry);
		set_cursor_waiting(TRUE);
		if (seqargs->seqEntry && seqargs->seqEntry[0] == '\0')
			seqargs->seqEntry = "dec_";
		seqargs->seq = &com.seq;

		apply_deconvolve_to_sequence(seqargs);
	} else {
		the_fit = &gfit;
		start_in_new_thread(deconvolve, NULL);
	}
}

void on_bdeconv_estimate_clicked(GtkButton *button, gpointer user_data) {
	sequence_is_running = 0;
	control_window_switch_to_tab(OUTPUT_LOGS);
	if(!sequence_is_loaded())
		the_fit = &gfit;
	if(!com.headless)
		set_estimate_params(); // Do this before entering the thread as it contains GTK functions
	start_in_new_thread(estimate_only, NULL);
}

// Actual drawing function
void drawing_the_PSF(GtkWidget *widget, cairo_t *cr) {
	if (!com.kernel || !com.kernelsize) return;
	int width =  gtk_widget_get_allocated_width(widget);
	int height = gtk_widget_get_allocated_height(widget);

	int stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, com.kernelsize);

	float minval = FLT_MAX, maxval = -FLT_MAX;
	for (int i = 0; i < com.kernelsize * com.kernelsize; i++) {
		if (maxval < com.kernel[i]) maxval = com.kernel[i];
		if (minval > com.kernel[i]) minval = com.kernel[i];
	}
	float invrange = (maxval == minval) ? 1.f : 1.f / (maxval - minval);

	guchar *buf = calloc(com.kernelsize * com.kernelsize * 4, sizeof(guchar));
	for (int i = 0; i <com.kernelsize; i++) {
		for (int j = 0; j < com.kernelsize; j++) {
			float val = pow((com.kernel[(com.kernelsize - i - 1) * com.kernelsize + j] - minval) * invrange, 0.5f);
			buf[i * stride + 4 * j + 0] = float_to_uchar_range(val);
			buf[i * stride + 4 * j + 1] = buf[i * stride + 4 * j];
			buf[i * stride + 4 * j + 2] = buf[i * stride + 4 * j];
		}
	}

	if (surface)
		cairo_surface_destroy(surface);

	surface = cairo_image_surface_create_for_data(buf,
				CAIRO_FORMAT_RGB24, com.kernelsize, com.kernelsize, stride);

	cairo_set_source_surface(cr, surface, 0, 0);
	cairo_surface_set_device_scale(surface, (double)com.kernelsize / (double)width, (double)com.kernelsize / (double)height);
	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_FAST);
	cairo_paint(cr);
	free(buf);
}

// PSF drawing callback
gboolean on_PSFkernel_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	drawing_the_PSF(widget, cr);
	return FALSE;
}

// caller
static void DrawPSF() {
	if (!drawingPSF) {
		drawingPSF = lookup_widget("bdeconv_drawingarea");
	}
	if (!com.kernelsize)
		return;
	gtk_widget_queue_draw(drawingPSF);
}

///////// ****** SEQUENCE PROCESSING ****** //////////

static int deconvolution_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
// The deconvolution routines use FFTW which optimises its planning on the basis of the virtual
// memory model and available threads. It is therefore pragmatic to allow FFTW to optimise itself
// and run the sequence sequentially rather than in parallel.
	return 1;
}
int deconvolution_finalize_hook(struct generic_seq_args *seqargs) {
	int retval = 0;
//	deconvolution_sequence_data *seqdata = (deconvolution_sequence_data *) args->user;
	g_assert(seqargs->has_output); // don't call this hook otherwise
	if ((seqargs->force_ser_output || (seqargs->seq->type == SEQ_SER && !seqargs->force_fitseq_output)) && seqargs->new_ser) {
		retval = ser_write_and_close(seqargs->new_ser);
		free(seqargs->new_ser);
	}
	else if ((seqargs->force_fitseq_output || (seqargs->seq->type == SEQ_FITSEQ && !seqargs->force_ser_output)) && seqargs->new_fitseq) {
		retval = fitseq_close_file(seqargs->new_fitseq);
		free(seqargs->new_fitseq);
	}
	args.psftype = args.oldpsftype; // Restore consistency
	sequence_is_running = 0;
	return retval;
}

int deconvolution_image_hook(struct generic_seq_args *seqargs, int o, int i, fits *fit, rectangle *_, int threads) {
	int ret = 0;
//	deconvolution_sequence_data *seqdata = (deconvolution_sequence_data *) args->user;
	the_fit = fit;
	deconvolve(NULL);
	args.oldpsftype = args.psftype; // Need to store the previous psf type so we can restore
	// it later and avoid inconsistency between the GTK widget and the parameter.
	args.psftype = PSF_PREVIOUS; // For all but the first image in the sequence we will reuse the kernel calculated for the first image.

	return ret;
}

void apply_deconvolve_to_sequence(struct deconvolution_sequence_data *seqdata) {
	sequence_is_running = 1;
	struct generic_seq_args *seqargs = create_default_seqargs(seqdata->seq);
	seqargs->seq = seqdata->seq;
	seqargs->filtering_criterion = seq_filter_included;
	seqargs->nb_filtered_images = seqdata->seq->selnum;
	seqargs->compute_mem_limits_hook = deconvolution_compute_mem_limits;
	seqargs->prepare_hook = seq_prepare_hook;
	seqargs->finalize_hook = deconvolution_finalize_hook;
	seqargs->image_hook = deconvolution_image_hook;
	seqargs->description = _("Deconvolution");
	seqargs->has_output = TRUE;
	seqargs->output_type = get_data_type(seqargs->seq->bitpix);
	seqargs->new_seq_prefix = seqdata->seqEntry;
	seqargs->load_new_sequence = TRUE;
	seqargs->user = seqdata;

	start_in_new_thread(generic_sequence_worker, seqargs);
}

gpointer deconvolve_sequence_command(gpointer p, sequence* seqname) {
	int retval = 0;
	if (p != NULL) {
		estk_data *command_data = (estk_data *) p;
		memcpy(&args, command_data, sizeof(estk_data));
		free(command_data);
	}
	deconvolution_sequence_data* seqargs = malloc(sizeof(deconvolution_sequence_data));
	seqargs->seqEntry = "dec_";
	seqargs->seq = seqname;
	apply_deconvolve_to_sequence(seqargs);
	return GINT_TO_POINTER(retval);

}
