/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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
#include <math.h>
#include <locale.h>
#include <gdk/gdk.h>
#include "core/siril.h"
#include "core/siril_date.h"
#include "core/command.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/registration_preview.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "gui/utils.h"
#include "gui/siril_preview.h"
#include "gui/progress_and_log.h"
#include "gui/newdeconv.h"
#include "filters/deconvolution/deconvolution.h"
#include "algos/PSF.h"
#include "io/sequence.h"

extern gboolean aperture_warning_given;
extern gboolean bad_load;
extern orientation_t imageorientation;
extern gboolean next_psf_is_previous;
static gboolean first_time = TRUE;
extern fits* the_fit;
extern estk_data args;
static cairo_surface_t *surface = NULL;
extern int sequence_is_running;

// Static widget pointers
static GtkButton *bdeconv_close = NULL, *bdeconv_roi_preview = NULL, *bdeconv_apply = NULL, *bdeconv_estimate = NULL, *bdeconv_savekernel = NULL, *bdeconv_advice_button = NULL;
static GtkComboBox *bdeconv_blindtype = NULL, *bdeconv_profile = NULL, *bdeconv_nonblindtype = NULL, *bdeconv_rl_regularization = NULL, *bdeconv_rl_method = NULL;
static GtkDialog *bdeconv_dialog = NULL;
static GtkDrawingArea *bdeconv_drawingarea = NULL;
static GtkEntry *bdeconv_seq_prefix = NULL;
static GtkExpander *bdeconv_expander = NULL;
static GtkFileChooser *bdeconv_filechooser = NULL;
static GtkFrame *bdeconv_psfcontrols = NULL, *bdeconv_blindcontrols = NULL, *bdeconv_l0controls = NULL, *bdeconv_gfcontrols = NULL, *bdeconv_starpsf_details = NULL;
static GtkGrid *bdeconv_manual_stars = NULL, *bdeconv_manual_airy = NULL;
static GtkLabel *bdeconv_starprofile_text = NULL, *bdeconv_starfwhm_text = NULL, *bdeconv_starratio_text = NULL, *bdeconv_starangle_text = NULL, *bdeconv_starbeta_text = NULL, *bdeconv_iterlabel = NULL, *regul_label = NULL, *bdeconv_steplabel = NULL, *algo_method_label = NULL;
static GtkToggleButton *bdeconv_psfblind = NULL, *bdeconv_psfstars = NULL, *bdeconv_psfmanual = NULL, *bdeconv_psfprevious = NULL, *bdeconv_psfselection = NULL;
static GtkSpinButton *bdeconv_ks = NULL, *bdeconv_psfratio = NULL, *bdeconv_psfwhm = NULL, *bdeconv_psfangle = NULL, *bdeconv_psfbeta = NULL, *airy_pixelsize = NULL, *airy_wl = NULL, *airy_fl = NULL, *airy_diameter = NULL, *airy_obstruction = NULL, *bdeconv_gflambda = NULL, *bdeconv_gamma = NULL, *bdeconv_iters = NULL, *bdeconv_lambdaratio = NULL, *bdeconv_lambdamin = NULL, *bdeconv_scalefactor = NULL, *bdeconv_upsampleblur = NULL, *bdeconv_downsampleblur = NULL, *bdeconv_kthresh = NULL, *bdeconv_kl1 = NULL, *bdeconv_ncomp = NULL, *bdeconv_ntries = NULL, *bdeconv_nouter = NULL, *bdeconv_ninner = NULL, *bdeconv_finaliters = NULL, *bdeconv_alpha = NULL, *bdeconv_stopcriterion = NULL, *bdeconv_stepsize = NULL;
static GtkToggleButton *bdeconv_multiscale = NULL, *bdeconv_betterkernel = NULL, *bdeconv_symkern = NULL, *bdeconv_stopping_toggle = NULL, *bdeconv_seqapply = NULL;

void bdeconv_dialog_init_statics() {
	if (bdeconv_close == NULL) {
		// GtkButton
		bdeconv_close = GTK_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_close"));
		bdeconv_roi_preview = GTK_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_roi_preview"));
		bdeconv_apply = GTK_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_apply"));
		bdeconv_estimate = GTK_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_estimate"));
		bdeconv_savekernel = GTK_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_savekernel"));
		bdeconv_advice_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_advice_button"));
		// GtkComboBoxText
		bdeconv_blindtype = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "bdeconv_blindtype"));
		bdeconv_profile = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "bdeconv_profile"));
		bdeconv_nonblindtype = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "bdeconv_nonblindtype"));
		bdeconv_rl_regularization = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "bdeconv_rl_regularization"));
		bdeconv_rl_method = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "bdeconv_rl_method"));
		// GtkDialog
		bdeconv_dialog = GTK_DIALOG(gtk_builder_get_object(gui.builder, "bdeconv_dialog"));
		// GtkDrawingArea
		bdeconv_drawingarea = GTK_DRAWING_AREA(gtk_builder_get_object(gui.builder, "bdeconv_drawingarea"));
		// GtkEntry
		bdeconv_seq_prefix = GTK_ENTRY(gtk_builder_get_object(gui.builder, "bdeconv_seq_prefix"));
		// GtkExpander
		bdeconv_expander = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "bdeconv_expander"));
		// GtkFileChooserButton (cast as a GtkFileChooser*)
		bdeconv_filechooser = GTK_FILE_CHOOSER(gtk_builder_get_object(gui.builder, "bdeconv_filechooser"));
		// GtkFrame
		bdeconv_psfcontrols = GTK_FRAME(gtk_builder_get_object(gui.builder, "bdeconv_psfcontrols"));
		bdeconv_blindcontrols = GTK_FRAME(gtk_builder_get_object(gui.builder, "bdeconv_blindcontrols"));
		bdeconv_l0controls = GTK_FRAME(gtk_builder_get_object(gui.builder, "bdeconv_l0controls"));
		bdeconv_gfcontrols = GTK_FRAME(gtk_builder_get_object(gui.builder, "bdeconv_gfcontrols"));
		bdeconv_starpsf_details = GTK_FRAME(gtk_builder_get_object(gui.builder, "bdeconv_starpsf_details"));
		// GtkGrid
		bdeconv_manual_stars = GTK_GRID(gtk_builder_get_object(gui.builder, "bdeconv_manual_stars"));
		bdeconv_manual_airy = GTK_GRID(gtk_builder_get_object(gui.builder, "bdeconv_manual_airy"));
		// GtkLabel
		bdeconv_starprofile_text = GTK_LABEL(gtk_builder_get_object(gui.builder, "bdeconv_starprofile_text"));
		bdeconv_starfwhm_text = GTK_LABEL(gtk_builder_get_object(gui.builder, "bdeconv_starfwhm_text"));
		bdeconv_starratio_text = GTK_LABEL(gtk_builder_get_object(gui.builder, "bdeconv_starratio_text"));
		bdeconv_starangle_text = GTK_LABEL(gtk_builder_get_object(gui.builder, "bdeconv_starangle_text"));
		bdeconv_starbeta_text = GTK_LABEL(gtk_builder_get_object(gui.builder, "bdeconv_starbeta_text"));
		bdeconv_iterlabel = GTK_LABEL(gtk_builder_get_object(gui.builder, "bdeconv_iterlabel"));
		regul_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "regul_label"));
		bdeconv_steplabel = GTK_LABEL(gtk_builder_get_object(gui.builder, "bdeconv_steplabel"));
		algo_method_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "algo_method_label"));
		// GtkRadioButton (treat these as toggle buttons)
		bdeconv_psfblind = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfblind"));
		bdeconv_psfstars = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfstars"));
		bdeconv_psfmanual = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfmanual"));
		bdeconv_psfprevious = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfprevious"));
		bdeconv_psfselection = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfselection"));
		// GtkSpinButton
		bdeconv_ks = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_ks"));
		bdeconv_psfratio = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfratio"));
		bdeconv_psfwhm = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfwhm"));
		bdeconv_psfangle = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfangle"));
		bdeconv_psfbeta = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_psfbeta"));
		airy_pixelsize = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "airy_pixelsize"));
		airy_wl = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "airy_wl"));
		airy_fl = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "airy_fl"));
		airy_diameter = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "airy_diameter"));
		airy_obstruction = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "airy_obstruction"));
		bdeconv_gflambda = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_gflambda"));
		bdeconv_gamma = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_gamma"));
		bdeconv_iters = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_iters"));
		bdeconv_lambdaratio = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_lambdaratio"));
		bdeconv_lambdamin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_lambdamin"));
		bdeconv_scalefactor = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_scalefactor"));
		bdeconv_upsampleblur = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_upsampleblur"));
		bdeconv_downsampleblur = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_downsampleblur"));
		bdeconv_kthresh = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_kthresh"));
		bdeconv_kl1 = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_kl1"));
		bdeconv_ncomp = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_ncomp"));
		bdeconv_ntries = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_ntries"));
		bdeconv_nouter = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_nouter"));
		bdeconv_ninner = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_ninner"));
		bdeconv_finaliters = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_finaliters"));
		bdeconv_alpha = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_alpha"));
		bdeconv_stopcriterion = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_stopcriterion"));
		bdeconv_stepsize = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_stepsize"));
		// GtkToggleButton
		bdeconv_multiscale = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_multiscale"));
		bdeconv_betterkernel = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_betterkernel"));
		bdeconv_symkern = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_symkern"));
		bdeconv_stopping_toggle = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_stopping_toggle"));
		bdeconv_seqapply = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "bdeconv_seqapply"));
	}
}

//set below flag to zero to avoid kernel printout to stdout
#define DEBUG_PSF 0

gboolean reset_conv_controls(gpointer user_data) {
	gtk_combo_box_set_active(bdeconv_profile, args.profile);
	gtk_combo_box_set_active(bdeconv_blindtype, args.blindtype);
	gtk_combo_box_set_active(bdeconv_nonblindtype, args.nonblindtype);
	gtk_combo_box_set_active(bdeconv_rl_regularization, args.regtype);
	gtk_combo_box_set_active(bdeconv_rl_method, args.rl_method);
	gtk_toggle_button_set_active(bdeconv_multiscale, args.multiscale);
	gtk_toggle_button_set_active(bdeconv_psfblind, TRUE);
	gtk_spin_button_set_value(bdeconv_lambdaratio, args.lambda_ratio);
	gtk_spin_button_set_value(bdeconv_lambdamin, args.lambda_min);
	gtk_spin_button_set_value(bdeconv_gamma, args.gamma);
	gtk_spin_button_set_value(bdeconv_iters, args.iterations);
	gtk_spin_button_set_value(bdeconv_scalefactor, args.scalefactor);
	gtk_spin_button_set_value(bdeconv_kthresh, args.kernel_threshold_max);
	gtk_toggle_button_set_active(bdeconv_betterkernel, args.better_kernel);
	gtk_spin_button_set_value(bdeconv_upsampleblur, args.upscaleblur);
	gtk_spin_button_set_value(bdeconv_downsampleblur, args.downscaleblur);
	gtk_spin_button_set_value(bdeconv_kl1, args.k_l1);
	gtk_spin_button_set_value(bdeconv_alpha, 1.f / args.alpha);
	gtk_spin_button_set_value(bdeconv_psfwhm, args.psf_fwhm);
	gtk_spin_button_set_value(bdeconv_psfbeta, args.psf_beta);
	gtk_spin_button_set_value(bdeconv_psfratio, args.psf_ratio);
	gtk_spin_button_set_value(airy_diameter, args.airy_diameter);
	gtk_spin_button_set_value(airy_fl, args.airy_fl);
	gtk_spin_button_set_value(airy_wl, args.airy_wl);
	gtk_spin_button_set_value(airy_pixelsize, args.airy_pixelsize);
	gtk_spin_button_set_value(airy_obstruction, args.airy_obstruction);
	gtk_spin_button_set_value(bdeconv_ninner, args.ninner);
	gtk_spin_button_set_value(bdeconv_ntries, args.ntries);
	gtk_spin_button_set_value(bdeconv_nouter, args.nouter);
	gtk_spin_button_set_value(bdeconv_finaliters, args.finaliters);
	gtk_spin_button_set_value(bdeconv_stopcriterion, args.stopcriterion);
	gtk_spin_button_set_value(bdeconv_stepsize, args.stepsize);
	gtk_spin_button_set_value(bdeconv_ncomp, args.compensationfactor);
	return FALSE;
}

void on_bdeconv_psfblind_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_BLIND;
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_blindcontrols), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_psfcontrols), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_starpsf_details), FALSE);
	if (args.blindtype == BLIND_SI) {
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_gfcontrols), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_l0controls), FALSE);
	} else if (args.blindtype == BLIND_L0) {
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_gfcontrols), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_l0controls), TRUE);
	}
}

void on_bdeconv_ks_value_changed(GtkSpinButton *button, gpointer user_data) {
	args.ks = gtk_spin_button_get_value(button);
	// Prevent setting ks to an even value.
	if(!(args.ks % 2)) {
		args.ks++;
		gtk_spin_button_set_value(button, args.ks);
	}
	reset_conv_kernel();
	DrawPSF(NULL);
}

void on_bdeconv_advice_button_clicked(GtkButton *button, gpointer user_data) {
	// Copypasta from documentation.c but with a specific URL to point to the deconvolution tips page
	#define GET_DOCUMENTATION_URL "https://siril.readthedocs.io"
	#define DECONVOLUTION_TIPS_URL "processing/deconvolution.html#deconvolution-usage-tips"

	gboolean ret;
	const char *locale;
	const char *supported_languages[] = { "de", "fr", NULL }; // en is NULL: default language
	gchar *lang = NULL;
	int i = 0;

	if (!com.pref.lang || !g_strcmp0(com.pref.lang, "")) {
		locale = setlocale(LC_MESSAGES, NULL);
	} else {
		locale = com.pref.lang;
	}

	if (locale) {
		while (supported_languages[i]) {
			if (!strncmp(locale, supported_languages[i], 2)) {
				lang = g_strndup(locale, 2);
				break;
			}
			i++;
		}
	}
	if (!lang) {
		lang = g_strdup_printf("en"); // Last gasp fallback in case there is an error with the locale
	}
	/* Use the tag when documentation will be tagged */
	const gchar *version = NULL;
#ifdef SIRIL_UNSTABLE
	version = "latest";
#else
	version = "stable";
#endif

	gchar *url = g_strdup_printf("%s/%s/%s/%s", GET_DOCUMENTATION_URL, lang, version, DECONVOLUTION_TIPS_URL);
	siril_log_message(_("Deconvolution usage hints and tips URL: %s\n"), url);
#if GTK_CHECK_VERSION(3, 22, 0)
	GtkWidget* win = lookup_widget("control_window");
	ret = gtk_show_uri_on_window(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)), url,
			gtk_get_current_event_time(), NULL);
#else
	ret = gtk_show_uri(gdk_screen_get_default(), url,
			gtk_get_current_event_time(), NULL);
#endif
	if (!ret) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not show link"),
				_("Please go to <a href=\""GET_DOCUMENTATION_URL"\">"GET_DOCUMENTATION_URL"</a> "
								"by copying the link."));
	}
	g_free(url);
	g_free(lang);
}

void on_bdeconv_blindtype_changed(GtkComboBox *combo, gpointer user_data) {
	args.blindtype = gtk_combo_box_get_active(combo);
	if (args.psftype == PSF_BLIND) {
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_psfcontrols), FALSE);
		switch (args.blindtype) {
			case BLIND_SI:
				gtk_widget_set_visible(GTK_WIDGET(bdeconv_l0controls), FALSE);
				gtk_widget_set_visible(GTK_WIDGET(bdeconv_gfcontrols), TRUE);
				args.intermediatedeconvolutionweight = 15.f;
				args.finaldeconvolutionweight = 15.f;
				args.lambda = 1.f / 15.f;
				gtk_spin_button_set_value(bdeconv_gflambda, 15.);
				break;
			case BLIND_L0:
				gtk_widget_set_visible(GTK_WIDGET(bdeconv_l0controls), TRUE);
				gtk_widget_set_visible(GTK_WIDGET(bdeconv_gfcontrols), FALSE);
				args.lambda = 1.f / 3000.f;
				gtk_spin_button_set_value(bdeconv_gflambda, 3000.);
				break;
		}
	}
}

void on_bdeconv_nonblindtype_changed(GtkComboBox *combo, gpointer user_data) {
	args.nonblindtype = gtk_combo_box_get_active(combo);
	switch (args.nonblindtype) {
		case DECONV_SB:
			args.finaliters = 1; // Default niters for Split Bregman
			gtk_spin_button_set_value(bdeconv_finaliters, (double) args.finaliters);
			args.alpha = 3000.f;
			gtk_spin_button_set_value(bdeconv_alpha, (double) args.alpha);
			gtk_widget_set_visible(GTK_WIDGET(regul_label), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_rl_regularization), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(algo_method_label), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_rl_method), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stopping_toggle), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stopcriterion), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stepsize), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_steplabel), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_finaliters), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_iterlabel), TRUE);
			break;
		case DECONV_RL:
			args.finaliters = 10; // Default niters for RL. Reduce this if using the multiplicative
			// algorithm in order to avoid burning holes round your stars!
			gtk_spin_button_set_value(bdeconv_finaliters, (double) args.finaliters);
			args.alpha = 3000.f;
			gtk_spin_button_set_value(bdeconv_alpha, (double) args.alpha);
			gtk_widget_set_visible(GTK_WIDGET(regul_label), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_rl_regularization), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(algo_method_label), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_rl_method), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stopping_toggle), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stopcriterion), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stepsize), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_steplabel), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_finaliters), TRUE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_iterlabel), TRUE);
			break;
		case DECONV_WIENER:
			args.alpha = 500.f;
			gtk_spin_button_set_value(bdeconv_alpha, (double) args.alpha);
			gtk_widget_set_visible(GTK_WIDGET(regul_label), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_rl_regularization), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(algo_method_label), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_rl_method), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stopping_toggle), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stopcriterion), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_stepsize), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_steplabel), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_finaliters), FALSE);
			gtk_widget_set_visible(GTK_WIDGET(bdeconv_iterlabel), FALSE);
			break;
	}
	gtk_spin_button_set_value(bdeconv_finaliters, args.finaliters);
}

void on_bdeconv_expander_activate(GtkExpander *expander, gpointer user_data) {
	gtk_window_resize(GTK_WINDOW(bdeconv_dialog), 1, 1);
	// This callback is just to prevent excessive blank space after shrinking the expander
}

void on_airy_diameter_value_changed(GtkWidget *button, gpointer user_data) {
	GtkCssProvider *css = gtk_css_provider_new();
	gtk_css_provider_load_from_data(css, "* { background-image:none; color:@theme_color;}",-1,NULL);
	GtkStyleContext * context = gtk_widget_get_style_context(button);

	gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
	g_object_unref(css);
}

void on_airy_fl_value_changed(GtkWidget *button, gpointer user_data) {
	GtkCssProvider *css = gtk_css_provider_new();
	gtk_css_provider_load_from_data(css, "* { background-image:none; color:@theme_color;}",-1,NULL);
	GtkStyleContext * context = gtk_widget_get_style_context(button);

	gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
	g_object_unref(css);
}

void on_airy_pixelsize_value_changed(GtkWidget *button, gpointer user_data) {
	GtkCssProvider *css = gtk_css_provider_new();
	gtk_css_provider_load_from_data(css, "* { background-image:none; color:@theme_color;}",-1,NULL);
	GtkStyleContext * context = gtk_widget_get_style_context(button);

	gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
	g_object_unref(css);
}

static void initialize_airy_parameters() {
	// Get initial stab at parameters for Airy function from FITS header. Not essential they be correct as they are user-editable in the UI.
	args.airy_fl = the_fit->keywords.focal_length;
	gtk_spin_button_set_value(airy_fl, args.airy_fl);
	if (the_fit->keywords.focal_length < 10.) {
		GtkWidget *button = GTK_WIDGET(airy_fl);
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:salmon;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(button);
		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
	}
	args.airy_diameter = the_fit->keywords.aperture;
	gtk_spin_button_set_value(airy_diameter, args.airy_diameter);
	if (the_fit->keywords.aperture < 10.) {
		GtkWidget *button = GTK_WIDGET(airy_diameter);
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:salmon;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(button);
		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
	}
	args.airy_pixelsize = the_fit->keywords.pixel_size_x;
	gtk_spin_button_set_value(airy_pixelsize, args.airy_pixelsize);
	if (the_fit->keywords.pixel_size_x <=1.) {
		GtkWidget *button = GTK_WIDGET(airy_pixelsize);
		GtkCssProvider *css = gtk_css_provider_new();
		gtk_css_provider_load_from_data(css, "* { background-image:none; color:salmon;}",-1,NULL);
		GtkStyleContext * context = gtk_widget_get_style_context(button);
		gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(css),GTK_STYLE_PROVIDER_PRIORITY_USER);
		g_object_unref(css);
	}
	args.airy_wl = 525.f; // Approximation of the average wavelength of light accepted by an OSC camera covering 400-656nm. User editable to give precise values e.g. for narrowband.
	gtk_spin_button_set_value(airy_wl, args.airy_wl);
}

void on_bdeconv_psfmanual_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_MANUAL;
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_psfcontrols), TRUE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_l0controls), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_gfcontrols), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_blindcontrols), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_starpsf_details), FALSE);
}

void on_bdeconv_profile_changed(GtkComboBox *combo, gpointer user_data) {
	int profile = gtk_combo_box_get_active(combo);
	if (profile != PROFILE_AIRY) {
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_manual_stars), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_manual_airy), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(bdeconv_psfangle), (profile == PROFILE_MOFFAT || profile == PROFILE_GAUSSIAN));
		gtk_widget_set_sensitive(GTK_WIDGET(bdeconv_psfratio), (profile == PROFILE_MOFFAT || profile == PROFILE_GAUSSIAN));
		gtk_widget_set_sensitive(GTK_WIDGET(bdeconv_psfbeta), profile == PROFILE_MOFFAT);
	} else {
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_manual_stars), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(bdeconv_manual_airy), TRUE);
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
	gtk_widget_set_sensitive(GTK_WIDGET(bdeconv_stepsize), gtk_combo_box_get_active(combo) == 1);
}

void on_bdeconv_psfprevious_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_PREVIOUS;
}

void on_bdeconv_reset_clicked(GtkButton *button, gpointer user_data) {
	reset_conv_controls_and_args();
}

void calculate_parameters() {
	if (com.stars && com.stars[0]) {
		int i = 0;
		double FWHMx = 0.0, FWHMy = 0.0, beta = 0.0, angle = 0.0;
		gboolean unit_is_arcsec = FALSE;
		int n = 0, layer = 0;
		starprofile profiletype = PSF_GAUSSIAN;
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
				siril_log_color_message(_("Error: all stars' FWHM must have the same units.\n"), "red");
				return;
			}
			else if (layer != com.stars[i]->layer ) {
				siril_log_color_message(_("Error: stars' properties must all be computed on the same layer"), "red");
				return;
			}
			else if (profiletype != com.stars[i]->profile) {
				siril_log_color_message(_("Error: stars must all be modeled with the same profile type"), "red");
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
		if (i <= 0 || n <= 0)
			return;
		/* compute average */
		args.psf_fwhm = (float) FWHMx / (float) n;
		FWHMy = (float) FWHMy / (float) n;
		args.psf_beta = (float) beta / (float) n;
		args.psf_ratio = args.symkern ? 1.f : args.psf_fwhm / FWHMy;
		if (unit_is_arcsec) {
			double bin_X = com.pref.binning_update ? (double) the_fit->keywords.binning_x : 1.0;
			double conversionfactor = (((3600.0 * 180.0) / G_PI) / 1.0E3 * (double)the_fit->keywords.pixel_size_x / the_fit->keywords.focal_length) * bin_X;
			args.psf_fwhm /= (float) conversionfactor;
		}
		args.psf_angle = (float) angle / (float) n;
		if (!com.headless && !com.script) {
			if (com.stars[0]->profile == PSF_GAUSSIAN)
				gtk_label_set_text(bdeconv_starprofile_text, _("Gaussian"));
			else
				gtk_label_set_text(bdeconv_starprofile_text, _("Moffat"));
			char temp[10];
			sprintf(temp, "%.2f", args.psf_fwhm);
			gtk_label_set_text(bdeconv_starfwhm_text, temp);
			sprintf(temp, "%.2f", args.psf_ratio);
			gtk_label_set_text(bdeconv_starratio_text, temp);
			sprintf(temp, "%.2f", args.psf_angle);
			gtk_label_set_text(bdeconv_starangle_text, temp);
			sprintf(temp, "%.2f", args.psf_beta);
			gtk_label_set_text(bdeconv_starbeta_text, temp);
			return;
		}
	} else if (!com.headless && !com.script) {
		gtk_label_set_text(bdeconv_starprofile_text, _("No stars selected"));
	}
	return;
}

void on_bdeconv_psfstars_toggled(GtkToggleButton *button, gpointer user_data) {
	args.psftype = PSF_STARS;
	calculate_parameters();
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_psfcontrols), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_l0controls), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_gfcontrols), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_blindcontrols), FALSE);
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_starpsf_details), TRUE);
}

void deconv_roi_callback() {
	gui.roi.operation_supports_roi = TRUE;
	gtk_widget_set_visible(GTK_WIDGET(bdeconv_roi_preview), gui.roi.active);
	the_fit = gui.roi.active ? &gui.roi.fit : &gfit;
	copy_backup_to_gfit();
	notify_gfit_modified();
}

void close_deconv() {
	roi_supported(FALSE);
	siril_preview_hide();
	remove_roi_callback(deconv_roi_callback);
	siril_close_dialog("bdeconv_dialog");
}

void on_bdeconv_close_clicked(GtkButton *button, gpointer user_data) {
	close_deconv();
}

void on_bdeconv_dialog_show(GtkWidget *widget, gpointer user_data) {
	bdeconv_dialog_init_statics();
	the_fit = gui.roi.active ? &gui.roi.fit : &gfit;
	roi_supported(TRUE);
	deconv_roi_callback();
	add_roi_callback(deconv_roi_callback);
	copy_gfit_to_backup();
	if (first_time) {
		reset_conv_controls_and_args();
		first_time = FALSE;
	}
	if (com.kernel && com.kernelsize > 0) {
		args.psftype = PSF_PREVIOUS;
			gtk_toggle_button_set_active(bdeconv_psfprevious, TRUE);
	}
	calculate_parameters();
	initialize_airy_parameters();
	control_window_switch_to_tab(OUTPUT_LOGS);
}

void on_bdeconv_savekernel_clicked(GtkButton *button, gpointer user_data) {
	// Only allocate as much space for filenames as required - we determine the max pathlength
	long pathmax = get_pathmax();
	gchar *filename = NULL;
	gchar *temp1 = NULL;
	gchar *temp2 = NULL;
	gchar *temp3 = NULL;
	gchar *temp4 = NULL;
	gchar *temp5 = NULL;
	gchar *temp6 = NULL;
	gchar kernelsuffix[10] = "_PSF";
	// Set up paths and filenames
	if (single_image_is_loaded())
		temp1 = g_path_get_basename(com.uniq->filename);
	else if (sequence_is_loaded())
		temp1 = g_strdup(com.seq.seqname);
	else
		temp1 = g_strdup_printf("deconvolution");
	temp2 = build_timestamp_filename();
	temp3 = g_strdup_printf("%s_%s", temp2, temp1);
	temp4 = g_build_filename(com.wd, temp3, NULL);
	temp5 = remove_ext_from_filename(temp4);
	temp6 = g_strdup_printf("%s%s", temp5, kernelsuffix);
#ifdef HAVE_LIBTIFF
	filename = g_strdup_printf("%s.tif", temp6);
#else
	filename = g_strdup_printf("%s.fit", temp6);
#endif
	if (strlen(filename) > pathmax) {
		siril_log_color_message(_("Error: file path too long!\n"), "red");
	} else {
		save_kernel(filename);
	}
	g_free(temp1);
	g_free(temp2);
	g_free(temp3);
	g_free(temp4);
	g_free(temp5);
	g_free(temp6);
	g_free(filename);
	return;
}

gboolean set_kernel_size_in_gui(gpointer user_data) {
	gtk_spin_button_set_value(bdeconv_ks, com.kernelsize);
	return FALSE;
}

void on_bdeconv_filechooser_file_set(GtkFileChooser *filechooser, gpointer user_data) {
	gchar* filename = siril_file_chooser_get_filename(filechooser);
	if (filename == NULL) {
		siril_log_color_message(_("No PSF file selected.\n"), "red");
		gtk_file_chooser_unselect_all(filechooser);
	} else {
		reset_conv_kernel();
		load_kernel(filename);
	}
	if (filename)
		g_free(filename);
	if (bad_load) {
		gtk_file_chooser_unselect_all(filechooser);
		bad_load = FALSE;
	} else {
		args.psftype = PSF_PREVIOUS; // Set to use previous kernel
		gtk_toggle_button_set_active(bdeconv_psfprevious, TRUE);
	}
	return;
}

void set_estimate_params() {
	if (com.headless)
		return;
	args.ks = gtk_spin_button_get_value(bdeconv_ks);
	args.gamma = gtk_spin_button_get_value(bdeconv_gamma);
	args.iterations = gtk_spin_button_get_value(bdeconv_iters);
	args.lambda_ratio = gtk_spin_button_get_value(bdeconv_lambdaratio);
	args.lambda_min = gtk_spin_button_get_value(bdeconv_lambdamin);
	args.scalefactor = gtk_spin_button_get_value(bdeconv_scalefactor);
	args.upscaleblur = gtk_spin_button_get_value(bdeconv_upsampleblur);
	args.downscaleblur = gtk_spin_button_get_value(bdeconv_downsampleblur);
	args.kernel_threshold_max = gtk_spin_button_get_value(bdeconv_kthresh);
	args.k_l1 = gtk_spin_button_get_value(bdeconv_kl1);
	args.ninner = gtk_spin_button_get_value(bdeconv_ninner);
	args.ntries = gtk_spin_button_get_value(bdeconv_ntries);
	args.nouter = gtk_spin_button_get_value(bdeconv_ntries);
	args.compensationfactor = gtk_spin_button_get_value(bdeconv_ntries);
	args.intermediatedeconvolutionweight = gtk_spin_button_get_value(bdeconv_gflambda);
	args.finaldeconvolutionweight = args.intermediatedeconvolutionweight;
	args.lambda = 1.f / args.intermediatedeconvolutionweight;
	args.profile = gtk_combo_box_get_active(bdeconv_profile);
	args.psf_fwhm = gtk_spin_button_get_value(bdeconv_psfwhm);
	args.psf_ratio = gtk_spin_button_get_value(bdeconv_psfratio);
	args.psf_beta = gtk_spin_button_get_value(bdeconv_psfbeta);
	args.psf_angle = gtk_spin_button_get_value(bdeconv_psfangle);
	args.airy_diameter = gtk_spin_button_get_value(airy_diameter);
	args.airy_fl = gtk_spin_button_get_value(airy_fl);
	args.airy_wl = gtk_spin_button_get_value(airy_wl);
	args.airy_pixelsize = gtk_spin_button_get_value(airy_pixelsize);
	args.airy_obstruction = gtk_spin_button_get_value(airy_obstruction) / 100.f;
	args.better_kernel = gtk_toggle_button_get_active(bdeconv_betterkernel);
	args.multiscale = gtk_toggle_button_get_active(bdeconv_multiscale);
}

gboolean estimate_idle(gpointer arg) {
	if (args.psftype == PSF_BLIND || args.psftype == PSF_STARS) {
		args.psftype = PSF_PREVIOUS; // If a blind estimate is done, switch to previous kernel in order to avoid recalculating it when Apply is clicked
		gtk_toggle_button_set_active(bdeconv_psfprevious, TRUE);
	}
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	set_cursor_waiting(FALSE);
	siril_debug_print("Estimate idle stopping processing thread\n");
	stop_processing_thread();
	siril_free(args.fdata);
	args.fdata = NULL;
	DrawPSF(NULL);
	return FALSE;
}

void set_deconvolve_params() {
	if (com.headless)
		return;
	args.finaliters = gtk_spin_button_get_value(bdeconv_finaliters);
	args.alpha = 1.f / gtk_spin_button_get_value(bdeconv_alpha);
	args.stopcriterion = gtk_spin_button_get_value(bdeconv_stopcriterion);
	args.stopcriterion_active = gtk_toggle_button_get_active(bdeconv_stopping_toggle) ? 1 : 0;
	args.rl_method = gtk_combo_box_get_active(bdeconv_rl_method);
	args.stepsize = gtk_spin_button_get_value(bdeconv_stepsize);
	args.regtype = gtk_combo_box_get_active(bdeconv_rl_regularization);
}

gboolean deconvolve_idle(gpointer arg) {
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	siril_free(args.fdata);
	args.fdata = NULL;
	if (!args.previewing) {
		copy_gfit_to_backup();
		populate_roi();
	}
	notify_gfit_modified(); // Also stops the thread and updates the cursor
	if (next_psf_is_previous && !com.headless && !com.script) {
		args.psftype = PSF_PREVIOUS;
		gtk_toggle_button_set_active(bdeconv_psfprevious, TRUE);
	}
	siril_debug_print("Deconvolve idle stopping processing thread\n");
	return FALSE;
}

void on_bdeconv_symkern_toggled(GtkToggleButton *button, gpointer user_data) {
	args.symkern = gtk_toggle_button_get_active(button);
	start_in_new_thread(estimate_only, NULL);
}

void on_bdeconv_roi_preview_clicked(GtkButton *button, gpointer user_data) {
	sequence_is_running = 0;
	control_window_switch_to_tab(OUTPUT_LOGS);
	set_estimate_params(); // Do this before entering the thread as it contains GTK functions
	set_deconvolve_params();
	if (!check_ok_if_cfa())
		return;
	set_cursor_waiting(TRUE);
	if (gtk_toggle_button_get_active(bdeconv_seqapply) && sequence_is_loaded()) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Sequence selected"), _("Preview cannot be used with \"Apply to Sequence\" selected"));
	} else {
		copy_backup_to_gfit();
		args.previewing = TRUE;
		the_fit = (!com.headless && gui.roi.active) ? &gui.roi.fit : &gfit;
		start_in_new_thread(deconvolve, NULL);
	}
}

void on_bdeconv_apply_clicked(GtkButton *button, gpointer user_data) {
	sequence_is_running = 0;
	control_window_switch_to_tab(OUTPUT_LOGS);
	deconvolution_sequence_data* seqargs = NULL;
	set_estimate_params(); // Do this before entering the thread as it contains GTK functions
	set_deconvolve_params();
	args.previewing = FALSE;
	if (!check_ok_if_cfa())
		return;
	set_cursor_waiting(TRUE);
	if (gtk_toggle_button_get_active(bdeconv_seqapply) && sequence_is_loaded()) {
		seqargs = siril_calloc(1, sizeof(deconvolution_sequence_data));
		seqargs->seq = &com.seq;
		seqargs->from_command = FALSE;
		seqargs->seqEntry = strdup(gtk_entry_get_text(bdeconv_seq_prefix));
		if (seqargs->seqEntry && seqargs->seqEntry[0] == '\0') {
			siril_free(seqargs->seqEntry);
			seqargs->seqEntry = strdup("dec_");
		}
		apply_deconvolve_to_sequence(seqargs);
	} else {
		copy_backup_to_gfit();
		the_fit = &gfit;
		start_in_new_thread(deconvolve, NULL);
	}
}

void on_bdeconv_estimate_clicked(GtkButton *button, gpointer user_data) {
	sequence_is_running = 0;
	gboolean abort = FALSE;
	control_window_switch_to_tab(OUTPUT_LOGS);
	gtk_file_chooser_unselect_all(bdeconv_filechooser);
	if(!sequence_is_loaded())
		the_fit = &gfit; // The blind estimate is still always done on the whole image.
		// TODO: consider if this should be done on the ROI if active...
	if(!com.headless)
		set_estimate_params(); // Do this before entering the thread as it contains GTK functions
	if (args.psftype == PSF_STARS || args.psftype == PSF_BLIND)
		abort = !check_ok_if_cfa();
	if (!abort) {
		set_cursor_waiting(TRUE);
		start_in_new_thread(estimate_only, NULL);
	}
}

// Actual drawing function
void drawing_the_PSF(GtkWidget *widget, cairo_t *cr) {
	static GMutex psf_preview_mutex;
	if (!com.kernel || !com.kernelsize) return;
	g_mutex_lock(&psf_preview_mutex);
	int width =  gtk_widget_get_allocated_width(widget);
	int height = gtk_widget_get_allocated_height(widget);

	int stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, com.kernelsize);

	float minval = FLT_MAX, maxval = -FLT_MAX;
	for (int i = 0; i < com.kernelsize * com.kernelsize * com.kernelchannels; i++) {
		if (maxval < com.kernel[i]) maxval = com.kernel[i];
		if (minval > com.kernel[i]) minval = com.kernel[i];
	}
	float invrange = (maxval == minval) ? 1.f : 1.f / (maxval - minval);
	int kpixels = com.kernelsize * com.kernelsize;

	guchar *buf = siril_calloc(com.kernelsize * com.kernelsize * 4, sizeof(guchar));
	for (int i = 0; i < com.kernelsize; i++) {
		for (int j = 0; j < com.kernelsize; j++) {
			float val[3] = { 0.f };
			if (com.kernelchannels == 1) {
				val[0] = pow((com.kernel[(com.kernelsize - i - 1) * com.kernelsize + j] - minval) * invrange, 0.5f);
				val[1] = val[0];
				val[2] = val[0];
			} else if (com.kernelchannels == 3) {
				val[0] = pow((com.kernel[2 * kpixels + (com.kernelsize - i - 1) * com.kernelsize + j] - minval) * invrange, 0.5f);
				val[1] = pow((com.kernel[kpixels + (com.kernelsize - i - 1) * com.kernelsize + j] - minval) * invrange, 0.5f);
				val[2] = pow((com.kernel[(com.kernelsize - i - 1) * com.kernelsize + j] - minval) * invrange, 0.5f);
			}
			buf[i * stride + 4 * j + 0] = float_to_uchar_range(val[0]);
			buf[i * stride + 4 * j + 1] = float_to_uchar_range(val[1]);
			buf[i * stride + 4 * j + 2] = float_to_uchar_range(val[2]);
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
	siril_free(buf);
	g_mutex_unlock(&psf_preview_mutex);
	return;
}

// PSF drawing callback
gboolean on_PSFkernel_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	if (get_imageorientation() != args.kernelorientation) {
		gtk_widget_set_tooltip_text(GTK_WIDGET(bdeconv_drawingarea), _("Image row order has changed since PSF was generated. The PSF orientation will be updated to match before applying to this image."));
	} else {
		gtk_widget_set_tooltip_text(GTK_WIDGET(bdeconv_drawingarea), "");
	}
	drawing_the_PSF(widget, cr);
	return FALSE;
}

// caller
gboolean DrawPSF(gpointer user_data) {
	if (GTK_IS_WIDGET(bdeconv_drawingarea))
		gtk_widget_queue_draw(GTK_WIDGET(bdeconv_drawingarea));
	return FALSE;
}
