/*
 * This file is part of Siril, an astronomy image processor.
 *
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#include <math.h>

#include "algos/star_finder.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "io/single_image.h"
#include "io/sequence.h"


void on_toggle_relax_checks_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.pref.starfinder_conf.relax_checks = gtk_toggle_button_get_active(togglebutton);
}

void on_spin_sf_radius_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.pref.starfinder_conf.radius = (int)gtk_spin_button_get_value(spinbutton);
}

void on_spin_sf_threshold_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.pref.starfinder_conf.sigma = gtk_spin_button_get_value(spinbutton);
}

void on_spin_sf_roundness_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	gboolean range = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("psf_roundness_range_check_button")));
	double minr = gtk_spin_button_get_value(spinbutton);
	if (range) {
		double maxr = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_maxr")));
		if (minr >= maxr - 0.01)
			gtk_spin_button_set_value(spinbutton, maxr - 0.01);
	}
	com.pref.starfinder_conf.roundness = minr;
}

void on_spin_sf_minbeta_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.pref.starfinder_conf.min_beta = gtk_spin_button_get_value(spinbutton);
}

void on_spin_sf_convergence_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.pref.starfinder_conf.convergence = (int)gtk_spin_button_get_value(spinbutton);
}

void on_spin_sf_minA_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	double maxA = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_maxA")));
	double minA = gtk_spin_button_get_value(spinbutton);
	if (minA >= maxA - 0.01) {
		gtk_spin_button_set_value(spinbutton, maxA - 0.01);
	}
	com.pref.starfinder_conf.min_A = gtk_spin_button_get_value(spinbutton);
	siril_debug_print("minA = %f, maxA = %f\n", com.pref.starfinder_conf.min_A, com.pref.starfinder_conf.max_A);
}

void on_spin_sf_maxA_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	double minA = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_minA")));
	double maxA = gtk_spin_button_get_value(spinbutton);
	if (maxA <= minA + 0.01) {
		gtk_spin_button_set_value(spinbutton, minA + 0.01);
	}
	com.pref.starfinder_conf.max_A = gtk_spin_button_get_value(spinbutton);
	siril_debug_print("minA = %f, maxA = %f\n", com.pref.starfinder_conf.min_A, com.pref.starfinder_conf.max_A);
}

void on_spin_sf_maxr_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	double minr = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round")));
	double maxr = gtk_spin_button_get_value(spinbutton);
	if (maxr <= minr + 0.01) {
		gtk_spin_button_set_value(spinbutton, minr + 0.01);
	}
	com.pref.starfinder_conf.max_r = gtk_spin_button_get_value(spinbutton);
	siril_debug_print("minr = %f, maxr = %f\n", com.pref.starfinder_conf.roundness, com.pref.starfinder_conf.max_r);
}

void on_psf_amplitude_range_check_button_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *minw = lookup_widget("spinstarfinder_minA"), *maxw = lookup_widget("spinstarfinder_maxA");
	gboolean enabled = gtk_toggle_button_get_active(togglebutton);
	if (enabled) {
		com.pref.starfinder_conf.min_A = gtk_spin_button_get_value(GTK_SPIN_BUTTON(minw));
		com.pref.starfinder_conf.max_A = gtk_spin_button_get_value(GTK_SPIN_BUTTON(maxw));
	} else {
		com.pref.starfinder_conf.min_A = 0.0;
		com.pref.starfinder_conf.max_A = 0.0;
	}
	gtk_widget_set_sensitive(minw, enabled);
	gtk_widget_set_sensitive(maxw, enabled);
	siril_debug_print("minA = %f, maxA = %f\n", com.pref.starfinder_conf.min_A, com.pref.starfinder_conf.max_A);
}

void on_psf_roundness_range_check_button_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *maxw = lookup_widget("spinstarfinder_maxr");
	gboolean enabled = gtk_toggle_button_get_active(togglebutton);
	if (enabled) {
		double minr = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round")));
		double maxr = gtk_spin_button_get_value(GTK_SPIN_BUTTON(maxw));
		if (maxr <= minr + 0.01) {
			maxr = minr + 0.01;
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(maxw), maxr);
		}
		com.pref.starfinder_conf.max_r = maxr;
	}
	else com.pref.starfinder_conf.max_r = 1.0;
	gtk_widget_set_sensitive(maxw, enabled);
	siril_debug_print("minr = %f, maxr = %f\n", com.pref.starfinder_conf.roundness, com.pref.starfinder_conf.max_r);
}

void on_combostarfinder_profile_changed(GtkComboBox *combo, gpointer user_data) {
	GtkWidget *beta_control_label = GTK_WIDGET(lookup_widget("beta_control_label"));
	GtkWidget *beta_control_spin = GTK_WIDGET(lookup_widget("spin_minbeta"));
	com.pref.starfinder_conf.profile = gtk_combo_box_get_active(combo);

	gtk_widget_set_visible(beta_control_label, com.pref.starfinder_conf.profile != PSF_GAUSSIAN);
	gtk_widget_set_visible(beta_control_spin, com.pref.starfinder_conf.profile != PSF_GAUSSIAN);
}

void on_reset_findstar_button_clicked(GtkButton *button, gpointer user_data) {
	com.pref.starfinder_conf = (star_finder_params) {
		.radius = DEF_BOX_RADIUS, .sigma = 1., .roundness = 0.5,
			.convergence = 1, .relax_checks = FALSE, .profile = PSF_GAUSSIAN,
			.min_beta = 1.5, .min_A = 0.0, .max_A = 0.0, .max_r = 1.0 };
	update_peaker_GUI();
}

void update_peaker_GUI() {
	static GtkSpinButton *spin_radius = NULL, *spin_sigma = NULL,
			*spin_roundness = NULL, *spin_convergence = NULL, *spin_minbeta = NULL,
			*spin_minA = NULL, *spin_maxA = NULL, *spin_maxR = NULL;
	static GtkToggleButton *toggle_checks = NULL, *toggle_A = NULL, *toggle_R;
	static GtkComboBox *combostarfinder_profile;

	if (spin_radius == NULL) {
		spin_radius = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_radius"));
		spin_sigma = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_threshold"));
		spin_roundness = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round"));
		spin_convergence = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_convergence"));
		spin_minbeta = GTK_SPIN_BUTTON(lookup_widget("spin_minbeta"));
		spin_minA = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_minA"));
		spin_maxA = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_maxA"));
		spin_maxR = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_maxr"));
		toggle_A = GTK_TOGGLE_BUTTON(lookup_widget("psf_amplitude_range_check_button"));
		toggle_R = GTK_TOGGLE_BUTTON(lookup_widget("psf_roundness_range_check_button"));
		toggle_checks = GTK_TOGGLE_BUTTON(lookup_widget("toggle_relax_checks"));
		combostarfinder_profile = GTK_COMBO_BOX(lookup_widget("combostarfinder_profile"));
	}
	gtk_spin_button_set_value(spin_radius, (double) com.pref.starfinder_conf.radius);
	gtk_spin_button_set_value(spin_sigma, com.pref.starfinder_conf.sigma);
	gtk_spin_button_set_value(spin_roundness, com.pref.starfinder_conf.roundness);
	gtk_spin_button_set_value(spin_convergence, com.pref.starfinder_conf.convergence);
	gtk_spin_button_set_value(spin_minbeta, com.pref.starfinder_conf.min_beta);
	//gtk_spin_button_set_value(spin_minA, com.pref.starfinder_conf.min_A);
	//gtk_spin_button_set_value(spin_maxA, com.pref.starfinder_conf.max_A);
	gtk_spin_button_set_value(spin_maxR, com.pref.starfinder_conf.max_r);
	gtk_toggle_button_set_active(toggle_checks, com.pref.starfinder_conf.relax_checks);
	gtk_combo_box_set_active(combostarfinder_profile, com.pref.starfinder_conf.profile);
	// in case the toggles are being reset as well
	if (com.pref.starfinder_conf.max_r == 1.0) {
		gtk_toggle_button_set_active(toggle_R, FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(spin_maxR), FALSE);
	}
	g_signal_handlers_block_by_func(toggle_A, on_psf_amplitude_range_check_button_toggled, NULL);
	if (com.pref.starfinder_conf.min_A == 0.0 && com.pref.starfinder_conf.max_A == 0.0) {
		gtk_toggle_button_set_active(toggle_A, FALSE);
		gtk_spin_button_set_value(spin_maxA, 1.0);
		gtk_widget_set_sensitive(GTK_WIDGET(spin_minA), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(spin_maxA), FALSE);
	} else {
		gtk_toggle_button_set_active(toggle_A, TRUE);
		gtk_spin_button_set_value(spin_minA, com.pref.starfinder_conf.min_A);
		gtk_spin_button_set_value(spin_maxA, com.pref.starfinder_conf.max_A);
		gtk_widget_set_sensitive(GTK_WIDGET(spin_minA), TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(spin_maxA), TRUE);
	}
	g_signal_handlers_unblock_by_func(toggle_A, on_psf_amplitude_range_check_button_toggled, NULL);
}

void confirm_peaker_GUI() {
	static GtkSpinButton *spin_radius = NULL, *spin_sigma = NULL,
			*spin_roundness = NULL, *spin_convergence = NULL, *spin_minbeta = NULL,
			*spin_minA = NULL, *spin_maxA = NULL, *spin_maxR = NULL;

	if (spin_radius == NULL) {
		spin_radius = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_radius"));
		spin_sigma = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_threshold"));
		spin_roundness = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round"));
		spin_convergence = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_convergence"));
		spin_minbeta = GTK_SPIN_BUTTON(lookup_widget("spin_minbeta"));
		spin_minA = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_minA"));
		spin_maxA = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_maxA"));
		spin_maxR = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_maxr"));
	}
	gtk_spin_button_update(spin_radius);
	gtk_spin_button_update(spin_sigma);
	gtk_spin_button_update(spin_roundness);
	gtk_spin_button_update(spin_convergence);
	gtk_spin_button_update(spin_minbeta);
	gtk_spin_button_update(spin_minA);
	gtk_spin_button_update(spin_maxA);
	gtk_spin_button_update(spin_maxR);
}

void on_process_starfinder_button_clicked(GtkButton *button, gpointer user_data) {
	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		siril_log_color_message(_("Load an image first, aborted.\n"), "red");
		return;
	}
	confirm_peaker_GUI(); //making sure the spin buttons values are read even without confirmation

	struct starfinder_data *args = calloc(1, sizeof(struct starfinder_data));
	args->im.fit = &gfit;
	if (sequence_is_loaded() && com.seq.current >= 0) {
		args->im.from_seq = &com.seq;
		args->im.index_in_seq = com.seq.current;
	} else {
		args->im.from_seq = NULL;
		args->im.index_in_seq = -1;
	}
	args->layer = select_vport(gui.cvport);
	args->max_stars_fitted = 0;
	args->starfile = NULL;
	args->threading = MULTI_THREADED;
	args->update_GUI = TRUE;
	if (com.selection.w != 0 && com.selection.h != 0) {
		args->selection = com.selection;
	}

	if (!start_in_new_thread(findstar_worker, args)) {
		free(args);
	}
}

/* Function to add star one by one, from the selection rectangle, the
 * minimization is run and the star is detected and added to the list of stars.
 *
 * IF A STAR IS FOUND and not already present in com.stars, the return value is
 * the new star and index is set to the index of the new star in com.stars.
 * IF NO NEW STAR WAS FOUND, either because it was already in the list, or a
 * star failed to be detected in the selection, or any other error, the return
 * value is NULL and index is set to -1.
 */
psf_star *add_star(fits *fit, int layer, int *index) {
	int i = 0;
	gboolean already_found = FALSE;
	starprofile profile;
	if (com.stars && com.stars[0])
		// Add star depending on fit profile used for existing stars in the array
		profile = (com.stars[0]->profile == PSF_GAUSSIAN ? PSF_GAUSSIAN : PSF_MOFFAT_BFREE);
	else
		// Default to Gaussian (or get this from a parameter in the GUI, tbd)
		profile = com.pref.starfinder_conf.profile;

	*index = -1;
	psf_star *result = psf_get_minimisation(&gfit, layer, &com.selection, FALSE, FALSE, NULL, TRUE, profile, NULL);
	if (!result) // we don't check for errors as we assume the user has selected a star
		return NULL;
	result->angle = -result->angle; // we need to invert the angle because of the way the matrix is passed to minimizer
	/* We do not check if it's matching with the "reject_star()" criteria.
	 * Indeed, in this case the user can add manually stars missed by star_finder */

	if (com.stars && !com.star_is_seqdata) {
		// check if the star was already detected/peaked
		while (com.stars[i]) {
			if (fabs(result->x0 + com.selection.x - com.stars[i]->xpos) < 0.9
					&& fabs(com.selection.y + com.selection.h - result->y0
									- com.stars[i]->ypos) < 0.9)
				already_found = TRUE;
			i++;
		}
	} else {
		if (com.star_is_seqdata) {
			/* com.stars was allocated with a size of 2, we need to free it before reallocating */
			clear_stars_list(TRUE);
		}
		com.stars = new_fitted_stars(MAX_STARS);
		if (!com.stars) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		com.star_is_seqdata = FALSE;
	}

	if (already_found) {
		free_psf(result);
		result = NULL;
		char *msg = siril_log_message(_("This star has already been picked !\n"));
		siril_message_dialog(GTK_MESSAGE_INFO, _("Peaker"), msg);
	} else {
		if (i < MAX_STARS) {
			result->xpos = result->x0 + com.selection.x;
			result->ypos = com.selection.y + com.selection.h - result->y0;
			psf_star **newstars = realloc(com.stars, (i + 2) * sizeof(psf_star *));
			if (!newstars)
				PRINT_ALLOC_ERR;
			else {
				com.stars = newstars;
				com.stars[i] = result;
				com.stars[i + 1] = NULL;
				*index = i;
			}
		} else {
			free_psf(result);
			result = NULL;
		}
	}
	return result;
}

static int get_comstar_count() {
	int i = 0;
	while (com.stars[i])
		i++;
	return i;
}

/* Remove a star from com.stars, at index index. The star is freed. */
int remove_star(int index) {
	if (index < 0 || !com.stars || !com.stars[index])
		return 1;

	int N = get_comstar_count() + 1;

	free_psf(com.stars[index]);
	memmove(&com.stars[index], &com.stars[index + 1],
			(N - index - 1) * sizeof(*com.stars));
	redraw(REDRAW_OVERLAY);
	return 0;
}

gboolean end_findstar(gpointer p) {
	struct starfinder_data *args = (struct starfinder_data *) p;
	stop_processing_thread();
	set_cursor_waiting(FALSE);
	free(args);
	return FALSE;
}

