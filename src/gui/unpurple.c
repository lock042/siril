 /*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril_free-astro.org/index.php/Siril
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
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "io/image_format_fits.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"
#include "opencv/opencv.h"
#include "filters/synthstar.h"
#include "filters/unpurple.h"

static double mod_b = 1.0, thresh = 0.0, old_thresh = 0.0;
static fits starmask = {0};
static gboolean is_roi = FALSE;

int generate_binary_starmask(fits *fit, fits **star_mask, double threshold) {
	gboolean stars_needs_freeing = FALSE;
	psf_star **stars = NULL;
	int channel = 1;

	int nb_stars = starcount(com.stars);
	int dimx = fit->naxes[0];
	int dimy = fit->naxes[1];
	int count = dimx * dimy;

	// Do we have stars from Dynamic PSF or not?
	if (nb_stars < 1) {
		image *input_image = NULL;
		input_image = siril_calloc(1, sizeof(image));
		input_image->fit = fit;
		input_image->from_seq = NULL;
		input_image->index_in_seq = -1;
		stars = peaker(input_image, channel, &com.pref.starfinder_conf, &nb_stars,
						NULL, FALSE, FALSE, 0, com.pref.starfinder_conf.profile, com.max_thread);
		siril_free(input_image);
		stars_needs_freeing = TRUE;
	} else {
		stars = com.stars;
		stars_needs_freeing = FALSE;
	}

	if (starcount(stars) < 1) {
		siril_log_color_message(_("No stars detected in the image.\n"), "red");
		return -1;
	}

	siril_log_message(_("Creating binary star mask for %d stars...\n"), nb_stars);
	if (new_fit_image(star_mask, dimx, dimy, 1, DATA_USHORT)) {
		return -1;
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread) if (com.max_thread > 1)
#endif
	for (size_t i = 0; i < dimx * dimy; i++) {
		(*star_mask)->pdata[0][i] = 0;
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread)
#endif
	for (int n = 0; n < nb_stars; n++) {
		int size = (int) 2 * max(stars[n]->fwhmx, stars[n]->fwhmy);

		// The fringe factor tries to adjust the star size to allow for the purple fringe
		// The threshold slider decides the scaling factor
		double base_fringe = 10.0;
		double fringe = base_fringe + threshold * pow(size, 1.5);
		size = sqrt(pow(size, 2) + pow(fringe, 2));

		if (size % 2 == 0)
			size++;

		int x0 = (int)(stars[n]->xpos - size / 2);
		int y0 = (int)((dimy - stars[n]->ypos) - size / 2);
		for (int y = 0; y < size; y++) {
			for (int x = 0; x < size; x++) {
				int px = x0 + x;
				int py = y0 + y;
				if (px >= 0 && px < dimx && py >= 0 && py < dimy) {
					int idx = py * dimx + px;
					if (idx >= 0 && idx < count) {
						double dx = x - (size / 2.0);
						double dy = y - (size / 2.0);
						double distance = sqrt(dx * dx + dy * dy);
						int is_star = (distance <= (size / 2.0)) ? 1 : 0;
						(*star_mask)->pdata[0][idx] = is_star ? USHRT_MAX : (*star_mask)->pdata[0][idx];
					}
				}
			}
		}
	}

	if (stars_needs_freeing)
		free_fitted_stars(stars);

	return 0;
}

static int unpurple_update_preview() {
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview"))))
		copy_backup_to_gfit();

	gboolean withstarmask = TRUE;
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_stars"))))
		withstarmask = FALSE;

	fits *fit = gui.roi.active ? &gui.roi.fit : &gfit;

	// We need to create a star mask if one doesn't already exist
	// and we need to recreate it if we change roi
	if (withstarmask) {
		if (starmask.naxis == 0 || gui.roi.active != is_roi || old_thresh != thresh) {
			is_roi = gui.roi.active;
			old_thresh = thresh;
			fits *starmask_ptr = &starmask;
			generate_binary_starmask(fit, &starmask_ptr, thresh);
		}
	}

	struct unpurpleargs *args = siril_calloc(1, sizeof(struct unpurpleargs));
	*args = (struct unpurpleargs){.fit = fit, .starmask = &starmask, .withstarmask = withstarmask, .thresh = thresh, .mod_b = mod_b, .verbose = FALSE, .for_final = FALSE};
	set_cursor_waiting(TRUE);
	// we call the unpurple_filter directly here because update_preview already handles the ROI mutex lock
	if (!start_in_new_thread(unpurple_filter, args)) {
		siril_free(args);
	}
	return 0;
}

void unpurple_change_between_roi_and_image() {
	// If we are showing the preview, update it after the ROI change.
	gui.roi.operation_supports_roi = TRUE;
	roi_supported(TRUE);
	update_image *param = siril_malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

static void unpurple_startup() {
	add_roi_callback(unpurple_change_between_roi_and_image);
	roi_supported(TRUE);
	copy_gfit_to_backup();
}

static void unpurple_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(), _("Unpurple filter: (thresh=%2.2lf, mod_b=%2.2lf)"), thresh, mod_b);
	}
	roi_supported(FALSE);
	remove_roi_callback(unpurple_change_between_roi_and_image);
	clearfits(&starmask);
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int unpurple_process_all() {
	set_cursor_waiting(TRUE);
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview"))))
		copy_backup_to_gfit();
	fits *fit = &gfit;

	gboolean withstarmask = TRUE;
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_stars"))))
		withstarmask = FALSE;

	//TODO: Optimization: Can we reuse the starmask we already have?
	if (withstarmask) {
		fits *starmask_ptr = &starmask;
		generate_binary_starmask(fit, &starmask_ptr, thresh);
	}

	struct unpurpleargs *args = siril_calloc(1, sizeof(struct unpurpleargs));
	*args = (struct unpurpleargs){.fit = fit, .starmask = &starmask, .withstarmask = withstarmask, .thresh = thresh, .mod_b = mod_b, .verbose = FALSE, .for_final = TRUE};

	// We call the unpurple handler here because we don't have update_preview to handle the ROI mutex for us
	if (!start_in_new_thread(unpurple_handler, args)) {
		siril_free(args);
	}

	return 0;
}

static void apply_unpurple_changes() {
	gboolean status = (mod_b != 1.0) || (thresh != 0.0);
	unpurple_close(!status);
}

void apply_unpurple_cancel() {
	unpurple_close(TRUE);
	siril_close_dialog("unpurple_dialog");
}

/*** callbacks **/

void on_unpurple_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_unpurple_mod_b = GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_mod_b"));
	GtkSpinButton *spin_unpurple_thresh = GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_thresh"));

	unpurple_startup();
	mod_b = 1.0, thresh = 0.0;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_unpurple_mod_b, mod_b);
	gtk_spin_button_set_value(spin_unpurple_thresh, thresh);

	set_notify_block(FALSE);

	// Default parameters transform the image, so update the preview if toggle is active
	update_image *param = siril_malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

void on_unpurple_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_unpurple_cancel();
}

void on_unpurple_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview"))) || gui.roi.active) {
		unpurple_process_all();
	}

	apply_unpurple_changes();
	siril_close_dialog("unpurple_dialog");
}

void on_unpurple_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_unpurple_changes();
}

void on_unpurple_undo_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_unpurple_mod_b = GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_mod_b"));
	GtkSpinButton *spin_unpurple_thresh = GTK_SPIN_BUTTON(lookup_widget("spin_unpurple_thresh"));

	mod_b = 1.0, thresh = 0.0;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(spin_unpurple_mod_b, mod_b);
	gtk_spin_button_set_value(spin_unpurple_thresh, thresh);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = siril_malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

/*** adjusters **/
void on_spin_unpurple_mod_b_value_changed(GtkSpinButton *button, gpointer user_data) {
	mod_b = gtk_spin_button_get_value(button);
	update_image *param = siril_malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

void on_spin_unpurple_thresh_value_changed(GtkSpinButton *button, gpointer user_data) {
	thresh = gtk_spin_button_get_value(button);
	update_image *param = siril_malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}

void on_unpurple_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	cancel_pending_update();
	if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")))) {
		/* if user click very fast */
		waiting_for_thread();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();
		update_image *param = siril_malloc(sizeof(update_image));
		param->update_preview_fn = unpurple_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer)param);
	}
}

void on_unpurple_stars_toggled(GtkToggleButton *button, gpointer user_data) {
	update_image *param = siril_malloc(sizeof(update_image));
	param->update_preview_fn = unpurple_update_preview;
	param->show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("unpurple_preview")));
	notify_update((gpointer)param);
}
