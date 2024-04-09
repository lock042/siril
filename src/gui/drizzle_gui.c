/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/command.h"
#include "drizzle/driz_portability.h"
#include "drizzle/cdrizzleutil.h"
#include "drizzle/cdrizzlemap.h"
#include "drizzle/cdrizzlebox.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "gui/message_dialog.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"

void on_apply_drizzle_clicked(GtkButton *button, gpointer user_data) {
	struct driz_args_t *driz = calloc(1, sizeof(struct driz_args_t));
	driz->seq = &com.seq;
	driz->reference_image = sequence_find_refimage(&com.seq);
	driz->keep_counts = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("generate_pixcnt")));
//	driz->use_wcs = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_driz_method"))) == 1;
	driz->use_wcs = FALSE;
	driz->load_new_sequence = TRUE;
	driz->use_flats = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("driz_use_flats")));
	driz->scale = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_driz_scale")));
	driz->weight_scale = 1.f; // Not used for now
	driz->kernel = (enum e_kernel_t) gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_driz_kernel")));
	driz->pixel_fraction = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_driz_dropsize")));
	driz->framing = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboreg_framing")));
	driz->prefix = strdup(gtk_entry_get_text(GTK_ENTRY(lookup_widget("regseqname_entry"))));
	get_reg_sequence_filtering_from_gui(&driz->filtering_criterion, &driz->filtering_parameter, -1);
	if (driz->use_flats) {
		fits reffit = { 0 };
		GtkEntry *entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		const gchar *flat_filename = gtk_entry_get_text(entry);
		gchar *error = NULL;
		int status;
		gchar *expression = path_parse(&reffit, flat_filename, PATHPARSE_MODE_READ, &status);
		if (status) {
			error = _("NOT USING FLAT: could not parse the expression");
			driz->use_flats = FALSE;
		} else {
			free(expression);
			if (flat_filename[0] == '\0') {
				siril_log_message(_("Error: no master flat specified in the preprocessing tab.\n"));
				free(driz->prefix);
				free(driz);
				return;
			} else {
				set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
				driz->flat = calloc(1, sizeof(fits));
				if (!readfits(flat_filename, driz->flat, NULL, TRUE)) {
					if (driz->flat->naxes[2] != driz->seq->nb_layers) {
						error = _("NOT USING FLAT: number of channels is different");
					} else if (driz->flat->naxes[0] != driz->seq->rx ||
							driz->flat->naxes[1] != driz->seq->ry) {
						error = _("NOT USING FLAT: image dimensions are different");
					} else {
						// no need to deal with bitdepth conversion as readfits has already forced conversion to float
						siril_log_message(_("Master flat read for use as initial pixel weight\n"));
					}
				} else error = _("NOT USING FLAT: cannot open the file");
				if (error) {
					siril_log_color_message("%s\n", "red", error);
					set_progress_bar_data(error, PROGRESS_DONE);
					if (driz->flat) {
						clearfits(driz->flat);
						free(driz->flat);
					}
					free(driz->prefix);
					free(driz);
					return;
				}
			}
		}
	}

	apply_drizzle(driz);
}

void on_drizzleCheckButton_toggled(GtkToggleButton* button, gpointer user_data) {
	gboolean state = gtk_toggle_button_get_active(button);
	gtk_widget_set_visible(lookup_widget("box_drizzle_controls"), state);
	if (state) {
		gtk_notebook_set_current_page(GTK_NOTEBOOK(lookup_widget("notebook_registration")), 4);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("upscaleCheckButton")), FALSE);
	}
	gtk_widget_set_visible(lookup_widget("interp_box"), !state);
	gtk_widget_set_visible(lookup_widget("toggle_reg_clamp"), !state);
	gtk_widget_set_visible(lookup_widget("regNoOutput"), !state);

}

void on_upscaleCheckButton_toggled(GtkToggleButton* button, gpointer user_data) {
	gboolean state = gtk_toggle_button_get_active(button);
	if (state) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("drizzleCheckButton")), FALSE);
	}
}
