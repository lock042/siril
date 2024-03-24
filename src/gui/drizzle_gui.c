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
#include "gui/message_dialog.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"

void on_drizzleTab_show(GtkWidget *widget, gpointer user_data) {
	if (com.seq.cfa_opened_monochrome) {
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_driz_bayer")), "active");
	} else {
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_driz_bayer")), "inactive");
	}
}

void on_apply_drizzle_clicked(GtkButton *button, gpointer user_data) {
	struct driz_args_t *driz = calloc(1, sizeof(struct driz_args_t));
	driz->seq = &com.seq;
	driz->reference_image = sequence_find_refimage(&com.seq);
	driz->keep_counts = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("generate_pixcnt")));
	driz->use_wcs = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_driz_method"))) == 1;
	driz->use_flats = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("driz_use_flats")));
	driz->use_bias = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("driz_use_bias")));
	driz->scale = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_driz_scale")));
	driz->weight_scale = 1.f; // Not used for now
	driz->kernel = (enum e_kernel_t) gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_driz_kernel")));
	driz->pixel_fraction = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_driz_dropsize")));
	driz->filtering_criterion = seq_filter_all; // Needs a GUI element
	driz->filtering_parameter = 1.0; // Needs a GUI element
	driz->framing = FRAMING_MAX; // for testing purposes, this is probably the one that best shows any shifts / rotation between frames. Needs a GUI element
	apply_drizzle(driz);
}


