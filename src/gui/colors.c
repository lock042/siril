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

#include <assert.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/colors.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "algos/statistics.h"

#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/registration_preview.h"
#include "gui/message_dialog.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "gui/colors.h"

void on_button_bkg_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };

	if ((!com.selection.h) || (!com.selection.w)) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_h"));
	}

	gtk_spin_button_set_value(selection_black_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_black_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_black_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_black_value[3], com.selection.h);
}

void initialize_calibration_interface() {
	static GtkAdjustment *selection_black_adjustment[4] = { NULL, NULL, NULL,
		NULL };
	static GtkAdjustment *selection_white_adjustment[4] = { NULL, NULL, NULL,
		NULL };

	if (!selection_black_adjustment[0]) {
		selection_black_adjustment[0] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_bkg_x"));
		selection_black_adjustment[1] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_bkg_y"));
		selection_black_adjustment[2] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_bkg_w"));
		selection_black_adjustment[3] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_bkg_h"));
	}
	if (!selection_white_adjustment[0]) {
		selection_white_adjustment[0] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_white_x"));
		selection_white_adjustment[1] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_white_y"));
		selection_white_adjustment[2] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_white_w"));
		selection_white_adjustment[3] = GTK_ADJUSTMENT(
				gtk_builder_get_object(gui.builder, "adjustment_white_h"));
	}
	gtk_adjustment_set_upper(selection_black_adjustment[0], gfit->rx);
	gtk_adjustment_set_upper(selection_black_adjustment[1], gfit->ry);
	gtk_adjustment_set_upper(selection_black_adjustment[2], gfit->rx);
	gtk_adjustment_set_upper(selection_black_adjustment[3], gfit->ry);
	gtk_adjustment_set_value(selection_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_black_adjustment[3], 0);

	gtk_adjustment_set_upper(selection_white_adjustment[0], gfit->rx);
	gtk_adjustment_set_upper(selection_white_adjustment[1], gfit->ry);
	gtk_adjustment_set_upper(selection_white_adjustment[2], gfit->rx);
	gtk_adjustment_set_upper(selection_white_adjustment[3], gfit->ry);
	gtk_adjustment_set_value(selection_white_adjustment[0], 0);
	gtk_adjustment_set_value(selection_white_adjustment[1], 0);
	gtk_adjustment_set_value(selection_white_adjustment[2], 0);
	gtk_adjustment_set_value(selection_white_adjustment[3], 0);
}

void on_button_bkg_neutralization_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	rectangle black_selection;
	int width, height;

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(gui.builder, "spin_bkg_h"));
	}
	width = (int) gtk_spin_button_get_value(selection_black_value[2]);
	height = (int) gtk_spin_button_get_value(selection_black_value[3]);

	if ((!width) || (!height)) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}
	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);

	undo_save_state(gfit, _("Background neutralization"));

	set_cursor_waiting(TRUE);
	background_neutralize(gfit, black_selection);
	populate_roi();
	delete_selected_area();

	update_gfit_histogram_if_needed();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
}

void on_button_white_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_white_value[4] = { NULL, NULL, NULL, NULL };

	if (!selection_white_value[0]) {
		selection_white_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_white_x"));
		selection_white_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_white_y"));
		selection_white_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_white_w"));
		selection_white_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_white_h"));
	}

	if ((!com.selection.h) || (!com.selection.w)) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the white reference area"));
		return;
	}

	gtk_spin_button_set_value(selection_white_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_white_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_white_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_white_value[3], com.selection.h);
}

static void white_balance(fits *fit, gboolean is_manual, rectangle white_selection,
		rectangle black_selection) {
	int chan;
	double norm, low, high;
	double kw[3] = { 0.0, 0.0, 0.0 };
	double bg[3] = { 0.0, 0.0, 0.0 };
	static GtkRange *scale_white_balance[3] = { NULL, NULL, NULL };
	static GtkRange *scaleLimit[2] = { NULL, NULL };

	if (scale_white_balance[RLAYER] == NULL) {
		scale_white_balance[RLAYER] = GTK_RANGE(lookup_widget("scale_r"));
		scale_white_balance[GLAYER] = GTK_RANGE(lookup_widget("scale_g"));
		scale_white_balance[BLAYER] = GTK_RANGE(lookup_widget("scale_b"));

		scaleLimit[0] = GTK_RANGE(lookup_widget("lowWhiteColorCalibScale"));
		scaleLimit[1] = GTK_RANGE(lookup_widget("upWhiteColorCalibScale"));
	}

	assert(fit->naxes[2] == 3);
	norm = get_normalized_value(fit);

	if (is_manual) {
		kw[RLAYER] = gtk_range_get_value(scale_white_balance[RLAYER]);
		kw[GLAYER] = gtk_range_get_value(scale_white_balance[GLAYER]);
		kw[BLAYER] = gtk_range_get_value(scale_white_balance[BLAYER]);

	} else {
		low = gtk_range_get_value(scaleLimit[0]);
		high = gtk_range_get_value(scaleLimit[1]);
		get_coeff_for_wb(fit, white_selection, black_selection, kw, bg, norm, low, high);
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(chan) schedule(static)
#endif
	for (chan = 0; chan < 3; chan++) {
		if (kw[chan] == 1.0) continue;
		calibrate(fit, chan, kw[chan], bg[chan], norm);
	}

	invalidate_stats_from_fit(fit);
	invalidate_gfit_histogram();
}

void on_calibration_apply_button_clicked(GtkButton *button, gpointer user_data) {
	rectangle black_selection, white_selection;
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	static GtkSpinButton *selection_white_value[4] = { NULL, NULL, NULL, NULL };
	struct timeval t_start, t_end;

	siril_log_color_message(_("Color Calibration: processing...\n"), "green");
	gettimeofday(&t_start, NULL);

	GtkToggleButton *manual = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_manual_calibration"));
	gboolean is_manual = gtk_toggle_button_get_active(manual);

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_bkg_h"));
	}

	if (!selection_white_value[0]) {
		selection_white_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_white_x"));
		selection_white_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_white_y"));
		selection_white_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_white_w"));
		selection_white_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_white_h"));
	}

	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);

	if ((!black_selection.w || !black_selection.h) && !is_manual) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}

	white_selection.x = gtk_spin_button_get_value(selection_white_value[0]);
	white_selection.y = gtk_spin_button_get_value(selection_white_value[1]);
	white_selection.w = gtk_spin_button_get_value(selection_white_value[2]);
	white_selection.h = gtk_spin_button_get_value(selection_white_value[3]);

	if ((!white_selection.w || !white_selection.h) && !is_manual) {
		siril_message_dialog( GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the white reference area"));
		return;
	}

	set_cursor_waiting(TRUE);
	undo_save_state(gfit, _("Color Calibration"));
	white_balance(gfit, is_manual, white_selection, black_selection);

	gettimeofday(&t_end, NULL);

	show_time(t_start, t_end);

	populate_roi();
	delete_selected_area();

	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	update_gfit_histogram_if_needed();
	set_cursor_waiting(FALSE);
}

void on_calibration_close_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("color_calibration");
}

gboolean calibration_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("color_calibration");
	return TRUE;
}

void on_checkbutton_manual_calibration_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	GtkWidget *cc_box_red = lookup_widget("cc_box_red");
	GtkWidget *scale_r = lookup_widget("scale_r");
	GtkWidget *cc_box_green = lookup_widget("cc_box_green");
	GtkWidget *scale_g = lookup_widget("scale_g");
	GtkWidget *cc_box_blue = lookup_widget("cc_box_blue");
	GtkWidget *scale_b = lookup_widget("scale_b");
	gtk_widget_set_sensitive(cc_box_red, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(scale_r, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(cc_box_green, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(scale_g, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(cc_box_blue, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(scale_b, gtk_toggle_button_get_active(togglebutton));
}

void negative_processing() {
	set_cursor_waiting(TRUE);
	undo_save_state(gfit, _("Negative Transformation"));
	siril_log_color_message(_("Negative Transformation\n"), "green");
	pos_to_neg(gfit);
	invalidate_stats_from_fit(gfit);
	invalidate_gfit_histogram();
	update_gfit_histogram_if_needed();
	redraw(REMAP_ALL);
	gui_function(redraw_previews, NULL);
	set_cursor_waiting(FALSE);
}
/**********************************************************************/

void on_extract_channel_button_close_clicked(GtkButton *button,
		gpointer user_data) {
	siril_close_dialog("extract_channel_dialog");
}

void on_combo_extract_colors_changed(GtkComboBox *box, gpointer user_data) {
	switch(gtk_combo_box_get_active(box)) {
		default:
		case 0: // RGB
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c1")), _("Red: "));
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c2")), _("Green: "));
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c3")), _("Blue: "));
			break;
		case 1: // HSL
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c1")), _("Hue: "));
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c2")), _("Saturation: "));
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c3")), _("Lightness: "));
			break;
		case 2: // HSV
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c1")), _("Hue: "));
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c2")), _("Saturation: "));
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c3")), _("Value: "));
			break;
		case 3: // CIE L*a*b*
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c1")), "L*: ");
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c2")), "a*: ");
			gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c3")), "b*: ");
	}
}

void on_extract_channel_button_ok_clicked(GtkButton *button, gpointer user_data) {
	static GtkEntry *channel_extract_entry[3] = { NULL, NULL, NULL };
	static GtkComboBox *combo_extract_channel = NULL;

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	struct extract_channels_data *args = calloc(1, sizeof(struct extract_channels_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}

	if (combo_extract_channel == NULL) {
		combo_extract_channel = GTK_COMBO_BOX(lookup_widget("combo_extract_colors"));
		channel_extract_entry[0] = GTK_ENTRY(lookup_widget("Ch1_extract_channel_entry"));
		channel_extract_entry[1] = GTK_ENTRY(lookup_widget("Ch2_extract_channel_entry"));
		channel_extract_entry[2] = GTK_ENTRY(lookup_widget("Ch3_extract_channel_entry"));
	}

	args->type = gtk_combo_box_get_active(combo_extract_channel);
	args->str_type = gtk_combo_box_get_active_id(combo_extract_channel);

	if (args->type != EXTRACT_RGB) {
		// Not RGB, so we need to value_check the image to avoid out-of-range pixels
		if (!value_check(gfit)) {
			siril_log_color_message(_("Error in value_check(). This should not happen...\n"), "red");
			return;
		}
	}

	args->channel[0] = args->channel[1] = args->channel[2] = NULL;

	for (int i = 0; i < 3; i++) {
	    const gchar *text = gtk_entry_get_text(channel_extract_entry[i]);
	    if (text && *text) {
	        args->channel[i] = g_strdup_printf("%s%s", text, com.pref.ext);
	    }
	}

	args->fit = calloc(1, sizeof(fits));
	set_cursor_waiting(TRUE);
	if (copyfits(gfit, args->fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
		siril_log_message(_("Could not copy the input image, aborting.\n"));
		clearfits(args->fit);
		free(args->fit);
		free(args->channel[0]);
		free(args->channel[1]);
		free(args->channel[2]);
		free(args);
	} else {
		copy_fits_metadata(gfit, args->fit);
		if (!start_in_new_thread(extract_channels, args)) {
			clearfits(args->fit);
			free(args->fit);
			free(args->channel[0]);
			free(args->channel[1]);
			free(args->channel[2]);
			free(args);
		}
	}
}

void update_button_sensitivity(GtkWidget *entry, gpointer user_data) {
    GtkWidget *button = GTK_WIDGET(user_data);
    GtkEntry *channel_extract_entry[3] = {
        GTK_ENTRY(lookup_widget("Ch1_extract_channel_entry")),
        GTK_ENTRY(lookup_widget("Ch2_extract_channel_entry")),
        GTK_ENTRY(lookup_widget("Ch3_extract_channel_entry"))
    };

    gboolean has_text = FALSE;

    for (int i = 0; i < 3; i++) {
        const gchar *text = gtk_entry_get_text(channel_extract_entry[i]);
        if (text && *text) {
            has_text = TRUE;
            break;
        }
    }

    gtk_widget_set_sensitive(button, has_text);
}


/* Helper function to read matrix from GUI */
static void get_ccm_values(ccm matrix, float *power) {
	matrix[0][0] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m00"))), NULL);
	matrix[0][1] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m01"))), NULL);
	matrix[0][2] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m02"))), NULL);
	matrix[1][0] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m10"))), NULL);
	matrix[1][1] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m11"))), NULL);
	matrix[1][2] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m12"))), NULL);
	matrix[2][0] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m20"))), NULL);
	matrix[2][1] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m21"))), NULL);
	matrix[2][2] = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(lookup_widget("entry_m22"))), NULL);

	if (power) {
		*power = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spin_ccm_power")));
	}
}

void on_ccm_apply_clicked(GtkButton* button, gpointer user_data) {
	struct ccm_data *args = new_ccm_data();
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}

	get_ccm_values(args->matrix, &args->power);

	GtkToggleButton *btn = GTK_TOGGLE_BUTTON(lookup_widget("check_apply_seq_ccm"));
	gboolean seq_toggle = gtk_toggle_button_get_active(btn);

	if (seq_toggle && sequence_is_loaded()) {
		GtkEntry *ccmSeqEntry = GTK_ENTRY(lookup_widget("entryCCMSeq"));
		args->seqEntry = strdup(gtk_entry_get_text(ccmSeqEntry));
		args->seq = &com.seq;
		apply_ccm_to_sequence(args);
	} else {
		if (seq_toggle) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No sequence is loaded"));
			free_ccm_data(args);
			free(args);
			return;
		}

		// Check for ICC profile warning
		if (gfit->icc_profile && gfit->color_managed) {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("ICC Profile"),
				_("This image has an attached ICC profile. Applying the CCM will invalidate the "
				"ICC profile therefore color management will be disabled. When you have completed low-level color manipulation and returned the image "
				"to the color space described by its ICC profile you can re-enable it using the button at the bottom of this dialog."));
			color_manage(gfit, FALSE);
			gtk_widget_set_sensitive(lookup_widget("ccm_restore_icc"), TRUE);
		}

		// Free the args structure as we're using the worker
		ccm temp_matrix;
		float temp_power = args->power;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				temp_matrix[i][j] = args->matrix[i][j];
			}
		}
		free_ccm_data(args);

		// Process using worker
		set_cursor_waiting(TRUE);
		if (ccm_process_with_worker(temp_matrix, temp_power) == 0) {
			invalidate_stats_from_fit(gfit);
		}
		set_cursor_waiting(FALSE);
	}
}

void on_ccm_restore_icc_clicked(GtkButton *button, gpointer user_data) {
	if (gfit->icc_profile) {
		color_manage(gfit, TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(button), FALSE);
	}
}

static void update_ccm_matrix(ccm matrix) {
	gchar buf[G_ASCII_DTOSTR_BUF_SIZE+1];
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m00")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[0][0]));
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m01")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[0][1]));
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m02")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[0][2]));
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m10")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[1][0]));
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m11")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[1][1]));
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m12")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[1][2]));
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m20")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[2][0]));
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m21")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[2][1]));
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_m22")), g_ascii_dtostr(buf, G_ASCII_DTOSTR_BUF_SIZE, matrix[2][2]));
}

void on_ccm_reset_clicked(GtkButton* button, gpointer user_data) {
	ccm matrix = { { 0.f } };
	matrix[0][0] = matrix[1][1] = matrix[2][2] = 1.0f;
	update_ccm_matrix(matrix);
}

void on_combo_ccm_preset_changed(GtkComboBox *combo, gpointer user_data) {
	ccm matrix;
	float power;
	int index = gtk_combo_box_get_active(combo);
	switch (index) {
		case 0:
			// Custom - does nothing
			return;
		case 1: // Linear Rec.709 to XYZ
			matrix[0][0] = 0.4124564f;
			matrix[0][1] = 0.3575761f;
			matrix[0][2] = 0.1804375f;
			matrix[1][0] = 0.2126729f;
			matrix[1][1] = 0.7151522f;
			matrix[1][2] = 0.0721750f;
			matrix[2][0] = 0.0193339f;
			matrix[2][1] = 0.1191920f;
			matrix[2][2] = 0.9503041f;
			power = 1.f;
			break;
		case 3:
			matrix[0][0] = 0.4124564f;
			matrix[0][1] = 0.3575761f;
			matrix[0][2] = 0.1804375f;
			matrix[1][0] = 0.2126729f;
			matrix[1][1] = 0.7151522f;
			matrix[1][2] = 0.0721750f;
			matrix[2][0] = 0.0193339f;
			matrix[2][1] = 0.1191920f;
			matrix[2][2] = 0.9503041f;
			power = 2.2f;
			break;
		case 2: // XYZ to Linear Rec.709
			matrix[0][0] = 3.2404542f;
			matrix[0][1] = -1.5371385f;
			matrix[0][2] = -0.4985314f;
			matrix[1][0] = -0.9692660f;
			matrix[1][1] = 1.8760108f;
			matrix[1][2] = 0.0415560f;
			matrix[2][0] = 0.0556434f;
			matrix[2][1] = -0.2040259f;
			matrix[2][2] = 1.0572253f;
			power = 1.f;
			break;
		case 4:
			matrix[0][0] = 3.2404542f;
			matrix[0][1] = -1.5371385f;
			matrix[0][2] = -0.4985314f;
			matrix[1][0] = -0.9692660f;
			matrix[1][1] = 1.8760108f;
			matrix[1][2] = 0.0415560f;
			matrix[2][0] = 0.0556434f;
			matrix[2][1] = -0.2040259f;
			matrix[2][2] = 1.0572253f;
			power = 1.f / 2.2f;
			break;
		default:
			siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), _("This case is not handled yet"));
			return;
	}
	update_ccm_matrix(matrix);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_ccm_power")), power);
}

void on_ccm_close_clicked(GtkButton* button, gpointer user_data) {
	siril_close_dialog("ccm_dialog");
}

gboolean ccm_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("ccm_dialog");
	return TRUE;
}
