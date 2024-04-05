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

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/siril_wcs.h"
#include "algos/photometry.h"
#include "algos/comparison_stars.h"
#include "algos/occult_time.h"
#include "occultation.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "io/sequence.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *delay_cam = NULL;
static GtkWidget *sep = NULL;
static GtkWidget *apply_offset = NULL;
static double delay = 0.0;
//static GtkWidget *delta_vmag_entry = NULL;
//static GtkWidget *delta_bv_entry = NULL;
//static GtkWidget *emag_entry = NULL;
//static GtkWidget *target_entry = NULL;
//static GtkWidget *apass_radio = NULL;
//static GtkWidget *check_narrow = NULL;

static void on_occult_response(GtkDialog* self, gint response_id, gpointer user_data);
static void on_find_clicked(GtkDialog* self, gint response_id, gpointer user_data);

static void build_the_dialog() {
	dialog = gtk_dialog_new_with_buttons(_("Calibrate Timestamps over PPS"), NULL,
			0, _("_OK"), GTK_RESPONSE_ACCEPT, _("_Close"), GTK_RESPONSE_REJECT, NULL);


	// If the user clicks one of these dialog buttons, GtkDialog will emit
	// the GtkDialog::response signal with the corresponding response ID
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 400);
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
	g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
	g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_occult_response), NULL);

	GtkWidget *label = gtk_label_new(_("Calibration of a camera setup with a 1PPS signal trigger"));
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	g_object_set(G_OBJECT(label), "margin", 15, NULL);


	GtkWidget *label2 = gtk_label_new(_("The first step aimes at finding the necessary delay to re-time\nand synchronize the timestamps to a 1PPS signal from a GPS module.\n\nTo get this parameter, you need:\n -a loaded sequence (from a .ser file)\n -this sequence must have a pseudo-star blinking at a beat of 1 pulse per second\n -make a selection around this pseudo-star, including enough background\n\nThen you can either click on Find delay or write your own delay if you know it\n"));
	gtk_label_set_line_wrap(GTK_LABEL(label2), TRUE);
	g_object_set(G_OBJECT(label2), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(label2), "margin-top", 15, NULL);


//	target_entry = gtk_entry_new();
//	gtk_entry_set_placeholder_text(GTK_ENTRY(target_entry), "Target star name");
//	gtk_widget_set_tooltip_text(target_entry, _("Enter the target star name to search in cataloges"));
//	g_object_set(G_OBJECT(target_entry), "margin", 15, NULL);
	GtkWidget *label_delay_cam = gtk_label_new(_("Camera delay to handle (ms):"));
	gtk_widget_set_halign(label_delay_cam, GTK_ALIGN_CENTER);
//	g_object_set(G_OBJECT(label_delay_cam), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(label_delay_cam), "margin-top", 5, NULL);
	g_object_set(G_OBJECT(label_delay_cam), "margin-bottom", 0, NULL);

	delay_cam = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(delay_cam), "0.0");
//	delay_cam.props.xalign = 0.5;
	gtk_entry_set_alignment (GTK_ENTRY (delay_cam), 0.5);
	gtk_widget_set_tooltip_text(delay_cam, _("Camera delay"));
	gtk_widget_set_halign(delay_cam, GTK_ALIGN_CENTER);
//	g_object_set(G_OBJECT(delay_cam), "margin-left", 35, NULL);
//	g_object_set(G_OBJECT(delay_cam), "margin-right", 40, NULL);
//	g_object_set(G_OBJECT(delay_cam), "margin-top", 0, NULL);

	GtkWidget *find_delay = gtk_button_new_with_mnemonic (_("Find Delay"));
	g_signal_connect (find_delay, "clicked", G_CALLBACK (on_find_clicked), NULL);

	sep = gtk_separator_new (GTK_ORIENTATION_HORIZONTAL);
	gtk_widget_set_size_request (sep, 1, 5);
	g_object_set(G_OBJECT(sep), "margin-top", 15, NULL);

	GtkWidget *label3 = gtk_label_new(_("Once the delay is set, you can choose to apply it or not \n"));
//	gtk_label_set_line_wrap(GTK_LABEL(label3), TRUE);
	gtk_widget_set_halign(label3, GTK_ALIGN_START);
	g_object_set(G_OBJECT(label3), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(label3), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(label3), "margin-top", 15, NULL);

	apply_offset = gtk_check_button_new_with_label(_("Apply Offset"));
	gtk_widget_set_halign(apply_offset, GTK_ALIGN_CENTER);
	g_object_set(G_OBJECT(apply_offset), "margin-bottom", 35, NULL);
	g_object_set(G_OBJECT(apply_offset), "margin-top", 5, NULL);

/*	check_narrow = gtk_check_button_new_with_label(_("Narrow field of view"));
	gtk_widget_set_tooltip_text(check_narrow, _("Tick this box to use a narrow field of view centered about the target star"));
	gtk_widget_set_halign(check_narrow, GTK_ALIGN_START);
	g_object_set(G_OBJECT(check_narrow), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(check_narrow), "margin-top", 5, NULL);
	g_object_set(G_OBJECT(check_narrow), "margin-bottom", 0, NULL);

	GtkWidget *labelvmag = gtk_label_new(_("Allowed visual magnitude range:"));
	gtk_widget_set_halign(labelvmag, GTK_ALIGN_START);
	g_object_set(G_OBJECT(labelvmag), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(labelvmag), "margin-top", 10, NULL);
	g_object_set(G_OBJECT(labelvmag), "margin-bottom", 0, NULL);

	delta_vmag_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(delta_vmag_entry), "3.0");
	gtk_widget_set_tooltip_text(delta_vmag_entry, _("Allowed range of visual magnitude between the target star and the comparison stars"));
	g_object_set(G_OBJECT(delta_vmag_entry), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(delta_vmag_entry), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(delta_vmag_entry), "margin-top", 0, NULL);

	GtkWidget *labelbv = gtk_label_new(_("Allowed B-V index range:"));
	gtk_widget_set_halign(labelbv, GTK_ALIGN_START);
	g_object_set(G_OBJECT(labelbv), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(labelbv), "margin-top", 10, NULL);
	g_object_set(G_OBJECT(labelbv), "margin-bottom", 0, NULL);

	delta_bv_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(delta_bv_entry), "0.5");
	gtk_widget_set_tooltip_text(delta_bv_entry, _("Allowed range of B-V index (color) between the target star and the comparison stars"));
	g_object_set(G_OBJECT(delta_bv_entry), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(delta_bv_entry), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(delta_bv_entry), "margin-top", 0, NULL);

	GtkWidget *labelemag = gtk_label_new(_("Allowed magnitude error:"));
	gtk_widget_set_halign(labelemag, GTK_ALIGN_START);
	g_object_set(G_OBJECT(labelemag), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(labelemag), "margin-top", 10, NULL);
	g_object_set(G_OBJECT(labelemag), "margin-bottom", 0, NULL);

	emag_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(emag_entry), "0.03");
	gtk_widget_set_tooltip_text(emag_entry, _("Allowed catalogue magnitude error for comparison stars"));
	g_object_set(G_OBJECT(emag_entry), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(emag_entry), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(emag_entry), "margin-top", 0, NULL);
*/
	/* catalogue choice */
/*	GtkWidget *radio2, *radiobox;
	radiobox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(radiobox), TRUE);
	gtk_widget_set_tooltip_text(radiobox, _("Recommended catalogue for this feature is APASS"));

	apass_radio = gtk_radio_button_new_with_label_from_widget(NULL, _("APASS catalogue"));
	radio2 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(apass_radio),
			_("NOMAD catalogue"));
	gtk_container_add(GTK_CONTAINER(radiobox), apass_radio);
	gtk_container_add(GTK_CONTAINER(radiobox), radio2);
	g_object_set(G_OBJECT(radiobox), "margin", 15, NULL);

	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_box_set_spacing(GTK_BOX(content_area), 15);
	gtk_container_add(GTK_CONTAINER(content_area), label);
	gtk_container_add(GTK_CONTAINER(content_area), target_entry);
	gtk_container_add(GTK_CONTAINER(content_area), check_narrow);
	gtk_container_add(GTK_CONTAINER(content_area), labelvmag);
	gtk_container_add(GTK_CONTAINER(content_area), delta_vmag_entry);
	gtk_container_add(GTK_CONTAINER(content_area), labelbv);
	gtk_container_add(GTK_CONTAINER(content_area), delta_bv_entry);
	gtk_container_add(GTK_CONTAINER(content_area), labelemag);
	gtk_container_add(GTK_CONTAINER(content_area), emag_entry);
	gtk_container_add(GTK_CONTAINER(content_area), radiobox);
*/
	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), label);
	gtk_container_add(GTK_CONTAINER(content_area), label2);
	gtk_container_add(GTK_CONTAINER(content_area), label_delay_cam);
	gtk_container_add(GTK_CONTAINER(content_area), delay_cam);
	gtk_container_add(GTK_CONTAINER(content_area), find_delay);
	gtk_container_add(GTK_CONTAINER(content_area), sep);
	gtk_container_add(GTK_CONTAINER(content_area), label3);
	gtk_container_add(GTK_CONTAINER(content_area), apply_offset);
	gtk_widget_show_all(GTK_WIDGET(content_area));

}

// the public getter
GtkWidget *get_occult_dialog() {
	if (!dialog)
		build_the_dialog();
	
	return dialog;
}

static void on_find_clicked(GtkDialog* self, gint response_id, gpointer user_data){
//	struct phot_config *args = calloc(1, sizeof(struct phot_config));
	control_window_switch_to_tab(OUTPUT_LOGS);
	siril_log_message(_("Find button clicked \n"));		//c'est ici que la procedure de calcul doit etre lancée

	int layer = -1;
	if (com.seq.regparam) {
		for (int i = 0; i < com.seq.nb_layers; i++) {
			if (com.seq.regparam[i]) {
				layer = i;
				break;
			}
		}
	}
	if (layer == -1) {
		siril_debug_print("unregistered sequence\n");
		//siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The sequence must be registered"));
		if (com.seq.nb_layers == 3)
			layer = 1;
		else layer = 0;
	}

// Tedst if area selected
// Test if seq loaded

//	struct phot_config *args = calloc(1, sizeof(struct phot_config));
	struct seqpsf_args *spsfargs = malloc(sizeof(struct seqpsf_args));
	spsfargs->for_photometry = TRUE;

	struct light_curve_args *args = calloc(1, sizeof(struct light_curve_args));
	args->seq = &com.seq;
	args->layer = layer;
	start_in_new_thread(occultation_worker, args);


	gchar *end;
	const gchar *text = gtk_entry_get_text(GTK_ENTRY(delay_cam));
	delay = g_ascii_strtod(text, &end);
	if (text == end || delay <= -10.0 || delay > 30.7) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("DELAY range not accepted (should be ]0, 0.7]"));
		return;
	}
	siril_log_message(_("Delay value: %lf \n"), delay);
	delay += 1.0;

	gtk_entry_set_text(GTK_ENTRY(delay_cam), g_strdup_printf("%0.2lf", delay));
}


static void on_occult_response(GtkDialog* self, gint response_id, gpointer user_data) {
	siril_debug_print("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		gtk_widget_hide(dialog);
		return;
	}


//	if (!has_wcs(&gfit)) {
//		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The currently loaded image must be plate solved"));
//		gtk_widget_hide(dialog);
//		return;
//	}

//	const gchar *entered_target_name = gtk_entry_get_text(GTK_ENTRY(target_entry));
//	gchar *target_name = g_strdup(entered_target_name);
//	g_strstrip(target_name);
//	if (target_name[0] == '\0') {
//		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Enter the name of the target star"));
//		return;
//	}

/*	gchar *end;
	const gchar *text = gtk_entry_get_text(GTK_ENTRY(delta_vmag_entry));
	double delta_Vmag = g_ascii_strtod(text, &end);
	if (text == end || delta_Vmag <= 0.0 || delta_Vmag > 6.0) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Vmag range not accepted (should be ]0, 6])"));
		return;
	}
	text = gtk_entry_get_text(GTK_ENTRY(delta_bv_entry));
	double delta_BV = g_ascii_strtod(text, &end);
	if (text == end || delta_BV <= 0.0 || delta_BV > 0.7) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("BV range not accepted (should be ]0, 0.7]"));
		return;
	}
	text = gtk_entry_get_text(GTK_ENTRY(emag_entry));
	double emag = g_ascii_strtod(text, &end);
	if (text == end || emag <= 0.0 || emag > 0.1) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Magnitude error not accepted (should be ]0, 0.1["));
		return;
	}
*/
//	gboolean use_apass = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(apass_radio));
//	gboolean narrow = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_narrow));
	control_window_switch_to_tab(OUTPUT_LOGS);
	gboolean use_offset = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(apply_offset));
//	struct phot_config *args = calloc(1, sizeof(struct phot_config));
	struct light_curve_args *args = calloc(1, sizeof(struct light_curve_args));
//	args->seq->psf->occult_is_valid = TRUE;
	if (use_offset) {
		siril_log_message(_("Applied offset: %lf \n"), delay);
		args->time_offset = TRUE;
		args->JD_offset = delay;
	}
	else {
		siril_log_message(_("No offset applied \n"));
		args->time_offset = FALSE;
		args->JD_offset = 0.0;
	}
	gtk_widget_hide(dialog);
	

//	struct compstars_arg *args = calloc(1, sizeof(struct compstars_arg));
/*	args->fit = &gfit;
	args->target_name = g_strdup(target_name);
	args->narrow_fov = narrow;
	args->cat = use_apass ? CAT_APASS : CAT_NOMAD;
	args->delta_Vmag = delta_Vmag;
	args->delta_BV = delta_BV;
	args->max_emag = emag;
	args->nina_file = g_strdup("auto");
*/

//	start_in_new_thread(occultation_worker, args);
}
