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
#include "gui/image_display.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *delay_cam = NULL;
static GtkWidget *sep = NULL;
static GtkWidget *apply_offset = NULL;
static double delay = 0.0;

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
	if (com.seq.photometry[0] != NULL) free_photometry_set(&com.seq, 0);
	control_window_switch_to_tab(OUTPUT_LOGS);

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
	if (com.selection.w == 0 || com.selection.h  == 0) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Make a selection aroud the blinking star first"));
		return;
	}

	if (com.selection.w < com.pref.phot_set.outer || com.selection.h < com.pref.phot_set.outer) {
		siril_log_color_message(_("The selection has benen resized \n"), "salmon");
		com.selection.w = 2.0 * com.pref.phot_set.outer;
		com.selection.h = 2.0 * com.pref.phot_set.outer;
	}

	struct light_curve_args *args = calloc(1, sizeof(struct light_curve_args));
	args->seq = &com.seq;
	start_in_new_thread(occultation_worker, args);

	delay = 0.0;
	gtk_entry_set_text(GTK_ENTRY(delay_cam), g_strdup_printf("%0.3lf", delay));
	control_window_switch_to_tab(OUTPUT_LOGS);
}

gboolean end_occultation_worker(gpointer p) {
	if (!com.script) {
		struct light_curve_args *args = (struct light_curve_args *)p;
		args->seq = &com.seq;
		delay = args->JD_offset;
		gtk_entry_set_text(GTK_ENTRY(delay_cam), g_strdup_printf("%0.3lf", delay));
		control_window_switch_to_tab(OUTPUT_LOGS);
		free_light_curve_args(args);
	}
	if (sequence_is_loaded()) {
		drawPlot();
		notify_new_photometry();	// switch to and update plot tab
		redraw(REDRAW_OVERLAY);
	}
	return end_generic(NULL);
}

static void on_occult_response(GtkDialog* self, gint response_id, gpointer user_data) {
	siril_debug_print("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		if (com.seq.photometry[0] != NULL) free_photometry_set(&com.seq, 0);
		gtk_widget_hide(dialog);
		return;
	}

	control_window_switch_to_tab(OUTPUT_LOGS);
	gboolean use_offset = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(apply_offset));
	struct light_curve_args *args = calloc(1, sizeof(struct light_curve_args));
	args->seq = &com.seq;
	if (use_offset) {
		args->time_offset = TRUE;
		args->JD_offset = delay;
		siril_log_message(_("Applied offset: %0.3lf (ms) \n"), args->JD_offset);
	}
	else {
		siril_log_message(_("No offset applied \n"));
		args->time_offset = FALSE;
		args->JD_offset = 0.0;
	}

	if (com.seq.photometry[0] != NULL) free_photometry_set(&com.seq, 0);
	free_light_curve_args(args);

	gtk_widget_hide(dialog);


}
