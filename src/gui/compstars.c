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

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/siril_wcs.h"
#include "algos/comparison_stars.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *delta_vmag_entry = NULL;
static GtkWidget *delta_bv_entry = NULL;
static GtkWidget *target_entry = NULL;
static GtkWidget *apass_radio = NULL;

static void on_compstars_response(GtkDialog* self, gint response_id, gpointer user_data);

static void build_the_dialog() {
	dialog = gtk_dialog_new_with_buttons(_("Create a comparison stars list"), NULL,
			0, _("_OK"), GTK_RESPONSE_ACCEPT, _("_Close"), GTK_RESPONSE_REJECT, NULL);
	// If the user clicks one of these dialog buttons, GtkDialog will emit
	// the GtkDialog::response signal with the corresponding response ID
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
	g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_compstars_response), NULL);

	GtkWidget *label = gtk_label_new(_("Find comparison stars for the currently loaded image and the given star name"));
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	g_object_set(G_OBJECT(label), "margin", 15, NULL);

	target_entry = gtk_entry_new();
	gtk_entry_set_placeholder_text(GTK_ENTRY(target_entry), "Target star name");
	gtk_widget_set_tooltip_text(target_entry, _("Enter the target star name to search in cataloges"));
	g_object_set(G_OBJECT(target_entry), "margin", 15, NULL);

	GtkWidget *labelvmag = gtk_label_new(_("Allowed visual magnitude range:"));
	gtk_widget_set_halign(labelvmag, GTK_ALIGN_START);
	g_object_set(G_OBJECT(labelvmag), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(labelvmag), "margin-top", 10, NULL);
	g_object_set(G_OBJECT(labelvmag), "margin-bottom", 0, NULL);

	delta_vmag_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(delta_vmag_entry), "3.0");
	gtk_widget_set_tooltip_text(delta_vmag_entry, _("Allowed range of visual magnitude between the target star and the comparison stars"));
	g_object_set(G_OBJECT(delta_vmag_entry), "margin-left", 15, NULL);
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
	g_object_set(G_OBJECT(delta_bv_entry), "margin-top", 0, NULL);

	/* catalogue choice */
	GtkWidget *radio2, *radiobox;
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
	gtk_container_add(GTK_CONTAINER(content_area), labelvmag);
	gtk_container_add(GTK_CONTAINER(content_area), delta_vmag_entry);
	gtk_container_add(GTK_CONTAINER(content_area), labelbv);
	gtk_container_add(GTK_CONTAINER(content_area), delta_bv_entry);
	gtk_container_add(GTK_CONTAINER(content_area), radiobox);
	gtk_widget_show_all(GTK_WIDGET(content_area));
}

// the public getter
GtkWidget *get_compstars_dialog() {
	if (!dialog)
		build_the_dialog();
	return dialog;
}

static void on_compstars_response(GtkDialog* self, gint response_id, gpointer user_data) {
	siril_debug_print("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		gtk_widget_hide(dialog);
		return;
	}
	if (!has_wcs(&gfit)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The currently loaded image must be plate solved"));
		gtk_widget_hide(dialog);
		return;
	}

	const gchar *entered_target_name = gtk_entry_get_text(GTK_ENTRY(target_entry));
	gchar *target_name = g_strdup(entered_target_name);
	g_strstrip(target_name);
	if (target_name[0] == '\0') {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Enter the name of the target star"));
		return;
	}

	gchar *end;
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

	gboolean use_apass = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(apass_radio));
	control_window_switch_to_tab(OUTPUT_LOGS);

	struct compstars_arg *args = calloc(1, sizeof(struct compstars_arg));
	args->target_name = g_strdup(target_name);
	args->narrow_fov = TRUE;
	args->cat = use_apass ? CAT_APASS : CAT_NOMAD;
	args->delta_Vmag = delta_Vmag;
	args->delta_BV = delta_BV;
	args->nina_file = g_strdup("auto");

	start_in_new_thread(compstars_worker, args);
}
