/*
 * This file is part of Siril, an astronomy image processor.
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

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/comparison_stars.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/PSF_list.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *delta_vmag_entry = NULL;
static GtkWidget *delta_bv_entry = NULL;
static GtkWidget *emag_entry = NULL;
static GtkWidget *target_entry = NULL;
static GtkWidget *manu_target_entry = NULL;
static GtkWidget *apass_radio = NULL;
static GtkWidget *check_narrow = NULL;
static GtkWidget *auto_mode, *mode_grp, *manual_mode, *sub_manu_box;
static GtkWidget *auto_data_grp;

static void on_compstars_response(GtkDialog* self, gint response_id, gpointer user_data);

static void output_state(GtkToggleButton *source, gpointer user_data) {
    gtk_widget_set_sensitive(auto_data_grp, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(auto_mode)));
	gtk_widget_set_sensitive(sub_manu_box, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(manual_mode)));
}

static void build_the_dialog() {
	dialog = gtk_dialog_new_with_buttons(_("Create a comparison stars list"), NULL,
			0, _("_Close"), GTK_RESPONSE_REJECT, _("_OK"), GTK_RESPONSE_ACCEPT, NULL);
	// If the user clicks one of these dialog buttons, GtkDialog will emit
	// the GtkDialog::response signal with the corresponding response ID
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
	g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(siril_widget_hide_on_delete), NULL);
	g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_compstars_response), NULL);

	// Sets the "OK" button as default one
	GtkWidget *OK_button = gtk_dialog_get_widget_for_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
	gtk_widget_set_can_default(OK_button, TRUE);
	gtk_widget_grab_default(OK_button);
	gtk_widget_grab_focus(OK_button);
	gtk_style_context_add_class(gtk_widget_get_style_context(OK_button), "suggested-action");


	/* Mode (Auto/Manu) choice */
	mode_grp = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(mode_grp), TRUE);
	gtk_widget_set_tooltip_text(mode_grp, _("Toggle Manual mode or Automatic mode for Comparison stars list"));

	manual_mode = gtk_radio_button_new_with_label_from_widget(NULL, _("Use the stars selected in the currently loaded image"));
	g_signal_connect (manual_mode, "toggled",G_CALLBACK (output_state), NULL);
	gtk_container_add(GTK_CONTAINER(mode_grp), manual_mode);

	// Name of the output file in manu mode
	// Definition of the 3 elements horizontal sub-box
	sub_manu_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 3);
	gtk_box_set_homogeneous(GTK_BOX(sub_manu_box), TRUE);
	gtk_widget_set_tooltip_text(sub_manu_box, _("Enter your own target name"));
	// 1st element, a label
	GtkWidget *label1_user_name = gtk_label_new(_("Output file name: "));
	gtk_widget_set_halign(label1_user_name, GTK_ALIGN_END);
	gtk_container_add(GTK_CONTAINER(sub_manu_box), label1_user_name);
	// 2nd element, the target user name
	manu_target_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(manu_target_entry), "V_SirilstarList_user");
	gtk_widget_set_tooltip_text(manu_target_entry, _("Change the default file name if needed"));
	gtk_widget_set_halign(manu_target_entry, GTK_ALIGN_CENTER);
	gtk_entry_set_alignment(GTK_ENTRY (manu_target_entry), 0.5);
	g_object_set(G_OBJECT(manu_target_entry), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(manu_target_entry), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(sub_manu_box), manu_target_entry);
	// 3rd element, another label
	GtkWidget *label2_user_name = gtk_label_new(_(".csv"));
	gtk_widget_set_halign(label2_user_name, GTK_ALIGN_START);
	gtk_container_add(GTK_CONTAINER(sub_manu_box), label2_user_name);
	// and finally include that box to the upper level box
	gtk_container_add(GTK_CONTAINER(mode_grp), sub_manu_box);

	auto_mode = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(manual_mode), _("Find comparison stars from catalogue requests"));
	g_signal_connect (auto_mode, "toggled",G_CALLBACK (output_state), NULL);
	gtk_container_add(GTK_CONTAINER(mode_grp), auto_mode);

	gtk_widget_set_halign(mode_grp, GTK_ALIGN_START);
	g_object_set(G_OBJECT(mode_grp), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(mode_grp), "margin-top", 20, NULL);
	g_object_set(G_OBJECT(mode_grp), "margin-bottom", 20, NULL);

	// Defines the group for the auto mode parameters
	auto_data_grp = gtk_box_new(GTK_ORIENTATION_VERTICAL, 9);
	gtk_box_set_homogeneous(GTK_BOX(auto_data_grp), TRUE);
	gtk_widget_set_tooltip_text(auto_data_grp, _("Variable star data and sorting parameters for catalogue request"));

	// Defines the parameters of the automatic mode
	target_entry = gtk_entry_new();
	gtk_entry_set_placeholder_text(GTK_ENTRY(target_entry), "Target star name");
	gtk_widget_set_tooltip_text(target_entry, _("Enter the target star name to search in catalogues"));
	g_object_set(G_OBJECT(target_entry), "margin", 15, NULL);
	g_object_set(G_OBJECT(target_entry), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(target_entry), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), target_entry);

	check_narrow = gtk_check_button_new_with_label(_("Narrow field of view"));
	gtk_widget_set_tooltip_text(check_narrow, _("Tick this box to use a narrow field of view centered about the target star"));
	gtk_widget_set_halign(check_narrow, GTK_ALIGN_START);
	g_object_set(G_OBJECT(check_narrow), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(check_narrow), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(check_narrow), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), check_narrow);

	GtkWidget *labelvmag = gtk_label_new(_("Allowed visual magnitude range:"));
	gtk_widget_set_halign(labelvmag, GTK_ALIGN_START);
	g_object_set(G_OBJECT(labelvmag), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(labelvmag), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(labelvmag), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), labelvmag);

	delta_vmag_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(delta_vmag_entry), "3.0");
	gtk_widget_set_tooltip_text(delta_vmag_entry, _("Allowed range of visual magnitude between the target star and the comparison stars"));
	g_object_set(G_OBJECT(delta_vmag_entry), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(delta_vmag_entry), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(delta_vmag_entry), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(delta_vmag_entry), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), delta_vmag_entry);

	GtkWidget *labelbv = gtk_label_new(_("Allowed B-V index range:"));
	gtk_widget_set_halign(labelbv, GTK_ALIGN_START);
	g_object_set(G_OBJECT(labelbv), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(labelbv), "margin-top", 10, NULL);
	g_object_set(G_OBJECT(labelbv), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), labelbv);

	delta_bv_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(delta_bv_entry), "0.5");
	gtk_widget_set_tooltip_text(delta_bv_entry, _("Allowed range of B-V index (color) between the target star and the comparison stars"));
	g_object_set(G_OBJECT(delta_bv_entry), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(delta_bv_entry), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(delta_bv_entry), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(delta_bv_entry), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), delta_bv_entry);

	GtkWidget *labelemag = gtk_label_new(_("Allowed magnitude error:"));
	gtk_widget_set_halign(labelemag, GTK_ALIGN_START);
	g_object_set(G_OBJECT(labelemag), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(labelemag), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(labelemag), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), labelemag);

	emag_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(emag_entry), "0.03");
	gtk_widget_set_tooltip_text(emag_entry, _("Allowed catalogue magnitude error for comparison stars"));
	g_object_set(G_OBJECT(emag_entry), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(emag_entry), "margin-right", 15, NULL);
	g_object_set(G_OBJECT(emag_entry), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(emag_entry), "margin-bottom", 0, NULL);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), emag_entry);

	/* catalogue choice */
	GtkWidget *nomad_radio, *cat_choice_box;
	cat_choice_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(cat_choice_box), TRUE);
	gtk_widget_set_tooltip_text(cat_choice_box, _("Recommended catalogue for this feature is APASS"));

	apass_radio = gtk_radio_button_new_with_label_from_widget(NULL, _("APASS catalogue"));
	nomad_radio = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(apass_radio), _("NOMAD catalogue"));
	gtk_container_add(GTK_CONTAINER(cat_choice_box), apass_radio);
	gtk_container_add(GTK_CONTAINER(cat_choice_box), nomad_radio);
	gtk_container_add(GTK_CONTAINER(auto_data_grp), cat_choice_box);
	g_object_set(G_OBJECT(cat_choice_box), "margin-left", 15, NULL);
	g_object_set(G_OBJECT(cat_choice_box), "margin-top", 0, NULL);
	g_object_set(G_OBJECT(cat_choice_box), "margin-bottom", 0, NULL);

	// Gather the graphic items
	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_box_set_spacing(GTK_BOX(content_area), 15);
	gtk_container_add(GTK_CONTAINER(content_area), mode_grp);
	gtk_container_add(GTK_CONTAINER(content_area), auto_data_grp);
	gtk_widget_set_sensitive (auto_data_grp, FALSE);
	gtk_widget_show_all(GTK_WIDGET(content_area));
}

// The process to perform a **Manual** Compstar List
static void manual_photometry_data (sequence *seq) {
	gchar *entered_target_name = g_strdup(gtk_entry_get_text(GTK_ENTRY(manu_target_entry)));
	if (entered_target_name [0] == '\0') {
		g_free(entered_target_name);
		entered_target_name = g_strdup("V_SirilstarList_user");
		gtk_entry_set_text(GTK_ENTRY(manu_target_entry), "V_SirilstarList_user");
	}

	gchar *temp_name = g_strdup(entered_target_name);
	g_strstrip(temp_name);
	gchar *target_name = g_strdup_printf("%s.csv", temp_name);
	g_free(temp_name);

	double ra, dec;
	// Gather the selected stars by hand
	int nb_ref_stars = 0;
	if (!seq->photometry[0] || !seq->photometry[1]) {
		g_free(target_name);
		siril_log_color_message(_("One Variable star and one comparison star at least are required. Cannot create any file\n"), "salmon");
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("One Variable star and one comparison star at least are required. Cannot create any file"));
		g_free(entered_target_name);
		return;
	}
	point sel_item[MAX_SEQPSF];

	for (int r = 0; r < MAX_SEQPSF && seq->photometry[r]; r++) {
		if (get_ra_and_dec_from_star_pos(seq->photometry[r][seq->current], &ra, &dec)) {
			siril_log_color_message(_("Problem with conversion\n"), "red"); // PB in the conversion pix->wcs
			g_free(entered_target_name);
			g_free(target_name);
			return;
		}
		sel_item[r].x = ra;
		sel_item[r].y = dec;
		nb_ref_stars++;
	}

	control_window_switch_to_tab(OUTPUT_LOGS);

	siril_catalogue *comp_sta = siril_catalog_new(CAT_COMPSTARS);

	// Header for the console display
	siril_log_message(_("-> %i comparison stars selected\n"), nb_ref_stars - 1);
	siril_log_message("Star type        RA      DEC\n");
	// Allocating final sorted list to the required size
	cat_item *result = calloc(nb_ref_stars, sizeof(cat_item));
	// Write the target star
	fill_compstar_item(&result[0], sel_item[0].x, sel_item[0].y, 0.0, "V", "Target");
	siril_log_message(_("Target star  : %4.3lf, %+4.3lf\n"), sel_item[0].x, sel_item[0].y);
	// And write the selected comparison stars
	for (int i = 1; i < nb_ref_stars; i++) {
		gchar *name = g_strdup_printf("%d", i);
		fill_compstar_item(&result[i], sel_item[i].x, sel_item[i].y, 0.0, name, "Comp1");
		g_free(name);
		siril_log_message(_("Comp star %3d: %4.3lf, %+4.3lf\n"), i, sel_item[i].x, sel_item[i].y);
	}

	// Fill the catalogue structure
	comp_sta->cat_items = result;
	comp_sta->nbitems = nb_ref_stars;
	comp_sta->nbincluded = nb_ref_stars;
	// Fill the other catalogue  structure
	struct compstars_arg *args = calloc(1, sizeof(struct compstars_arg));

	args->comp_stars = comp_sta;
	args->nina_file = g_strdup(target_name);
	args->target_star = &result[0];
	args->delta_Vmag = 0.0;		// Explicitely set these three variables
	args->delta_BV = 0.0;
	args->max_emag = 0.0;
	args->cat = CAT_COMPSTARS;
	// Finally create the csv file
	write_nina_file(args);
	// and free the stuff
	siril_catalog_free(comp_sta);
	g_free(args->nina_file);
	g_free(target_name);
	g_free(entered_target_name);
	free(args);
}

// The process to perform an **Automatic** Compstar List
static void auto_photometry_data () {
	const gchar *entered_target_name = gtk_entry_get_text(GTK_ENTRY(target_entry));
	gchar *target_name = g_strdup(entered_target_name);
	g_strstrip(target_name);
	if (target_name[0] == '\0') {
		g_free(target_name);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Enter the name of the target star"));
		return;
	}

	gchar *end;
	const gchar *text = gtk_entry_get_text(GTK_ENTRY(delta_vmag_entry));
	double delta_Vmag = g_ascii_strtod(text, &end);
	if (text == end || delta_Vmag <= 0.0 || delta_Vmag > 6.0) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Vmag range not accepted (should be ]0, 6])"));
		g_free(target_name);
		return;
	}
	text = gtk_entry_get_text(GTK_ENTRY(delta_bv_entry));
	double delta_BV = g_ascii_strtod(text, &end);
	if (text == end || delta_BV <= 0.0 || delta_BV > 0.7) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("BV range not accepted (should be ]0, 0.7]"));
		g_free(target_name);
		return;
	}
	text = gtk_entry_get_text(GTK_ENTRY(emag_entry));
	double emag = g_ascii_strtod(text, &end);
	if (text == end || emag <= 0.0 || emag > 0.1) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Magnitude error not accepted (should be ]0, 0.1["));
		g_free(target_name);
		return;
	}

	gboolean use_apass = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(apass_radio));
	gboolean narrow = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_narrow));
	control_window_switch_to_tab(OUTPUT_LOGS);

	struct compstars_arg *args = calloc(1, sizeof(struct compstars_arg));
	args->fit = &gfit;
	args->target_name = g_strdup(target_name);
	g_free(target_name);
	args->narrow_fov = narrow;
	args->cat = use_apass ? CAT_APASS : CAT_NOMAD;
	args->delta_Vmag = delta_Vmag;
	args->delta_BV = delta_BV;
	args->max_emag = emag;
	args->nina_file = g_strdup("auto");

	if(!start_in_new_thread(compstars_worker, args)) {
		g_free(args->target_name);
		free(args);
	}
}

// the public getter
GtkWidget *get_compstars_dialog() {
	if (!dialog)
		build_the_dialog();
	return dialog;
}

static void on_compstars_response(GtkDialog* self, gint response_id, gpointer user_data) {
//	auto_manu =gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(auto_mode))
	siril_debug_print("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		gtk_widget_grab_focus(gtk_dialog_get_widget_for_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT));
		gtk_widget_hide(dialog);
		return;
	}

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(manual_mode)))
		manual_photometry_data(&com.seq);

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(auto_mode)))
		auto_photometry_data();

}
