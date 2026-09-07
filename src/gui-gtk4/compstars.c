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

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/comparison_stars.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/PSF_list.h"
#include "io/annotation_catalogues.h"
#include "io/siril_catalogues.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *delta_vmag_entry = NULL;
static GtkWidget *delta_bv_entry = NULL;
static GtkWidget *emag_entry = NULL;
static GtkWidget *target_entry = NULL;
static GtkWidget *manu_target_entry = NULL;
static GtkWidget *apass_radio = NULL;
static GtkWidget *band_combo = NULL;
static GtkWidget *check_narrow = NULL;
static GtkWidget *labelmag = NULL, *labelcolor = NULL;
static GtkWidget *auto_mode, *mode_grp, *manual_mode, *sub_manu_box;
static GtkWidget *auto_data_grp;

static void on_compstars_response(GtkWindow *self, gint response_id, gpointer user_data);

static GtkWidget *compstars_ok_button = NULL;

/* Phase 14G.4: GtkDialog → GtkWindow.  Bridge per-button "clicked" signals
 * back to the legacy on_compstars_response shape. */
static void compstars_btn_ok_clicked(GtkButton *btn, gpointer ud)    { (void)btn; on_compstars_response(GTK_WINDOW(dialog), GTK_RESPONSE_ACCEPT, ud); }
static void compstars_btn_close_clicked(GtkButton *btn, gpointer ud) { (void)btn; on_compstars_response(GTK_WINDOW(dialog), GTK_RESPONSE_REJECT, ud); }

/* Idle callback invoked on the GUI thread when compstars_worker finishes. */
static gboolean end_compstars(gpointer p) {
	siril_log_debug("end_compstars\n");
	struct compstars_arg *args = (struct compstars_arg *) p;

	clear_stars_list(args->has_GUI);
	if (args->has_GUI && !args->retval) {
		purge_user_catalogue(CAT_AN_USER_TEMP);
		if (!load_siril_cat_to_temp(args->comp_stars)) {
			/* the annotate button is bound to the win.annotate-object action;
			 * unlike GTK3, setting the toggle active in GTK4 does not activate
			 * the action, so drive the action state directly */
			GActionMap *map = G_ACTION_MAP(gtk_builder_get_object(gui.builder, "control_window"));
			GAction *annotate = g_action_map_lookup_action(map, "annotate-object");
			refresh_found_objects();
			GVariant *state = g_action_get_state(annotate);
			if (!g_variant_get_boolean(state)) {
				g_action_change_state(annotate, g_variant_new_boolean(TRUE));
			} else {
				refresh_found_objects();
				redraw(REDRAW_OVERLAY);
			}
			g_variant_unref(state);
		}
	} else {
		siril_catalog_free(args->comp_stars);
	}
	redraw(REDRAW_OVERLAY);
	free_compstars_arg(args);
	return end_generic(NULL);
}

static void output_state(GtkToggleButton *source, gpointer user_data) {
    gtk_widget_set_sensitive(auto_data_grp, siril_toggle_get_active(GTK_WIDGET(auto_mode)));
	gtk_widget_set_sensitive(sub_manu_box, siril_toggle_get_active(GTK_WIDGET(manual_mode)));
}

/* The band dropdown only lists what the selected catalogue can supply, and the
 * two range labels name the band and the colour index that go with it. */
static void fill_band_combo(siril_cat_index cat, phot_band selected) {
	int index = 0, to_select = 0;
	siril_drop_down_clear_strings(GTK_DROP_DOWN(band_combo));
	for (int band = 0; band < PHOT_NB_BANDS; band++) {
		if (!catalogue_has_band(cat, band))
			continue;
		gchar *label = g_strdup_printf("%s (%s)", phot_band_to_str(band), phot_band_description(band));
		siril_drop_down_append_text(GTK_DROP_DOWN(band_combo), label);
		g_free(label);
		if (band == selected)
			to_select = index;
		index++;
	}
	gtk_drop_down_set_selected(GTK_DROP_DOWN(band_combo), to_select);
}

// the dropdown only holds the bands of the current catalogue, map back to phot_band
static phot_band get_selected_band(siril_cat_index cat) {
	guint selected = gtk_drop_down_get_selected(GTK_DROP_DOWN(band_combo));
	guint index = 0;
	for (int band = 0; band < PHOT_NB_BANDS; band++) {
		if (!catalogue_has_band(cat, band))
			continue;
		if (index == selected)
			return (phot_band)band;
		index++;
	}
	return PHOT_BAND_V;
}

static void update_band_labels() {
	siril_cat_index cat = siril_toggle_get_active(GTK_WIDGET(apass_radio)) ? CAT_APASS : CAT_NOMAD;
	phot_band band = get_selected_band(cat);
	gchar *text = g_strdup_printf(_("Allowed %s magnitude range:"), phot_band_to_str(band));
	gtk_label_set_text(GTK_LABEL(labelmag), text);
	g_free(text);
	text = g_strdup_printf(_("Allowed %s index range:"), phot_band_color_to_str(band));
	gtk_label_set_text(GTK_LABEL(labelcolor), text);
	g_free(text);
}

static void on_band_changed(GObject *self, GParamSpec *pspec, gpointer user_data) {
	update_band_labels();
}

static void on_catalogue_changed(GtkToggleButton *source, gpointer user_data) {
	siril_cat_index cat = siril_toggle_get_active(GTK_WIDGET(apass_radio)) ? CAT_APASS : CAT_NOMAD;
	// keep the band across the change when the new catalogue also provides it
	phot_band band = get_selected_band(cat == CAT_APASS ? CAT_NOMAD : CAT_APASS);
	fill_band_combo(cat, catalogue_has_band(cat, band) ? band : PHOT_BAND_V);
	update_band_labels();
}

static void build_the_dialog() {
	/* Phase 14G.4: GtkDialog → GtkWindow. */
	dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Create a comparison stars list"));
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	/* A modal GtkWindow with no transient parent grabs input with no anchor
	 * and can wedge the whole session on Wayland/Xorg. */
	gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(lookup_widget("control_window")));
	gtk_window_set_hide_on_close(GTK_WINDOW(dialog), TRUE);
	g_signal_connect(G_OBJECT(dialog), "close-request", G_CALLBACK(siril_widget_hide_on_delete), NULL);


	/* Mode (Auto/Manu) choice */
	mode_grp = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(mode_grp), TRUE);
	gtk_widget_set_tooltip_text(mode_grp, _("Toggle Manual mode or Automatic mode for Comparison stars list"));

	/* Phase 17.6: GtkRadioButton removed in GTK4; use GtkCheckButton + group. */
	manual_mode = gtk_check_button_new_with_label(_("Use the stars selected in the currently loaded image"));
	g_signal_connect (manual_mode, "toggled",G_CALLBACK (output_state), NULL);
	gtk_box_append(GTK_BOX(mode_grp), manual_mode);

	// Name of the output file in manu mode
	// Definition of the 3 elements horizontal sub-box
	sub_manu_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 3);
	gtk_box_set_homogeneous(GTK_BOX(sub_manu_box), TRUE);
	gtk_widget_set_tooltip_text(sub_manu_box, _("Enter your own target name"));
	// 1st element, a label
	GtkWidget *label1_user_name = gtk_label_new(_("Output file name: "));
	gtk_widget_set_halign(label1_user_name, GTK_ALIGN_END);
	gtk_box_append(GTK_BOX(sub_manu_box), label1_user_name);
	// 2nd element, the target user name
	manu_target_entry = gtk_entry_new();
	gtk_editable_set_text(GTK_EDITABLE(manu_target_entry), "V_SirilstarList_user");
	gtk_widget_set_tooltip_text(manu_target_entry, _("Change the default file name if needed"));
	gtk_widget_set_halign(manu_target_entry, GTK_ALIGN_CENTER);
	gtk_entry_set_alignment(GTK_ENTRY (manu_target_entry), 0.5);
	gtk_widget_set_margin_top(GTK_WIDGET(manu_target_entry), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(manu_target_entry), 0);
	gtk_box_append(GTK_BOX(sub_manu_box), manu_target_entry);
	// 3rd element, another label
	GtkWidget *label2_user_name = gtk_label_new(_(".csv"));
	gtk_widget_set_halign(label2_user_name, GTK_ALIGN_START);
	gtk_box_append(GTK_BOX(sub_manu_box), label2_user_name);
	// and finally include that box to the upper level box
	gtk_box_append(GTK_BOX(mode_grp), sub_manu_box);

	auto_mode = gtk_check_button_new_with_label(_("Find comparison stars from catalogue requests"));
	gtk_check_button_set_group(GTK_CHECK_BUTTON(auto_mode), GTK_CHECK_BUTTON(manual_mode));
	g_signal_connect (auto_mode, "toggled",G_CALLBACK (output_state), NULL);
	gtk_box_append(GTK_BOX(mode_grp), auto_mode);

	gtk_widget_set_halign(mode_grp, GTK_ALIGN_START);
	gtk_widget_set_margin_start(GTK_WIDGET(mode_grp), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(mode_grp), 20);
	gtk_widget_set_margin_bottom(GTK_WIDGET(mode_grp), 20);
	// Defines the group for the auto mode parameters
	auto_data_grp = gtk_box_new(GTK_ORIENTATION_VERTICAL, 9);
	gtk_box_set_homogeneous(GTK_BOX(auto_data_grp), TRUE);
	gtk_widget_set_tooltip_text(auto_data_grp, _("Variable star data and sorting parameters for catalogue request"));

	// Defines the parameters of the automatic mode
	target_entry = gtk_entry_new();
	gtk_entry_set_placeholder_text(GTK_ENTRY(target_entry), "Target star name");
	gtk_widget_set_tooltip_text(target_entry, _("Enter the target star name to search in catalogues"));
	gtk_widget_set_margin_start(GTK_WIDGET(target_entry), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(target_entry), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(target_entry), 15);
	gtk_widget_set_margin_bottom(GTK_WIDGET(target_entry), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(target_entry), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(target_entry), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), target_entry);

	check_narrow = gtk_check_button_new_with_label(_("Narrow field of view"));
	gtk_widget_set_tooltip_text(check_narrow, _("Tick this box to use a narrow field of view centered about the target star"));
	gtk_widget_set_halign(check_narrow, GTK_ALIGN_START);
	gtk_widget_set_margin_start(GTK_WIDGET(check_narrow), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(check_narrow), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(check_narrow), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), check_narrow);

	GtkWidget *band_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_tooltip_text(band_box, _("Photometric band the comparison stars are selected in, "
				"named after the AAVSO filter designations. Choose the one matching the filter "
				"used for the observation"));
	GtkWidget *labelband = gtk_label_new(_("Photometric band:"));
	gtk_widget_set_halign(labelband, GTK_ALIGN_START);
	gtk_box_append(GTK_BOX(band_box), labelband);
	band_combo = gtk_drop_down_new(NULL, NULL);
	gtk_widget_set_hexpand(band_combo, TRUE);
	g_signal_connect(band_combo, "notify::selected", G_CALLBACK(on_band_changed), NULL);
	gtk_box_append(GTK_BOX(band_box), band_combo);
	gtk_widget_set_margin_start(GTK_WIDGET(band_box), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(band_box), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(band_box), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(band_box), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), band_box);

	labelmag = gtk_label_new(NULL);	// text set by update_band_labels()
	gtk_widget_set_halign(labelmag, GTK_ALIGN_START);
	gtk_widget_set_margin_start(GTK_WIDGET(labelmag), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(labelmag), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(labelmag), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), labelmag);

	delta_vmag_entry = gtk_entry_new();
	gtk_editable_set_text(GTK_EDITABLE(delta_vmag_entry), "3.0");
	gtk_widget_set_tooltip_text(delta_vmag_entry, _("Allowed range of magnitude, in the selected band, between the target star and the comparison stars"));
	gtk_widget_set_margin_start(GTK_WIDGET(delta_vmag_entry), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(delta_vmag_entry), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(delta_vmag_entry), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(delta_vmag_entry), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), delta_vmag_entry);

	labelcolor = gtk_label_new(NULL);	// text set by update_band_labels()
	gtk_widget_set_halign(labelcolor, GTK_ALIGN_START);
	gtk_widget_set_margin_start(GTK_WIDGET(labelcolor), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(labelcolor), 10);
	gtk_widget_set_margin_bottom(GTK_WIDGET(labelcolor), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), labelcolor);

	delta_bv_entry = gtk_entry_new();
	gtk_editable_set_text(GTK_EDITABLE(delta_bv_entry), "0.5");
	gtk_widget_set_tooltip_text(delta_bv_entry, _("Allowed range of color index between the target star and the comparison stars. "
				"The index is the one of the standard transformation equation for the selected band"));
	gtk_widget_set_margin_start(GTK_WIDGET(delta_bv_entry), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(delta_bv_entry), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(delta_bv_entry), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(delta_bv_entry), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), delta_bv_entry);

	GtkWidget *labelemag = gtk_label_new(_("Allowed magnitude error:"));
	gtk_widget_set_halign(labelemag, GTK_ALIGN_START);
	gtk_widget_set_margin_start(GTK_WIDGET(labelemag), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(labelemag), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(labelemag), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), labelemag);

	emag_entry = gtk_entry_new();
	gtk_editable_set_text(GTK_EDITABLE(emag_entry), "0.03");
	gtk_widget_set_tooltip_text(emag_entry, _("Allowed catalogue magnitude error for comparison stars, for catalogues supplying it"));
	gtk_widget_set_margin_start(GTK_WIDGET(emag_entry), 15);
	gtk_widget_set_margin_end(GTK_WIDGET(emag_entry), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(emag_entry), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(emag_entry), 0);
	gtk_box_append(GTK_BOX(auto_data_grp), emag_entry);

	/* catalogue choice */
	GtkWidget *nomad_radio, *cat_choice_box;
	cat_choice_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(cat_choice_box), TRUE);
	gtk_widget_set_tooltip_text(cat_choice_box, _("Recommended catalogue for this feature is APASS"));

	apass_radio = gtk_check_button_new_with_label(_("APASS catalogue"));
	nomad_radio = gtk_check_button_new_with_label(_("NOMAD catalogue"));
	gtk_check_button_set_group(GTK_CHECK_BUTTON(nomad_radio), GTK_CHECK_BUTTON(apass_radio));
	gtk_check_button_set_active(GTK_CHECK_BUTTON(apass_radio), TRUE);
	g_signal_connect(apass_radio, "toggled", G_CALLBACK(on_catalogue_changed), NULL);
	gtk_box_append(GTK_BOX(cat_choice_box), apass_radio);
	gtk_box_append(GTK_BOX(cat_choice_box), nomad_radio);
	gtk_box_append(GTK_BOX(auto_data_grp), cat_choice_box);
	gtk_widget_set_margin_start(GTK_WIDGET(cat_choice_box), 15);
	gtk_widget_set_margin_top(GTK_WIDGET(cat_choice_box), 0);
	gtk_widget_set_margin_bottom(GTK_WIDGET(cat_choice_box), 0);
	// Gather the graphic items
	GtkWidget *content_area = gtk_box_new(GTK_ORIENTATION_VERTICAL, 15);
	gtk_widget_set_margin_start(content_area, 12);
	gtk_widget_set_margin_end(content_area, 12);
	gtk_widget_set_margin_top(content_area, 12);
	gtk_widget_set_margin_bottom(content_area, 12);
	gtk_box_append(GTK_BOX(content_area), mode_grp);
	gtk_box_append(GTK_BOX(content_area), auto_data_grp);
	gtk_widget_set_sensitive (auto_data_grp, FALSE);

	/* Action area */
	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *btn_close = gtk_button_new_with_mnemonic(_("_Close"));
	g_signal_connect(btn_close, "clicked", G_CALLBACK(compstars_btn_close_clicked), NULL);
	gtk_box_append(GTK_BOX(bbox), btn_close);
	compstars_ok_button = gtk_button_new_with_mnemonic(_("_OK"));
	gtk_widget_add_css_class(compstars_ok_button, "suggested-action");
	g_signal_connect(compstars_ok_button, "clicked", G_CALLBACK(compstars_btn_ok_clicked), NULL);
	gtk_box_append(GTK_BOX(bbox), compstars_ok_button);
	gtk_box_append(GTK_BOX(content_area), bbox);
	gtk_window_set_default_widget(GTK_WINDOW(dialog), compstars_ok_button);

	gtk_window_set_child(GTK_WINDOW(dialog), content_area);

	fill_band_combo(CAT_APASS, PHOT_BAND_V);
	update_band_labels();
}

// The process to perform a **Manual** Compstar List
static void manual_photometry_data (sequence *seq) {
	gchar *entered_target_name = g_strdup(gtk_editable_get_text(GTK_EDITABLE(manu_target_entry)));
	if (entered_target_name [0] == '\0') {
		g_free(entered_target_name);
		entered_target_name = g_strdup("V_SirilstarList_user");
		gtk_editable_set_text(GTK_EDITABLE(manu_target_entry), "V_SirilstarList_user");
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
		siril_log_warning(_("One Variable star and one comparison star at least are required. Cannot create any file\n"));
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("One Variable star and one comparison star at least are required. Cannot create any file"));
		g_free(entered_target_name);
		return;
	}
	point sel_item[MAX_SEQPSF];

	for (int r = 0; r < MAX_SEQPSF && seq->photometry[r]; r++) {
		if (get_ra_and_dec_from_star_pos(seq->photometry[r][seq->current], &ra, &dec)) {
			siril_log_error(_("Problem with conversion\n")); // PB in the conversion pix->wcs
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
	args->delta_mag = 0.0;		// Explicitely set these three variables
	args->delta_color = 0.0;
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
	const gchar *entered_target_name = gtk_editable_get_text(GTK_EDITABLE(target_entry));
	gchar *target_name = g_strdup(entered_target_name);
	g_strstrip(target_name);
	if (target_name[0] == '\0') {
		g_free(target_name);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Enter the name of the target star"));
		return;
	}

	gchar *end;
	const gchar *text = gtk_editable_get_text(GTK_EDITABLE(delta_vmag_entry));
	double delta_Vmag = g_ascii_strtod(text, &end);
	if (text == end || delta_Vmag <= 0.0 || delta_Vmag > 6.0) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Vmag range not accepted (should be ]0, 6])"));
		g_free(target_name);
		return;
	}
	text = gtk_editable_get_text(GTK_EDITABLE(delta_bv_entry));
	double delta_BV = g_ascii_strtod(text, &end);
	if (text == end || delta_BV <= 0.0 || delta_BV > 0.7) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("BV range not accepted (should be ]0, 0.7]"));
		g_free(target_name);
		return;
	}
	text = gtk_editable_get_text(GTK_EDITABLE(emag_entry));
	double emag = g_ascii_strtod(text, &end);
	if (text == end || emag <= 0.0 || emag > 0.1) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Magnitude error not accepted (should be ]0, 0.1["));
		g_free(target_name);
		return;
	}

	gboolean use_apass = siril_toggle_get_active(GTK_WIDGET(apass_radio));
	siril_cat_index cat = use_apass ? CAT_APASS : CAT_NOMAD;
	gboolean narrow = siril_toggle_get_active(GTK_WIDGET(check_narrow));
	control_window_switch_to_tab(OUTPUT_LOGS);

	struct compstars_arg *args = calloc(1, sizeof(struct compstars_arg));
	args->fit = gfit;
	args->target_name = g_strdup(target_name);
	g_free(target_name);
	args->narrow_fov = narrow;
	args->cat = cat;
	args->band = get_selected_band(cat);
	args->delta_mag = delta_Vmag;
	args->delta_color = delta_BV;
	args->max_emag = emag;
	args->nina_file = g_strdup("auto");
	args->notify_done = end_compstars;

	if(!start_in_new_thread(compstars_worker, args)) {
		g_free(args->target_name);
		g_free(args->nina_file);
		free(args);
	}
}

// the public getter
GtkWidget *get_compstars_dialog() {
	if (!dialog)
		build_the_dialog();
	return dialog;
}

static void on_compstars_response(GtkWindow *self, gint response_id, gpointer user_data) {
	(void)self; (void)user_data;
	siril_log_debug("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		if (compstars_ok_button)
			gtk_widget_grab_focus(compstars_ok_button);
		gtk_widget_set_visible(dialog, FALSE);
		reactivate_parent(dialog);
		return;
	}

	if (siril_toggle_get_active(GTK_WIDGET(manual_mode)))
		manual_photometry_data(&com.seq);

	if (siril_toggle_get_active(GTK_WIDGET(auto_mode)))
		auto_photometry_data();

}
