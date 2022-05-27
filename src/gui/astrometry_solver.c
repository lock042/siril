/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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

#include <glib.h>
#include <gtk/gtk.h>
#include "algos/astrometry_solver.h"
#include "core/processing.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/PSF_list.h"
#include "gui/photometric_cc.h"
#include "io/single_image.h"
#include "io/sequence.h"

enum {
	COLUMN_RESOLVER,// string
	COLUMN_NAME,	// string
	N_COLUMNS
};

static GtkListStore *list_IPS = NULL;
extern struct sky_object platedObject[RESOLVER_NUMBER];

static void unselect_all_items();
void on_GtkTreeViewIPS_cursor_changed(GtkTreeView *tree_view, gpointer user_data);

static void initialize_ips_dialog() {
	GtkWidget *button_ips_ok, *button_cc_ok, *catalog_label, *catalog_box_ips,
			*catalog_box_pcc, *catalog_auto, *frame_cc_bkg, *frame_cc_norm,
			*catalog_label_pcc;
	GtkWindow *parent;

	button_ips_ok = lookup_widget("buttonIPS_ok");
	button_cc_ok = lookup_widget("button_cc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_ips = lookup_widget("ComboBoxIPSCatalog");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_auto = lookup_widget("GtkCheckButton_OnlineCat");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	frame_cc_norm = lookup_widget("frame_cc_norm");

	parent = GTK_WINDOW(lookup_widget("ImagePlateSolver_Dial"));

	gtk_widget_set_visible(button_ips_ok, TRUE);
	gtk_widget_set_visible(button_cc_ok, FALSE);
	gtk_widget_set_visible(catalog_label, TRUE);
	gtk_widget_set_visible(catalog_label_pcc, FALSE);
	gtk_widget_set_visible(catalog_box_ips, TRUE);
	gtk_widget_set_visible(catalog_box_pcc, FALSE);
	gtk_widget_set_visible(catalog_auto, TRUE);
	gtk_widget_set_visible(frame_cc_bkg, FALSE);
	gtk_widget_set_visible(frame_cc_norm, FALSE);

	gtk_window_set_title(parent, _("Image Plate Solver"));
}

static void get_mag_from_GUI(struct astrometry_data *args) {
	GtkToggleButton *autobutton = GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButton_Mag_Limit"));
	args->auto_magnitude = gtk_toggle_button_get_active(autobutton);
	if (!args->auto_magnitude) {
		GtkSpinButton *magButton = GTK_SPIN_BUTTON(
				lookup_widget("GtkSpinIPS_Mag_Limit"));
		args->forced_magnitude = gtk_spin_button_get_value(magButton);
	}
}

static double get_focal() {
	GtkEntry *focal_entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	const gchar *value = gtk_entry_get_text(focal_entry);

	return g_ascii_strtod(value, NULL);
}

/* get pixel in µm */
static double get_pixel() {
	GtkEntry *pixel_entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));
	const gchar *value = gtk_entry_get_text(pixel_entry);

	return g_ascii_strtod(value, NULL);
}

static int get_server_from_combobox() {
	GtkComboBoxText *box = GTK_COMBO_BOX_TEXT(lookup_widget("combo_server_ips"));
	return gtk_combo_box_get_active(GTK_COMBO_BOX(box));
}

static online_catalog get_online_catalog(double fov, double mag) {
	GtkToggleButton *auto_button;
	int ret;

	auto_button = GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButton_OnlineCat"));
	if (gtk_toggle_button_get_active(auto_button)) {
		if (mag <= 6.5) {
			ret = BRIGHT_STARS;
		} else if (fov > 180.0) {
			ret = NOMAD;
		} else if (fov < 30.0){
			ret = GAIAEDR3;
		} else {
			ret = PPMXL;
		}
		return ret;
	} else {
		GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("ComboBoxIPSCatalog"));
		ret = gtk_combo_box_get_active(box);
		return (ret < 0 ? NOMAD : ret);
	}
}

static SirilWorldCS *get_center_of_catalog() {
	GtkSpinButton *GtkSpinIPS_RA_h, *GtkSpinIPS_RA_m;
	GtkSpinButton *GtkSpinIPS_Dec_deg, *GtkSpinIPS_Dec_m;
	GtkEntry *GtkEntryIPS_RA_s, *GtkEntryIPS_Dec_s;
	GtkToggleButton *GtkCheckButtonIPS_S;

	/* get alpha center */
	GtkSpinIPS_RA_h = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_h"));
	GtkSpinIPS_RA_m = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_m"));
	GtkEntryIPS_RA_s = GTK_ENTRY(lookup_widget("GtkEntryIPS_RA_s"));

	gdouble hour = gtk_spin_button_get_value_as_int(GtkSpinIPS_RA_h);
	gdouble min = gtk_spin_button_get_value_as_int(GtkSpinIPS_RA_m);
	gdouble sec = g_ascii_strtod(gtk_entry_get_text(GtkEntryIPS_RA_s), NULL);

	/* get Dec center */
	GtkSpinIPS_Dec_deg = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_deg"));
	GtkSpinIPS_Dec_m = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_m"));
	GtkEntryIPS_Dec_s = GTK_ENTRY(lookup_widget("GtkEntryIPS_Dec_s"));
	GtkCheckButtonIPS_S = GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButtonIPS_S"));

	gdouble deg = gtk_spin_button_get_value_as_int(GtkSpinIPS_Dec_deg);
	gdouble m = gtk_spin_button_get_value_as_int(GtkSpinIPS_Dec_m);
	gdouble s = g_ascii_strtod(gtk_entry_get_text(GtkEntryIPS_Dec_s), NULL);
	if (gtk_toggle_button_get_active(GtkCheckButtonIPS_S)) {
		deg = -deg;
	}

	return siril_world_cs_new_from_ra_dec(hour, min, sec, deg, m, s);;
}

static gboolean is_detection_manual() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("checkButton_IPS_manual"));
	return gtk_toggle_button_get_active(button);
}

static gboolean flip_image_after_ps() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("checkButton_IPS_flip"));
	return gtk_toggle_button_get_active(button);
}

static gboolean is_downsample_activated() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("downsample_ips_button"));
	return gtk_toggle_button_get_active(button);
}

static gboolean is_autocrop_activated() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("autocrop_ips_button"));
	return gtk_toggle_button_get_active(button);
}

static void update_pixel_size() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));
	float pixel;

	pixel = gfit.pixel_size_x > gfit.pixel_size_y ? gfit.pixel_size_x : gfit.pixel_size_y;

	if (pixel > 0.f) {
		gchar *cpixels = g_strdup_printf("%.2lf", (double) pixel);
		gtk_entry_set_text(entry, cpixels);
		g_free(cpixels);
	}
}

static void update_focal() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	double focal;

	focal = gfit.focal_length;

	if (focal > 0.0) {
		gchar *cfocal = g_strdup_printf("%.1lf", focal);
		gtk_entry_set_text(entry, cfocal);
		g_free(cfocal);
	}
}

static void update_resolution_field() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_resolution"));
	double res = get_resolution(get_focal(), get_pixel());
	gchar *cres;

	cres = g_strdup_printf("%1.3lf", res);
	gtk_entry_set_text(entry, cres);
	g_free(cres);
}

static void update_coordinates(SirilWorldCS *world_cs) {
	gchar *RA_sec, *Dec_sec;
	gint ra_h, ra_m;
	gint dec_deg, dec_m;
	gdouble ra_s, dec_s;

	siril_world_cs_get_ra_hour_min_sec(world_cs, &ra_h, &ra_m, &ra_s);
	siril_world_cs_get_dec_deg_min_sec(world_cs, &dec_deg, &dec_m, &dec_s);

	RA_sec = g_strdup_printf("%6.4lf", ra_s);
	Dec_sec = g_strdup_printf("%6.4lf", dec_s);

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButtonIPS_S")), dec_deg < 0);

	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_h")), ra_h);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_m")), ra_m);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_RA_s")), RA_sec);

	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_deg")), abs(dec_deg));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_m")), dec_m);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_Dec_s")), Dec_sec);

	g_free(RA_sec);
	g_free(Dec_sec);
}

void update_coords() {
	SirilWorldCS *world_cs = NULL;
	fits *fit = &gfit;

	if (fit->wcsdata.ra != 0.0 && fit->wcsdata.dec != 0.0) {

		world_cs = siril_world_cs_new_from_a_d(fit->wcsdata.ra, fit->wcsdata.dec);

		update_coordinates(world_cs);
		unselect_all_items();
	} else if (fit->wcsdata.objctra[0] != '\0' && fit->wcsdata.objctdec[0] != '\0') {

		world_cs = siril_world_cs_new_from_objct_ra_dec(fit->wcsdata.objctra, fit->wcsdata.objctdec);

		update_coordinates(world_cs);
		unselect_all_items();
	} else if (fit->wcsdata.crval[0] != 0.0 && fit->wcsdata.crval[1] != 0.0) {

		// first transform coords to alpha and delta
		world_cs = siril_world_cs_new_from_a_d(fit->wcsdata.crval[0], fit->wcsdata.crval[1]);

		update_coordinates(world_cs);
		unselect_all_items();
	}
	if (world_cs)
		siril_world_cs_unref(world_cs);
}

static void update_image_parameters_GUI() {
	/* update all fields. Resolution is updated as well
	 thanks to the Entry and combo changed signal
	  */
	update_focal();
	update_pixel_size();
	update_coords();
}

gboolean end_plate_solver(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *) p;
	stop_processing_thread();

	set_cursor_waiting(FALSE);

	if (args->ret) {
		char *title = siril_log_color_message(_("Plate Solving failed. "
				"The image could not be aligned with the reference stars.\n"), "red");
		if (!args->message) {
			args->message = g_strdup(_("This is usually because the initial parameters (pixel size, focal length, initial coordinates) "
					"are too far from the real metadata of the image.\n\n"
					"You could also try to look into another catalogue, or try to click on the \"Downsampling\" button, especially for image done with Drizzle.\n\n"
					"Finally, keep in mind that plate solving algorithm should only be applied on linear image."));
		}
		siril_message_dialog(GTK_MESSAGE_ERROR, title, args->message);
		g_free(args->message);
	} else {
		image_solved *solution = args->solution;

		/* update UI */
		update_image_parameters_GUI();
		set_GUI_CAMERA();
		update_coordinates(solution->image_center);
		delete_selected_area();
		/* ****** */

		control_window_switch_to_tab(OUTPUT_LOGS);
		if (args->flip_image || args->for_photometry_cc)
			redraw(REMAP_ALL);
		else redraw(REDRAW_OVERLAY);

		siril_world_cs_unref(solution->px_cat_center);
		siril_world_cs_unref(solution->image_center);
		free(solution);
	}
	if (args->image_flipped)
		clear_stars_list(TRUE);
	update_MenuItem();
	free(args);
	return FALSE;
}

static void start_image_plate_solve() {
	struct astrometry_data *args = calloc(1, sizeof(struct astrometry_data));

	args->for_photometry_cc = FALSE;
	if (!fill_plate_solver_structure_from_GUI(args)) {
		set_cursor_waiting(TRUE);
		start_in_new_thread(match_catalog, args);
	}
}

static void get_list_IPS() {
	if (list_IPS == NULL)
		list_IPS = GTK_LIST_STORE(gtk_builder_get_object(gui.builder, "liststoreIPS"));
}

static void clear_all_objects() {
	gtk_list_store_clear(list_IPS);
}

static void add_object_to_list() {
	GtkTreeIter iter;

	get_list_IPS();
	clear_all_objects();

	if (platedObject[RESOLVER_NED].name) {
		gtk_list_store_append(list_IPS, &iter);
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "NED", COLUMN_NAME,
				platedObject[RESOLVER_NED].name, -1);
	}

	if (platedObject[RESOLVER_SIMBAD].name) {
		gtk_list_store_append(list_IPS, &iter);
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "Simbad",
				COLUMN_NAME, platedObject[RESOLVER_SIMBAD].name, -1);
	}

	if (platedObject[RESOLVER_VIZIER].name) {
		gtk_list_store_append(list_IPS, &iter);
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "VizieR",
				COLUMN_NAME, platedObject[RESOLVER_VIZIER].name, -1);
	}

}

static void unselect_all_items() {
	GtkTreeSelection *selection;

	selection = GTK_TREE_SELECTION(gtk_builder_get_object(gui.builder, "gtkselectionIPS"));
	gtk_tree_selection_unselect_all(selection);
}

static void add_object_in_tree_view(const gchar *object) {
	struct sky_object obj;
	GtkTreeView *GtkTreeViewIPS;
	GtkTreeViewIPS = GTK_TREE_VIEW(lookup_widget("GtkTreeViewIPS"));
	if (!object || object[0] == '\0')
		return;

	set_cursor_waiting(TRUE);

	query_server server = get_server_from_combobox();
	gchar *result = search_in_catalogs(object, server);
	if (result) {
		free_Platedobject();
		parse_content_buffer(result, &obj);
		if (!has_nonzero_coords()) {
			g_free(result);
			set_cursor_waiting(FALSE);
			siril_log_color_message(_("No catalog\n"), "red");
			return;
		}
		g_signal_handlers_block_by_func(GtkTreeViewIPS, on_GtkTreeViewIPS_cursor_changed, NULL);
		add_object_to_list();
		g_signal_handlers_unblock_by_func(GtkTreeViewIPS, on_GtkTreeViewIPS_cursor_changed, NULL);
		g_free(result);

		/* select first object found in the list */
		GtkTreeIter iter;
		GtkTreeSelection *selection = gtk_tree_view_get_selection(GtkTreeViewIPS);
		if (gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_IPS), &iter)) {
			gtk_tree_selection_select_iter(selection, &iter);
			g_signal_emit_by_name(GTK_TREE_VIEW(GtkTreeViewIPS), "cursor-changed");
		}

	}
	set_cursor_waiting(FALSE);
}

/*****
 * CALLBACKS FUNCTIONS
 */

void on_GtkEntry_IPS_focal_changed(GtkEditable *editable, gpointer user_data) {
	update_resolution_field();
	com.pref.focal = g_ascii_strtod(gtk_editable_get_chars(editable, 0, -1), NULL);
}

void on_GtkEntry_IPS_pixels_changed(GtkEditable *editable, gpointer user_data) {
	update_resolution_field();
	com.pref.pitch = g_ascii_strtod(gtk_editable_get_chars(editable, 0, -1), NULL);
}

void on_GtkEntry_IPS_insert_text(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);
	int i, count = 0;

	gchar *result = g_strndup(text, length);

	for (i = 0; i < length; i++) {
		if (!g_ascii_isdigit(text[i]) && text[i] != '.')
			continue;
		result[count++] = text[i];
	}

	if (count > 0) {
		g_signal_handlers_block_by_func(G_OBJECT (editable),
				G_CALLBACK (on_GtkEntry_IPS_insert_text), data);
		gtk_editable_insert_text(editable, result, count, position);
		g_signal_handlers_unblock_by_func(G_OBJECT (editable),
				G_CALLBACK (on_GtkEntry_IPS_insert_text), data);
	}
	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");

	g_free(result);
}

void on_buttonIPS_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("ImagePlateSolver_Dial");
}

void on_GtkTreeViewIPS_cursor_changed(GtkTreeView *tree_view, gpointer user_data) {

	GtkTreeModel *treeModel = gtk_tree_view_get_model(tree_view);
	GtkTreeSelection *selection = gtk_tree_view_get_selection (tree_view);
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;

	if (gtk_tree_model_get_iter_first(treeModel, &iter) == FALSE)
		return;	//The tree is empty
	if (gtk_tree_selection_get_selected(selection, &treeModel, &iter)) { //get selected item
		gtk_tree_model_get_value(treeModel, &iter, COLUMN_RESOLVER, &value);
		const gchar *res = g_value_get_string(&value);
		int selected_item;

		if (!g_strcmp0(res, "NED")) {
			selected_item = 0;
		} else if (!g_strcmp0(res, "Simbad")) {
			selected_item = 1;
		} else if (!g_strcmp0(res, "VizieR")) {
			selected_item = 2;
		} else {
			selected_item = -1;
		}

		if (selected_item >= 0) {
			if (platedObject[selected_item].world_cs)
				update_coordinates(platedObject[selected_item].world_cs);
			else {
				char *msg = siril_log_message(_("There are no available coordinates with this name, try with another name\n"));
				siril_message_dialog(GTK_MESSAGE_WARNING, _("No coordinates"), msg);
			}
		}

		g_value_unset(&value);
	}
}

void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data) {
	if (!has_any_keywords()) {
		char *msg = siril_log_message(_("There are no keywords stored in the FITS header.\n"));
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No metadata"), msg);
	} else {
		update_image_parameters_GUI();
	}
}

void on_GtkButtonIPS_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;

	entry = GTK_ENTRY(lookup_widget("GtkSearchIPS"));
	add_object_in_tree_view(gtk_entry_get_text(GTK_ENTRY(entry)));
}

void on_buttonIPS_ok_clicked(GtkButton *button, gpointer user_data) {
	start_image_plate_solve();
}

void on_GtkSearchIPS_activate(GtkEntry *entry, gpointer user_data) {
	add_object_in_tree_view(gtk_entry_get_text(GTK_ENTRY(entry)));
}

void on_GtkCheckButton_Mag_Limit_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *spinmag;

	spinmag = lookup_widget("GtkSpinIPS_Mag_Limit");
	gtk_widget_set_sensitive(spinmag, !gtk_toggle_button_get_active(button));
}

void on_GtkCheckButton_OnlineCat_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *combobox;

	combobox = lookup_widget("ComboBoxIPSCatalog");
	gtk_widget_set_sensitive(combobox, !gtk_toggle_button_get_active(button));
}


void open_astrometry_dialog() {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		initialize_ips_dialog();
		siril_open_dialog("ImagePlateSolver_Dial");
	}
}

int fill_plate_solver_structure_from_GUI(struct astrometry_data *args) {
	args->fit = &gfit;
	args->pixel_size = get_pixel();
	args->focal_length = get_focal();
	// args->onlineCatalog requires processed arguments in automatic mode, setting later

	SirilWorldCS *catalog_center = get_center_of_catalog();
	if (siril_world_cs_get_alpha(catalog_center) == 0.0 && siril_world_cs_get_delta(catalog_center) == 0.0) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No coordinates"), _("Please enter object coordinates."));
		set_cursor_waiting(FALSE);
		return 1;
	}
	args->cat_center = catalog_center;

	args->downsample = is_downsample_activated();
	args->autocrop = is_autocrop_activated();
	args->flip_image = flip_image_after_ps();
	args->manual = is_detection_manual();
	get_mag_from_GUI(args);

	set_cursor_waiting(TRUE);

	process_plate_solver_input(args);

	args->onlineCatalog = args->for_photometry_cc ? get_photometry_catalog() :
		get_online_catalog(args->used_fov * CROP_ALLOWANCE, args->limit_mag);

	/* currently the GUI version downloads the catalog here, because
	 * siril_message_dialog() doesn't use idle function, we could change that */
	GFile *catalog_file = download_catalog(args->onlineCatalog, catalog_center, args->used_fov * CROP_ALLOWANCE, args->limit_mag);
	if (!catalog_file) {
		siril_world_cs_unref(catalog_center);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No catalog"), _("Cannot download the online star catalog."));
		set_cursor_waiting(FALSE);
		return 1;
	}
	args->catalog_file = catalog_file;

	set_cursor_waiting(FALSE);
	return 0;
}

gboolean confirm_delete_wcs_keywords(fits *fit) {
	gboolean erase = TRUE;

	if (fit->wcsdata.equinox > 0.0) {
		erase = siril_confirm_dialog(_("Astrometric solution detected"),
				_("The astrometric solution contained in "
				"the image will be erased by the geometric transformation and no undo "
				"will be possible."), _("Process"));
	}
	return erase;
}

void set_focal_and_pixel_pitch() {
	GtkEntry *focal, *pitch;
	gchar *f_str, *p_str;

	focal = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	pitch = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));

	f_str = g_strdup_printf("%.1lf", com.pref.focal);
	p_str = g_strdup_printf("%.2lf", com.pref.pitch);

	gtk_entry_set_text(focal, f_str);
	gtk_entry_set_text(pitch, p_str);

	g_free(f_str);
	g_free(p_str);
}

