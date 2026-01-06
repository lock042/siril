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

#include <ctype.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "algos/search_objects.h"
#include "algos/siril_wcs.h"
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "io/local_catalogues.h"
#include "io/siril_catalogues.h"


enum {
	CONESEARCH_PAGE,
	SHOW_PAGE,
	SEARCH_OBJECT_PAGE
};

// Statics declarations

GtkButton *annotate_clear = NULL, *annotate_close = NULL, *annotate_apply = NULL, *annotate_save_as_button = NULL, *show_button = NULL, *show_button_save_to_DSO = NULL;
GtkSearchEntry *search_objects_entry = NULL;
GtkWindow *annotate_dialog = NULL;
GtkEntry *conesearch_save_entry = NULL, *annotate_obscode_entry = NULL, *show_file_entry = NULL, *show_ra_entry = NULL, *show_dec_entry = NULL, *show_name_entry = NULL;
GtkComboBoxText *conesearch_combo = NULL;
GtkAdjustment *adj_mag_limit = NULL;
GtkSpinButton *conesearch_maglimit = NULL;
GtkStack *stack_show = NULL;
GtkNotebook *notebook_annotate = NULL;
GtkToggleButton *conesearch_photometric = NULL, *conesearch_tag = NULL, *conesearch_log = NULL, *show_tag = NULL, *show_log = NULL;
GtkBox *annotate_obscode_box = NULL;

static int local_cat = BOOL_NOT_SET;

// Statics init

void annotate_dialog_init_statics() {
	if (adj_mag_limit == NULL) {
		// GtkNotebook
		notebook_annotate = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook_annotate"));
		// GtkEntry
		conesearch_save_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "conesearch_save_entry"));
		annotate_obscode_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "annotate_obscode_entry"));
		show_file_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "show_file_entry"));
		show_ra_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "show_ra_entry"));
		show_dec_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "show_dec_entry"));
		show_name_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "show_name_entry"));
		// GtkDialog
		annotate_dialog = GTK_WINDOW(gtk_builder_get_object(gui.builder, "annotate_dialog"));
		// GtkSearchEntry
		search_objects_entry = GTK_SEARCH_ENTRY(gtk_builder_get_object(gui.builder, "search_objects_entry"));
		// GtkButton
		annotate_clear = GTK_BUTTON(gtk_builder_get_object(gui.builder, "annotate_clear"));
		annotate_close = GTK_BUTTON(gtk_builder_get_object(gui.builder, "annotate_close"));
		annotate_apply = GTK_BUTTON(gtk_builder_get_object(gui.builder, "annotate_apply"));
		annotate_save_as_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "annotate_save_as_button"));
		show_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "show_button"));
		show_button_save_to_DSO = GTK_BUTTON(gtk_builder_get_object(gui.builder, "show_button_save_to_DSO"));
		// GtkSpinButton
		conesearch_maglimit = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "conesearch_maglimit"));
		// GtkStack
		stack_show = GTK_STACK(gtk_builder_get_object(gui.builder, "stack_show"));
		// GtkComboBoxText
		conesearch_combo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(gui.builder, "conesearch_combo"));
		// GtkAdjustment
		adj_mag_limit = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adj_mag_limit"));
		// GtkToggleButton
		conesearch_photometric = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "conesearch_photometric"));
		conesearch_tag = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "conesearch_tag"));
		conesearch_log = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "conesearch_log"));
		show_tag = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "show_tag"));
		show_log = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "show_log"));
		// GtkBox
		annotate_obscode_box = GTK_BOX(gtk_builder_get_object(gui.builder, "annotate_obscode_box"));
	}
}

static siril_cat_index get_cat_index_from_combo() {
	siril_cat_index cat = CAT_AUTO;
	const gchar *cat_char = gtk_combo_box_text_get_active_text(conesearch_combo);
	if (g_str_has_prefix(cat_char, "tycho2"))
		cat = CAT_TYCHO2;
	else if (g_str_has_prefix(cat_char, "nomad"))
		cat = CAT_NOMAD;
	else if (g_str_has_prefix(cat_char, "gaia"))
		cat = CAT_GAIADR3;
	else if (g_str_has_prefix(cat_char, "ppmxl"))
		cat = CAT_PPMXL;
	else if (g_str_has_prefix(cat_char, "bsc"))
		cat = CAT_BSC;
	else if (g_str_has_prefix(cat_char, "apass"))
		cat = CAT_APASS;
	else if (g_str_has_prefix(cat_char, "gcvs"))
		cat = CAT_GCVS;
	else if (g_str_has_prefix(cat_char, "vsx"))
		cat = CAT_VSX;
	else if (g_str_has_prefix(cat_char, "varisum"))
		cat = CAT_VARISUM;
	else if (g_str_has_prefix(cat_char, "simbad"))
		cat = CAT_SIMBAD;
	else if (g_str_has_prefix(cat_char, "exo"))
		cat = CAT_EXOPLANETARCHIVE;
	else if (g_str_has_prefix(cat_char, "pgc"))
		cat = CAT_PGC;
	else if (g_str_has_prefix(cat_char, "aavso_chart"))
		cat = CAT_AAVSO_CHART;
	else if (g_str_has_prefix(cat_char, "solsys")) {
		cat = CAT_IMCCE;
		if (!gfit->keywords.date_obs) {
			siril_log_color_message(_("This option only works on images that have observation date information\n"), "red");
			return CAT_AUTO;
		}
	} else {
		cat = CAT_AUTO;
	}
	if (local_cat == BOOL_NOT_SET)
		local_cat = (int)local_catalogues_available();
	if (cat == CAT_AUTO) {
		cat = (local_cat) ? get_local_catalogue_index() : CAT_NOMAD;
	}

	return cat;
}

void on_conesearch_combo_changed(GtkComboBox *widget, gpointer user_data) {
	siril_cat_index cat = get_cat_index_from_combo();
	uint32_t columns = siril_catalog_columns(cat);
	gboolean has_phot = has_field_from_columns(columns, BMAG);
	gboolean has_names = has_field_from_columns(columns, NAME);
	gboolean logtag = display_names_for_catalogue(cat);

	float mag = siril_catalog_get_default_limit_mag(cat);
	gtk_widget_set_sensitive(GTK_WIDGET(conesearch_tag), has_names);
	gtk_widget_set_sensitive(GTK_WIDGET(conesearch_log), has_names);
	gtk_toggle_button_set_active(conesearch_tag, logtag);
	gtk_toggle_button_set_active(conesearch_log, logtag);
	gtk_widget_set_sensitive(GTK_WIDGET(conesearch_photometric), has_phot);
	gtk_widget_set_visible(GTK_WIDGET(annotate_obscode_box), cat == CAT_IMCCE);
	gtk_widget_set_sensitive(GTK_WIDGET(conesearch_maglimit), mag != 0.f);
	if (mag != 0) {
		gtk_spin_button_set_value(conesearch_maglimit, mag);
	}
}

static void set_filter(GtkFileChooser *dialog) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("csv files (*.csv)"));
	gtk_file_filter_add_pattern(f, "*.csv");
	gtk_file_chooser_add_filter(dialog, f);
	gtk_file_chooser_set_filter(dialog, f);
}

void on_annotate_dialog_show(GtkWidget *widget, gpointer user_data) {
	annotate_dialog_init_statics();
	if (com.pref.astrometry.default_obscode != NULL) {
		gtk_entry_set_text(annotate_obscode_entry, com.pref.astrometry.default_obscode);
	}
	on_conesearch_combo_changed(NULL, NULL);
}

void on_annotate_save_as_button_clicked(GtkButton *button, gpointer user_data) {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	gint res;
	gchar *filename;

	filename = g_strdup(".csv");

	widgetdialog = siril_file_chooser_save(annotate_dialog, GTK_FILE_CHOOSER_ACTION_SAVE);
	dialog = GTK_FILE_CHOOSER(widgetdialog);
	gtk_file_chooser_set_current_folder(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, FALSE); // the overwrite is checked when applied!
	gtk_file_chooser_set_current_name(dialog, filename);
	gtk_file_chooser_set_local_only(dialog, FALSE);
	set_filter(dialog);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_file_chooser_get_filename(dialog);
		gtk_entry_set_text(conesearch_save_entry, file);
		gtk_editable_set_position(GTK_EDITABLE(conesearch_save_entry), -1);
		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
	g_free(filename);
}

void on_show_button_clicked(GtkButton *button, gpointer user_data) {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	gint res;

	widgetdialog = siril_file_chooser_open(annotate_dialog, GTK_FILE_CHOOSER_ACTION_OPEN);
	dialog = GTK_FILE_CHOOSER(widgetdialog);
	gtk_file_chooser_set_current_folder(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	gtk_file_chooser_set_local_only(dialog, FALSE);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_file_chooser_get_filename(dialog);
		gtk_entry_set_text(show_file_entry, file);
		gtk_editable_set_position(GTK_EDITABLE(show_file_entry), -1);
		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
}

void on_show_button_get_coords_clicked(GtkButton *button, gpointer user_data) {
	if (has_wcs(gfit) && (com.selection.h && com.selection.w)) {
		psf_error error = PSF_NO_ERR;
		psf_star *result = psf_get_minimisation(gfit, select_vport(gui.cvport), &com.selection, FALSE, FALSE, NULL, FALSE, com.pref.starfinder_conf.profile, &error);
		if (result && error == PSF_NO_ERR) {
			double world_x, world_y;
			gchar *ra, *dec;

			result->xpos = result->x0 + com.selection.x;
			if (gfit->top_down)
				result->ypos = result->y0 + com.selection.y;
			else
				result->ypos = com.selection.y + com.selection.h - result->y0;

			pix2wcs(gfit, result->xpos, (double) gfit->ry - result->ypos - 1.0, &world_x, &world_y);
			if (world_x >= 0.0 && !isnan(world_x) && !isnan(world_y)) {
				SirilWorldCS *world_cs = siril_world_cs_new_from_a_d(world_x, world_y);
				if (world_cs) {
					ra = siril_world_cs_alpha_format(world_cs, "%02d %02d %.3lf");
					dec = siril_world_cs_delta_format(world_cs, "%c%02d %02d %.3lf");

					gtk_entry_set_text(GTK_ENTRY(lookup_widget("show_ra_entry")), ra);
					gtk_entry_set_text(GTK_ENTRY(lookup_widget("show_dec_entry")), dec);

					g_free(ra), g_free(dec);
					siril_world_cs_unref(world_cs);
				}
			}
		}
		free_psf(result);
	}
}

static int collect_single_coords_and_name(double *ra, double *dec, gchar **name) {
	const gchar *ra_str = gtk_entry_get_text(show_ra_entry);
	const gchar *dec_str = gtk_entry_get_text(show_dec_entry);

	SirilWorldCS *coords = siril_world_cs_new_from_objct_ra_dec((gchar *)ra_str, (gchar *)dec_str);
	if (!coords)
		return 1;
	*ra = siril_world_cs_get_alpha(coords);
	*dec = siril_world_cs_get_delta(coords);
	siril_world_cs_unref(coords);

	const gchar *name_entry = gtk_entry_get_text(show_name_entry);
	if (name_entry != NULL && strlen(name_entry) != 0) {
		*name = g_strdup(name_entry);
	} else {
		*name = NULL;
	}
	return 0;
}

void on_show_button_save_to_DSO_clicked(GtkButton *button, gpointer user_data) {
	double ra, dec;
	gchar *name = NULL;
	if (collect_single_coords_and_name(&ra, &dec, &name) || !name) {
		siril_log_color_message(_("Could not parse coordinates or name, aborting\n"), "red");
		return;
	}
	cat_item *item = calloc(1, sizeof(cat_item));
	item->name = name;
	item->ra = ra;
	item->dec = dec;
	add_item_in_catalogue(item, CAT_AN_USER_DSO, TRUE);
	set_annotation_visibility(CAT_AN_USER_DSO, TRUE);	// and display it
	siril_catalog_free_item(item);
	free(item);
	refresh_found_objects();
}

static conesearch_params* parse_conesearch_ui() {
	conesearch_params *params = init_conesearch_params();

	params->cat = get_cat_index_from_combo();
	params->display_log = gtk_toggle_button_get_active(conesearch_log) ? BOOL_TRUE : BOOL_FALSE;
	params->display_tag = gtk_toggle_button_get_active(conesearch_tag) ? BOOL_TRUE : BOOL_FALSE;
	params->photometric = gtk_toggle_button_get_active(conesearch_photometric);
	params->limit_mag = gtk_spin_button_get_value(conesearch_maglimit);
	const gchar *output = gtk_entry_get_text(conesearch_save_entry);

	if (output != NULL && strlen(output) != 0) {
		params->outfilename = g_strdup(output);
	}
	const gchar *obscode = gtk_entry_get_text(annotate_obscode_entry);
	if (obscode != NULL && strlen(obscode) != 0) {
		params->default_obscode_used = FALSE;
		params->obscode = g_strdup(obscode);
	}
	return params;
}

static show_params* parse_show_ui() {
	show_params *params = g_new0(show_params, 1);
	params->display_log = BOOL_NOT_SET;
	params->display_tag = BOOL_NOT_SET;
	params->name = NULL;

	const gchar *visible_child_name = gtk_stack_get_visible_child_name(stack_show);

	if (!g_strcmp0("page0", visible_child_name)) {
		const gchar *input = gtk_entry_get_text(show_file_entry);
		if (input != NULL && strlen(input) != 0) {
			params->file = g_strdup(input);
		} else {
			g_free(params);
			return NULL;
		}
		params->coords = NULL;
		params->name = NULL;
		params->display_log = gtk_toggle_button_get_active(show_log);
		params->display_tag = gtk_toggle_button_get_active(show_tag);
	} else {
		params->file = NULL;
		double ra = 0., dec = 0.;
		if (collect_single_coords_and_name(&ra, &dec, &params->name)) {
			g_free(params);
			return NULL;
		}
		params->coords = siril_world_cs_new_from_a_d(ra, dec);
	}
	return params;
}

void on_annotate_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry = NULL;
	conesearch_params *params_cone = NULL;
	show_params *params_show = NULL;

	set_cursor_waiting(TRUE);
	int page = gtk_notebook_get_current_page(notebook_annotate);
	switch (page) {
	case CONESEARCH_PAGE:
		params_cone = parse_conesearch_ui();
		if (!params_cone) {
			set_cursor_waiting(FALSE);
			return;
		}
		if (params_cone->outfilename && g_file_test(params_cone->outfilename, G_FILE_TEST_EXISTS)) {
			gchar *basename = g_path_get_basename(params_cone->outfilename);
			gchar *dir_path = g_path_get_dirname(params_cone->outfilename);
			gchar *last_dir = g_path_get_basename(dir_path);

			gchar *title = g_strdup_printf("A file named \"%s\" already exists. Do you want to replace it?", basename);
			gchar *txt = g_strdup_printf("The file already exists in \"%s\". Replacing it will overwrite its contents.", last_dir);

			if (!siril_confirm_dialog(N_(title), N_(txt), _("Replace"))) {
				set_cursor_waiting(FALSE);
				g_free(basename); g_free(dir_path); g_free(last_dir); g_free(title); g_free(txt);
				g_free(params_cone);
				return;
			}
			g_free(basename); g_free(dir_path); g_free(last_dir); g_free(title); g_free(txt);
		}
		// Allocate and initialize generic_img_args
		struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
		if (!args) {
			free_conesearch_params(params_cone);
			PRINT_ALLOC_ERR;
			return;
		}
		args->fit = gfit;
		args->mem_ratio = 1.0f;
		args->image_hook = conesearch_image_hook;
		args->description = _("Cone search");
		args->verbose = TRUE;
		args->command_updates_gfit = FALSE;
		args->command = TRUE;
		args->user = params_cone;
		args->log_hook = NULL;

		if (!start_in_new_thread(generic_image_worker, args)) {
			siril_log_color_message(_("Error: failed to start conesearch image worker\n"), "red");
			free_generic_img_args(args);
			return;
		}
		break;
	case SHOW_PAGE:
		params_show = parse_show_ui();
		if (!params_show) {
			set_cursor_waiting(FALSE);
			return;
		}
		execute_show_command(params_show);
		g_free(params_show->name);
		g_free(params_show->file);
		g_free(params_show);
		break;
	case SEARCH_OBJECT_PAGE:
		entry = GTK_ENTRY(search_objects_entry);
		search_object(entry);
		break;
	}
	set_cursor_waiting(FALSE);
}

void on_annotate_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("annotate_dialog");
}

void on_annotate_clear_clicked(GtkButton *button, gpointer user_data) {
	purge_user_catalogue(CAT_AN_USER_TEMP);
	redraw(REDRAW_OVERLAY);
}
