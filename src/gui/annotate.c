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
#include "gui/utils.h"
#include "io/local_catalogues.h"
#include "io/siril_catalogues.h"

#include "annotate.h"

enum {
	CONESEARCH_PAGE,
	SHOW_PAGE,
	SEARCH_OBJECT_PAGE
};

// Statics declarations

GtkButton *annotate_close = NULL, *annotate_apply = NULL, *annotate_save_as_button = NULL, *show_button = NULL;
GtkSearchEntry *search_objects_entry = NULL;
GtkWindow *annotate_dialog = NULL;
GtkLabel *show_ra = NULL, *show_dec = NULL;
GtkEntry *conesearch_save_entry = NULL, *annotate_obscode_entry = NULL, *show_file_entry = NULL, *show_ra_entry = NULL, *show_dec_entry = NULL, *show_name_entry = NULL;
GtkComboBoxText *conesearch_combo = NULL;
GtkAdjustment *adj_mag_limit = NULL;
GtkSpinButton *conesearch_maglimit = NULL;
GtkStack *stack_show = NULL;
GtkNotebook *notebook_annotate = NULL;
GtkToggleButton *conesearch_photometric = NULL, *conesearch_tag = NULL, *conesearch_log = NULL, *show_tag = NULL, *show_log = NULL, *show_clear = NULL;


// Statics init

static void annotate_dialog_init_statics() {
	if (adj_mag_limit == NULL) {
		// GtkButton
		annotate_close = GTK_BUTTON(lookup_widget("annotate_close"));
		annotate_apply = GTK_BUTTON(lookup_widget("annotate_apply"));
		annotate_save_as_button = GTK_BUTTON(lookup_widget("annotate_save_as_button"));
		show_button = GTK_BUTTON(lookup_widget("show_button"));
		// GtkSearchEntry
		search_objects_entry = GTK_SEARCH_ENTRY(lookup_widget("search_objects_entry"));
		// GtkWindow
		annotate_dialog = GTK_WINDOW(lookup_widget("annotate_dialog"));
		// GtkLabel
		show_ra = GTK_LABEL(lookup_widget("show_ra"));
		show_dec = GTK_LABEL(lookup_widget("show_dec"));
		// GtkEntry
		conesearch_save_entry = GTK_ENTRY(lookup_widget("conesearch_save_entry"));
		annotate_obscode_entry = GTK_ENTRY(lookup_widget("annotate_obscode_entry"));
		show_file_entry = GTK_ENTRY(lookup_widget("show_file_entry"));
		show_ra_entry = GTK_ENTRY(lookup_widget("show_ra_entry"));
		show_dec_entry = GTK_ENTRY(lookup_widget("show_dec_entry"));
		show_name_entry = GTK_ENTRY(lookup_widget("show_name_entry"));
		// GtkComboBoxText
		conesearch_combo = GTK_COMBO_BOX_TEXT(lookup_widget("conesearch_combo"));
		// GtkAdjustment
		adj_mag_limit = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adj_mag_limit"));
		// GtkSpinButton
		conesearch_maglimit = GTK_SPIN_BUTTON(lookup_widget("conesearch_maglimit"));
		// GtkStack
		stack_show = GTK_STACK(lookup_widget("stack_show"));
		// GtkNotebook
		notebook_annotate = GTK_NOTEBOOK(lookup_widget("notebook_annotate"));
		// GtkToggleButton
		conesearch_photometric = GTK_TOGGLE_BUTTON(lookup_widget("conesearch_photometric"));
		conesearch_tag = GTK_TOGGLE_BUTTON(lookup_widget("conesearch_tag"));
		conesearch_log = GTK_TOGGLE_BUTTON(lookup_widget("conesearch_log"));
		show_tag = GTK_TOGGLE_BUTTON(lookup_widget("show_tag"));
		show_log = GTK_TOGGLE_BUTTON(lookup_widget("show_log"));
		show_clear = GTK_TOGGLE_BUTTON(lookup_widget("show_clear"));
	}
}


static siril_cat_index get_cat_index_from_combo() {
	siril_cat_index cat = CAT_AUTO;
	const gchar *cat_char = gtk_combo_box_text_get_active_text(conesearch_combo);
	if (!g_strcmp0(cat_char, "tycho2"))
		cat = CAT_TYCHO2;
	else if (!g_strcmp0(cat_char, "nomad"))
		cat = CAT_NOMAD;
	else if (!g_strcmp0(cat_char, "gaia"))
		cat = CAT_GAIADR3;
	else if (!g_strcmp0(cat_char, "ppmxl"))
		cat = CAT_PPMXL;
	else if (!g_strcmp0(cat_char, "bsc"))
		cat = CAT_BSC;
	else if (!g_strcmp0(cat_char, "apass"))
		cat = CAT_APASS;
	else if (!g_strcmp0(cat_char, "gcvs"))
		cat = CAT_GCVS;
	else if (!g_strcmp0(cat_char, "vsx"))
		cat = CAT_VSX;
	else if (!g_strcmp0(cat_char, "varisum"))
		cat = CAT_VARISUM;
	else if (!g_strcmp0(cat_char, "simbad"))
		cat = CAT_SIMBAD;
	else if (!g_strcmp0(cat_char, "exo"))
		cat = CAT_EXOPLANETARCHIVE;
	else if (!g_strcmp0(cat_char, "pgc"))
		cat = CAT_PGC;
	else if (!g_strcmp0(cat_char, "aavso_chart"))
		cat = CAT_AAVSO_CHART;
	else if (!g_strcmp0(cat_char, "solsys")) {
		cat = CAT_IMCCE;
		if (!gfit.keywords.date_obs) {
			siril_log_color_message(_("This option only works on images that have observation date information\n"), "red");
			return CAT_AUTO;
		}
	} else {
		cat = CAT_AUTO;
	}

	return cat;
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
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
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

static conesearch_params* parse_conesearch_ui() {
	conesearch_params *params = g_new0(conesearch_params, 1);
	params->limit_mag = -1.0f;
	params->photometric = FALSE;
	params->display_tag = BOOL_NOT_SET;
	params->display_log = BOOL_NOT_SET;
	params->cat = CAT_AUTO;
	params->obscode = NULL;
	params->default_obscode_used = FALSE;
	params->trixel = -1;
	params->outfilename = NULL;

	gboolean local_cat = local_catalogues_available();

	if (com.pref.astrometry.default_obscode != NULL) {
		params->obscode = g_strdup(com.pref.astrometry.default_obscode);
		params->default_obscode_used = TRUE;
	}
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

	if (params->cat == CAT_AUTO) {
		params->cat = (local_cat) ? CAT_LOCAL : CAT_NOMAD;
		if (params->trixel >= 0 && params->cat == CAT_LOCAL)
			params->cat = CAT_LOCAL_TRIX;
	}

	return params;
}

int execute_conesearch(conesearch_params *params) {
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}

	// Preparing the catalogue query
	siril_catalogue *siril_cat = siril_catalog_fill_from_fit(&gfit, params->cat, params->limit_mag);
	siril_cat->phot = params->photometric;
	if (params->cat == CAT_IMCCE) {
		if (params->obscode) {
			siril_cat->IAUcode = params->obscode;
			if (params->default_obscode_used) {
				siril_log_message(_("Using default observatory code %s\n"), params->obscode);
			}
		} else {
			siril_cat->IAUcode = g_strdup("500");
			siril_log_color_message(_("Did not specify an observatory code, using geocentric by default, positions may not be accurate\n"), "salmon");
		}
	} else if (params->obscode) {
		g_free(params->obscode);
		params->obscode = NULL;
	}
	if (params->cat == CAT_LOCAL_TRIX)
		siril_cat->trixel = params->trixel;

	siril_debug_print("centre coords: %f, %f, radius: %f arcmin\n", siril_cat->center_ra, siril_cat->center_dec, siril_cat->radius);
	conesearch_args *args = init_conesearch();
	args->fit = &gfit;
	args->siril_cat = siril_cat;
	args->has_GUI = !com.script;
	args->display_log = (params->display_log == BOOL_NOT_SET) ? display_names_for_catalogue(params->cat) : (gboolean) params->display_log;
	args->display_tag = (params->display_tag == BOOL_NOT_SET) ? display_names_for_catalogue(params->cat) : (gboolean) params->display_tag;
	args->outfilename = params->outfilename;
	if (check_conesearch_args(args)) { // can't fail for now
		free_conesearch(args);
		return CMD_GENERIC_ERROR;
	}
	start_in_new_thread(conesearch_worker, args);
	return CMD_OK;
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
		}
		params->coords = NULL;
		params->name = NULL;
	} else {
		params->file = NULL;
		const gchar *ra = gtk_entry_get_text(show_ra_entry);
		const gchar *dec = gtk_entry_get_text(show_dec_entry);

		if (is_string_numeric(ra) && is_string_numeric(dec)) {
			params->coords = siril_world_cs_new_from_objct_ra_dec((gchar *)ra, (gchar *)dec);
		} else {
			params->coords = NULL;
		}
	    const gchar *name_entry = gtk_entry_get_text(show_name_entry);
		if (name_entry != NULL && strlen(name_entry) != 0) {
			params->name = g_strdup(name_entry);
		}
	}

    params->clear = gtk_toggle_button_get_active(show_clear);
    params->display_log = gtk_toggle_button_get_active(show_log);
    params->display_tag = gtk_toggle_button_get_active(show_tag);


	return params;
}

int execute_show_command(show_params *params) {
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}

	if (params->clear) {
		purge_user_catalogue(CAT_AN_USER_TEMP);
		redraw(REDRAW_OVERLAY);
	}

	siril_catalogue *siril_cat = calloc(1, sizeof(siril_catalogue));
	siril_cat->cat_index = CAT_SHOW;
	siril_cat->columns = siril_catalog_columns(siril_cat->cat_index);
	conesearch_args *args = init_conesearch();
	args->siril_cat = siril_cat;
	args->has_GUI = TRUE;
	args->fit = &gfit;

	if (params->file) {
		int check = siril_catalog_load_from_file(siril_cat, params->file);
		if (check > 0) {
			free_conesearch(args);
			return CMD_ARG_ERROR;
		}
		if (check == -1) {
			free_conesearch(args);
			return CMD_OK;
		}
		args->display_log =
				(params->display_log == BOOL_NOT_SET) ?
						(gboolean) has_field(siril_cat, NAME) :
						(gboolean) params->display_log;
		args->display_tag =
				(params->display_tag == BOOL_NOT_SET) ?
						(gboolean) has_field(siril_cat, NAME) :
						(gboolean) params->display_tag;
		start_in_new_thread(conesearch_worker, args);
		return CMD_OK;
	}

	if (params->coords) {
		cat_item *item = calloc(1, sizeof(cat_item));
		item->ra = siril_world_cs_get_alpha(params->coords);
		item->dec = siril_world_cs_get_delta(params->coords);
		siril_world_cs_unref(params->coords);
		item->name = g_strdup(params->name);
		siril_cat->columns |= (1 << CAT_FIELD_NAME);
		args->display_log = TRUE;
		args->display_tag = TRUE;

		siril_catalog_append_item(siril_cat, item);
		siril_catalog_free_item(item);
		free(item);
		start_in_new_thread(conesearch_worker, args);
		return CMD_OK;
	}

	free_conesearch(args);
	return CMD_ARG_ERROR;
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
			return;
		}
		execute_conesearch(params_cone);
		g_free(params_cone);
		break;
	case SHOW_PAGE:
		params_show = parse_show_ui();
		if (!params_show) {
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
