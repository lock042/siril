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
#include "core/siril_networking.h"
#include "algos/colors.h"
//#include "algos/siril_wcs.h"
#include "io/local_catalogues.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/dialogs.h"
#include "gui/PSF_list.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *tempbox = NULL;
static GtkWidget *dtempbox = NULL;
static GtkWidget *check_filter = NULL;
static GtkWidget *temp_entry = NULL;
static GtkWidget *dtemp_entry = NULL;

static struct catmag_data args;

static void on_catmag_response(GtkDialog* self, gint response_id, gpointer user_data);

static void filter_checked(GtkToggleButton *source, gpointer user_data) {
	gtk_widget_set_sensitive(tempbox, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_filter)));
	gtk_widget_set_sensitive(dtempbox, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_filter)));
}

static void build_the_dialog() {
	dialog = gtk_dialog_new_with_buttons(_("Calibrate image magnitudes"), NULL,
			0, _("_Close"), GTK_RESPONSE_REJECT, _("_OK"), GTK_RESPONSE_ACCEPT, NULL);
	// If the user clicks one of these dialog buttons, GtkDialog will emit
	// the GtkDialog::response signal with the corresponding response ID
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
	g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(siril_widget_hide_on_delete), NULL);
	g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_catmag_response), NULL);

	/********************************************
	 *    -Image must be plate solved first-    *
	 *      Using catalogue: <detected>         *
	 *                                          *
	 *      [x] Filter stars                    *
	 *          Temperature: [ 5550 ] K         *
	 *          Allowed range: [ 500 ] K        *
	 *                                          *
	 *                                   OK     *
	 ********************************************/

	// Sets the "OK" button as default one
	GtkWidget *OK_button = gtk_dialog_get_widget_for_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
	gtk_widget_set_can_default(OK_button, TRUE);
	gtk_widget_grab_default(OK_button);
	gtk_widget_grab_focus(OK_button);
	gtk_style_context_add_class(gtk_widget_get_style_context(OK_button), "suggested-action");

	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12); // the main box

	/*if (!has_wcs(gfit)) {
		// unfortunately this is not called when the window is opened but later, so
		// this test doesn't work
		GtkWidget *label_platesolve = gtk_label_new("");
		gtk_label_set_markup(GTK_LABEL(label_platesolve), _("<b>Image must be plate solved first</b>"));
		gtk_container_add(GTK_CONTAINER(box), label_platesolve);
		g_object_set(G_OBJECT(label_platesolve), "margin-top", 20, NULL);
		g_object_set(G_OBJECT(label_platesolve), "margin-bottom", 10, NULL);
	}*/

	GtkWidget *catbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	GtkWidget *label_cat = gtk_label_new(_("Using catalogue:"));
	gtk_container_add(GTK_CONTAINER(catbox), label_cat);
	GtkWidget *label_cat_name = gtk_label_new("");
	gtk_container_add(GTK_CONTAINER(catbox), label_cat_name);
	gtk_container_add(GTK_CONTAINER(box), catbox);
	g_object_set(G_OBJECT(catbox), "margin-top", 20, NULL);
	g_object_set(G_OBJECT(catbox), "margin-left", 15, NULL);

	check_filter = gtk_check_button_new_with_label(_("Filter stars"));
	gtk_container_add(GTK_CONTAINER(box), check_filter);
	g_signal_connect (check_filter, "toggled", G_CALLBACK (filter_checked), NULL);
	g_object_set(G_OBJECT(check_filter), "margin-left", 12, NULL);

	tempbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	GtkWidget *label_temp = gtk_label_new(_("Temperature:"));
	temp_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(temp_entry), "5550.0");
	gtk_widget_set_tooltip_text(temp_entry, _("Stars with temperatures around this value will be used"));
	GtkWidget *label_K1 = gtk_label_new("K");
	gtk_container_add(GTK_CONTAINER(tempbox), label_temp);
	gtk_container_add(GTK_CONTAINER(tempbox), temp_entry);
	gtk_container_add(GTK_CONTAINER(tempbox), label_K1);
	gtk_container_add(GTK_CONTAINER(box), tempbox);
	gtk_widget_set_sensitive (tempbox, FALSE);
	g_object_set(G_OBJECT(tempbox), "margin-left", 35, NULL);

	dtempbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	GtkWidget *label_dtemp = gtk_label_new(_("Allowed range:"));
	dtemp_entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(dtemp_entry), "500.0");
	gtk_widget_set_tooltip_text(dtemp_entry, _("Stars with temperatures within the reference plus or minus this value will be used"));
	GtkWidget *label_K2 = gtk_label_new("K");
	gtk_container_add(GTK_CONTAINER(dtempbox), label_dtemp);
	gtk_container_add(GTK_CONTAINER(dtempbox), dtemp_entry);
	gtk_container_add(GTK_CONTAINER(dtempbox), label_K2);
	gtk_container_add(GTK_CONTAINER(box), dtempbox);
	gtk_widget_set_sensitive (dtempbox, FALSE);
	g_object_set(G_OBJECT(dtempbox), "margin-left", 35, NULL);

	/* dynamic displays */
	//gtk_widget_set_visible(label_platesolve, has_wcs(gfit));

	if (local_gaia_available()) {
		args.catalogue = CAT_LOCAL_GAIA_ASTRO;
	}
	else if (local_kstars_available()) {
		args.catalogue = CAT_LOCAL_KSTARS;
	}
	else {
		if (!is_online()) {
			gtk_label_set_text(GTK_LABEL(label_cat_name), _("local catalogues not found and offline mode is enabled"));
			args.catalogue = CAT_UNDEF;
		}
		args.catalogue = CAT_GAIADR3;
	}
	if (args.catalogue != CAT_UNDEF) {
		gchar *str = g_strdup_printf("%s (%s)", catalog_to_str(args.catalogue),
				args.catalogue == CAT_GAIADR3 ? "remote" : "local");
		gtk_label_set_text(GTK_LABEL(label_cat_name), str);
		g_free(str);
	}

	// Gather the graphic items
	GtkWidget *content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_box_set_spacing(GTK_BOX(content_area), 15);
	gtk_container_add(GTK_CONTAINER(content_area), box);
	gtk_widget_show_all(GTK_WIDGET(content_area));
}

// the public getter
GtkWidget *get_catmag_dialog() {
	if (!dialog)
		build_the_dialog();
	return dialog;
}

static void on_catmag_response(GtkDialog* self, gint response_id, gpointer user_data) {
	//siril_debug_print("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		gtk_widget_grab_focus(gtk_dialog_get_widget_for_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT));
		gtk_widget_hide(dialog);
		return;
	}

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_filter))) {
		gchar *end;
		const gchar *text = gtk_entry_get_text(GTK_ENTRY(temp_entry));
		double reftemp = g_ascii_strtod(text, &end);
		if (text == end || reftemp < 3000.0 || reftemp > 50000.0) { // range to be defined
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
					_("Reference temperature should be [3000, 50000]"));
			return;
		}

		text = gtk_entry_get_text(GTK_ENTRY(dtemp_entry));
		double dtemp = g_ascii_strtod(text, &end);
		if (text == end || reftemp < 1.0 || reftemp > 20000.0) { // range to be defined
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
					_("Reference temperature should be [1, 20000]"));
			return;
		}

		if (args.catalogue != CAT_LOCAL_KSTARS) {
			args.limit_temperature = TRUE;
			args.limit_BV = FALSE;
			args.refT = reftemp;
			args.dT = dtemp;
		} else {
			args.limit_BV = TRUE;
			args.limit_temperature = FALSE;
			args.refBV = T_to_BV(reftemp);
			args.dBV = (T_to_BV(reftemp - dtemp) - T_to_BV(reftemp + dtemp)) / 2.0;
		}
	} else {
		args.limit_temperature = FALSE;
		args.limit_BV = FALSE;
	}

	args.fit = gfit;
	control_window_switch_to_tab(OUTPUT_LOGS);
	struct catmag_data *worker_args = malloc(sizeof(struct catmag_data));
	memcpy(worker_args, &args, sizeof(struct catmag_data));
	start_in_new_thread(catmag_mono_worker, worker_args);
}
