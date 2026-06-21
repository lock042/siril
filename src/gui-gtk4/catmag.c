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
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/PSF_list.h"

static GtkWidget *dialog = NULL;	// the window, a GtkDialog
static GtkWidget *tempbox = NULL;
static GtkWidget *dtempbox = NULL;
static GtkWidget *check_filter = NULL;
static GtkWidget *temp_entry = NULL;
static GtkWidget *dtemp_entry = NULL;

static struct catmag_data args;

static void on_catmag_response(GtkWindow *self, gint response_id, gpointer user_data);

static GtkWidget *catmag_ok_button = NULL;

/* Phase 14G.4: GtkDialog → GtkWindow.  Bridge per-button "clicked" signals
 * back to the legacy on_catmag_response(self, response_id, …) shape. */
static void catmag_btn_ok_clicked(GtkButton *btn, gpointer ud)    { (void)btn; on_catmag_response(GTK_WINDOW(dialog), GTK_RESPONSE_ACCEPT, ud); }
static void catmag_btn_close_clicked(GtkButton *btn, gpointer ud) { (void)btn; on_catmag_response(GTK_WINDOW(dialog), GTK_RESPONSE_REJECT, ud); }

static void filter_checked(GtkToggleButton *source, gpointer user_data) {
	gtk_widget_set_sensitive(tempbox, siril_toggle_get_active(GTK_WIDGET(check_filter)));
	gtk_widget_set_sensitive(dtempbox, siril_toggle_get_active(GTK_WIDGET(check_filter)));
}

static void build_the_dialog() {
	/* Phase 14G.4: GtkDialog → GtkWindow.  We host the previous content
	 * area in a vertical GtkBox under the window's child slot, and append
	 * an explicit button row at the bottom for Close/OK. */
	dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Calibrate image magnitudes"));
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 200);
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	/* A modal GtkWindow with no transient parent grabs input with no anchor
	 * and can wedge the whole session on Wayland/Xorg. */
	gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(lookup_widget("control_window")));
	gtk_window_set_hide_on_close(GTK_WINDOW(dialog), TRUE);
	g_signal_connect(G_OBJECT(dialog), "close-request", G_CALLBACK(siril_widget_hide_on_delete), NULL);

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

	GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 12); // the main box

	/*if (!has_wcs(gfit)) {
		// unfortunately this is not called when the window is opened but later, so
		// this test doesn't work
		GtkWidget *label_platesolve = gtk_label_new("");
		gtk_label_set_markup(GTK_LABEL(label_platesolve), _("<b>Image must be plate solved first</b>"));
		gtk_box_append(GTK_BOX(box), label_platesolve);
		gtk_widget_set_margin_top(GTK_WIDGET(label_platesolve), 20);
		gtk_widget_set_margin_bottom(GTK_WIDGET(label_platesolve), 10);
	}*/

	GtkWidget *catbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	GtkWidget *label_cat = gtk_label_new(_("Using catalogue:"));
	gtk_box_append(GTK_BOX(catbox), label_cat);
	GtkWidget *label_cat_name = gtk_label_new("");
	gtk_box_append(GTK_BOX(catbox), label_cat_name);
	gtk_box_append(GTK_BOX(box), catbox);
	gtk_widget_set_margin_top(GTK_WIDGET(catbox), 20);
	gtk_widget_set_margin_start(GTK_WIDGET(catbox), 15);
	check_filter = gtk_check_button_new_with_label(_("Filter stars"));
	gtk_box_append(GTK_BOX(box), check_filter);
	g_signal_connect (check_filter, "toggled", G_CALLBACK (filter_checked), NULL);
	gtk_widget_set_margin_start(GTK_WIDGET(check_filter), 12);
	tempbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	GtkWidget *label_temp = gtk_label_new(_("Temperature:"));
	temp_entry = gtk_entry_new();
	gtk_editable_set_text(GTK_EDITABLE(temp_entry), "5550.0");
	gtk_widget_set_tooltip_text(temp_entry, _("Stars with temperatures around this value will be used"));
	GtkWidget *label_K1 = gtk_label_new("K");
	gtk_box_append(GTK_BOX(tempbox), label_temp);
	gtk_box_append(GTK_BOX(tempbox), temp_entry);
	gtk_box_append(GTK_BOX(tempbox), label_K1);
	gtk_box_append(GTK_BOX(box), tempbox);
	gtk_widget_set_sensitive (tempbox, FALSE);
	gtk_widget_set_margin_start(GTK_WIDGET(tempbox), 35);
	dtempbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	GtkWidget *label_dtemp = gtk_label_new(_("Allowed range:"));
	dtemp_entry = gtk_entry_new();
	gtk_editable_set_text(GTK_EDITABLE(dtemp_entry), "500.0");
	gtk_widget_set_tooltip_text(dtemp_entry, _("Stars with temperatures within the reference plus or minus this value will be used"));
	GtkWidget *label_K2 = gtk_label_new("K");
	gtk_box_append(GTK_BOX(dtempbox), label_dtemp);
	gtk_box_append(GTK_BOX(dtempbox), dtemp_entry);
	gtk_box_append(GTK_BOX(dtempbox), label_K2);
	gtk_box_append(GTK_BOX(box), dtempbox);
	gtk_widget_set_sensitive (dtempbox, FALSE);
	gtk_widget_set_margin_start(GTK_WIDGET(dtempbox), 35);
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
	GtkWidget *content_area = gtk_box_new(GTK_ORIENTATION_VERTICAL, 15);
	gtk_widget_set_margin_start(content_area, 12);
	gtk_widget_set_margin_end(content_area, 12);
	gtk_widget_set_margin_top(content_area, 12);
	gtk_widget_set_margin_bottom(content_area, 12);
	gtk_box_append(GTK_BOX(content_area), box);

	/* Action area */
	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	GtkWidget *btn_close = gtk_button_new_with_mnemonic(_("_Close"));
	g_signal_connect(btn_close, "clicked", G_CALLBACK(catmag_btn_close_clicked), NULL);
	gtk_box_append(GTK_BOX(bbox), btn_close);
	catmag_ok_button = gtk_button_new_with_mnemonic(_("_OK"));
	gtk_widget_add_css_class(catmag_ok_button, "suggested-action");
	g_signal_connect(catmag_ok_button, "clicked", G_CALLBACK(catmag_btn_ok_clicked), NULL);
	gtk_box_append(GTK_BOX(bbox), catmag_ok_button);
	gtk_box_append(GTK_BOX(content_area), bbox);
	gtk_window_set_default_widget(GTK_WINDOW(dialog), catmag_ok_button);

	gtk_window_set_child(GTK_WINDOW(dialog), content_area);
}

// the public getter
GtkWidget *get_catmag_dialog() {
	if (!dialog)
		build_the_dialog();
	return dialog;
}

static void on_catmag_response(GtkWindow *self, gint response_id, gpointer user_data) {
	(void)self; (void)user_data;
	if (response_id != GTK_RESPONSE_ACCEPT) {
		if (catmag_ok_button)
			gtk_widget_grab_focus(catmag_ok_button);
		gtk_widget_set_visible(dialog, FALSE);
		reactivate_parent(dialog);
		return;
	}

	if (siril_toggle_get_active(GTK_WIDGET(check_filter))) {
		gchar *end;
		const gchar *text = gtk_editable_get_text(GTK_EDITABLE(temp_entry));
		double reftemp = g_ascii_strtod(text, &end);
		if (text == end || reftemp < 3000.0 || reftemp > 50000.0) { // range to be defined
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
					_("Reference temperature should be [3000, 50000]"));
			return;
		}

		text = gtk_editable_get_text(GTK_EDITABLE(dtemp_entry));
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
