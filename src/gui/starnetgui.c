/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#include <math.h>

#include <string.h>
#include <math.h>
#include <float.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_app_dirs.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/remixer.h"
#include "gui/siril_preview.h"
#include "gui/registration_preview.h"
#include "core/undo.h"

#include "filters/starnet.h"
static gboolean sgui_customstride;
static gboolean sgui_starmask;
static gboolean sgui_upscale;
static gboolean sgui_linear;
static gboolean sgui_follow_on;
static double sgui_starnet_stride;

static void starnet_startup() {
	sgui_linear = FALSE;
	sgui_customstride = FALSE;
	sgui_upscale = FALSE;
	sgui_starmask = TRUE;
	sgui_follow_on = FALSE;
	sgui_starnet_stride = 256.0;
}

/*** callbacks **/

void on_starnet_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_starnet_stride = GTK_SPIN_BUTTON(lookup_widget("spin_starnet_stride"));
	GtkToggleButton *toggle_starnet_stretch = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_stretch"));
	GtkToggleButton *toggle_starnet_followon = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_postremix"));
	GtkToggleButton *toggle_starnet_upsample = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_upsample"));
	GtkToggleButton *toggle_starnet_starmask = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_starmask"));
	GtkToggleButton *toggle_starnet_customstride = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_customstride"));
	GtkLabel *label_starnetinfo = GTK_LABEL(lookup_widget("label_starnetinfo"));

	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("stride_control")), FALSE);
	if (!com.pref.starnet_dir || g_access(com.pref.starnet_dir, R_OK)) {
		gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), FALSE);
		gtk_label_set_text(label_starnetinfo, "Starnet++ installation directory not set.\nMust be configured in Preferences / Miscellaneous.");
	} else
		gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), TRUE);

#ifndef HAVE_LIBTIFF
	gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), FALSE);
	gtk_label_set_text(label_starnetinfo, "Starnet++ unavailable: requires Siril to be compiled with libtiff support.");
#endif
#ifdef HAVE_LIBTIFF
	if (!starnet_executablecheck()) {
		gtk_label_set_text(label_starnetinfo, "No valid Starnet++ executable found in the configured Starnet++ installation directory.\nCheck your Starnet++ installation and Siril configuration.");
		gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), FALSE);
	} else {
		gtk_label_set_text(label_starnetinfo, "");
		gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), TRUE);
	}
#endif
	starnet_startup();

	set_notify_block(TRUE);
	gtk_toggle_button_set_active(toggle_starnet_followon, sgui_follow_on);
	gtk_toggle_button_set_active(toggle_starnet_stretch, sgui_linear);
	gtk_toggle_button_set_active(toggle_starnet_upsample, sgui_upscale);
	gtk_toggle_button_set_active(toggle_starnet_starmask, sgui_starmask);
	gtk_toggle_button_set_active(toggle_starnet_customstride, sgui_customstride);
	gtk_spin_button_set_value(spin_starnet_stride, sgui_starnet_stride);
	set_notify_block(FALSE);
}

void on_starnet_cancel_clicked(GtkButton *button, gpointer user_data) {
	sgui_linear = FALSE;
	sgui_customstride = FALSE;
	sgui_upscale = FALSE;
	sgui_starmask = TRUE;
	siril_close_dialog("starnet_dialog");
}

#ifdef HAVE_LIBTIFF
void on_starnet_execute_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_starnet_stride = GTK_SPIN_BUTTON(lookup_widget("spin_starnet_stride"));
	GtkToggleButton *toggle_starnet_stretch = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_stretch"));
	GtkToggleButton *toggle_starnet_upsample = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_upsample"));
	GtkToggleButton *toggle_starnet_starmask = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_starmask"));
	GtkToggleButton *toggle_starnet_customstride = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_customstride"));

	sgui_customstride = gtk_toggle_button_get_active(toggle_starnet_customstride);
	sgui_starmask = gtk_toggle_button_get_active(toggle_starnet_starmask);
	sgui_upscale = gtk_toggle_button_get_active(toggle_starnet_upsample);
	sgui_linear = gtk_toggle_button_get_active(toggle_starnet_stretch);
	sgui_starnet_stride = gtk_spin_button_get_value(spin_starnet_stride);
	if (fmod(sgui_starnet_stride, 2))
		sgui_starnet_stride = 2 * round(sgui_starnet_stride / 2);
	if (sgui_starnet_stride < 2.0)
		sgui_starnet_stride = 2.0;
	if (sgui_starnet_stride > 256.0)
		sgui_starnet_stride = 256.0;
	starnet_data *starnet_args;
	starnet_args = malloc(sizeof(starnet_data));
	memset(starnet_args->stride, 0, sizeof(starnet_args->stride));
	starnet_args->customstride = sgui_customstride;
	starnet_args->upscale = sgui_upscale;
	starnet_args->linear = sgui_linear;
	starnet_args->starmask = sgui_starmask;
	sprintf(starnet_args->stride, "%d", (int) sgui_starnet_stride);
	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);
	starnet_args->follow_on = sgui_follow_on;
	start_in_new_thread(do_starnet, starnet_args);
	siril_close_dialog("starnet_dialog");
}
#endif

/*** adjusters **/
void on_spin_starnet_stride_value_changed(GtkSpinButton *button, gpointer user_data) {
	sgui_starnet_stride = gtk_spin_button_get_value(button);
	if (fmod(sgui_starnet_stride, 2))
		sgui_starnet_stride = 2 * round(sgui_starnet_stride / 2);
	if (sgui_starnet_stride < 2.0)
		sgui_starnet_stride = 2.0;
	if (sgui_starnet_stride > 256.0)
		sgui_starnet_stride = 256.0;
}

void on_toggle_starnet_customstride_toggled(GtkToggleButton *button, gpointer user_data) {
	sgui_customstride = gtk_toggle_button_get_active(button);
	if (sgui_customstride) {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("stride_control")), TRUE);
	} else {
		gtk_widget_set_visible(GTK_WIDGET(lookup_widget("stride_control")), FALSE);
	}
}

void on_toggle_starnet_stretch_toggled(GtkToggleButton *button, gpointer user_data) {
	sgui_linear = gtk_toggle_button_get_active(button);
}
void on_toggle_starnet_postremix_toggled(GtkToggleButton *button, gpointer user_data) {
	sgui_follow_on = gtk_toggle_button_get_active(button);
	if (sgui_follow_on && !sgui_starmask) {
		sgui_starmask = TRUE;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_starmask")), sgui_starmask);
	}
}

void on_toggle_starnet_upsample_toggled(GtkToggleButton *button, gpointer user_data) {
	sgui_upscale = gtk_toggle_button_get_active(button);
}

void on_toggle_starnet_starmask_toggled(GtkToggleButton *button, gpointer user_data) {
	sgui_starmask = gtk_toggle_button_get_active(button);
		if (sgui_follow_on && !sgui_starmask) {
		sgui_follow_on = FALSE;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_postremix")), sgui_follow_on);
	}


}
