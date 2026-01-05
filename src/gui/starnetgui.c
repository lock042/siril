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

#include <math.h>

#include <string.h>
#include <float.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "filters/starnet.h"

static gboolean sgui_customstride;
static gboolean sgui_starmask;
static gboolean sgui_upscale;
static gboolean sgui_linear;
static gboolean sgui_follow_on;
static int sgui_starnet_stride;
static int sgui_next_seq_index;

static void starnet_startup() {
	sgui_linear = TRUE;
	sgui_customstride = FALSE;
	sgui_upscale = FALSE;
	sgui_starmask = TRUE;
	sgui_follow_on = FALSE;
	sgui_starnet_stride = 256;
	sgui_next_seq_index = 0;
}

/*** callbacks **/

void on_starnet_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_starnet_stride = GTK_SPIN_BUTTON(lookup_widget("spin_starnet_stride"));
	GtkToggleButton *toggle_starnet_stretch = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_stretch"));
	GtkToggleButton *toggle_starnet_followon = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_postremix"));
	GtkToggleButton *toggle_starnet_upsample = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_upsample"));
	GtkToggleButton *toggle_starnet_starmask = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_starmask"));
	GtkToggleButton *toggle_starnet_customstride = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_customstride"));
	GtkToggleButton *toggle_starnet_sequence = GTK_TOGGLE_BUTTON(lookup_widget("starnet_sequence_toggle"));
	GtkWidget *starnet_next_sequence_controls = lookup_widget("starnet_next_sequence_controls");
	GtkLabel *label_starnetinfo = GTK_LABEL(lookup_widget("label_starnetinfo"));
	gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), TRUE);

	gtk_widget_set_visible(GTK_WIDGET(lookup_widget("stride_control")), FALSE);
	if (!com.pref.starnet_exe || g_access(com.pref.starnet_exe, R_OK)) {
		gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), FALSE);
		gtk_label_set_text(label_starnetinfo, "StarNet installation directory not set.\nMust be configured in Preferences / Miscellaneous.");
	} else
		gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), TRUE);

#ifndef HAVE_LIBTIFF
	gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), FALSE);
	gtk_label_set_text(label_starnetinfo, "StarNet unavailable: requires Siril to be compiled with libtiff support.");
#endif
#ifdef HAVE_LIBTIFF
	gchar *temp = NULL, *winext = NULL, *starnetcommand = NULL;
	gchar *bothtext = g_strdup(_("Valid StarNet v1 mono and RGB executables found in the configured StarNet installation directory.\nStarNet can process mono and RGB images."));
	switch (starnet_executablecheck(com.pref.starnet_exe)) {
		case NIL:
			gtk_label_set_text(label_starnetinfo, _("No valid StarNet executable found in the configured StarNet installation directory.\nCheck your StarNet installation and Siril configuration."));
			gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), FALSE);
			break;
		case TORCH:
			gtk_label_set_text(label_starnetinfo, _("Valid StarNet v2-Torch executable found in the configured StarNet installation directory.\nStarNet can process mono and RGB images."));
			break;
		case V2:
			gtk_label_set_text(label_starnetinfo, _("Valid StarNet v2 executable found in the configured StarNet installation directory.\nStarNet can process mono and RGB images."));
			break;
		case V1MONO:
			temp = g_path_get_dirname(com.pref.starnet_exe);
			winext = g_str_has_suffix(com.pref.starnet_exe, ".exe") ? g_strdup(".exe") : g_strdup("\0");
			starnetcommand = g_strdup_printf("%s/rgb_starnet++%s", temp, winext);
			if (starnet_executablecheck(starnetcommand) == V1RGB)
				gtk_label_set_text(label_starnetinfo, bothtext);
			else {
				gtk_label_set_text(label_starnetinfo, _("Only the StarNet v1 mono executable was found in the configured StarNet installation directory.\nStarNet can process mono images only."));
				if (gfit->naxes[2] == 3)
					gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), FALSE);
			}
			g_free(temp);
			temp = NULL;
			g_free(winext);
			winext = NULL;
			g_free(starnetcommand);
			starnetcommand = NULL;
			break;
		case V1RGB:
			temp = g_path_get_dirname(com.pref.starnet_exe);
			winext = g_str_has_suffix(com.pref.starnet_exe, ".exe") ? g_strdup(".exe") : g_strdup("\0");
			starnetcommand = g_strdup_printf("%s/mono_starnet++%s", temp, winext);
			if (starnet_executablecheck(starnetcommand) == V1MONO)
				gtk_label_set_text(label_starnetinfo, bothtext);
			else {
				gtk_label_set_text(label_starnetinfo, _("Only the StarNet v1 RGB executable was found in the configured StarNet installation directory.\nStarNet can process RGB images only."));
				if (gfit->naxes[2] == 1)
					gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("starnet_apply")), FALSE);
			}
			g_free(temp);
			temp = NULL;
			g_free(winext);
			winext = NULL;
			g_free(starnetcommand);
			starnetcommand = NULL;
			break;
	}
	g_free(bothtext);
#endif
	starnet_startup();

	set_notify_block(TRUE);
	gtk_toggle_button_set_active(toggle_starnet_followon, sgui_follow_on);
	gtk_toggle_button_set_active(toggle_starnet_sequence, FALSE);
	gtk_toggle_button_set_active(toggle_starnet_stretch, sgui_linear);
	gtk_toggle_button_set_active(toggle_starnet_upsample, sgui_upscale);
	gtk_toggle_button_set_active(toggle_starnet_starmask, sgui_starmask);
	gtk_toggle_button_set_active(toggle_starnet_customstride, sgui_customstride);
	gtk_spin_button_set_value(spin_starnet_stride, sgui_starnet_stride);
	gtk_widget_set_sensitive(starnet_next_sequence_controls, sgui_starmask && gtk_toggle_button_get_active(toggle_starnet_sequence));
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
	if (!check_ok_if_cfa())
		return;

	GtkSpinButton *spin_starnet_stride = GTK_SPIN_BUTTON(lookup_widget("spin_starnet_stride"));
	GtkToggleButton *toggle_starnet_stretch = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_stretch"));
	GtkToggleButton *toggle_starnet_sequence = GTK_TOGGLE_BUTTON(lookup_widget("starnet_sequence_toggle"));
	GtkToggleButton *toggle_starnet_upsample = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_upsample"));
	GtkToggleButton *toggle_starnet_starmask = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_starmask"));
	GtkToggleButton *toggle_starnet_customstride = GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_customstride"));
	GtkComboBox *combo_starnet_next_sequence = GTK_COMBO_BOX(lookup_widget("combo_starnet_next_sequence"));

	sgui_customstride = gtk_toggle_button_get_active(toggle_starnet_customstride);
	sgui_starmask = gtk_toggle_button_get_active(toggle_starnet_starmask);
	sgui_upscale = gtk_toggle_button_get_active(toggle_starnet_upsample);
	sgui_linear = gtk_toggle_button_get_active(toggle_starnet_stretch);
	sgui_starnet_stride = (int) gtk_spin_button_get_value(spin_starnet_stride);
	if (sgui_starnet_stride % 2)
		sgui_starnet_stride++;
	if (sgui_starnet_stride < 2)
		sgui_starnet_stride = 2;
	if (sgui_starnet_stride > 256)
		sgui_starnet_stride = 256;

	// Allocate parameters using the allocator
	starnet_data *starnet_params = new_starnet_args();
	if (!starnet_params) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		return;
	}

	starnet_params->starnet_fit = gfit;
	starnet_params->imgnumber = -1;
	starnet_params->customstride = sgui_customstride;
	starnet_params->upscale = sgui_upscale;
	starnet_params->linear = sgui_linear;
	starnet_params->starmask = sgui_starmask;
	starnet_params->stride = g_strdup_printf("%d", sgui_starnet_stride);
	starnet_params->follow_on = sgui_follow_on;

	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);

	if (gtk_toggle_button_get_active(toggle_starnet_sequence) == FALSE) {
		if (single_image_is_loaded()) {
			// Save backup for undo before processing
			copy_gfit_to_backup();

			// Allocate generic_img_args
			struct generic_img_args *args = calloc(1, sizeof(struct generic_img_args));
			if (!args) {
				PRINT_ALLOC_ERR;
				free_starnet_args(starnet_params);
				set_cursor_waiting(FALSE);
				return;
			}

			// Set up generic_img_args
			args->fit = gfit;
			args->mem_ratio = 3.0f;
			args->image_hook = starnet_single_image_hook;
			args->idle_function = starnet_single_image_idle;
			args->description = _("StarNet");
			args->verbose = TRUE;
			args->user = starnet_params;
			args->max_threads = com.max_thread;
			args->for_preview = FALSE;
			args->for_roi = FALSE;

			if (!start_in_new_thread(generic_image_worker, args)) {
				free_generic_img_args(args);
			}
			siril_close_dialog("starnet_dialog");
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Not in single image mode"),
				_("Unable to apply StarNet to a single image as no single image is loaded. Did you mean to apply to sequence?"));
			free_starnet_args(starnet_params);
		}
		set_cursor_waiting(FALSE);
	} else {
		if (sequence_is_loaded()) {
			starnet_params->follow_on = FALSE;
			struct multi_output_data *multi_args = calloc(1, sizeof(struct multi_output_data));
			multi_args->user_data = (gpointer) starnet_params;
			starnet_params->multi_args = multi_args;
			multi_args->seq = &com.seq;
			multi_args->n = sgui_starmask ? 2 : 1;
			multi_args->new_seq_index = gtk_combo_box_get_active(combo_starnet_next_sequence);
			multi_args->prefixes = calloc(multi_args->n, sizeof(char*));
			multi_args->prefixes[0] = g_strdup("starless_");
			if (sgui_starmask) {
				multi_args->prefixes[1] = g_strdup("starmask_");
			}
			multi_args->seqEntry = strdup(multi_args->prefixes[0]);
			apply_starnet_to_sequence(multi_args);
			siril_close_dialog("starnet_dialog");
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("No sequence loaded"),
				_("Unable to apply StarNet to a sequence as no sequence is loaded. Did you mean to uncheck the apply to sequence option?"));
			free_starnet_args(starnet_params);
		}
		set_cursor_waiting(FALSE);
	}
}
#endif

/*** adjusters **/
void on_spin_starnet_stride_value_changed(GtkSpinButton *button, gpointer user_data) {
	sgui_starnet_stride = (int) gtk_spin_button_get_value(button);
	if (sgui_starnet_stride % 2)
		sgui_starnet_stride++;
	if (sgui_starnet_stride < 2)
		sgui_starnet_stride = 2;
	if (sgui_starnet_stride > 256)
		sgui_starnet_stride = 256;
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
	GtkToggleButton *toggle_starnet_sequence = GTK_TOGGLE_BUTTON(lookup_widget("starnet_sequence_toggle"));
	sgui_starmask = gtk_toggle_button_get_active(button);
		if (sgui_follow_on && !sgui_starmask) {
		sgui_follow_on = FALSE;
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_starnet_postremix")), sgui_follow_on);
	}
	gtk_widget_set_sensitive(lookup_widget("starnet_next_sequence_controls"), sgui_starmask && gtk_toggle_button_get_active(toggle_starnet_sequence));
}

void on_starnet_sequence_toggle_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("starnet_next_sequence_controls"), sgui_starmask && gtk_toggle_button_get_active(button));
}

void on_combo_starnet_next_sequence_changed(GtkComboBox *combo, gpointer user_data) {
	sgui_next_seq_index = gtk_combo_box_get_active(combo);
}
