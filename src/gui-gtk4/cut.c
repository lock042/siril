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

#include "algos/demosaicing.h"
#include "algos/PSF.h"
#include "core/cut.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/utils.h"
#include "io/sequence.h"

gboolean reset_cut_gui_filedependent(gpointer user_data) { // Separated out to avoid having to repeat too much after opening a new file
	GtkWidget *colorbutton = (GtkWidget*) lookup_widget("cut_radio_color");
	GtkWidget *cfabutton = (GtkWidget*) lookup_widget("cut_cfa");
	gtk_widget_set_sensitive(colorbutton, (gfit->naxes[2] == 3));
	sensor_pattern pattern = get_cfa_pattern_index_from_string(gfit->keywords.bayer_pattern);
	gboolean cfa_disabled = (gfit->naxes[2] > 1 || pattern < BAYER_FILTER_MIN || pattern > BAYER_FILTER_MAX);
	gtk_widget_set_sensitive(cfabutton, !cfa_disabled);
	GtkToggleButton* as = (GtkToggleButton*) lookup_widget("cut_dist_pref_as");
	siril_toggle_set_active(GTK_WIDGET(as), gfit->keywords.wcsdata.pltsolvd);
	return FALSE;
}

gboolean reset_cut_gui(gpointer user_data) {
	GtkToggleButton *radio_mono = (GtkToggleButton*) lookup_widget("cut_radio_mono");
	siril_toggle_set_active(GTK_WIDGET(radio_mono), TRUE);
	GtkToggleButton *save_dat = (GtkToggleButton*) lookup_widget("cut_save_checkbutton");
	siril_toggle_set_active(GTK_WIDGET(save_dat), FALSE);
	GtkToggleButton *save_png = (GtkToggleButton*) lookup_widget("cut_save_png");
	siril_toggle_set_active(GTK_WIDGET(save_png), FALSE);
	GtkSpinButton *cut_startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
	GtkSpinButton *cut_starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
	GtkSpinButton *cut_finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
	GtkSpinButton *cut_finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
	gtk_spin_button_set_value(cut_startx, -1);
	gtk_spin_button_set_value(cut_starty, -1);
	gtk_spin_button_set_value(cut_finishx, -1);
	gtk_spin_button_set_value(cut_finishy, -1);
	GtkSpinButton *cut_width = (GtkSpinButton*) lookup_widget("cut_spin_width");
	GtkSpinButton *cut_wn1 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber1");
	GtkSpinButton *cut_wn2 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber2");
	gtk_spin_button_set_value(cut_width, 1);
	gtk_spin_button_set_value(cut_wn1, -1);
	gtk_spin_button_set_value(cut_wn2, -1);
	GtkLabel *wn1x = (GtkLabel*) lookup_widget("label_wn1_x");
	GtkLabel *wn1y = (GtkLabel*) lookup_widget("label_wn1_y");
	GtkLabel *wn2x = (GtkLabel*) lookup_widget("label_wn2_x");
	GtkLabel *wn2y = (GtkLabel*) lookup_widget("label_wn2_y");
	gtk_label_set_text(wn1x, "");
	gtk_label_set_text(wn1y, "");
	gtk_label_set_text(wn2x, "");
	gtk_label_set_text(wn2y, "");
	GtkSpinButton *bgpoly = GTK_SPIN_BUTTON(lookup_widget("spin_spectro_bgpoly"));
	gtk_spin_button_set_value(bgpoly, 3);
	GtkCheckButton *plot_spectro_bg = GTK_CHECK_BUTTON(lookup_widget("cut_spectro_plot_bg"));
	siril_toggle_set_active(GTK_WIDGET(plot_spectro_bg), FALSE);
	GtkEntry *title = (GtkEntry*) lookup_widget("cut_title");
	gtk_editable_set_text(GTK_EDITABLE(title), _("Intensity Profile"));
	reset_cut_gui_filedependent(NULL);
	return FALSE;
}



void measure_line(fits *fit, point start, point finish, gboolean pref_as) {
	int deg = -1;
	static const gchar *label_selection[] = { "labelselection_red", "labelselection_green", "labelselection_blue", "labelselection_rgb" };
	static gchar measurement_buffer[256] = { 0 };
	point delta = { finish.x - start.x, finish.y - start.y };
	double pixdist = sqrt(delta.x * delta.x + delta.y * delta.y);
	if (pixdist == 0.) {
		measurement_buffer[0] = '\0';
	} else {
		double conversionfactor = get_conversion_factor(fit);
		if (conversionfactor != -DBL_MAX && pref_as) {
			double asdist = pixdist * conversionfactor;
			if (asdist < 60.0) {
				g_sprintf(measurement_buffer, _("Measurement: %.1f\""), asdist);
			} else {
				int min = (int) asdist / 60;
				double sec = asdist - (min * 60);
				if (asdist < 3600) {
					g_sprintf(measurement_buffer, _("Measurement: %d\' %.1f\""), min, sec);
				} else {
					deg = (int) asdist / 3600;
					min -= (deg * 60);
					g_sprintf(measurement_buffer, _("Measurement: %dº %d\' %.0f\""), deg, min, sec);
				}
			}
		} else {
			g_sprintf(measurement_buffer, _("Measurement: %.1f px"), pixdist);
		}
		gtk_label_set_text(GTK_LABEL(lookup_widget(label_selection[gui.cvport])), measurement_buffer);
		if (deg > 10.0) {
			control_window_switch_to_tab(OUTPUT_LOGS);
			siril_log_warning(_("Warning: angular measurement > 10º. Error is > 1%\n"));
		}
	}
}

static void update_spectro_coords() {
	GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
	GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
	GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
	GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
	gtk_spin_button_set_value(startx, gui.cut.cut_start.x);
	gtk_spin_button_set_value(starty, gui.cut.cut_start.y);
	gtk_spin_button_set_value(finishx, gui.cut.cut_end.x);
	gtk_spin_button_set_value(finishy, gui.cut.cut_end.y);
}

//// GUI callbacks ////
void on_cut_sequence_apply_from_gui() {
	GtkToggleButton* cut_color = (GtkToggleButton*)lookup_widget("cut_radio_color");
	cut_struct *arg = calloc(1, sizeof(cut_struct));
	memcpy(arg, &gui.cut, sizeof(cut_struct));
	arg->title = g_strdup(gui.cut.title);
	arg->user_title = g_strdup(gui.cut.user_title);
	arg->filename = g_strdup(gui.cut.filename);
	arg->save_png_too = FALSE;
	arg->fit = NULL;
	arg->seq = &com.seq;
	if (siril_toggle_get_active(GTK_WIDGET(cut_color)))
		arg->mode = CUT_COLOR;
	else
		arg->mode = CUT_MONO;
	arg->display_graph = FALSE;
	arg->cut_measure = FALSE;
	control_window_switch_to_tab(OUTPUT_LOGS);
	// Check args are cromulent
	if (cut_struct_is_valid(arg))
		apply_cut_to_sequence(arg);
	else
		free_cut_args(arg);
}

void on_cut_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry* entry = (GtkEntry*) lookup_widget("cut_title");
	if (gui.cut.user_title)
		g_free(gui.cut.user_title);
	gui.cut.user_title = g_strdup(gtk_editable_get_text(GTK_EDITABLE(entry)));
	GtkToggleButton* cut_color = (GtkToggleButton*)lookup_widget("cut_radio_color");
	if (siril_toggle_get_active(GTK_WIDGET(cut_color))) {
		gui.cut.mode = CUT_COLOR;
	} else {
		gui.cut.mode = CUT_MONO;
	}
	GtkToggleButton* apply_to_sequence = (GtkToggleButton*)lookup_widget("cut_apply_to_sequence");
	gui.cut.vport = gui.cvport;
	if (siril_toggle_get_active(GTK_WIDGET(apply_to_sequence))) {
		if (sequence_is_loaded())
			on_cut_sequence_apply_from_gui();
		else
			siril_message_dialog(GTK_MESSAGE_ERROR,
					_("No sequence is loaded"),
					_("The Apply to sequence option is checked, but no sequence is loaded."));

	} else {
		gui.cut.fit = gfit;
		gui.cut.seq = NULL;
		gui.cut.display_graph = TRUE;
		// We have to pass a dynamically allocated copy of gui.cut
		// otherwise start_in_new_thread() can try to free gui.cut
		// if the processing thread is already running
		cut_struct *p = malloc(sizeof(cut_struct));
		memcpy(p, &gui.cut, sizeof(cut_struct));
		if (p->tri) {
			siril_log_debug("Tri-profile\n");
			if (!start_in_new_thread(tri_cut, p))
				free(p);
		} else if (p->cfa) {
			siril_log_debug("CFA profiling\n");
			if (!start_in_new_thread(cfa_cut, p))
				free(p);
		} else {
			siril_log_debug("Single profile\n");
			if (!start_in_new_thread(cut_profile, p))
				free(p);
		}
	}
}

void on_cut_close_button_clicked(GtkButton *button, gpointer user_data) {
	/* Drive the win.cut stateful action to FALSE so the toolbar toggle
	 * tracks the dialog being dismissed.  cut_state handles the close
	 * and the mouse_status reset, and the GtkToggleButton bound to
	 * action-name="win.cut" follows the state automatically.  Setting
	 * the state explicitly (rather than activating) avoids opening the
	 * dialog if the action ever got out of sync with the toggle. */
	GtkRoot *root = gtk_widget_get_root(GTK_WIDGET(button));
	if (root && G_IS_ACTION_GROUP(root))
		g_action_group_change_action_state(G_ACTION_GROUP(root), "cut",
			g_variant_new_boolean(FALSE));
	else
		siril_close_dialog("cut_dialog");
}

void match_adjustments_to_gfit() {
	GtkAdjustment *sxa = (GtkAdjustment*) lookup_adjustment("adj_cut_xstart");
	GtkAdjustment *fxa = (GtkAdjustment*) lookup_adjustment("adj_cut_xfinish");
	GtkAdjustment *sya = (GtkAdjustment*) lookup_adjustment("adj_cut_ystart");
	GtkAdjustment *fya = (GtkAdjustment*) lookup_adjustment("adj_cut_yfinish");
	gtk_adjustment_set_upper(sxa, gui.cut.fit->rx);
	gtk_adjustment_set_upper(fxa, gui.cut.fit->rx);
	gtk_adjustment_set_upper(sya, gui.cut.fit->ry);
	gtk_adjustment_set_upper(fya, gui.cut.fit->ry);
}

void on_cut_manual_coords_button_clicked(GtkButton* button, gpointer user_data) {
	match_adjustments_to_gfit();
	g_signal_handlers_block_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
	GtkWidget *cut_coords_dialog = lookup_widget("cut_coords_dialog");
	if (gui.cut.cut_start.x != -1) // If there is a cut line already made, show the
							   // endpoint coordinates in the dialog
		update_spectro_coords();
	if (!gtk_widget_is_visible(cut_coords_dialog))
		siril_open_dialog("cut_coords_dialog");
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
}

void on_cut_spectroscopic_button_clicked(GtkButton* button, gpointer user_data) {
	g_signal_handlers_block_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
	GtkWidget *cut_spectroscopy_dialog = lookup_widget("cut_spectroscopy_dialog");
	if (!gtk_widget_is_visible(cut_spectroscopy_dialog))
		siril_open_dialog("cut_spectroscopy_dialog");
	mouse_status = MOUSE_ACTION_NONE;

}

void update_spectro_labels() {
	GtkWidget* monobutton = lookup_widget("cut_radio_mono");
	GtkWidget* colorbutton = lookup_widget("cut_radio_color");
	GtkWidget* tributton = lookup_widget("cut_tri_cut");
	GtkWidget* cfabutton = lookup_widget("cut_cfa");
	GtkWidget* cut_offset_label = lookup_widget("cut_offset_label");
	GtkWidget* pixels = lookup_widget("cut_dist_pref_px");
	GtkWidget* arcsec = lookup_widget("cut_dist_pref_as");
	if (spectroscopy_selections_are_valid(&gui.cut)) {
		/* cut_radio_mono and cut_tri_cut are GtkCheckButton in cut_dialog.ui.
		 * In GTK4 GtkCheckButton no longer derives from GtkButton (unlike
		 * GTK3) so GTK_BUTTON() asserts and gtk_button_set_label crashes.
		 * Use the GtkCheckButton API instead. */
		gtk_check_button_set_label(GTK_CHECK_BUTTON(monobutton), _("Spectroscopic"));
		gtk_widget_set_tooltip_text(monobutton, _("Reduces a spectrum without background removal. This is suitable when the entire image represents a calibrated spectrum"));
		gtk_check_button_set_label(GTK_CHECK_BUTTON(tributton), _("Spectro w/ bg removal"));
		gtk_label_set_text(GTK_LABEL(cut_offset_label), _("Spectro bg offset (px)"));
		gtk_widget_set_tooltip_text(tributton, _("Reduces a spectrum with background removal. This is suitable when background removal is required: the background is computed along parallel lines equidistant from the central spectral profile line"));
		gtk_widget_set_visible(colorbutton, FALSE);
		gtk_widget_set_visible(cfabutton, FALSE);
		gtk_widget_set_visible(pixels, FALSE);
		gtk_widget_set_visible(arcsec, FALSE);
	} else {
		gtk_check_button_set_label(GTK_CHECK_BUTTON(monobutton), _("Mono"));
		gtk_widget_set_tooltip_text(monobutton, _("Generates a single luminance profile along the profile line"));
		gtk_check_button_set_label(GTK_CHECK_BUTTON(tributton), _("Tri-profile (mono)"));
		gtk_widget_set_tooltip_text(tributton, _("Generates 3 parallel intensity profiles separated by a given number of pixels. Tri-profiles always plot luminance along each profile"));
		gtk_label_set_text(GTK_LABEL(cut_offset_label), _("Tri-profile offset (px)"));
		gtk_widget_set_visible(colorbutton, TRUE);
		gtk_widget_set_visible(cfabutton, TRUE);
		gtk_widget_set_visible(pixels, TRUE);
		gtk_widget_set_visible(arcsec, TRUE);
	}
}

void on_cut_dialog_show(GtkWindow *dialog, gpointer user_data) {
	GtkWidget* colorbutton = lookup_widget("cut_radio_color");
	GtkWidget* cfabutton = lookup_widget("cut_cfa");
	GtkCheckButton *plot_bg = GTK_CHECK_BUTTON(lookup_widget("cut_spectro_plot_bg"));
	GtkToggleButton* seqbutton = (GtkToggleButton*) lookup_widget("cut_apply_to_sequence");
	GtkToggleButton* pngbutton = (GtkToggleButton*) lookup_widget("cut_save_png");
	gtk_widget_set_sensitive(colorbutton, (gfit->naxes[2] == 3));
	sensor_pattern pattern = get_cfa_pattern_index_from_string(gfit->keywords.bayer_pattern);
	gboolean cfa_disabled = ((gfit->naxes[2] > 1) || ((!(pattern == BAYER_FILTER_RGGB || pattern == BAYER_FILTER_GRBG || pattern == BAYER_FILTER_BGGR || pattern == BAYER_FILTER_GBRG))));
	gtk_widget_set_sensitive(cfabutton, !cfa_disabled);
	if (siril_toggle_get_active(GTK_WIDGET(seqbutton)))
		siril_toggle_set_active(GTK_WIDGET(pngbutton), TRUE);
	GtkToggleButton *save_dat = (GtkToggleButton*) lookup_widget("cut_save_checkbutton");
	siril_toggle_set_active(GTK_WIDGET(save_dat), FALSE);
	gui.cut.save_dat = siril_toggle_get_active(GTK_WIDGET(save_dat));
	update_spectro_labels();
	gui.cut.plot_spectro_bg = siril_toggle_get_active(GTK_WIDGET(plot_bg));
}

void on_cut_spectro_cancel_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("cut_spectroscopy_dialog");
}

void on_cut_coords_cancel_button_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
	GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
	GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
	GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
	siril_close_dialog("cut_coords_dialog");
	// Reset coords widgets to match the values in the struct gui.cut
	// Do this after the dialog is closed in order to avoid potential
	// momentary flickering of the widget values
	gtk_spin_button_set_value(startx, gui.cut.cut_start.x);
	gtk_spin_button_set_value(starty, gui.cut.cut_start.y);
	gtk_spin_button_set_value(finishx, gui.cut.cut_end.x);
	gtk_spin_button_set_value(finishy, gui.cut.cut_end.y);
}

void on_cut_coords_dialog_hide(GtkWindow *window, gpointer user_data) {
	siril_open_dialog("cut_dialog");
	mouse_status = MOUSE_ACTION_CUT_SELECT;
	g_signal_handlers_unblock_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
}

void on_cut_spectroscopy_dialog_hide(GtkWindow *window, gpointer user_data) {
	siril_open_dialog("cut_dialog");
	mouse_status = MOUSE_ACTION_CUT_SELECT;
	g_signal_handlers_unblock_by_func(GTK_WINDOW(lookup_widget("cut_dialog")), on_cut_close_button_clicked, NULL);
}

void on_cut_coords_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
	GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
	GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
	GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
	int sx = (int) gtk_spin_button_get_value(startx);
	int sy = (int) gtk_spin_button_get_value(starty);
	int fx = (int) gtk_spin_button_get_value(finishx);
	int fy = (int) gtk_spin_button_get_value(finishy);
	siril_log_debug("start (%d, %d) finish (%d, %d)\n", sx, sy, fx, fy);
	// Update the struct gui.cut with the entered values
	// This is all done at once when Apply is clicked rather than individually in the
	// GtkSpinButton callbacks in order to avoid the endpoints jumping about and the
	// line looking like it's in the wrong place.
	gui.cut.cut_start.x = sx;
	gui.cut.cut_start.y = sy;
	gui.cut.cut_end.x = fx;
	gui.cut.cut_end.y = fy;
	measure_line(gfit, gui.cut.cut_start, gui.cut.cut_end, gui.cut.pref_as);
	redraw(REDRAW_OVERLAY);
	siril_close_dialog("cut_coords_dialog");
}

void on_cut_wavenumber1_clicked(GtkButton *button, gpointer user_data) {
	set_cursor("crosshair");
	mouse_status = MOUSE_ACTION_CUT_WN1;
}

void on_cut_wavenumber2_clicked(GtkButton *button, gpointer user_data) {
	set_cursor("crosshair");
	mouse_status = MOUSE_ACTION_CUT_WN2;
}

void on_cut_spin_width_value_changed(GtkSpinButton *button, gpointer user_data) {
	int n = (int) gtk_spin_button_get_value(button);
	if (!(n % 2)) {
		n++;
		gtk_spin_button_set_value(button, n);
	}
	gui.cut.width = n;
}

void on_cut_spectro_apply_button_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton* cut_spin_wavenumber1 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber1");
	GtkSpinButton* cut_spin_wavenumber2 = (GtkSpinButton*) lookup_widget("cut_spin_wavenumber2");
	GtkSpinButton* cut_spin_width = (GtkSpinButton*) lookup_widget("cut_spin_width");
	gui.cut.wavenumber1 = gtk_spin_button_get_value(cut_spin_wavenumber1);
	gui.cut.wavenumber2 = gtk_spin_button_get_value(cut_spin_wavenumber2);
	gui.cut.width = (int) gtk_spin_button_get_value(cut_spin_width);
	siril_close_dialog("cut_spectroscopy_dialog");
}

void on_select_from_star_clicked(GtkButton *button, gpointer user_data) {
	psf_star *result = NULL;
	int layer = 0; // Detect stars in layer 0 (mono or red) as this will always be present
	const gchar *caller = gtk_buildable_get_buildable_id(GTK_BUILDABLE(button));
	gboolean is_start = !g_strcmp0("start_select_from_star", caller);

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		psf_error error = PSF_NO_ERR;
		result = psf_get_minimisation(gui.cut.fit, layer, &com.selection, FALSE, FALSE, NULL, FALSE, com.pref.starfinder_conf.profile, &error);
		set_cursor_waiting(FALSE);
		if (result && error == PSF_NO_ERR) {
			if (is_start) {
				gui.cut.cut_start.x = result->x0 + com.selection.x;
				gui.cut.cut_start.y = com.selection.y + com.selection.h - result->y0;
				GtkSpinButton* startx = (GtkSpinButton*) lookup_widget("cut_xstart_spin");
				GtkSpinButton* starty = (GtkSpinButton*) lookup_widget("cut_ystart_spin");
				gtk_spin_button_set_value(startx, gui.cut.cut_start.x);
				gtk_spin_button_set_value(starty, gui.cut.cut_start.y);
			} else {
				gui.cut.cut_end.x = result->x0 + com.selection.x;
				gui.cut.cut_end.y = com.selection.y + com.selection.h - result->y0;
				GtkSpinButton* finishx = (GtkSpinButton*) lookup_widget("cut_xfinish_spin");
				GtkSpinButton* finishy = (GtkSpinButton*) lookup_widget("cut_yfinish_spin");
				gtk_spin_button_set_value(finishx, gui.cut.cut_end.x);
				gtk_spin_button_set_value(finishy, gui.cut.cut_end.y);
			}
			redraw(REDRAW_OVERLAY);
			measure_line(gfit, gui.cut.cut_start, gui.cut.cut_end, gui.cut.pref_as);
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR,
						_("No star detected"),
						_("Siril cannot set the star coordinate as no star has been detected in the selection"));
		}
		free_psf(result);
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR,
					_("No selection"),
					_("Siril cannot set the star coordinate as no selection is made"));
	}
}

void on_cut_save_png_toggled(GtkCheckButton *button, gpointer user_data) {
	GtkToggleButton* seqbutton = (GtkToggleButton*) lookup_widget("cut_apply_to_sequence");
	gui.cut.save_png_too = siril_toggle_get_active(GTK_WIDGET(button));
	if (!gui.cut.save_png_too)
		siril_toggle_set_active(GTK_WIDGET(seqbutton), FALSE);
}

void on_cut_apply_to_sequence_toggled(GtkCheckButton *button, gpointer user_data) {
	// The sequence mode will always save PNGs, so we set the option at the same time.
	// It's not really necessary as the structure member isn't used, but it keeps
	// things consistent for the user.
	GtkToggleButton *pngbutton = (GtkToggleButton*) lookup_widget("cut_save_png");
	if (siril_toggle_get_active(GTK_WIDGET(button))) {
		siril_toggle_set_active(GTK_WIDGET(pngbutton), TRUE);
		gui.cut.save_png_too = TRUE;
	}
}

void on_cut_tricut_step_value_changed(GtkSpinButton *button, gpointer user_data) {
	int n = (int) gtk_spin_button_get_value(button);
	if (!(n % 2)) {
		n++;
		gtk_spin_button_set_value(button, n);
	}
	gui.cut.step = n;
	redraw(REDRAW_OVERLAY);
}

void on_cut_tri_cut_toggled(GtkCheckButton *button, gpointer user_data) {
	gui.cut.tri = siril_toggle_get_active(GTK_WIDGET(button));
	if (gui.cut.tri) {
		gui.cut.cfa = FALSE;
	}
	gtk_widget_set_sensitive(GTK_WIDGET(user_data), gui.cut.tri);
	redraw(REDRAW_OVERLAY);
}

void on_spectro_x_axis_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *combo = GTK_DROP_DOWN(obj);
	(void)pspec;
	gui.cut.plot_as_wavenumber = gtk_drop_down_get_selected(combo);
}

void on_cut_cfa_toggled(GtkCheckButton *button, gpointer user_data) {
	gui.cut.cfa = siril_toggle_get_active(GTK_WIDGET(button));
	if (gui.cut.cfa) {
		gui.cut.tri = FALSE;
	}
	redraw(REDRAW_OVERLAY);
}

void on_cut_save_checkbutton_toggled(GtkCheckButton *button, gpointer user_data) {
	gui.cut.save_dat = siril_toggle_get_active(GTK_WIDGET(button));
}


void on_cut_spin_point1_value_changed(GtkSpinButton* button, gpointer user_data) {
	gboolean wl_changed = ((GtkWidget*) button == lookup_widget("cut_spin_wavelength1"));
	GtkSpinButton* wn1 = GTK_SPIN_BUTTON(lookup_widget("cut_spin_wavenumber1"));
	GtkSpinButton* wl1 = GTK_SPIN_BUTTON(lookup_widget("cut_spin_wavelength1"));
	double val = 10000000. / gtk_spin_button_get_value(button);
	if (wl_changed)
		gtk_spin_button_set_value(wn1, val);
	else
		gtk_spin_button_set_value(wl1, val);
}

void on_cut_spin_point2_value_changed(GtkSpinButton* button, gpointer user_data) {
	gboolean wl_changed = ((GtkWidget*) button == lookup_widget("cut_spin_wavelength2"));
	GtkSpinButton* wn2 = GTK_SPIN_BUTTON(lookup_widget("cut_spin_wavenumber2"));
	GtkSpinButton* wl2 = GTK_SPIN_BUTTON(lookup_widget("cut_spin_wavelength2"));
	double val = 10000000. / gtk_spin_button_get_value(button);
	if (wl_changed)
		gtk_spin_button_set_value(wn2, val);
	else
		gtk_spin_button_set_value(wl2, val);
}

void on_cut_dist_pref_as_group_changed(GtkCheckButton *button, gpointer user_data) {
	GtkToggleButton *as = (GtkToggleButton*) lookup_widget("cut_dist_pref_as");
	gui.cut.pref_as = siril_toggle_get_active(GTK_WIDGET(as));
}

void on_cut_spectro_polyorder_changed(GtkSpinButton* button, gpointer user_data) {
	gui.cut.bg_poly_order = gtk_spin_button_get_value(button);
}

void on_cut_spectro_plot_bg_toggled(GtkCheckButton *button, gpointer user_data) {
	gui.cut.plot_spectro_bg = siril_toggle_get_active(GTK_WIDGET(button));
}

//// Sequence Processing ////

