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

#include "core/siril.h"
#include "core/proto.h"
#include "core/command.h"
#include "core/undo.h"
#include "core/icc_profile.h"
#include "core/command.h"
#include "core/siril_update.h"
#include "core/siril_cmd_help.h"
#include "core/siril_log.h"
#include "core/initfile.h"
#include "compositing/align_rgb.h"
#include "io/annotation_catalogues.h"
#include "algos/astrometry_solver.h"
#include "algos/noise.h"
#include "algos/geometry.h"
#include "algos/photometry.h"
#include "algos/siril_wcs.h"
#include "algos/ccd-inspector.h"
#include "compositing/compositing.h"
#include "filters/synthstar.h"
#include "gui/about_dialog.h"
#include "gui/utils.h"
#include "gui/colors.h"
#include "gui/callbacks.h"
#include "gui/documentation.h"
#include "gui/histogram.h"
#include "gui/icc_profile.h"
#include "gui/open_dialog.h"
#include "gui/message_dialog.h"
#include "gui/PSF_list.h"
#include "gui/save_dialog.h"
#include "gui/sequence_list.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/script_menu.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/photometric_cc.h"
#include "gui/menu_gray_geometry.h"
#include "gui/registration_preview.h"
#include "gui/remixer.h"
#include "livestacking/livestacking.h"
#include "registration/registration.h"
#include "io/siril_catalogues.h"

#include "siril_actions.h"

void open_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	header_open_button_clicked();
	launch_clipboard_survey();
}

void cwd_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	cwd_btton_clicked();
}

void livestacking_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWidget *w = lookup_widget("livestacking_player");

	gtk_widget_show(w);
	gtk_window_set_keep_above(GTK_WINDOW(w), TRUE);
}

void save_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_save_button_clicked();
}

void save_as_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_save_as_button_clicked();
}

void snapshot_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_snapshot_button_clicked(FALSE);
}

void clipboard_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_snapshot_button_clicked(TRUE);
}

void undo_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	set_cursor_waiting(TRUE);
	undo_display_data(UNDO);
	set_cursor_waiting(FALSE);
}

void redo_action_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	set_cursor_waiting(TRUE);
	undo_display_data(REDO);
	set_cursor_waiting(FALSE);
}

void quit_action_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	siril_quit();
}

void about_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_show_about_dialog();
}

void preferences_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("settings_window");
}

void close_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	process_close(0);
}

void scripts_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("settings_window");
	gtk_stack_set_visible_child((GtkStack*) lookup_widget("stack_pref"), lookup_widget("scripts_page"));
//	siril_get_on_script_pages();
}

void updates_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
#if defined(HAVE_LIBCURL)
	siril_check_updates(TRUE);
#else
	siril_log_message(_("Cannot check for updates with this version, missing dependency\n"));
#endif
}

void doc_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_documentation();
}

static gboolean is_extended = FALSE;

void full_screen_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWindow *window;
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	GtkWidget *control_center_box = lookup_widget("control_center_box");
	GtkButton *button = GTK_BUTTON(lookup_widget("button_paned"));
	gboolean is_control_box_visible;
	gboolean is_fullscreen;

	window = GTK_WINDOW(GTK_APPLICATION_WINDOW(user_data));

	GdkWindow *gdk_window = gtk_widget_get_window(GTK_WIDGET(window));
	is_fullscreen = gdk_window_get_state(gdk_window) & GDK_WINDOW_STATE_FULLSCREEN;
	is_control_box_visible = gtk_widget_get_visible(control_center_box);

	if (is_fullscreen) {
		gtk_window_unfullscreen(window);
		if (is_extended)
			gtk_button_clicked(button);
	} else {
		gtk_window_fullscreen(window);
		if (is_control_box_visible) {
			gtk_button_clicked(button);
		}
		is_extended = is_control_box_visible;
	}
	gtk_widget_set_visible(toolbarbox, is_fullscreen);
}

void panel_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
	GtkImage *image = GTK_IMAGE(gtk_bin_get_child(GTK_BIN(GTK_BUTTON(lookup_widget("button_paned")))));
	GtkWidget *widget = gtk_paned_get_child2(paned);

	gboolean is_visible = gtk_widget_is_visible(widget);

	gtk_widget_set_visible(widget, !is_visible);

	if (!is_visible) {
		gtk_image_set_from_icon_name(image, "pan-end-symbolic", GTK_ICON_SIZE_BUTTON);
		if (gui.icc.iso12646)
			disable_iso12646_conditions(TRUE, FALSE, TRUE);
	} else {
		gtk_image_set_from_icon_name(image, "pan-start-symbolic", GTK_ICON_SIZE_BUTTON);
	}
	if (com.pref.gui.remember_windows) {
		com.pref.gui.is_extended = !is_visible;
		writeinitfile();
	}
}

void keyboard_shortcuts_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWindow *window;

	window = GTK_WINDOW(GTK_APPLICATION_WINDOW(user_data));

	siril_cmd_help_keyboard_shortcuts(window);
}

void tab_conversion_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(FILE_CONVERSION);
}

void tab_sequence_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(IMAGE_SEQ);
}

void tab_prepro_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PRE_PROC);
}

void tab_registration_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(REGISTRATION);
}

void tab_plot_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PLOT);
}

void tab_stacking_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

void tab_logs_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(OUTPUT_LOGS);
}

void toolbar_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWidget *w = lookup_widget("toolbarbox");
	gtk_widget_set_visible(w, !gtk_widget_get_visible(w));
}

void change_zoom_fit_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	if (g_variant_get_boolean(state)) {
		gui.zoom_value = ZOOM_FIT;
		if (gui.icc.iso12646)
			disable_iso12646_conditions(FALSE, TRUE, TRUE);
		reset_display_offset();
		redraw(REDRAW_IMAGE);
	} else {
		gui.zoom_value = get_zoom_val();
	}
	update_zoom_label();
	g_simple_action_set_state(action, state);
}

void zoom_fit_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void zoom_in_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	point center = get_center_of_vport();
	update_zoom(center.x, center.y, ZOOM_IN);
	if (gui.icc.iso12646)
		disable_iso12646_conditions(FALSE, TRUE, TRUE);
}

void zoom_out_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	point center = get_center_of_vport();
	update_zoom(center.x, center.y, ZOOM_OUT);
	if (gui.icc.iso12646)
		disable_iso12646_conditions(FALSE, TRUE, TRUE);
}

void zoom_one_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	update_zoom_fit_button();
	gui.zoom_value = ZOOM_NONE;
	reset_display_offset();
	update_zoom_label();
	redraw(REDRAW_IMAGE);
}

void negative_view_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	g_simple_action_set_state(action, state);
	set_cursor_waiting(TRUE);
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void negative_view_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void photometry_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	mouse_status = g_variant_get_boolean(state) ? MOUSE_ACTION_PHOTOMETRY : MOUSE_ACTION_SELECT_REG_AREA;
	g_simple_action_set_state(action, state);
	free(gui.qphot);
	gui.qphot = NULL;
	redraw(REDRAW_OVERLAY);
}

void photometry_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void color_map_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	g_simple_action_set_state(action, state);
	set_cursor_waiting(TRUE);
	redraw(REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void color_map_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void chain_channels_state_change(GSimpleAction *action, GVariant *state, gpointer user_data) {
	g_simple_action_set_state(action, state);
	set_unlink_channels(!g_variant_get_boolean(state));
}

void chain_channels_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void astrometry_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	open_astrometry_dialog();
}

void dyn_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	int confirm = TRUE;
	if (gfit.naxes[2] == 1 && gfit.keywords.bayer_pattern[0] != '\0') {
		confirm = siril_confirm_dialog(_("Undebayered CFA image loaded"),
				_("The star detection functions in the Dynamic PSF dialog may produce results for this CFA image but will not perform optimally and star parameters may be inaccurate. Are you sure you wish to proceed?"), _("Proceed"));
	}
	if (confirm)
		siril_open_dialog("stars_list_window");
}

void pick_star_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	pick_a_star();
}

void psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	psf_star *result = NULL;
	int layer = select_vport(gui.cvport);

	if (layer == -1)
		return;
	if (!(com.selection.h && com.selection.w))
		return;
	if (com.selection.w > 300 || com.selection.h > 300) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Current selection is too large"),
				_("To determine the PSF, please make a selection around a star."));
		return;
	}
	struct phot_config *ps = phot_set_adjusted_for_image(&gfit);
	result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE, ps, TRUE, com.pref.starfinder_conf.profile, NULL);
	free(ps);
	if (!result)
		return;

	popup_psf_result(result, &com.selection, &gfit);
	free_psf(result);
}

void seq_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	process_seq_psf(0);
}

void crop_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_crop();
}

void seq_crop_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("crop_dialog");
}

void search_object_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (has_wcs(&gfit))
		siril_open_dialog("search_objects");
}

void search_object_solar_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (has_wcs(&gfit)) {
		siril_catalogue *siril_cat = siril_catalog_fill_from_fit(&gfit, CAT_IMCCE, -1.f);
		if (com.pref.astrometry.default_obscode != NULL) {
			if (strlen(com.pref.astrometry.default_obscode) == 3) {
				siril_cat->IAUcode = g_strdup(com.pref.astrometry.default_obscode);
				siril_log_message(_("Using preferred observatory code %s\n"), siril_cat->IAUcode);
			} else {
				g_free(com.pref.astrometry.default_obscode);
				com.pref.astrometry.default_obscode = NULL;
				siril_cat->IAUcode = g_strdup("500");
				siril_log_color_message(_("Invalid observatory code found in preferences. Code must be 3 characters. Resetting preference and using default geocentric observer code.\n"), "red");
			}
		} else {
			siril_cat->IAUcode = g_strdup("500");
			siril_log_message(_("Using default geocentric observer code\n"));
		}
		conesearch_args *args = init_conesearch();
		args->fit = &gfit;
		args->has_GUI = TRUE;
		args->siril_cat = siril_cat;
		args->display_log = TRUE;
		args->display_tag = TRUE;
		start_in_new_thread(conesearch_worker, args);
	}
}

void annotate_object_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	if (g_variant_get_boolean(state)) {
		if (has_wcs(&gfit)) {
			com.found_object = find_objects_in_field(&gfit);
		}
	} else {
		g_slist_free(com.found_object);
		com.found_object = NULL;
		purge_user_catalogue(CAT_AN_USER_TEMP);
		refresh_annotation_visibility();
	}
	g_simple_action_set_state(action, state);
	redraw(REDRAW_OVERLAY);
}

void wcs_grid_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	gui.show_wcs_grid = g_variant_get_boolean(state);
	g_simple_action_set_state(action, state);
	redraw(REDRAW_OVERLAY);
}

void annotate_object_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void wcs_grid_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void regframe_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	GtkToggleButton *drawframe;
	drawframe = GTK_TOGGLE_BUTTON(lookup_widget("drawframe_check"));
	gtk_toggle_button_set_active(drawframe, g_variant_get_boolean(state));
	g_simple_action_set_state(action, state);
	redraw(REDRAW_OVERLAY);
}

void regframe_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;
	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void seq_list_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_widget_get_visible(lookup_widget("seqlist_dialog"))) {
		siril_close_dialog("seqlist_dialog");
	} else {
		gboolean confirm = TRUE;
		if (com.seq.current == RESULT_IMAGE) {
			confirm = siril_confirm_dialog(_("Save your changes before loading a frame of the sequence."),
					_("The image currently displayed is the result of the previous stack. "
							"If you load an image from the sequence, you might lose the entire process you performed on the image, "
							"but not the image itself. You need to save your data before doing this."),
					_("Load another image"));
		}
		if (confirm) {
			int layer = get_registration_layer(&com.seq);
			update_seqlist(layer);
			siril_open_dialog("seqlist_dialog");
		}
	}
}

void statistics_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	set_cursor_waiting(TRUE);
	computeStat();
	siril_open_dialog("StatWindow");
	set_cursor_waiting(FALSE);
}

void noise_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	evaluate_noise_in_image();
}

void ccd_inspector_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	compute_aberration_inspector();
}

void show_tilt_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void show_tilt_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	set_cursor_waiting(TRUE);
	if (g_variant_get_boolean(state)) {
		draw_sensor_tilt(&gfit);

	} else {
		clear_sensor_tilt();
		redraw(REDRAW_OVERLAY);
	}
	g_simple_action_set_state(action, state);
	set_cursor_waiting(FALSE);
}

void image_information_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("file_information");
}

void image_fits_header_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("keywords_dialog");
}

/******* processing menu **************/

void remove_green_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("SCNR_dialog");
}

void saturation_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("satu_dialog");
}

void color_calib_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	initialize_calibration_interface();
	siril_open_dialog("color_calibration");
}

void pcc_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	initialize_photometric_cc_dialog();
	siril_open_dialog("s_pcc_dialog");
}

void spcc_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	initialize_spectrophotometric_cc_dialog();
	siril_open_dialog("s_pcc_dialog");
}

void split_channel_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	siril_open_dialog("extract_channel_dialog");
}

void negative_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	negative_processing();
}

void histo_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	toggle_histogram_window_visibility(1);
}

void fix_banding_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("canon_fixbanding_dialog");
}

void cosmetic_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("cosmetic_dialog");
}

void background_extr_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("background_extraction_dialog");
}

void asinh_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("asinh_dialog");
}

void epf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("epf_dialog");
}
void starnet_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("starnet_dialog");
}
void deconvolution_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("bdeconv_dialog");
}

void payne_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	toggle_histogram_window_visibility(2);
}

void binning_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("binxy_dialog");
}

void resample_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("resample_dialog");
}

void rotation_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (com.selection.w == 0 || com.selection.h == 0) {
		com.selection = (rectangle){ 0, 0, gfit.rx, gfit.ry };
	}
	siril_open_dialog("rotation_dialog");
	redraw(REDRAW_OVERLAY);
}

void rotation90_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_rotate90();
}

void rotation270_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_rotate270();
}

void mirrorx_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	mirrorx_gui(&gfit);
}

void mirrory_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	mirrory_gui(&gfit);
}

void wavelets_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	siril_open_dialog("wavelets_dialog");
}

void split_wavelets_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("extract_wavelets_layers_dialog");
}

void medianfilter_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("Median_dialog");
}

void rgradient_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("rgradient_dialog");
}

void clahe_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("CLAHE_dialog");
}

void linearmatch_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_set_file_filter("reference_filechooser_linearmatch", "filefilter_fits");
	siril_open_dialog("linearmatch_dialog");
}

void fft_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkFileChooserButton *magbutton, *phasebutton;

	magbutton = GTK_FILE_CHOOSER_BUTTON(lookup_widget("filechooser_mag"));
	phasebutton = GTK_FILE_CHOOSER_BUTTON(lookup_widget("filechooser_phase"));
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(magbutton), com.wd);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(phasebutton), com.wd);
	siril_set_file_filter("filechooser_mag", "filefilter_fits");
	siril_set_file_filter("filechooser_phase", "filefilter_fits");
	siril_open_dialog("dialog_FFT");
}

void rgb_compositing_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	open_compositing_window();
}

void star_remix_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	toggle_remixer_window_visibility(CALL_FROM_MENU, NULL, NULL);
}

void pixel_math_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("pixel_math_dialog");
}

void split_cfa_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("split_cfa_dialog");
}

void nina_lc_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("nina_light_curve");
}

void compstars_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("compstars");
}

void denoise_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("denoise_dialog");
}

void merge_cfa_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("merge_cfa_dialog");
}

void star_desaturate_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	undo_save_state(&gfit, "Synthetic stars: desaturate clipped stars");
	control_window_switch_to_tab(OUTPUT_LOGS);
	start_in_new_thread(fix_saturated_stars, NULL);
}

void star_synthetic_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (!check_ok_if_cfa())
		return;
	undo_save_state(&gfit, "Synthetic stars: full replacement");
	control_window_switch_to_tab(OUTPUT_LOGS);
	start_in_new_thread(do_synthstar, NULL);
}

void align_dft_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	undo_save_state(&gfit, _("RGB alignment (DFT)"));
	rgb_align(1);
}

void align_global_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	undo_save_state(&gfit, _("RGB alignment (Global stars)"));
	rgb_align(2);
}

void align_kombat_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	undo_save_state(&gfit, _("RGB alignment (KOMBAT)"));
	rgb_align(3);
}

void align_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	undo_save_state(&gfit, _("RGB alignment (PSF)"));
	if (com.selection.w > 300 || com.selection.h > 300) {
		siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
		return;
	}
	rgb_align(0);
}

void icc_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("icc_dialog");
}

void cut_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkToggleToolButton *button = (GtkToggleToolButton*) lookup_widget("cut_button");
	if (gtk_toggle_tool_button_get_active(button)) {
		mouse_status = MOUSE_ACTION_CUT_SELECT;
		siril_open_dialog("cut_dialog");
	} else {
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		siril_close_dialog("cut_coords_dialog");
		siril_close_dialog("cut_spectroscopy_dialog");
		siril_close_dialog("cut_dialog");
	}
}

void clear_roi(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_clear_roi();
}

void set_roi(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_set_roi();
}

void ccm_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("ccm_dialog");
}
