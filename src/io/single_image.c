/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#include <string.h>

#include "core/siril.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "algos/statistics.h"
#include "io/annotation_catalogues.h"
#include "algos/ccd-inspector.h"
#include "algos/background_extraction.h"
#include "algos/astrometry_solver.h"
#include "algos/demosaicing.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/utils.h"
#include "gui/cut.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/icc_profile.h"
#include "gui/message_dialog.h"
#include "gui/plot.h"
#include "gui/registration_preview.h"
#include "gui/registration.h"
#include "gui/user_polygons.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "core/undo.h"
#include "core/processing.h"

/* Closes and frees resources attached to the single image opened in gfit.
 * If a sequence is loaded and one of its images is displayed, nothing is done.
 */
void close_single_image() {
	if (sequence_is_loaded() && com.seq.current >= 0)
		return;
	/* We need to clear display and soft proofing transforms and a few other
	 * color management data */

	siril_debug_print("MODE: closing single image\n");
	undo_flush();
	/* we need to close all dialogs in order to avoid bugs
	 * with previews
	 */
	on_clear_roi();
	free_image_data();
}

static gboolean free_image_data_gui(gpointer p) {
	disable_iso12646_conditions(TRUE, FALSE, FALSE);
	//reset_compositing_module();
	clear_user_polygons(); // clear list of user polygons
	delete_selected_area(); // this triggers a redraw
	reset_plot(); // clear existing plot if any
	siril_close_preview_dialogs();
	/* It is better to close all other dialog. Indeed, some dialog are not compatible with all images */
	display_filename();
	update_zoom_label();
	update_display_fwhm();
	adjust_sellabel();
	gui_function(update_MenuItem, NULL);
	reset_3stars();
	gui_function(close_tab, NULL);	// close Green and Blue tabs
	free_cut_args(&gui.cut);
	initialize_cut_struct(&gui.cut);

	GtkComboBox *binning = GTK_COMBO_BOX(gtk_builder_get_object(gui.builder, "combobinning"));
	GtkEntry* focal_entry = GTK_ENTRY(lookup_widget("focal_entry"));
	GtkEntry* pitchX_entry = GTK_ENTRY(lookup_widget("pitchX_entry"));
	GtkEntry* pitchY_entry = GTK_ENTRY(lookup_widget("pitchY_entry"));
	// avoid redrawing plot while com.seq has not been updated
	g_signal_handlers_block_by_func(focal_entry, on_focal_entry_changed, NULL);
	g_signal_handlers_block_by_func(pitchX_entry, on_pitchX_entry_changed, NULL);
	g_signal_handlers_block_by_func(pitchY_entry, on_pitchY_entry_changed, NULL);
	g_signal_handlers_block_by_func(binning, on_combobinning_changed, NULL);
	clear_stars_list(TRUE);
	clear_sampling_setting_box();	// clear focal and pixel pitch info
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	cleanup_annotation_catalogues(TRUE);
	reset_display_offset();
	reset_menu_toggle_button();
	reset_zoom_default();
	free(gui.qphot);
	gui.qphot = NULL;
	gui.show_wcs_disto = FALSE;
	clear_sensor_tilt();
	g_signal_handlers_unblock_by_func(focal_entry, on_focal_entry_changed, NULL);
	g_signal_handlers_unblock_by_func(pitchX_entry, on_pitchX_entry_changed, NULL);
	g_signal_handlers_unblock_by_func(pitchY_entry, on_pitchY_entry_changed, NULL);
	g_signal_handlers_unblock_by_func(binning, on_combobinning_changed, NULL);
	siril_debug_print("free_image_data_idle() complete\n");

	/* free display image data */
	for (int vport = 0; vport < MAXVPORT; vport++) {
		struct image_view *view = &gui.view[vport];
		if (view->buf) {
			free(view->buf);
			view->buf = NULL;
		}
		if (view->full_surface) {
			cairo_surface_destroy(view->full_surface);
			view->full_surface = NULL;
		}
		view->full_surface_stride = 0;
		view->full_surface_height = 0;

		if (view->disp_surface) {
			cairo_surface_destroy(view->disp_surface);
			view->disp_surface = NULL;
		}
		view->view_width = -1;
		view->view_height= -1;
	}
	clear_previews();
	free_reference_image();
	siril_debug_print("free_image_data_gui() complete\n");
	return FALSE;
}

/* frees resources when changing sequence or closing a single image */
void free_image_data() {
	siril_debug_print("free_image_data() called, clearing loaded image\n");
	/* WARNING: single_image.fit references the actual fits image,
	 * shouldn't it be used here instead of gfit? */
	cmsCloseProfile(gfit.icc_profile);
	gfit.icc_profile = NULL;
	reset_icc_transforms();
	if (!single_image_is_loaded() && sequence_is_loaded())
		save_stats_from_fit(&gfit, &com.seq, com.seq.current);

	invalidate_gfit_histogram();

	if (com.uniq) {
		free(com.uniq->filename);
		com.uniq->fileexist = FALSE;
		free(com.uniq->comment);
		free(com.uniq);
		com.uniq = NULL;
	}
	/* this function frees resources used in the GUI, some of these resources
	 * need to be handled in the GTK+ main thread, so we use an idle function
	 * to deal with them */

	if (!com.headless) {
		if (com.script || com.python_command) {
			execute_idle_and_wait_for_it(free_image_data_gui, NULL);
		} else if (!g_main_context_is_owner(g_main_context_default())) {
			siril_add_idle(free_image_data_gui, NULL);
		} else {
			free_image_data_gui(NULL);
		}
	}
	clearfits(&gfit);
	siril_debug_print("free_image_data() complete\n");
}

static gboolean end_read_single_image(gpointer p) {
	set_GUI_CAMERA();
	return FALSE;
}

/**
 * Reads an image from disk and stores it in the user allocated destination
 * fits.
 * @param filename
 * @param dest
 * @param realname_out argument can be NULL, and if not, it is set to the
 * real file name of the loaded file, since the given filename can be without
 * extension.
 * @param is_sequence is set to TRUE if the loaded image is in fact a SER or AVI sequence. Can be NULL
 * @return 0 on success
 */
int read_single_image(const char *filename, fits *dest, char **realname_out,
		gboolean allow_sequences, gboolean *is_sequence, gboolean allow_dialogs,
		gboolean force_float) {
	int retval;
	image_type imagetype;
	char *realname = NULL;
	gboolean single_sequence = FALSE;

	retval = stat_file(filename, &imagetype, &realname);
	if (retval) {
		siril_log_color_message(_("Error opening image %s: file not found or not supported.\n"), "red", filename);
		free(realname);
		return 1;
	}
	if (imagetype == TYPESER || imagetype == TYPEAVI ||
				(imagetype == TYPEFITS && fitseq_is_fitseq(realname, NULL))) {
		if (allow_sequences) {
			retval = read_single_sequence(realname, imagetype);
			single_sequence = TRUE;
		} else {
			siril_log_color_message(_("Cannot open a sequence from here\n"), "red");
			free(realname);
			return 1;
		}
	} else {
		retval = any_to_fits(imagetype, realname, dest, allow_dialogs, force_float);
		if (!retval)
			debayer_if_needed(imagetype, dest, FALSE);
		if (com.pref.debayer.open_debayer || imagetype != TYPEFITS)
			update_fits_header(dest);
	}
	if (is_sequence) {
		*is_sequence = single_sequence;
	}
	if (retval && retval != OPEN_IMAGE_CANCEL)
		siril_log_color_message(_("Opening %s failed.\n"), "red", realname);
	if (realname_out)
		*realname_out = realname;
	else
		free(realname);
	gui.file_ext_filter = (int) imagetype;
	update_gain_from_gfit();
	siril_add_idle(end_read_single_image, NULL);
	return retval;
}

gboolean end_open_single_image(gpointer arg) {
	com.icc.srgb_hint = FALSE;
	if (!com.headless)
		open_single_image_from_gfit(NULL);
	return FALSE;
}

/* filename will be freed when the unique file is closed */
int create_uniq_from_gfit(char *filename, gboolean exists) {
	com.uniq = calloc(1, sizeof(single));
	if (!com.uniq) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	com.uniq->filename = filename;
	com.uniq->fileexist = exists;
	com.uniq->nb_layers = gfit.naxes[2];
	com.uniq->fit = &gfit;
	return 0;
}

/* This function is used to load a single image, meaning outside a sequence,
 * whether a sequence is loaded or not, whether an image was already loaded or
 * not. The opened file is available in the usual global variable for current
 * image, gfit.
 */
int open_single_image(const char* filename) {
	int retval = 0;
	char *realname = NULL;
	gboolean is_single_sequence;

	/* Check we aren't running a processing thread otherwise it will clobber gfit
	 * when it finishes and cause a segfault.
	 */
	if ((retval = get_thread_run())) {
		siril_log_message(_("Cannot open another file while the processing thread is still operating on the current one!\n"));
	}

	/* first, close everything */
	if (!retval) {
		close_sequence(FALSE);	// closing a sequence if loaded
		close_single_image();	// close the previous image and free resources

		/* open the new file */
		retval = read_single_image(filename, &gfit, &realname, TRUE, &is_single_sequence, TRUE, FALSE);
	}
	if (retval) {
		queue_message_dialog(GTK_MESSAGE_ERROR, _("Error opening file"),
				_("There was an error when opening this image. "
						"See the log for more information."));
		free(realname);
		return 1;
	}

	if (!is_single_sequence) {
		siril_debug_print("Loading image OK, now displaying\n");

		/* Now initializing com struct */
		com.seq.current = UNRELATED_IMAGE;
		create_uniq_from_gfit(realname, get_type_from_filename(realname) == TYPEFITS);
		if (!com.headless)
			execute_idle_and_wait_for_it(end_open_single_image, NULL);
	} else {
		free(realname);
	}
	gui_function(reset_cut_gui_filedependent, NULL);
	check_gfit_profile_identical_to_monitor();
	return retval;
}

/* updates the GUI to reflect the opening of a single image, found in gfit and com.uniq */
gboolean open_single_image_from_gfit(gpointer user_data) {
	siril_debug_print("gui_function(open_single_image_from_gfit, NULL)\n");
	/* now initializing everything
	 * code based on seq_load_image or set_seq (sequence.c) */

	initialize_display_mode();

	update_zoom_label();

	init_layers_hi_and_lo_values(MIPSLOHI); // If MIPS-LO/HI exist we load these values. If not it is min/max

	sliders_mode_set_state(gui.sliders);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();

	set_display_mode();
	update_prepro_interface(TRUE);
	adjust_sellabel();

	display_filename();	// display filename in gray window
	set_precision_switch(NULL); // set precision on screen

	/* update menus */
	update_MenuItem(NULL);

	close_tab(NULL);
	init_right_tab(NULL);

	update_gfit_histogram_if_needed();
	redraw(REMAP_ALL);
	return FALSE;
}

gboolean update_single_image_from_gfit(gpointer user_data) {
	/* a variation on open_single_image_from_gfit that only
	 does the things necessary when key aspects may have
	 changed (eg changed number of channels, bitpix etc.)*/

	init_layers_hi_and_lo_values(MIPSLOHI); // If MIPS-LO/HI exist we load these values. If not it is min/max

	sliders_mode_set_state(gui.sliders);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();

	set_precision_switch(NULL); // set precision on screen

	close_tab(NULL);
	init_right_tab(NULL);

	update_gfit_histogram_if_needed();
	redraw(REMAP_ALL);
	return FALSE;
}

/* searches the image for minimum and maximum pixel value, on each layer
 * the values are stored in fit->min[layer] and fit->max[layer] */
int image_find_minmax(fits *fit) {
	int layer;
	if (fit->maxi > 0.0)
		return 0;
	fit->mini = DBL_MAX;
	fit->maxi = -DBL_MAX;
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		// calling statistics() saves stats in the fit already, we don't need
		// to use the returned handle
		free_stats(statistics(NULL, -1, fit, layer, NULL, STATS_MINMAX, SINGLE_THREADED));
		if (!fit->stats || !fit->stats[layer])
			return -1;
		fit->maxi = max(fit->maxi, fit->stats[layer]->max);
		fit->mini = min(fit->mini, fit->stats[layer]->min);
	}
	return 0;
}

static void fit_lohi_to_layers(fits *fit, double lo, double hi) {
	if (fit->type == DATA_USHORT) {
		gui.lo = (WORD)lo;
		gui.hi = (WORD)hi;
	}
	else if (fit->type == DATA_FLOAT) {
		gui.lo = float_to_ushort_range((float)lo);
		gui.hi = float_to_ushort_range((float)hi);
	}
}

/* gfit has been loaded, now we copy the hi and lo values into the com.uniq or com.seq layers.
 * gfit.hi and gfit.lo may only be available in some FITS files; if they are not available, the
 * min and max value for the layer is used.
 * If gfit changed, its hi and lo values need to be updated, and they are taken from min and
 * max.
 */
void init_layers_hi_and_lo_values(sliders_mode force_minmax) {
	if (force_minmax == USER) return;
	if (gfit.keywords.hi == 0 || force_minmax == MINMAX) {
		gui.sliders = MINMAX;
		image_find_minmax(&gfit);
		fit_lohi_to_layers(&gfit, gfit.mini, gfit.maxi);
	} else {
		gui.sliders = MIPSLOHI;
		gui.hi = gfit.keywords.hi;
		gui.lo = gfit.keywords.lo;
	}
}

int single_image_is_loaded() {
	return (com.uniq != NULL && com.uniq->nb_layers > 0);
}

/**************** updating the single image *******************/

/* generic idle function for end of operation on gfit */
gboolean end_gfit_operation(gpointer data) {
	// Remove unused argument warnings
	(void) data;
	// this function should not contain anything required by the execution
	// of the operation because it won't be run in headless

	siril_debug_print("end of gfit operation - idle function\n");
	stop_processing_thread();

	update_gfit_histogram_if_needed();

	/* update bit depth selector */
	gui_function(set_precision_switch, NULL);

	/* update display of gfit name (useful if it changes) */
	adjust_sellabel();
	display_filename();

	// compute new min and max if needed for display and update sliders
	init_layers_hi_and_lo_values(gui.sliders);
	set_cutoff_sliders_values();

	if (com.python_command) // must be synchronous to prevent a crash where this is still running while the next command runs
		redraw(REMAP_ALL);
	else
		queue_redraw(REMAP_ALL);	// queues a redraw if !com.script

	gui_function(redraw_previews, NULL);	// queues redraws of the registration previews if !com.script

	set_cursor_waiting(FALSE); // called from current thread if !com.script, idle else
	return FALSE;
}

/* to be called after each operation that modifies the content of gfit, at the
 * end of a processing operation, not for previews */
void notify_gfit_modified() {
	siril_debug_print("end of gfit operation\n");
	invalidate_stats_from_fit(&gfit);
	invalidate_gfit_histogram();

	gui_function(end_gfit_operation, NULL);
}
