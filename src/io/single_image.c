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

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "algos/statistics.h"
#include "io/annotation_catalogues.h"
#include "algos/ccd-inspector.h"
#include "algos/background_extraction.h"
#include "algos/astrometry_solver.h"
#include "algos/demosaicing.h"
#include "core/gui_iface.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "core/undo.h"
#include "core/processing.h"

/* Closes and frees resources attached to the single image opened in gfit->
 * If a sequence is loaded and one of its images is displayed, nothing is done.
 */
void close_single_image() {
	if (sequence_is_loaded() && com.seq.current >= 0)
		return;
	/* We need to clear display and soft proofing transforms and a few other
	 * color management data */

	siril_log_debug("MODE: closing single image\n");
	undo_flush();
	/* we need to close all dialogs in order to avoid bugs
	 * with previews
	 */
	gui_iface.clear_roi();
	free_image_data();
}

/* free_image_data_gui moved to gui/gui_iface_impl.c */
/* frees resources when changing sequence or closing a single image */
void free_image_data() {
	siril_log_debug("free_image_data() called, clearing loaded image\n");
	/* WARNING: single_image.fit references the actual fits image,
	 * shouldn't it be used here instead of gfit? */
	cmsCloseProfile(gfit->icc_profile);
	gfit->icc_profile = NULL;
	reset_icc_transforms();
	if (!single_image_is_loaded() && sequence_is_loaded())
		save_stats_from_fit(gfit, &com.seq, com.seq.current);

	gui_iface.invalidate_histogram();

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

	if (!com.headless)
		gui_iface.on_image_closed();
	clearfits(gfit);
	siril_log_debug("free_image_data() complete\n");
}

static gboolean end_read_single_image(gpointer p) {
	gui_iface.set_GUI_CAMERA();
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
		siril_log_error(_("Error opening image %s: file not found or not supported.\n"), filename);
		free(realname);
		return 1;
	}
	if (imagetype == TYPESER || imagetype == TYPEAVI ||
				(imagetype == TYPEFITS && fitseq_is_fitseq(realname, NULL))) {
		if (allow_sequences) {
			retval = read_single_sequence(realname, imagetype);
			single_sequence = TRUE;
		} else {
			siril_log_error(_("Cannot open a sequence from here\n"));
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
		siril_log_error(_("Opening %s failed.\n"), realname);
	if (realname_out)
		*realname_out = realname;
	else
		free(realname);
	gui_iface.set_last_opened_filetype((int)imagetype);
	update_gain_from_gfit();
	siril_add_idle(end_read_single_image, NULL);
	return retval;
}

gboolean end_open_single_image(gpointer arg) {
	com.icc.srgb_hint = FALSE;
	if (!com.headless) {
		g_rw_lock_reader_lock(&gfit->rwlock);
		gui_iface.on_image_loaded();
		g_rw_lock_reader_unlock(&gfit->rwlock);
	}
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
	com.uniq->nb_layers = gfit->naxes[2];
	com.uniq->fit = gfit;
	return 0;
}

/* This function is used to load a single image, meaning outside a sequence,
 * whether a sequence is loaded or not, whether an image was already loaded or
 * not. The opened file is available in the usual global variable for current
 * image, gfit->
 */
int open_single_image(const char* filename) {
	int retval = 0;
	char *realname = NULL;
	gboolean is_single_sequence;

	/* Check we aren't running a processing thread otherwise it will clobber gfit
	 * when it finishes and cause a segfault.
	 */
	if ((retval = processing_is_job_active())) {
		siril_log_error(_("Cannot open another file while the processing thread is still operating on the current one!\n"));
	}

	/* first, close everything */
	if (!retval) {
		close_sequence(FALSE);	// closing a sequence if loaded
		close_single_image();	// close the previous image and free resources

		/* open the new file */
		retval = read_single_image(filename, gfit, &realname, TRUE, &is_single_sequence, TRUE, FALSE);
	}
	if (retval) {
		gui_iface.message_dialog(SIRIL_MSG_ERROR, _("Error opening file"),
				_("There was an error when opening this image. "
						"See the log for more information."));
		free(realname);
		return 1;
	}

	if (!is_single_sequence) {
		siril_log_debug("Loading image OK, now displaying\n");

		/* Now initializing com struct */
		com.seq.current = UNRELATED_IMAGE;
		create_uniq_from_gfit(realname, get_type_from_filename(realname) == TYPEFITS);
		if (!com.headless)
			gui_iface.execute_idle_sync(end_open_single_image, NULL);
	} else {
		free(realname);
	}
	gui_iface.reset_cut_gui_filedependent(NULL);
	gui_iface.check_icc_identical_to_monitor();
	return retval;
}

/* updates the GUI to reflect the opening of a single image, found in gfit and com.uniq */
gboolean open_single_image_from_gfit(gpointer user_data) {
	siril_log_debug("gui_function(open_single_image_from_gfit, NULL)\n");
	/* now initializing everything
	 * code based on seq_load_image or set_seq (sequence.c) */

	gui_iface.initialize_display_mode();

	gui_iface.update_zoom_label();

	init_layers_hi_and_lo_values(MIPSLOHI); // If MIPS-LO/HI exist we load these values. If not it is min/max

	gui_iface.sliders_mode_set_state(gui_iface.get_sliders_mode());
	gui_iface.set_cutoff_sliders_max_values();
	gui_iface.set_cutoff_sliders_values();

	gui_iface.update_display_mode_state();
	gui_iface.update_prepro_interface(TRUE);
	gui_iface.adjust_sellabel();

	gui_iface.display_filename();	// display filename in gray window
	gui_iface.set_precision_switch(); // set precision on screen

	/* update menus */
	gui_iface.update_menu_item();

	gui_iface.close_tab();
	gui_iface.init_right_tab();

	gui_iface.remap_all_vports();
	gui_iface.update_histogram();
	gui_iface.redraw_image(REMAP_ALL);
	return FALSE;
}

gboolean update_single_image_from_gfit(gpointer user_data) {
	/* a variation on open_single_image_from_gfit that only
	 does the things necessary when key aspects may have
	 changed (eg changed number of channels, bitpix etc.)*/

	g_rw_lock_reader_lock(&gfit->rwlock);
	init_layers_hi_and_lo_values(MIPSLOHI); // If MIPS-LO/HI exist we load these values. If not it is min/max

	gui_iface.sliders_mode_set_state(gui_iface.get_sliders_mode());
	gui_iface.set_cutoff_sliders_max_values();
	gui_iface.set_cutoff_sliders_values();

	gui_iface.set_precision_switch(); // set precision on screen

	gui_iface.close_tab();
	gui_iface.init_right_tab();

	gui_iface.remap_all_vports();
	gui_iface.update_histogram();
	g_rw_lock_reader_unlock(&gfit->rwlock);
	gui_iface.redraw_image(REMAP_ALL);
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

static void fit_lohi_to_layers(fits *fit, double lo_in, double hi_in, WORD *lo_out, WORD *hi_out) {
	if (fit->type == DATA_USHORT) {
		*lo_out = (WORD)lo_in;
		*hi_out = (WORD)hi_in;
	} else if (fit->type == DATA_FLOAT) {
		*lo_out = float_to_ushort_range((float)lo_in);
		*hi_out = float_to_ushort_range((float)hi_in);
	}
}

/* gfit has been loaded, now we copy the hi and lo values into the com.uniq or com.seq layers.
 * gfit->hi and gfit->lo may only be available in some FITS files; if they are not available, the
 * min and max value for the layer is used.
 * If gfit changed, its hi and lo values need to be updated, and they are taken from min and
 * max.
 */
void init_layers_hi_and_lo_values(sliders_mode force_minmax) {
	if (force_minmax == USER) return;
	WORD lo = 0, hi = 0xFFFF;
	sliders_mode sliders;
	if (gfit->keywords.hi == 0 || force_minmax == MINMAX) {
		sliders = MINMAX;
		image_find_minmax(gfit);
		fit_lohi_to_layers(gfit, gfit->mini, gfit->maxi, &lo, &hi);
	} else {
		sliders = MIPSLOHI;
		hi = gfit->keywords.hi;
		lo = gfit->keywords.lo;
	}
	gui_iface.update_display_range_after_load((int)sliders, (int)lo, (int)hi);
}

int single_image_is_loaded() {
	return (com.uniq != NULL && com.uniq->nb_layers > 0);
}

/**************** updating the single image *******************/

/* generic idle function for end of operation on gfit */
gboolean end_gfit_operation(gpointer data G_GNUC_UNUSED) {
	// this function should not contain anything required by the execution
	// of the operation because it won't be run in headless

	siril_log_debug("end of gfit operation - idle function\n");
	stop_processing_thread();

	// Check the mask tab visibility is correct
	gui_iface.show_or_hide_mask_tab_async();

	gui_iface.refresh_histogram_if_visible(); // histogram data already computed in notify_gfit_data_modified()

	/* update bit depth selector */
	gui_iface.on_precision_changed();

	/* update display of gfit name (useful if it changes) */
	gui_iface.adjust_sellabel();
	gui_iface.display_filename();

	// compute new min and max if needed for display and update sliders
	gui_iface.set_cutoff_sliders_values();

	/* re-enable the display-mode menu disabled at the start of single-image ops */
	gui_iface.enable_display_mode_menu();

	if (com.python_command) // must be synchronous to prevent a crash where this is still running while the next command runs
		gui_iface.redraw_image(REMAP_ALL);
	else
		gui_iface.redraw_image_async(REMAP_ALL);	// queues a redraw if !com.script

	gui_iface.redraw_previews();

	gui_iface.set_busy(FALSE); // called from current thread if !com.script, idle else
	return FALSE;
}

/* to be called after each operation that modifies the content of gfit, at the
 * end of a processing operation, not for previews */
void gfit_modified_update_gui() {
	gui_iface.execute_idle_sync(end_gfit_operation, NULL);
}

/* Must be called on the data-processing thread (worker or script thread) after
 * gfit data is modified, before the thread ends.  Handles all aspects relating
 * to gfit itself: invalidates cached statistics and histogram, computes fresh
 * histogram data, remaps the Cairo display buffers, and recomputes the display
 * range.  The pixel work is thread-safe; any GTK widget-state updates it
 * triggers are dispatched as idle callbacks so they run on the main thread.
 *
 * In non-Python script mode the expensive work (histogram computation, Cairo
 * remap, display-range recalculation) is skipped: there is no point updating
 * the display buffers for every intermediate command in a script when the
 * result will never be shown until the script ends.  The cache-invalidation
 * calls are always made so that subsequent commands see stale stats/histograms
 * as dirty and recompute them on demand.  execute_script() calls this function
 * again after clearing com.script so the final result is displayed correctly. */
void notify_gfit_data_modified() {
	invalidate_stats_from_fit(gfit);
	// The following are only required in GUI mode
	if (!com.headless) {
		/* Hold histogram_mutex across the invalidate+recompute pair so the
		 * main thread never sees a partially-nullified layers_hist[] array.
		 * update_histo_mtf() on the main thread acquires the same mutex. */
		g_mutex_lock(&com.histogram_mutex);
		gui_iface.invalidate_histogram();
		// Skip expensive pixel work mid-script; display is flushed at script end.
		if (com.script && !com.python_script) {
			g_mutex_unlock(&com.histogram_mutex);
			return;
		}
		/* Do not remap if a job on the worker thread is currently writing gfit.
		 * The worker calls this function itself from within generic_image_worker,
		 * so we must allow that path through — hence the processing_in_worker_thread()
		 * exemption.  Any other caller that arrives while a job is active (e.g. a
		 * GTK slider callback during a Python script command) would race against the
		 * worker reading/writing gfit pixel data.  The job's own call, and the
		 * subsequent end_gfit_operation idle, will update the display correctly. */
		if (!processing_in_worker_thread() && processing_is_job_active()) {
			g_mutex_unlock(&com.histogram_mutex);
			return;
		}
		/* If a ROI is active and contains processed data, merge it back into
		 * gfit now — before computing the histogram and before remap_all()
		 * builds the Cairo display buffers — so that both operations see the
		 * fully-updated pixel data.  This is the correct point to do this:
		 * gui_iface.redraw_image() must remain a pure "repaint from Cairo buffers" function
		 * and must not write gfit. */
		fits *roi_fit = (fits*)gui_iface.get_roi_fit();
		if (gui_iface.roi_is_active() && gui_iface.roi_operation_supports() &&
				roi_fit &&
				((gfit->type == DATA_FLOAT && roi_fit->fdata) ||
				 (gfit->type == DATA_USHORT && roi_fit->data)))
			gui_iface.copy_roi_into_gfit();
		gui_iface.compute_histo_for_fit(gfit); // reads gfit pixel data; GTK toggle update deferred to idle
		g_mutex_unlock(&com.histogram_mutex);
		gui_iface.remap_all_vports(); // Updates the Cairo image buffers based on applying the remap LUT to gfit
		/* gui.hi / gui.lo are read on the GTK main thread (set_cutoff_sliders_values);
		 * protect the write with com.mutex to prevent a data race. */
		g_mutex_lock(&com.mutex);
		init_layers_hi_and_lo_values((sliders_mode)gui_iface.get_sliders_mode());
		g_mutex_unlock(&com.mutex);
	}
}

gboolean enforce_area_in_fits(fits *fit, rectangle *area) {
        gboolean has_crossed = FALSE;
        int rx = fit->rx, ry = fit->ry;
        if (area->w > rx) {
                area->w = rx;
                has_crossed = TRUE;
        }
        if (area->h > ry) {
                area->h = ry;
                has_crossed = TRUE;
        }
        if (area->x < 0) {
                area->x = 0;
                has_crossed = TRUE;
        }
        if (area->y < 0) {
                area->y = 0;
                has_crossed = TRUE;
        }
        if (area->x + area->w > rx) {
                area->x = rx - area->w;
                has_crossed = TRUE;
        }
        if (area->y + area->h > ry) {
                area->y = ry - area->h;
                has_crossed = TRUE;
        }
        return has_crossed;
}

