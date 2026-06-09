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

#include "siril_preview.h"

#include "core/siril.h"
#include "gui-gtk4/gui_state.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/siril_log.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/callbacks.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "algos/statistics.h"


static gboolean notify_is_blocked;
static gboolean preview_is_active;

/* ── Preview scheduling state ───────────────────────────────────────────────
 *
 * All fields are accessed exclusively from the GTK main thread.
 *
 * pending_preview  — the most recent notify_update() request that arrived
 *                    while a preview was already in flight.  At most one is
 *                    kept; later arrivals replace earlier ones (coalescing).
 *
 * preview_in_flight — TRUE from the moment we start dispatching a preview
 *                     until the last pending request has been serviced.
 *                     Used to decide whether to coalesce new requests.
 *
 * preview_job_active — TRUE when the current preview dispatched a job to
 *                      the processing thread (async path).  Only set when
 *                      we know a thread job is running for us, so we can
 *                      safely call processing_request_cancel() without
 *                      accidentally aborting an unrelated job.
 */
static update_image *pending_preview   = NULL;
static gboolean      preview_in_flight = FALSE;
static gboolean      preview_job_active = FALSE;
static cmsHPROFILE preview_icc_backup = NULL;
static fits preview_roi_backup;
static fits preview_gfit_backup = { 0 };

/* Forward declaration. */
static void dispatch_preview(update_image *im);

/* Called (on the GTK main thread, via g_idle_add) after every processing-
 * thread job completes.  Also called directly from dispatch_preview() when
 * update_preview_fn ran synchronously (no thread job was submitted).      */
static void on_preview_done(void) {
	preview_job_active = FALSE;

	if (pending_preview) {
		update_image *next = pending_preview;
		pending_preview = NULL;
		/* preview_in_flight stays TRUE — we are immediately starting another. */
		dispatch_preview(next);
	} else {
		preview_in_flight = FALSE;
	}
}

/* GSourceFunc registered with processing_set_post_job_callback().
 * Fires after every processing-thread job, not only preview jobs.
 * The preview_in_flight / preview_job_active guards make it a no-op when
 * an unrelated job finishes while no preview is in progress.              */
static gboolean preview_post_job_callback(gpointer data G_GNUC_UNUSED) {
	if (preview_in_flight && preview_job_active)
		on_preview_done();
	return G_SOURCE_REMOVE;
}

/* Execute one preview update.  Always called from the GTK main thread.
 * Consumes (frees) im.                                                    */
static void dispatch_preview(update_image *im) {
	lock_roi_mutex();

	if (notify_is_blocked) {
		unlock_roi_mutex();
		free(im);
		preview_in_flight = FALSE;
		return;
	}

	if (im->show_preview) {
		siril_log_debug("update preview\n");
		set_cursor_waiting(TRUE);
		im->update_preview_fn();
	}

	waiting_for_thread(); /* no-op on GTK main thread; kept for safety */
	set_progress_bar_data(NULL, PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	/* Don't gfit_modified_update_gui() here; callers are responsible. */
	unlock_roi_mutex();
	free(im);

	if (processing_is_job_active()) {
		/* update_preview_fn submitted an async job.  The post-job callback
		 * will call on_preview_done() when the worker finishes.           */
		preview_job_active = TRUE;
	} else {
		/* update_preview_fn ran synchronously (no thread job submitted).
		 * Drive the next pending update immediately.                      */
		on_preview_done();
	}
}

void copy_gfit_icc_to_backup() {
	if (!gfit->icc_profile)
		return;
	if (preview_icc_backup)
		cmsCloseProfile(preview_icc_backup);
	preview_icc_backup = copyICCProfile(gfit->icc_profile);
}

static void copy_backup_icc_to_gfit() {
	if (gfit->icc_profile)
		cmsCloseProfile(gfit->icc_profile);
	gfit->icc_profile = copyICCProfile(preview_icc_backup);
}

static void clear_backup_icc() {
	if (preview_icc_backup) {
		cmsCloseProfile(preview_icc_backup);
		preview_icc_backup = NULL;
	}
}

int backup_roi() {
	int retval;
	if ((retval = copyfits(&gui.roi.fit, &preview_roi_backup, CP_ALLOC | CP_COPYA | CP_FORMAT | CP_COPYMASK, -1)))
		siril_log_debug("Image copy error in ROI\n");

	return retval;
}

int restore_roi() {
	int retval;
	if ((retval = copyfits(&preview_roi_backup, &gui.roi.fit, CP_ALLOC | CP_COPYA | CP_FORMAT | CP_COPYMASK, -1)))
		siril_log_debug("Image copy error in ROI\n");

	return retval;
}

void copy_gfit_to_backup() {
	guint64 gfit_size = gfit->rx * gfit->ry * gfit->naxes[2] * gfit->type == DATA_FLOAT ? 4 : 2;
	if (!preview_is_active && (get_available_memory() < (gfit_size * 2))) {
		siril_log_warning(_("Warning: insufficient memory available to create a preview.\n"));
		return;
	}
	// We need the backup to have the mask state copied to it, because image operations start from the backup if a preview is active
	if (copyfits(gfit, &preview_gfit_backup, CP_ALLOC | CP_COPYA | CP_COPYMASK | CP_FORMAT, -1)) {
		siril_log_debug("Image copy error in previews\n");
		return;
	}
	copy_fits_metadata(gfit, &preview_gfit_backup);
	if (!com.script)
		copy_gfit_icc_to_backup();
	if (gui.roi.active && backup_roi()) {
		siril_log_debug("Image copy error in ROI\n");
		return;
	}
	preview_is_active = TRUE;
}

int copy_backup_to_gfit() {
	int retval = 0;
	if (!gfit->data && !gfit->fdata)
		retval = 1;
	else {
		// Restore the mask state too
		if (copyfits(&preview_gfit_backup, gfit, CP_COPYA | CP_COPYMASK, -1)) {
			siril_log_debug("Image copy error in previews\n");
			retval = 1;
		} else if (!com.script) {
			copy_backup_icc_to_gfit();
		}
		if (retval == 0) {
			copy_fits_metadata(&preview_gfit_backup, gfit);
			/* Stats were computed on preview-modified pixels; invalidate them
			 * so stale values are never saved back to the input sequence. */
			invalidate_stats_from_fit(gfit);
		}
		if (gui.roi.active && restore_roi()) {
			siril_log_debug("Image copy error in ROI\n");
			retval = 1;
		}
	}
	return retval;
}

fits *get_preview_gfit_backup() {
	return (is_preview_active()) ? &preview_gfit_backup : gfit;
}

fits *get_roi_backup() {
	return (is_preview_active()) ? &preview_roi_backup : &gui.roi.fit;
}

gboolean is_preview_active() {
	return preview_is_active;
}

void clear_backup() {
	clearfits(&preview_gfit_backup);
	clear_backup_icc();
	preview_is_active = FALSE;
}

void set_notify_block(gboolean value) {
	notify_is_blocked = value;
}

void cancel_pending_update(void) {
	if (pending_preview) {
		free(pending_preview);
		pending_preview = NULL;
	}
}

void siril_preview_hide() {
	copy_backup_to_gfit();
	clear_backup();
	notify_gfit_data_modified();  // remap Cairo buffers to the restored image
	gfit_modified_update_gui();
}

void notify_update(gpointer user_data) {
	/* One-time registration of the post-job callback with the processing
	 * thread.  Done here rather than at application startup so that the
	 * preview module is self-contained and requires no external init call. */
	static gboolean callback_registered = FALSE;
	if (!callback_registered) {
		processing_set_post_job_callback(preview_post_job_callback);
		callback_registered = TRUE;
	}

	if (notify_is_blocked) {
		free(user_data);
		return;
	}

	if (preview_in_flight) {
		/* A preview is already in progress.  Coalesce: keep only the latest
		 * request and cancel the running thread job (if any) so the worker
		 * finishes sooner and we can start the new one.                    */
		if (pending_preview)
			free(pending_preview);
		pending_preview = (update_image *) user_data;
		if (preview_job_active)
			processing_request_cancel();
		return;
	}

	/* No preview in flight: start immediately. */
	preview_in_flight = TRUE;
	dispatch_preview((update_image *) user_data);
}
