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

#include "siril_preview.h"

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"


#define PREVIEW_DELAY 200

static guint timer_id = 0;
static gboolean notify_is_blocked;
static gboolean preview_is_active;
static cmsHPROFILE preview_icc_backup = NULL;
static fits preview_roi_backup;
static fits preview_gfit_backup = { 0 };

static gboolean update_preview(gpointer user_data) {
	lock_roi_mutex();
	if (notify_is_blocked)
		return FALSE;
	update_image *im = (update_image*) user_data;

	if (im->show_preview) {
		siril_debug_print("update preview\n");
		set_cursor_waiting(TRUE);
		im->update_preview_fn();
	}

	waiting_for_thread(); // in case function is run in another thread
	set_progress_bar_data(NULL, PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	// Don't notify_gfit_modified() here, it must be done by the callers
	unlock_roi_mutex();
	return FALSE;
}

static void free_struct(gpointer user_data) {
	update_image *im = (update_image*) user_data;

	timer_id = 0;
	free(im);
}

void copy_gfit_icc_to_backup() {
	if (!gfit.icc_profile)
		return;
	if (preview_icc_backup)
		cmsCloseProfile(preview_icc_backup);
	preview_icc_backup = copyICCProfile(gfit.icc_profile);
}

static void copy_backup_icc_to_gfit() {
	if (gfit.icc_profile)
		cmsCloseProfile(gfit.icc_profile);
	gfit.icc_profile = copyICCProfile(preview_icc_backup);
}

int backup_roi() {
	int retval;
	if ((retval = copyfits(&gui.roi.fit, &preview_roi_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)))
		siril_debug_print("Image copy error in ROI\n");

	return retval;
}

int restore_roi() {
	int retval;
	if ((retval = copyfits(&preview_roi_backup, &gui.roi.fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)))
		siril_debug_print("Image copy error in ROI\n");

	return retval;
}

void copy_gfit_to_backup() {
	guint64 gfit_size = gfit.rx * gfit.ry * gfit.naxes[2] * gfit.type == DATA_FLOAT ? 4 : 2;
	if (!preview_is_active && (get_available_memory() < (gfit_size * 2))) {
		siril_log_color_message(_("Warning: insufficient memory available to create a preview.\n"), "salmon");
		return;
	}
	if (copyfits(&gfit, &preview_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
		siril_debug_print("Image copy error in previews\n");
		return;
	}
	copy_fits_metadata(&gfit, &preview_gfit_backup);
	if (!com.script)
		copy_gfit_icc_to_backup();
	if (gui.roi.active && backup_roi()) {
		siril_debug_print("Image copy error in ROI\n");
		return;
	}
	preview_is_active = TRUE;
}

int copy_backup_to_gfit() {
	int retval = 0;
	if (!gfit.data && !gfit.fdata)
		retval = 1;
	else {
		if (copyfits(&preview_gfit_backup, &gfit, CP_COPYA, -1)) {
			siril_debug_print("Image copy error in previews\n");
			retval = 1;
		} else if (!com.script) {
			copy_backup_icc_to_gfit();
		}
		if (retval == 0) copy_fits_metadata(&preview_gfit_backup, &gfit);
		if (gui.roi.active && restore_roi()) {
			siril_debug_print("Image copy error in ROI\n");
			retval = 1;
		}
	}
	return retval;
}

fits *get_preview_gfit_backup() {
	return (is_preview_active()) ? &preview_gfit_backup : &gfit;
}

fits *get_roi_backup() {
	return (is_preview_active()) ? &preview_roi_backup : &gui.roi.fit;
}

gboolean is_preview_active() {
	return preview_is_active;
}

void clear_backup() {
	clearfits(&preview_gfit_backup);
	preview_is_active = FALSE;
}

void set_notify_block(gboolean value) {
	notify_is_blocked = value;
}

void cancel_pending_update() {
    if (timer_id != 0) {
        g_source_remove(timer_id);
        timer_id = 0;
    }
}

void siril_preview_hide() {
	copy_backup_to_gfit();
	clear_backup();
	notify_gfit_modified();
}

void notify_update(gpointer user_data) {
	if (timer_id != 0) {
		g_source_remove(timer_id);
	}
	timer_id = g_timeout_add_full(G_PRIORITY_DEFAULT_IDLE,
			PREVIEW_DELAY, (GSourceFunc) update_preview, user_data,
			(GDestroyNotify) free_struct);
}
