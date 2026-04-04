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

/* ---- FLIS group-member backup ----------------------------------------
 * Stores pixel-data copies for the non-active members of a group that is
 * being processed with live preview.  The active layer (gfit) is handled
 * by the existing preview_gfit_backup; these entries cover the rest.
 * Uses only fits* — no FLIS type knowledge required here.
 * -------------------------------------------------------------------- */
typedef struct {
	fits *target;   /* lay->fit pointer — NOT owned by this module */
	fits  backup;   /* owned pixel copy */
} flis_member_backup_t;

static flis_member_backup_t *flis_member_backups  = NULL;
static int                   n_flis_member_backups = 0;

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

int copy_flis_group_members_to_backup(fits **member_fits, int n_members) {
	if (!member_fits || n_members <= 0)
		return 0;
	flis_member_backup_t *arr = calloc(n_members, sizeof(flis_member_backup_t));
	if (!arr) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	int saved = 0;
	for (int i = 0; i < n_members; i++) {
		if (!member_fits[i])
			continue;
		arr[saved].target = member_fits[i];
		memset(&arr[saved].backup, 0, sizeof(fits));
		if (copyfits(member_fits[i], &arr[saved].backup,
		             CP_ALLOC | CP_FORMAT | CP_COPYA | CP_COPYMASK, -1)) {
			siril_debug_print("copy_flis_group_members_to_backup: copy error for member %d\n", i);
			/* Free what we have so far and abort */
			for (int j = 0; j < saved; j++)
				clearfits(&arr[j].backup);
			free(arr);
			return 1;
		}
		copy_fits_metadata(member_fits[i], &arr[saved].backup);
		saved++;
	}
	/* Free any pre-existing group backup before replacing */
	free_flis_group_backups();
	flis_member_backups  = arr;
	n_flis_member_backups = saved;
	return 0;
}

void restore_flis_group_members_from_backup(void) {
	if (!flis_member_backups || n_flis_member_backups <= 0)
		return;
	for (int i = 0; i < n_flis_member_backups; i++) {
		fits *tgt = flis_member_backups[i].target;
		if (!tgt) continue;
		if (copyfits(&flis_member_backups[i].backup, tgt, CP_COPYA | CP_COPYMASK, -1))
			siril_debug_print("restore_flis_group_members_from_backup: copy error for member %d\n", i);
		else
			copy_fits_metadata(&flis_member_backups[i].backup, tgt);
	}
}

void free_flis_group_backups(void) {
	if (!flis_member_backups)
		return;
	for (int i = 0; i < n_flis_member_backups; i++)
		clearfits(&flis_member_backups[i].backup);
	free(flis_member_backups);
	flis_member_backups   = NULL;
	n_flis_member_backups = 0;
}

gboolean is_flis_group_backup_active(void) {
	return flis_member_backups != NULL;
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
		siril_debug_print("Image copy error in ROI\n");

	return retval;
}

int restore_roi() {
	int retval;
	if ((retval = copyfits(&preview_roi_backup, &gui.roi.fit, CP_ALLOC | CP_COPYA | CP_FORMAT | CP_COPYMASK, -1)))
		siril_debug_print("Image copy error in ROI\n");

	return retval;
}

void copy_gfit_to_backup() {
	guint64 gfit_size = gfit->rx * gfit->ry * gfit->naxes[2] * gfit->type == DATA_FLOAT ? 4 : 2;
	if (!preview_is_active && (get_available_memory() < (gfit_size * 2))) {
		siril_log_color_message(_("Warning: insufficient memory available to create a preview.\n"), "salmon");
		return;
	}
	// We need the backup to have the mask state copied to it, because image operations start from the backup if a preview is active
	if (copyfits(gfit, &preview_gfit_backup, CP_ALLOC | CP_COPYA | CP_COPYMASK | CP_FORMAT, -1)) {
		siril_debug_print("Image copy error in previews\n");
		return;
	}
	copy_fits_metadata(gfit, &preview_gfit_backup);
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
	if (!gfit->data && !gfit->fdata)
		retval = 1;
	else {
		// Restore the mask state too
		if (copyfits(&preview_gfit_backup, gfit, CP_COPYA | CP_COPYMASK, -1)) {
			siril_debug_print("Image copy error in previews\n");
			retval = 1;
		} else if (!com.script) {
			copy_backup_icc_to_gfit();
		}
		if (retval == 0) copy_fits_metadata(&preview_gfit_backup, gfit);
		if (gui.roi.active && restore_roi()) {
			siril_debug_print("Image copy error in ROI\n");
			retval = 1;
		}
	}
	/* Restore non-active FLIS group members if a group preview is active.
	 * This is a no-op for plain images and single-layer FLIS operations. */
	restore_flis_group_members_from_backup();
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
	/* Release any FLIS group-member backups.  The restore was already done
	 * (by copy_backup_to_gfit or by the final-apply worker path), so we just
	 * free the memory here without restoring. */
	free_flis_group_backups();
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
