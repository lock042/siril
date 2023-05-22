/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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

#include "siril_preview.h"

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/undo.h"
#include "gui/progress_and_log.h"
#include "gui/image_display.h"
#include "gui/registration_preview.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"


#define PREVIEW_DELAY 200

static guint timer_id;
static gboolean notify_is_blocked;
static gboolean preview_is_active;
static fits preview_gfit_backup = { 0 };

static gboolean update_preview(gpointer user_data) {
	update_image *im = (update_image*) user_data;

	if (notify_is_blocked) return FALSE;

	if (im->show_preview) {
		siril_debug_print("update preview\n");
		set_cursor_waiting(TRUE);
		im->update_preview_fn();
	}

	waiting_for_thread(); // in case function is run in another thread
	set_progress_bar_data(NULL, PROGRESS_DONE);
	if (im->show_preview) {
		notify_gfit_modified();
	}
	return FALSE;
}

static void free_struct(gpointer user_data) {
	update_image *im = (update_image*) user_data;

	timer_id = 0;
	free(im);
}

void copy_gfit_to_backup() {
	if (copyfits(&gfit, &preview_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
		siril_log_message(_("Image copy error in previews\n"));
		return;
	}
//	preview_gfit_backup.icc_profile = copyICCProfile(gfit.icc_profile);
	preview_is_active = TRUE;
}

int copy_backup_to_gfit() {
	int retval = 0;
	if (!gfit.data && !gfit.fdata)
		retval = 1;
	else if (copyfits(&preview_gfit_backup, &gfit, CP_COPYA, -1)) {
		siril_log_message(_("Image copy error in previews\n"));
		retval = 1;
	}
//	gfit.icc_profile = copyICCProfile(preview_gfit_backup.icc_profile);
	return retval;
}

fits *get_preview_gfit_backup() {
	return (is_preview_active()) ? &preview_gfit_backup : &gfit;
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
