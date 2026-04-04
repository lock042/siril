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
#ifndef SRC_GUI_SIRIL_PREVIEW_H_
#define SRC_GUI_SIRIL_PREVIEW_H_

#include "core/siril.h"

typedef struct update_preview_struct {
	gboolean show_preview;
	int (*update_preview_fn)(void);
} update_image;

void copy_gfit_icc_to_backup();
int backup_roi();
int restore_roi();
void copy_gfit_to_backup();
int copy_backup_to_gfit();
fits *get_roi_backup();
fits *get_preview_gfit_backup();
gboolean is_preview_active();
void clear_backup();

/* Group-member backup for FLIS layer-group preview.
 *
 * When a group operation runs with for_preview=TRUE, the non-active group
 * members need the same backup/restore treatment that copy_gfit_to_backup /
 * copy_backup_to_gfit gives the active layer (gfit).  These functions store
 * pixel copies for an array of fits* targets (the lay->fit of each non-active
 * member) using only fits* — no FLIS type awareness is required here.
 *
 * copy_backup_to_gfit() calls restore_flis_group_members_from_backup() so
 * the restore is always paired with the gfit restore.  clear_backup() calls
 * free_flis_group_backups() to release memory once the backup is consumed.
 */
int      copy_flis_group_members_to_backup(fits **member_fits, int n_members);
void     restore_flis_group_members_from_backup(void);
void     free_flis_group_backups(void);
gboolean is_flis_group_backup_active(void);

void siril_preview_hide();

void cancel_pending_update();
void set_notify_block(gboolean value);
void notify_update(gpointer user_data);

#endif /* SRC_GUI_SIRIL_PREVIEW_H_ */
