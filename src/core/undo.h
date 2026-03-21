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
#ifndef UNDO_H_
#define UNDO_H_

#define UNDO		-1
#define REDO		 1
#define HISTORY_SIZE 30		// size of the operations history

/* Sentinel value for historic.flis_layer_id indicating that the state was
 * saved from a plain FITS image (not a FLIS layer). */
#define FLIS_UNDO_LAYER_NONE  (-1)

/*
 * IMPORTANT — siril.h change required
 * =====================================
 * The `historic` struct in siril.h must have the following field added:
 *
 *   gint flis_layer_id;
 *
 * This field stores the item_id of the FLIS layer that was active when
 * the undo state was saved.  It is set to FLIS_UNDO_LAYER_NONE (-1) for
 * states saved from plain FITS images.  Because item_ids are stable and
 * never reused within a session, the undo mechanism can locate the correct
 * layer even after the stack has been reordered.
 *
 * Suggested placement immediately after mask_bitpix:
 *
 *   typedef struct {
 *       char  *filename;
 *       char  *mask_filename;
 *       int    mask_bitpix;
 *       gint   flis_layer_id;   ← ADD THIS
 *       int    rx, ry;
 *       ...
 *   } historic;
 */

gboolean is_undo_available(void);
gboolean is_redo_available(void);
int      undo_display_data(int dir);
int      undo_save_state(fits *fit, const char *message, ...);
int      undo_flush(void);
void     set_undo_redo_tooltip(void);
gboolean undo_in_thread(gpointer user_data);
gboolean redo_in_thread(gpointer user_data);

/**
 * flis_undo_purge_layer:
 * @item_id: the stable item_id of the layer being removed.
 *
 * Purges all undo/redo states belonging to the specified layer from the
 * history ring buffer and frees their swap files.  All states for other
 * layers are preserved.
 *
 * Call this from flis_gui.c (on_flis_remove_layer_clicked) immediately
 * BEFORE the layer is destroyed, so the item_id is still valid.
 *
 * Add, duplicate, and reorder operations do NOT need to call this:
 *   - Add / duplicate: the new layer has no prior states.
 *   - Reorder: layer lookup uses item_id not position, so reordering is
 *     transparent to the undo mechanism.
 */
void flis_undo_purge_layer(gint item_id);

#endif /* UNDO_H_ */
