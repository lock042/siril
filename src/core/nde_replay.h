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
#ifndef _NDE_REPLAY_H_
#define _NDE_REPLAY_H_

/**
 * \file nde_replay.h
 * \brief NDE replay engine (phase 2, nde-phase2-3-plan.md P2.D).
 *
 * Rebuilds an item's pixels headlessly: baseline checkpoint → each Tier-A
 * record deserialized and applied through generic_image_worker with the
 * nde_replay flag, always into a PRIVATE fits (never gfit).  Phase 2 uses
 * it for verification (flis_replay_check); phase 3's amend-and-replay
 * feeds an edited chain through the same driver.
 *
 * Threading: nde_chain_replay() MUST run inside a processing job (the
 * caller owns the slot — one continuous occupancy for the whole chain, so
 * python cannot interleave).  nde_chain_build() only snapshots and
 * validates; it is safe anywhere.
 */

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ffit;

typedef struct {
	gint       item_id;    /* layer item_id; -1 = plain single image */
	GPtrArray *records;    /* chain members: nde_record* views into ->snapshot */
	GPtrArray *snapshot;   /* owns the records (live-prefix deep copy) */
	gboolean   replayable; /* TRUE when reasons is empty */
	GPtrArray *reasons;    /* gchar*: why the chain cannot be replayed */
} nde_chain;

/**
 * Snapshot the live history and extract + validate the replay chain for
 * @item_id.  Membership: LAYER-scope records targeting the item; CANVAS-
 * scope records on plain images (the image IS the layer).  Blockers:
 * Tier-B records, active-mask records, unknown ops / missing
 * deserializers / newer-version params / unparsable params, canvas-scope
 * records targeting the item on a FLIS document (their layer-offset side
 * effects cannot be verified on a scratch fits yet), merge-down/flatten
 * that replaced the item's pixels, and a missing baseline checkpoint.
 * DOCUMENT-scope property/structure records are ignored (they do not
 * touch the item's pixels; layer.add is embodied in the baseline).
 * Never returns NULL.
 */
nde_chain *nde_chain_build(gint item_id);

void nde_chain_free(nde_chain *chain);

/**
 * Replay the chain from the item's baseline into a freshly allocated fits
 * (caller clearfits()+free()s).  NULL on failure with a heap message in
 * @err (caller g_free()s).  Honours cancellation between records
 * (processing_should_continue()).  Job-slot contract in the file header.
 */
struct ffit *nde_chain_replay(const nde_chain *chain, gchar **err);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_REPLAY_H_ */
