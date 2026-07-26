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

#ifndef SIRIL_CORE_NDE_SCRIPT_SCOPE_H
#define SIRIL_CORE_NDE_SCRIPT_SCOPE_H

#include <glib.h>

/*
 * NDE provenance scope for a Python script run (phase 5, Tier-C capture half —
 * nde-phase5-plan.md "Engine: the script provenance scope").
 *
 * Brackets a script's run and coalesces its whole NET EFFECT on the loaded
 * image into AT MOST ONE nde record, decided at scope end:
 *   1. a valid replay declaration (record_args called, script hash known)
 *        -> Tier-C replayable record (script path + sha256 + args), post-script
 *           pixels as the output checkpoint;
 *   2. else a net change to the image (pixels differ from scope start, or a
 *        non-pixel mutation was flagged) -> Tier-B opaque barrier;
 *   3. else (read-only, or previewed-then-reverted to the start) -> nothing.
 *
 * Decision is by net effect, never per set_pixeldata call: a script that
 * previews via set_image_pixeldata then reverts on cancel ends at its starting
 * pixels and leaves no record.  While a scope is active, per-op NDE capture in
 * generic_image_worker is suppressed — the script record subsumes it.
 */

/* Open a scope for a launching script.  @script_path is the entry file (NULL
 * for an inline -c script: no hash, so Tier-C is impossible — only Tier-B or
 * nothing).  Snapshots the pre-script pixel hash and ensures the chain
 * baseline.  Nested/concurrent begins flatten onto the first (depth-counted). */
void     nde_script_scope_begin(const char *script_path);

/* Close the scope: run the net-effect decision, capture the record (if any),
 * and restore normal per-op capture.  Balances nde_script_scope_begin. */
void     nde_script_scope_end(void);

/* TRUE while a scope is open — the gate generic_image_worker and the
 * set_pixeldata handler consult to suppress their own per-op capture. */
gboolean nde_script_scope_active(void);

/* Mutation entry points flag the scope so the end-of-scope compare runs.
 * _pixels_dirty: a pixel-writing op ran (the pixel hash will confirm net
 * effect).  _nonpixel_dirty: a mutation the pixel hash cannot see (mask, and
 * later metadata) — forces a record regardless of the pixel comparison. */
void     nde_script_scope_mark_pixels_dirty(void);
void     nde_script_scope_mark_nonpixel_dirty(void);

/* Replay declaration (wired to CMD_DECLARE_REPLAYABLE / _RECORD_REPLAY_ARGS in
 * phase 5 step 2).  declare_replayable stashes the argparse schema; record_args
 * marks the script's committed intent and stashes the recorded arg values (a kv
 * blob).  A well-behaved script calls record_args only on Apply, never on a
 * preview or cancel — so the declaration IS the commit signal. */
void     nde_script_scope_declare_replayable(const char *schema);
void     nde_script_scope_record_args(const char *args_kv);

#endif /* SIRIL_CORE_NDE_SCRIPT_SCOPE_H */
