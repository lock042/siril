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

/*
 * NDE script provenance scope — see nde_script_scope.h for the contract and
 * nde-phase5-plan.md "Engine: the script provenance scope".
 */

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/nde_history.h"
#include "core/nde_checkpoint.h"
#include "core/nde_script_scope.h"

/* A single global scope guarded by scope_mutex.  Mutations arrive from the
 * python-comm thread (set_pixeldata), the processing worker (cmd ops) and the
 * pyscript / main thread (begin/end), so every field access takes the lock.
 * Concurrent scripts flatten onto the first scope via depth (rare, and the
 * processing-slot reservation already serialises their actual mutations). */
static GMutex scope_mutex;
static struct {
	gboolean active;
	gint     depth;
	gint     item_id;
	gchar   *script_path;
	gchar   *script_sha256;
	gchar   *basename;         /* for the record summary */
	gchar   *baseline_hash;    /* pixel sha256 at scope begin, or NULL */
	gboolean wrote_pixels;     /* a pixel-writing op ran under the scope */
	gboolean nonpixel_dirty;   /* a mutation the pixel hash cannot see (mask) */
	gboolean declared;         /* record_args was called (commit intent) */
	gchar   *decl_args;        /* recorded argparse values (kv blob) */
	gchar   *decl_schema;      /* argparse schema (informational) */
} sc;

/* SHA-256 of gfit's pixel buffer — mirrors flis_pixel_sha256 (static there).
 * Caller need not hold the lock on gfit; we take a brief reader lock. */
static gchar *scope_pixel_hash(void) {
	if (!gfit)
		return NULL;
	g_rw_lock_reader_lock(&gfit->rwlock);
	gchar *hex = NULL;
	size_t npix = (size_t)gfit->rx * gfit->ry *
	              (gfit->naxes[2] ? gfit->naxes[2] : 1);
	const guchar *buf = NULL;
	gsize nbytes = 0;
	if (gfit->type == DATA_FLOAT) {
		buf = (const guchar *)gfit->fdata;
		nbytes = npix * sizeof(float);
	} else {
		buf = (const guchar *)gfit->data;
		nbytes = npix * sizeof(WORD);
	}
	if (buf && npix) {
		GChecksum *ck = g_checksum_new(G_CHECKSUM_SHA256);
		g_checksum_update(ck, buf, nbytes);
		hex = g_strdup(g_checksum_get_string(ck));
		g_checksum_free(ck);
	}
	g_rw_lock_reader_unlock(&gfit->rwlock);
	return hex;
}

static void scope_reset_locked(void) {
	sc.active = FALSE;
	sc.depth = 0;
	sc.item_id = -1;
	g_clear_pointer(&sc.script_path, g_free);
	g_clear_pointer(&sc.script_sha256, g_free);
	g_clear_pointer(&sc.basename, g_free);
	g_clear_pointer(&sc.baseline_hash, g_free);
	sc.wrote_pixels = FALSE;
	sc.nonpixel_dirty = FALSE;
	sc.declared = FALSE;
	g_clear_pointer(&sc.decl_args, g_free);
	g_clear_pointer(&sc.decl_schema, g_free);
}

void nde_script_scope_begin(const char *script_path) {
	/* Snapshot pixels and hash the entry file BEFORE taking the lock (both may
	 * be slow-ish; nothing else can open a scope meanwhile because only script
	 * launches call begin, and they are serialised at the launch site). */
	gchar *sha256 = script_path ? nde_file_sha256(script_path, NULL) : NULL;
	gchar *base = script_path ? g_path_get_basename(script_path) : NULL;
	gchar *baseline = scope_pixel_hash();
	gint item = nde_checkpoint_active_item_id();

	g_mutex_lock(&scope_mutex);
	if (sc.active) {
		/* Nested / concurrent begin: flatten onto the outermost scope. */
		sc.depth++;
		g_mutex_unlock(&scope_mutex);
		g_free(sha256);
		g_free(base);
		g_free(baseline);
		return;
	}
	scope_reset_locked();
	sc.active = TRUE;
	sc.depth = 1;
	sc.item_id = item;
	sc.script_path = g_strdup(script_path);
	sc.script_sha256 = sha256;    /* ownership transferred */
	sc.basename = base;           /* ownership transferred */
	sc.baseline_hash = baseline;  /* ownership transferred */
	g_mutex_unlock(&scope_mutex);

	/* Ensure the chain baseline exists for the pre-script state so records
	 * that precede the script (if any) remain replayable.  Idempotent. */
	if (gfit && gfit->rx) {
		g_rw_lock_reader_lock(&gfit->rwlock);
		nde_checkpoint_baseline_ensure(gfit, item);
		g_rw_lock_reader_unlock(&gfit->rwlock);
	}
}

gboolean nde_script_scope_active(void) {
	g_mutex_lock(&scope_mutex);
	gboolean r = sc.active;
	g_mutex_unlock(&scope_mutex);
	return r;
}

void nde_script_scope_mark_pixels_dirty(void) {
	g_mutex_lock(&scope_mutex);
	if (sc.active)
		sc.wrote_pixels = TRUE;
	g_mutex_unlock(&scope_mutex);
}

void nde_script_scope_mark_nonpixel_dirty(void) {
	g_mutex_lock(&scope_mutex);
	if (sc.active)
		sc.nonpixel_dirty = TRUE;
	g_mutex_unlock(&scope_mutex);
}

void nde_script_scope_declare_replayable(const char *schema) {
	g_mutex_lock(&scope_mutex);
	if (sc.active) {
		g_free(sc.decl_schema);
		sc.decl_schema = g_strdup(schema);
	}
	g_mutex_unlock(&scope_mutex);
}

void nde_script_scope_record_args(const char *args_kv) {
	g_mutex_lock(&scope_mutex);
	if (sc.active) {
		g_free(sc.decl_args);
		sc.decl_args = g_strdup(args_kv);
		sc.declared = TRUE;   /* commit intent */
	}
	g_mutex_unlock(&scope_mutex);
}

void nde_script_scope_end(void) {
	g_mutex_lock(&scope_mutex);
	if (!sc.active) {
		g_mutex_unlock(&scope_mutex);
		return;
	}
	if (--sc.depth > 0) {          /* still inside a nested begin */
		g_mutex_unlock(&scope_mutex);
		return;
	}
	/* Take ownership of the scope state and close it under the lock. */
	gint     item          = sc.item_id;
	gchar   *baseline_hash  = sc.baseline_hash;   sc.baseline_hash = NULL;
	gboolean wrote_pixels   = sc.wrote_pixels;
	gboolean nonpixel_dirty = sc.nonpixel_dirty;
	gboolean declared_intent = sc.declared;
	gboolean declared       = sc.declared && sc.script_sha256 != NULL;
	gchar   *script_path    = sc.script_path;     sc.script_path = NULL;
	gchar   *script_sha256  = sc.script_sha256;   sc.script_sha256 = NULL;
	gchar   *decl_args      = sc.decl_args;        sc.decl_args = NULL;
	gchar   *basename       = sc.basename;         sc.basename = NULL;
	scope_reset_locked();
	g_mutex_unlock(&scope_mutex);

	/* Perf gate: a purely read-only script mutated nothing → no record, and no
	 * need to hash the (possibly large) final image. */
	if (!wrote_pixels && !nonpixel_dirty)
		goto done;

	/* Net-effect test: compare the post-script pixels to the scope-start hash.
	 * A previewed-then-reverted script writes pixels but ends where it began,
	 * so its hash matches → no record. */
	gboolean pixels_changed = FALSE;
	if (wrote_pixels) {
		gchar *final_hash = scope_pixel_hash();
		pixels_changed = baseline_hash && final_hash &&
		                 g_strcmp0(baseline_hash, final_hash) != 0;
		g_free(final_hash);
	}
	gboolean changed = pixels_changed || nonpixel_dirty;
	if (!changed)
		goto done;

	const gchar *summary = basename ? basename : _("Python script");

	/* The script declared replay support but the launch gave us no stable file
	 * to hash (inline -c, or an unsaved editor buffer run via a temp file) —
	 * the recipe would be unreproducible, so the declaration is dropped.  Say
	 * so, or the silent Tier-B downgrade is baffling. */
	if (declared_intent && !declared)
		siril_log_message(_("Script '%s' declared replay support but was not run from a saved script file — recording an opaque history step\n"),
		                  summary);

	if (declared) {
		/* Tier-C: a replayable record carrying script path + sha256 + args. */
		GString *kv = nde_kv_start();
		nde_kv_add_str(kv, "script", script_path ? script_path : "");
		nde_kv_add_str(kv, "sha256", script_sha256 ? script_sha256 : "");
		if (decl_args)
			nde_kv_add_str(kv, "args", decl_args);
		gchar *params = nde_kv_end(kv);   /* ownership passed to capture */
		g_rw_lock_reader_lock(&gfit->rwlock);
		nde_capture_script("python.script", NDE_SCOPE_LAYER, item,
		                   params, summary, gfit);
		g_rw_lock_reader_unlock(&gfit->rwlock);
	} else {
		/* Tier-B opaque barrier — as python.set_pixeldata was per call, but now
		 * once for the whole script's net effect. */
		g_rw_lock_reader_lock(&gfit->rwlock);
		nde_capture_opaque("python.script", NDE_SCOPE_LAYER, item,
		                   summary, gfit);
		g_rw_lock_reader_unlock(&gfit->rwlock);
	}

done:
	g_free(baseline_hash);
	g_free(script_path);
	g_free(script_sha256);
	g_free(decl_args);
	g_free(basename);
}
