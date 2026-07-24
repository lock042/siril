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
 * NDE baseline checkpoint store — see nde_checkpoint.h for the contract and
 * flis-nde-sketch.md §1/§7 for the design.  Baselines are pre-op pixel
 * snapshots kept in delete-on-close swap files, keyed by layer item_id.
 */

#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#ifdef _WIN32
#include <windows.h>
#include <io.h>  /* _get_osfhandle */
#endif

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/nde_checkpoint.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "algos/siril_wcs.h"

/* ======================================================================= */
/* Factored swap-file storage (flis-nde-sketch.md §7 decision 13)          */
/*                                                                         */
/* A clean pair writing/reading a fits' pixel payload + geometry to/from a */
/* delete-on-close swap file, kept deliberately independent of the         */
/* checkpoint registry so it can later become the shared snapshot store.   */
/* Mirrors undo.c's undo_build_swapfile / undo_get_data_* (the statics     */
/* there are not exported; this is the factored equivalent).               */
/* ======================================================================= */

/* One stored pixel snapshot: an open, delete-on-close fd plus the geometry
 * and metadata needed to reconstruct an equivalent fits on read. */
typedef struct {
	int        fd;            /* open swap-file fd; delete-on-close */
	guint      rx, ry;
	guint      nchans;
	data_type  type;
	wcs_info   wcsdata;
	wcsprm_t  *wcslib;        /* deep copy; freed with the snapshot */
	double     focal_length;
} swap_snapshot;

/* Mark @fd delete-on-close so the swap file vanishes when the fd closes,
 * even on abnormal termination (mirrors undo.c's swap_mark_delete_on_close). */
static void swap_mark_delete_on_close(int fd, const gchar *path) {
#ifndef _WIN32
	(void)fd;
	(void) g_unlink(path);
#else
	(void)path;
	HANDLE h = (HANDLE)_get_osfhandle(fd);
	if (h != INVALID_HANDLE_VALUE) {
		FILE_DISPOSITION_INFO fdi = { .DeleteFile = TRUE };
		(void) SetFileInformationByHandle(h, FileDispositionInfo, &fdi, sizeof(fdi));
	}
#endif
}

/* Deep-copy @src's pixels + geometry into a fresh delete-on-close swap file.
 * Returns a heap swap_snapshot (free with swap_snapshot_free) or NULL on
 * error.  Performs file I/O only — no lock is taken here. */
static swap_snapshot *swap_snapshot_write(const fits *src) {
	if (!src)
		return NULL;
	if (src->type != DATA_USHORT && src->type != DATA_FLOAT)
		return NULL;
	if ((src->type == DATA_USHORT && !src->data) ||
	    (src->type == DATA_FLOAT && !src->fdata))
		return NULL;

	gchar *nameBuff = g_build_filename(com.pref.swap_dir, "siril_ndeb-XXXXXX", NULL);
#ifndef _WIN32
	int fd = g_mkstemp(nameBuff);
#else
	int fd = g_mkstemp_full(nameBuff, _O_RDWR | _O_CREAT | _O_EXCL | _O_BINARY | _O_TEMPORARY,
	    _S_IREAD | _S_IWRITE);
#endif
	if (fd < 0) {
		siril_log_error(_("File I/O Error: Unable to create checkpoint swap file in %s: [%s]\n"),
		                com.pref.swap_dir, strerror(errno));
		g_free(nameBuff);
		return NULL;
	}
	swap_mark_delete_on_close(fd, nameBuff);
	g_free(nameBuff);

	size_t n = (size_t)src->naxes[0] * src->naxes[1] * src->naxes[2];
	errno = 0;
	ssize_t written = -1;
	if (src->type == DATA_USHORT)
		written = write(fd, src->data, n * sizeof(WORD));
	else
		written = write(fd, src->fdata, n * sizeof(float));
	if (written < 0) {
		siril_log_error(_("File I/O Error: Unable to write checkpoint swap file: [%s]\n"),
		                strerror(errno));
		g_close(fd, NULL);
		return NULL;
	}

	swap_snapshot *s = g_new0(swap_snapshot, 1);
	s->fd           = fd;
	s->rx           = src->rx;
	s->ry           = src->ry;
	s->nchans       = (guint)src->naxes[2];
	s->type         = src->type;
	s->wcsdata      = src->keywords.wcsdata;
	s->wcslib       = NULL;
	if (src->keywords.wcslib) {
		int wstatus = -1;
		s->wcslib = wcs_deepcopy(src->keywords.wcslib, &wstatus);
		if (wstatus)
			siril_log_debug("nde checkpoint: could not copy wcslib struct\n");
	}
	s->focal_length = src->keywords.focal_length;
	return s;
}

/* Reconstruct a fully-owned fits from @s (caller clearfits()+frees).  Reads
 * the swap file back bit-exactly.  File I/O only — no lock is taken. */
static fits *swap_snapshot_read(const swap_snapshot *s) {
	if (!s)
		return NULL;
	fits *out = NULL;
	if (new_fit_image(&out, (int)s->rx, (int)s->ry, (int)s->nchans, s->type))
		return NULL;

	size_t n = (size_t)s->rx * s->ry * s->nchans;
	size_t bytes = n * (s->type == DATA_USHORT ? sizeof(WORD) : sizeof(float));
	void *dst = (s->type == DATA_USHORT) ? (void *)out->data : (void *)out->fdata;

	if (lseek(s->fd, 0, SEEK_SET) == (off_t)-1) {
		siril_log_error(_("nde checkpoint: error seeking swap file: [%s]\n"), strerror(errno));
		clearfits(out);
		free(out);
		return NULL;
	}
	errno = 0;
	if (read(s->fd, dst, bytes) < (ssize_t)bytes) {
		siril_log_error(_("nde checkpoint: read failed with error [%s]\n"), strerror(errno));
		clearfits(out);
		free(out);
		return NULL;
	}

	memcpy(&out->keywords.wcsdata, &s->wcsdata, sizeof(wcs_info));
	free_wcs(out);
	if (s->wcslib) {
		int wstatus = -1;
		out->keywords.wcslib = wcs_deepcopy(s->wcslib, &wstatus);
		if (wstatus)
			siril_log_debug("nde checkpoint: could not copy wcslib struct\n");
	} else {
		reset_wcsdata(out);
	}
	out->keywords.focal_length = s->focal_length;
	return out;
}

static void swap_snapshot_free(swap_snapshot *s) {
	if (!s)
		return;
	if (s->fd >= 0)
		g_close(s->fd, NULL);   /* delete-on-close removes the file */
	if (s->wcslib) {
		wcsfree(s->wcslib);
		free(s->wcslib);
	}
	g_free(s);
}

/* ======================================================================= */
/* Registry (LEAF lock — never held across the swap I/O above)             */
/* ======================================================================= */

/* item_id (GINT_TO_POINTER) -> swap_snapshot*.  Guarded by cp_mutex, which
 * is a strict leaf: swap_snapshot_write/read/free (file I/O, fits copies)
 * all run OUTSIDE the lock; only the hash-table pointer op is inside. */
static GMutex      cp_mutex;
static GHashTable *cp_table;   /* created lazily */

static void cp_table_ensure_locked(void) {
	if (!cp_table)
		cp_table = g_hash_table_new_full(g_direct_hash, g_direct_equal,
		                                 NULL, (GDestroyNotify)swap_snapshot_free);
}

/* Output-checkpoint table (phase-4 barriers): record_id -> post-op pixels.
 * Declared here because drop/purge below cover both tables. */
typedef struct {
	swap_snapshot *snap;
	gint item_id;      /* owning layer, for drop-by-layer */
} output_ckpt;

/* gint64* (owned) -> output_ckpt*.  Shares cp_mutex (still a strict leaf:
 * swap I/O runs outside it, exactly like the baseline table). */
static GHashTable *out_table;

static void output_ckpt_free(gpointer p) {
	output_ckpt *c = p;
	swap_snapshot_free(c->snap);
	g_free(c);
}

static void out_table_ensure_locked(void) {
	if (!out_table)
		out_table = g_hash_table_new_full(g_int64_hash, g_int64_equal,
		                                  g_free, output_ckpt_free);
}


void nde_checkpoint_baseline_ensure(const fits *pre, gint item_id) {
	if (!pre)
		return;
	/* Fast path: bail without any copy if a baseline already exists.  The
	 * table lookup takes the leaf lock only. */
	g_mutex_lock(&cp_mutex);
	gboolean exists = cp_table &&
	    g_hash_table_contains(cp_table, GINT_TO_POINTER(item_id));
	g_mutex_unlock(&cp_mutex);
	if (exists)
		return;

	/* Prepare the (expensive) swap copy OUTSIDE the lock. */
	swap_snapshot *s = swap_snapshot_write(pre);
	if (!s)
		return;

	/* Re-check under the lock: another thread may have won the race between
	 * our exists check and here.  First writer keeps its snapshot; the loser
	 * frees its own — the baseline must reflect the pre-FIRST-op pixels. */
	g_mutex_lock(&cp_mutex);
	cp_table_ensure_locked();
	if (g_hash_table_contains(cp_table, GINT_TO_POINTER(item_id))) {
		g_mutex_unlock(&cp_mutex);
		swap_snapshot_free(s);
		return;
	}
	g_hash_table_insert(cp_table, GINT_TO_POINTER(item_id), s);
	g_mutex_unlock(&cp_mutex);
}

void nde_checkpoint_baseline_adopt(const fits *src, gint item_id) {
	if (!src)
		return;
	/* Prepare the copy outside the lock. */
	swap_snapshot *s = swap_snapshot_write(src);
	if (!s)
		return;
	g_mutex_lock(&cp_mutex);
	cp_table_ensure_locked();
	/* Overwrite: on load the on-disk baseline is authoritative.  The
	 * hash table's destroy-func frees any previous snapshot. */
	g_hash_table_insert(cp_table, GINT_TO_POINTER(item_id), s);
	g_mutex_unlock(&cp_mutex);
}

fits *nde_checkpoint_baseline_get(gint item_id) {
	/* Take a private handle to the snapshot under the lock, then do the
	 * read OUTSIDE it.  The snapshot cannot be freed underneath us because
	 * purge/drop only run on the same threads that mutate the document (no
	 * concurrent free with a live get in the single-slot model), and the
	 * read only touches the fd + immutable geometry. */
	g_mutex_lock(&cp_mutex);
	swap_snapshot *s = cp_table ?
	    g_hash_table_lookup(cp_table, GINT_TO_POINTER(item_id)) : NULL;
	g_mutex_unlock(&cp_mutex);
	if (!s)
		return NULL;
	return swap_snapshot_read(s);
}

gboolean nde_checkpoint_baseline_exists(gint item_id) {
	g_mutex_lock(&cp_mutex);
	gboolean e = cp_table &&
	    g_hash_table_contains(cp_table, GINT_TO_POINTER(item_id));
	g_mutex_unlock(&cp_mutex);
	return e;
}

void nde_checkpoint_drop(gint item_id) {
	g_mutex_lock(&cp_mutex);
	if (cp_table)
		g_hash_table_remove(cp_table, GINT_TO_POINTER(item_id));
	if (out_table) {
		/* output checkpoints of a dying layer die with it */
		GHashTableIter it;
		gpointer k, v;
		g_hash_table_iter_init(&it, out_table);
		while (g_hash_table_iter_next(&it, &k, &v)) {
			output_ckpt *c = v;
			if (c->item_id == item_id)
				g_hash_table_iter_remove(&it);
		}
	}
	g_mutex_unlock(&cp_mutex);
}

void nde_checkpoint_purge(void) {
	GHashTable *doomed = NULL, *doomed_out = NULL;
	g_mutex_lock(&cp_mutex);
	doomed = cp_table;
	cp_table = NULL;
	doomed_out = out_table;
	out_table = NULL;
	g_mutex_unlock(&cp_mutex);
	/* Destroy (which closes fds -> deletes files) outside the lock. */
	if (doomed)
		g_hash_table_destroy(doomed);
	if (doomed_out)
		g_hash_table_destroy(doomed_out);
}

/* ======================================================================= */
/* Output checkpoints (phase-4 barriers) — keyed by record_id              */
/* ======================================================================= */

void nde_checkpoint_output_store(const fits *post, gint64 record_id, gint item_id) {
	if (!post || record_id <= 0)
		return;
	swap_snapshot *s = swap_snapshot_write(post);   /* outside the lock */
	if (!s)
		return;
	output_ckpt *c = g_new0(output_ckpt, 1);
	c->snap = s;
	c->item_id = item_id;
	gint64 *key = g_new(gint64, 1);
	*key = record_id;
	g_mutex_lock(&cp_mutex);
	out_table_ensure_locked();
	g_hash_table_insert(out_table, key, c);   /* overwrite = refresh */
	g_mutex_unlock(&cp_mutex);
}

fits *nde_checkpoint_output_get(gint64 record_id) {
	g_mutex_lock(&cp_mutex);
	output_ckpt *c = out_table ?
	    g_hash_table_lookup(out_table, &record_id) : NULL;
	swap_snapshot *s = c ? c->snap : NULL;
	g_mutex_unlock(&cp_mutex);
	if (!s)
		return NULL;
	return swap_snapshot_read(s);   /* outside the lock, same rationale as
	                                 * nde_checkpoint_baseline_get */
}

gboolean nde_checkpoint_output_exists(gint64 record_id) {
	g_mutex_lock(&cp_mutex);
	gboolean e = out_table && g_hash_table_contains(out_table, &record_id);
	g_mutex_unlock(&cp_mutex);
	return e;
}

void nde_checkpoint_output_drop(gint64 record_id) {
	g_mutex_lock(&cp_mutex);
	if (out_table)
		g_hash_table_remove(out_table, &record_id);
	g_mutex_unlock(&cp_mutex);
}

gint nde_checkpoint_active_item_id(void) {
	if (is_current_image_flis()) {
		flis_layer_t *lay = flis_active_layer();
		if (lay)
			return lay->item_id;
	}
	return -1;
}
