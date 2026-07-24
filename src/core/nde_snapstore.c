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
 * Shared refcounted pixel-snapshot store — see nde_snapstore.h for the
 * contract and nde-convergence-plan.md C1 for the design.  The swap-file
 * payload code was factored out of nde_checkpoint.c (which now stores
 * nde_snap references), which in turn mirrored undo.c's swap machinery.
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
#include "core/nde_snapstore.h"
#include "io/image_format_fits.h"
#include "algos/statistics.h"
#include "algos/siril_wcs.h"

typedef struct wcsprm wcsprm_t;

typedef struct {
	gint     item_id;
	gint64   record_id;
	gboolean post;
} snap_tag;

struct nde_snap {
	gint       refcount;      /* atomic */
	int        fd;            /* open swap-file fd; delete-on-close */
	guint      rx, ry;
	guint      nchans;
	data_type  type;
	gint64     bytes;         /* pixel payload size on disk */
	wcs_info   wcsdata;
	wcsprm_t  *wcslib;        /* deep copy; freed with the snapshot */
	double     focal_length;
	gboolean   tagged;
	snap_tag   tag;
};

/* ======================================================================= */
/* Store state (LEAF lock — never held across swap I/O or fits copies)     */
/* ======================================================================= */

static GMutex store_mutex;
static GHashTable *registry;       /* snap_tag* (owned) -> nde_snap* (NO ref) */
static GQueue pool = G_QUEUE_INIT; /* nde_snap* (pool holds a ref); head = MRU */
static gint64 pool_bytes;
static nde_snapstore_stats_t stats;

static guint tag_hash(gconstpointer p) {
	const snap_tag *t = p;
	return g_int_hash(&t->item_id) ^ g_int64_hash(&t->record_id) ^ (t->post ? 0x9e3779b9u : 0);
}

static gboolean tag_equal(gconstpointer a, gconstpointer b) {
	const snap_tag *x = a, *y = b;
	return x->item_id == y->item_id && x->record_id == y->record_id && x->post == y->post;
}

static void registry_ensure_locked(void) {
	if (!registry)
		registry = g_hash_table_new_full(tag_hash, tag_equal, g_free, NULL);
}

/* Remove @s's registry entry if it is still the indexed snapshot for its
 * tag (a later same-tag snapshot may have replaced it).  Mutex held. */
static void deregister_locked(nde_snap *s) {
	if (!s->tagged || !registry)
		return;
	nde_snap *cur = g_hash_table_lookup(registry, &s->tag);
	if (cur == s)
		g_hash_table_remove(registry, &s->tag);
	s->tagged = FALSE;
}

/* ======================================================================= */
/* Swap payload (I/O only — no lock)                                       */
/* ======================================================================= */

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

nde_snap *nde_snap_create(const fits *src) {
	if (!src)
		return NULL;
	if (src->type != DATA_USHORT && src->type != DATA_FLOAT)
		return NULL;
	if ((src->type == DATA_USHORT && !src->data) ||
	    (src->type == DATA_FLOAT && !src->fdata))
		return NULL;

	gchar *nameBuff = g_build_filename(com.pref.swap_dir, "siril_ndes-XXXXXX", NULL);
#ifndef _WIN32
	int fd = g_mkstemp(nameBuff);
#else
	int fd = g_mkstemp_full(nameBuff, _O_RDWR | _O_CREAT | _O_EXCL | _O_BINARY | _O_TEMPORARY,
	    _S_IREAD | _S_IWRITE);
#endif
	if (fd < 0) {
		siril_log_error(_("File I/O Error: Unable to create snapshot swap file in %s: [%s]\n"),
		                com.pref.swap_dir, strerror(errno));
		g_free(nameBuff);
		return NULL;
	}
	swap_mark_delete_on_close(fd, nameBuff);
	g_free(nameBuff);

	size_t n = (size_t)src->naxes[0] * src->naxes[1] * src->naxes[2];
	size_t bytes = n * (src->type == DATA_USHORT ? sizeof(WORD) : sizeof(float));
	errno = 0;
	ssize_t written = (src->type == DATA_USHORT) ?
			write(fd, src->data, bytes) : write(fd, src->fdata, bytes);
	if (written < (ssize_t)bytes) {
		siril_log_error(_("File I/O Error: Unable to write snapshot swap file: [%s]\n"),
		                strerror(errno));
		g_close(fd, NULL);
		return NULL;
	}

	nde_snap *s = g_new0(nde_snap, 1);
	s->refcount     = 1;
	s->fd           = fd;
	s->rx           = src->rx;
	s->ry           = src->ry;
	s->nchans       = (guint)src->naxes[2];
	s->type         = src->type;
	s->bytes        = (gint64)bytes;
	s->wcsdata      = src->keywords.wcsdata;
	if (src->keywords.wcslib) {
		int wstatus = -1;
		s->wcslib = wcs_deepcopy(src->keywords.wcslib, &wstatus);
		if (wstatus)
			siril_log_debug("nde snapstore: could not copy wcslib struct\n");
	}
	s->focal_length = src->keywords.focal_length;
	return s;
}

nde_snap *nde_snap_ref(nde_snap *s) {
	if (s)
		g_atomic_int_inc(&s->refcount);
	return s;
}

void nde_snap_unref(nde_snap *s) {
	if (!s)
		return;
	if (!g_atomic_int_dec_and_test(&s->refcount))
		return;
	g_mutex_lock(&store_mutex);
	deregister_locked(s);
	g_mutex_unlock(&store_mutex);
	if (s->fd >= 0)
		g_close(s->fd, NULL);   /* delete-on-close removes the file */
	if (s->wcslib) {
		wcsfree(s->wcslib);
		free(s->wcslib);
	}
	g_free(s);
}

fits *nde_snap_read(const nde_snap *s) {
	if (!s)
		return NULL;
	fits *out = NULL;
	if (new_fit_image(&out, (int)s->rx, (int)s->ry, (int)s->nchans, s->type))
		return NULL;

	void *dst = (s->type == DATA_USHORT) ? (void *)out->data : (void *)out->fdata;
	if (lseek(s->fd, 0, SEEK_SET) == (off_t)-1) {
		siril_log_error(_("nde snapstore: error seeking swap file: [%s]\n"), strerror(errno));
		clearfits(out);
		free(out);
		return NULL;
	}
	errno = 0;
	if (read(s->fd, dst, s->bytes) < (ssize_t)s->bytes) {
		siril_log_error(_("nde snapstore: read failed with error [%s]\n"), strerror(errno));
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
			siril_log_debug("nde snapstore: could not copy wcslib struct\n");
	} else {
		reset_wcsdata(out);
	}
	out->keywords.focal_length = s->focal_length;
	return out;
}

int nde_snap_read_into(const nde_snap *s, fits *dst) {
	if (!s || !dst)
		return 1;
	if (lseek(s->fd, 0, SEEK_SET) == (off_t)-1) {
		siril_log_error(_("nde snapstore: error seeking swap file: [%s]\n"), strerror(errno));
		return 1;
	}
	size_t n = (size_t)s->rx * s->ry;
	dst->rx = dst->naxes[0] = s->rx;
	dst->ry = dst->naxes[1] = s->ry;
	dst->naxes[2] = s->nchans;
	dst->naxis = s->nchans == 1 ? 2 : 3;
	dst->type = s->type;

	errno = 0;
	if (s->type == DATA_USHORT) {
		WORD *newdata = realloc(dst->data, s->bytes);
		if (!newdata) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		dst->data = newdata;
		if (read(s->fd, dst->data, s->bytes) < (ssize_t)s->bytes) {
			siril_log_error(_("nde snapstore: read failed with error [%s]\n"), strerror(errno));
			return 1;
		}
		dst->pdata[RLAYER] = dst->data;
		if (dst->naxes[2] > 1) {
			dst->pdata[GLAYER] = dst->data + n;
			dst->pdata[BLAYER] = dst->data + n * 2;
		} else {
			dst->pdata[GLAYER] = dst->pdata[BLAYER] = dst->pdata[RLAYER];
		}
	} else {
		float *newdata = realloc(dst->fdata, s->bytes);
		if (!newdata) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		dst->fdata = newdata;
		if (read(s->fd, dst->fdata, s->bytes) < (ssize_t)s->bytes) {
			siril_log_error(_("nde snapstore: read failed with error [%s]\n"), strerror(errno));
			return 1;
		}
		dst->fpdata[RLAYER] = dst->fdata;
		if (dst->naxes[2] > 1) {
			dst->fpdata[GLAYER] = dst->fdata + n;
			dst->fpdata[BLAYER] = dst->fdata + n * 2;
		} else {
			dst->fpdata[GLAYER] = dst->fpdata[BLAYER] = dst->fpdata[RLAYER];
		}
	}
	full_stats_invalidation_from_fit(dst);
	return 0;
}

gint64 nde_snap_size(const nde_snap *s) {
	return s ? s->bytes : 0;
}

/* ======================================================================= */
/* Tag registry                                                            */
/* ======================================================================= */

void nde_snap_set_tag(nde_snap *s, gint item_id, gint64 record_id, gboolean post) {
	if (!s)
		return;
	g_mutex_lock(&store_mutex);
	registry_ensure_locked();
	deregister_locked(s);   /* drop any previous tag of this snap */
	s->tag.item_id = item_id;
	s->tag.record_id = record_id;
	s->tag.post = post;
	s->tagged = TRUE;
	snap_tag *key = g_new(snap_tag, 1);
	*key = s->tag;
	g_hash_table_replace(registry, key, s);   /* latest same-tag snap wins */
	g_mutex_unlock(&store_mutex);
}

gboolean nde_snap_tag_get(const nde_snap *s, gint *item_id, gint64 *record_id, gboolean *post) {
	if (!s || !s->tagged)
		return FALSE;
	if (item_id)   *item_id = s->tag.item_id;
	if (record_id) *record_id = s->tag.record_id;
	if (post)      *post = s->tag.post;
	return TRUE;
}

gboolean nde_snapstore_has(gint item_id, gint64 record_id, gboolean post) {
	snap_tag t = { item_id, record_id, post };
	g_mutex_lock(&store_mutex);
	gboolean found = registry && g_hash_table_contains(registry, &t);
	g_mutex_unlock(&store_mutex);
	return found;
}

/* Move @s to the pool head if it is a pool member.  Mutex held. */
static void pool_touch_locked(nde_snap *s) {
	GList *link = g_queue_find(&pool, s);
	if (link) {
		g_queue_unlink(&pool, link);
		g_queue_push_head_link(&pool, link);
	}
}

fits *nde_snapstore_lookup(gint item_id, gint64 record_id, gboolean post) {
	snap_tag t = { item_id, record_id, post };
	g_mutex_lock(&store_mutex);
	stats.lookups++;
	nde_snap *s = registry ? g_hash_table_lookup(registry, &t) : NULL;
	if (s) {
		stats.hits++;
		nde_snap_ref(s);
		pool_touch_locked(s);
	}
	g_mutex_unlock(&store_mutex);
	if (!s)
		return NULL;
	fits *out = nde_snap_read(s);   /* outside the lock */
	nde_snap_unref(s);
	return out;
}

/* ======================================================================= */
/* LRU cache pool                                                          */
/* ======================================================================= */

static gint64 pool_budget_bytes(void) {
	return (gint64)com.pref.nde_cache_mb * 1024 * 1024;
}

/* Unlink pool entries chosen by @pick into @doomed (which takes the pool's
 * refs).  Mutex held; the actual unrefs happen OUTSIDE the lock because
 * final unref re-enters it. */
static void pool_collect_locked(GPtrArray *doomed,
                                gboolean (*pick)(nde_snap *, gpointer),
                                gpointer user) {
	GList *l = pool.head;
	while (l) {
		GList *next = l->next;
		nde_snap *s = l->data;
		if (pick(s, user)) {
			g_queue_unlink(&pool, l);
			g_list_free_1(l);
			pool_bytes -= s->bytes;
			stats.evictions++;
			g_ptr_array_add(doomed, s);
		}
		l = next;
	}
}

static void pool_release(GPtrArray *doomed) {
	for (guint i = 0; i < doomed->len; i++)
		nde_snap_unref(g_ptr_array_index(doomed, i));
	g_ptr_array_unref(doomed);
}

static gboolean pick_same_tag(nde_snap *s, gpointer user) {
	return tag_equal(&s->tag, user);
}

static gboolean pick_all(nde_snap *s, gpointer user) {
	(void)s; (void)user;
	return TRUE;
}

static gboolean pick_invalidated(nde_snap *s, gpointer user) {
	const snap_tag *edit = user;
	if (s->tag.item_id != edit->item_id)
		return FALSE;
	/* POST(K >= M) and PRE(K > M) no longer describe the amended chain. */
	return s->tag.post ? (s->tag.record_id >= edit->record_id)
	                   : (s->tag.record_id >  edit->record_id);
}

void nde_snapstore_deposit(const fits *state, gint item_id, gint64 record_id) {
	if (!state || record_id <= 0 || com.pref.nde_cache_mb <= 0)
		return;
	nde_snap *s = nde_snap_create(state);   /* I/O outside the lock */
	if (!s)
		return;
	nde_snap_set_tag(s, item_id, record_id, TRUE);

	snap_tag t = { item_id, record_id, TRUE };
	GPtrArray *doomed = g_ptr_array_new();
	g_mutex_lock(&store_mutex);
	/* replace any same-tag pool entry, then insert at MRU */
	pool_collect_locked(doomed, pick_same_tag, &t);
	g_queue_push_head(&pool, s);   /* pool takes the initial ref */
	pool_bytes += s->bytes;
	stats.deposits++;
	/* evict LRU over budget (never the snap just inserted) */
	gint64 budget = pool_budget_bytes();
	while (pool_bytes > budget && pool.length > 1) {
		nde_snap *victim = g_queue_pop_tail(&pool);
		pool_bytes -= victim->bytes;
		stats.evictions++;
		g_ptr_array_add(doomed, victim);
	}
	g_mutex_unlock(&store_mutex);
	pool_release(doomed);
}

void nde_snapstore_invalidate_from(gint item_id, gint64 record_id) {
	snap_tag edit = { item_id, record_id, FALSE /* unused by pick */ };
	GPtrArray *doomed = g_ptr_array_new();
	g_mutex_lock(&store_mutex);
	pool_collect_locked(doomed, pick_invalidated, &edit);
	g_mutex_unlock(&store_mutex);
	pool_release(doomed);
}

static gboolean pick_record(nde_snap *s, gpointer user) {
	const gint64 *id = user;
	return s->tag.record_id == *id;
}

void nde_snapstore_evict_record(gint64 record_id) {
	GPtrArray *doomed = g_ptr_array_new();
	g_mutex_lock(&store_mutex);
	pool_collect_locked(doomed, pick_record, &record_id);
	g_mutex_unlock(&store_mutex);
	pool_release(doomed);
}

void nde_snapstore_pool_purge(void) {
	GPtrArray *doomed = g_ptr_array_new();
	g_mutex_lock(&store_mutex);
	pool_collect_locked(doomed, pick_all, NULL);
	g_mutex_unlock(&store_mutex);
	pool_release(doomed);
}

/* ======================================================================= */
/* Stats                                                                   */
/* ======================================================================= */

void nde_snapstore_stats(nde_snapstore_stats_t *out) {
	g_mutex_lock(&store_mutex);
	*out = stats;
	g_mutex_unlock(&store_mutex);
}

void nde_snapstore_stats_reset(void) {
	g_mutex_lock(&store_mutex);
	memset(&stats, 0, sizeof(stats));
	g_mutex_unlock(&store_mutex);
}
