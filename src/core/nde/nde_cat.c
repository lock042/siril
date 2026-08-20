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

/* See nde_cat.h for the model (stash → adopt → persist → replay). */

#include "core/siril.h"
#include "core/nde/nde_cat.h"
#include "core/siril_log.h"

static GMutex cat_mutex;        /* leaf: nothing acquired while held */
static GHashTable *cat_table;   /* gint64* -> siril_catalogue* (owned) */
static siril_catalogue *cat_pending;

static void cat_value_free(gpointer p) {
	siril_catalog_free(p);
}

siril_catalogue *nde_cat_copy(const siril_catalogue *cat) {
	if (!cat)
		return NULL;
	siril_catalogue *out = siril_catalog_new(cat->cat_index);
	out->center_ra = cat->center_ra;
	out->center_dec = cat->center_dec;
	out->radius = cat->radius;
	out->limitmag = cat->limitmag;
	out->phot = cat->phot;
	out->columns = cat->columns;
	out->nbitems = cat->nbitems;
	if (cat->nbitems > 0 && cat->cat_items) {
		out->cat_items = calloc(cat->nbitems, sizeof(cat_item));
		for (int i = 0; i < cat->nbitems; i++) {
			const cat_item *src = &cat->cat_items[i];
			cat_item *dst = &out->cat_items[i];
			dst->ra = src->ra;
			dst->dec = src->dec;
			dst->pmra = src->pmra;
			dst->pmdec = src->pmdec;
			dst->mag = src->mag;
			dst->bmag = src->bmag;
			dst->e_mag = src->e_mag;
			dst->e_bmag = src->e_bmag;
			dst->teff = src->teff;
			dst->gaiasourceid = src->gaiasourceid;
			dst->index = src->index;
			dst->name = g_strdup(src->name);
			dst->included = TRUE;   /* projection re-derives exclusions */
			if (src->xp_sampled) {
				dst->xp_sampled = malloc(XPSAMPLED_LEN * sizeof(double));
				memcpy(dst->xp_sampled, src->xp_sampled,
				       XPSAMPLED_LEN * sizeof(double));
			}
		}
	}
	return out;
}

void nde_cat_stash_pending(const siril_catalogue *cat) {
	siril_catalogue *copy = nde_cat_copy(cat);
	g_mutex_lock(&cat_mutex);
	if (cat_pending)
		siril_catalog_free(cat_pending);
	cat_pending = copy;
	g_mutex_unlock(&cat_mutex);
}

/* The photometric ops whose captures may claim a pending catalogue: the
 * single-image record and the group joint record (nde_joint.h). */
static gboolean op_uses_catalogue(const char *op_id) {
	return op_id && (!g_strcmp0(op_id, "color.photometric_cc") ||
	                 !g_strcmp0(op_id, "flis.group_calibration"));
}

static void register_locked(gint64 record_id, siril_catalogue *cat) {
	if (!cat_table)
		cat_table = g_hash_table_new_full(g_int64_hash, g_int64_equal,
		                                  g_free, cat_value_free);
	gint64 *key = g_new(gint64, 1);
	*key = record_id;
	g_hash_table_replace(cat_table, key, cat);
}

void nde_cat_adopt_pending(gint64 record_id, const char *op_id) {
	g_mutex_lock(&cat_mutex);
	if (cat_pending) {
		if (record_id > 0 && op_uses_catalogue(op_id)) {
			register_locked(record_id, cat_pending);
		} else {
			/* Any other capture means the pipeline's run was not the one
			 * being recorded — the stash is stale. */
			siril_catalog_free(cat_pending);
		}
		cat_pending = NULL;
	}
	g_mutex_unlock(&cat_mutex);
}

void nde_cat_register(gint64 record_id, siril_catalogue *cat) {
	if (!cat)
		return;
	g_mutex_lock(&cat_mutex);
	register_locked(record_id, cat);
	g_mutex_unlock(&cat_mutex);
}

siril_catalogue *nde_cat_get_copy(gint64 record_id) {
	g_mutex_lock(&cat_mutex);
	siril_catalogue *found = cat_table ?
			g_hash_table_lookup(cat_table, &record_id) : NULL;
	siril_catalogue *copy = nde_cat_copy(found);
	g_mutex_unlock(&cat_mutex);
	return copy;
}

siril_catalogue *nde_cat_peek(gint64 record_id) {
	g_mutex_lock(&cat_mutex);
	siril_catalogue *found = cat_table ?
			g_hash_table_lookup(cat_table, &record_id) : NULL;
	g_mutex_unlock(&cat_mutex);
	return found;
}

gboolean nde_cat_has(gint64 record_id) {
	return nde_cat_peek(record_id) != NULL;
}

void nde_cat_drop(gint64 record_id) {
	g_mutex_lock(&cat_mutex);
	if (cat_table)
		g_hash_table_remove(cat_table, &record_id);
	g_mutex_unlock(&cat_mutex);
}

void nde_cat_purge(void) {
	g_mutex_lock(&cat_mutex);
	if (cat_table)
		g_hash_table_remove_all(cat_table);
	if (cat_pending) {
		siril_catalog_free(cat_pending);
		cat_pending = NULL;
	}
	g_mutex_unlock(&cat_mutex);
}
