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

#include "core/siril.h"
#include "opencv/opencv.h"
#include "clahe.h"
#include "core/op_descriptors.h"
#include "core/nde_history.h"

/* NDE serializers for filters.clahe.  The hook reads clip and tileSize
 * (cvClahe(fit, clip, tileSize)). */
static gchar *clahe_serialize(gconstpointer user) {
	const clahe_params *d = user;
	GString *kv = nde_kv_start();
	nde_kv_add_double(kv, "clip", d->clip);
	nde_kv_add_int(kv, "tileSize", d->tileSize);
	return nde_kv_end(kv);
}

static gpointer clahe_deserialize(const gchar *blob, int version) {
	if (version > op_desc_clahe.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	double clip;
	gint64 tileSize;
	NDE_KV_REQUIRE(kv, nde_kv_get_double(kv, "clip", &clip) &&
	                   nde_kv_get_int(kv, "tileSize", &tileSize));
	clahe_params *d = calloc(1, sizeof(*d));
	if (d) {
		d->destroy_fn = free;
		d->clip = clip;
		d->tileSize = (int)tileSize;
	}
	g_hash_table_unref(kv);
	return d;
}

/* Op descriptor — single source of truth for this operation (op_descriptor.h) */
const op_descriptor op_desc_clahe = {
	.id = "filters.clahe", .version = 1,
	.image_hook = clahe_image_hook,
	.log_hook = clahe_log_hook,
	.description = N_("CLAHE"),
	.mem_ratio = 2.0f,
	.flags = OP_MASK_CAPABLE,
	.serialize = clahe_serialize, .deserialize = clahe_deserialize,
};

/* The actual CLAHE processing hook */
int clahe_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	clahe_params *params = (clahe_params *)args->user;
	if (!params)
		return 1;

	return cvClahe(fit, params->clip, params->tileSize);
}

gchar *clahe_log_hook(gpointer p, log_hook_detail detail) {
	clahe_params *params = (clahe_params *) p;
	if (!params) return NULL;
	gchar *message = g_strdup_printf(_("CLAHE (size=%d, clip=%.2f)"), params->tileSize, params->clip);
	return message;
}

/* GUI callbacks moved to src/gui/clahe.c */
