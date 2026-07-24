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

#include <string.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "algos/fitting.h"
#include "algos/statistics.h"
#include "io/image_format_fits.h"
#include "linear_match.h"
#include "core/op_descriptors.h"
#include "core/nde_history.h"

/* NDE serializers (phase 4.5 Convention 1 — file operands).  low/high are POD;
 * the reference image is pinned by path (+hash/size at capture) and loaded by
 * replay_pre into the embedded `ref` via the new_linear_match_data ownership
 * dance (which zeroes the caller's fits so clearfits frees each buffer once). */
static gchar *linear_match_serialize(gconstpointer user) {
	const struct linear_match_data *d = user;
	GString *kv = nde_kv_start();
	nde_kv_add_double(kv, "low", d->low);
	nde_kv_add_double(kv, "high", d->high);
	nde_kv_add_bool(kv, "force_to_float", d->force_to_float);
	nde_kv_add_str(kv, "operand_path", d->operand_path ? d->operand_path : "");
	gint64 size = 0;
	gchar *sha = nde_file_sha256(d->operand_path, &size);
	if (sha) {
		nde_kv_add_str(kv, "operand_sha256", sha);
		nde_kv_add_int(kv, "operand_size", size);
		g_free(sha);
	}
	return nde_kv_end(kv);
}

static gpointer linear_match_deserialize(const gchar *blob, int version) {
	if (version > op_desc_linear_match.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	double low, high;
	gboolean force_to_float;
	const char *path = nde_kv_get_str(kv, "operand_path");
	if (!nde_kv_get_double(kv, "low", &low) ||
	    !nde_kv_get_double(kv, "high", &high) ||
	    !nde_kv_get_bool(kv, "force_to_float", &force_to_float) ||
	    !path || !*path) {
		g_hash_table_unref(kv);
		return NULL;
	}
	struct linear_match_data *d = calloc(1, sizeof(*d));
	if (d) {
		d->destroy_fn = free_linear_match_data;
		d->low = low;
		d->high = high;
		d->force_to_float = force_to_float;
		d->operand_path = strdup(path);
	}
	g_hash_table_unref(kv);
	return d;
}

/* replay_pre: verify + readfits the reference into the embedded `ref`.  The
 * struct owns `ref` (freed by free_linear_match_data's clearfits(&data->ref)),
 * so we readfits directly into it. */
static int linear_match_replay_pre(gpointer user, GHashTable *kv, fits *target) {
	(void)target;
	struct linear_match_data *d = user;
	const char *path = nde_kv_get_str(kv, "operand_path");
	gint64 expect_size = 0;
	nde_kv_get_int(kv, "operand_size", &expect_size);
	const char *sha = nde_kv_get_str(kv, "operand_sha256");
	if (!nde_operand_verify(path, expect_size, sha))
		return 1;
	clearfits(&d->ref);
	if (readfits(path, &d->ref, NULL, d->force_to_float))
		return 1;
	return 0;
}

/* Op descriptor — single source of truth for this operation (op_descriptor.h) */
const op_descriptor op_desc_linear_match = {
	.id = "color.linear_match", .version = 1,
	.image_hook = linear_match_image_hook,
	.log_hook = linear_match_log_hook,
	.description = N_("Linear Match"),
	.mem_ratio = 0.0f,
	.flags = 0,
	.serialize = linear_match_serialize, .deserialize = linear_match_deserialize,
	.replay_pre = linear_match_replay_pre,
};

static void apply_linear_to_fits_ushort(fits *fit, const double *a, const double *b) {
	size_t size = fit->rx * fit->ry;

	invalidate_stats_from_fit(fit);
	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		for (size_t i = 0; i < size; i++) {
			fit->pdata[channel][i] = round_to_WORD(fit->pdata[channel][i] * a[channel] + b[channel] * USHRT_MAX_DOUBLE);
		}
	}
}

static void apply_linear_to_fits_float(fits *fit, const double *a, const double *b) {
	size_t size = fit->rx * fit->ry;

	invalidate_stats_from_fit(fit);
	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		for (size_t i = 0; i < size; i++) {
			fit->fpdata[channel][i] = fit->fpdata[channel][i] * a[channel] + b[channel];
		}
	}
}

void apply_linear_to_fits(fits *fit, double *a, double *b) {
	if (fit->type == DATA_USHORT) {
		apply_linear_to_fits_ushort(fit, a, b);
	} else if (fit->type == DATA_FLOAT) {
		apply_linear_to_fits_float(fit, a, b);
	}
}

struct linear_match_data *new_linear_match_data(fits *ref_fit, double low, double high,
                                                const char *ref_path) {
	struct linear_match_data *data = calloc(1, sizeof(struct linear_match_data));
	if (!data)
		return NULL;
	data->destroy_fn = free_linear_match_data;
	data->low = low;
	data->high = high;
	/* pin the reference for NDE replay (phase 4.5 Convention 1) */
	data->operand_path = ref_path ? strdup(ref_path) : NULL;
	data->force_to_float = (ref_fit->type == DATA_FLOAT);
	/* Caller must have already loaded ref_fit; we take ownership of its data */
	data->ref = *ref_fit;
	/* Zero out the caller's copy so clearfits() won't free the same buffers */
	memset(ref_fit, 0, sizeof(fits));
	return data;
}

void free_linear_match_data(void *p) {
	struct linear_match_data *data = (struct linear_match_data *)p;
	if (!data)
		return;
	clearfits(&data->ref);
	free(data->operand_path);
	free(data);
}

int linear_match_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct linear_match_data *data = (struct linear_match_data *)args->user;
	double a[3] = { 0.0 }, b[3] = { 0.0 };
	if (find_linear_coeff(fit, &data->ref, data->low, data->high, a, b, NULL))
		return 1;
	apply_linear_to_fits(fit, a, b);
	return 0;
}

gchar *linear_match_log_hook(gpointer p, log_hook_detail detail) {
	struct linear_match_data *data = (struct linear_match_data *)p;
	return g_strdup_printf(_("Linear match (low: %.3f, high: %.3f)"), data->low, data->high);
}
