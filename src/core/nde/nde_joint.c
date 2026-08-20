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

/* See nde_joint.h for the model (one record, N participating chains). */

#include <math.h>

#include "core/siril.h"
#include "core/nde/nde_joint.h"
#include "core/nde/nde_history.h"
#include "core/nde/nde_op_class.h"
#include "core/nde/nde_checkpoint.h"     /* baseline offsets (register L2 anchor) */
#include "core/nde/nde_compositing.h"
#include "core/nde/nde_replay.h"
#include "core/op_descriptors.h"
#include "core/siril_log.h"
#include "core/proto.h"          /* get_normalized_value */
#include "algos/statistics.h"
#include "algos/colors.h"        /* get_coeff_for_wb */
#include "algos/photometric_cc.h"    /* photometric_cc_image_hook (group L2) */
#include "io/image_format_fits.h"    /* copy_fits_metadata (WCS donor) */
#include "io/image_format_flis.h"
#include "io/flis_compose.h"     /* flis_render_layers_ctx */
#include "opencv/opencv.h"       /* cvTransformImage (register L1) */
#include "registration/flis_register.h"  /* the register re-solve (L2) */

#define OP_LAYERS_MATCH      "flis.layers_match"
#define OP_GROUP_CALIBRATION "flis.group_calibration"
#define OP_REGISTER          "flis.register"

/* A corrupt blob must not size an allocation: no real composition has more
 * layers than this. */
#define JOINT_MAX_PARTICIPANTS 256

/* Nor size an IMAGE: a registration record's output dimensions come straight
 * out of the params, and a corrupt one must not ask cvTransformImage for a
 * terabyte.  Well above any real astronomical frame. */
#define NDE_JOINT_MAX_DIM 200000

gboolean nde_joint_is_op(const char *op_id) {
	return nde_op_class_for(op_id)->family == NDE_OPC_JOINT;
}

/* Registration warps its participants AND moves them on the canvas; the other
 * two joint ops only rescale pixels.  The trait comes from the descriptor's
 * OP_GEOMETRY_CHANGING, so this stays in step with the flag that
 * nde_joint_geometry_signature already reads. */
gboolean nde_joint_is_geometric_op(const char *op_id) {
	const nde_op_class *cls = nde_op_class_for(op_id);
	return cls->family == NDE_OPC_JOINT && (cls->traits & NDE_OPT_GEOMETRIC);
}

/* ======================================================================= */
/* Participant block codec (shared by both joint ops)                      */
/* ======================================================================= */

struct nde_joint_layers_match_data *nde_joint_layers_match_data_new(guint n) {
	struct nde_joint_layers_match_data *p = g_new0(struct nde_joint_layers_match_data, 1);
	p->destroy_fn = nde_joint_layers_match_data_free;
	p->n = n;
	p->parts = g_new0(nde_joint_participant, n);
	return p;
}

void nde_joint_layers_match_data_free(void *ptr) {
	struct nde_joint_layers_match_data *p = ptr;
	if (!p)
		return;
	for (guint i = 0; i < p->n; i++)
		g_free(p->parts[i].name);
	g_free(p->parts);
	g_free(p);
}

/* Append participant @k of @parts to @kv under the i<k>_* keys. */
static void participant_serialize(GString *kv, guint k,
                                  const nde_joint_participant *part) {
	gchar key[32];
	g_snprintf(key, sizeof(key), "i%u_item", k);
	nde_kv_add_int(kv, key, part->item_id);
	g_snprintf(key, sizeof(key), "i%u_name", k);
	nde_kv_add_str(kv, key, part->name ? part->name : "");
	g_snprintf(key, sizeof(key), "i%u_tint", k);
	nde_kv_add_bool(kv, key, part->tinted);
	if (part->tinted) {
		g_snprintf(key, sizeof(key), "i%u_tr", k);
		nde_kv_add_double(kv, key, part->tint[0]);
		g_snprintf(key, sizeof(key), "i%u_tg", k);
		nde_kv_add_double(kv, key, part->tint[1]);
		g_snprintf(key, sizeof(key), "i%u_tb", k);
		nde_kv_add_double(kv, key, part->tint[2]);
	}
	g_snprintf(key, sizeof(key), "i%u_s", k);
	nde_kv_add_double(kv, key, part->diag_scale);
	if (part->diag_offset != 0.0) {
		g_snprintf(key, sizeof(key), "i%u_o", k);
		nde_kv_add_double(kv, key, part->diag_offset);
	}
}

/* Parse participant @k from @kv into @part.  FALSE when the entry is absent
 * or invalid (the whole blob is then rejected — a joint record with a
 * half-readable participant list must not replay against the wrong layers). */
static gboolean participant_parse(GHashTable *kv, guint k,
                                  nde_joint_participant *part) {
	gchar key[32];
	gint64 item;
	g_snprintf(key, sizeof(key), "i%u_item", k);
	if (!nde_kv_get_int(kv, key, &item) || item <= 0)
		return FALSE;
	part->item_id = (gint)item;
	g_snprintf(key, sizeof(key), "i%u_name", k);
	part->name = g_strdup(nde_kv_get_str(kv, key) ? nde_kv_get_str(kv, key) : "");
	g_snprintf(key, sizeof(key), "i%u_tint", k);
	if (!nde_kv_get_bool(kv, key, &part->tinted))
		part->tinted = FALSE;
	part->tint[0] = part->tint[1] = part->tint[2] = 1.0;
	if (part->tinted) {
		double r, g, b;
		gchar kr[32], kg[32], kb[32];
		g_snprintf(kr, sizeof(kr), "i%u_tr", k);
		g_snprintf(kg, sizeof(kg), "i%u_tg", k);
		g_snprintf(kb, sizeof(kb), "i%u_tb", k);
		if (!nde_kv_get_double(kv, kr, &r) ||
		    !nde_kv_get_double(kv, kg, &g) ||
		    !nde_kv_get_double(kv, kb, &b) ||
		    !(r >= 0. && r <= 1. && g >= 0. && g <= 1. && b >= 0. && b <= 1.))
			return FALSE;
		part->tint[0] = r; part->tint[1] = g; part->tint[2] = b;
	}
	part->diag_scale = 1.0;
	g_snprintf(key, sizeof(key), "i%u_s", k);
	nde_kv_get_double(kv, key, &part->diag_scale);
	if (part->diag_scale < 0.0)
		return FALSE;
	part->diag_offset = 0.0;
	g_snprintf(key, sizeof(key), "i%u_o", k);
	nde_kv_get_double(kv, key, &part->diag_offset);
	return TRUE;
}

/* The participant count of a parsed blob: n present and in range.  Layers
 * match needs at least 2 (matching fewer is not a match); group calibration
 * accepts a single-layer group. */
static gboolean participant_count(GHashTable *kv, guint min, guint *n_out) {
	gint64 n;
	if (!nde_kv_get_int(kv, "n", &n) || n < (gint64)min || n > JOINT_MAX_PARTICIPANTS)
		return FALSE;
	*n_out = (guint)n;
	return TRUE;
}

gint *nde_joint_params_participants(const char *params, guint *n_out) {
	if (n_out)
		*n_out = 0;
	if (!params)
		return NULL;
	GHashTable *kv = nde_kv_parse(params);
	guint n;
	gint *items = NULL;
	if (participant_count(kv, 1, &n)) {
		items = g_new0(gint, n);
		for (guint k = 0; k < n; k++) {
			gchar key[32];
			gint64 item;
			g_snprintf(key, sizeof(key), "i%u_item", k);
			if (!nde_kv_get_int(kv, key, &item) || item <= 0) {
				g_free(items);
				items = NULL;
				n = 0;
				break;
			}
			items[k] = (gint)item;
		}
		if (items && n_out)
			*n_out = n;
	}
	g_hash_table_unref(kv);
	return items;
}

gint *nde_joint_record_participants(const struct nde_record *rec, guint *n_out) {
	if (n_out)
		*n_out = 0;
	if (!rec || !nde_joint_is_op(rec->op_id))
		return NULL;
	return nde_joint_params_participants(rec->params, n_out);
}

GPtrArray *nde_joint_output_pins(const char *op_id, const char *params,
                                 gint target_item_id) {
	if (!nde_joint_is_op(op_id))
		return NULL;
	guint n = 0;
	gint *items = nde_joint_params_participants(params, &n);
	if (!items)
		return NULL;
	/* A bare record to collect into, so the "output 0 is the target" invariant
	 * is maintained by the one function that knows it (nde_history.h) rather
	 * than restated here. */
	nde_record tmp = { .target_item_id = target_item_id };
	for (guint k = 0; k < n; k++) {
		gchar role[16];
		g_snprintf(role, sizeof(role), "in%u", k);
		nde_record_add_output(&tmp, role, items[k]);
	}
	g_free(items);
	return tmp.outputs;
}

void nde_joint_sync_outputs(struct nde_record *rec) {
	if (!rec)
		return;
	GPtrArray *outs = nde_joint_output_pins(rec->op_id, rec->params,
	                                        rec->target_item_id);
	/* Unreadable params, or not a joint record at all: keep whatever output
	 * list is already there.  One loaded from a file wrote its outputs down,
	 * and a blob this build cannot parse is no reason to forget them. */
	if (!outs)
		return;
	if (rec->outputs)
		g_ptr_array_unref(rec->outputs);
	rec->outputs = outs;
}

gboolean nde_joint_params_same_participants(const char *a, const char *b) {
	if (!a || !b)
		return FALSE;
	GHashTable *ka = nde_kv_parse(a), *kb = nde_kv_parse(b);
	guint na = 0, nb = 0;
	gboolean same = participant_count(ka, 1, &na) &&
	                participant_count(kb, 1, &nb) && na == nb;
	for (guint k = 0; same && k < na; k++) {
		gchar key[32];
		gint64 ia, ib;
		g_snprintf(key, sizeof(key), "i%u_item", k);
		same = nde_kv_get_int(ka, key, &ia) && nde_kv_get_int(kb, key, &ib) &&
		       ia == ib;
	}
	g_hash_table_unref(ka);
	g_hash_table_unref(kb);
	return same;
}

/* ======================================================================= */
/* Geometry signature (the L1/L2 discriminator for flis.register)          */
/* ======================================================================= */

gchar *nde_joint_geometry_signature(gint item_id, gint64 before_record_id) {
	GString *s = g_string_new(NULL);
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		const nde_record *rec = g_ptr_array_index(snap, i);
		/* Positional, not by id: after a committed reorder the ids are no
		 * longer monotonic in log order, and it is the ORDER the layer's
		 * pixels see that decides what geometry it reaches the record with. */
		if (before_record_id && rec->record_id == before_record_id)
			break;
		if (!nde_record_writes_item(rec, item_id))
			continue;
		const op_descriptor *op = op_descriptor_by_id(rec->op_id);
		if (!op || !(op->flags & OP_GEOMETRY_CHANGING))
			continue;
		/* op id AND params: amending a rotation's angle leaves the id alone
		 * but changes the geometry the layer arrives with, which is exactly
		 * the case the signature exists to catch. */
		g_string_append(s, rec->op_id ? rec->op_id : "?");
		g_string_append_c(s, '\x1f');
		g_string_append(s, rec->params ? rec->params : "");
		g_string_append_c(s, '\x1e');
	}
	if (snap)
		g_ptr_array_unref(snap);
	if (!s->len)
		return g_string_free(s, TRUE), g_strdup("");
	gchar *sig = g_compute_checksum_for_string(G_CHECKSUM_SHA256, s->str, -1);
	g_string_free(s, TRUE);
	return sig;
}

/* ======================================================================= */
/* Generation-keyed memo: record_id -> one opaque blob                     */
/* ======================================================================= */

/* Each joint op caches its analysis against the history generation it was
 * computed under, because one edit must cost ONE analysis however many
 * participants replay it: without this, an edit upstream of a five-layer
 * registration would re-solve the registration five times.  A generation bump
 * invalidates everything — there is no finer-grained invalidation, and the
 * analyses are cheap enough relative to a replay that there need not be.
 *
 * The value is an opaque blob; callers pack and unpack it.  Storing bytes
 * rather than a typed value is what lets one memo serve both the scalar ops'
 * pair of double arrays and the register op's array of solutions. */
typedef struct {
	guint64  generation;
	gsize    size;
	gpointer blob;
} nde_memo_entry;

typedef struct {
	GMutex      mutex;   /* leaf: nothing acquired while held */
	GHashTable *table;   /* gint64* -> nde_memo_entry*, created on first put */
} nde_memo;

static void nde_memo_entry_free(gpointer ptr) {
	nde_memo_entry *e = ptr;
	if (!e)
		return;
	g_free(e->blob);
	g_free(e);
}

/* Copies the blob out on a generation hit, NULL on a miss.  Caller g_free()s.
 * *size_out is the stored size, which a caller that knows what it expects
 * should check: a same-generation entry of the wrong length is a bug it can
 * catch for free. */
static gpointer nde_memo_get(nde_memo *m, gint64 record_id, guint64 generation,
                             gsize *size_out) {
	gpointer out = NULL;
	g_mutex_lock(&m->mutex);
	if (m->table) {
		nde_memo_entry *e = g_hash_table_lookup(m->table, &record_id);
		if (e && e->generation == generation) {
			out = g_memdup2(e->blob, e->size);
			if (size_out)
				*size_out = e->size;
		}
	}
	g_mutex_unlock(&m->mutex);
	return out;
}

static void nde_memo_put(nde_memo *m, gint64 record_id, guint64 generation,
                         gconstpointer blob, gsize size) {
	nde_memo_entry *e = g_new0(nde_memo_entry, 1);
	e->generation = generation;
	e->size = size;
	e->blob = g_memdup2(blob, size);
	gint64 *key = g_new(gint64, 1);
	*key = record_id;
	g_mutex_lock(&m->mutex);
	if (!m->table)
		m->table = g_hash_table_new_full(g_int64_hash, g_int64_equal,
		                                 g_free, nde_memo_entry_free);
	g_hash_table_replace(m->table, key, e);
	g_mutex_unlock(&m->mutex);
}

static void nde_memo_invalidate(nde_memo *m) {
	g_mutex_lock(&m->mutex);
	if (m->table)
		g_hash_table_remove_all(m->table);
	g_mutex_unlock(&m->mutex);
}

/* ======================================================================= */
/* Factor memo: record_id -> the per-participant factors {a[n], b[n]}      */
/* ======================================================================= */

static nde_memo joint_factors_memo;
static nde_memo register_solution_memo;

void nde_joint_cache_invalidate(void) {
	nde_memo_invalidate(&joint_factors_memo);
	nde_memo_invalidate(&register_solution_memo);
}

/* The two arrays go in as one 2n-double blob, a then b, so n comes back out of
 * the stored size — which matters because the lookup happens before the
 * record's params are parsed, i.e. before n is known. */
static gboolean joint_factors_get(gint64 record_id, guint64 generation,
                                  double **a_out, double **b_out) {
	gsize size = 0;
	double *blob = nde_memo_get(&joint_factors_memo, record_id, generation, &size);
	if (!blob)
		return FALSE;
	gsize n = size / (2 * sizeof(double));
	*a_out = g_memdup2(blob, n * sizeof(double));
	*b_out = g_memdup2(blob + n, n * sizeof(double));
	g_free(blob);
	return TRUE;
}

static void joint_factors_put(gint64 record_id, guint64 generation, guint n,
                              const double *a, const double *b) {
	double *blob = g_new(double, 2 * (gsize)n);
	memcpy(blob, a, n * sizeof(double));
	memcpy(blob + n, b, n * sizeof(double));
	nde_memo_put(&joint_factors_memo, record_id, generation,
	             blob, 2 * (gsize)n * sizeof(double));
	g_free(blob);
}

/* ======================================================================= */
/* flis.layers_match descriptor                                            */
/* ======================================================================= */

static gchar *layers_match_serialize(gconstpointer user) {
	const struct nde_joint_layers_match_data *p = user;
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "n", p->n);
	for (guint k = 0; k < p->n; k++)
		participant_serialize(kv, k, &p->parts[k]);
	return nde_kv_end(kv);
}

static gpointer layers_match_deserialize(const gchar *blob, int version) {
	if (version > op_desc_flis_layers_match.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	guint n;
	struct nde_joint_layers_match_data *p = NULL;
	if (participant_count(kv, 2, &n)) {
		p = nde_joint_layers_match_data_new(n);
		for (guint k = 0; p && k < n; k++) {
			if (!participant_parse(kv, k, &p->parts[k])) {
				nde_joint_layers_match_data_free(p);
				p = NULL;
			}
		}
		/* A layer listed twice would receive its factor twice. */
		for (guint i = 0; p && i < n; i++)
			for (guint j = i + 1; p && j < n; j++)
				if (p->parts[i].item_id == p->parts[j].item_id) {
					nde_joint_layers_match_data_free(p);
					p = NULL;
				}
	}
	g_hash_table_unref(kv);
	return p;
}

static gchar *layers_match_log_hook(gpointer user, log_hook_detail detail) {
	const struct nde_joint_layers_match_data *p = user;
	(void)detail;
	if (!p)
		return g_strdup(_("Layers match"));
	GString *s = g_string_new(_("Layers match"));
	const char *sep = " — ";
	for (guint i = 0; i < p->n; i++) {
		if (p->parts[i].name && *p->parts[i].name) {
			g_string_append(s, sep);
			g_string_append(s, p->parts[i].name);
			sep = ", ";
		}
	}
	if (!g_strcmp0(sep, " — "))
		g_string_append_printf(s, _(" (%u layers)"), p->n);
	return g_string_free(s, FALSE);
}

/* Never reached: the replay driver special-cases joint ops before the
 * generic dispatch (replay_apply_records), and capture applies the factors
 * itself.  Reaching this means a call site handed the op to a generic
 * worker it must not go through. */
static int layers_match_image_hook(struct generic_img_args *args, fits *fit,
                                   int nb_threads) {
	(void)args; (void)fit; (void)nb_threads;
	siril_log_message(_("flis.layers_match cannot run as a generic image operation\n"));
	return 1;
}

const op_descriptor op_desc_flis_layers_match = {
	.id = OP_LAYERS_MATCH, .version = 1,
	.image_hook = layers_match_image_hook,
	.log_hook = layers_match_log_hook,
	.description = N_("Layers match"),
	.mem_ratio = 1.0f,
	.flags = 0,
	.serialize = layers_match_serialize,
	.deserialize = layers_match_deserialize,
};

/* ======================================================================= */
/* flis.group_calibration descriptor                                       */
/* ======================================================================= */

struct nde_joint_group_calib_data *nde_joint_group_calib_data_new(guint n) {
	struct nde_joint_group_calib_data *p = g_new0(struct nde_joint_group_calib_data, 1);
	p->destroy_fn = nde_joint_group_calib_data_free;
	p->n = n;
	p->parts = g_new0(nde_joint_participant, n);
	p->manual_kw[0] = p->manual_kw[1] = p->manual_kw[2] = 1.0;
	p->diag_K[0] = p->diag_K[1] = p->diag_K[2] = 1.0;
	return p;
}

void nde_joint_group_calib_data_free(void *ptr) {
	struct nde_joint_group_calib_data *p = ptr;
	if (!p)
		return;
	for (guint i = 0; i < p->n; i++)
		g_free(p->parts[i].name);
	g_free(p->parts);
	g_free(p->pcc_blob);
	g_free(p);
}

static gchar *group_calib_serialize(gconstpointer user) {
	const struct nde_joint_group_calib_data *p = user;
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "n", p->n);
	for (guint k = 0; k < p->n; k++)
		participant_serialize(kv, k, &p->parts[k]);
	nde_kv_add_int(kv, "kind", p->kind);
	nde_kv_add_int(kv, "group", p->group_id);
	if (p->kind == 0 /* CC */) {
		nde_kv_add_bool(kv, "manual", p->manual);
		if (p->manual) {
			nde_kv_add_double(kv, "mkw0", p->manual_kw[0]);
			nde_kv_add_double(kv, "mkw1", p->manual_kw[1]);
			nde_kv_add_double(kv, "mkw2", p->manual_kw[2]);
		} else {
			nde_kv_add_int(kv, "wx", p->white_sel.x);
			nde_kv_add_int(kv, "wy", p->white_sel.y);
			nde_kv_add_int(kv, "ww", p->white_sel.w);
			nde_kv_add_int(kv, "wh", p->white_sel.h);
			nde_kv_add_int(kv, "bx", p->black_sel.x);
			nde_kv_add_int(kv, "by", p->black_sel.y);
			nde_kv_add_int(kv, "bw", p->black_sel.w);
			nde_kv_add_int(kv, "bh", p->black_sel.h);
			nde_kv_add_double(kv, "low", p->low);
			nde_kv_add_double(kv, "high", p->high);
		}
	}
	if (p->pcc_blob)
		nde_kv_add_str(kv, "pcc", p->pcc_blob);   /* nested, codec-escaped */
	for (int c = 0; c < 3; c++) {
		gchar key[8];
		g_snprintf(key, sizeof(key), "dK%d", c);
		nde_kv_add_double(kv, key, p->diag_K[c]);
		g_snprintf(key, sizeof(key), "dO%d", c);
		nde_kv_add_double(kv, key, p->diag_O[c]);
	}
	return nde_kv_end(kv);
}

static gpointer group_calib_deserialize(const gchar *blob, int version) {
	if (version > op_desc_flis_group_calibration.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	guint n;
	struct nde_joint_group_calib_data *p = NULL;
	gint64 kind;
	if (participant_count(kv, 1, &n) &&
	    nde_kv_get_int(kv, "kind", &kind) &&
	    kind >= 0 && kind <= NDE_JOINT_GROUP_DIRECT) {
		p = nde_joint_group_calib_data_new(n);
		p->kind = (gint)kind;
		for (guint k = 0; p && k < n; k++) {
			if (!participant_parse(kv, k, &p->parts[k])) {
				nde_joint_group_calib_data_free(p);
				p = NULL;
			}
		}
		for (guint i = 0; p && i < n; i++)
			for (guint j = i + 1; p && j < n; j++)
				if (p->parts[i].item_id == p->parts[j].item_id) {
					nde_joint_group_calib_data_free(p);
					p = NULL;
				}
	}
	if (p) {
		gint64 gi;
		if (nde_kv_get_int(kv, "group", &gi))
			p->group_id = (gint)gi;
		nde_kv_get_bool(kv, "manual", &p->manual);
		nde_kv_get_double(kv, "mkw0", &p->manual_kw[0]);
		nde_kv_get_double(kv, "mkw1", &p->manual_kw[1]);
		nde_kv_get_double(kv, "mkw2", &p->manual_kw[2]);
		gint64 v;
		if (nde_kv_get_int(kv, "wx", &v)) p->white_sel.x = (int)v;
		if (nde_kv_get_int(kv, "wy", &v)) p->white_sel.y = (int)v;
		if (nde_kv_get_int(kv, "ww", &v)) p->white_sel.w = (int)v;
		if (nde_kv_get_int(kv, "wh", &v)) p->white_sel.h = (int)v;
		if (nde_kv_get_int(kv, "bx", &v)) p->black_sel.x = (int)v;
		if (nde_kv_get_int(kv, "by", &v)) p->black_sel.y = (int)v;
		if (nde_kv_get_int(kv, "bw", &v)) p->black_sel.w = (int)v;
		if (nde_kv_get_int(kv, "bh", &v)) p->black_sel.h = (int)v;
		nde_kv_get_double(kv, "low", &p->low);
		nde_kv_get_double(kv, "high", &p->high);
		const char *pcc = nde_kv_get_str(kv, "pcc");
		p->pcc_blob = pcc && *pcc ? g_strdup(pcc) : NULL;
		gboolean diag_ok = TRUE;
		for (int c = 0; c < 3; c++) {
			gchar key[8];
			g_snprintf(key, sizeof(key), "dK%d", c);
			diag_ok = nde_kv_get_double(kv, key, &p->diag_K[c]) && diag_ok;
			g_snprintf(key, sizeof(key), "dO%d", c);
			diag_ok = nde_kv_get_double(kv, key, &p->diag_O[c]) && diag_ok;
		}
		/* The stored composite affine is the L1 fallback AND the whole
		 * parameter for DIRECT records — a blob without it is unusable. */
		if (!diag_ok) {
			nde_joint_group_calib_data_free(p);
			p = NULL;
		}
	}
	g_hash_table_unref(kv);
	return p;
}

static gchar *group_calib_log_hook(gpointer user, log_hook_detail detail) {
	const struct nde_joint_group_calib_data *p = user;
	(void)detail;
	const char *what = _("Group calibration");
	if (p) {
		what = p->kind == 0 ? _("Colour calibration (group)") :
		       p->kind == 1 ? _("Photometric CC (group)") :
		       p->kind == 2 ? _("Spectrophotometric CC (group)") :
		                      _("Group calibration");
	}
	if (!p)
		return g_strdup(what);
	GString *s = g_string_new(what);
	const char *sep = " — ";
	for (guint i = 0; i < p->n; i++) {
		if (p->parts[i].name && *p->parts[i].name) {
			g_string_append(s, sep);
			g_string_append(s, p->parts[i].name);
			sep = ", ";
		}
	}
	return g_string_free(s, FALSE);
}

static int group_calib_image_hook(struct generic_img_args *args, fits *fit,
                                  int nb_threads) {
	(void)args; (void)fit; (void)nb_threads;
	siril_log_message(_("flis.group_calibration cannot run as a generic image operation\n"));
	return 1;
}

const op_descriptor op_desc_flis_group_calibration = {
	.id = OP_GROUP_CALIBRATION, .version = 1,
	.image_hook = group_calib_image_hook,
	.log_hook = group_calib_log_hook,
	.description = N_("Group calibration"),
	.mem_ratio = 1.0f,
	.flags = 0,
	.serialize = group_calib_serialize,
	.deserialize = group_calib_deserialize,
};

/* ======================================================================= */
/* flis.register descriptor                                                */
/* ======================================================================= */

struct nde_joint_register_data *nde_joint_register_data_new(guint n) {
	struct nde_joint_register_data *p = g_new0(struct nde_joint_register_data, 1);
	p->destroy_fn = nde_joint_register_data_free;
	p->n = n;
	p->parts = g_new0(nde_joint_register_participant, n);
	for (guint k = 0; k < n; k++) {
		p->parts[k].H[0] = p->parts[k].H[4] = p->parts[k].H[8] = 1.0;
		p->parts[k].geom_sig = g_strdup("");
	}
	return p;
}

void nde_joint_register_data_free(void *ptr) {
	struct nde_joint_register_data *p = ptr;
	if (!p)
		return;
	for (guint i = 0; i < p->n; i++) {
		g_free(p->parts[i].name);
		g_free(p->parts[i].geom_sig);
	}
	g_free(p->parts);
	g_free(p);
}

static gchar *register_serialize(gconstpointer user) {
	const struct nde_joint_register_data *p = user;
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "n", p->n);
	nde_kv_add_int(kv, "method", p->method);
	nde_kv_add_int(kv, "tx", p->tx_type);
	nde_kv_add_int(kv, "interp", p->interpolation);
	nde_kv_add_bool(kv, "clamp", p->clamp);
	nde_kv_add_int(kv, "ref", p->ref_item);
	nde_kv_add_int(kv, "sel_x", p->selection.x);
	nde_kv_add_int(kv, "sel_y", p->selection.y);
	nde_kv_add_int(kv, "sel_w", p->selection.w);
	nde_kv_add_int(kv, "sel_h", p->selection.h);
	nde_kv_add_int(kv, "cw", p->canvas_w);
	nde_kv_add_int(kv, "ch", p->canvas_h);
	for (guint k = 0; k < p->n; k++) {
		const nde_joint_register_participant *part = &p->parts[k];
		gchar key[32];
		g_snprintf(key, sizeof(key), "i%u_item", k);
		nde_kv_add_int(kv, key, part->item_id);
		g_snprintf(key, sizeof(key), "i%u_name", k);
		nde_kv_add_str(kv, key, part->name ? part->name : "");
		for (int c = 0; c < 9; c++) {
			g_snprintf(key, sizeof(key), "i%u_h%d", k, c);
			nde_kv_add_double(kv, key, part->H[c]);
		}
		g_snprintf(key, sizeof(key), "i%u_x", k);
		nde_kv_add_int(kv, key, part->pos_x);
		g_snprintf(key, sizeof(key), "i%u_y", k);
		nde_kv_add_int(kv, key, part->pos_y);
		g_snprintf(key, sizeof(key), "i%u_rx", k);
		nde_kv_add_int(kv, key, part->out_rx);
		g_snprintf(key, sizeof(key), "i%u_ry", k);
		nde_kv_add_int(kv, key, part->out_ry);
		g_snprintf(key, sizeof(key), "i%u_geom", k);
		nde_kv_add_str(kv, key, part->geom_sig ? part->geom_sig : "");
	}
	return nde_kv_end(kv);
}

/* Parse participant @k.  FALSE rejects the WHOLE blob, same rule as the
 * scalar joint ops: a half-readable participant list would replay a warp
 * against the wrong layer, which destroys pixels rather than mis-scaling
 * them.  The output size is the one hard requirement — without it there is
 * nothing to warp into and no way to tell a truncated blob from a valid one. */
static gboolean register_participant_parse(GHashTable *kv, guint k,
                                           nde_joint_register_participant *part) {
	gchar key[32];
	gint64 v;
	g_snprintf(key, sizeof(key), "i%u_item", k);
	if (!nde_kv_get_int(kv, key, &v) || v <= 0)
		return FALSE;
	part->item_id = (gint)v;
	g_snprintf(key, sizeof(key), "i%u_name", k);
	part->name = g_strdup(nde_kv_get_str(kv, key) ? nde_kv_get_str(kv, key) : "");
	for (int c = 0; c < 9; c++) {
		g_snprintf(key, sizeof(key), "i%u_h%d", k, c);
		if (!nde_kv_get_double(kv, key, &part->H[c]) ||
		    !isfinite(part->H[c]))
			return FALSE;
	}
	g_snprintf(key, sizeof(key), "i%u_x", k);
	if (!nde_kv_get_int(kv, key, &v))
		return FALSE;
	part->pos_x = (gint)v;
	g_snprintf(key, sizeof(key), "i%u_y", k);
	if (!nde_kv_get_int(kv, key, &v))
		return FALSE;
	part->pos_y = (gint)v;
	g_snprintf(key, sizeof(key), "i%u_rx", k);
	if (!nde_kv_get_int(kv, key, &v) || v <= 0 || v > NDE_JOINT_MAX_DIM)
		return FALSE;
	part->out_rx = (gint)v;
	g_snprintf(key, sizeof(key), "i%u_ry", k);
	if (!nde_kv_get_int(kv, key, &v) || v <= 0 || v > NDE_JOINT_MAX_DIM)
		return FALSE;
	part->out_ry = (gint)v;
	g_snprintf(key, sizeof(key), "i%u_geom", k);
	const char *sig = nde_kv_get_str(kv, key);
	part->geom_sig = g_strdup(sig ? sig : "");
	return TRUE;
}

static gpointer register_deserialize(const gchar *blob, int version) {
	if (version > op_desc_flis_register.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	guint n;
	struct nde_joint_register_data *p = NULL;
	/* Registering fewer than two layers is not a registration. */
	if (participant_count(kv, 2, &n)) {
		p = nde_joint_register_data_new(n);
		for (guint k = 0; p && k < n; k++) {
			g_free(p->parts[k].geom_sig);
			p->parts[k].geom_sig = NULL;
			if (!register_participant_parse(kv, k, &p->parts[k])) {
				nde_joint_register_data_free(p);
				p = NULL;
			}
		}
		/* A layer listed twice would be warped twice. */
		for (guint i = 0; p && i < n; i++)
			for (guint j = i + 1; p && j < n; j++)
				if (p->parts[i].item_id == p->parts[j].item_id) {
					nde_joint_register_data_free(p);
					p = NULL;
				}
	}
	if (p) {
		gint64 v;
		if (nde_kv_get_int(kv, "method", &v)) p->method = (gint)v;
		if (nde_kv_get_int(kv, "tx", &v))     p->tx_type = (gint)v;
		if (nde_kv_get_int(kv, "interp", &v)) p->interpolation = (gint)v;
		nde_kv_get_bool(kv, "clamp", &p->clamp);
		if (nde_kv_get_int(kv, "ref", &v))    p->ref_item = (gint)v;
		if (nde_kv_get_int(kv, "sel_x", &v))  p->selection.x = (int)v;
		if (nde_kv_get_int(kv, "sel_y", &v))  p->selection.y = (int)v;
		if (nde_kv_get_int(kv, "sel_w", &v))  p->selection.w = (int)v;
		if (nde_kv_get_int(kv, "sel_h", &v))  p->selection.h = (int)v;
		if (nde_kv_get_int(kv, "cw", &v))     p->canvas_w = (gint)v;
		if (nde_kv_get_int(kv, "ch", &v))     p->canvas_h = (gint)v;
		/* The reference must be one of the participants: the re-run builds
		 * its internal sequence around it, and a reference outside the list
		 * has no pixels to align to. */
		gboolean ref_ok = FALSE;
		for (guint k = 0; k < p->n; k++)
			ref_ok = ref_ok || p->parts[k].item_id == p->ref_item;
		if (!ref_ok) {
			nde_joint_register_data_free(p);
			p = NULL;
		}
	}
	g_hash_table_unref(kv);
	return p;
}

static gchar *register_log_hook(gpointer user, log_hook_detail detail) {
	const struct nde_joint_register_data *p = user;
	(void)detail;
	if (!p)
		return g_strdup(_("Register layers"));
	GString *s = g_string_new(_("Register layers"));
	const char *sep = " — ";
	for (guint i = 0; i < p->n; i++) {
		if (p->parts[i].name && *p->parts[i].name) {
			g_string_append(s, sep);
			g_string_append(s, p->parts[i].name);
			sep = ", ";
		}
	}
	if (!g_strcmp0(sep, " — "))
		g_string_append_printf(s, _(" (%u layers)"), p->n);
	return g_string_free(s, FALSE);
}

/* Never reached — see layers_match_image_hook. */
static int register_image_hook(struct generic_img_args *args, fits *fit,
                               int nb_threads) {
	(void)args; (void)fit; (void)nb_threads;
	siril_log_message(_("flis.register cannot run as a generic image operation\n"));
	return 1;
}

const op_descriptor op_desc_flis_register = {
	.id = OP_REGISTER, .version = 1,
	.image_hook = register_image_hook,
	.log_hook = register_log_hook,
	.description = N_("Register layers"),
	.mem_ratio = 2.0f,
	/* Registration resamples every participant into a new bounding box, so
	 * a SECOND registration upstream of this one invalidates its stored
	 * transforms exactly as a rotation or a crop would.  The flag is what
	 * nde_joint_geometry_signature reads to notice that. */
	.flags = OP_GEOMETRY_CHANGING,
	.serialize = register_serialize,
	.deserialize = register_deserialize,
};

gboolean nde_joint_register_apply_settings(struct nde_joint_register_data *p,
                                           gint method, gint tx_type,
                                           gint interpolation, gboolean clamp,
                                           gint ref_item, gchar **err) {
	if (err) *err = NULL;
	if (!p)
		return FALSE;
	/* The reference has to stay inside the participant list for the same
	 * reason the deserializer insists on it: the re-solve builds its internal
	 * sequence around the reference, and one outside the list has no pixels to
	 * align to.  Refusing here keeps a rejected amend from writing a blob its
	 * own deserializer would then throw away. */
	gboolean ref_ok = FALSE;
	for (guint k = 0; k < p->n; k++)
		ref_ok = ref_ok || p->parts[k].item_id == ref_item;
	if (!ref_ok) {
		if (err)
			*err = g_strdup(_("the reference layer must be one of the "
			                  "registered layers"));
		return FALSE;
	}

	/* The selection is a SETTING of this record, not live state — replay can
	 * happen much later, headless, with com.selection empty or somewhere else
	 * entirely.  So a method change is the moment to capture it.
	 *
	 * Which selection to store takes some care.  The record may well not have
	 * one: a method that requires none (the global star alignments) stores
	 * {0,0,0,0}, so switching such a record to KOMBAT or DFT MUST pick up the
	 * selection the user has drawn — otherwise the re-solve runs with an empty
	 * rectangle, and an empty template is an ASSERT inside cv::matchTemplate,
	 * i.e. the process aborts.  But re-opening a KOMBAT record merely to
	 * change the interpolation must not silently drag its selection off to
	 * wherever the cursor last left one.  Hence: prefer a valid live
	 * selection, else keep a valid stored one, else refuse. */
	selection_type sel_req = REQUIRES_NO_SELECTION;
	rectangle new_sel = p->selection;
	if (flis_register_resolve_method((flis_reg_method_id)method, &sel_req, NULL)
	    && sel_req != REQUIRES_NO_SELECTION) {
		if (flis_register_selection_ok(sel_req, &com.selection, 0, 0, NULL)) {
			new_sel = com.selection;
		} else if (!flis_register_selection_ok(sel_req, &p->selection, 0, 0, NULL)) {
			gchar *why = NULL;
			flis_register_selection_ok(sel_req, &com.selection, 0, 0, &why);
			if (err)
				*err = g_strdup_printf(_("%s — the previous method needed no "
				                         "selection, so this record has none "
				                         "stored"), why ? why : "?");
			g_free(why);
			return FALSE;
		}
	}

	const gboolean changed = p->method        != method ||
	                         p->tx_type       != tx_type ||
	                         p->interpolation != interpolation ||
	                         (!p->clamp)      != (!clamp) ||
	                         p->ref_item      != ref_item ||
	                         memcmp(&p->selection, &new_sel, sizeof new_sel) != 0;
	p->selection     = new_sel;
	p->method        = method;
	p->tx_type       = tx_type;
	p->interpolation = interpolation;
	p->clamp         = clamp;
	p->ref_item      = ref_item;
	if (!changed)
		return FALSE;   /* byte-identical re-serialization: OK is a no-op */

	/* The stored transforms ARE the old settings' answer, so once the settings
	 * move they are stale by definition — no upstream geometry had to change
	 * for that to be true.  Replay picks L1 vs L2 by comparing each
	 * participant's geom_sig against a freshly computed signature
	 * (register_geometry_unchanged), so poisoning the signatures is how this
	 * layer of the code says "do not reuse these": the comparison cannot
	 * match, replay takes L2, and the registration is honestly re-solved with
	 * what the user just chose.  The H values are left in place as the L3
	 * fallback for when the re-solve itself cannot run.
	 *
	 * The marker is deliberately NOT the empty string.  "" is a REAL signature
	 * — the one a layer with no geometry upstream of the record computes — and
	 * that is the common case for a first registration, so clearing to ""
	 * would compare equal and take L1, the exact opposite of the intent.  Any
	 * value a signature cannot take will do; nde_joint_geometry_signature
	 * returns "" or SHA256 hex, so a non-hex token is safe forever. */
	for (guint k = 0; k < p->n; k++) {
		g_free(p->parts[k].geom_sig);
		p->parts[k].geom_sig = g_strdup(NDE_JOINT_GEOM_SIG_STALE);
	}
	return TRUE;
}

/* ======================================================================= */
/* flis.register replay                                                    */
/* ======================================================================= */

/* Shared with the scalar joint ops' factor computation below. */
static gboolean resolve_participants_ids(const struct nde_record *rec,
                                         const gint *item_ids,
                                         gchar * const *names, guint n,
                                         fits *self_scratch, gint self_item,
                                         gboolean require_mono,
                                         fits **pixels, gboolean *owned,
                                         gchar **err);

/* How many times an analysis has actually run (cache misses), across all three
 * joint ops.  A test counter: each generation cache's whole point is that one
 * edit costs one analysis however many participants replay, and only this can
 * prove it.  Shared deliberately — the tests measure deltas around a single
 * op, so one counter reads the same as three. */
static guint analysis_runs;

/* The solutions memo (declared with the factor memo above) exists for the same
 * reason: the L2 re-solve runs a whole registration, so N participant replays
 * of one edit must share ONE solve.  Unlike the factors, n is known before the
 * lookup, so the stored size is checked against it. */
static flis_reg_solution *register_solution_get(gint64 record_id, guint64 generation,
                                                guint n) {
	gsize size = 0;
	flis_reg_solution *sol = nde_memo_get(&register_solution_memo, record_id,
	                                      generation, &size);
	if (sol && size != n * sizeof(*sol)) {
		g_free(sol);
		sol = NULL;
	}
	return sol;
}

static void register_solution_put(gint64 record_id, guint64 generation, guint n,
                                  const flis_reg_solution *sol) {
	nde_memo_put(&register_solution_memo, record_id, generation,
	             sol, n * sizeof(*sol));
}

/* TRUE when EVERY participant reaches the record with the geometry it had at
 * capture, i.e. the stored transforms are still the right answer (L1).  A
 * single mismatch condemns them all: the transforms are defined against ONE
 * shared reference frame, so a participant whose geometry moved invalidates
 * the whole solve, not just its own row of it. */
static gboolean register_geometry_unchanged(const struct nde_record *rec,
                                            const struct nde_joint_register_data *p) {
	for (guint k = 0; k < p->n; k++) {
		gchar *now = nde_joint_geometry_signature(p->parts[k].item_id,
		                                          rec->record_id);
		gboolean same = !g_strcmp0(now, p->parts[k].geom_sig ? p->parts[k].geom_sig : "");
		g_free(now);
		if (!same)
			return FALSE;
	}
	return TRUE;
}

/* The stored answer, verbatim — L1, and the L3 fallback when a re-solve that
 * SHOULD have run could not. */
static void register_solutions_from_record(const struct nde_joint_register_data *p,
                                           flis_reg_solution *sol) {
	for (guint k = 0; k < p->n; k++) {
		memcpy(sol[k].H, p->parts[k].H, sizeof(sol[k].H));
		sol[k].out_rx = p->parts[k].out_rx;
		sol[k].out_ry = p->parts[k].out_ry;
		sol[k].pos_x  = p->parts[k].pos_x;
		sol[k].pos_y  = p->parts[k].pos_y;
	}
}

/* L2: re-run the recorded registration against the participants as they now
 * stand.  FALSE leaves @sol untouched and the caller falls back to L1. */
static gboolean register_resolve(const struct nde_record *rec,
                                 const struct nde_joint_register_data *p,
                                 fits *self_scratch, gint self_item,
                                 flis_reg_solution *sol) {
	selection_type sel = REQUIRES_NO_SELECTION;
	transformation_type tx = HOMOGRAPHY_TRANSFORMATION;
	registration_function method =
			flis_register_resolve_method((flis_reg_method_id)p->method, &sel, &tx);
	if (!method) {
		siril_log_warning(_("Register layers replay: the recorded method is not "
		                    "available — replaying the recorded transforms\n"));
		return FALSE;
	}
	/* The transformation type is a RECORDED setting, not a property of the
	 * method: the user could have solved a shift-only alignment with a
	 * homography-capable method.  Honour what the record says. */
	tx = (transformation_type)p->tx_type;

	guint n = p->n;
	fits **pixels = g_new0(fits *, n);
	gboolean *owned = g_new0(gboolean, n);
	gint *ids = g_new0(gint, n);
	gchar **names = g_new0(gchar *, n);
	for (guint k = 0; k < n; k++) {
		ids[k] = p->parts[k].item_id;
		names[k] = p->parts[k].name;
	}
	gboolean ok = FALSE;
	gchar *inner = NULL;
	if (!resolve_participants_ids(rec, ids, names, n, self_scratch, self_item,
	                              FALSE, pixels, owned, &inner)) {
		siril_log_warning(_("Register layers replay: %s — replaying the "
		                    "recorded transforms\n"), inner ? inner : "?");
		g_free(inner);
		goto out;
	}

	gint ref_index = 0, canvas_index = -1;
	for (guint k = 0; k < n; k++)
		if (p->parts[k].item_id == p->ref_item)
			ref_index = (gint)k;
	/* The canvas (base) layer is the document's FIRST layer; it only fixes
	 * the canvas origin when it is one of the participants. */
	flis_layer_t *base = com.uniq && com.uniq->layers ?
			(flis_layer_t *)com.uniq->layers->data : NULL;
	if (base)
		for (guint k = 0; k < n; k++)
			if (p->parts[k].item_id == base->item_id)
				canvas_index = (gint)k;
	/* Where the reference sat before the registration, for the no-base-layer
	 * branch of the origin choice.  The record has it: the reference's own
	 * capture-time position minus nothing, because a reference that is not
	 * the base layer keeps whatever offset it already had. */
	flis_layer_t *ref_lay = flis_layer_get_by_id(p->ref_item);
	gint ref_pos_x = 0, ref_pos_y = 0;
	if (nde_checkpoint_baseline_position(p->ref_item, &ref_pos_x, &ref_pos_y)) {
		/* recorded starting position — the anchor the live path used */
	} else if (ref_lay) {
		ref_pos_x = ref_lay->position_x;
		ref_pos_y = ref_lay->position_y;
	}

	analysis_runs++;
	if (flis_register_solve(pixels, (int)n, ref_index, method, sel, tx,
	                        p->selection, canvas_index, ref_pos_x, ref_pos_y,
	                        sol)) {
		siril_log_warning(_("Register layers replay: the registration could not "
		                    "be re-run — replaying the recorded transforms\n"));
		goto out;
	}
	ok = TRUE;
out:
	for (guint k = 0; k < n; k++) {
		if (owned[k] && pixels[k]) {
			clearfits(pixels[k]);
			free(pixels[k]);
		}
	}
	g_free(pixels);
	g_free(owned);
	g_free(ids);
	g_free(names);
	return ok;
}

gboolean nde_joint_register_apply(const struct nde_record *rec, nde_state *state,
                                  gint item_id, gchar **err) {
	g_return_val_if_fail(rec != NULL && err != NULL, FALSE);
	g_return_val_if_fail(state != NULL && state->pix != NULL, FALSE);
	fits *scratch = state->pix;
	struct nde_joint_register_data *p =
			register_deserialize(rec->params, rec->op_version);
	if (!p) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
		return FALSE;
	}
	gint idx = -1;
	for (guint k = 0; k < p->n; k++)
		if (p->parts[k].item_id == item_id)
			idx = (gint)k;
	if (idx < 0) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) does not name this layer as a participant"),
		                       rec->record_id, rec->op_id);
		nde_joint_register_data_free(p);
		return FALSE;
	}

	guint64 generation = nde_history_generation();
	flis_reg_solution *sol = register_solution_get(rec->record_id, generation, p->n);
	if (!sol) {
		sol = g_new0(flis_reg_solution, p->n);
		register_solutions_from_record(p, sol);
		/* The intent ladder.  L1 first because it is both cheaper AND more
		 * faithful: re-solving when nothing moved would re-fit the stars and
		 * shift every layer by a fraction of a pixel for no reason. */
		if (!register_geometry_unchanged(rec, p)) {
			flis_reg_solution *fresh = g_new0(flis_reg_solution, p->n);
			if (register_resolve(rec, p, scratch, item_id, fresh)) {
				g_free(sol);
				sol = fresh;
			} else {
				g_free(fresh);   /* L3: the recorded transforms, warned about */
			}
		}
		register_solution_put(rec->record_id, generation, p->n, sol);
	}

	const flis_reg_solution *s = &sol[idx];
	Homography H = { 0 };
	H.h00 = s->H[0]; H.h01 = s->H[1]; H.h02 = s->H[2];
	H.h10 = s->H[3]; H.h11 = s->H[4]; H.h12 = s->H[5];
	H.h20 = s->H[6]; H.h21 = s->H[7]; H.h22 = s->H[8];
	int rc = cvTransformImage(scratch, (unsigned int)s->out_rx,
	                          (unsigned int)s->out_ry, H, 1.f,
	                          p->interpolation, p->clamp, NULL);
	if (rc) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): the layer "
		                         "could not be warped"),
		                       rec->record_id, rec->op_id);
		g_free(sol);
		nde_joint_register_data_free(p);
		return FALSE;
	}
	invalidate_stats_from_fit(scratch);
	/* Registration moves the layer as well as warping it, so the value it
	 * produces is the warped pixels AT the new position — only if this replay
	 * is carrying one, though; the verification path is not (nde_state.h). */
	if (state->has_pos) {
		state->pos_x = s->pos_x;
		state->pos_y = s->pos_y;
	}
	g_free(sol);
	nde_joint_register_data_free(p);
	return TRUE;
}

/* ======================================================================= */
/* Factor computation                                                      */
/* ======================================================================= */

/* Resolve every participant's pixels at the record's position.  The item
 * whose chain is being replayed contributes its accumulated scratch (the
 * same "the base is already resolved" rule composite_apply follows); every
 * other participant is a live positional edge, replayed strictly up to just
 * before the record.  On success fills @pixels (borrowed for self, owned
 * otherwise, @owned says which); on failure everything owned is freed. */
/* @item_ids / @names are the participant list flattened out of whichever
 * params struct the caller has (the two scalar joint ops share
 * nde_joint_participant; flis.register has its own).  @require_mono is the
 * colour-model check the SCALAR ops need — their analysis is defined over
 * tinted mono layers — and which registration must not apply, because
 * warping an RGB layer is perfectly well defined. */
static gboolean resolve_participants_ids(const struct nde_record *rec,
                                         const gint *item_ids,
                                         gchar * const *names, guint n,
                                         fits *self_scratch, gint self_item,
                                         gboolean require_mono,
                                         fits **pixels, gboolean *owned,
                                         gchar **err) {
	for (guint k = 0; k < n; k++) {
		const char *name = names && names[k] ? names[k] : "?";
		if (item_ids[k] == self_item && self_scratch) {
			pixels[k] = self_scratch;
		} else {
			pixels[k] = nde_replay_resolve_before(item_ids[k],
			                                      rec->record_id, err);
			if (!pixels[k]) {
				gchar *inner = *err;
				*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): "
				                         "layer '%s' cannot be rebuilt: %s"),
				                       rec->record_id, rec->op_id, name,
				                       inner ? inner : "?");
				g_free(inner);
				goto fail;
			}
			owned[k] = TRUE;
		}
		if (require_mono && pixels[k]->naxes[2] != 1) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): "
			                         "layer '%s' is no longer a mono layer"),
			                       rec->record_id, rec->op_id, name);
			goto fail;
		}
	}
	return TRUE;
fail:
	for (guint k = 0; k < n; k++) {
		if (owned[k] && pixels[k]) {
			clearfits(pixels[k]);
			free(pixels[k]);
			pixels[k] = NULL;
			owned[k] = FALSE;
		}
	}
	return FALSE;
}

static gboolean resolve_participants(const struct nde_record *rec,
                                     const nde_joint_participant *parts, guint n,
                                     fits *self_scratch, gint self_item,
                                     fits **pixels, gboolean *owned, gchar **err) {
	gint *ids = g_new0(gint, n);
	gchar **names = g_new0(gchar *, n);
	for (guint k = 0; k < n; k++) {
		ids[k] = parts[k].item_id;
		names[k] = parts[k].name;
	}
	gboolean ok = resolve_participants_ids(rec, ids, names, n, self_scratch,
	                                       self_item, TRUE, pixels, owned, err);
	g_free(ids);
	g_free(names);
	return ok;
}

/* The tint each participant carries AS OF the record's log position: the
 * positional fold when the log says anything about the layer's tint before
 * the record, the tint recorded in the joint params otherwise (a document
 * from before tint capture existed).  Fills @tints as n RGB triples with
 * 1/1/1 for an untinted layer. */
static void fold_participant_tints(const struct nde_record *rec,
                                   const nde_joint_participant *parts, guint n,
                                   double *tints) {
	for (guint k = 0; k < n; k++) {
		gboolean tinted;
		double tint[3];
		if (nde_compositing_has_tint_record_upto(parts[k].item_id, rec->record_id)) {
			nde_compositing_fold_upto(parts[k].item_id, rec->record_id,
			                          NULL, NULL, NULL, &tinted, tint);
		} else {
			tinted = parts[k].tinted;
			tint[0] = parts[k].tint[0];
			tint[1] = parts[k].tint[1];
			tint[2] = parts[k].tint[2];
		}
		tints[k * 3 + 0] = tinted ? tint[0] : 1.0;
		tints[k * 3 + 1] = tinted ? tint[1] : 1.0;
		tints[k * 3 + 2] = tinted ? tint[2] : 1.0;
	}
}

/* Background median of each resolved participant, normalised to [0,1].
 * SINGLE_THREADED like the live analysis, so replay is deterministic. */
static gboolean participant_medians(const struct nde_record *rec,
                                    const nde_joint_participant *parts, guint n,
                                    fits **pixels, double *medians, gchar **err) {
	for (guint k = 0; k < n; k++) {
		imstats *st = statistics(NULL, -1, pixels[k], 0, NULL, STATS_BASIC,
		                         SINGLE_THREADED);
		if (!st) {
			*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): "
			                         "statistics failed on layer '%s'"),
			                       rec->record_id, rec->op_id,
			                       parts[k].name ? parts[k].name : "?");
			return FALSE;
		}
		double med = st->median;
		free_stats(st);
		if (pixels[k]->type == DATA_USHORT)
			med /= USHRT_MAX_DOUBLE;
		medians[k] = med;
	}
	return TRUE;
}

guint nde_joint_analysis_runs(void) {
	return analysis_runs;
}

/* The group composite as of the record's position: synthetic layer shells
 * over the resolved pixels (the borrowed-struct pattern nde_composite.c's
 * replay established), compositing state from the positional fold.  NULL
 * when the render fails or does not composite to colour. */
static fits *render_group_composite(const struct nde_record *rec,
                                    const struct nde_joint_group_calib_data *p,
                                    fits **pixels) {
	flis_layer_t *lays = g_new0(flis_layer_t, p->n);
	GSList *list = NULL;
	for (guint i = 0; i < p->n; i++) {
		gfloat opacity;
		gint blend;
		gboolean visible, tinted;
		double tint[3];
		nde_compositing_fold_upto(p->parts[i].item_id, rec->record_id,
		                          &opacity, &blend, &visible, &tinted, tint);
		if (!nde_compositing_has_tint_record_upto(p->parts[i].item_id,
		                                          rec->record_id)) {
			tinted = p->parts[i].tinted;
			tint[0] = p->parts[i].tint[0];
			tint[1] = p->parts[i].tint[1];
			tint[2] = p->parts[i].tint[2];
		}
		flis_layer_t *live = flis_layer_get_by_id(p->parts[i].item_id);
		lays[i].fit        = pixels[i];
		lays[i].item_id    = p->parts[i].item_id;
		lays[i].blend_mode = (flis_blend_mode_t)blend;
		lays[i].opacity    = opacity;
		lays[i].visible    = visible;
		lays[i].has_tint   = tinted;
		lays[i].layer_tint = (flis_tint_t){ tint[0], tint[1], tint[2] };
		lays[i].position_x = live ? live->position_x : 0;
		lays[i].position_y = live ? live->position_y : 0;
		lays[i].group_id   = live ? live->group_id : p->group_id;
		list = g_slist_append(list, &lays[i]);
	}
	flis_render_ctx ctx = { 0 };
	ctx.canvas_w = com.uniq ? com.uniq->canvas_w : 0;
	ctx.canvas_h = com.uniq ? com.uniq->canvas_h : 0;
	ctx.groups   = com.uniq ? com.uniq->groups : NULL;
	/* Same sub-composite contract as the live analysis: screen stack over
	 * black, no canvas background (flis_group_render_calib_composite). */
	fits *out = flis_render_layers_ctx(list, &ctx, TRUE, FALSE);
	g_slist_free(list);
	g_free(lays);
	if (out && out->naxes[2] != 3) {
		clearfits(out);
		free(out);
		out = NULL;
	}
	/* WCS donor for the photometric pipeline, exactly as the live analysis
	 * transplants one (flis_group_render_calib_composite): a solved
	 * participant first, then any solved document layer with canvas
	 * geometry. */
	if (out) {
		flis_layer_t *donor = NULL;
		for (guint i = 0; i < p->n && !donor; i++) {
			flis_layer_t *lay = flis_layer_get_by_id(p->parts[i].item_id);
			if (lay && lay->fit && lay->fit->keywords.wcslib &&
			    lay->fit->rx == out->rx && lay->fit->ry == out->ry &&
			    lay->position_x == 0 && lay->position_y == 0)
				donor = lay;
		}
		for (GSList *l = com.uniq ? com.uniq->layers : NULL; l && !donor; l = l->next) {
			flis_layer_t *lay = l->data;
			if (lay && lay->fit && lay->fit->keywords.wcslib &&
			    lay->fit->rx == out->rx && lay->fit->ry == out->ry &&
			    lay->position_x == 0 && lay->position_y == 0)
				donor = lay;
		}
		if (donor)
			copy_fits_metadata(donor->fit, out);
	}
	return out;
}

/* The composite-level affine (K, O) the record's parameters describe.  The
 * intent ladder: L2 re-runs the recorded analysis against the resolved
 * composite; when it cannot run, L1 falls back to the stored diag_K/diag_O
 * — the DISTRIBUTION below still recomputes against current tints and
 * medians, so even L1 adapts to the narrowband tint edit. */
static void group_calib_derive_KO(const struct nde_record *rec,
                                  const struct nde_joint_group_calib_data *p,
                                  fits **pixels, double K[3], double O[3]) {
	memcpy(K, p->diag_K, 3 * sizeof(double));
	memcpy(O, p->diag_O, 3 * sizeof(double));
	if (p->kind == NDE_JOINT_GROUP_DIRECT)
		return;   /* the affine IS the operation parameter */
	if (p->kind == 0 && p->manual) {
		for (int c = 0; c < 3; c++) {
			K[c] = p->manual_kw[c];
			O[c] = 0.0;
		}
		return;
	}
	if (p->kind == 0) {
		fits *composite = render_group_composite(rec, p, pixels);
		if (!composite) {
			siril_log_warning(_("Group calibration replay: the composite could not be re-rendered — using the stored calibration\n"));
			return;
		}
		double kw[3] = { 0 }, bg[3] = { 0 };
		get_coeff_for_wb(composite, p->white_sel, p->black_sel, kw, bg,
		                 get_normalized_value(composite), p->low, p->high);
		clearfits(composite);
		free(composite);
		if (kw[0] > 0.0 && kw[1] > 0.0 && kw[2] > 0.0) {
			/* calibrate(): C' = (C - bg)*kw + bg  ==  kw*C + bg*(1-kw) */
			for (int c = 0; c < 3; c++) {
				K[c] = kw[c];
				O[c] = bg[c] * (1.0 - kw[c]);
			}
		} else {
			siril_log_warning(_("Group calibration replay: white balance could not be recomputed from the selections — using the stored calibration\n"));
		}
		return;
	}
	/* PCC / SPCC: re-run the full photometric pipeline on the re-rendered
	 * composite.  The star data comes from the record's embedded catalogue
	 * (nde_cat.h) via the current-record id — offline and deterministic. */
	if (p->pcc_blob) {
		const op_descriptor *pop = op_descriptor_by_id("color.photometric_cc");
		struct photometric_cc_data *pcc = (pop && pop->deserialize) ?
				pop->deserialize(p->pcc_blob, pop->version) : NULL;
		gboolean derived = FALSE;
		if (pcc && !pcc->have_effective) {
			fits *composite = render_group_composite(rec, p, pixels);
			if (composite) {
				struct generic_img_args ia = { 0 };
				ia.user = pcc;
				ia.nde_replay = TRUE;   /* no capture, no stash */
				int rv = photometric_cc_image_hook(&ia, composite, com.max_thread);
				pcc->fit = NULL;
				if (!rv && pcc->have_effective) {
					/* apply_photometric_color_correction():
					 * C' = kw*C + (bg_mean - bg*kw) */
					double bg_mean = (pcc->eff_bg[0] + pcc->eff_bg[1] +
					                  pcc->eff_bg[2]) / 3.0;
					for (int c = 0; c < 3; c++) {
						K[c] = pcc->eff_kw[c];
						O[c] = bg_mean - pcc->eff_bg[c] * K[c];
					}
					derived = TRUE;
				}
				clearfits(composite);
				free(composite);
			}
		}
		if (pcc && pcc->destroy_fn)
			pcc->destroy_fn(pcc);
		if (derived)
			return;
	}
	siril_log_message(_("Group calibration replay: the photometric analysis could not be re-run — using the stored calibration\n"));
}

static gboolean group_factors_compute(const struct nde_record *rec,
                                      fits *self_scratch, gint self_item,
                                      guint64 generation,
                                      double **a_out, double **b_out,
                                      gchar **err) {
	struct nde_joint_group_calib_data *p =
			group_calib_deserialize(rec->params, rec->op_version);
	if (!p) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
		                       rec->record_id, rec->op_id);
		return FALSE;
	}
	guint n = p->n;
	fits **pixels = g_new0(fits *, n);
	gboolean *owned = g_new0(gboolean, n);
	double *tints = g_new0(double, n * 3);
	double *medians = g_new0(double, n);
	double *a = g_new0(double, n);
	double *b = g_new0(double, n);
	gboolean ok = FALSE;

	if (!resolve_participants(rec, p->parts, n, self_scratch, self_item,
	                          pixels, owned, err))
		goto out;
	fold_participant_tints(rec, p->parts, n, tints);
	if (!participant_medians(rec, p->parts, n, pixels, medians, err))
		goto out;

	analysis_runs++;
	double K[3], O[3];
	/* The sibling resolutions above recursed through the replay driver and
	 * cleared its current-record id; the analysis below (the PCC pipeline)
	 * reads THIS record's embedded catalogue under that id — restore it. */
	nde_replay_set_current_record(rec->record_id);
	group_calib_derive_KO(rec, p, pixels, K, O);
	gchar *why = NULL;
	if (flis_group_distribute(tints, medians, (int)n, K, O, a, b, &why)) {
		/* L0: the recorded per-layer factors, verbatim, with a warning —
		 * intent cannot be recomputed but the pixels stay reproducible. */
		siril_log_warning(_("Group calibration replay: %s — replaying the recorded per-layer factors\n"),
		                  why ? why : "?");
		g_free(why);
		for (guint k = 0; k < n; k++) {
			a[k] = p->parts[k].diag_scale;
			b[k] = p->parts[k].diag_offset;
		}
	}
	joint_factors_put(rec->record_id, generation, n, a, b);
	*a_out = a;
	*b_out = b;
	a = b = NULL;
	ok = TRUE;
out:
	for (guint k = 0; k < n; k++) {
		if (owned[k] && pixels[k]) {
			clearfits(pixels[k]);
			free(pixels[k]);
		}
	}
	g_free(pixels);
	g_free(owned);
	g_free(tints);
	g_free(medians);
	g_free(a);
	g_free(b);
	nde_joint_group_calib_data_free(p);
	return ok;
}

gboolean nde_joint_factors(const struct nde_record *rec, fits *self_scratch,
                           gint self_item, double **a_out, double **b_out,
                           gchar **err) {
	g_return_val_if_fail(rec != NULL && a_out != NULL && b_out != NULL && err != NULL, FALSE);
	*err = NULL;
	guint64 generation = nde_history_generation();
	if (joint_factors_get(rec->record_id, generation, a_out, b_out))
		return TRUE;

	if (!g_strcmp0(rec->op_id, OP_GROUP_CALIBRATION))
		return group_factors_compute(rec, self_scratch, self_item, generation,
		                             a_out, b_out, err);
	if (g_strcmp0(rec->op_id, OP_LAYERS_MATCH)) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): not a joint operation"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
		return FALSE;
	}

	struct nde_joint_layers_match_data *p =
			layers_match_deserialize(rec->params, rec->op_version);
	if (!p) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
		                       rec->record_id, rec->op_id);
		return FALSE;
	}
	guint n = p->n;
	fits **pixels = g_new0(fits *, n);
	gboolean *owned = g_new0(gboolean, n);
	double *tints = g_new0(double, n * 3);
	double *medians = g_new0(double, n);
	double *a = g_new0(double, n);
	double *b = g_new0(double, n);   /* layers match: offsets stay 0 */
	guint8 *infeasible = g_new0(guint8, n);
	gboolean ok = FALSE;

	if (!resolve_participants(rec, p->parts, n, self_scratch, self_item,
	                          pixels, owned, err))
		goto out;
	fold_participant_tints(rec, p->parts, n, tints);
	if (!participant_medians(rec, p->parts, n, pixels, medians, err))
		goto out;

	analysis_runs++;
	gchar *why = NULL;
	int rc = flis_layers_match_solve(tints, medians, (int)n, a, infeasible, &why);
	if (rc) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): %s"),
		                       rec->record_id, rec->op_id, why ? why : "?");
		g_free(why);
		goto out;
	}
	for (guint k = 0; k < n; k++) {
		if (infeasible[k])
			siril_log_warning(_("Layers match: layer '%s' left unscaled — tint infeasible\n"),
			                  p->parts[k].name ? p->parts[k].name : "?");
	}
	joint_factors_put(rec->record_id, generation, n, a, b);
	*a_out = a;
	*b_out = b;
	a = b = NULL;   /* transferred */
	ok = TRUE;
out:
	for (guint k = 0; k < n; k++) {
		if (owned[k] && pixels[k]) {
			clearfits(pixels[k]);
			free(pixels[k]);
		}
	}
	g_free(pixels);
	g_free(owned);
	g_free(tints);
	g_free(medians);
	g_free(a);
	g_free(b);
	g_free(infeasible);
	nde_joint_layers_match_data_free(p);
	return ok;
}

gboolean nde_joint_factor_for_item(const struct nde_record *rec,
                                   fits *self_scratch, gint self_item,
                                   double *scale_out, double *offset_out,
                                   gchar **err) {
	g_return_val_if_fail(scale_out != NULL && offset_out != NULL && err != NULL, FALSE);
	guint n = 0;
	gint *items = nde_joint_record_participants(rec, &n);
	if (!items) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s): parameters failed to parse"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
		return FALSE;
	}
	gint idx = -1;
	for (guint k = 0; k < n; k++)
		if (items[k] == self_item)
			idx = (gint)k;
	g_free(items);
	if (idx < 0) {
		*err = g_strdup_printf(_("record %" G_GINT64_FORMAT " (%s) does not name this layer as a participant"),
		                       rec->record_id, rec->op_id ? rec->op_id : "?");
		return FALSE;
	}
	double *a = NULL, *b = NULL;
	if (!nde_joint_factors(rec, self_scratch, self_item, &a, &b, err))
		return FALSE;
	*scale_out = a[idx];
	*offset_out = b[idx];
	g_free(a);
	g_free(b);
	return TRUE;
}
