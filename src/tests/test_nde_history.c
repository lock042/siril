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
 * test_nde_history — NDE provenance core (phase 1, sketch §18):
 *   - by-id registry lookups over the op_descriptor list
 *   - key=value codec escaping round-trips
 *   - append / live_count / on_undo / on_redo / truncate-on-append semantics
 *   - snapshot deep-copy isolation
 *   - worker-thread append vs snapshot smoke test
 */

#include <criterion/criterion.h>
#include <string.h>
#include "core/siril.h"
#include "core/op_descriptor.h"
#include "core/op_descriptors.h"
#include "core/nde_history.h"
#include "algos/geometry.h"

cominfo com;	// the core data struct
fits *gfit;	// currently loaded image

static void setup(void) {
	memset(&com, 0, sizeof(com));
	com.uniq = g_new0(single, 1);
}

static void teardown(void) {
	nde_history_attach(NULL);
	g_free(com.uniq);
	com.uniq = NULL;
}

TestSuite(nde_history, .init = setup, .fini = teardown);

/* Convenience: append a minimal record with the given op id. */
static gint64 append_op(const char *op_id) {
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup(op_id);
	rec->op_version = 1;
	rec->tier = NDE_TIER_B;
	rec->summary = g_strdup(op_id);
	return nde_history_append(rec);
}

/* ---------------- registry ---------------- */

Test(nde_history, registry_by_id) {
	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	cr_assert(n > 0);
	/* every descriptor must be findable by its own id */
	for (size_t i = 0; i < n; i++)
		cr_assert_eq(op_descriptor_by_id(all[i]->id), all[i],
		             "lookup of '%s' returned the wrong descriptor", all[i]->id);
	cr_assert_null(op_descriptor_by_id("no.such_op"));
	cr_assert_null(op_descriptor_by_id(""));
	cr_assert_null(op_descriptor_by_id(NULL));
}

/* ---------------- kv codec ---------------- */

Test(nde_history, kv_roundtrip_plain) {
	GString *kv = nde_kv_start();
	nde_kv_add_str(kv, "name", "value");
	nde_kv_add_int(kv, "count", -42);
	nde_kv_add_bool(kv, "flag", TRUE);
	nde_kv_add_float(kv, "ratio", 1.5f);
	gchar *blob = nde_kv_end(kv);
	cr_assert_str_eq(blob, "name=value;count=-42;flag=1;ratio=1.5");

	GHashTable *parsed = nde_kv_parse(blob);
	cr_assert_str_eq(nde_kv_get_str(parsed, "name"), "value");
	gint64 i64;
	cr_assert(nde_kv_get_int(parsed, "count", &i64));
	cr_assert_eq(i64, -42);
	gboolean b;
	cr_assert(nde_kv_get_bool(parsed, "flag", &b));
	cr_assert(b);
	float f;
	cr_assert(nde_kv_get_float(parsed, "ratio", &f));
	cr_assert_float_eq(f, 1.5f, 1e-9);
	cr_assert(!nde_kv_get_int(parsed, "absent", &i64));
	cr_assert_null(nde_kv_get_str(parsed, "absent"));
	g_hash_table_unref(parsed);
	g_free(blob);
}

Test(nde_history, kv_roundtrip_escaping) {
	static const char *nasty[] = {
		"semi;colon", "eq=uals", "back\\slash", "new\nline",
		";;==\\\\\n\n", "trailing\\", "", "plain",
	};
	GString *kv = nde_kv_start();
	for (size_t i = 0; i < G_N_ELEMENTS(nasty); i++) {
		gchar *key = g_strdup_printf("k%zu", i);
		nde_kv_add_str(kv, key, nasty[i]);
		g_free(key);
	}
	gchar *blob = nde_kv_end(kv);
	GHashTable *parsed = nde_kv_parse(blob);
	for (size_t i = 0; i < G_N_ELEMENTS(nasty); i++) {
		gchar *key = g_strdup_printf("k%zu", i);
		cr_assert_str_eq(nde_kv_get_str(parsed, key), nasty[i],
		                 "value %zu did not round-trip", i);
		g_free(key);
	}
	g_hash_table_unref(parsed);
	g_free(blob);
}

Test(nde_history, kv_float_precision) {
	/* %.9g must round-trip any float32 exactly */
	const float values[] = { 0.1f, 1.0f/3.0f, 3.4028235e38f, 1.17549435e-38f, -0.0f };
	GString *kv = nde_kv_start();
	for (size_t i = 0; i < G_N_ELEMENTS(values); i++) {
		gchar *key = g_strdup_printf("f%zu", i);
		nde_kv_add_float(kv, key, values[i]);
		g_free(key);
	}
	gchar *blob = nde_kv_end(kv);
	GHashTable *parsed = nde_kv_parse(blob);
	for (size_t i = 0; i < G_N_ELEMENTS(values); i++) {
		gchar *key = g_strdup_printf("f%zu", i);
		float f;
		cr_assert(nde_kv_get_float(parsed, key, &f));
		cr_assert(memcmp(&f, &values[i], sizeof(float)) == 0,
		          "float %zu did not round-trip bit-exactly", i);
		g_free(key);
	}
	g_hash_table_unref(parsed);
	g_free(blob);
}

Test(nde_history, kv_parse_edge_cases) {
	GHashTable *parsed = nde_kv_parse(NULL);
	cr_assert_eq(g_hash_table_size(parsed), 0);
	g_hash_table_unref(parsed);

	parsed = nde_kv_parse("");
	cr_assert_eq(g_hash_table_size(parsed), 0);
	g_hash_table_unref(parsed);

	/* malformed pair (no '=') is skipped, valid neighbours survive */
	parsed = nde_kv_parse("orphan;a=1;;b=2");
	cr_assert_eq(g_hash_table_size(parsed), 2);
	cr_assert_str_eq(nde_kv_get_str(parsed, "a"), "1");
	cr_assert_str_eq(nde_kv_get_str(parsed, "b"), "2");
	g_hash_table_unref(parsed);

	/* empty value is a value, not a malformed pair */
	parsed = nde_kv_parse("empty=");
	cr_assert_str_eq(nde_kv_get_str(parsed, "empty"), "");
	g_hash_table_unref(parsed);
}

/* ---------------- append / live_count ---------------- */

Test(nde_history, append_assigns_monotonic_ids) {
	cr_assert_eq(nde_history_live_count(), 0);
	gint64 id1 = append_op("stretch.mtf");
	gint64 id2 = append_op("filters.median");
	cr_assert_eq(id1, 1);
	cr_assert_eq(id2, 2);
	cr_assert_eq(nde_history_live_count(), 2);
}

Test(nde_history, append_without_single_image_refuses) {
	teardown();  /* drop com.uniq */
	setup();
	g_free(com.uniq);
	com.uniq = NULL;
	cr_assert_eq(append_op("stretch.mtf"), 0);
	com.uniq = g_new0(single, 1);  /* teardown expects it back */
}

Test(nde_history, undo_redo_move_live_count) {
	gint64 a = append_op("a.a"), b = append_op("b.b"), c = append_op("c.c");
	cr_assert_eq(nde_history_live_count(), 3);
	nde_history_on_undo(c);
	cr_assert_eq(nde_history_live_count(), 2);
	nde_history_on_undo(b);
	cr_assert_eq(nde_history_live_count(), 1);
	nde_history_on_redo(b);
	cr_assert_eq(nde_history_live_count(), 2);
	nde_history_on_redo(c);
	cr_assert_eq(nde_history_live_count(), 3);
	nde_history_on_undo(a);
	cr_assert_eq(nde_history_live_count(), 0);
	/* unknown / zero ids are ignored */
	nde_history_on_undo(0);
	nde_history_on_redo(999);
	cr_assert_eq(nde_history_live_count(), 0);
}

Test(nde_history, append_truncates_dead_tail) {
	append_op("a.a");
	gint64 b = append_op("b.b");
	append_op("c.c");
	nde_history_on_undo(b);   /* live: [a]; dead: b, c */
	cr_assert_eq(nde_history_live_count(), 1);
	gint64 d = append_op("d.d");
	cr_assert_eq(nde_history_live_count(), 2);
	/* ids stay monotonic — never reused even after truncation */
	cr_assert_eq(d, 4);
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 2);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(snap, 0))->op_id, "a.a");
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(snap, 1))->op_id, "d.d");
	g_ptr_array_unref(snap);
}

Test(nde_history, snapshot_is_live_prefix_deep_copy) {
	append_op("a.a");
	gint64 b = append_op("b.b");
	nde_history_on_undo(b);
	gint64 next_id = 0;
	GPtrArray *snap = nde_history_snapshot(&next_id);
	cr_assert_not_null(snap);
	cr_assert_eq(snap->len, 1, "snapshot must contain only the live prefix");
	cr_assert_eq(next_id, 3);
	nde_record *copy = g_ptr_array_index(snap, 0);
	cr_assert_str_eq(copy->op_id, "a.a");
	/* deep copy: mutating the snapshot must not touch the canonical log */
	g_free(copy->summary);
	copy->summary = g_strdup("mutated");
	g_ptr_array_unref(snap);
	snap = nde_history_snapshot(NULL);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(snap, 0))->summary, "a.a");
	g_ptr_array_unref(snap);

	/* empty history snapshots as NULL */
	nde_history_attach(NULL);
	cr_assert_null(nde_history_snapshot(NULL));
}

Test(nde_history, attach_replaces_and_stale_flag) {
	append_op("a.a");
	nde_history *h = g_new0(nde_history, 1);
	h->records = g_ptr_array_new_with_free_func((GDestroyNotify)nde_record_free);
	nde_record *rec = nde_record_new();
	rec->record_id = 7;
	rec->op_id = g_strdup("loaded.op");
	g_ptr_array_add(h->records, rec);
	h->live_count = 1;
	h->next_record_id = 8;
	nde_history_attach(h);

	cr_assert_eq(nde_history_live_count(), 1);
	gint64 next_id = 0;
	GPtrArray *snap = nde_history_snapshot(&next_id);
	cr_assert_eq(next_id, 8);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(snap, 0))->op_id, "loaded.op");
	g_ptr_array_unref(snap);

	cr_assert(!nde_history_is_stale());
	nde_history_set_stale(TRUE);
	cr_assert(nde_history_is_stale());
}

/* ---------------- capture helpers ---------------- */

Test(nde_history, capture_structural_fills_fields_and_monotonic) {
	gchar *params = g_strdup("w=10;h=20");
	gint64 id1 = nde_capture_structural("canvas.resize", NDE_SCOPE_CANVAS, -1,
	                                    params, "Resize canvas");
	cr_assert_eq(id1, 1, "first structural record gets id 1");
	gint64 id2 = nde_capture_structural("layer.add", NDE_SCOPE_DOCUMENT, 5,
	                                    g_strdup("name=Ha"), "Add layer");
	cr_assert_eq(id2, 2, "ids are monotonic");

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap);
	cr_assert_eq(snap->len, 2);
	nde_record *r0 = g_ptr_array_index(snap, 0);
	cr_assert_str_eq(r0->op_id, "canvas.resize");
	cr_assert_eq(r0->scope, NDE_SCOPE_CANVAS);
	cr_assert_eq(r0->target_item_id, -1);
	cr_assert_eq(r0->tier, NDE_TIER_A);
	cr_assert_str_eq(r0->params, "w=10;h=20");
	cr_assert_str_eq(r0->summary, "Resize canvas");
	cr_assert_not_null(r0->timestamp);
	cr_assert_not_null(r0->impl);
	nde_record *r1 = g_ptr_array_index(snap, 1);
	cr_assert_str_eq(r1->op_id, "layer.add");
	cr_assert_eq(r1->target_item_id, 5);
	g_ptr_array_unref(snap);
}

Test(nde_history, capture_opaque_is_tier_b_null_params) {
	gint64 id = nde_capture_opaque("python.set_pixeldata", NDE_SCOPE_LAYER, 3,
	                               "Python script pixel update", NULL);
	cr_assert_eq(id, 1);
	GPtrArray *snap = nde_history_snapshot(NULL);
	nde_record *r = g_ptr_array_index(snap, 0);
	cr_assert_eq(r->tier, NDE_TIER_B);
	cr_assert_null(r->params);
	cr_assert_str_eq(r->op_id, "python.set_pixeldata");
	cr_assert_eq(r->target_item_id, 3);
	g_ptr_array_unref(snap);
}

/* ---------------- thread smoke ---------------- */

#define SMOKE_APPENDS 500

static gpointer smoke_appender(gpointer data) {
	for (int i = 0; i < SMOKE_APPENDS; i++)
		append_op("smoke.op");
	return NULL;
}

Test(nde_history, append_vs_snapshot_thread_smoke) {
	GThread *worker = g_thread_new("nde-smoke", smoke_appender, NULL);
	for (int i = 0; i < 200; i++) {
		GPtrArray *snap = nde_history_snapshot(NULL);
		if (snap) {
			/* records must be complete and ordered whenever we look */
			gint64 prev = 0;
			for (guint j = 0; j < snap->len; j++) {
				nde_record *rec = g_ptr_array_index(snap, j);
				cr_assert_str_eq(rec->op_id, "smoke.op");
				cr_assert(rec->record_id > prev);
				prev = rec->record_id;
			}
			g_ptr_array_unref(snap);
		}
	}
	g_thread_join(worker);
	cr_assert_eq(nde_history_live_count(), SMOKE_APPENDS);
}

Test(nde_history, snapshot_all_includes_dead_tail) {
	append_op("a.a");
	gint64 b = append_op("b.b");
	nde_history_on_undo(b);      /* live: [a]; dead: [b] */

	guint live = 0;
	GPtrArray *all = nde_history_snapshot_all(&live);
	cr_assert_not_null(all);
	cr_assert_eq(all->len, 2, "snapshot_all must include the dead tail");
	cr_assert_eq(live, 1);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(all, 1))->op_id, "b.b");
	g_ptr_array_unref(all);

	/* the persistence snapshot stays live-only */
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 1);
	g_ptr_array_unref(snap);

	/* empty history → NULL, live 0 */
	nde_history_attach(NULL);
	live = 99;
	cr_assert_null(nde_history_snapshot_all(&live));
	cr_assert_eq(live, 0);
}

Test(nde_history, capture_from_descriptor) {
	/* Tier A: a descriptor with a serializer (crop is POD) */
	struct crop_args ca = { 0 };
	ca.area.x = 1; ca.area.y = 2; ca.area.w = 30; ca.area.h = 40;
	gint64 id = nde_capture_from_descriptor(&op_desc_crop, &ca, "crop it", NULL);
	cr_assert(id > 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	nde_record *rec = g_ptr_array_index(snap, snap->len - 1);
	cr_assert_str_eq(rec->op_id, op_desc_crop.id);
	cr_assert_eq(rec->tier, NDE_TIER_A);
	cr_assert_not_null(rec->params);
	cr_assert(strstr(rec->params, "w=30") != NULL, "params blob missing field: %s", rec->params);
	cr_assert_eq(rec->scope, NDE_SCOPE_CANVAS, "crop carries OP_GEOMETRY_CHANGING");
	cr_assert_str_eq(rec->summary, "crop it");
	cr_assert_not_null(rec->timestamp);
	cr_assert_not_null(rec->impl);
	g_ptr_array_unref(snap);

	/* Tier B: first descriptor without a serializer */
	size_t n = 0;
	const op_descriptor *const *all = op_descriptor_all(&n);
	const op_descriptor *opaque = NULL;
	for (size_t i = 0; i < n && !opaque; i++)
		if (!all[i]->serialize)
			opaque = all[i];
	cr_assert_not_null(opaque, "expected at least one serializer-less descriptor");
	id = nde_capture_from_descriptor(opaque, NULL, "opaque op", NULL);
	cr_assert(id > 0);
	snap = nde_history_snapshot(NULL);
	rec = g_ptr_array_index(snap, snap->len - 1);
	cr_assert_eq(rec->tier, NDE_TIER_B);
	cr_assert_null(rec->params);
	g_ptr_array_unref(snap);
}

/* ---------------- P3.A: amend / delete ---------------- */

static gint64 append_tier_a(const char *op_id, const char *params) {
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup(op_id);
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->params = g_strdup(params);
	rec->summary = g_strdup(op_id);
	rec->timestamp = g_strdup("2026-01-01T00:00:00Z");
	rec->impl = g_strdup("test");
	return nde_history_append(rec);
}

Test(nde_history, amend_updates_params_and_truncates_dead_tail) {
	gint64 a = append_tier_a("geometry.mirrorx", "x_axis=1");
	gint64 c = append_tier_a("geometry.binning", "factor=2;mean=1");
	nde_history_on_undo(c);              /* dead tail: [c] */
	cr_assert_eq(nde_history_live_count(), 1);

	gchar *err = NULL;
	cr_assert(nde_history_amend(a, "x_axis=0", &err), "amend failed: %s", err ? err : "?");
	cr_assert_null(err);

	guint live = 0;
	GPtrArray *all = nde_history_snapshot_all(&live);
	cr_assert_eq(all->len, 1, "the dead tail must be truncated by an amend");
	cr_assert_eq(live, 1);
	nde_record *rec = g_ptr_array_index(all, 0);
	cr_assert_eq(rec->record_id, a, "record id must be preserved");
	cr_assert_str_eq(rec->params, "x_axis=0");
	cr_assert_str_neq(rec->timestamp, "2026-01-01T00:00:00Z", "timestamp must refresh");
	cr_assert_str_neq(rec->impl, "test", "impl must refresh");
	g_ptr_array_unref(all);
}

Test(nde_history, amend_rejections) {
	gchar *err = NULL;

	/* no history at all */
	cr_assert(!nde_history_amend(1, "x_axis=0", &err));
	cr_assert_not_null(err);
	g_clear_pointer(&err, g_free);

	gint64 a = append_tier_a("geometry.mirrorx", "x_axis=1");
	gint64 b = append_op("python.set_pixeldata");          /* Tier B */
	gint64 c = append_tier_a("future.unknown_op", "k=1");  /* Tier A, unknown op */
	gint64 d = append_tier_a("geometry.crop", "x=0;y=0;w=2;h=2");

	/* unknown id */
	cr_assert(!nde_history_amend(999, "x_axis=0", &err));
	g_clear_pointer(&err, g_free);

	/* Tier B */
	cr_assert(!nde_history_amend(b, "x=1", &err));
	cr_assert(strstr(err, "opaque") != NULL);
	g_clear_pointer(&err, g_free);

	/* unknown op */
	cr_assert(!nde_history_amend(c, "k=2", &err));
	g_clear_pointer(&err, g_free);

	/* unparsable params (crop requires its geometry keys) */
	cr_assert(!nde_history_amend(d, "nonsense=only", &err));
	cr_assert_not_null(err);
	g_clear_pointer(&err, g_free);

	/* dead record */
	nde_history_on_undo(d);
	cr_assert(!nde_history_amend(d, "x=0;y=0;w=1;h=1", &err));
	cr_assert(strstr(err, "undone") != NULL);
	g_clear_pointer(&err, g_free);

	/* nothing was modified by the failures */
	guint live = 0;
	GPtrArray *all = nde_history_snapshot_all(&live);
	cr_assert_eq(all->len, 4);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(all, 0))->params, "x_axis=1");
	g_ptr_array_unref(all);
	(void)a;
}

Test(nde_history, delete_removes_live_records) {
	gint64 a = append_tier_a("geometry.mirrorx", "x_axis=1");
	gint64 b = append_op("python.set_pixeldata");           /* Tier B */
	gint64 c = append_tier_a("geometry.binning", "factor=2;mean=0");
	gchar *err = NULL;

	/* Tier B IS deletable at the log level (deleting an opaque step is
	 * well-defined — policy checks live in nde_delete_execute), but this
	 * test keeps b as a survivor, so only assert the unknown-id refusal. */
	cr_assert(!nde_history_delete(999, &err));
	g_clear_pointer(&err, g_free);

	/* delete the first record; ids of the others are untouched */
	cr_assert(nde_history_delete(a, &err), "delete failed: %s", err ? err : "?");
	cr_assert_eq(nde_history_live_count(), 2);
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(((nde_record *)g_ptr_array_index(snap, 0))->record_id, b);
	cr_assert_eq(((nde_record *)g_ptr_array_index(snap, 1))->record_id, c);
	g_ptr_array_unref(snap);

	/* delete with a dead tail: undo c, delete b?  b is Tier B — use c:
	 * append d, undo d, then delete c must drop the dead d too */
	gint64 d = append_tier_a("geometry.mirrorx", "x_axis=1");
	nde_history_on_undo(d);
	cr_assert(nde_history_delete(c, &err), "delete failed: %s", err ? err : "?");
	guint live = 0;
	GPtrArray *all = nde_history_snapshot_all(&live);
	cr_assert_eq(all->len, 1, "dead tail must be truncated by a delete");
	cr_assert_eq(live, 1);
	cr_assert_eq(((nde_record *)g_ptr_array_index(all, 0))->record_id, b);
	g_ptr_array_unref(all);
}

/* ---------------- graph step 2: edit-at insertion point ---------------- */

/* Append a Tier-A record with an explicit scope/target, so the qualification
 * rules can be exercised without a document. */
static gint64 append_scoped(const char *op_id, gint scope, gint target) {
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup(op_id);
	rec->op_version = 1;
	rec->tier = NDE_TIER_A;
	rec->scope = scope;
	rec->target_item_id = target;
	rec->params = g_strdup("k=1");
	rec->summary = g_strdup(op_id);
	return nde_history_append(rec);
}

/* Ids in log order. */
static GArray *order_ids(void) {
	GArray *out = g_array_new(FALSE, FALSE, sizeof(gint64));
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; snap && i < snap->len; i++) {
		gint64 id = ((nde_record *)g_ptr_array_index(snap, i))->record_id;
		g_array_append_val(out, id);
	}
	if (snap)
		g_ptr_array_unref(snap);
	return out;
}

#define ASSERT_ORDER(...) do { \
	gint64 want[] = { __VA_ARGS__ }; \
	GArray *got = order_ids(); \
	cr_assert_eq(got->len, G_N_ELEMENTS(want), "expected %zu records, got %u", \
	             (size_t)G_N_ELEMENTS(want), got->len); \
	for (guint _i = 0; _i < got->len; _i++) \
		cr_assert_eq(g_array_index(got, gint64, _i), want[_i], \
		             "position %u: expected record %" G_GINT64_FORMAT ", got %" G_GINT64_FORMAT, \
		             _i, want[_i], g_array_index(got, gint64, _i)); \
	g_array_unref(got); \
} while (0)

Test(nde_history, insert_point_places_records_before_the_anchor) {
	gint64 a = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 b = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 c = append_scoped("filters.unsharp", NDE_SCOPE_LAYER, -1);

	gchar *err = NULL;
	cr_assert(nde_history_insert_point_set(c, -1, &err), "arm failed: %s", err ? err : "?");
	cr_assert_eq(nde_history_insert_point(), c);

	gint64 x = append_scoped("stretch.asinh", NDE_SCOPE_LAYER, -1);
	gint64 y = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	/* successive inserts accumulate in order, all still before the anchor */
	ASSERT_ORDER(a, b, x, y, c);
	cr_assert_eq(nde_history_live_count(), 5);

	GArray *ins = nde_history_insert_point_clear();
	cr_assert_eq(ins->len, 2);
	cr_assert_eq(g_array_index(ins, gint64, 0), x);
	cr_assert_eq(g_array_index(ins, gint64, 1), y);
	g_array_unref(ins);
	cr_assert_eq(nde_history_insert_point(), 0);

	/* disarmed: back to appending */
	gint64 z = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	ASSERT_ORDER(a, b, x, y, c, z);
}

Test(nde_history, insert_point_arming_requires_a_live_anchor) {
	gchar *err = NULL;
	gint64 a = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 b = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);

	cr_assert(!nde_history_insert_point_set(999, -1, &err));
	cr_assert_not_null(err);
	g_clear_pointer(&err, g_free);
	cr_assert_eq(nde_history_insert_point(), 0);

	/* an undone record is not a valid anchor */
	nde_history_on_undo(b);
	cr_assert(!nde_history_insert_point_set(b, -1, &err));
	cr_assert(strstr(err, "undone") != NULL);
	g_clear_pointer(&err, g_free);
	(void)a;
}

Test(nde_history, insert_point_truncates_the_dead_tail_when_armed) {
	gint64 a = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 b = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 c = append_scoped("filters.unsharp", NDE_SCOPE_LAYER, -1);
	nde_history_on_undo(c);                     /* dead tail: [c] */
	guint live = 0;
	GPtrArray *all = nde_history_snapshot_all(&live);
	cr_assert_eq(all->len, 3);
	cr_assert_eq(live, 2);
	g_ptr_array_unref(all);

	gchar *err = NULL;
	cr_assert(nde_history_insert_point_set(b, -1, &err), "arm failed: %s", err ? err : "?");
	all = nde_history_snapshot_all(&live);
	cr_assert_eq(all->len, 2, "arming must drop the dead tail");
	cr_assert_eq(live, 2);
	g_ptr_array_unref(all);
	(void)a;
}

Test(nde_history, insert_point_qualification_rules) {
	gint64 anchor = append_scoped("filters.unsharp", NDE_SCOPE_LAYER, 7);
	gchar *err = NULL;
	cr_assert(nde_history_insert_point_set(anchor, 7, &err), "arm failed: %s", err ? err : "?");

	/* LAYER scope on the insertion's item: inserted */
	gint64 mine = append_scoped("filters.gauss", NDE_SCOPE_LAYER, 7);
	/* LAYER scope on ANOTHER layer: appended, and harmless */
	gint64 other = append_scoped("filters.gauss", NDE_SCOPE_LAYER, 9);
	/* non-destructive structural: appended, and harmless */
	gint64 add = append_scoped("layer.add", NDE_SCOPE_DOCUMENT, 9);
	ASSERT_ORDER(mine, anchor, other, add);
	cr_assert(!nde_history_insert_point_disturbed());

	GArray *ins = nde_history_insert_point_clear();
	cr_assert_eq(ins->len, 1);
	cr_assert_eq(g_array_index(ins, gint64, 0), mine);
	g_array_unref(ins);
}

Test(nde_history, insert_point_canvas_qualifies_only_on_a_plain_image) {
	gint64 anchor = append_scoped("filters.unsharp", NDE_SCOPE_LAYER, -1);
	gchar *err = NULL;
	cr_assert(nde_history_insert_point_set(anchor, -1, &err), "arm failed: %s", err ? err : "?");
	/* plain image (item -1): the image IS the layer, so geometry belongs to
	 * the lineage and is inserted */
	gint64 crop = append_scoped("geometry.crop", NDE_SCOPE_CANVAS, -1);
	ASSERT_ORDER(crop, anchor);
	cr_assert(!nde_history_insert_point_disturbed());
	g_array_unref(nde_history_insert_point_clear());
}

Test(nde_history, insert_point_disturbing_records_raise_the_flag) {
	gint64 anchor = append_scoped("filters.unsharp", NDE_SCOPE_LAYER, 7);
	gchar *err = NULL;
	cr_assert(nde_history_insert_point_set(anchor, 7, &err), "arm failed: %s", err ? err : "?");
	cr_assert(!nde_history_insert_point_disturbed());

	/* canvas geometry on a LAYERED document is not this layer's lineage */
	append_scoped("geometry.crop", NDE_SCOPE_CANVAS, 7);
	cr_assert(nde_history_insert_point_disturbed());
	g_array_unref(nde_history_insert_point_clear());

	cr_assert(nde_history_insert_point_set(anchor, 7, &err));
	append_scoped("document.flatten", NDE_SCOPE_DOCUMENT, 7);
	cr_assert(nde_history_insert_point_disturbed());
	g_array_unref(nde_history_insert_point_clear());

	cr_assert(nde_history_insert_point_set(anchor, 7, &err));
	append_scoped("layer.merge_down", NDE_SCOPE_DOCUMENT, 7);
	cr_assert(nde_history_insert_point_disturbed());
	g_array_unref(nde_history_insert_point_clear());
}

Test(nde_history, insert_point_undo_lifts_the_record_out_and_redo_puts_it_back) {
	gint64 a = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 c = append_scoped("filters.unsharp", NDE_SCOPE_LAYER, -1);
	gchar *err = NULL;
	cr_assert(nde_history_insert_point_set(c, -1, &err), "arm failed: %s", err ? err : "?");

	gint64 x = append_scoped("stretch.asinh", NDE_SCOPE_LAYER, -1);
	gint64 y = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	ASSERT_ORDER(a, x, y, c);

	/* Undo must NOT declare the anchor and its successors dead: it lifts the
	 * inserted record out of the log entirely. */
	nde_history_on_undo(y);
	ASSERT_ORDER(a, x, c);
	cr_assert_eq(nde_history_live_count(), 3);
	guint live = 0;
	GPtrArray *all = nde_history_snapshot_all(&live);
	cr_assert_eq(all->len, 3, "an undone insert leaves no dead tail behind");
	g_ptr_array_unref(all);

	nde_history_on_redo(y);
	ASSERT_ORDER(a, x, y, c);

	/* undo both, then finish: nothing was inserted */
	nde_history_on_undo(y);
	nde_history_on_undo(x);
	ASSERT_ORDER(a, c);
	GArray *ins = nde_history_insert_point_clear();
	cr_assert_eq(ins->len, 0);
	g_array_unref(ins);
}

Test(nde_history, insert_point_capture_after_an_undo_supersedes_it) {
	gint64 c = append_scoped("filters.unsharp", NDE_SCOPE_LAYER, -1);
	gchar *err = NULL;
	cr_assert(nde_history_insert_point_set(c, -1, &err), "arm failed: %s", err ? err : "?");
	gint64 x = append_scoped("stretch.asinh", NDE_SCOPE_LAYER, -1);
	nde_history_on_undo(x);
	gint64 z = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	ASSERT_ORDER(z, c);
	/* x is gone for good: redoing it now must not resurrect it */
	nde_history_on_redo(x);
	ASSERT_ORDER(z, c);
	GArray *ins = nde_history_insert_point_clear();
	cr_assert_eq(ins->len, 1);
	cr_assert_eq(g_array_index(ins, gint64, 0), z);
	g_array_unref(ins);
}

Test(nde_history, drop_records_removes_only_the_named_live_records) {
	gint64 a = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 b = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 c = append_scoped("filters.unsharp", NDE_SCOPE_LAYER, -1);

	GArray *ids = g_array_new(FALSE, FALSE, sizeof(gint64));
	g_array_append_val(ids, b);
	gint64 missing = 4242;
	g_array_append_val(ids, missing);
	nde_history_drop_records(ids);
	g_array_unref(ids);

	ASSERT_ORDER(a, c);
	cr_assert_eq(nde_history_live_count(), 2);
}

Test(nde_history, truncate_dead_drops_the_undone_tail) {
	gint64 a = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	gint64 b = append_scoped("filters.gauss", NDE_SCOPE_LAYER, -1);
	nde_history_on_undo(b);
	guint live = 0;
	GPtrArray *all = nde_history_snapshot_all(&live);
	cr_assert_eq(all->len, 2);
	g_ptr_array_unref(all);

	nde_history_truncate_dead();
	all = nde_history_snapshot_all(&live);
	cr_assert_eq(all->len, 1);
	cr_assert_eq(live, 1);
	cr_assert_eq(((nde_record *)g_ptr_array_index(all, 0))->record_id, a);
	g_ptr_array_unref(all);
}

/* ---------------- graph step 4: input pins ---------------- */

Test(nde_history, input_pins_add_lookup_and_replace) {
	nde_record *rec = nde_record_new();
	cr_assert_null(rec->inputs, "no pins until one is added");
	cr_assert_null(nde_record_input(rec, "mask"));

	nde_record_add_input(rec, "mask", 3, 7);
	nde_record_add_input(rec, "overlay", 5, 0);
	cr_assert_eq(rec->inputs->len, 2);

	const nde_input_pin *m = nde_record_input(rec, "mask");
	cr_assert_not_null(m);
	cr_assert_eq(m->src_item_id, 3);
	cr_assert_eq(m->src_record_id, 7);

	/* a record has at most one input per role: re-adding rewires it */
	nde_record_add_input(rec, "mask", 9, 11);
	cr_assert_eq(rec->inputs->len, 2, "re-adding a role must not append");
	m = nde_record_input(rec, "mask");
	cr_assert_eq(m->src_item_id, 9);
	cr_assert_eq(m->src_record_id, 11);
	/* the other role is untouched */
	cr_assert_eq(nde_record_input(rec, "overlay")->src_item_id, 5);
	nde_record_free(rec);
}

Test(nde_history, input_pins_survive_a_record_copy) {
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("filters.gauss");
	nde_record_add_input(rec, "mask", 3, 7);
	nde_record *copy = nde_record_copy(rec);
	/* deep: freeing the original must not disturb the copy */
	nde_record_free(rec);
	cr_assert_not_null(copy->inputs);
	cr_assert_eq(copy->inputs->len, 1);
	const nde_input_pin *m = nde_record_input(copy, "mask");
	cr_assert_str_eq(m->role, "mask");
	cr_assert_eq(m->src_item_id, 3);
	cr_assert_eq(m->src_record_id, 7);
	nde_record_free(copy);
}

Test(nde_history, pins_codec_roundtrip) {
	nde_record *rec = nde_record_new();
	nde_record_add_input(rec, "mask", 3, 7);
	nde_record_add_input(rec, "base", -2, 0);
	gchar *blob = nde_pins_serialize(rec->inputs);
	cr_assert_not_null(blob);

	GPtrArray *back = nde_pins_parse(blob);
	cr_assert_not_null(back);
	cr_assert_eq(back->len, 2);
	const nde_input_pin *a = g_ptr_array_index(back, 0);
	cr_assert_str_eq(a->role, "mask");
	cr_assert_eq(a->src_item_id, 3);
	cr_assert_eq(a->src_record_id, 7);
	const nde_input_pin *b = g_ptr_array_index(back, 1);
	cr_assert_str_eq(b->role, "base");
	cr_assert_eq(b->src_item_id, -2);
	cr_assert_eq(b->src_record_id, 0);
	g_ptr_array_unref(back);
	g_free(blob);
	nde_record_free(rec);
}

Test(nde_history, pins_codec_edge_cases) {
	cr_assert_null(nde_pins_serialize(NULL), "nothing to persist");
	GPtrArray *empty = g_ptr_array_new_with_free_func((GDestroyNotify)nde_input_pin_free);
	cr_assert_null(nde_pins_serialize(empty));
	g_ptr_array_unref(empty);

	cr_assert_null(nde_pins_parse(NULL));
	cr_assert_null(nde_pins_parse(""));
	cr_assert_null(nde_pins_parse("n=0"));
	/* an incomplete pin points nowhere and is skipped rather than invented */
	cr_assert_null(nde_pins_parse("n=1;item0=3;rec0=7"), "no role");
	cr_assert_null(nde_pins_parse("n=1;role0=mask;rec0=7"), "no source item");
	/* an absent rec means "that item's baseline", which IS complete */
	GPtrArray *base = nde_pins_parse("n=1;role0=mask;item0=4");
	cr_assert_not_null(base);
	cr_assert_eq(((nde_input_pin *)g_ptr_array_index(base, 0))->src_record_id, 0);
	g_ptr_array_unref(base);
}

Test(nde_history, pins_codec_escapes_awkward_roles) {
	nde_record *rec = nde_record_new();
	/* roles are programmer constants, but the codec must not be the weak
	 * point if one ever carries a separator */
	nde_record_add_input(rec, "od;d=role\\x", 1, 2);
	gchar *blob = nde_pins_serialize(rec->inputs);
	GPtrArray *back = nde_pins_parse(blob);
	cr_assert_not_null(back);
	cr_assert_eq(back->len, 1);
	cr_assert_str_eq(((nde_input_pin *)g_ptr_array_index(back, 0))->role, "od;d=role\\x");
	g_ptr_array_unref(back);
	g_free(blob);
	nde_record_free(rec);
}

Test(nde_history, pins_to_string_names_the_baseline) {
	nde_record *rec = nde_record_new();
	nde_record_add_input(rec, "mask", 3, 7);
	nde_record_add_input(rec, "base", 1, 0);
	gchar *s = nde_pins_to_string(rec->inputs);
	cr_assert_str_eq(s, "mask@3:7, base@1:baseline");
	g_free(s);
	nde_record_free(rec);
	cr_assert_null(nde_pins_to_string(NULL));
}

Test(nde_history, last_record_for_item_tracks_the_end_of_a_chain) {
	cr_assert_eq(nde_history_last_record_for_item(5), 0, "no history yet");
	gint64 a = append_scoped("mask.from_stars", NDE_SCOPE_LAYER, 5);
	cr_assert_eq(nde_history_last_record_for_item(5), a);
	append_scoped("filters.gauss", NDE_SCOPE_LAYER, 9);   /* another item */
	cr_assert_eq(nde_history_last_record_for_item(5), a, "other items do not count");
	gint64 c = append_scoped("mask.blur", NDE_SCOPE_LAYER, 5);
	cr_assert_eq(nde_history_last_record_for_item(5), c);

	/* an undone record is not the item's current state */
	nde_history_on_undo(c);
	cr_assert_eq(nde_history_last_record_for_item(5), a);
	cr_assert_eq(nde_history_last_record_for_item(404), 0);
}
