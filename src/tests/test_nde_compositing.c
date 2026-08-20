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
 */

/*
 * test_nde_compositing — live compositing-state records (graph step 1,
 * nde-graph-design-note.md §6.1).  layer.set_opacity / set_blend /
 * set_visible are inputs to the COMPOSITOR, not pixel operations: editing one
 * re-folds the log onto the layer and re-renders, with no chain build and no
 * replay.  These tests cover the fold (last-value-wins, reset points,
 * properties the log does not describe), the params validator, and the
 * amend/delete policy predicates that used to refuse these records outright.
 */

#include <criterion/criterion.h>
#include <glib.h>
#include <string.h>
#include "flis_test_helpers.h"
#include "core/siril.h"
#include "core/nde/nde_history.h"
#include "core/nde/nde_compositing.h"
#include "core/nde/nde_replay.h"
#include "io/image_format_fits.h"

cominfo com;
fits *gfit;

static flis_layer_t *layer;

static void comp_test_init(void) {
	flis_test_init_com();
	fits *f = flis_test_make_mono_fits(8, 8, 0.25f);
	cr_assert_not_null(f);
	layer = flis_test_add_layer(f, "base");
	cr_assert_not_null(layer);
	gfit = layer->fit;
}

static void comp_test_fini(void) {
	layer = NULL;
	flis_test_cleanup_com();
}

TestSuite(nde_compositing, .init = comp_test_init, .fini = comp_test_fini);

/* Append a compositing record for the test layer and return its id. */
static gint64 capture_opacity(float v) {
	GString *kv = nde_kv_start();
	nde_kv_add_float(kv, "opacity", v);
	return nde_capture_structural("layer.set_opacity", NDE_SCOPE_LAYER,
	                              layer->item_id, nde_kv_end(kv), "Set opacity");
}

static gint64 capture_blend(gint64 v) {
	GString *kv = nde_kv_start();
	nde_kv_add_int(kv, "blend", v);
	return nde_capture_structural("layer.set_blend", NDE_SCOPE_LAYER,
	                              layer->item_id, nde_kv_end(kv), "Set blend mode");
}

static gint64 capture_visible(gboolean v) {
	GString *kv = nde_kv_start();
	nde_kv_add_bool(kv, "visible", v);
	return nde_capture_structural("layer.set_visible", NDE_SCOPE_LAYER,
	                              layer->item_id, nde_kv_end(kv), "Set visibility");
}

/* ---- op identification -------------------------------------------------- */

Test(nde_compositing, identifies_the_three_ops) {
	cr_assert(nde_compositing_is_op("layer.set_opacity"));
	cr_assert(nde_compositing_is_op("layer.set_blend"));
	cr_assert(nde_compositing_is_op("layer.set_visible"));
	cr_assert(!nde_compositing_is_op("geometry.mirrorx"));
	cr_assert(!nde_compositing_is_op("document.flatten"));
	cr_assert(!nde_compositing_is_op(NULL));
}

/* ---- validation --------------------------------------------------------- */

Test(nde_compositing, validate_accepts_good_params) {
	gchar *err = NULL;
	cr_assert(nde_compositing_validate("layer.set_opacity", "opacity=0.5", &err), "%s", err);
	cr_assert_null(err);
	cr_assert(nde_compositing_validate("layer.set_blend", "blend=2", &err), "%s", err);
	cr_assert(nde_compositing_validate("layer.set_visible", "visible=0", &err), "%s", err);
}

Test(nde_compositing, validate_rejects_out_of_range_opacity) {
	gchar *err = NULL;
	cr_assert(!nde_compositing_validate("layer.set_opacity", "opacity=1.5", &err));
	cr_assert_not_null(err);
	g_free(err);
	err = NULL;
	cr_assert(!nde_compositing_validate("layer.set_opacity", "opacity=-0.1", &err));
	g_free(err);
	err = NULL;
	cr_assert(!nde_compositing_validate("layer.set_opacity", "blend=2", &err),
	          "a blob without the opacity key must be rejected");
	g_free(err);
}

Test(nde_compositing, validate_rejects_bad_blend_mode) {
	gchar *err = NULL;
	cr_assert(!nde_compositing_validate("layer.set_blend", "blend=3", &err),
	          "3 is not a power-of-two blend flag");
	g_free(err);
	err = NULL;
	/* PASS_THROUGH is a GROUP mode — a layer must not be set to it. */
	cr_assert(!nde_compositing_validate("layer.set_blend", "blend=131072", &err));
	g_free(err);
}

/* ---- the fold ----------------------------------------------------------- */

Test(nde_compositing, fold_takes_the_last_value) {
	capture_opacity(0.5f);
	capture_opacity(0.25f);
	layer->opacity = 1.0f;   /* clobber, as an undo restore would */

	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 0.25f, 1e-6,
	                   "the fold must land on the LAST recorded opacity");
}

Test(nde_compositing, fold_returns_undescribed_properties_to_defaults) {
	/* The log is authoritative: a property with no record is at its default,
	 * whatever the live layer happens to hold.  Anything else would let the
	 * layer sit in a state the history does not describe. */
	capture_opacity(0.5f);
	layer->blend_mode = FLIS_BLEND_SCREEN;
	layer->visible = FALSE;

	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 0.5f, 1e-6);
	cr_assert_eq(layer->blend_mode, FLIS_BLEND_NORMAL,
	             "no blend record exists — the property folds to its default");
	cr_assert_eq(layer->visible, TRUE,
	             "no visibility record exists — the property folds to its default");
}

Test(nde_compositing, fold_covers_all_three_properties) {
	capture_opacity(0.75f);
	capture_blend(FLIS_BLEND_LUMINOSITY);
	capture_visible(FALSE);
	layer->opacity = 0.f;
	layer->blend_mode = FLIS_BLEND_NORMAL;
	layer->visible = TRUE;

	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 0.75f, 1e-6);
	cr_assert_eq(layer->blend_mode, FLIS_BLEND_LUMINOSITY);
	cr_assert_eq(layer->visible, FALSE);
}

Test(nde_compositing, flatten_is_a_reset_point) {
	capture_opacity(0.2f);
	capture_blend(FLIS_BLEND_SCREEN);
	/* flatten resets the surviving base to the compositing defaults */
	nde_capture_structural("document.flatten", NDE_SCOPE_DOCUMENT,
	                       layer->item_id, g_strdup("n_layers=3"), "Flatten image");
	layer->opacity = 0.f;
	layer->blend_mode = FLIS_BLEND_DARKEN;

	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 1.0f, 1e-6,
	                   "flatten resets opacity, superseding the earlier record");
	cr_assert_eq(layer->blend_mode, FLIS_BLEND_NORMAL);
}

Test(nde_compositing, records_after_a_reset_point_still_apply) {
	nde_capture_structural("document.flatten", NDE_SCOPE_DOCUMENT,
	                       layer->item_id, g_strdup("n_layers=2"), "Flatten image");
	capture_opacity(0.4f);

	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 0.4f, 1e-6);
}

Test(nde_compositing, merge_down_is_not_a_reset_point) {
	/* merge_down preserves the surviving bottom layer's blend and opacity —
	 * unlike flatten, it must NOT restore the defaults. */
	capture_opacity(0.3f);
	nde_capture_structural("layer.merge_down", NDE_SCOPE_DOCUMENT,
	                       layer->item_id, g_strdup("top_item=9"), "Merge down");
	layer->opacity = 1.0f;

	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 0.3f, 1e-6,
	                   "merge down must not reset the bottom layer's opacity");
}

Test(nde_compositing, other_items_records_are_not_adopted) {
	GString *kv = nde_kv_start();
	nde_kv_add_float(kv, "opacity", 0.1f);
	nde_capture_structural("layer.set_opacity", NDE_SCOPE_LAYER,
	                       layer->item_id + 50, nde_kv_end(kv), "Set opacity");

	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 1.0f, 1e-6,
	                   "another layer's opacity record must not be adopted; with "
	                   "no record of its own this layer folds to the default");
}

Test(nde_compositing, plain_image_is_a_no_op) {
	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(-1, &err));
	cr_assert_null(err);
}

Test(nde_compositing, missing_layer_reports_an_error) {
	gchar *err = NULL;
	cr_assert(!nde_compositing_recompute(layer->item_id + 500, &err));
	cr_assert_not_null(err);
	g_free(err);
}

/* ---- policy predicates -------------------------------------------------- */

Test(nde_compositing, records_are_amendable_and_deletable) {
	capture_opacity(0.5f);
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap);
	cr_assert_eq(snap->len, 1);
	const nde_record *rec = g_ptr_array_index(snap, 0);
	cr_assert(nde_record_amendable(rec),
	          "compositing records have editable params (validated by range, "
	          "not by an op descriptor)");
	cr_assert(nde_record_deletable(rec),
	          "deleting a compositing record is well-defined — the fold simply "
	          "falls back to the previous value");
	g_ptr_array_unref(snap);
}

/* ---- log edits re-fold -------------------------------------------------- */

Test(nde_compositing, amend_updates_the_layer) {
	gint64 id = capture_opacity(0.5f);
	gchar *err = NULL;
	cr_assert(nde_history_amend(id, "opacity=0.125", &err), "%s", err ? err : "");
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 0.125f, 1e-6,
	                   "amending the record must change the live layer opacity");
}

Test(nde_compositing, amend_rejects_invalid_params) {
	gint64 id = capture_opacity(0.5f);
	gchar *err = NULL;
	cr_assert(!nde_history_amend(id, "opacity=42", &err),
	          "the log commit must refuse an out-of-range opacity");
	cr_assert_not_null(err);
	g_free(err);
}

Test(nde_compositing, delete_falls_back_to_the_previous_value) {
	capture_opacity(0.5f);
	gint64 second = capture_opacity(0.25f);

	gchar *err = NULL;
	cr_assert(nde_history_delete(second, &err), "%s", err ? err : "");
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 0.5f, 1e-6,
	                   "with the last record gone the fold lands on its predecessor");
}

Test(nde_compositing, delete_of_the_only_record_restores_the_default) {
	gint64 id = capture_opacity(0.25f);
	gchar *err = NULL;
	cr_assert(nde_history_delete(id, &err), "%s", err ? err : "");
	layer->opacity = 0.25f;   /* the delete has not re-folded yet */
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 1.0f, 1e-6,
	                   "deleting the step that set the opacity must undo it, "
	                   "returning the layer to the default");
}

Test(nde_compositing, other_items_records_do_not_reset_this_one) {
	/* The fold is per item: a record for another layer must not drag this
	 * layer's own recorded opacity back to the default. */
	capture_opacity(0.6f);
	GString *kv = nde_kv_start();
	nde_kv_add_float(kv, "opacity", 0.1f);
	nde_capture_structural("layer.set_opacity", NDE_SCOPE_LAYER,
	                       layer->item_id + 50, nde_kv_end(kv), "Set opacity");

	gchar *err = NULL;
	cr_assert(nde_compositing_recompute(layer->item_id, &err), "%s", err ? err : "");
	cr_assert_float_eq(layer->opacity, 0.6f, 1e-6);
}
