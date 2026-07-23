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
 * test_nde_persist — FLIS_HIST persistence round-trips (nde sketch §14/§18):
 *   - save with history → load → record field equality, ids + FLISHSEQ kept
 *   - unknown op_id preservation
 *   - only the live prefix persists (dead tail dropped)
 *   - no history → no FLIS_HIST, empty on reload
 *   - PIXHASH: match → not stale; tampered layer HDU → stale;
 *     compressed save (hashes skipped) → not stale
 *   - FLISNDE=T with the table deleted → warn, empty history
 */

#include <criterion/criterion.h>
#include <unistd.h>
#include <fitsio.h>
#include "flis_test_helpers.h"
#include "core/nde_history.h"

cominfo com;
fits *gfit;

static char *tmpdir = NULL;
static char *tmppath = NULL;

static void setup(void) {
	flis_test_init_com();
	tmpdir = g_dir_make_tmp("nde-persist-XXXXXX", NULL);
	tmppath = g_build_filename(tmpdir, "test.flis", NULL);
}

static void teardown(void) {
	nde_history_attach(NULL);
	if (tmppath) { g_unlink(tmppath); g_free(tmppath); tmppath = NULL; }
	if (tmpdir)  { g_rmdir(tmpdir);   g_free(tmpdir);   tmpdir  = NULL; }
	flis_test_cleanup_com();
}

TestSuite(nde_persist, .init = setup, .fini = teardown);

/* Append a fully-populated record; returns its id. */
static gint64 append_full(const char *op_id, int version, const char *params,
                          const char *summary, gint tier, gint scope,
                          gint target, gboolean mask_active) {
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup(op_id);
	rec->op_version = version;
	rec->params = g_strdup(params);
	rec->summary = g_strdup(summary);
	rec->tier = tier;
	rec->scope = scope;
	rec->target_item_id = target;
	rec->mask_active = mask_active;
	rec->timestamp = nde_iso8601_now();
	rec->impl = nde_impl_string();
	return nde_history_append(rec);
}

Test(nde_persist, history_roundtrip_field_equality) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;
	gint64 id1 = append_full("stretch.mtf", 1, "linked=1;lo=0.001;mid=0.25;hi=1",
	                         "Midtone Transfer Function", NDE_TIER_A,
	                         NDE_SCOPE_LAYER, lid, TRUE);
	/* unknown op id + escaped params must survive the trip verbatim */
	gint64 id2 = append_full("future.shiny_op", 99, "text=semi\\;colon;w=\\=x\\\\y",
	                         "an op from the future", NDE_TIER_B,
	                         NDE_SCOPE_CANVAS, -1, FALSE);
	cr_assert_eq(id1, 1);
	cr_assert_eq(id2, 2);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert_eq(nde_history_live_count(), 2);
	cr_assert(!nde_history_is_stale(), "clean reload must not be stale");
	gint64 next_id = 0;
	GPtrArray *snap = nde_history_snapshot(&next_id);
	cr_assert_eq(next_id, 3, "FLISHSEQ must preserve the id sequence");

	nde_record *r1 = g_ptr_array_index(snap, 0);
	cr_assert_eq(r1->record_id, 1);
	cr_assert_str_eq(r1->op_id, "stretch.mtf");
	cr_assert_eq(r1->op_version, 1);
	cr_assert_str_eq(r1->params, "linked=1;lo=0.001;mid=0.25;hi=1");
	cr_assert_str_eq(r1->summary, "Midtone Transfer Function");
	cr_assert_eq(r1->tier, NDE_TIER_A);
	cr_assert_eq(r1->scope, NDE_SCOPE_LAYER);
	cr_assert_eq(r1->target_item_id, lid);
	cr_assert(r1->mask_active);
	cr_assert_not_null(r1->timestamp);
	cr_assert_not_null(r1->impl);
	cr_assert_null(r1->mask_ref);

	nde_record *r2 = g_ptr_array_index(snap, 1);
	cr_assert_eq(r2->record_id, 2);
	cr_assert_str_eq(r2->op_id, "future.shiny_op");
	cr_assert_eq(r2->op_version, 99);
	cr_assert_str_eq(r2->params, "text=semi\\;colon;w=\\=x\\\\y");
	cr_assert_eq(r2->tier, NDE_TIER_B);
	cr_assert_eq(r2->scope, NDE_SCOPE_CANVAS);
	cr_assert_eq(r2->target_item_id, -1);
	cr_assert(!r2->mask_active);
	g_ptr_array_unref(snap);

	/* new records continue the sequence — ids are never reused */
	cr_assert_eq(append_full("a.b", 1, NULL, "x", NDE_TIER_B,
	                         NDE_SCOPE_LAYER, -1, FALSE), 3);
}

Test(nde_persist, only_live_prefix_persists) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	append_full("a.a", 1, NULL, "a", NDE_TIER_B, NDE_SCOPE_LAYER, -1, FALSE);
	gint64 id2 = append_full("b.b", 1, NULL, "b", NDE_TIER_B, NDE_SCOPE_LAYER, -1, FALSE);
	nde_history_on_undo(id2);
	cr_assert_eq(nde_history_live_count(), 1);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert_eq(nde_history_live_count(), 1, "dead tail must not be persisted");
	gint64 next_id = 0;
	GPtrArray *snap = nde_history_snapshot(&next_id);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(snap, 0))->op_id, "a.a");
	cr_assert_eq(next_id, 3, "the undone record's id must stay burned");
	g_ptr_array_unref(snap);
}

Test(nde_persist, no_history_saves_and_loads_empty) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	cr_assert_eq(load_flis(tmppath), 0);
	cr_assert_eq(nde_history_live_count(), 0);
	cr_assert(!nde_history_is_stale());
}

Test(nde_persist, tampered_layer_pixels_set_stale) {
	flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	append_full("a.a", 1, NULL, "a", NDE_TIER_B, NDE_SCOPE_LAYER, -1, FALSE);
	cr_assert_eq(save_flis(tmppath), 0);

	/* Rewrite one pixel of the layer HDU behind the history's back.  With
	 * no ICC HDU the layout is: 1 thumbnail, 2 FLIS_META, 3 layer. */
	int status = 0;
	fitsfile *fptr = NULL;
	fits_open_diskfile(&fptr, tmppath, READWRITE, &status);
	cr_assert_eq(status, 0);
	fits_movabs_hdu(fptr, 3, NULL, &status);
	float poked = 0.75f;
	fits_write_img(fptr, TFLOAT, 1, 1, &poked, &status);
	fits_close_file(fptr, &status);
	cr_assert_eq(status, 0);

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);
	cr_assert_eq(nde_history_live_count(), 1, "records must load even when stale");
	cr_assert(nde_history_is_stale(), "externally modified pixels must mark the history stale");
}

Test(nde_persist, compressed_save_skips_hashes_not_stale) {
	flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.25f), "base");
	append_full("a.a", 1, NULL, "a", NDE_TIER_B, NDE_SCOPE_LAYER, -1, FALSE);
	com.pref.comp.fits_enabled = TRUE;
	com.pref.comp.fits_method = 0;          /* RICE */
	com.pref.comp.fits_quantization = 16.0;
	int rc = save_flis(tmppath);
	com.pref.comp.fits_enabled = FALSE;
	cr_assert_eq(rc, 0);

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);
	cr_assert_eq(nde_history_live_count(), 1);
	cr_assert(!nde_history_is_stale(),
	          "hash-less (compressed) files must not present as stale");
}

Test(nde_persist, declared_but_missing_table_loads_empty) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	append_full("a.a", 1, NULL, "a", NDE_TIER_B, NDE_SCOPE_LAYER, -1, FALSE);
	cr_assert_eq(save_flis(tmppath), 0);

	/* Delete the FLIS_HIST HDU, leaving FLISNDE=T behind. */
	int status = 0;
	fitsfile *fptr = NULL;
	fits_open_diskfile(&fptr, tmppath, READWRITE, &status);
	cr_assert_eq(status, 0);
	int nhdus = 0;
	fits_get_num_hdus(fptr, &nhdus, &status);
	gboolean deleted = FALSE;
	for (int h = 2; h <= nhdus && !deleted; h++) {
		char extname[FLEN_VALUE] = { 0 };
		fits_movabs_hdu(fptr, h, NULL, &status); status = 0;
		fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status); status = 0;
		if (!g_ascii_strcasecmp(extname, "FLIS_HIST")) {
			fits_delete_hdu(fptr, NULL, &status);
			deleted = TRUE;
		}
	}
	fits_close_file(fptr, &status);
	cr_assert(deleted, "fixture: FLIS_HIST HDU not found to delete");
	cr_assert_eq(status, 0);

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0, "missing table must degrade gracefully");
	cr_assert_eq(nde_history_live_count(), 0);
}
