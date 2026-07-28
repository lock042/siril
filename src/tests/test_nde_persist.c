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
#include "core/nde_checkpoint.h"

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

/* ---- NDE_BASE baseline persistence (plan P2.C) ------------------------ */

/* Build a float fits with a per-pixel ramp so bit-exactness is meaningful. */
static fits *ramp_float(int rx, int ry, int nch, float base) {
	fits *f = NULL;
	cr_assert_eq(new_fit_image(&f, rx, ry, nch, DATA_FLOAT), 0);
	size_t n = (size_t)rx * ry * nch;
	for (size_t i = 0; i < n; i++)
		f->fdata[i] = base + (float)i * 3.1e-4f;
	return f;
}

static fits *ramp_ushort(int rx, int ry, int nch, WORD base) {
	fits *f = NULL;
	cr_assert_eq(new_fit_image(&f, rx, ry, nch, DATA_USHORT), 0);
	size_t n = (size_t)rx * ry * nch;
	for (size_t i = 0; i < n; i++)
		f->data[i] = (WORD)(base + (WORD)(i * 13));
	return f;
}

/* Assert two float fits are pixel-identical. */
static void assert_float_exact(const fits *a, const fits *b) {
	cr_assert_not_null(a);
	cr_assert_not_null(b);
	cr_assert_eq(a->rx, b->rx);
	cr_assert_eq(a->ry, b->ry);
	cr_assert_eq(a->naxes[2], b->naxes[2]);
	cr_assert_eq(a->type, DATA_FLOAT);
	cr_assert_eq(b->type, DATA_FLOAT);
	size_t n = (size_t)a->rx * a->ry * a->naxes[2];
	for (size_t i = 0; i < n; i++)
		cr_assert_eq(a->fdata[i], b->fdata[i],
		             "float baseline pixel %zu: %g != %g", i, a->fdata[i], b->fdata[i]);
}

Test(nde_persist, baseline_roundtrip_float_exact) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(12, 10, 0.25f), "base");
	gint lid = l->item_id;
	/* Seed a distinctive baseline for this item and a Tier-A record for it. */
	fits *seed = ramp_float(12, 10, 1, 0.1f);
	nde_checkpoint_baseline_ensure(seed, lid);
	append_full("stretch.asinh", 1, "beta=5;offset=0", "Asinh",
	            NDE_TIER_A, NDE_SCOPE_LAYER, lid, FALSE);

	cr_assert_eq(save_flis(tmppath), 0);

	/* Reload from scratch — attach purges the store, then the loader adopts. */
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert(!nde_checkpoint_baseline_exists(lid), "purge must clear the store");
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert(nde_checkpoint_baseline_exists(lid), "baseline must be adopted on load");
	fits *got = nde_checkpoint_baseline_get(lid);
	assert_float_exact(got, seed);
	clearfits(got); free(got);
	clearfits(seed); free(seed);
}

Test(nde_persist, baseline_roundtrip_ushort_exact) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "base");
	gint lid = l->item_id;
	fits *seed = ramp_ushort(8, 8, 1, 1000);
	nde_checkpoint_baseline_ensure(seed, lid);
	append_full("filters.median", 1, "ksize=3", "Median",
	            NDE_TIER_A, NDE_SCOPE_LAYER, lid, FALSE);

	cr_assert_eq(save_flis(tmppath), 0);

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert(nde_checkpoint_baseline_exists(lid));
	fits *got = nde_checkpoint_baseline_get(lid);
	cr_assert_not_null(got);
	cr_assert_eq(got->type, DATA_USHORT);
	size_t n = (size_t)seed->rx * seed->ry * seed->naxes[2];
	for (size_t i = 0; i < n; i++)
		cr_assert_eq(got->data[i], seed->data[i],
		             "ushort baseline pixel %zu: %u != %u", i, got->data[i], seed->data[i]);
	clearfits(got); free(got);
	clearfits(seed); free(seed);
}

/* Decision 4: the baseline must round-trip LOSSLESSLY even when the user has
 * enabled (lossy) tile compression for ordinary layer HDUs. */
Test(nde_persist, baseline_lossless_under_user_compression) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.25f), "base");
	gint lid = l->item_id;
	fits *seed = ramp_float(16, 16, 1, 0.42f);
	nde_checkpoint_baseline_ensure(seed, lid);
	append_full("stretch.asinh", 1, "beta=5;offset=0", "Asinh",
	            NDE_TIER_A, NDE_SCOPE_LAYER, lid, FALSE);

	/* Lossy compression for the layer HDUs; the baseline must ignore it. */
	com.pref.comp.fits_enabled = TRUE;
	com.pref.comp.fits_method = 0;          /* RICE */
	com.pref.comp.fits_quantization = 16.0;
	int rc = save_flis(tmppath);
	com.pref.comp.fits_enabled = FALSE;
	cr_assert_eq(rc, 0);

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert(nde_checkpoint_baseline_exists(lid));
	fits *got = nde_checkpoint_baseline_get(lid);
	/* Bit-exact despite fits_enabled — GZIP + quantize 0.0 on NDE_BASE. */
	assert_float_exact(got, seed);
	clearfits(got); free(got);
	clearfits(seed); free(seed);
}

/* A file with no history has no NDE_BASE HDUs; a fresh load adopts nothing. */
Test(nde_persist, no_history_no_baselines) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	gint lid = l->item_id;
	cr_assert_eq(save_flis(tmppath), 0);

	/* Scan the saved file: there must be zero NDE_BASE HDUs. */
	int status = 0, nhdus = 0, found = 0;
	fitsfile *fptr = NULL;
	fits_open_diskfile(&fptr, tmppath, READONLY, &status);
	cr_assert_eq(status, 0);
	fits_get_num_hdus(fptr, &nhdus, &status);
	for (int h = 2; h <= nhdus; h++) {
		char extname[FLEN_VALUE] = { 0 };
		fits_movabs_hdu(fptr, h, NULL, &status); status = 0;
		fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status); status = 0;
		if (!g_ascii_strcasecmp(extname, "NDE_BASE")) found++;
	}
	fits_close_file(fptr, &status);
	cr_assert_eq(found, 0, "a history-less file must carry no NDE_BASE HDUs");

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);
	cr_assert_not(nde_checkpoint_baseline_exists(lid));
}

/* Pre-phase-2 file (FLIS_HIST present, NDE_BASE HDUs deleted): loads fine and
 * adopts no baseline — the chain is simply not replayable (decision 5). */
Test(nde_persist, pre_phase2_no_baseline_loads) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;
	fits *seed = ramp_float(8, 8, 1, 0.3f);
	nde_checkpoint_baseline_ensure(seed, lid);
	append_full("stretch.asinh", 1, "beta=5", "Asinh",
	            NDE_TIER_A, NDE_SCOPE_LAYER, lid, FALSE);
	cr_assert_eq(save_flis(tmppath), 0);

	/* Delete every NDE_BASE HDU, mimicking a pre-phase-2 writer. */
	int status = 0, nhdus = 0;
	fitsfile *fptr = NULL;
	fits_open_diskfile(&fptr, tmppath, READWRITE, &status);
	cr_assert_eq(status, 0);
	fits_get_num_hdus(fptr, &nhdus, &status);
	gboolean deleted = FALSE;
	for (int h = 2; h <= nhdus; ) {
		char extname[FLEN_VALUE] = { 0 };
		fits_movabs_hdu(fptr, h, NULL, &status); status = 0;
		fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status); status = 0;
		if (!g_ascii_strcasecmp(extname, "NDE_BASE")) {
			fits_delete_hdu(fptr, NULL, &status); status = 0;
			deleted = TRUE;
			fits_get_num_hdus(fptr, &nhdus, &status);
		} else {
			h++;
		}
	}
	fits_close_file(fptr, &status);
	cr_assert(deleted, "fixture: no NDE_BASE HDU found to delete");

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0, "missing baselines must load gracefully");
	cr_assert_eq(nde_history_live_count(), 1, "history still loads");
	cr_assert_not(nde_checkpoint_baseline_exists(lid),
	              "no baseline HDU → not replayable");
	clearfits(seed); free(seed);
}

/* ---- NDE_CKPT output-checkpoint persistence (plan P4.3) --------------- */

/* Count NDE_CKPT HDUs in the saved file. */
static int count_ckpt_hdus(const char *path) {
	int status = 0, nhdus = 0, found = 0;
	fitsfile *fptr = NULL;
	fits_open_diskfile(&fptr, path, READONLY, &status);
	cr_assert_eq(status, 0);
	fits_get_num_hdus(fptr, &nhdus, &status);
	for (int h = 2; h <= nhdus; h++) {
		char extname[FLEN_VALUE] = { 0 };
		fits_movabs_hdu(fptr, h, NULL, &status); status = 0;
		fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status); status = 0;
		if (!g_ascii_strcasecmp(extname, "NDE_CKPT")) found++;
	}
	fits_close_file(fptr, &status);
	return found;
}

/* A checkpointed barrier round-trips its output checkpoint bit-exactly
 * (float), keyed by record_id. */
Test(nde_persist, ckpt_roundtrip_float_exact) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(12, 10, 0.25f), "base");
	gint lid = l->item_id;
	/* A Tier-B barrier record + its output checkpoint. */
	gint64 bid = append_full("python.set_pixeldata", 1, NULL, "freehand",
	                         NDE_TIER_B, NDE_SCOPE_LAYER, lid, FALSE);
	fits *post = ramp_float(12, 10, 1, 0.7f);
	nde_checkpoint_output_store(post, bid, lid);

	cr_assert_eq(save_flis(tmppath), 0);
	cr_assert_eq(count_ckpt_hdus(tmppath), 1, "one barrier checkpoint HDU expected");

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert(!nde_checkpoint_output_exists(bid), "purge must clear the store");
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert(nde_checkpoint_output_exists(bid), "checkpoint must be adopted on load");
	fits *got = nde_checkpoint_output_get(bid);
	assert_float_exact(got, post);
	clearfits(got); free(got);
	clearfits(post); free(post);
}

/* Same, ushort — cheap enough to check both types. */
Test(nde_persist, ckpt_roundtrip_ushort_exact) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "base");
	gint lid = l->item_id;
	gint64 bid = append_full("color.rgb_align", 1, NULL, "RGB align",
	                         NDE_TIER_B, NDE_SCOPE_LAYER, lid, FALSE);
	fits *post = ramp_ushort(8, 8, 1, 2000);
	nde_checkpoint_output_store(post, bid, lid);

	cr_assert_eq(save_flis(tmppath), 0);
	cr_assert_eq(count_ckpt_hdus(tmppath), 1);

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert(nde_checkpoint_output_exists(bid));
	fits *got = nde_checkpoint_output_get(bid);
	cr_assert_not_null(got);
	cr_assert_eq(got->type, DATA_USHORT);
	size_t n = (size_t)post->rx * post->ry * post->naxes[2];
	for (size_t i = 0; i < n; i++)
		cr_assert_eq(got->data[i], post->data[i],
		             "ushort checkpoint pixel %zu: %u != %u", i, got->data[i], post->data[i]);
	clearfits(got); free(got);
	clearfits(post); free(post);
}

/* The checkpoint must round-trip LOSSLESSLY under the user's lossy tile
 * compression, exactly like NDE_BASE (decision 4). */
Test(nde_persist, ckpt_lossless_under_user_compression) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(16, 16, 0.25f), "base");
	gint lid = l->item_id;
	gint64 bid = append_full("python.set_pixeldata", 1, NULL, "freehand",
	                         NDE_TIER_B, NDE_SCOPE_LAYER, lid, FALSE);
	fits *post = ramp_float(16, 16, 1, 0.33f);
	nde_checkpoint_output_store(post, bid, lid);

	com.pref.comp.fits_enabled = TRUE;
	com.pref.comp.fits_method = 0;          /* RICE */
	com.pref.comp.fits_quantization = 16.0;
	int rc = save_flis(tmppath);
	com.pref.comp.fits_enabled = FALSE;
	cr_assert_eq(rc, 0);

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	cr_assert(nde_checkpoint_output_exists(bid));
	fits *got = nde_checkpoint_output_get(bid);
	assert_float_exact(got, post);   /* bit-exact despite fits_enabled */
	clearfits(got); free(got);
	clearfits(post); free(post);
}

/* A barrier WITHOUT a checkpoint (pre-phase-4 capture) writes no HDU; a
 * non-barrier Tier-A record never writes one either. */
Test(nde_persist, barrier_without_checkpoint_writes_no_hdu) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;
	/* barrier, but no checkpoint stored */
	append_full("python.set_pixeldata", 1, NULL, "freehand",
	            NDE_TIER_B, NDE_SCOPE_LAYER, lid, FALSE);
	/* a Tier-A record — not a barrier, never checkpointed */
	append_full("stretch.asinh", 1, "beta=5", "Asinh",
	            NDE_TIER_A, NDE_SCOPE_LAYER, lid, FALSE);

	cr_assert_eq(save_flis(tmppath), 0);
	cr_assert_eq(count_ckpt_hdus(tmppath), 0,
	             "no in-session checkpoint → no NDE_CKPT HDU");
}

/* A deleted/truncated barrier's checkpoint is dropped from the store, so it
 * is not written on the next save (only LIVE, retained checkpoints persist). */
Test(nde_persist, truncated_barrier_checkpoint_not_written) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;
	gint64 bid = append_full("python.set_pixeldata", 1, NULL, "freehand",
	                         NDE_TIER_B, NDE_SCOPE_LAYER, lid, FALSE);
	fits *post = ramp_float(8, 8, 1, 0.5f);
	nde_checkpoint_output_store(post, bid, lid);
	cr_assert(nde_checkpoint_output_exists(bid));

	/* Undo the barrier then append: the dead-tail truncation drops its
	 * checkpoint (P4.2 lifecycle) — the surviving record is not a barrier. */
	nde_history_on_undo(bid);
	append_full("stretch.asinh", 1, "beta=5", "Asinh",
	            NDE_TIER_A, NDE_SCOPE_LAYER, lid, FALSE);
	cr_assert(!nde_checkpoint_output_exists(bid),
	          "truncation must drop the barrier's checkpoint");

	cr_assert_eq(save_flis(tmppath), 0);
	cr_assert_eq(count_ckpt_hdus(tmppath), 0,
	             "a truncated barrier's checkpoint must not be persisted");
	clearfits(post); free(post);
}

/* Pre-phase-4 file (NDE_CKPT HDUs deleted): loads fine, adopts no checkpoint
 * — the barrier stays a full blocker (decision 5 analogue). */
Test(nde_persist, pre_phase4_no_ckpt_loads) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;
	gint64 bid = append_full("python.set_pixeldata", 1, NULL, "freehand",
	                         NDE_TIER_B, NDE_SCOPE_LAYER, lid, FALSE);
	fits *post = ramp_float(8, 8, 1, 0.6f);
	nde_checkpoint_output_store(post, bid, lid);
	cr_assert_eq(save_flis(tmppath), 0);

	/* Delete every NDE_CKPT HDU, mimicking a pre-phase-4 writer. */
	int status = 0, nhdus = 0;
	fitsfile *fptr = NULL;
	fits_open_diskfile(&fptr, tmppath, READWRITE, &status);
	cr_assert_eq(status, 0);
	fits_get_num_hdus(fptr, &nhdus, &status);
	gboolean deleted = FALSE;
	for (int h = 2; h <= nhdus; ) {
		char extname[FLEN_VALUE] = { 0 };
		fits_movabs_hdu(fptr, h, NULL, &status); status = 0;
		fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status); status = 0;
		if (!g_ascii_strcasecmp(extname, "NDE_CKPT")) {
			fits_delete_hdu(fptr, NULL, &status); status = 0;
			deleted = TRUE;
			fits_get_num_hdus(fptr, &nhdus, &status);
		} else {
			h++;
		}
	}
	fits_close_file(fptr, &status);
	cr_assert(deleted, "fixture: no NDE_CKPT HDU found to delete");

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0, "missing checkpoints must load gracefully");
	cr_assert_eq(nde_history_live_count(), 1, "history still loads");
	cr_assert_not(nde_checkpoint_output_exists(bid),
	              "no NDE_CKPT HDU → barrier stays a full blocker");
	clearfits(post); free(post);
}

Test(nde_persist, reordered_history_persists_in_new_order) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	append_full("a.a", 1, "k=1", "a", NDE_TIER_A, NDE_SCOPE_LAYER, -1, FALSE);
	append_full("b.b", 1, "k=2", "b", NDE_TIER_A, NDE_SCOPE_LAYER, -1, FALSE);
	append_full("c.c", 1, "k=3", "c", NDE_TIER_A, NDE_SCOPE_LAYER, -1, FALSE);

	/* log-level move: record 3 before record 1 → order [3, 1, 2] */
	gchar *err = NULL;
	cr_assert(nde_history_reorder(3, 1, &err), "reorder failed: %s", err ? err : "?");

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 3);
	cr_assert_eq(((nde_record *)g_ptr_array_index(snap, 0))->record_id, 3);
	cr_assert_eq(((nde_record *)g_ptr_array_index(snap, 1))->record_id, 1);
	cr_assert_eq(((nde_record *)g_ptr_array_index(snap, 2))->record_id, 2);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(snap, 0))->op_id, "c.c");
	g_ptr_array_unref(snap);
}

/* ---- input pins (graph step 4) ---------------------------------------- */

Test(nde_persist, input_pins_round_trip) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	gint64 id = append_full("filters.gauss", 1, "sigma=2", "blur",
	                        NDE_TIER_A, NDE_SCOPE_LAYER, 1, TRUE);
	cr_assert_eq(id, 1);
	/* Reach into the live record: capture sites add pins before append, but
	 * the fixture appends first, so add through a snapshot-free path. */
	GPtrArray *live = nde_history_snapshot(NULL);
	cr_assert_eq(live->len, 1);
	g_ptr_array_unref(live);

	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("filters.unsharp");
	rec->op_version = 1;
	rec->params = g_strdup("sigma=1;amount=0.5");
	rec->summary = g_strdup("unsharp");
	rec->tier = NDE_TIER_A;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = 1;
	rec->timestamp = nde_iso8601_now();
	rec->impl = nde_impl_string();
	nde_record_add_input(rec, "mask", 7, 3);
	nde_record_add_input(rec, "overlay", -2, 0);   /* 0 = that item's baseline */
	cr_assert(nde_history_append(rec) > 0);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_eq(snap->len, 2);
	const nde_record *plain = g_ptr_array_index(snap, 0);
	cr_assert_null(plain->inputs, "a record with no pins must load none");
	const nde_record *pinned = g_ptr_array_index(snap, 1);
	cr_assert_not_null(pinned->inputs);
	cr_assert_eq(pinned->inputs->len, 2);
	const nde_input_pin *m = nde_record_input(pinned, "mask");
	cr_assert_not_null(m);
	cr_assert_eq(m->src_item_id, 7);
	cr_assert_eq(m->src_record_id, 3);
	const nde_input_pin *o = nde_record_input(pinned, "overlay");
	cr_assert_not_null(o);
	cr_assert_eq(o->src_item_id, -2);
	cr_assert_eq(o->src_record_id, 0);
	cr_assert_null(nde_record_input(pinned, "nosuchrole"));
	g_ptr_array_unref(snap);
}

/* A file written before the INPUTS column existed is OLDER, not malformed:
 * the column lookup must tolerate its absence rather than discard the whole
 * history.  Simulated by deleting the column from a file we just wrote. */
Test(nde_persist, history_without_the_inputs_column_still_loads) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	append_full("stretch.mtf", 1, "lo=0.1", "mtf", NDE_TIER_A,
	            NDE_SCOPE_LAYER, 1, FALSE);
	nde_record *rec = nde_record_new();
	rec->op_id = g_strdup("filters.gauss");
	rec->op_version = 1;
	rec->params = g_strdup("sigma=2");
	rec->summary = g_strdup("blur");
	rec->tier = NDE_TIER_A;
	rec->scope = NDE_SCOPE_LAYER;
	rec->target_item_id = 1;
	rec->timestamp = nde_iso8601_now();
	rec->impl = nde_impl_string();
	nde_record_add_input(rec, "mask", 7, 1);
	cr_assert(nde_history_append(rec) > 0);
	cr_assert_eq(save_flis(tmppath), 0);

	int status = 0;
	fitsfile *fptr = NULL;
	fits_open_diskfile(&fptr, tmppath, READWRITE, &status);
	cr_assert_eq(status, 0);
	int nhdus = 0;
	fits_get_num_hdus(fptr, &nhdus, &status);
	gboolean dropped = FALSE;
	for (int h = 2; h <= nhdus && !dropped; h++) {
		char extname[FLEN_VALUE] = { 0 };
		fits_movabs_hdu(fptr, h, NULL, &status); status = 0;
		fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status); status = 0;
		if (!g_ascii_strcasecmp(extname, "FLIS_HIST")) {
			int col = 0;
			fits_get_colnum(fptr, CASEINSEN, "INPUTS", &col, &status);
			cr_assert_eq(status, 0, "fixture: INPUTS column not found");
			fits_delete_col(fptr, col, &status);
			cr_assert_eq(status, 0);
			/* and say it is a v1 table, as an old writer would have */
			int v1 = 1;
			fits_update_key(fptr, TINT, "FLISHVER", &v1, NULL, &status);
			dropped = TRUE;
		}
	}
	fits_close_file(fptr, &status);
	cr_assert(dropped, "fixture: FLIS_HIST not found");
	cr_assert_eq(status, 0);

	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0, "an older table must load, not be discarded");
	GPtrArray *snap = nde_history_snapshot(NULL);
	cr_assert_not_null(snap, "the history must survive, minus the pins");
	cr_assert_eq(snap->len, 2);
	cr_assert_str_eq(((nde_record *)g_ptr_array_index(snap, 0))->op_id, "stretch.mtf");
	cr_assert_null(((nde_record *)g_ptr_array_index(snap, 1))->inputs);
	g_ptr_array_unref(snap);
}
