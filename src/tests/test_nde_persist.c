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
#include "core/nde/nde_history.h"
#include "core/nde/nde_checkpoint.h"
#include "core/nde/nde_cat.h"
#include "core/op_descriptor.h"
#include "core/processing.h"

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

/* Build a small catalogue the way a PCC/SPCC run leaves one in memory. */
static siril_catalogue *make_test_catalogue(gboolean with_spectra) {
	siril_catalogue *cat = siril_catalog_new(CAT_LOCAL_GAIA_XPSAMP);
	cat->center_ra = 123.456;
	cat->center_dec = -54.321;
	cat->radius = 42.5;
	cat->limitmag = 17.25;
	cat->phot = TRUE;
	cat->columns = (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) |
	               (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) |
	               (1 << CAT_FIELD_TEFF);
	cat->nbitems = 2;
	cat->cat_items = calloc(2, sizeof(cat_item));
	for (int i = 0; i < 2; i++) {
		cat_item *it = &cat->cat_items[i];
		it->ra = 123.4 + i * 0.01;
		it->dec = -54.3 - i * 0.01;
		it->pmra = 1.5 * (i + 1);
		it->pmdec = -0.75 * (i + 1);
		it->mag = 12.5f + i;
		it->bmag = 13.25f + i;
		it->e_mag = 0.01f;
		it->e_bmag = 0.02f;
		it->teff = 5500.0f + 100 * i;
		it->gaiasourceid = 4000000000000000000ULL + i;
		it->name = g_strdup_printf("star%d", i);
		it->included = TRUE;
	}
	if (with_spectra) {
		cat->columns |= 1 << CAT_FIELD_XPSAMP;
		cat->cat_items[0].xp_sampled = malloc(XPSAMPLED_LEN * sizeof(double));
		for (int k = 0; k < XPSAMPLED_LEN; k++)
			cat->cat_items[0].xp_sampled[k] = 0.5 * k + 0.125;
	}
	return cat;
}

/* An SPCC record's embedded star catalogue survives save → load: the file
 * is self-contained, so the recompute runs offline on any machine. */
Test(nde_persist, embedded_catalogue_roundtrip) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint64 rid = append_full("color.photometric_cc", 2,
	                         "catalog=102;spcc=1;auto_replay=1", "SPCC",
	                         NDE_TIER_A, NDE_SCOPE_LAYER, l->item_id, FALSE);
	cr_assert(rid > 0);
	nde_cat_register(rid, make_test_catalogue(TRUE));
	cr_assert(nde_cat_has(rid));

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);   /* purges the registry too */
	cr_assert(!nde_cat_has(rid), "attach must purge the catalogue registry");
	cr_assert_eq(load_flis(tmppath), 0);

	siril_catalogue *back = nde_cat_get_copy(rid);
	cr_assert_not_null(back, "the embedded catalogue must reload");
	cr_assert_eq(back->cat_index, CAT_LOCAL_GAIA_XPSAMP);
	cr_assert_float_eq(back->center_ra, 123.456, 1e-9);
	cr_assert_float_eq(back->center_dec, -54.321, 1e-9);
	cr_assert_float_eq(back->radius, 42.5, 1e-9);
	cr_assert_float_eq(back->limitmag, 17.25, 1e-9);
	cr_assert(back->phot);
	cr_assert(has_field(back, XPSAMP));
	cr_assert_eq(back->nbitems, 2);
	const cat_item *it = &back->cat_items[0];
	cr_assert_float_eq(it->ra, 123.4, 1e-9);
	cr_assert_float_eq(it->pmdec, -0.75, 1e-9);
	cr_assert_float_eq(it->mag, 12.5f, 1e-6);
	cr_assert_float_eq(it->bmag, 13.25f, 1e-6);
	cr_assert_float_eq(it->teff, 5500.0f, 1e-3);
	cr_assert_eq(it->gaiasourceid, 4000000000000000000ULL);
	cr_assert_str_eq(it->name, "star0");
	cr_assert_not_null(it->xp_sampled, "the flux spectrum must survive the trip");
	double expect[XPSAMPLED_LEN];
	for (int k = 0; k < XPSAMPLED_LEN; k++)
		expect[k] = 0.5 * k + 0.125;
	cr_assert(memcmp(it->xp_sampled, expect, sizeof(expect)) == 0,
	          "spectra must round-trip bit-exactly");
	cr_assert_null(back->cat_items[1].xp_sampled,
	               "a star without a spectrum must stay without one");
	siril_catalog_free(back);
}

/* The stash → adopt handshake: only a photometric capture claims a pending
 * catalogue; anything else discards it as stale. */
Test(nde_persist, catalogue_stash_adopt_semantics) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	siril_catalogue *cat = make_test_catalogue(FALSE);

	nde_cat_stash_pending(cat);
	gint64 other = nde_capture_structural("layer.set_opacity", NDE_SCOPE_LAYER,
	                                      1, g_strdup("opacity=0.5"), "opacity");
	cr_assert(other > 0);
	cr_assert(!nde_cat_has(other), "a non-photometric capture must not adopt");

	nde_cat_stash_pending(cat);
	gint64 pcc = nde_capture_structural("color.photometric_cc", NDE_SCOPE_LAYER,
	                                    1, g_strdup("spcc=0"), "PCC");
	cr_assert(pcc > 0);
	cr_assert(!nde_cat_has(other));
	cr_assert(nde_cat_has(pcc), "the photometric capture must adopt the stash");

	siril_catalogue *back = nde_cat_get_copy(pcc);
	cr_assert_not_null(back);
	cr_assert_eq(back->nbitems, 2);
	siril_catalog_free(back);
	siril_catalog_free(cat);
}

/* The stand-in for op_desc_photometric_cc.  The capture reads only op->id,
 * op->version and op->serialize, so a descriptor that shares the id is enough
 * to exercise the handshake — and it keeps the test off the network and away
 * from the WCS the real hook needs. */
static int noop_image_hook(struct generic_img_args *args, fits *fit, int nb) {
	(void)args; (void)fit; (void)nb;
	return 0;
}

static const op_descriptor op_desc_fake_pcc = {
	.id = "color.photometric_cc",
	.version = 2,
	.image_hook = noop_image_hook,
	.description = "PCC",
};

/* The test above passes while the real thing is broken, because it captures
 * through nde_capture_structural — a capture_finish path, and capture_finish is
 * where the adopt lives.  Single-image PCC does NOT go that way: it runs
 * through generic_image_worker (gui-gtk4/photometric_cc.c), whose capture block
 * is a fourth, separate implementation that never adopted anything.  So the
 * stash photometric_cc.c takes "so the record can recompute offline forever"
 * was dropped on the floor for every non-group run, and replay silently fell
 * back to a fresh network query. */
Test(nde_persist, a_worker_captured_photometric_op_also_adopts_the_stash) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	gfit = flis_active_layer_fit();

	nde_cat_stash_pending(make_test_catalogue(FALSE));

	struct generic_img_args *args = calloc(1, sizeof(*args));
	args->fit = gfit;
	args->op = &op_desc_fake_pcc;
	args->user = NULL;
	args->command = TRUE;
	args->command_updates_gfit = TRUE;
	args->max_threads = 1;
	args->mem_ratio = -1.0f;
	gboolean prev = com.headless;
	com.headless = TRUE;
	cr_assert_eq(GPOINTER_TO_INT(generic_image_worker(args)), 0);
	com.headless = prev;

	gint64 rid = 0;
	GPtrArray *snap = nde_history_snapshot(NULL);
	for (guint i = 0; i < snap->len; i++) {
		const nde_record *r = g_ptr_array_index(snap, i);
		if (!g_strcmp0(r->op_id, "color.photometric_cc"))
			rid = r->record_id;
	}
	g_ptr_array_unref(snap);
	cr_assert(rid > 0, "the worker must have captured a record");

	cr_assert(nde_cat_has(rid),
	          "a photometric record captured by the worker must adopt the "
	          "pending catalogue, exactly as a dialog-captured one does");
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
	fits *got = nde_state_release(nde_checkpoint_baseline_get(lid));
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
	fits *got = nde_state_release(nde_checkpoint_baseline_get(lid));
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
	fits *got = nde_state_release(nde_checkpoint_baseline_get(lid));
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
	fits *got = nde_state_release(nde_checkpoint_output_get(bid));
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
	fits *got = nde_state_release(nde_checkpoint_output_get(bid));
	cr_assert_not_null(got);
	cr_assert_eq(got->type, DATA_USHORT);
	size_t n = (size_t)post->rx * post->ry * post->naxes[2];
	for (size_t i = 0; i < n; i++)
		cr_assert_eq(got->data[i], post->data[i],
		             "ushort checkpoint pixel %zu: %u != %u", i, got->data[i], post->data[i]);
	clearfits(got); free(got);
	clearfits(post); free(post);
}

/* A pinned 8-bit mask travels as a 16-bit checkpoint, and only orig_bitpix
 * records that its content is really 8-bit — BITPIX cannot say so.  Without
 * FLIS_OBPX the reader inferred 16-bit and fits_to_mask rebuilt the wrong
 * depth, so a saved mask came back a different kind of mask from the one
 * saved. */
Test(nde_persist, an_eight_bit_checkpoint_reloads_as_eight_bit) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.5f), "base");
	gint lid = l->item_id;
	gint64 bid = append_full("mask.blur", 1, NULL, "Blur mask",
	                         NDE_TIER_A, NDE_SCOPE_LAYER, lid, FALSE);
	fits *post = ramp_ushort(8, 8, 1, 2000);
	post->orig_bitpix = BYTE_IMG;
	nde_checkpoint_output_store(post, bid, lid);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	cr_assert_eq(load_flis(tmppath), 0);

	fits *got = nde_state_release(nde_checkpoint_output_get(bid));
	cr_assert_not_null(got);
	cr_assert_eq(got->orig_bitpix, BYTE_IMG,
	             "an 8-bit mask must not come back as a 16-bit one");
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
	fits *got = nde_state_release(nde_checkpoint_output_get(bid));
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
	const nde_pin *m = nde_record_input(pinned, "mask");
	cr_assert_not_null(m);
	cr_assert_eq(m->item_id, 7);
	cr_assert_eq(m->record_id, 3);
	const nde_pin *o = nde_record_input(pinned, "overlay");
	cr_assert_not_null(o);
	cr_assert_eq(o->item_id, -2);
	cr_assert_eq(o->record_id, 0);
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

/* Output checkpoints used to be written only for barrier records, because a
 * barrier was the only reason to keep one.  A mask referenced by another
 * record's input pin is kept under its own (non-barrier) record id, and
 * dropping it on save would silently turn its consumer back into a barrier on
 * reload.  The rule is now "has one" rather than "is a barrier". */
Test(nde_persist, checkpoints_of_non_barrier_records_persist) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "base");
	gint lid = l->item_id;
	/* A perfectly ordinary Tier-A record with no mask: not a barrier. */
	gint64 id = append_full("filters.gauss", 1, "sigma=2", "blur", NDE_TIER_A,
	                        NDE_SCOPE_LAYER, lid, FALSE);
	cr_assert(id > 0);
	fits *state = ramp_float(4, 4, 1, 0.125f);
	nde_checkpoint_output_store(state, id, lid);
	cr_assert(nde_checkpoint_output_exists(id));

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	nde_checkpoint_purge();
	cr_assert(!nde_checkpoint_output_exists(id), "fixture: store must be empty");

	cr_assert_eq(load_flis(tmppath), 0);
	cr_assert(nde_checkpoint_output_exists(id),
	          "a stored checkpoint must survive the save, barrier or not");
	fits *got = nde_state_release(nde_checkpoint_output_get(id));
	cr_assert_not_null(got);
	assert_float_exact(got, state);
	clearfits(got); free(got);
	clearfits(state); free(state);
}

/* A layer's replay value is pixels AND its position, so the position has to
 * travel with the checkpoint — otherwise a geometry chain becomes a blocker
 * the moment the file is reopened. */
Test(nde_persist, baseline_layer_offset_round_trips) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;
	append_full("geometry.crop", 1, "x=1;y=2;w=4;h=4", "crop", NDE_TIER_A,
	            NDE_SCOPE_CANVAS, lid, FALSE);
	fits *seed = ramp_float(8, 8, 1, 0.5f);
	/* Where the layer stands when the baseline is taken IS the baseline's
	 * position — there is nothing separate to record. */
	l->position_x = 17;
	l->position_y = 23;
	nde_checkpoint_baseline_ensure(seed, lid);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	nde_checkpoint_purge();
	cr_assert(!nde_checkpoint_baseline_has_position(lid));

	cr_assert_eq(load_flis(tmppath), 0);
	gint x = 0, y = 0;
	cr_assert(nde_checkpoint_baseline_position(lid, &x, &y),
	          "the baseline's position must survive the save");
	cr_assert_eq(x, 17);
	cr_assert_eq(y, 23);
	clearfits(seed); free(seed);
}

/* A restart point is a layer value too. */
Test(nde_persist, checkpoint_layer_offset_round_trips) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;
	gint64 id = append_full("filters.gauss", 1, "sigma=2", "blur", NDE_TIER_A,
	                        NDE_SCOPE_LAYER, lid, FALSE);
	fits *state = ramp_float(8, 8, 1, 0.125f);
	/* Post-op pixels go with wherever the op left the layer. */
	l->position_x = 5;
	l->position_y = 11;
	nde_checkpoint_output_store(state, id, lid);

	cr_assert_eq(save_flis(tmppath), 0);
	flis_free_layers(com.uniq);
	nde_history_attach(NULL);
	nde_checkpoint_purge();

	cr_assert_eq(load_flis(tmppath), 0);
	nde_state *back = nde_checkpoint_output_get(id);
	cr_assert_not_null(back);
	cr_assert(back->has_pos, "the restart point's position must survive the save");
	cr_assert_eq(back->pos_x, 5);
	cr_assert_eq(back->pos_y, 11);
	nde_state_free(back);
	clearfits(state); free(state);
}

/* A position is part of a stored value, so the two cannot come apart: there is
 * no way to record one for a baseline that does not exist, and dropping the
 * pixels takes the position with them.  A position outliving what it describes
 * would silently mis-anchor a later replay. */
Test(nde_persist, a_position_cannot_outlive_its_baseline) {
	flis_layer_t *l = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.25f), "base");
	gint lid = l->item_id;
	l->position_x = 1;
	l->position_y = 2;
	cr_assert(!nde_checkpoint_baseline_has_position(lid),
	          "no baseline for this item, so no position either");

	fits *seed = ramp_float(8, 8, 1, 0.5f);
	nde_checkpoint_baseline_ensure(seed, lid);
	gint x = 0, y = 0;
	cr_assert(nde_checkpoint_baseline_position(lid, &x, &y));
	cr_assert_eq(x, 1);
	cr_assert_eq(y, 2);
	nde_checkpoint_drop(lid);
	cr_assert(!nde_checkpoint_baseline_has_position(lid),
	          "dropping the pixels must drop the position with them");
	clearfits(seed); free(seed);
}
