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
 * test_flis_lifecycle — basic helper exercise for FLIS data-model
 * primitives.  Covers flis_layer_new, flis_layer_free, layermask_free,
 * flis_layer_add / _remove / _count, flis_group_add / _remove / _count.
 *
 * Pure data-model tests — no save/load or compositing.  Round-trip and
 * pixel-correctness tests live in dedicated files.
 */

#include <criterion/criterion.h>
#include "flis_test_helpers.h"

cominfo com;
fits *gfit;

TestSuite(flis_lifecycle, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

Test(flis_lifecycle, layer_new_with_name_sets_defaults) {
	fits *f = flis_test_make_mono_fits(8, 8, 0.5f);
	cr_assert_not_null(f, "fits creation failed");

	flis_layer_t *lay = flis_layer_new(f, "my layer");
	cr_assert_not_null(lay, "flis_layer_new returned NULL");
	cr_assert_eq(lay->fit, f, "layer->fit should equal the input fits");
	cr_assert_str_eq(lay->layer_name, "my layer");
	cr_assert_eq(lay->blend_mode, FLIS_BLEND_NORMAL);
	cr_assert_float_eq(lay->opacity, 1.0f, 1e-6);
	cr_assert(lay->visible);
	cr_assert(!lay->locked);
	cr_assert(!lay->has_tint);
	cr_assert(lay->lmask_active);
	cr_assert_null(lay->lmask);
	cr_assert_eq(lay->group_id, 0);
	cr_assert_eq(lay->position_x, 0);
	cr_assert_eq(lay->position_y, 0);

	flis_layer_free(lay);
}

Test(flis_lifecycle, layer_new_with_null_name_defaults_to_Layer) {
	fits *f = flis_test_make_mono_fits(8, 8, 0.0f);
	flis_layer_t *lay = flis_layer_new(f, NULL);
	cr_assert_not_null(lay);
	cr_assert_str_eq(lay->layer_name, "Layer");
	flis_layer_free(lay);
}

Test(flis_lifecycle, layer_free_handles_null) {
	flis_layer_free(NULL);  /* must not crash */
}

Test(flis_lifecycle, layermask_free_handles_null) {
	layermask_free(NULL);  /* must not crash */
}

Test(flis_lifecycle, layermask_free_frees_data) {
	layermask_t *m = flis_test_make_const_lmask(4, 4, 8, 1.0);
	cr_assert_not_null(m);
	cr_assert_not_null(m->data);
	cr_assert_eq(m->w, 4);
	cr_assert_eq(m->h, 4);
	cr_assert_eq(m->bitpix, 8);
	layermask_free(m);  /* must not leak — verified by ASAN if enabled */
}

Test(flis_lifecycle, layer_add_assigns_item_id_and_order) {
	flis_layer_t *base = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 0.0f), "base");
	cr_assert_not_null(base);
	cr_assert_eq(flis_layer_count(), 1);
	cr_assert_eq(base->item_id, 1, "first item should get item_id=1");

	flis_layer_t *top = flis_test_add_layer(flis_test_make_mono_fits(8, 8, 1.0f), "top");
	cr_assert_not_null(top);
	cr_assert_eq(flis_layer_count(), 2);
	cr_assert_eq(top->item_id, 2, "second item should get item_id=2");
	cr_assert_gt(top->layer_order, base->layer_order,
	             "newly added layer should stack above existing");
}

Test(flis_lifecycle, layer_remove_drops_layer_and_keeps_item_ids_unique) {
	flis_layer_t *a = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "a");
	flis_layer_t *b = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "b");
	cr_assert_eq(flis_layer_count(), 2);
	gint b_id = b->item_id;

	cr_assert_eq(flis_layer_remove(b), 0);
	cr_assert_eq(flis_layer_count(), 1);
	cr_assert_eq(a->item_id, 1, "remaining layer keeps its id");

	/* A new layer must not reuse the freed id. */
	flis_layer_t *c = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.5f), "c");
	cr_assert_neq(c->item_id, b_id, "item_id must never be reused");
	cr_assert_eq(c->item_id, 3, "next assignment is monotonic");
}

Test(flis_lifecycle, layer_remove_refuses_last_layer) {
	flis_layer_t *only = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "only");
	cr_assert_neq(flis_layer_remove(only), 0,
	              "removing the last layer should fail (a FLIS must keep ≥1 layer)");
	cr_assert_eq(flis_layer_count(), 1);
}

Test(flis_lifecycle, lookup_by_id_and_name) {
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "Background");
	flis_layer_t *ha = flis_test_add_layer(flis_test_make_mono_fits(4, 4, 1.0f), "Hα");

	cr_assert_eq(flis_layer_get_by_id(ha->item_id), ha);
	cr_assert_eq(flis_layer_get_by_name("Hα"), ha);
	cr_assert_null(flis_layer_get_by_name("does not exist"));
	cr_assert_null(flis_layer_get_by_id(99));
}

Test(flis_lifecycle, group_add_and_remove) {
	/* Need a layer first — group ops require a FLIS */
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "base");

	cr_assert_eq(flis_group_count(), 0);
	flis_group_t *g = flis_group_add("My Group");
	cr_assert_not_null(g);
	cr_assert_eq(flis_group_count(), 1);
	cr_assert_str_eq(g->name, "My Group");

	cr_assert_eq(flis_group_remove(g), 0);
	cr_assert_eq(flis_group_count(), 0);
}

Test(flis_lifecycle, is_current_image_flis_reflects_layers_state) {
	cr_assert(!is_current_image_flis(), "no layers → not FLIS");
	flis_test_add_layer(flis_test_make_mono_fits(4, 4, 0.0f), "base");
	cr_assert(is_current_image_flis(), "one layer → FLIS");
}
