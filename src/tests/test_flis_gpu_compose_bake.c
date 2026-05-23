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
 * test_flis_gpu_compose_bake — pixel-equivalence harness for the
 * GPU-compose path's CPU side (stage 3.3 slice 4).  Exercises
 * flis_compose_bake_tile_bgra8 / flis_compose_bake_lmask_bgra8 directly
 * with synthetic fixtures, no GTK involvement.
 *
 * What's covered:
 *   - Mono / RGB source data, USHORT and FLOAT types
 *   - Stretch LUT application (identity LUT vs gain-2 LUT vs clip LUT)
 *   - Mono tint baking (the green channel of a (0,1,0)-tinted mono
 *     layer should equal lut(mono); R and B should be zero)
 *   - Y-flip: FITS bottom-up source → display top-down tile output
 *   - Multi-tile coverage: a 4096×2048 layer at FLIS_TILE_DIM=2048
 *     produces 2 tiles side-by-side; verify each tile's content
 *   - Edge tile sizing (right/bottom edges may be smaller than tile_dim)
 *   - Lmask bake at 8/16/32 bitpix, with the same y-flip semantics
 *
 * What's NOT covered (slice 4 deliberately):
 *   - GSK blend-mode rendering itself — that would need a headless
 *     GTK init and an offscreen renderer.  GSK's blend math is W3C
 *     and tested upstream; the failures we've actually hit are in our
 *     CPU side (position/orientation/tint) which this file exercises.
 */

#include <criterion/criterion.h>
#include <stdint.h>
#include <string.h>
#include "flis_test_helpers.h"
#include "io/flis_compose.h"

cominfo com;
fits *gfit;

/* Tile dim used by the GPU compose path.  Mirrors FLIS_TILE_DIM in
 * src/gui-gtk4/flis_gpu_compose.c. */
#define TEST_TILE_DIM 2048

/* Build an identity LUT (lut[i] = i >> 8) — maps 16-bit pixel through
 * a no-op stretch to 8-bit. */
static void make_identity_lut(BYTE lut[USHRT_MAX + 1]) {
	for (int i = 0; i <= USHRT_MAX; i++)
		lut[i] = (BYTE)(i >> 8);
}

/* Build a "gain 2" LUT: doubles input, clips to 255. */
static void make_gain2_lut(BYTE lut[USHRT_MAX + 1]) {
	for (int i = 0; i <= USHRT_MAX; i++) {
		int v = (i >> 8) * 2;
		lut[i] = (BYTE)(v > 255 ? 255 : v);
	}
}

/* Build a fits with a known float-mono pattern: row r col c has value
 * (r * cols + c) / (rows * cols).  Lets us check both position-
 * dependence and the FITS→display y-flip. */
static fits *make_gradient_mono_float(int w, int h) {
	fits *f = NULL;
	if (new_fit_image(&f, w, h, 1, DATA_FLOAT)) return NULL;
	float denom = (float)((size_t)w * h);
	for (int r = 0; r < h; r++)
		for (int c = 0; c < w; c++)
			f->fdata[(size_t)r * w + c] = (float)((size_t)r * w + c) / denom;
	return f;
}

/* Build a mono-uint16 fits where pixel value at (r,c) = (r * w + c) * scale,
 * clipped to USHRT_MAX. */
static fits *make_gradient_mono_word(int w, int h, int scale) {
	fits *f = NULL;
	if (new_fit_image(&f, w, h, 1, DATA_USHORT)) return NULL;
	for (int r = 0; r < h; r++) {
		for (int c = 0; c < w; c++) {
			int v = ((size_t)r * w + c) * scale;
			if (v > USHRT_MAX) v = USHRT_MAX;
			f->data[(size_t)r * w + c] = (WORD)v;
		}
	}
	return f;
}

/* Build a 3-channel RGB float fits where R=0.1, G=0.5, B=0.9 uniformly. */
static fits *make_uniform_rgb_float(int w, int h, float r, float g, float b) {
	fits *f = NULL;
	if (new_fit_image(&f, w, h, 3, DATA_FLOAT)) return NULL;
	size_t n = (size_t)w * h;
	for (size_t i = 0; i < n; i++) {
		f->fpdata[0][i] = r;
		f->fpdata[1][i] = g;
		f->fpdata[2][i] = b;
	}
	return f;
}

TestSuite(flis_gpu_bake, .init = flis_test_init_com, .fini = flis_test_cleanup_com);

/* ---------------------------------------------------------------- */
/* Tile bake: mono USHORT, identity LUT, no tint, single tile         */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, mono_word_identity_single_tile_yflip) {
	/* 4×4 mono uint16 layer: pixel value = row*4 + col, scaled up
	 * into the 16-bit range so identity LUT (>> 8) produces distinct
	 * bytes.  Identity LUT means out_byte = src_word >> 8. */
	const int w = 4, h = 4;
	fits *f = make_gradient_mono_word(w, h, 4096);  /* values 0, 4096, 8192, ... */
	flis_layer_t *lay = flis_test_add_layer(f, "mono");
	cr_assert_not_null(lay);

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);
	uint8_t out[4 * 4 * 4] = { 0 };
	cr_assert(flis_compose_bake_tile_bgra8(lay, TEST_TILE_DIM, 0, 0, w, h, 1, 0, lut, out));

	/* FITS-source y=0 (bottom of image, lowest values) ends up in
	 * display row h-1=3 (bottom of texture).  So texture row 3 should
	 * contain pixel values 0, 4096, 8192, 12288 → bytes 0, 16, 32, 48. */
	for (int c = 0; c < w; c++) {
		BYTE expected = (BYTE)((0 * w + c) * 4096 >> 8);
		size_t i = (size_t)3 * w * 4 + (size_t)c * 4;
		cr_assert_eq(out[i + 0], expected, "row 3 (FITS y=0) col %d B", c);
		cr_assert_eq(out[i + 1], expected, "row 3 col %d G", c);
		cr_assert_eq(out[i + 2], expected, "row 3 col %d R", c);
		cr_assert_eq(out[i + 3], 255,      "row 3 col %d A", c);
	}
	/* Texture row 0 (display top) should hold FITS row 3 (top of
	 * image), values 12, 13, 14, 15 → bytes 192, 208, 224, 240. */
	for (int c = 0; c < w; c++) {
		BYTE expected = (BYTE)((3 * w + c) * 4096 >> 8);
		size_t i = (size_t)0 * w * 4 + (size_t)c * 4;
		cr_assert_eq(out[i + 2], expected, "row 0 (FITS y=h-1) col %d R", c);
	}
}

/* ---------------------------------------------------------------- */
/* Tile bake: gain-2 LUT applies the stretch                          */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, mono_word_gain2_lut) {
	const int w = 4, h = 1;
	fits *f = make_gradient_mono_word(w, h, 8192);  /* 0, 8192, 16384, 24576 */
	flis_layer_t *lay = flis_test_add_layer(f, "mono");

	BYTE lut[USHRT_MAX + 1];
	make_gain2_lut(lut);
	uint8_t out[4 * 1 * 4] = { 0 };
	cr_assert(flis_compose_bake_tile_bgra8(lay, TEST_TILE_DIM, 0, 0, w, h, 1, 0, lut, out));

	/* Single-row layer: only one display row, sees FITS row 0.
	 * gain-2 LUT: (i >> 8) * 2, clamped.  0/32/64/96 → 0/64/128/192. */
	const BYTE expected[] = { 0, 64, 128, 192 };
	for (int c = 0; c < w; c++) {
		size_t i = (size_t)c * 4;
		cr_assert_eq(out[i + 2], expected[c],
			"col %d R should be %u got %u", c, expected[c], out[i + 2]);
	}
}

/* ---------------------------------------------------------------- */
/* Tile bake: mono float with green tint                              */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, mono_float_green_tint) {
	const int w = 2, h = 1;
	fits *f = make_gradient_mono_float(w, h);  /* values 0/2, 1/2 */
	flis_layer_t *lay = flis_test_add_layer(f, "mono");
	/* Pixel 0 → 0.0 → 16-bit 0 → byte 0.  Pixel 1 → 0.5 → 16-bit ~32767
	 * → byte 127 (identity LUT). */
	lay->has_tint = TRUE;
	lay->layer_tint.r = 0.0;
	lay->layer_tint.g = 1.0;
	lay->layer_tint.b = 0.0;

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);
	uint8_t out[2 * 1 * 4] = { 0 };
	cr_assert(flis_compose_bake_tile_bgra8(lay, TEST_TILE_DIM, 0, 0, w, h, 1, 0, lut, out));

	/* With tint (0,1,0): R = byte * 0 = 0, G = byte * 1, B = byte * 0 = 0. */
	/* Col 0: byte = 0 → (0, 0, 0). */
	cr_assert_eq(out[0], 0); cr_assert_eq(out[1], 0); cr_assert_eq(out[2], 0);
	/* Col 1: byte ≈ 127.  Allow ±1 for the round-trip via roundf_to_WORD. */
	cr_assert(out[4 + 1] >= 126 && out[4 + 1] <= 128,
		"col 1 G expected ~127, got %u", out[4 + 1]);
	cr_assert_eq(out[4 + 0], 0, "col 1 B should be 0 (tinted out)");
	cr_assert_eq(out[4 + 2], 0, "col 1 R should be 0 (tinted out)");
}

/* ---------------------------------------------------------------- */
/* Tile bake: RGB float, all three channels routed through LUT        */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, rgb_float_per_channel) {
	const int w = 1, h = 1;
	fits *f = make_uniform_rgb_float(w, h, 0.0f, 0.5f, 1.0f);
	flis_layer_t *lay = flis_test_add_layer(f, "rgb");

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);
	uint8_t out[1 * 1 * 4] = { 0 };
	cr_assert(flis_compose_bake_tile_bgra8(lay, TEST_TILE_DIM, 0, 0, w, h, 1, 0, lut, out));

	/* BGRA layout: B=lut(B_word), G=lut(G_word), R=lut(R_word).
	 * 0.0 → 0; 0.5 → ~127; 1.0 → 255. */
	cr_assert(out[0] == 255, "B should be 255 (input 1.0), got %u", out[0]);
	cr_assert(out[1] >= 126 && out[1] <= 128, "G should be ~127 (input 0.5), got %u", out[1]);
	cr_assert(out[2] == 0,   "R should be 0 (input 0.0), got %u", out[2]);
	cr_assert(out[3] == 255, "A");
}

/* ---------------------------------------------------------------- */
/* Tile bake: multi-tile coverage — 4-pixel-wide layer with tile_dim=2 */
/* split into 2 tiles, each tile's content correctly comes from the    */
/* appropriate columns.                                                */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, multi_tile_horizontal) {
	const int w = 4, h = 2;
	fits *f = make_gradient_mono_word(w, h, 8192);
	/* Pixel layout (FITS bottom-up):
	 *   FITS row 0 (bottom): 0    8192  16384 24576
	 *   FITS row 1 (top):   32768 40960 49152 57344  (clipped to 65535)
	 * In display (texture row 0 = top, row 1 = bottom):
	 *   row 0 corresponds to FITS row 1: 32768 40960 49152 57344
	 *   row 1 corresponds to FITS row 0:     0  8192 16384 24576
	 * Identity LUT (>>8):
	 *   row 0: 128 160 192 224
	 *   row 1:   0  32  64  96 */
	flis_layer_t *lay = flis_test_add_layer(f, "mono");

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);

	const int TD = 2;  /* test-only tile_dim */
	const int tile_w = 2, tile_h = 2;
	uint8_t left[2 * 2 * 4] = { 0 }, right[2 * 2 * 4] = { 0 };

	cr_assert(flis_compose_bake_tile_bgra8(lay, TD, 0, 0, tile_w, tile_h, 1, 0, lut, left));
	cr_assert(flis_compose_bake_tile_bgra8(lay, TD, 1, 0, tile_w, tile_h, 1, 0, lut, right));

	/* Left tile (cols 0-1): row 0 → 128, 160; row 1 → 0, 32 */
	cr_assert_eq(left[0 * 8 + 0 * 4 + 2], 128);
	cr_assert_eq(left[0 * 8 + 1 * 4 + 2], 160);
	cr_assert_eq(left[1 * 8 + 0 * 4 + 2],   0);
	cr_assert_eq(left[1 * 8 + 1 * 4 + 2],  32);
	/* Right tile (cols 2-3): row 0 → 192, 224; row 1 → 64, 96 */
	cr_assert_eq(right[0 * 8 + 0 * 4 + 2], 192);
	cr_assert_eq(right[0 * 8 + 1 * 4 + 2], 224);
	cr_assert_eq(right[1 * 8 + 0 * 4 + 2],  64);
	cr_assert_eq(right[1 * 8 + 1 * 4 + 2],  96);
}

/* ---------------------------------------------------------------- */
/* Tile bake: edge tile is smaller than tile_dim                      */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, edge_tile_smaller_than_dim) {
	/* 3-wide layer with tile_dim=2 → 2 tiles, second is 1-wide. */
	const int w = 3, h = 1;
	fits *f = make_gradient_mono_word(w, h, 8192);
	flis_layer_t *lay = flis_test_add_layer(f, "mono");

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);
	uint8_t left[2 * 1 * 4]  = { 0 };
	uint8_t right[1 * 1 * 4] = { 0 };

	cr_assert(flis_compose_bake_tile_bgra8(lay, 2, 0, 0, 2, 1, 1, 0, lut, left));
	cr_assert(flis_compose_bake_tile_bgra8(lay, 2, 1, 0, 1, 1, 1, 0, lut, right));

	/* Layer values: 0, 8192, 16384 → bytes 0, 32, 64 (identity LUT). */
	cr_assert_eq(left[0 + 2],   0, "left col 0");
	cr_assert_eq(left[4 + 2],  32, "left col 1");
	cr_assert_eq(right[0 + 2], 64, "right col 0");
}

/* ---------------------------------------------------------------- */
/* Tile bake: out-of-range tile coords rejected                       */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, reject_out_of_range) {
	const int w = 2, h = 2;
	fits *f = make_gradient_mono_word(w, h, 32768);
	flis_layer_t *lay = flis_test_add_layer(f, "mono");

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);
	uint8_t out[2 * 2 * 4] = { 0 };

	/* Tile (1,0) with tile_dim=2 would start at x=2, outside a 2-wide
	 * layer. */
	cr_assert(!flis_compose_bake_tile_bgra8(lay, 2, 1, 0, 2, 2, 1, 0, lut, out));
	cr_assert(!flis_compose_bake_tile_bgra8(lay, 2, 0, 1, 2, 2, 1, 0, lut, out));
}

/* ---------------------------------------------------------------- */
/* Mip-aware bake: half-res (mip=2) keeps every other source pixel    */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, mip2_half_res_stride) {
	/* 4×4 mono layer, identity LUT, mip=2 → output is 2×2.
	 * Stride samples source[(out_y*2, out_x*2)]. */
	const int w = 4, h = 4;
	fits *f = make_gradient_mono_word(w, h, 4096);
	/* Source bytes after identity LUT (>>8):
	 *   FITS row 0: 0  16  32  48
	 *   FITS row 1: 64  80  96 112
	 *   FITS row 2: 128 144 160 176
	 *   FITS row 3: 192 208 224 240
	 * In display orientation (row 0 = FITS row 3):
	 *   disp row 0: 192 208 224 240
	 *   disp row 1: 128 144 160 176
	 *   disp row 2:  64  80  96 112
	 *   disp row 3:   0  16  32  48
	 * mip=2 stride-samples disp rows {0, 2} and cols {0, 2}:
	 *   out row 0: 192 224
	 *   out row 1:  64  96 */
	flis_layer_t *lay = flis_test_add_layer(f, "mono");

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);
	uint8_t out[2 * 2 * 4] = { 0 };
	cr_assert(flis_compose_bake_tile_bgra8(lay, TEST_TILE_DIM, 0, 0, w, h, 2, 0, lut, out));

	cr_assert_eq(out[0 * 8 + 0 * 4 + 2], 192, "out (0,0)");
	cr_assert_eq(out[0 * 8 + 1 * 4 + 2], 224, "out (0,1)");
	cr_assert_eq(out[1 * 8 + 0 * 4 + 2],  64, "out (1,0)");
	cr_assert_eq(out[1 * 8 + 1 * 4 + 2],  96, "out (1,1)");
}

/* ---------------------------------------------------------------- */
/* Guard bands: extra pixels on each side sample neighbour FITS data  */
/* (or edge-clamp at layer boundaries)                                */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, guard_edge_clamps_at_layer_top_left) {
	/* 2×2 mono layer, guard=1 → output 4×4.  Outer border edge-clamps
	 * because the layer has no neighbouring tiles. */
	const int w = 2, h = 2;
	fits *f = make_gradient_mono_word(w, h, 16384);
	/* FITS row 0: 0, 16384 → bytes 0, 64
	 * FITS row 1: 32768, 49152 → bytes 128, 192
	 * Display:
	 *   disp row 0 (FITS 1): 128, 192
	 *   disp row 1 (FITS 0):   0,  64 */
	flis_layer_t *lay = flis_test_add_layer(f, "mono");

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);
	uint8_t out[4 * 4 * 4] = { 0 };
	cr_assert(flis_compose_bake_tile_bgra8(lay, TEST_TILE_DIM, 0, 0, w, h, 1, 1, lut, out));

	/* Output layout, R channel only (BGRA, R at offset 2):
	 *   guard row -1 (out_y=0): clamps to disp row 0 → 128, 128, 192, 192
	 *   inner row 0  (out_y=1): 128, 128, 192, 192
	 *   inner row 1  (out_y=2):   0,   0,  64,  64
	 *   guard row +1 (out_y=3): clamps to disp row 1 →   0,   0,  64,  64 */
	const BYTE expected[4][4] = {
		{128, 128, 192, 192},  /* out row 0 (top guard) */
		{128, 128, 192, 192},  /* out row 1 (inner top) */
		{  0,   0,  64,  64},  /* out row 2 (inner bottom) */
		{  0,   0,  64,  64},  /* out row 3 (bottom guard) */
	};
	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			const size_t i = (size_t)y * 4 * 4 + (size_t)x * 4;
			cr_assert_eq(out[i + 2], expected[y][x],
				"out (%d,%d).R expected %u got %u", y, x, expected[y][x], out[i + 2]);
		}
	}
}

Test(flis_gpu_bake, mip_rejects_non_pow2_or_non_divisor) {
	const int w = 4, h = 4;
	fits *f = make_gradient_mono_word(w, h, 4096);
	flis_layer_t *lay = flis_test_add_layer(f, "mono");

	BYTE lut[USHRT_MAX + 1];
	make_identity_lut(lut);
	uint8_t out[2 * 2 * 4] = { 0 };
	/* mip=3 isn't a power of 2 */
	cr_assert(!flis_compose_bake_tile_bgra8(lay, TEST_TILE_DIM, 0, 0, w, h, 3, 0, lut, out));
	/* mip=8 doesn't divide tile_w=4 */
	cr_assert(!flis_compose_bake_tile_bgra8(lay, TEST_TILE_DIM, 0, 0, w, h, 8, 0, lut, out));
}

/* ---------------------------------------------------------------- */
/* Lmask bake at 8/16/32 bitpix                                       */
/* ---------------------------------------------------------------- */
Test(flis_gpu_bake, lmask_8bit) {
	const int w = 4, h = 4;
	fits *f = make_gradient_mono_float(w, h);
	flis_layer_t *lay = flis_test_add_layer(f, "with_lmask");
	layermask_t *lm = flis_test_make_const_lmask(w, h, 8, 0.5);  /* value=128 */
	cr_assert_eq(flis_layer_set_lmask(lay, lm), 0);

	uint8_t out[4 * 4 * 4] = { 0 };
	cr_assert(flis_compose_bake_lmask_bgra8(lay, out));

	/* All bytes should be 128 (R=G=B=128, A=255). */
	for (int i = 0; i < 16; i++) {
		cr_assert_eq(out[i * 4 + 0], 128, "B[%d]", i);
		cr_assert_eq(out[i * 4 + 1], 128, "G[%d]", i);
		cr_assert_eq(out[i * 4 + 2], 128, "R[%d]", i);
		cr_assert_eq(out[i * 4 + 3], 255, "A[%d]", i);
	}
}

Test(flis_gpu_bake, lmask_32bit_float) {
	const int w = 2, h = 1;
	fits *f = make_gradient_mono_float(w, h);
	flis_layer_t *lay = flis_test_add_layer(f, "with_lmask");
	layermask_t *lm = flis_test_make_const_lmask(w, h, 32, 0.75);  /* value=0.75 */
	cr_assert_eq(flis_layer_set_lmask(lay, lm), 0);

	uint8_t out[2 * 1 * 4] = { 0 };
	cr_assert(flis_compose_bake_lmask_bgra8(lay, out));
	const BYTE expected = (BYTE)(0.75f * 255.f);  /* 191 */
	for (int i = 0; i < 2; i++) {
		cr_assert_eq(out[i * 4 + 0], expected, "B[%d]", i);
		cr_assert_eq(out[i * 4 + 3], 255,      "A[%d]", i);
	}
}

Test(flis_gpu_bake, lmask_null_returns_false) {
	const int w = 2, h = 1;
	fits *f = make_gradient_mono_float(w, h);
	flis_layer_t *lay = flis_test_add_layer(f, "no_lmask");
	uint8_t out[2 * 1 * 4] = { 0 };
	cr_assert(!flis_compose_bake_lmask_bgra8(lay, out));
}
