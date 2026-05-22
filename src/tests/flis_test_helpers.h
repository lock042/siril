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
 * flis_test_helpers.h — small inline helpers shared across FLIS test files.
 *
 * Each test file declares `cominfo com` and `fits *gfit` itself (link-time
 * required globals), then includes this header for the shared boilerplate
 * for building fits structs and FLIS layer stacks.
 */

#ifndef FLIS_TEST_HELPERS_H
#define FLIS_TEST_HELPERS_H

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <glib/gstdio.h>

#include "core/siril.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/flis_compose.h"

/* Minimal com initialisation: set up enough state so single_image_is_loaded
 * works and com.uniq is fresh.  Tests that exercise actual FLIS code paths
 * call this once per test. */
static inline void flis_test_init_com(void) {
	memset(&com, 0, sizeof(com));
	com.uniq = calloc(1, sizeof(single));
	com.uniq->next_item_id = 1;
	com.uniq->active_layer = 0;
	com.uniq->chans = 0;  /* set when first layer added */
	com.max_thread = 1;
	com.pref.swap_dir = g_strdup(g_get_tmp_dir());
	com.script = TRUE;  /* tests run headless; ensures undo paths short-circuit */
}

/* Tear down everything allocated by flis_test_init_com plus any layers /
 * groups that were added. */
static inline void flis_test_cleanup_com(void) {
	if (com.uniq) {
		if (com.uniq->layers) flis_free_layers(com.uniq);
		if (com.uniq->groups) flis_free_groups(com.uniq);
		g_free(com.uniq->filename);
		g_free(com.uniq->comment);
		free(com.uniq);
		com.uniq = NULL;
	}
	if (com.pref.swap_dir) {
		g_free(com.pref.swap_dir);
		com.pref.swap_dir = NULL;
	}
	gfit = NULL;
}

/* Build a fresh constant-colour float fits.  For mono pass chans=1; for RGB
 * pass chans=3 and r/g/b values via the array form. */
static inline fits *flis_test_make_mono_fits(int rx, int ry, float value) {
	fits *f = NULL;
	if (new_fit_image(&f, rx, ry, 1, DATA_FLOAT)) return NULL;
	size_t n = (size_t)rx * ry;
	for (size_t i = 0; i < n; i++) f->fdata[i] = value;
	return f;
}

static inline fits *flis_test_make_rgb_fits(int rx, int ry,
                                            float r, float g, float b) {
	fits *f = NULL;
	if (new_fit_image(&f, rx, ry, 3, DATA_FLOAT)) return NULL;
	size_t n = (size_t)rx * ry;
	float channels[3] = { r, g, b };
	for (int c = 0; c < 3; c++)
		for (size_t i = 0; i < n; i++)
			f->fpdata[c][i] = channels[c];
	return f;
}

/* Allocate-and-fill a layermask of constant value.  bitpix must be 8, 16,
 * or 32; value is interpreted in the natural range for that bit depth. */
static inline layermask_t *flis_test_make_const_lmask(size_t w, size_t h,
                                                     guint8 bitpix, double value) {
	layermask_t *lm = calloc(1, sizeof(layermask_t));
	if (!lm) return NULL;
	lm->w = w; lm->h = h; lm->bitpix = bitpix;
	size_t n = w * h;
	switch (bitpix) {
		case 8: {
			lm->data = malloc(n);
			uint8_t v = (uint8_t)(value * 255.0 + 0.5);
			memset(lm->data, v, n);
			break;
		}
		case 16: {
			lm->data = malloc(n * sizeof(uint16_t));
			uint16_t v = (uint16_t)(value * 65535.0 + 0.5);
			uint16_t *p = (uint16_t *)lm->data;
			for (size_t i = 0; i < n; i++) p[i] = v;
			break;
		}
		case 32: {
			lm->data = malloc(n * sizeof(float));
			float *p = (float *)lm->data;
			float v = (float)value;
			for (size_t i = 0; i < n; i++) p[i] = v;
			break;
		}
		default:
			free(lm);
			return NULL;
	}
	return lm;
}

/* Add a layer to com.uniq backed by the given fits.  Returns the layer or
 * NULL.  Ownership of fit is transferred to the layer. */
static inline flis_layer_t *flis_test_add_layer(fits *f, const char *name) {
	return flis_layer_add(f, name);
}

/* Approximate equality for float compares; tolerance is per-channel for
 * BGRA8/float-RGB results. */
static inline gboolean flis_test_approx_eq(float a, float b, float eps) {
	return fabsf(a - b) <= eps;
}

#endif /* FLIS_TEST_HELPERS_H */
