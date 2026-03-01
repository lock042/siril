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
 * along with Siril. If not, see <http://www.gnu.org/licenses/ >.
 */

#include <criterion/criterion.h>
#include <math.h>
#include "core/siril.h"
#include "core/arithm.h"
#include "io/image_format_fits.h"

cominfo com;	// the core data struct
guiinfo gui;	// the gui data struct
fits *gfit;	// currently loaded image (now a pointer)

// Setup function to allocate global gfit as in the main program
static void setup(void) {
	gfit = calloc(1, sizeof(fits));
	cr_assert(gfit != NULL, "Failed to allocate global gfit");
}

// Teardown function to free global gfit after tests
static void teardown(void) {
	if (gfit) {
		clearfits(gfit);
		free(gfit);
		gfit = NULL;
	}
}

// Define test suite with setup/teardown
TestSuite(arithmetics, .init = setup, .fini = teardown);

static void set_ushort_data(fits *fit, WORD *data, int length) {
	memcpy(fit->data, data, length * sizeof(WORD));
}

static void set_float_data(fits *fit, float *data, int length) {
	memcpy(fit->fdata, data, length * sizeof(float));
}

static WORD *alloc_data(const WORD *from, int length) {
	WORD *data = malloc(5 * sizeof(WORD));
	memcpy(data, from, length * sizeof(WORD));
	return data;
}

static float *alloc_fdata(const float *from, int length) {
	float *data = malloc(5 * sizeof(float));
	memcpy(data, from, length * sizeof(float));
	return data;
}

/* image a is ushort, */
void test_ushort() {
	fits *a = NULL;
	int size = 5;

	WORD origa[] = { 0, 1, 2, 1000, 65535 };

	WORD *dataa = alloc_data(origa, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);

	int retval = soper(a, 2.f, OPER_MUL, FALSE);
	cr_assert(!retval, "soper MUL failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 2);
	cr_expect_eq(a->data[2], 4);
	cr_expect_eq(a->data[3], 2000);
	cr_expect_eq(a->data[4], 65535);

	set_ushort_data(a, origa, size);
	retval = soper(a, 2.f, OPER_DIV, FALSE);
	cr_assert(!retval, "soper DIV failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 1);
	cr_expect_eq(a->data[2], 1);
	cr_expect_eq(a->data[3], 500);
	cr_expect_eq(a->data[4], 32768);

	set_ushort_data(a, origa, size);
	retval = soper(a, 2.f * INV_USHRT_MAX_SINGLE, OPER_ADD, FALSE);
	cr_assert(!retval, "soper ADD failed");
	cr_expect_eq(a->data[0], 2);
	cr_expect_eq(a->data[1], 3);
	cr_expect_eq(a->data[2], 4);
	cr_expect_eq(a->data[3], 1002);
	cr_expect_eq(a->data[4], 65535);

	set_ushort_data(a, origa, size);
	retval = soper(a, 2.f * INV_USHRT_MAX_SINGLE, OPER_SUB, FALSE);
	cr_assert(!retval, "soper SUB failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 0);
	cr_expect_eq(a->data[2], 0);
	cr_expect_eq(a->data[3], 998);
	cr_expect_eq(a->data[4], 65533);
}

void test_ushort_force_float() {
	fits *a = NULL;
	int size = 5;

	WORD origa[] = { 0, 1, 2, 1000, 65535 };
	/*
	 * 0     ---> 0.f
	 * 1     ---> 0.0000152590
	 * 2     ---> 0.0000305180
	 * 1000  ---> 0.0152590219
	 */

	WORD *dataa = alloc_data(origa, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);

	int retval = soper(a, 2.f, OPER_MUL, TRUE);
	cr_assert(!retval, "soper MUL failed");
	cr_expect_float_eq(a->fdata[0], 0.0f, 1e-7, "expected 0.0, got %.7f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], 0.0000305f, 1e-7, "expected 0.0000305, got %.7f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 0.0000610f, 1e-7, "expected 0.0000610, got %.7f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 0.0305180f, 1e-7, "expected 0.0305180, got %.7f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 2.0f, 1e-7, "expected 1.0, got %f", a->fdata[4]);

	clearfits(a);
	dataa = alloc_data(origa, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);

	retval = soper(a, 2.f, OPER_DIV, TRUE);
	cr_assert(!retval, "soper MUL failed");
	cr_expect_float_eq(a->fdata[0], 0.0f, 1e-7, "expected 0.0, got %.7f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], 0.0000076f, 1e-7, "expected 0.0000076, got %.7f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 0.0000153f, 1e-7, "expected 0.0000153, got %.7f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 0.0076295f, 1e-7, "expected 0.0076295, got %.7f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 0.5f, 1e-7, "expected 0.5, got %f", a->fdata[4]);

	clearfits(a);
	dataa = alloc_data(origa, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);

	retval = soper(a, 0.1f, OPER_ADD, TRUE);
	cr_assert(!retval, "soper MUL failed");
	cr_expect_float_eq(a->fdata[0], 0.1f, 1e-7, "expected 0.1, got %.7f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], 0.1000153f, 1e-7, "expected 0.1000153, got %.7f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 0.1000305f, 1e-7, "expected 0.1000305, got %.7f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 0.1152590f, 1e-7, "expected 0.1152590, got %.7f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 1.1f, 1e-7, "expected 1.1, got %f", a->fdata[4]);

	clearfits(a);
	dataa = alloc_data(origa, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);

	retval = soper(a, 0.1f, OPER_SUB, TRUE);
	cr_assert(!retval, "soper MUL failed");
	cr_expect_float_eq(a->fdata[0], -0.1f, 1e-7, "expected -0.1, got %.7f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], -0.0999847f, 1e-7, "expected -0.0999847, got %.7f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], -0.0999695f, 1e-7, "expected -0.0999695, got %.7f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], -0.0847410f, 1e-7, "expected -0.084741, got %.7f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 0.9f, 1e-7, "expected 0.9, got %f", a->fdata[4]);

}

void test_float() {
	fits *a = NULL;
	int size = 5;

	float origa[] = { 0.0f, 0.01f, 0.1f, 0.3f, 1.0f };

	float *dataa = alloc_fdata(origa, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_FLOAT, dataa);

	int retval = soper(a, 2.f, OPER_MUL, FALSE);
	cr_assert(!retval, "soper MUL failed");
	cr_expect_float_eq(a->fdata[0], 0.0f, 1e-7, "expected 0.0, got %f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], 0.02f, 1e-7, "expected 0.02, got %f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 0.2f, 1e-7, "expected 0.2, got %f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 0.6f, 1e-7, "expected 0.6, got %f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 2.0f, 1e-7, "expected 1.0, got %f", a->fdata[4]);

	set_float_data(a, origa, size);
	retval = soper(a, 2.f, OPER_DIV, FALSE);
	cr_assert(!retval, "soper DIV failed");
	cr_expect_float_eq(a->fdata[0], 0.0f, 1e-7, "expected 0.0, got %f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], 0.005f, 1e-7, "expected 0.005, got %f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 0.05f, 1e-7, "expected 0.05, got %f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 0.15f, 1e-7, "expected 0.15, got %f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 0.5f, 1e-7, "expected 0.5, got %f", a->fdata[4]);

	set_float_data(a, origa, size);
	retval = soper(a, 0.08, OPER_ADD, FALSE);
	cr_assert(!retval, "soper ADD failed");
	cr_expect_float_eq(a->fdata[0], 0.08f, 1e-7, "expected 0.08, got %f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], 0.09f, 1e-7, "expected 0.09, got %f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 0.18f, 1e-7, "expected 0.18, got %f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 0.38f, 1e-7, "expected 0.38, got %f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 1.08f, 1e-7, "expected 1.08, got %f", a->fdata[4]);

	set_float_data(a, origa, size);
	retval = soper(a, 0.3, OPER_SUB, FALSE);
	cr_assert(!retval, "soper SUB failed");
	cr_expect_float_eq(a->fdata[0], -0.3f, 1e-7, "expected -0.3, got %f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], -0.29f, 1e-7, "expected -0.29, got %f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], -0.2f, 1e-7, "expected -0.2, got %f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 0.0f, 1e-7, "expected 0.0, got %f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 0.7f, 1e-7, "expected 0.7, got %f", a->fdata[4]);
}

Test(arithmetics, ushort_test) { test_ushort(); }
Test(arithmetics, ushort_test_force) { test_ushort_force_float(); }
Test(arithmetics, float_test) { test_float(); }
