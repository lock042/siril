// SPDX-License-Identifier: GPL-3.0-or-later

#include <criterion/criterion.h>

#include "core/siril_spawn.h"

extern GPtrArray* siril_flatpak_argv(gchar** argv);

void test_siril_flatpak_argv_null() {
	gchar* nullargs[] = {NULL};

	g_autoptr(GPtrArray) argv_null = siril_flatpak_argv(nullargs);
	cr_assert_eq(4, argv_null->len);

	const gchar* const* newargs = (const gchar* const*)argv_null->pdata;
	cr_expect_str_eq("/bin/flatpak-spawn", (gchar*)newargs[0]);
	cr_expect_str_eq("--host", (gchar*)newargs[1]);
	cr_expect_str_eq("--", (gchar*)newargs[2]);
	cr_expect_eq(NULL, newargs[3]);
}

void test_siril_flatpak_argv_populated() {
	gchar* args[] = { "echo", "hello", NULL };

	g_autoptr(GPtrArray) argv_populated = siril_flatpak_argv(args);
	cr_assert_eq(6, argv_populated->len);

	const gchar* const* newargs = (const gchar* const*)argv_populated->pdata;
	cr_expect_str_eq("/bin/flatpak-spawn", newargs[0]);
	cr_expect_str_eq("--host", newargs[1]);
	cr_expect_str_eq("--", newargs[2]);
	cr_expect_str_eq("echo", newargs[3]);
	cr_expect_str_eq("hello", newargs[4]);
	cr_expect_eq(NULL, newargs[5]);
}

Test(siril, flatpak_argv_null) { test_siril_flatpak_argv_null(); }
Test(siril, flatpak_argv_populated) { test_siril_flatpak_argv_populated(); }
