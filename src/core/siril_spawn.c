// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "core/siril_spawn.h"

/*
 * Copy non-owned pointer args to a temporary array,
 * prepending the wrapper program which must be used
 * on flatpak.
 */
GPtrArray*
siril_flatpak_argv(gchar** argv) {
	GPtrArray* argv_new = g_ptr_array_new_full(32, NULL);
	g_ptr_array_add(argv_new, (gchar*)"/bin/flatpak-spawn");
	g_ptr_array_add(argv_new, (gchar*)"--host");
	g_ptr_array_add(argv_new, (gchar*)"--");

	for (gchar** arg = argv; *arg != NULL; ++arg) {
		g_ptr_array_add(argv_new, *arg);
	}

	g_ptr_array_add(argv_new, NULL);

	return argv_new;
}

gboolean
siril_spawn_host_sync(
	const gchar* working_directory,
	gchar** argv,
	gchar** envp,
	GSpawnFlags flags,
	GSpawnChildSetupFunc child_setup,
	gpointer user_data,
	gchar** standard_output,
	gchar** standard_error,
	gint* wait_status,
	GError** error
) {
	g_return_val_if_fail(argv != NULL, FALSE);
	g_return_val_if_fail(argv[0] != NULL, FALSE);
	g_return_val_if_fail(flags & G_SPAWN_SEARCH_PATH, FALSE);

#ifdef FLATPAK_ID
	// prepend flatpak-spawn command
	g_autoptr(GPtrArray) argv_new = siril_flatpak_argv(argv);
	argv = (gchar **)argv_new->pdata;

	// use absolute path only
	flags &= ~G_SPAWN_SEARCH_PATH;
#endif

	return g_spawn_sync(
		working_directory,
		argv,
		envp,
		flags,
		child_setup,
		user_data,
		standard_output,
		standard_error,
		wait_status,
		error
	);
}

gboolean
siril_spawn_host_async_with_pipes(
	const gchar* working_directory,
	gchar** argv,
	gchar** envp,
	GSpawnFlags flags,
	GSpawnChildSetupFunc child_setup,
	gpointer user_data,
	GPid* child_pid,
	gint* standard_input,
	gint* standard_output,
	gint* standard_error,
	GError** error
) {
	g_return_val_if_fail(argv != NULL, FALSE);
	g_return_val_if_fail(argv[0] != NULL, FALSE);
	g_return_val_if_fail(flags & G_SPAWN_SEARCH_PATH, FALSE);

#ifdef FLATPAK_ID
	// prepend flatpak-spawn command
	g_autoptr(GPtrArray) argv_new = siril_flatpak_argv(argv);
	argv = (gchar**)argv_new->pdata;

	// use absolute path only
	flags &= ~G_SPAWN_SEARCH_PATH;
#endif

	return g_spawn_async_with_pipes(
		working_directory,
		argv,
		envp,
		flags,
		child_setup,
		user_data,
		child_pid,
		standard_input,
		standard_output,
		standard_error,
		error
	);
}
