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
#ifndef SRC_CORE_SIRIL_SPAWN_H_
#define SRC_CORE_SIRIL_SPAWN_H_

#include <glib.h>

/**
 * @brief Wrapper around g_spawn_sync
 *
 * This method permits Siril to call external executables that
 * reside on the host. To harmonize behavior among non-flatpak
 * and flatpak platforms, this method requires the caller to
 * set the `G_SPAWN_SEARCH_PATH` flag. If a path search is not
 * desired, the caller should provide an absolute executable path.
 *
 * - On non-flatpak platforms, this function is aliased to
 *   `g_spawn_sync()` and behaves identically—except that
 *   `G_SPAWN_SEARCH_PATH` must be set.
 *
 * - On flatpak, this function requests that the flatpak runtime
 *   start the specified executable (by `argv`) outside of the
 *   flatpak sandbox, directly on the host.
 *
 * The arguments and behavior are otherwise as per
 * \ref g_spawn_sync()
 */
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
);

/**
 * @brief Wrapper around g_spawn_async_with_pipes
 *
 * This method permits Siril to call external executables that
 * reside on the host. To harmonize behavior among non-flatpak
 * and flatpak platforms, this method requires the caller to
 * set the `G_SPAWN_SEARCH_PATH` flag. If a path search is not
 * desired, the caller should provide an absolute executable path.
 *
 * - On non-flatpak platforms, this function is aliased to
 *   `g_spawn_async_with_pipes()` and behaves identically—except
 *   that `G_SPAWN_SEARCH_PATH` must be set.
 *
 * - On flatpak, this function requests that the flatpak runtime
 *   start the specified executable (by `argv`) outside of the
 *   flatpak sandbox, directly on the host.
 *
 * The arguments and behavior are otherwise as per
 * \ref g_spawn_async_with_pipes()
 */
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
);

/**
 * @def siril_spawn_check_wait_status
 * @brief Checks waitpid() status
 *
 * On GLIB >= 2.70, this is an alias for
 * `g_spawn_check_wait_status()`. On previous versions of GLIB,
 * it is an alias to a now-deprecated function with the same
 * name.
 *
 * This method checks the result of a `waitpid()` or
 * `g_spawn_sync()` to determine if the process exited
 * successfully.
 */

#if GLIB_CHECK_VERSION(2, 70, 0)
#define siril_spawn_check_wait_status g_spawn_check_wait_status
#else
#define siril_spawn_check_wait_status g_spawn_check_exit_status
#endif

#endif /* SRC_CORE_SIRIL_SPAWN_H_ */
