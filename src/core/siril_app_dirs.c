/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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

#include "core/siril.h"
#include "core/proto.h"

#include "siril_app_dirs.h"

static GUserDirectory sdir[] = { G_USER_DIRECTORY_PICTURES,
		G_USER_DIRECTORY_DOCUMENTS };

static const gchar *siril_share_dir = NULL;
static const gchar *siril_config_dir = NULL;
static const gchar *siril_startup_dir = NULL;
static const gchar *siril_locale_dir = NULL;
static const gchar *siril_scripts_repo_dir = NULL;
static const gchar *siril_spcc_repo_dir = NULL;
static const gchar *siril_user_data_dir = NULL;

/* To set the data dir we are looking for the glade file */
static void search_for_data_dir() {
	/* First we are looking for in the package_data_dir */
#ifdef _WIN32 // in the case where the app is started with double click on seq file
	gchar *execname = g_win32_get_package_installation_directory_of_module(NULL);
	gchar *path = g_build_filename(execname, "share", PACKAGE, NULL);
	g_free(execname);
	if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
		siril_share_dir = g_strdup(path);
	}
	g_free(path);
#elif defined(ENABLE_RELOCATABLE_RESOURCES) && defined(OS_OSX)
	const gchar *relocated_path = g_getenv("SIRIL_RELOCATED_RES_DIR");
	if (relocated_path != NULL) {
		gchar *path = g_build_filename(relocated_path, "share", PACKAGE, NULL);
		if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
			siril_share_dir = g_strdup(path);
		}
		g_free(path);
	}
#elif defined(ENABLE_RELOCATABLE_RESOURCES) && defined(__linux__) // appimage
	/* For AppImage */
	const gchar *relocated_path = g_getenv("APPDIR");
	if (relocated_path != NULL) {
		siril_share_dir = g_strdup(g_build_filename(relocated_path, "usr", "share", PACKAGE, NULL));
	}
#else
	gchar *path = g_build_filename(PACKAGE_DATA_DIR, NULL);
	if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
		siril_share_dir = g_strdup(path);
	}
	g_free(path);
#endif
}

static void search_for_config_dir() {
	siril_config_dir = g_get_user_config_dir();
}

/** This function tries to set a startup directory. It first looks at the "Pictures" directory,
 *  then if it does not exist, the "Document" one, Finally, if it fails on some UNIX systems
 *  the dir is set to the home directory.
 *  @return a working directory path if success, NULL if error
 */

static void search_for_startup_dir() {
	const gchar *dir = NULL;
	gint i = 0;
	size_t size;

	size = sizeof(sdir) / sizeof(GUserDirectory);

	while (!dir && i < size) {
		dir = g_get_user_special_dir(sdir[i]);
		i++;
	}
	/* Not every platform has a directory for these logical id */
	if (!dir) {
		dir = g_get_home_dir();
	}
	if (dir)
		siril_startup_dir = g_strdup(dir);
}

/**
 * This function search for the locale dir
 * @return the locale dir
 */
static void search_for_locale_dir() {
#ifdef _WIN32
	gchar *win32_dir;

	win32_dir = g_win32_get_package_installation_directory_of_module(NULL);
	gchar *locale_dir = g_build_filename(win32_dir, "share", "locale", NULL);

	g_free(win32_dir);

	siril_locale_dir = locale_dir;
#elif defined(ENABLE_RELOCATABLE_RESOURCES) && defined(OS_OSX)
	const gchar *relocated_path = g_getenv("SIRIL_RELOCATED_RES_DIR");
	if (relocated_path != NULL) {
		siril_locale_dir = g_build_filename(relocated_path, "share", "locale", NULL);
	} else {
		gchar *path = g_build_filename(LOCALEDIR, NULL);
		if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
			siril_locale_dir = g_strdup(path);
		} else {
			g_warning("Locale directory %s not found in OSX, using fallback", path);
			siril_locale_dir = g_strdup("/usr/share/locale");
		}
		g_free(path);
	}
#elif defined(ENABLE_RELOCATABLE_RESOURCES) && defined(__linux__) // appimage
	const gchar *relocated_path = g_getenv("APPDIR");
	if (relocated_path != NULL) {
		siril_locale_dir = g_build_filename(relocated_path, "usr", "share", "locale", NULL);
	} else {
		g_warning("APPDIR environment variable not set in AppImage");
		gchar *path = g_build_filename(LOCALEDIR, NULL);
		if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
			siril_locale_dir = g_strdup(path);
		} else {
			g_warning("Locale directory %s not found in AppImage, using fallback", path);
			siril_locale_dir = g_strdup("/usr/share/locale");
		}
		g_free(path);
	}
#else
	gchar *path = g_build_filename(LOCALEDIR, NULL);

	if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
		siril_locale_dir = g_strdup(path);
	} else {
		g_warning("Locale directory %s not found, using fallback", path);
		siril_locale_dir = g_strdup("/usr/share/locale");
	}
	g_free(path);
#endif
}

static void search_for_scripts_repo_dir() {
	siril_scripts_repo_dir = g_build_filename(g_get_user_data_dir(), "siril-scripts", NULL);
}

static void search_for_spcc_repo_dir() {
	siril_spcc_repo_dir = g_build_filename(g_get_user_data_dir(), "siril-spcc-database", NULL);
}

static void search_for_user_data_dir() {
	siril_user_data_dir = g_build_filename(g_get_user_data_dir(), "siril", NULL);
}

/** Public functions **/

const gchar* siril_get_locale_dir() {
	return siril_locale_dir;
}
const gchar* siril_get_startup_dir() {
	return siril_startup_dir;
}

const gchar* siril_get_system_data_dir() {
	return siril_share_dir;
}

const gchar* siril_get_config_dir() {
	return siril_config_dir;
}

void initialize_siril_directories() {
	search_for_data_dir();
	search_for_locale_dir();
	search_for_startup_dir();
	search_for_config_dir();
	search_for_scripts_repo_dir();
	search_for_spcc_repo_dir();
	search_for_user_data_dir();
}

const gchar* siril_get_scripts_repo_path() {
	return siril_scripts_repo_dir;
}

const gchar* siril_get_spcc_repo_path() {
	return siril_spcc_repo_dir;
}

const gchar *siril_get_user_data_dir() {
	return siril_user_data_dir;

}
