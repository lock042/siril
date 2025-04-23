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
#include <glib.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/siril_app_dirs.h"
#include "algos/photometry.h"
#include "algos/star_finder.h"
#include "io/sequence.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "stacking/stacking.h"

#include "initfile.h"

#define STR_INDIR(x) #x
#define STR(x) STR_INDIR(x)
#define GLIB_CONFIG_FILE "config." STR(SIRIL_MAJOR_VERSION) "." STR(SIRIL_MINOR_VERSION) ".ini"
static char *configfiles[] = { GLIB_CONFIG_FILE, "config.1.3.ini", "config.ini" };

static int get_key_data(GKeyFile *kf, struct settings_access *desc) {
	GError *error = NULL;
	gboolean boolval;
	int intval;
	double doubleval;
	gchar *strval = NULL;
	gsize len;
	gchar **strs = NULL;
	switch (desc->type) {
		case STYPE_BOOL:
			boolval = g_key_file_get_boolean(kf, desc->group, desc->key, &error);
			if (error && error->code == G_KEY_FILE_ERROR_INVALID_VALUE) {
				gchar* keystring = g_key_file_get_string(kf, desc->group, desc->key, NULL);
				siril_log_message(_("error in config file for %s.%s: %s (value: %s)\n"),
						desc->group, desc->key, error->message, keystring);
				g_free(keystring);
				return 1;
			}
			*((gboolean*)desc->data) = boolval;
			break;
		case STYPE_INT:
			intval = g_key_file_get_integer(kf, desc->group, desc->key, &error);
			if (error && error->code == G_KEY_FILE_ERROR_INVALID_VALUE) {
				gchar* keystring = g_key_file_get_string(kf, desc->group, desc->key, NULL);
				siril_log_message(_("error in config file for %s.%s: %s (value: %s)\n"),
						desc->group, desc->key, error->message, keystring);
				g_free(keystring);
				return 1;
			}
			if (desc->range_int.min != 0 || desc->range_int.max != 0) {
				if (intval < desc->range_int.min || intval > desc->range_int.max) {
					siril_log_message(_("value %d is out of range [%d, %d] for %s.%s\n"),
							intval, desc->range_int.min, desc->range_int.max,
							desc->group, desc->key);
					return 1;
				}
			}

			*((int*)desc->data) = intval;
			break;
		case STYPE_DOUBLE:
			doubleval = g_key_file_get_double(kf, desc->group, desc->key, &error);
			if (error && error->code == G_KEY_FILE_ERROR_INVALID_VALUE) {
				gchar* keystring = g_key_file_get_string(kf, desc->group, desc->key, NULL);
				siril_log_message(_("error in config file for %s.%s: %s (value: %s)\n"),
						desc->group, desc->key, error->message, keystring);
				g_free(keystring);
				return 1;
			}
			if (desc->range_double.min != 0.0 || desc->range_double.max != 0.0) {
				if (doubleval < desc->range_double.min || doubleval > desc->range_double.max) {
					siril_log_message(_("value %f is out of range [%f, %f] for %s.%s\n"),
							doubleval, desc->range_double.min, desc->range_double.max,
							desc->group, desc->key);
					return 1;
				}
			}

			*((double*)desc->data) = doubleval;
			break;
		case STYPE_STR:
		case STYPE_STRDIR:
			strval = g_key_file_get_string(kf, desc->group, desc->key, NULL);
			if (!strval) {
				siril_log_message(_("unknown error in config file for %s.%s\n"),
						desc->group, desc->key);
				return 1;
			}
			if (strval[0] == '\0') {
				g_free(strval);
				return 1;
			}
			if (desc->type == STYPE_STRDIR && !g_file_test(strval, G_FILE_TEST_IS_DIR)) {
				siril_log_color_message(_("directory `%s' for config key %s.%s doesn't exist, not using it.\n"),
						"salmon", strval, desc->group, desc->key);
				g_free(strval);
				return 1;
			}
			gchar *old_value = *((gchar**)desc->data);
			if (old_value)
				g_free(old_value);
			*((gchar**)desc->data) = strval;
			break;
		case STYPE_STRLIST:
			strs = g_key_file_get_string_list(kf, desc->group, desc->key, &len, NULL);
			if (strs && len > 0) {
				GSList *list = NULL;
				for (gsize i = 0; i < len; i++)
					list = g_slist_prepend(list, strs[i]);
				list = g_slist_reverse(list); // Preserve order
				GSList *old_list = *((GSList**)desc->data);
				if (old_list)
					g_slist_free_full(old_list, g_free);
				*((GSList**)desc->data) = list;
			}
			g_free(strs);
			break;
	}
	return 0;
}

static gchar *get_initfile(const gchar *pathname) {
	for (int i = 0; i < G_N_ELEMENTS(configfiles); i++) {
		gchar *current = g_build_filename(pathname, configfiles[i], NULL);
		if (g_file_test(current, G_FILE_TEST_EXISTS)) {
			if (i > 0)
				siril_log_message(_("Converting previous initfile: %s\n"), current);
			return current;
		}
		g_free(current);
	}
	return NULL;
}

/**
 * Public functions
 */

int read_keyfile(GKeyFile *kf) {
	gsize nb_keys_read = 0;
	gsize nb_groups;
	gchar **groups = g_key_file_get_groups(kf, &nb_groups);
	for (gsize group = 0; group < nb_groups; group++) {
		gsize nb_keys;
		gchar **keys = g_key_file_get_keys(kf, groups[group], &nb_keys, NULL);

		for (gsize key = 0; key < nb_keys; key++) {
			struct settings_access *desc = get_key_settings(groups[group], keys[key]);
			if (!desc) {
				siril_log_message(_("unknown settings variable %s.%s\n"), groups[group], keys[key]);
				continue;
			}
			if (!get_key_data(kf, desc))
				nb_keys_read++;
		}
		g_strfreev(keys);
	}
	g_strfreev(groups);
	siril_debug_print("read %zd keys from key file\n", nb_keys_read);
	return nb_keys_read;
}

int readinitfile(gchar *fname) {
	GKeyFile *kf = g_key_file_new();
	GError *error = NULL;
	if (!g_key_file_load_from_file(kf, fname, G_KEY_FILE_NONE, &error)) {
		if (error != NULL) {
			siril_log_color_message(_("Settings could not be loaded from %s: %s\n"), "red", fname, error->message);
			g_clear_error(&error);
		}
		return 1;
	}
	com.pref.check_update = FALSE;
	if (read_keyfile(kf) == 0)
		siril_log_message(_("Warning: nothing could be read from the settings file\n"));
	g_key_file_free(kf);
	return 0;
}


int checkinitfile() {
	/* com.initfile will contain the path passed with -i if any, NULL else */
	if (com.initfile) {
		siril_log_message(_("Reading configuration file %s\n"), com.initfile);
		return readinitfile(com.initfile);
	}
	set_wisdom_file();
	int retval = 0;

	/* set com.initfile to default location */
	gchar *pathname = g_build_filename(siril_get_config_dir(), PACKAGE, NULL);
	gchar *config_file = g_build_filename(pathname, GLIB_CONFIG_FILE, NULL);
	gchar *existing_config_file = get_initfile(pathname);
	if (!existing_config_file) {
			/* neither files found, create the directory and load defaults */
			initialize_default_settings();

			if (g_mkdir_with_parents(pathname, 0755) == 0) {
				g_fprintf(stdout, "Created config dir %s\n", pathname);
				com.initfile = config_file;
			} else {
				siril_log_message(_("Failed to create config dir %s\n"), pathname);
				g_free(config_file);
				config_file = NULL;
				//config_dir_failed = TRUE;
				if (com.headless)
					siril_log_message(_("Continuing without configuration file\n"));
				else retval = 1;
			}
	}
	else {
		com.initfile = existing_config_file;
		retval = readinitfile(com.initfile);
		g_free(existing_config_file);
		com.initfile = config_file;
	}
	g_free(pathname);
	return retval;
}

int writeinitfile() {
	if (!com.initfile) {
		siril_debug_print("not saving settings, file not available\n");
		return 0;
	}
	siril_debug_print("saving ini file %s\n", com.initfile);

	GKeyFile *kf = g_key_file_new();
	struct settings_access *desc = get_all_settings();
	while (desc->group) {
		gchar *str;
		guint count;
		GSList *list;
		gchar **strs;
		guint i;
		switch (desc->type) {
			case STYPE_BOOL:
				g_key_file_set_boolean(kf, desc->group, desc->key, *((gboolean*)desc->data));
				break;
			case STYPE_INT:
				g_key_file_set_integer(kf, desc->group, desc->key, *((int*)desc->data));
				break;
			case STYPE_DOUBLE:
				g_key_file_set_double(kf, desc->group, desc->key, *((double*)desc->data));
				break;
			case STYPE_STR:
			case STYPE_STRDIR:
				str = *((gchar**)desc->data);
				if (!str) str = "";
				g_key_file_set_string(kf, desc->group, desc->key, str);
				break;
			case STYPE_STRLIST:
				list = *((GSList**)desc->data);
				count = g_slist_length(list);
				strs = malloc((count + 1) * sizeof(gchar *));
				for (i = 0; i < count; i++) {
					strs[i] = list->data;
					list = list->next;
				}
				strs[i] = NULL;
				g_key_file_set_string_list(kf, desc->group, desc->key,
						(const gchar * const*)strs, (gsize)count);
				free(strs);
				break;
		}
		desc++;
	}

	GError *error = NULL;
	if (!g_key_file_save_to_file(kf, com.initfile, &error)) {
		siril_log_color_message(_("Could not save the settings in %s: %s\n"), "salmon", com.initfile, error->message);
		g_free(com.initfile);
		com.initfile = NULL;
		g_key_file_free(kf);
		g_clear_error(&error);
		return 1;
	}
	g_clear_error(&error);
	g_key_file_free(kf);
	return 0;
}
