/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2015 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#ifdef HAVE_LIBCONFIG
#include <libconfig.h>
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

#define GLIB_CONFIG_FILE "config.ini"

#ifdef HAVE_LIBCONFIG
#define LIBCONFIG_FILE "siril.config"
static const char *keywords[] = { "working-directory",
		"debayer-settings", "prepro-settings", "registration-settings",
		"stacking-settings", "astrometry-settings", "photometry-settings",
		"misc-settings", "compression-settings" };

static int readinitfile_libconfig(gchar *path) {
	config_t config;
	const char *dir = NULL;
	GSList *list = NULL;
	gchar *file_path;

	config_init(&config);

#ifdef _WIN32
	/* in the case the filename is given as argument */
	file_path = g_win32_locale_filename_from_utf8(path);
#else
	file_path = path;
#endif

	if (config_read_file(&config, file_path) == CONFIG_FALSE) {
		siril_log_color_message(_("Cannot load initfile: %s\n"), "red", config_error_text(&config));
		config_destroy(&config);
		return 1;
	}
	siril_log_message(_("Loading old configuration file: '%s'\n"), file_path);

	/* Keeping the up-scaled files poses a few problems with sequence
	 * filtering changing and user comprehension, so for now it can only be
	 * enabled by uncommenting the following line. */
	//com.cache_upscaled = TRUE;

	/* Working directory */
	if (config_lookup_string(&config, keywords[WD], &dir)) {
		free(com.wd);
		com.wd = g_strdup(dir);
	}

	/* Debayer setting */
	config_setting_t *debayer_setting = config_lookup(&config, keywords[BAY]);
	if (debayer_setting) {
		config_setting_lookup_bool(debayer_setting, "use_bayer_header", &com.pref.debayer.use_bayer_header);
		config_setting_lookup_int(debayer_setting, "bayer_pattern", &com.pref.debayer.bayer_pattern);
		config_setting_lookup_bool(debayer_setting, "roworder_top_down", &com.pref.debayer.top_down);
		config_setting_lookup_int(debayer_setting, "debayer_algo", (int*)&com.pref.debayer.bayer_inter);
		config_setting_lookup_int(debayer_setting, "x_bayer_offset", &com.pref.debayer.xbayeroff);
		config_setting_lookup_int(debayer_setting, "y_bayer_offset", &com.pref.debayer.ybayeroff);
	}

	/* Preprocessing settings */
	config_setting_t *prepro_setting = config_lookup(&config, keywords[PRE]);
	if (prepro_setting) {
		const char *bias = NULL, *bias_synth = NULL, *dark = NULL, *flat = NULL;

		config_setting_lookup_bool(prepro_setting, "cfa", &com.pref.prepro.cfa);
		config_setting_lookup_bool(prepro_setting, "equalize_cfa", &com.pref.prepro.equalize_cfa);
		config_setting_lookup_bool(prepro_setting, "fix_xtrans", &com.pref.prepro.fix_xtrans);

		config_setting_lookup_string(prepro_setting, "bias_lib", &bias);
		com.pref.prepro.bias_lib = g_strdup(bias);
		config_setting_lookup_bool(prepro_setting, "use_bias_lib", &com.pref.prepro.use_bias_lib);


		config_setting_lookup_string(prepro_setting, "bias_synth", &bias_synth);
		com.pref.prepro.bias_synth = g_strdup(bias_synth);
		config_setting_lookup_bool(prepro_setting, "use_bias_synth", &com.pref.prepro.use_bias_synth);

		config_setting_lookup_string(prepro_setting, "dark_lib", &dark);
		com.pref.prepro.dark_lib = g_strdup(dark);
		config_setting_lookup_bool(prepro_setting, "use_dark_lib", &com.pref.prepro.use_dark_lib);

		config_setting_lookup_string(prepro_setting, "flat_lib", &flat);
		com.pref.prepro.flat_lib = g_strdup(flat);
		config_setting_lookup_bool(prepro_setting, "use_flat_lib", &com.pref.prepro.use_flat_lib);

		prepro_setting = config_lookup(&config, "prepro-settings.xtrans_af");
		if (prepro_setting != NULL) {
			com.pref.prepro.xtrans_af.x = config_setting_get_int_elem(prepro_setting, 0);
			com.pref.prepro.xtrans_af.y = config_setting_get_int_elem(prepro_setting, 1);
			com.pref.prepro.xtrans_af.w = config_setting_get_int_elem(prepro_setting, 2);
			com.pref.prepro.xtrans_af.h = config_setting_get_int_elem(prepro_setting, 3);
		}
		prepro_setting = config_lookup(&config, "prepro-settings.xtrans_sample");
		if (prepro_setting != NULL) {
			com.pref.prepro.xtrans_sample.x = config_setting_get_int_elem(prepro_setting, 0);
			com.pref.prepro.xtrans_sample.y = config_setting_get_int_elem(prepro_setting, 1);
			com.pref.prepro.xtrans_sample.w = config_setting_get_int_elem(prepro_setting, 2);
			com.pref.prepro.xtrans_sample.h = config_setting_get_int_elem(prepro_setting, 3);
		}
	}

	/* Registration setting */
	config_setting_t *reg_setting = config_lookup(&config, keywords[REG]);
	if (reg_setting) {
		config_setting_lookup_int(reg_setting, "method", &com.pref.gui.reg_settings);
	}

	/* Stacking setting */
	config_setting_t *stack_setting = config_lookup(&config, keywords[STK]);
	if (stack_setting) {
		config_setting_lookup_int(stack_setting, "method", &com.pref.stack.method);
		config_setting_lookup_int(stack_setting, "rejection", &com.pref.stack.rej_method);
		config_setting_lookup_int(stack_setting, "normalisation", &com.pref.stack.normalisation_method);
		config_setting_lookup_float(stack_setting, "sigma_low", &com.pref.stack.sigma_low);
		config_setting_lookup_float(stack_setting, "sigma_high", &com.pref.stack.sigma_high);
		config_setting_lookup_float(stack_setting, "linear_low", &com.pref.stack.linear_low);
		config_setting_lookup_float(stack_setting, "linear_high", &com.pref.stack.linear_high);
		config_setting_lookup_float(stack_setting, "percentile_low", &com.pref.stack.percentile_low);
		config_setting_lookup_float(stack_setting, "percentile_high", &com.pref.stack.percentile_high);


		config_setting_lookup_int(stack_setting, "mem_mode", (int*)&com.pref.mem_mode);
		config_setting_lookup_float(stack_setting, "maxmem", &com.pref.memory_ratio);
		config_setting_lookup_float(stack_setting, "maxmem_gb",	&com.pref.memory_amount);
	}
	if (com.pref.mem_mode < 0 || com.pref.mem_mode > 2)
		com.pref.mem_mode = RATIO;
	if (com.pref.memory_ratio <= 0.05)
		com.pref.memory_ratio = 0.9;

	/* FITS compression setting */
	config_setting_t *comp_setting = config_lookup(&config, keywords[CMP]);
	if (comp_setting) {
		config_setting_lookup_bool(comp_setting, "compress_enabled", &com.pref.comp.fits_enabled);
		config_setting_lookup_int(comp_setting, "compress_method", &com.pref.comp.fits_method);
		config_setting_lookup_float(comp_setting, "compress_quantization", &com.pref.comp.fits_quantization);
		config_setting_lookup_float(comp_setting, "compress_hcompress_scale", &com.pref.comp.fits_hcompress_scale);
	}

	/* Astrometry setting */
	config_setting_t *astrometry_setting = config_lookup(&config, keywords[AST]);
	if (astrometry_setting) {
		config_setting_lookup_bool(astrometry_setting, "messier", &com.pref.gui.catalog[0]);
		config_setting_lookup_bool(astrometry_setting, "ngc", &com.pref.gui.catalog[1]);
		config_setting_lookup_bool(astrometry_setting, "ic", &com.pref.gui.catalog[2]);
		config_setting_lookup_bool(astrometry_setting, "ldn", &com.pref.gui.catalog[3]);
		config_setting_lookup_bool(astrometry_setting, "sh2", &com.pref.gui.catalog[4]);
		config_setting_lookup_bool(astrometry_setting, "stars", &com.pref.gui.catalog[5]);
		config_setting_lookup_bool(astrometry_setting, "user", &com.pref.gui.catalog[6]);
		config_setting_lookup_int(astrometry_setting, "compass_position", &com.pref.gui.position_compass);
		config_setting_lookup_int(astrometry_setting, "wcs_formalism", &com.pref.wcs_formalism);

	} else {
		for (int i = 0; i < 6; i ++) {
			com.pref.gui.catalog[i] = TRUE;
		}
	}

	/* Photometry setting */
	config_setting_t *photometry_setting = config_lookup(&config, keywords[PTM]);
	if (photometry_setting) {
		config_setting_lookup_float(photometry_setting, "gain", &com.pref.phot_set.gain);
		config_setting_lookup_float(photometry_setting, "inner-radius", &com.pref.phot_set.inner);
		config_setting_lookup_float(photometry_setting, "outer-radius", &com.pref.phot_set.outer);
		config_setting_lookup_float(photometry_setting, "aperture-radius", &com.pref.phot_set.aperture);
		config_setting_lookup_bool(photometry_setting, "force-radius", &com.pref.phot_set.force_radius);
		//config_setting_lookup_float(photometry_setting, "minval", &com.pref.phot_set.minval);
		// we don't want to keep this old value ^ in the new version because it's now negative by default
		config_setting_lookup_float(photometry_setting, "maxval", &com.pref.phot_set.maxval);
		if (com.pref.phot_set.inner == 0.0 || com.pref.phot_set.outer == 0.0) {
			initialize_photometric_param();
		}
	}

	/* Misc setting */
	config_setting_t *misc_setting = config_lookup(&config, keywords[MISC]);
	if (misc_setting) {
		int type;
		const char *swap_dir = NULL, *starnet_dir = NULL, *extension = NULL, *lang = NULL, *copyright = NULL;

		config_setting_lookup_int(misc_setting, "pan_position", &com.pref.gui.pan_position);
		config_setting_lookup_int(misc_setting, "hd_bitdepth", &com.pref.hd_bitdepth);

		if (config_setting_lookup_bool(misc_setting, "is_extended", &com.pref.gui.is_extended) == CONFIG_FALSE) {
			com.pref.gui.is_extended = TRUE;
		}
		if (config_setting_lookup_bool(misc_setting, "first_start_1", &com.pref.gui.first_start) == CONFIG_FALSE) {
			com.pref.gui.first_start = TRUE;
		}
		if (config_setting_lookup_bool(misc_setting, "confirm_quit", &com.pref.gui.silent_quit) == CONFIG_FALSE) {
			com.pref.gui.silent_quit = FALSE;
		}
		if (config_setting_lookup_bool(misc_setting, "scripts_warning", &com.pref.gui.warn_script_run) == CONFIG_FALSE) {
			com.pref.gui.warn_script_run = TRUE;
		}
		if (config_setting_lookup_bool(misc_setting, "check_requires", &com.pref.script_check_requires) == CONFIG_FALSE) {
			com.pref.script_check_requires = TRUE;
		}
		if (config_setting_lookup_bool(misc_setting, "show_thumbnails", &com.pref.gui.show_thumbnails) == CONFIG_FALSE) {
			com.pref.gui.show_thumbnails = TRUE;
		}
		if (config_setting_lookup_bool(misc_setting, "remember_winpos", &com.pref.gui.remember_windows) == CONFIG_FALSE) {
			com.pref.gui.remember_windows = TRUE;
		}
#ifdef HAVE_JSON_GLIB
		if (config_setting_lookup_bool(misc_setting, "check_update_at_startup", &com.pref.check_update) == CONFIG_FALSE) {
			com.pref.check_update = TRUE;
		}
#else
		com.pref.check_update = FALSE;
#endif
		if (config_setting_lookup_float(misc_setting, "font_scale", &com.pref.gui.font_scale) == CONFIG_FALSE) {
			com.pref.gui.font_scale = 100.0;
		}
		if (config_setting_lookup_bool(misc_setting, "icon_symbolic", &com.pref.gui.icon_symbolic) == CONFIG_FALSE) {
			com.pref.gui.icon_symbolic = FALSE;
		}
		if (config_setting_lookup_bool(misc_setting, "rgb_aladin", &com.pref.rgb_aladin) == CONFIG_FALSE) {
			com.pref.rgb_aladin = FALSE;
		}
		if (config_setting_lookup_float(misc_setting, "focal", &com.pref.focal) == CONFIG_FALSE) {
			com.pref.focal = 1000.0;
		}
		if (config_setting_lookup_float(misc_setting, "pitch", &com.pref.pitch) == CONFIG_FALSE) {
			com.pref.pitch = 5.0;
		}
		config_setting_lookup_int(misc_setting, "thumbnail_size", &com.pref.gui.thumbnail_size);
		config_setting_lookup_int(misc_setting, "theme", &com.pref.gui.combo_theme);
		config_setting_lookup_string(misc_setting, "lang", &lang);
		com.pref.lang = g_strdup(lang);
		config_setting_lookup_bool(misc_setting, "is_maximized", &com.pref.gui.is_maximized);
		config_setting_lookup_string(misc_setting, "swap_directory", &swap_dir);
		com.pref.swap_dir = g_strdup(swap_dir);
		config_setting_lookup_string(misc_setting, "starnet_directory", &starnet_dir);
		com.pref.starnet_dir = g_strdup(starnet_dir);
		config_setting_lookup_string(misc_setting, "extension", &extension);
		com.pref.ext = g_strdup(extension);
		config_setting_lookup_int(misc_setting, "FITS_type", &type);
		com.pref.force_16bit = (type == 0);
		config_setting_lookup_int(misc_setting, "selection_guides", &com.pref.gui.selection_guides);
		config_setting_lookup_string(misc_setting, "copyright", &copyright);
		com.pref.copyright = g_strdup(copyright);

		misc_setting = config_lookup(&config, "misc-settings.scripts_paths");
		if (misc_setting != NULL) {
			unsigned int count = config_setting_length(misc_setting);
			unsigned int i;
			const char *tmp = NULL;

			for (i = 0; i < count; ++i) {
				tmp = config_setting_get_string_elem(misc_setting, i);
				list = g_slist_append(list, g_strdup(tmp));
			}
		}
		misc_setting = config_lookup(&config, "misc-settings.main_w_pos");
		if (misc_setting != NULL) {
			com.pref.gui.main_w_pos.x = config_setting_get_int_elem(misc_setting, 0);
			com.pref.gui.main_w_pos.y = config_setting_get_int_elem(misc_setting, 1);
			com.pref.gui.main_w_pos.w = config_setting_get_int_elem(misc_setting, 2);
			com.pref.gui.main_w_pos.h = config_setting_get_int_elem(misc_setting, 3);
		}

	}
	com.pref.gui.script_path = list;
	config_destroy(&config);
	return 0;
}

#endif

static int get_key_data(GKeyFile *kf, struct settings_access *desc) {
	GError *error = NULL;
	gboolean boolval;
	int intval;
	double doubleval;
	gchar *strval;
	gsize len;
	gchar **strs;
	switch (desc->type) {
		case STYPE_BOOL:
			boolval = g_key_file_get_boolean(kf, desc->group, desc->key, &error);
			if (error && error->code == G_KEY_FILE_ERROR_INVALID_VALUE) {
				siril_log_message(_("error in config file for %s.%s: %s (value: %s)\n"),
						desc->group, desc->key, error->message,
						g_key_file_get_string(kf, desc->group, desc->key, NULL));
				return 1;
			}
			*((gboolean*)desc->data) = boolval;
			break;
		case STYPE_INT:
			intval = g_key_file_get_integer(kf, desc->group, desc->key, &error);
			if (error && error->code == G_KEY_FILE_ERROR_INVALID_VALUE) {
				siril_log_message(_("error in config file for %s.%s: %s (value: %s)\n"),
						desc->group, desc->key, error->message,
						g_key_file_get_string(kf, desc->group, desc->key, NULL));
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
				siril_log_message(_("error in config file for %s.%s: %s (value: %s)\n"),
						desc->group, desc->key, error->message,
						g_key_file_get_string(kf, desc->group, desc->key, NULL));
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
			if (strval[0] == '\0')
				return 1;
			if (desc->type == STYPE_STRDIR && !g_file_test(strval, G_FILE_TEST_IS_DIR)) {
				siril_log_color_message(_("directory `%s' for config key %s.%s doesn't exist, not using it.\n"),
						"salmon", strval, desc->group, desc->key);
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
				GSList *old_list = *((GSList**)desc->data);
				if (old_list)
					g_slist_free_full(old_list, g_free);
				*((GSList**)desc->data) = list;
				g_free(strs);
			}
			break;
	}
	return 0;
}

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
	}
	siril_debug_print("read %zd keys from key file\n", nb_keys_read);
	return nb_keys_read == 0;
}

int readinitfile(char *path) {
	GKeyFile *kf = g_key_file_new();
	gchar *fname = get_locale_filename(path);
	GError *error = NULL;
	if (!g_key_file_load_from_file(kf, fname, G_KEY_FILE_NONE, &error)) {
		if (error != NULL) {
			siril_log_color_message(_("Settings could not be loaded from %s: %s\n"), "red", fname, error->message);
			g_clear_error(&error);
		}
		g_free(fname);
		return 1;
	}
	g_free(fname);
#ifndef HAVE_JSON_GLIB
	com.pref.check_update = FALSE;
#endif
	return read_keyfile(kf);
}

/**
 * Public functions
 */


int checkinitfile() {
	/* com.initfile will contain the path passed with -i if any, NULL else */
	if (com.initfile) {
		siril_log_message(_("Reading configuration file %s\n"), com.initfile);
		return readinitfile(com.initfile);
	}

	int retval = 0;

	/* set com.initfile to default location */
	gchar *pathname = g_build_filename(siril_get_config_dir(), PACKAGE, NULL);
	gchar *config_file = g_build_filename(pathname, GLIB_CONFIG_FILE, NULL);
	if (!g_file_test(config_file, G_FILE_TEST_EXISTS)) {
#ifdef HAVE_LIBCONFIG
		/* try the old config file */
		gchar *libconfig_file = g_build_filename(pathname, LIBCONFIG_FILE, NULL);
		if (g_file_test(libconfig_file, G_FILE_TEST_EXISTS) &&
				!readinitfile_libconfig(libconfig_file)) {
			siril_log_color_message(_("With this new version, the format of the configuration file changed. Settings were imported from the old format. This is a one-time message.\n"), "salmon");
			com.initfile = config_file;
			if (writeinitfile()) {
				siril_log_color_message(_("But the new file failed to be saved to %s, changes in settings will not be saved\n"), "red", config_file);
				com.initfile = NULL;
			}
			else siril_log_message(_("New settings saved in %s\n"), com.initfile);
		} else {
#endif
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
#ifdef HAVE_LIBCONFIG
		}
		g_free(libconfig_file);
#endif
	}
	else {
		com.initfile = config_file;
		retval = readinitfile(com.initfile);
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
		return 1;
	}
	g_key_file_free(kf);
	return 0;
}
