/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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
#else
#define SIRIL_UNSTABLE 1
#endif
#include "core/settings.h"
#include "core/siril.h"
#include "core/siril_log.h"

/* the settings as initialized in static.
 * the dynamic fields are set in initialize_default_settings() */
preferences pref_init = {
	.wd = NULL,
	.ext = NULL,
	.force_16bit = FALSE,
	.mem_mode = RATIO,
	.memory_ratio = 0.9,
	.memory_amount = 10,
	.script_check_requires = TRUE,
	.pipe_check_requires = FALSE,
	.check_update = !SIRIL_UNSTABLE,
	.lang = 0,
	.swap_dir = NULL,
	.focal = 1000,
	.pitch = 5,
	/* catalogs for astrometry */
	.catalog[0] = TRUE,
	.catalog[1] = TRUE,
	.catalog[2] = TRUE,
	.catalog[3] = TRUE,
	.catalog[4] = TRUE,
	.catalog[5] = TRUE,
	.catalog[6] = TRUE,
	.wcs_formalism = WCS_FORMALISM_1,
	.rgb_aladin = FALSE,
	.copyright = NULL,
	.starfinder_conf = { // starfinder_conf
		.radius = 10,
		.adjust = TRUE,
		.sigma = 1.0,
		.roundness = 0.5,
		.focal_length = 0.,
		.pixel_size_x = 0.,
	},
	.prepro = {
		.cfa = FALSE,
		.equalize_cfa = TRUE,
		.fix_xtrans = FALSE,
		{ // xtrans_af
				.x = 0,
				.y = 0,
				.w = 0,
				.h = 0
		},
		{ //xtrans_sample
				.x = 0,
				.y = 0,
				.w = 0,
				.h = 0
		},
		.bias_lib = NULL,
		.use_bias_lib = FALSE,
		.bias_synth = NULL,
		.use_bias_synth = FALSE,
		.dark_lib = NULL,
		.use_dark_lib = FALSE,
		.flat_lib = NULL,
		.use_flat_lib = FALSE,
	},
	.gui = {
		.first_start = TRUE,
		.confirm_quit = TRUE,
		.remember_windows = TRUE,
		.main_w_pos = {
				.x = 0,
				.y = 0,
				.w = 0,
				.h = 0
		},
		.pan_position = -1,
		.is_extended = TRUE,
		.is_maximized = FALSE,
		.combo_theme = 0,
		.font_scale = 100,
		.icon_symbolic = FALSE,
		.script_path = NULL,
		.warn_script_run = TRUE,
		.show_thumbnails = TRUE,
		.thumbnail_size = 256,
		.position_compass = 1,
		.selection_guides = 0,
		.reg_settings = 0
	},
	.debayer = {
		.open_debayer = FALSE,
		.use_bayer_header = TRUE,
		.bayer_pattern = BAYER_FILTER_RGGB,
		.bayer_inter = BAYER_RCD,
		.top_down = TRUE,
		.xbayeroff = 0,
		.ybayeroff = 0,
	},
	.phot_set = {
		.gain = 2.3,
		.inner = 20.0,
		.outer = 30.0,
		.aperture = 10.0,
		.force_radius = FALSE,
		.minval = 0,
		.maxval = 60000,
	},
	.stack = {
		.method = 0,
		.normalisation_method = ADDITIVE_SCALING,
		.rej_method = WINSORIZED,
		.sigma_low = 3.0, .sigma_high = 3.0,
		.linear_low = 5.0, .linear_high = 5.0,
		.percentile_low = 3.0, .percentile_high = 3.0,
	},
	.comp = {
		.fits_enabled = FALSE,
		.fits_method = 0,
		.fits_quantization = 16.0,
		.fits_hcompress_scale = 4.0,
	}
};

static void initialize_settings_to_default() {
	memcpy(&com.pref, &pref_init, sizeof(preferences));
}

void free_preferences(preferences *pref) {
	g_free(pref->ext);
	pref->ext = NULL;
	g_free(pref->swap_dir);
	pref->swap_dir = NULL;
	g_free(pref->copyright);
	pref->copyright = NULL;
	g_free(pref->lang);
	pref->lang = NULL;
	g_slist_free_full(pref->gui.script_path, g_free);
	pref->gui.script_path = NULL;
}

/* static + dynamic settings initialization */
void initialize_default_settings() {
	initialize_settings_to_default();
	com.pref.ext = g_strdup(".fit");
	com.pref.swap_dir = g_strdup(g_get_tmp_dir());
	//com.pref.wd = g_strdup(siril_get_startup_dir());
	// TODO check what's needed here
}

struct settings_access all_settings[] = {
	{ "core", "wd", STYPE_STRDIR, "current working directory", &com.pref.wd },
	{ "core", "extension", STYPE_STR, "FITS file extension", &com.pref.ext },
	{ "core", "force_16bit", STYPE_BOOL, "don't use 32 bits for pixel depth", &com.pref.force_16bit },
	{ "core", "mem_mode", STYPE_INT, "memory mode (0 ratio, 1 amount)", &com.pref.mem_mode, { .range_int = { 0, 1 } } },
	{ "core", "mem_ratio", STYPE_DOUBLE, "memory ratio of available", &com.pref.memory_ratio, { .range_double = { 0.05, 4.0 } }  },
	{ "core", "mem_amount", STYPE_DOUBLE, "amount of memory in GB", &com.pref.memory_amount },
	{ "core", "script_check_requires", STYPE_BOOL, "need requires cmd in pupe", &com.pref.script_check_requires },
	{ "core", "pipe_check_requires", STYPE_BOOL, "need requires cmd in pupe", &com.pref.pipe_check_requires },
	{ "core", "check_updates", STYPE_BOOL, "check update at start-up", &com.pref.check_update },
	{ "core", "lang", STYPE_STR, "active siril language", &com.pref.lang },
	{ "core", "swap_dir", STYPE_STRDIR, "swap directory", &com.pref.swap_dir },
	{ "core", "wcs_formalism", STYPE_INT, "WCS formalism used in FITS header", &com.pref.wcs_formalism },
	{ "core", "rgb_aladin", STYPE_BOOL, "add CTYPE3='RGB' in the FITS header", &com.pref.rgb_aladin },
	{ "core", "copyright", STYPE_STR, "user copyright to put in file header", &com.pref.copyright },

	// TODO: complete the list. If a field is not listed, it won't be read from or saved to ini file

	{ "gui", "first_start", STYPE_BOOL, "first start of siril", &com.pref.gui.first_start },
	{ "gui", "remember_windows", STYPE_BOOL, "remember window position", &com.pref.gui.remember_windows },
	{ "gui", "script_path", STYPE_STRLIST, "list of script directories", &com.pref.gui.script_path },

	{ NULL, NULL, STYPE_BOOL, NULL, NULL }
};

struct settings_access *get_key_settings(const char *group, const char *key) {
	int nb_settings = sizeof(all_settings) / sizeof(struct settings_access) - 1;
	for (int i = 0; i < nb_settings; i++) {
		if (!strcmp(all_settings[i].group, group) && !strcmp(all_settings[i].key, key))
			return all_settings + i;
	}
	return NULL;
}

struct settings_access *get_all_settings() {
	return all_settings;
}

static const char *settings_type_to_string(enum settings_type type) {
	switch (type) {
		case STYPE_BOOL:
			return "boolean";
		case STYPE_INT:
			return "integer";
		case STYPE_DOUBLE:
			return "double";
		case STYPE_STR:
			return "string";
		case STYPE_STRDIR:
			return "directory";
		case STYPE_STRLIST:
			return "list of strings";
		default:
			return "unknown";
	}
}

int print_settings_key(const char *group, const char *key, gboolean with_details) {
	struct settings_access *desc = get_key_settings(group, key);
	if (!desc) {
		siril_log_message(_("unknown settings variable %s.%s\n"), group, key);
		return 1;
	}
	GString *str = g_string_sized_new(120);
	g_string_printf(str, "%s.%s = ", desc->group, desc->key);
	GSList *list;
	switch (desc->type) {
		case STYPE_BOOL:
			g_string_append(str, *((gboolean*)desc->data) ? "true" : "false");
			break;
		case STYPE_INT:
			g_string_append_printf(str, "%d", *((int*)desc->data));
			break;
		case STYPE_DOUBLE:
			g_string_append_printf(str, "%f", *((double*)desc->data));
			break;
		case STYPE_STR:
		case STYPE_STRDIR:
			g_string_append(str, *((gchar**)desc->data));
			break;
		case STYPE_STRLIST:
			list = *((GSList**)desc->data);
			while (list) {
				if (list != *((GSList **)desc->data))
					g_string_append(str, " ; ");
				g_string_append(str, list->data);
				list = list->next;
			}
			break;
	}
	if (with_details) {
		g_string_append_printf(str, " (%s)", settings_type_to_string(desc->type));
		g_string_append_printf(str, ", %s", desc->desc);
		// TODO add the non working ranges
	}
	gchar *s = g_string_free(str, FALSE);
	siril_log_message("%s\n", s);
	return 0;
}

int print_all_settings(gboolean with_details) {
	int nb_settings = sizeof(all_settings) / sizeof(struct settings_access) - 1;
	for (int i = 0; i < nb_settings; i++) {
		print_settings_key(all_settings[i].group, all_settings[i].key, with_details);
	}
	return 0;
}
