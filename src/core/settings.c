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
#include "io/catalogues.h"

/* the settings as initialized in static.
 * the dynamic fields are set in initialize_default_settings() */
preferences pref_init = {
	.wd = NULL,
	.ext = NULL,
	.add_fz = FALSE,
	.force_16bit = FALSE,
	.allow_heterogeneous_fitseq = FALSE,
	.mem_mode = RATIO,
	.memory_ratio = 0.9,
	.memory_amount = 10,
	.hd_bitdepth = 20,
	.script_check_requires = TRUE,
	.pipe_check_requires = FALSE,
#ifdef HAVE_JSON_GLIB
	.check_update = !SIRIL_UNSTABLE,
#else
	.check_update = FALSE,
#endif
	.lang = 0,
	.swap_dir = NULL,
	.focal = 1000,
	.pitch = 5,
	.wcs_formalism = WCS_FORMALISM_1,
	.catalogue_paths[0] = NULL,
	.catalogue_paths[1] = NULL,
	.catalogue_paths[2] = NULL,
	.catalogue_paths[3] = NULL,
	.rgb_aladin = FALSE,
	.copyright = NULL,
	.starnet_dir = NULL,
	.starfinder_conf = { // starfinder_conf
		.radius = DEF_BOX_RADIUS,
		.sigma = 1.0,
		.roundness = 0.5,
		.focal_length = 0.,
		.pixel_size_x = 0.,
		.convergence = 1,
		.relax_checks = FALSE,
		.profile = PSF_GAUSSIAN,
		.min_beta = 1.5
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
		.dark_lib = NULL,
		.use_dark_lib = FALSE,
		.flat_lib = NULL,
		.use_flat_lib = FALSE,
		.stack_default = NULL,
		.use_stack_default = TRUE,
	},
	.gui = {
		.first_start = TRUE,
		.silent_quit = FALSE,
		.silent_linear = FALSE,
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
		.default_rendering_mode = LINEAR_DISPLAY,
		.display_histogram_mode = LINEAR_DISPLAY,
		.catalog[0] = TRUE,
		.catalog[1] = TRUE,
		.catalog[2] = TRUE,
		.catalog[3] = TRUE,
		.catalog[4] = TRUE,
		.catalog[5] = TRUE,
		.catalog[6] = TRUE,
		.catalog[7] = TRUE,
		.catalog[8] = TRUE,
		.position_compass = 1,
		.selection_guides = 0,
		.reg_settings = 0,
		.reg_interpolation = OPENCV_LANCZOS4,
		.reg_clamping = TRUE,
		.pm_presets = NULL
	},
	.debayer = {
		.open_debayer = FALSE,
		.use_bayer_header = TRUE,
		.bayer_pattern = BAYER_FILTER_RGGB,
		.bayer_inter = BAYER_RCD,
		.top_down = TRUE,
		.xbayeroff = 0,
		.ybayeroff = 0,
		.xtrans_passes = 1
	},
	.phot_set = {
		.gain = 2.3,
		.inner = 20.0,
		.outer = 30.0,
		.aperture = 10.0,
		.force_radius = FALSE,
		.minval = -500.0,
		.maxval = 60000.0,
	},
	.analysis = {
		.mosaic_panel = 256,
		.mosaic_window = 381,
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

void free_preferences(preferences *pref) {
	g_free(pref->ext);
	pref->ext = NULL;
	g_free(pref->swap_dir);
	pref->swap_dir = NULL;
	g_free(pref->copyright);
	pref->copyright = NULL;
	g_free(pref->starnet_dir);
	pref->starnet_dir = NULL;
	g_free(pref->lang);
	pref->lang = NULL;
	g_slist_free_full(pref->gui.script_path, g_free);
	pref->gui.script_path = NULL;
}

/* static + dynamic settings initialization */
void initialize_default_settings() {
	com.pref = pref_init;
	com.pref.ext = g_strdup(".fit");
	com.pref.prepro.stack_default = g_strdup("$seqname$stacked");
	com.pref.swap_dir = g_strdup(g_get_tmp_dir());
	initialize_local_catalogues_paths();
}

void update_gain_from_gfit() {
	if (gfit.cvf > 0.0)
		com.pref.phot_set.gain = gfit.cvf;
}

struct settings_access all_settings[] = {
	{ "core", "wd", STYPE_STRDIR, N_("current working directory"), &com.pref.wd },
	{ "core", "extension", STYPE_STR, N_("FITS file extension"), &com.pref.ext },
	{ "core", "fz_extension", STYPE_BOOL, N_("add fz to the FITS extension"), &com.pref.add_fz },
	{ "core", "force_16bit", STYPE_BOOL, N_("don't use 32 bits for pixel depth"), &com.pref.force_16bit },
	{ "core", "allow_heterogeneous_fitseq", STYPE_BOOL, N_("allow FITS cubes to have different sizes"), &com.pref.allow_heterogeneous_fitseq },
	{ "core", "mem_mode", STYPE_INT, N_("memory mode (0 ratio, 1 amount)"), &com.pref.mem_mode, { .range_int = { 0, 1 } } },
	{ "core", "mem_ratio", STYPE_DOUBLE, N_("memory ratio of available"), &com.pref.memory_ratio, { .range_double = { 0.05, 4.0 } } },
	{ "core", "mem_amount", STYPE_DOUBLE, N_("amount of memory in GB"), &com.pref.memory_amount, { .range_double = { 0.1, 1000000. } } },
	{ "core", "hd_bitdepth", STYPE_INT, N_("HD AutoStretch bit depth"), &com.pref.hd_bitdepth, { .range_int = { 17, 24 } } },
	{ "core", "script_check_requires", STYPE_BOOL, N_("need requires cmd in script"), &com.pref.script_check_requires },
	{ "core", "pipe_check_requires", STYPE_BOOL, N_("need requires cmd in pipe"), &com.pref.pipe_check_requires },
	{ "core", "check_updates", STYPE_BOOL, N_("check update at start-up"), &com.pref.check_update },
	{ "core", "lang", STYPE_STR, N_("active siril language"), &com.pref.lang },
	{ "core", "swap_dir", STYPE_STRDIR, N_("swap directory"), &com.pref.swap_dir },
	{ "core", "wcs_formalism", STYPE_INT, N_("WCS formalism used in FITS header"), &com.pref.wcs_formalism, { .range_int = { 0, 1 } } },
	{ "core", "catalogue_namedstars", STYPE_STR, N_("Path of the namedstars.dat catalogue"), &com.pref.catalogue_paths[0] },
	{ "core", "catalogue_unnamedstars", STYPE_STR, N_("Path of the unnamedstars.dat catalogue"), &com.pref.catalogue_paths[1] },
	{ "core", "catalogue_tycho2", STYPE_STR, N_("Path of the deepstars.dat catalogue"), &com.pref.catalogue_paths[2] },
	{ "core", "catalogue_nomad", STYPE_STR, N_("Path of the USNO-NOMAD-1e8.dat catalogue"), &com.pref.catalogue_paths[3] },
	{ "core", "rgb_aladin", STYPE_BOOL, N_("add CTYPE3='RGB' in the FITS header"), &com.pref.rgb_aladin },
	{ "core", "copyright", STYPE_STR, N_("user copyright to put in file header"), &com.pref.copyright },
	{ "core", "starnet_dir", STYPE_STR, N_("directory of the starnet++ installation"), &com.pref.starnet_dir },

	{ "starfinder", "focal_length", STYPE_DOUBLE, N_("focal length in mm for radius adjustment"), &com.pref.starfinder_conf.focal_length, { .range_double = { 0., 999999. } } },
	{ "starfinder", "pixel_size", STYPE_DOUBLE, N_("pixel size in µm for radius adjustment"), &com.pref.starfinder_conf.pixel_size_x, { .range_double = { 0., 99. } } },

	{ "debayer", "use_debayer_header", STYPE_BOOL, N_("use pattern from the file header"), &com.pref.debayer.use_bayer_header },
	{ "debayer", "pattern", STYPE_INT, N_("index of the Bayer pattern"), &com.pref.debayer.bayer_pattern, { .range_int = { 0, XTRANS_FILTER_4 } } },
	{ "debayer", "interpolation", STYPE_INT, N_("type of interpolation"), &com.pref.debayer.bayer_inter, { .range_int = { 0, XTRANS } } },
	{ "debayer", "top_down", STYPE_BOOL, N_("force debayer top-down"), &com.pref.debayer.top_down },
	{ "debayer", "offset_x", STYPE_INT, N_("Bayer matrix offset X"), &com.pref.debayer.xbayeroff, { .range_int = { 0, 1 } } },
	{ "debayer", "offset_y", STYPE_INT, N_("Bayer matrix offset Y"), &com.pref.debayer.ybayeroff, { .range_int = { 0, 1 } } },
	{ "debayer", "xtrans_passes", STYPE_INT, N_("Number of passes for the X-Trans Markesteijn algorithm"), &com.pref.debayer.xtrans_passes, { .range_int = { 1, 4 } } },

	{ "photometry", "gain", STYPE_DOUBLE, N_("electrons per ADU for noise estimation"), &com.pref.phot_set.gain, { .range_double = { 0., 10. } } },
	{ "photometry", "inner", STYPE_DOUBLE, N_("inner radius for background annulus"), &com.pref.phot_set.inner, { .range_double = { 2., 100. } } },
	{ "photometry", "outer", STYPE_DOUBLE, N_("outer radius for background annulus"), &com.pref.phot_set.outer, { .range_double = { 3., 200. } } },
	{ "photometry", "force_radius", STYPE_BOOL, N_("force flux aperture value"), &com.pref.phot_set.force_radius },
	{ "photometry", "aperture", STYPE_DOUBLE, N_("forced aperture for flux computation"), &com.pref.phot_set.aperture, { .range_double = { 1., 100. } } },
	{ "photometry", "minval", STYPE_DOUBLE, N_("minimum valid pixel value for photometry"), &com.pref.phot_set.minval, { .range_double = { -65536.0, 65534.0 } } },
	{ "photometry", "maxval", STYPE_DOUBLE, N_("maximum valid pixel value for photometry"), &com.pref.phot_set.maxval, { .range_double = { 1.0, 65535.0 } } },

	{ "analysis", "panel", STYPE_INT, N_("panel size of aberration inspector"), &com.pref.analysis.mosaic_panel, { .range_int = { 127, 1024 } } },
	{ "analysis", "window", STYPE_INT, N_("window size of aberration inspector"), &com.pref.analysis.mosaic_window, { .range_int = { 300, 1600 } } },

	{ "compression", "enabled", STYPE_BOOL, N_("FITS compression enabled"), &com.pref.comp.fits_enabled },
	{ "compression", "method", STYPE_INT, N_("FITS compression method"), &com.pref.comp.fits_method, { .range_int = { 0, 3 } } },
	{ "compression", "quantization", STYPE_DOUBLE, N_("quantization factor for 32-bit float"), &com.pref.comp.fits_quantization, { .range_double = { 8., 256. } }  },
	{ "compression", "hcompress_scale", STYPE_DOUBLE, N_("Hcompress scale factor"), &com.pref.comp.fits_hcompress_scale, { .range_double = { 0., 256. } }  },

	/* the GUI part, not as useful but still required to be listed in order to be saved in the ini file */
	{ "gui_prepro", "cfa", STYPE_BOOL, N_("type of sensor for cosmetic correction"), &com.pref.prepro.cfa },
	{ "gui_prepro", "equalize_cfa", STYPE_BOOL, N_("equalize flat channels"), &com.pref.prepro.equalize_cfa },
	{ "gui_prepro", "fix_xtrans", STYPE_BOOL, N_("enable correction for X-Trans sensor"), &com.pref.prepro.fix_xtrans },
	{ "gui_prepro", "xtrans_af_x", STYPE_INT, N_("if no X-Trans model found, use this"), &com.pref.prepro.xtrans_af.x },
	{ "gui_prepro", "xtrans_af_y", STYPE_INT, N_("if no X-Trans model found, use this"), &com.pref.prepro.xtrans_af.y },
	{ "gui_prepro", "xtrans_af_w", STYPE_INT, N_("if no X-Trans model found, use this"), &com.pref.prepro.xtrans_af.w },
	{ "gui_prepro", "xtrans_af_h", STYPE_INT, N_("if no X-Trans model found, use this"), &com.pref.prepro.xtrans_af.h },
	{ "gui_prepro", "xtrans_sample_x", STYPE_INT, N_("if no X-Trans model found, use this"), &com.pref.prepro.xtrans_sample.x },
	{ "gui_prepro", "xtrans_sample_y", STYPE_INT, N_("if no X-Trans model found, use this"), &com.pref.prepro.xtrans_sample.y },
	{ "gui_prepro", "xtrans_sample_w", STYPE_INT, N_("if no X-Trans model found, use this"), &com.pref.prepro.xtrans_sample.w },
	{ "gui_prepro", "xtrans_sample_h", STYPE_INT, N_("if no X-Trans model found, use this"), &com.pref.prepro.xtrans_sample.h },
	{ "gui_prepro", "bias_lib", STYPE_STR, N_("default master bias"), &com.pref.prepro.bias_lib },
	{ "gui_prepro", "use_bias_lib", STYPE_BOOL, N_("use default master bias"), &com.pref.prepro.use_bias_lib },
	{ "gui_prepro", "dark_lib", STYPE_STR, N_("default master dark"), &com.pref.prepro.dark_lib },
	{ "gui_prepro", "use_dark_lib", STYPE_BOOL, N_("use default master dark"), &com.pref.prepro.use_dark_lib },
	{ "gui_prepro", "flat_lib", STYPE_STR, N_("default master flat"), &com.pref.prepro.flat_lib },
	{ "gui_prepro", "use_flat_lib", STYPE_BOOL, N_("use default master flat"), &com.pref.prepro.use_flat_lib },
	{ "gui_prepro", "stack_default", STYPE_STR, N_("default stack name"), &com.pref.prepro.stack_default },
	{ "gui_prepro", "use_stack_default", STYPE_BOOL, N_("use preferred stack name"), &com.pref.prepro.use_stack_default },

	{ "gui_registration", "method", STYPE_INT, N_("index of the selected registration method"), &com.pref.gui.reg_settings, { .range_int = { 0, 7 } } },
	{ "gui_registration", "interpolation", STYPE_INT, N_("index of the selected interpolation method"), &com.pref.gui.reg_interpolation, { .range_int = { 0, 5 } } },
	{ "gui_registration", "clamping", STYPE_BOOL, N_("use clamping method with Lanczos and Cubic interpolation"), &com.pref.gui.reg_clamping },

	{ "gui_stack", "method", STYPE_INT, N_("index of the selected method"), &com.pref.stack.method, { .range_int = { 0, 4 } } },
	{ "gui_stack", "normalization", STYPE_INT, N_("index of the normalization method"), &com.pref.stack.normalisation_method, { .range_int = { 0, MULTIPLICATIVE_SCALING } } },
	{ "gui_stack", "rejection", STYPE_INT, N_("index of the rejection method"), &com.pref.stack.rej_method, { .range_int = { 0, GESDT } } },
	{ "gui_stack", "sigma_low", STYPE_DOUBLE, N_("sigma low value for rejection"), &com.pref.stack.sigma_low, { .range_double = { 0., 20. } } },
	{ "gui_stack", "sigma_high", STYPE_DOUBLE, N_("sigma high value for rejection"), &com.pref.stack.sigma_high, { .range_double = { 0., 20. } } },
	{ "gui_stack", "linear_low", STYPE_DOUBLE, N_("linear low value for rejection"), &com.pref.stack.linear_low, { .range_double = { 0., 20. } } },
	{ "gui_stack", "linear_high", STYPE_DOUBLE, N_("linear high value for rejection"), &com.pref.stack.linear_high, { .range_double = { 0., 20. } } },
	{ "gui_stack", "percentile_low", STYPE_DOUBLE, N_("percentile low value for rejection"), &com.pref.stack.percentile_low, { .range_double = { 0., 100. } } },
	{ "gui_stack", "percentile_high", STYPE_DOUBLE, N_("percentile high value for rejection"), &com.pref.stack.percentile_high, { .range_double = { 0., 100. } } },

	{ "gui", "first_start", STYPE_BOOL, N_("first start of siril"), &com.pref.gui.first_start },
	{ "gui", "silent_quit", STYPE_BOOL, N_("don't confirm quit when exiting"), &com.pref.gui.silent_quit },
	{ "gui", "silent_linear", STYPE_BOOL, N_("don't confirm save when non linear mode"), &com.pref.gui.silent_linear },
	{ "gui", "remember_windows", STYPE_BOOL, N_("remember window position"), &com.pref.gui.remember_windows },
	{ "gui", "main_win_pos_x", STYPE_INT, N_("main window position"), &com.pref.gui.main_w_pos.x },
	{ "gui", "main_win_pos_y", STYPE_INT, N_("main window position"), &com.pref.gui.main_w_pos.y },
	{ "gui", "main_win_pos_w", STYPE_INT, N_("main window position"), &com.pref.gui.main_w_pos.w },
	{ "gui", "main_win_pos_h", STYPE_INT, N_("main window position"), &com.pref.gui.main_w_pos.h },
	{ "gui", "pan_position", STYPE_INT, N_("position of the two sides separator"), &com.pref.gui.pan_position },
	{ "gui", "extended", STYPE_BOOL, N_("main window is extended"), &com.pref.gui.is_extended },
	{ "gui", "maximized", STYPE_BOOL, N_("main window is maximized"), &com.pref.gui.is_maximized },
	{ "gui", "theme", STYPE_INT, N_("index of the selected theme"), &com.pref.gui.combo_theme, { .range_int = { 0, 1 } } },
	{ "gui", "font_scale", STYPE_DOUBLE, N_("font scale in percent"), &com.pref.gui.font_scale },
	{ "gui", "icon_symbolic", STYPE_BOOL, N_("icon style"), &com.pref.gui.icon_symbolic },
	{ "gui", "script_path", STYPE_STRLIST, N_("list of script directories"), &com.pref.gui.script_path },
	{ "gui", "warn_script_run", STYPE_BOOL, N_("warn when launching a script"), &com.pref.gui.warn_script_run },
	{ "gui", "show_thumbnails", STYPE_BOOL, N_("show thumbnails in open dialog"), &com.pref.gui.show_thumbnails },
	{ "gui", "thumbnail_size", STYPE_INT, N_("size of the thumbnails"), &com.pref.gui.thumbnail_size },
	{ "gui", "selection_guides", STYPE_INT, N_("number of elements of the grid guides"), &com.pref.gui.selection_guides },
	{ "gui", "default_rendering_mode", STYPE_INT, N_("default display mode"), &com.pref.gui.default_rendering_mode, { .range_int = { 0, 6 } } },
	{ "gui", "display_histogram_mode", STYPE_INT, N_("default histogram display mode"), &com.pref.gui.display_histogram_mode, { .range_int = { 0, 1 } } },

	{ "gui_astrometry", "compass_position", STYPE_INT, N_("index of the compass position over grid"), &com.pref.gui.position_compass, { .range_int = { 0, 5 } } },
	{ "gui_astrometry", "cat_messier", STYPE_BOOL, N_("show Messier objects in annotations"), &com.pref.gui.catalog[0] },
	{ "gui_astrometry", "cat_ngc", STYPE_BOOL, N_("show NGC objects in annotations"), &com.pref.gui.catalog[1] },
	{ "gui_astrometry", "cat_ic", STYPE_BOOL, N_("show IC objects in annotations"), &com.pref.gui.catalog[2] },
	{ "gui_astrometry", "cat_ldn", STYPE_BOOL, N_("show LDN objects in annotations"), &com.pref.gui.catalog[3] },
	{ "gui_astrometry", "cat_sh2", STYPE_BOOL, N_("show SH2 objects in annotations"), &com.pref.gui.catalog[4] },
	{ "gui_astrometry", "cat_stars", STYPE_BOOL, N_("show stars in annotations"), &com.pref.gui.catalog[5] },
	{ "gui_astrometry", "cat_user", STYPE_BOOL, N_("show user objects in annotations"), &com.pref.gui.catalog[6] },

	{ "gui_pixelmath", "pm_presets", STYPE_STRLIST, N_("list of pixel math presets"), &com.pref.gui.pm_presets },

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
		siril_log_message(_("Unknown settings variable %s.%s\n"), group, key);
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
			g_string_append_printf(str, "%g", *((double*)desc->data));
			break;
		case STYPE_STR:
		case STYPE_STRDIR:
			g_string_append(str, *((gchar**)desc->data) ? *((gchar**)desc->data) : "(not set)");
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
		if (desc->type == STYPE_INT && (desc->range_int.min != 0 || desc->range_int.max != 0))
			g_string_append_printf(str, " [%d, %d]",
					desc->range_int.min, desc->range_int.max);
		else if (desc->type == STYPE_DOUBLE && (desc->range_double.min != 0.0 || desc->range_double.max != 0.0))
			g_string_append_printf(str, " [%g, %g]",
					desc->range_double.min, desc->range_double.max);
		g_string_append_printf(str, " (%s)", settings_type_to_string(desc->type));
		g_string_append_printf(str, ", %s", desc->desc);
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
