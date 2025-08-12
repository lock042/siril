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
#else
#define SIRIL_UNSTABLE 1
#endif
#include "core/settings.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "io/local_catalogues.h"
#include "stacking/stacking.h"

/* the settings as initialized in static.
 * the dynamic fields are set in initialize_default_settings() */
preferences pref_init = {
	.wd = NULL,
	.ext = NULL,
	.force_16bit = FALSE,
	.allow_heterogeneous_fitseq = FALSE,
	.mem_mode = RATIO,
	.memory_ratio = 0.9,
	.memory_amount = 10,
	.hd_bitdepth = 20,
	.script_check_requires = TRUE,
	.pipe_check_requires = FALSE,
 #ifdef SIRIL_UNSTABLE
	.check_update = FALSE,
 #else
	.check_update = TRUE,
 #endif
	.lang = 0,
	.swap_dir = NULL,
	.binning_update = TRUE,
	.wcs_formalism = WCS_FORMALISM_1,
	.catalogue_paths[0] = NULL,
	.catalogue_paths[1] = NULL,
	.catalogue_paths[2] = NULL,
	.catalogue_paths[3] = NULL,
	.catalogue_paths[4] = NULL,
	.catalogue_paths[5] = NULL,
	.rgb_aladin = FALSE,
	.use_checksum = FALSE,
	.copyright = NULL,
	.starnet_exe = NULL,
	.starnet_weights = NULL,
	.graxpert_path = NULL,
	.asnet_dir = NULL,
	.selected_scripts = NULL,
	.use_scripts_repository = TRUE,
	.auto_script_update = TRUE,
	.starfinder_conf = { // starfinder_conf
		.radius = DEF_BOX_RADIUS,
		.sigma = 1.0,
		.roundness = 0.5,
		.focal_length = 0.,
		.pixel_size_x = 0.,
		.convergence = 1,
		.relax_checks = FALSE,
		.profile = PSF_GAUSSIAN,
		.min_beta = 1.5,
		.min_A = 0.0,
		.max_A = 0.0,
		.max_r = 1.0	// min_r is .roundness
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
		.disto_lib = NULL,
		.use_disto_lib = FALSE,
		.stack_default = NULL,
		.use_stack_default = TRUE,
	},
	.gui = {
		.first_start = NULL,
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
		.warn_scripts_run = TRUE,
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
		.catalog[9] = TRUE,
		.catalog[10] = TRUE,
		.position_compass = 1,
		.selection_guides = 0,
		.show_deciasec = FALSE,
		.reg_settings = 0,
		.reg_interpolation = OPENCV_LANCZOS4,
		.reg_clamping = TRUE,
		.pm_presets = NULL,
		.roi_mode = ROI_MANUAL,
		.enable_roi_warning = TRUE,
		.config_colors.color_bkg_samples = NULL,
		.config_colors.color_std_annotations = NULL,
		.config_colors.color_dso_annotations = NULL,
		.config_colors.color_sso_annotations = NULL,
		.config_colors.color_tmp_annotations = NULL,
		.mmb_action = MMB_ZOOM_FIT,
		.mouse_speed_limit = 0.0,
		.mouse_cfg = {
			.mouse_actions_array = NULL,
			.scroll_actions_array = NULL
		},
		.editor_cfg = {
			.highlight_syntax = TRUE,
			.highlight_bracketmatch = TRUE,
			.rmargin = TRUE,
			.rmargin_pos = 80,
			.show_linenums = TRUE,
			.show_linemarks = FALSE,
			.highlight_currentline = TRUE,
			.autoindent = TRUE,
			.indentontab = TRUE,
			.smartbs = TRUE,
			.smarthomeend = TRUE,
			.showspaces = FALSE,
			.shownewlines = FALSE,
			.minimap = FALSE
		}
	},
	.debayer = {
		.open_debayer = FALSE,
		.use_bayer_header = TRUE,
		.bayer_pattern = BAYER_FILTER_RGGB,
		.bayer_inter = BAYER_RCD,
		.orientation = ROW_ORDER_HEADER_TOPDOWN,
		.xbayeroff = 0,
		.ybayeroff = 0,
		.xtrans_passes = 1
	},
	.phot_set = {
		.gain = 2.3,
		.inner = 20.0,
		.outer = 30.0,
		.auto_inner_factor = 4.2,
		.auto_outer_factor = 6.3,
		.aperture = 10.0,
		.auto_aperture_factor = 4.0,
		.force_radius = TRUE,
		.minval = -1500.0,
		.maxval = 60000.0,
		.discard_var_catalogues = 4,
	},
	.astrometry = {
		.update_default_scale = TRUE,
		.percent_scale_range = 20,
		.sip_correction_order = 3,
		.radius_degrees = 10.0,
		.keep_xyls_files = FALSE,
		.keep_wcs_files = FALSE,
		.max_seconds_run = 30,
		.show_asnet_output = FALSE,
		.default_obscode = NULL,
	},
	.analysis = {
		.mosaic_panel = 256,
		.mosaic_window = 381,
	},
	.stack = {
		.method = STACK_SUM,
		.normalisation_method = ADDITIVE_SCALING,
		.weighting_method = NO_WEIGHT,
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
	},
	.fftw_conf = {
		.timelimit = 60,
		.strategy = 0,
		.multithreaded = TRUE,
		.wisdom_file = NULL,
		.fft_cutoff = 15,
	},
	.max_slice_size = 32769,
	.fits_save_icc = TRUE,
	.icc = {
		.rendering_intent = INTENT_RELATIVE_COLORIMETRIC,
		.proofing_intent = INTENT_RELATIVE_COLORIMETRIC,
		.export_intent = INTENT_RELATIVE_COLORIMETRIC,
		.processing_intent = INTENT_PERCEPTUAL,
		.working_gamut = TYPE_SRGB,
		.icc_path_monitor = NULL,
		.icc_path_soft_proof = NULL,
		.custom_monitor_profile_active = FALSE,
		.soft_proofing_profile_active = FALSE,
		.custom_icc_trc = NULL,
		.custom_icc_gray = NULL,
		.export_8bit_method = EXPORT_SRGB,
		.export_16bit_method = EXPORT_WORKING,
		.default_to_srgb = TRUE,
		.rendering_bpc = TRUE,
		.autoconversion = ICC_NEVER_AUTOCONVERT,
		.autoassignment = ICC_ASSIGN_ON_STRETCH,
		.pedantic_linear = FALSE,
		.cmf = CMF_1931_2DEG
	},
	.spcc = {
		.use_spcc_repository = TRUE,
		.auto_spcc_update = TRUE,
		.redpref = NULL,
		.greenpref = NULL,
		.bluepref = NULL,
		.lpfpref = NULL,
		.oscfilterpref = NULL,
		.monosensorpref = NULL,
		.oscsensorpref = NULL,
		.is_mono = TRUE,
		.is_dslr = FALSE,
		.nb_mode = FALSE,
		.red_wl = 656.28,
		.green_wl = 500.70,
		.blue_wl = 500.70,
		.red_bw = 6.0,
		.green_bw = 6.0,
		.blue_bw = 6.0
	}
};

void free_preferences(preferences *pref) {
	g_free(pref->ext);
	pref->ext = NULL;
	g_free(pref->swap_dir);
	pref->swap_dir = NULL;
	g_free(pref->copyright);
	pref->copyright = NULL;
	g_free(pref->starnet_exe);
	pref->starnet_exe = NULL;
	g_free(pref->graxpert_path);
	pref->graxpert_path = NULL;
	g_free(pref->starnet_weights);
	pref->starnet_weights = NULL;
	g_free(pref->asnet_dir);
	pref->asnet_dir = NULL;
	g_free(pref->lang);
	pref->lang = NULL;
	g_slist_free_full(g_steal_pointer(&pref->gui.script_path), g_free);
	pref->gui.script_path = NULL;
	g_free(pref->fftw_conf.wisdom_file);
	pref->fftw_conf.wisdom_file = NULL;
	g_free(pref->astrometry.default_obscode);
	pref->astrometry.default_obscode = NULL;
}

void set_wisdom_file() {
	if (com.pref.fftw_conf.wisdom_file)
		g_free(com.pref.fftw_conf.wisdom_file);
	if (com.pref.fftw_conf.multithreaded)
		com.pref.fftw_conf.wisdom_file = g_build_filename(g_get_user_cache_dir(), "siril_fftw_threaded.wisdom", NULL);
	else
		com.pref.fftw_conf.wisdom_file = g_build_filename(g_get_user_cache_dir(), "siril_fftw.wisdom", NULL);
}

static void initialize_configurable_colors() {
	com.pref.gui.config_colors.color_bkg_samples = g_strdup("rgba(255, 51, 26, 1.0)");
	com.pref.gui.config_colors.color_std_annotations = g_strdup("rgba(128, 255, 77, 0.9)");
	com.pref.gui.config_colors.color_dso_annotations = g_strdup("rgba(255, 128, 0, 0.9)");
	com.pref.gui.config_colors.color_sso_annotations = g_strdup("rgba(255, 255, 0, 0.9)");
	com.pref.gui.config_colors.color_tmp_annotations = g_strdup("rgba(255, 0, 0, 0.9)");
}

/* static + dynamic settings initialization */
void initialize_default_settings() {
	com.pref = pref_init;
	com.pref.spcc.oscfilterpref = g_strdup("No filter");
	com.pref.ext = g_strdup(".fit");
	com.pref.prepro.stack_default = g_strdup("$seqname$stacked");
	com.pref.swap_dir = g_strdup(g_get_tmp_dir());
	initialize_local_catalogues_paths();
	initialize_configurable_colors();
}

void update_gain_from_gfit() {
	if (gfit.keywords.cvf > 0.0)
		com.pref.phot_set.gain = gfit.keywords.cvf;
}

struct settings_access all_settings[] = {
	{ "core", "wd", STYPE_STRDIR, N_("current working directory"), &com.pref.wd },
	{ "core", "extension", STYPE_STR, N_("FITS file extension"), &com.pref.ext },
	{ "core", "force_16bit", STYPE_BOOL, N_("don't use 32 bits for pixel depth"), &com.pref.force_16bit },
	{ "core", "fits_save_icc", STYPE_BOOL, N_("embed ICC profiles in FITS when saving"), &com.pref.fits_save_icc },
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
	{ "core", "binning_update", STYPE_BOOL, N_("update pixel size of binned images"), &com.pref.binning_update },
	{ "core", "wcs_formalism", STYPE_INT, N_("WCS formalism used in FITS header"), &com.pref.wcs_formalism, { .range_int = { 0, 1 } } },
	{ "core", "catalogue_namedstars", STYPE_STR, N_("Path of the namedstars.dat catalogue"), &com.pref.catalogue_paths[0] },
	{ "core", "catalogue_unnamedstars", STYPE_STR, N_("Path of the unnamedstars.dat catalogue"), &com.pref.catalogue_paths[1] },
	{ "core", "catalogue_tycho2", STYPE_STR, N_("Path of the deepstars.dat catalogue"), &com.pref.catalogue_paths[2] },
	{ "core", "catalogue_nomad", STYPE_STR, N_("Path of the USNO-NOMAD-1e8.dat catalogue"), &com.pref.catalogue_paths[3] },
	{ "core", "catalogue_gaia_astro", STYPE_STR, N_("Path of the local Gaia astrometric catalogue"), &com.pref.catalogue_paths[4] },
	{ "core", "catalogue_gaia_photo", STYPE_STR, N_("Path of the local Gaia photometric catalogue"), &com.pref.catalogue_paths[5] },
	{ "core", "rgb_aladin", STYPE_BOOL, N_("add CTYPE3='RGB' in the FITS header"), &com.pref.rgb_aladin },
	{ "core", "use_checksum", STYPE_BOOL, N_("Verify file checksums if they exist"), &com.pref.use_checksum },
	{ "core", "copyright", STYPE_STR, N_("user copyright to put in file header"), &com.pref.copyright },
	{ "core", "starnet_exe", STYPE_STR, N_("location of the StarNet executable"), &com.pref.starnet_exe },
	{ "core", "starnet_weights", STYPE_STR, N_("location of the StarNet-torch weights file"), &com.pref.starnet_weights },
	{ "core", "graxpert_path", STYPE_STR, N_("location of the GraXpert executable"), &com.pref.graxpert_path },
#ifdef _WIN32
	{ "core", "asnet_dir", STYPE_STR, N_("directory of the asnet_ansvr installation"), &com.pref.asnet_dir },
#else
	{ "core", "asnet_dir", STYPE_STR, N_("directory containing the solve-field executable"), &com.pref.asnet_dir },
#endif
	{ "core", "fftw_timelimit", STYPE_DOUBLE, N_("FFTW planning timelimit"), &com.pref.fftw_conf.timelimit },
	{ "core", "fftw_conv_fft_cutoff", STYPE_INT, N_("Convolution minimum kernel size to use FFTW"), &com.pref.fftw_conf.fft_cutoff },
	{ "core", "fftwf_strategy", STYPE_INT, N_("FFTW planning strategy"), &com.pref.fftw_conf.strategy },
	{ "core", "fftw_multithreaded", STYPE_BOOL, N_("multithreaded FFTW"), &com.pref.fftw_conf.multithreaded },
	{ "core", "max_slice_size", STYPE_INT, N_("Maximum slice size for automated slice processing"), &com.pref.max_slice_size, { .range_int = { 512, 32769 } } },

	{ "starfinder", "focal_length", STYPE_DOUBLE, N_("focal length in mm for radius adjustment"), &com.pref.starfinder_conf.focal_length, { .range_double = { 0., 999999. } } },
	{ "starfinder", "pixel_size", STYPE_DOUBLE, N_("pixel size in µm for radius adjustment"), &com.pref.starfinder_conf.pixel_size_x, { .range_double = { 0., 99. } } },

	{ "debayer", "use_bayer_header", STYPE_BOOL, N_("use pattern from the file header"), &com.pref.debayer.use_bayer_header },
	{ "debayer", "pattern", STYPE_INT, N_("index of the Bayer pattern"), &com.pref.debayer.bayer_pattern, { .range_int = { 0, XTRANS_FILTER_4 } } },
	{ "debayer", "interpolation", STYPE_INT, N_("type of interpolation"), &com.pref.debayer.bayer_inter, { .range_int = { 0, XTRANS } } },
	{ "debayer", "orientation", STYPE_INT, N_("row-order preference"), &com.pref.debayer.orientation, { .range_int = { 0, ROW_ORDER_FORCE_BOTTOMUP } } },
	{ "debayer", "offset_x", STYPE_INT, N_("Bayer matrix offset X"), &com.pref.debayer.xbayeroff, { .range_int = { 0, 1 } } },
	{ "debayer", "offset_y", STYPE_INT, N_("Bayer matrix offset Y"), &com.pref.debayer.ybayeroff, { .range_int = { 0, 1 } } },
	{ "debayer", "xtrans_passes", STYPE_INT, N_("Number of passes for the X-Trans Markesteijn algorithm"), &com.pref.debayer.xtrans_passes, { .range_int = { 1, 4 } } },

	{ "photometry", "gain", STYPE_DOUBLE, N_("electrons per ADU for noise estimation"), &com.pref.phot_set.gain, { .range_double = { 0., 10. } } },
	{ "photometry", "inner", STYPE_DOUBLE, N_("inner radius for background annulus"), &com.pref.phot_set.inner, { .range_double = { 2., 100. } } },
	{ "photometry", "outer", STYPE_DOUBLE, N_("outer radius for background annulus"), &com.pref.phot_set.outer, { .range_double = { 3., 200. } } },
	{ "photometry", "inner_factor", STYPE_DOUBLE, N_("factor for inner radius automatic computation"), &com.pref.phot_set.auto_inner_factor, { .range_double = { 2.0, 50.0 } } },
	{ "photometry", "outer_factor", STYPE_DOUBLE, N_("factor for outer radius automatic computation"), &com.pref.phot_set.auto_outer_factor, { .range_double = { 2.0, 50.0 } } },
	{ "photometry", "force_radius", STYPE_BOOL, N_("force flux aperture value"), &com.pref.phot_set.force_radius },
	{ "photometry", "auto_aperture_factor", STYPE_DOUBLE, N_("Radius/halfFWHM ratio"), &com.pref.phot_set.auto_aperture_factor, { .range_double = { 1., 5. } }  },
	{ "photometry", "aperture", STYPE_DOUBLE, N_("forced aperture for flux computation"), &com.pref.phot_set.aperture, { .range_double = { 1., 100. } } },
	{ "photometry", "minval", STYPE_DOUBLE, N_("minimum valid pixel value for photometry"), &com.pref.phot_set.minval, { .range_double = { -65536.0, 65534.0 } } },
	{ "photometry", "maxval", STYPE_DOUBLE, N_("maximum valid pixel value for photometry"), &com.pref.phot_set.maxval, { .range_double = { 1.0, 65535.0 } } },
	{ "photometry", "discard_var_catalogues", STYPE_INT, N_("catalogues to be used to discard the variable stars from the comparison stars list"), &com.pref.phot_set.discard_var_catalogues, { .range_int = { 0, 7 } } },
	{ "photometry", "redpref", STYPE_STR, N_("preferred SPCC red filter"), &com.pref.spcc.redpref },
	{ "photometry", "greenpref", STYPE_STR, N_("preferred SPCC green filter"), &com.pref.spcc.greenpref },
	{ "photometry", "bluepref", STYPE_STR, N_("preferred SPCC blue filter"), &com.pref.spcc.bluepref },
	{ "photometry", "lpfpref", STYPE_STR, N_("preferred SPCC DSLR LPF filter"), &com.pref.spcc.lpfpref },
	{ "photometry", "oscfilterpref", STYPE_STR, N_("preferred SPCC OSC filter"), &com.pref.spcc.oscfilterpref },
	{ "photometry", "monosensorpref", STYPE_STR, N_("preferred SPCC mono sensor"), &com.pref.spcc.monosensorpref },
	{ "photometry", "oscsensorpref", STYPE_STR, N_("preferred SPCC OSC sensor"), &com.pref.spcc.oscsensorpref },
	{ "photometry", "is_mono", STYPE_BOOL, N_("is the SPCC sensor mono?"), &com.pref.spcc.is_mono },
	{ "photometry", "is_dslr", STYPE_BOOL, N_("is the SPCC OSC sensor a DSLR?"), &com.pref.spcc.is_dslr },
	{ "photometry", "nb_mode", STYPE_BOOL, N_("are we in narrowband mode?"), &com.pref.spcc.nb_mode },
	{ "photometry", "r_wl", STYPE_DOUBLE, N_("red NB filter wavelength"), &com.pref.spcc.red_wl },
	{ "photometry", "r_bw", STYPE_DOUBLE, N_("red NB filter bandwidth"), &com.pref.spcc.red_bw },
	{ "photometry", "g_wl", STYPE_DOUBLE, N_("green NB filter wavelength"), &com.pref.spcc.green_wl },
	{ "photometry", "g_bw", STYPE_DOUBLE, N_("green NB filter bandwidth"), &com.pref.spcc.green_bw },
	{ "photometry", "b_wl", STYPE_DOUBLE, N_("blue NB filter wavelength"), &com.pref.spcc.blue_wl },
	{ "photometry", "b_bw", STYPE_DOUBLE, N_("blue NB filter bandwidth"), &com.pref.spcc.blue_bw },

	{ "astrometry", "asnet_keep_xyls", STYPE_BOOL, N_("do not delete .xyls FITS tables"), &com.pref.astrometry.keep_xyls_files },
	{ "astrometry", "asnet_keep_wcs", STYPE_BOOL, N_("do not delete .wcs result files"), &com.pref.astrometry.keep_wcs_files },
	{ "astrometry", "asnet_show_output", STYPE_BOOL, N_("show solve-field output in main log"), &com.pref.astrometry.show_asnet_output },

	{ "astrometry", "sip_order", STYPE_INT, N_("degrees of the polynomial correction"), &com.pref.astrometry.sip_correction_order, { .range_int = { 1, 5 } } },
	{ "astrometry", "radius", STYPE_DOUBLE, N_("radius around the target coordinates (degrees)"), &com.pref.astrometry.radius_degrees, { .range_double = { 0.01, 30.0 } } },
	{ "astrometry", "max_seconds_run", STYPE_INT, N_("maximum seconds to try solving"), &com.pref.astrometry.max_seconds_run, { .range_int = { 0, 100000 } } },
	{ "astrometry", "update_default_scale", STYPE_BOOL, N_("update default focal length and pixel size from the result"), &com.pref.astrometry.update_default_scale },
	{ "astrometry", "percent_scale_range", STYPE_INT, N_("percent below and above the expected sampling to allow"), &com.pref.astrometry.percent_scale_range, { .range_int = { 1, 50 } } },
	{ "astrometry", "default_obscode", STYPE_STR, N_("default IAU observatory code"), &com.pref.astrometry.default_obscode },

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
	{ "gui_prepro", "disto_lib", STYPE_STR, N_("default distortion master"), &com.pref.prepro.disto_lib },
	{ "gui_prepro", "use_disto_lib", STYPE_BOOL, N_("use default master distortion"), &com.pref.prepro.use_disto_lib },
	{ "gui_prepro", "stack_default", STYPE_STR, N_("default stack name"), &com.pref.prepro.stack_default },
	{ "gui_prepro", "use_stack_default", STYPE_BOOL, N_("use preferred stack name"), &com.pref.prepro.use_stack_default },

	{ "gui_registration", "method", STYPE_INT, N_("index of the selected registration method"), &com.pref.gui.reg_settings, { .range_int = { 0, 7 } } },
	{ "gui_registration", "interpolation", STYPE_INT, N_("index of the selected interpolation method"), &com.pref.gui.reg_interpolation, { .range_int = { 0, 5 } } },
	{ "gui_registration", "clamping", STYPE_BOOL, N_("use clamping method with Lanczos and Cubic interpolation"), &com.pref.gui.reg_clamping },

	{ "gui_stack", "method", STYPE_INT, N_("index of the selected method"), &com.pref.stack.method, { .range_int = { 0, STACK_MIN } } },
	{ "gui_stack", "normalization", STYPE_INT, N_("index of the normalization method"), &com.pref.stack.normalisation_method, { .range_int = { 0, MULTIPLICATIVE_SCALING } } },
	{ "gui_stack", "rejection", STYPE_INT, N_("index of the rejection method"), &com.pref.stack.rej_method, { .range_int = { 0, GESDT } } },
	{ "gui_stack", "weighting", STYPE_INT, N_("index of the weighting method"), &com.pref.stack.weighting_method, { .range_int = { 0, NBSTACK_WEIGHT } } },
	{ "gui_stack", "sigma_low", STYPE_DOUBLE, N_("sigma low value for rejection"), &com.pref.stack.sigma_low, { .range_double = { 0., 20. } } },
	{ "gui_stack", "sigma_high", STYPE_DOUBLE, N_("sigma high value for rejection"), &com.pref.stack.sigma_high, { .range_double = { 0., 20. } } },
	{ "gui_stack", "linear_low", STYPE_DOUBLE, N_("linear low value for rejection"), &com.pref.stack.linear_low, { .range_double = { 0., 20. } } },
	{ "gui_stack", "linear_high", STYPE_DOUBLE, N_("linear high value for rejection"), &com.pref.stack.linear_high, { .range_double = { 0., 20. } } },
	{ "gui_stack", "percentile_low", STYPE_DOUBLE, N_("percentile low value for rejection"), &com.pref.stack.percentile_low, { .range_double = { 0., 100. } } },
	{ "gui_stack", "percentile_high", STYPE_DOUBLE, N_("percentile high value for rejection"), &com.pref.stack.percentile_high, { .range_double = { 0., 100. } } },

	{ "gui", "first_start", STYPE_STR, N_("first start of siril"), &com.pref.gui.first_start },
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
	{ "gui", "use_scripts_repository", STYPE_BOOL, N_("use and sync online scripts repository"), &com.pref.use_scripts_repository },
	{ "gui", "use_spcc_repository", STYPE_BOOL, N_("use and sync spcc-database repository"), &com.pref.spcc.use_spcc_repository },
	{ "gui", "auto_update_scripts", STYPE_BOOL, N_("auto sync online scripts repository"), &com.pref.auto_script_update },
	{ "gui", "auto_update_spcc", STYPE_BOOL, N_("auto sync spcc-database repository"), &com.pref.spcc.auto_spcc_update },
	{ "gui", "selected_scripts", STYPE_STRLIST, N_("list of scripts selected from the repository"), &com.pref.selected_scripts },
	{ "gui", "warn_scripts_run", STYPE_BOOL, N_("warn when launching a script"), &com.pref.gui.warn_scripts_run },
	{ "gui", "show_thumbnails", STYPE_BOOL, N_("show thumbnails in open dialog"), &com.pref.gui.show_thumbnails },
	{ "gui", "thumbnail_size", STYPE_INT, N_("size of the thumbnails"), &com.pref.gui.thumbnail_size },
	{ "gui", "selection_guides", STYPE_INT, N_("number of elements of the grid guides"), &com.pref.gui.selection_guides },
	{ "gui", "show_deciasec", STYPE_BOOL, N_("show tenths of arcseconds on hover"), &com.pref.gui.show_deciasec },
	{ "gui", "default_rendering_mode", STYPE_INT, N_("default display mode"), &com.pref.gui.default_rendering_mode, { .range_int = { 0, 6 } } },
	{ "gui", "display_histogram_mode", STYPE_INT, N_("default histogram display mode"), &com.pref.gui.display_histogram_mode, { .range_int = { 0, 1 } } },
	{ "gui", "roi_mode", STYPE_INT, N_("ROI selection mode"), &com.pref.gui.roi_mode },
	{ "gui", "roi_warning", STYPE_BOOL, N_("enable ROI dialog warning"), &com.pref.gui.enable_roi_warning },
	{ "gui", "mmb_zoom_action", STYPE_INT, N_("Middle mouse button double click zoom action"), &com.pref.gui.mmb_action },
	{ "gui", "mouse_speed_limit", STYPE_DOUBLE, N_("Mouse smooth scroll speed limit"), &com.pref.gui.mouse_speed_limit },
	{ "gui", "color_bkg_samples", STYPE_STR, N_("configure background samples color"), &com.pref.gui.config_colors.color_bkg_samples },
	{ "gui", "color_std_annotations", STYPE_STR, N_("configure standard annotation color"), &com.pref.gui.config_colors.color_std_annotations },
	{ "gui", "color_dso_annotations", STYPE_STR, N_("configure dso annotation color"), &com.pref.gui.config_colors.color_dso_annotations },
	{ "gui", "color_sso_annotations", STYPE_STR, N_("configure sso annotation color"), &com.pref.gui.config_colors.color_sso_annotations },
	{ "gui", "color_tmp_annotations", STYPE_STR, N_("configure tmp annotation color"), &com.pref.gui.config_colors.color_tmp_annotations },
	{ "gui", "custom_monitor_profile", STYPE_STR, N_("path to custom monitor ICC profile"), &com.pref.icc.icc_path_monitor },
	{ "gui", "soft_proofing_profile", STYPE_STR, N_("path to soft proofing ICC profile"), &com.pref.icc.icc_path_soft_proof },
	{ "gui", "icc_custom_monitor_active", STYPE_BOOL, N_("custom monitor profile active"), &com.pref.icc.custom_monitor_profile_active },
	{ "gui", "icc_soft_proofing_active", STYPE_BOOL, N_("output proofing profile active"), &com.pref.icc.soft_proofing_profile_active },
	{ "gui", "custom_RGB_ICC_profile", STYPE_STR, N_("path to custom RGB ICC profile"), &com.pref.icc.custom_icc_trc },
	{ "gui", "custom_gray_ICC_profile", STYPE_STR, N_("path to custom gray ICC profile"), &com.pref.icc.custom_icc_gray },
	{ "gui", "rendering_intent", STYPE_INT, N_("color management rendering intent"), &com.pref.icc.rendering_intent },
	{ "gui", "proofing_intent", STYPE_INT, N_("color management soft proofing intent"), &com.pref.icc.proofing_intent },
	{ "gui", "export_intent", STYPE_INT, N_("color management export intent"), &com.pref.icc.export_intent },
	{ "gui", "default_to_srgb", STYPE_BOOL, N_("default to sRGB when exporting non color managed images"), &com.pref.icc.default_to_srgb },
	{ "gui", "working_gamut", STYPE_INT, N_("color management working gamut"), &com.pref.icc.working_gamut },
	{ "gui", "export_8bit_method", STYPE_INT, N_("color management export profile for 8bit files"), &com.pref.icc.export_8bit_method },
	{ "gui", "export_16bit_method", STYPE_INT, N_("color management export profile for 16bit files"), &com.pref.icc.export_16bit_method },
	{ "gui", "icc_autoconversion", STYPE_INT, N_("autoconvert images with an ICC profile to the working color space"), &com.pref.icc.autoconversion },
	{ "gui", "icc_autoassignment", STYPE_INT, N_("encodes ICC profile auto-assignment options"), &com.pref.icc.autoassignment },
	{ "gui", "icc_rendering_bpc", STYPE_BOOL, N_("enable rendering BPC"), &com.pref.icc.rendering_bpc },
	{ "gui", "icc_pedantic_linear", STYPE_BOOL, N_("pedantically assign linear ICC profiles"), &com.pref.icc.pedantic_linear },
	{ "gui", "mouse_actions", STYPE_STRLIST, N_("mouse actions config. Opaque data structure, edit this using the GUI"), &com.pref.gui.mouse_cfg.mouse_actions_array },
	{ "gui", "scroll_actions", STYPE_STRLIST, N_("scroll actions config. Opaque data structure, edit this using the GUI"), &com.pref.gui.mouse_cfg.scroll_actions_array },
	{ "gui_astrometry", "compass_position", STYPE_INT, N_("index of the compass position over grid"), &com.pref.gui.position_compass, { .range_int = { 0, 5 } } },
	{ "gui_astrometry", "cat_messier", STYPE_BOOL, N_("show Messier objects in annotations"), &com.pref.gui.catalog[0] },
	{ "gui_astrometry", "cat_ngc", STYPE_BOOL, N_("show NGC objects in annotations"), &com.pref.gui.catalog[1] },
	{ "gui_astrometry", "cat_ic", STYPE_BOOL, N_("show IC objects in annotations"), &com.pref.gui.catalog[2] },
	{ "gui_astrometry", "cat_ldn", STYPE_BOOL, N_("show LDN objects in annotations"), &com.pref.gui.catalog[3] },
	{ "gui_astrometry", "cat_sh2", STYPE_BOOL, N_("show SH2 objects in annotations"), &com.pref.gui.catalog[4] },
	{ "gui_astrometry", "cat_stars", STYPE_BOOL, N_("show stars in annotations"), &com.pref.gui.catalog[5] },
	{ "gui_astrometry", "cat_const", STYPE_BOOL, N_("show constellations in annotations"), &com.pref.gui.catalog[6] },
	{ "gui_astrometry", "cat_const_names", STYPE_BOOL, N_("show constellations names in annotations"), &com.pref.gui.catalog[7] },
	{ "gui_astrometry", "cat_user_dso", STYPE_BOOL, N_("show user DSO objects in annotations"), &com.pref.gui.catalog[8] },
	{ "gui_astrometry", "cat_user_sso", STYPE_BOOL, N_("show user SSO objects in annotations"), &com.pref.gui.catalog[9] },

	{ "gui_pixelmath", "pm_presets", STYPE_STRLIST, N_("list of pixel math presets"), &com.pref.gui.pm_presets },

	{ "script_editor", "highlight_syntax", STYPE_BOOL, N_("highlight syntax in the script editor"), &com.pref.gui.editor_cfg.highlight_syntax },
	{ "script_editor", "highlight_bracketmatch", STYPE_BOOL, N_("highlight matching brackets in the script editor"), &com.pref.gui.editor_cfg.highlight_bracketmatch },
	{ "script_editor", "rmargin", STYPE_BOOL, N_("show the right margin in the script editor"), &com.pref.gui.editor_cfg.rmargin },
	{ "script_editor", "rmargin_pos", STYPE_INT, N_("position of the right margin in the script editor"), &com.pref.gui.editor_cfg.rmargin_pos },
	{ "script_editor", "show_linenums", STYPE_BOOL, N_("show line numbers in the script editor"), &com.pref.gui.editor_cfg.show_linenums },
	{ "script_editor", "show_linemarks", STYPE_BOOL, N_("show line marks in the script editor"), &com.pref.gui.editor_cfg.show_linemarks },
	{ "script_editor", "highlight_currentline", STYPE_BOOL, N_("highlight the current line in the script editor"), &com.pref.gui.editor_cfg.highlight_currentline },
	{ "script_editor", "autoindent", STYPE_BOOL, N_("automatically indent new lines"), &com.pref.gui.editor_cfg.autoindent },
	{ "script_editor", "indentontab", STYPE_BOOL, N_("indent selected blocks of lines in the script editor using the tab key"), &com.pref.gui.editor_cfg.indentontab },
	{ "script_editor", "smartbs", STYPE_BOOL, N_("Smart Backspace behaviour in the script editor"), &com.pref.gui.editor_cfg.smartbs },
	{ "script_editor", "smarthomeend", STYPE_BOOL, N_("Smart Home / End behaviour in the script editor"), &com.pref.gui.editor_cfg.smarthomeend },
	{ "script_editor", "showspaces", STYPE_BOOL, N_("Show visible space and tab characters in the script editor"), &com.pref.gui.editor_cfg.showspaces },
	{ "script_editor", "shownewlines", STYPE_BOOL, N_("Show visible newline characters in the script editor"), &com.pref.gui.editor_cfg.shownewlines },
	{ "script_editor", "minimap", STYPE_BOOL, N_("Show a minimap in the script editor"), &com.pref.gui.editor_cfg.minimap },

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

gchar* get_settings_key(const char *group, const char *key, gboolean with_details) {
	struct settings_access *desc = get_key_settings(group, key);
	if (!desc) {
		siril_log_message(_("Unknown settings variable %s.%s\n"), group, key);
		return NULL;
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
		g_string_append_printf(str, ", %s", _(desc->desc));
	}
	gchar *s = g_string_free(str, FALSE);
	return s;
}

int print_settings_key(const char *group, const char *key, gboolean with_details) {
	gchar *s = get_settings_key(group, key, with_details);
	if (s) {
		siril_log_message("%s\n", s);
		g_free(s);
	}
	return 0;
}

int print_all_settings(gboolean with_details) {
	int nb_settings = sizeof(all_settings) / sizeof(struct settings_access) - 1;
	for (int i = 0; i < nb_settings; i++) {
		print_settings_key(all_settings[i].group, all_settings[i].key, with_details);
	}
	return 0;
}
