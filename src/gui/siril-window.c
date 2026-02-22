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
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include <gtk/gtk.h>
#include "core/siril_actions.h"

static GActionEntry win_entries[] = {
	{ "close", close_action_activate },
	{ "undo", undo_action_activate },
	{ "redo", redo_action_activate },
	{ "documentation", doc_action_activate },
	{ "updates", updates_action_activate },
	{ "full-screen", full_screen_activated},
	{ "hide-show-toolbar", toolbar_activate },
	{ "shortcuts", keyboard_shortcuts_activated},
	{ "cwd", cwd_action_activate },
	{ "livestacking", livestacking_action_activate },
	{ "panel", panel_activate },

	{ "conversion", tab_conversion_activate },
	{ "sequence", tab_sequence_activate },
	{ "registration", tab_registration_activate },
	{ "prepro", tab_prepro_activate },
	{ "plot", tab_plot_activate },
	{ "stacking", tab_stacking_activate },
	{ "logs", tab_logs_activate }
};

static GActionEntry image_entries[] = {
	{ "bit-depth", NULL },
	{ "zoom-out", zoom_out_activate },
	{ "zoom-in", zoom_in_activate },
	{ "histo_display", on_histogram_overlay_activate, NULL, "false", on_histogram_overlay_state },
	{ "zoom-fit", zoom_fit_activate, NULL, "true", change_zoom_fit_state },
	{ "zoom-one", zoom_one_activate },
	{ "negative-view", negative_view_activate, NULL, "false", negative_view_state },
	{ "color-map", color_map_activate, NULL, "false", color_map_state },
	{ "chain-chan", chain_channels_activate, NULL, "true", chain_channels_state_change },
	{ "snapshot", snapshot_action_activate },
	{ "clipboard", clipboard_action_activate },
	{ "fits-header", image_fits_header_activate },
	{ "statistics", statistics_activate },
	{ "evaluate-noise", noise_activate },
	{ "ccd-inspector", ccd_inspector_activate },
	{ "show-tilt", show_tilt_activate, NULL, "false", show_tilt_state },
	{ "show-disto", show_disto_activate, NULL, "false", show_disto_state },
	{ "astrometry", astrometry_activate },
	{ "photometry", photometry_activate, NULL, "false", photometry_state },
	{ "cut", cut_activate },
	{ "image-information", image_information_activate },
	{ "dyn-psf", dyn_psf_activate },
	{ "annotate-object", annotate_object_activate, NULL, "false", annotate_object_state },
	{ "wcs-grid", wcs_grid_activate, NULL, "false", wcs_grid_state },
	{ "annotate-dialog", annotate_dialog_activate },
	{ "seq-list", seq_list_activate },
	{ "regframe", regframe_activate , NULL, "true", regframe_state },
	{ "nina_light_curve", nina_lc_activate },
	{ "compstars", compstars_activate }
};

static GActionEntry selection_entries[] = {
	{ "pickstar", pick_star_activate },
	{ "psf", psf_activate },
	{ "crop", crop_activate },
	{ "set_roi", set_roi }
};

static GActionEntry selection_sequence_entries[] = {
	{ "seq-psf", seq_psf_activate },
	{ "seq-crop", seq_crop_activate }
};

static GActionEntry rgb_processing_entries[] = {
	{ "remove-green-processing", remove_green_activate },
	{ "saturation-processing", saturation_activate },
	{ "color-calib-processing", color_calib_activate },
	{ "pcc-processing", pcc_activate },
	{ "spcc-processing", spcc_activate },
	{ "align-dft", align_dft_activate },
	{ "align-psf", align_psf_activate },
	{ "align-global", align_global_activate },
	{ "align-kombat", align_kombat_activate },
	{ "split-channel-processing", split_channel_activate },
	{ "ccm-processing", ccm_activate },
	{ "unpurple-processing", unpurple_activate }
};

static GActionEntry any_processing_entries[] = {
	{ "negative-processing", negative_activate },
	{ "deconvolution-processing", deconvolution_activate },
	{ "histo-processing", histo_activate },
	{ "curves-processing", curves_activate },
	{ "payne-processing", payne_activate },
	{ "starnet-processing", starnet_activate },
	{ "fix-banding-processing", fix_banding_activate },
	{ "cosmetic-processing", cosmetic_activate },
	{ "background-extr-processing", background_extr_activate },
	{ "icc-tool", icc_activate },
	{ "clear_roi", clear_roi },
	{ "mask_from_color", mask_from_color_activate },
	{ "mask_from_image", mask_from_image_activate },
	{ "mask_from_stars", mask_from_stars_activate },
	{ "mask_from_file", mask_from_file_activate },
	{ "mask_from_gradient", mask_from_gradient_activate },
	{ "mask_add_from_poly", mask_add_from_poly_activate },
	{ "mask_clear_from_poly", mask_clear_from_poly_activate },
	{ "clear_mask", clear_mask_activate },
	{ "autostretch_mask", autostretch_mask_activate },
	{ "blur_mask", blur_mask_activate },
	{ "feather_mask", feather_mask_activate },
	{ "invert_mask", invert_mask_activate },
	{ "scale_mask", mask_scale_activate },
	{ "threshold_mask", threshold_mask_activate }
};

static GActionEntry any_mono_processing_entries[] = {
	{ "split-cfa-processing", split_cfa_activate }
};

static GActionEntry single_processing_entries[] = {
	{ "asinh-processing", asinh_activate },
	{ "epf-processing", epf_activate },
	{ "denoise-processing", denoise_activate },
	{ "binning-processing", binning_activate },
	{ "resample-processing", resample_activate },
	{ "rotation-processing", rotation_activate },
	{ "rotation90-processing", rotation90_activate },
	{ "rotation270-processing", rotation270_activate },
	{ "mirrorx-processing", mirrorx_activate },
	{ "mirrory-processing", mirrory_activate },
	{ "wavelets-processing", wavelets_activate },
	{ "split-wavelets-processing", split_wavelets_activate },
	{ "medianfilter-processing", medianfilter_activate },
	{ "rgradient-processing", rgradient_activate },
	{ "clahe-processing", clahe_activate },
	{ "linearmatch-processing", linearmatch_activate },
	{ "star-desaturate", star_desaturate_activate },
	{ "star-synthetic", star_synthetic_activate }
};

static GActionEntry none_processing_entries[] = {
	{ "fft-processing", fft_activate },
	{ "rgb-compositing-processing", rgb_compositing_activate },
	{ "star-remix-processing", star_remix_activate },
	{ "merge-cfa-processing", merge_cfa_activate },
	{ "pixel-math", pixel_math_activate }
};

static void _siril_window_enable_action_group(GActionMap *map,
		const gchar **group, gboolean enable) {
	GAction *action;

	for (const gchar **it = group; *it != NULL; it++) {
		action = g_action_map_lookup_action(map, *it);
		if (G_LIKELY(action))
			g_simple_action_set_enabled(G_SIMPLE_ACTION(action), enable);
		else
			g_warning("Action not found in action group: %s", *it);
	}
}

void siril_window_enable_image_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *image_actions[] = {
		"bit-depth",
		"zoom-out",
		"zoom-in",
		"zoom-fit",
		"zoom-one",
		"negative-view",
		"color-map",
		"snapshot",
		"clipboard",
		"statistics",
		"evaluate-noise",
		"ccd-inspector",
		"show-tilt",
		"astrometry",
		"photometry",
		"image-information",
		"dyn-psf",
		"seq-list",
		"regframe",
		"cut",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), image_actions, enable);
}

void siril_window_enable_wcs_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *wcs_processing_actions[] = {
		"annotate-object",
		"wcs-grid",
		"compstars",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), wcs_processing_actions, enable);
}

void siril_window_enable_wcs_disto_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *wcs_disto_processing_actions[] = {
		"show-disto",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), wcs_disto_processing_actions, enable);
}

void siril_window_autostretch_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *image_actions[] = {
		"chain-chan",
		NULL
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), image_actions, enable);
}

void siril_window_enable_rgb_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *rgb_processing_actions[] = {
		"remove-green-processing",
		"saturation-processing",
		"color-calib-processing",
		"split-channel-processing",
		"align-global",
		"unpurple-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), rgb_processing_actions, enable);
}

void siril_window_enable_any_rgb_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *any_rgb_processing_actions[] = {
		"ccm-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), any_rgb_processing_actions, enable);
}

void siril_window_enable_rgb_wcs_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *rgb_wcs_processing_actions[] = {
		"pcc-processing",
		"spcc-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), rgb_wcs_processing_actions, enable);
}

void siril_window_enable_any_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *any_processing_actions[] = {
		"negative-processing",
		"deconvolution-processing",
		"histo-processing",
		"payne-processing",
		"curves-processing",
		"fix-banding-processing",
		"starnet-processing",
		"cosmetic-processing",
		"background-extr-processing",
		"icc-tool",
		"clear_roi",
		"mask_from_image",
		"mask_from_stars",
		"mask_from_file",
		"clear_mask",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), any_processing_actions, enable);
}

void siril_window_enable_any_mono_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *any_mono_processing_actions[] = {
		"split-cfa-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), any_mono_processing_actions, enable);
}

void siril_window_enable_single_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *single_processing_actions[] = {
		"asinh-processing",
		"epf-processing",
		"denoise-processing",
		"resample-processing",
		"binning-processing",
		"rotation-processing",
		"rotation90-processing",
		"rotation270-processing",
		"mirrorx-processing",
		"mirrory-processing",
		"wavelets-processing",
		"split-wavelets-processing",
		"medianfilter-processing",
		"rgradient-processing",
		"clahe-processing",
		"linearmatch-processing",
		"star-desaturate",
		"star-synthetic",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), single_processing_actions, enable);
}

void siril_window_enable_none_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *none_processing_actions[] = {
		"fft-processing",
		"rgb-compositing-processing",
		"star-remix-processing",
		"merge-cfa-processing",
		"pixel-math",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), none_processing_actions, enable);
}

void siril_window_enable_if_selection_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *selection_actions[] = {
		"pickstar",
		"psf",
		"crop",
		"set_roi",

		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), selection_actions, enable);
}

void siril_window_enable_if_selection_rgb_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *selection_rgb_actions[] = {
		"align-dft",
		"align-psf",
		"align-kombat",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), selection_rgb_actions, enable);
}

void siril_window_enable_if_selection_sequence_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *selection_sequence_actions[] = {
		"seq-psf",
		"seq-crop",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), selection_sequence_actions, enable);
}

void siril_window_map_actions(GtkApplicationWindow *window) {
	g_action_map_add_action_entries(G_ACTION_MAP(window), win_entries, G_N_ELEMENTS(win_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), image_entries, G_N_ELEMENTS(image_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), rgb_processing_entries, G_N_ELEMENTS(rgb_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), any_processing_entries, G_N_ELEMENTS(any_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), any_mono_processing_entries, G_N_ELEMENTS(any_mono_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), single_processing_entries, G_N_ELEMENTS(single_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), none_processing_entries, G_N_ELEMENTS(none_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), selection_entries, G_N_ELEMENTS(selection_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), selection_sequence_entries, G_N_ELEMENTS(selection_sequence_entries), window);
}
