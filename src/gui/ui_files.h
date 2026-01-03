#ifndef UI_FILES_H
#define UI_FILES_H

/***************************************************************
 *
 * Add UI files here as needed. Widgets with dependencies
 * must appear after the GObjects they depend on, otherwise
 * the GtkBuilder will fail and the program will refuse to
 * start.
 * * Dependencies come first (file filters, menus, shortcuts)
 * * Then dialogs should mostly be independent of each other
 * and should be added preserving alphabetical order.
 * * Then comes the main application window.
 * * Finally there are a couple of dialogs that depend on the
 * application window, so they come last of all.
 *
 **************************************************************/

// This list contains most of the UI files. However any files with
// shortcuts that use the primary modifier will be written by Glade
// as GDK_CONTROL_MASK. These require rewriting on MacOS and are
// thus handled separately, and should be added to
// ui_files_with_primary_accelerator below.

const char* ui_files[] = {
// Dependencies first
	"/org/siril/ui/filefilters.ui",
	"/org/siril/ui/siril-shortcuts.ui",
// Dialogs (sorted)
	"/org/siril/ui/aavso_dialog.ui",
	"/org/siril/ui/aberration_inspector.ui",
	"/org/siril/ui/annotate_dialog.ui",
	"/org/siril/ui/asinh_dialog.ui",
	"/org/siril/ui/astrometry_dialog.ui",
	"/org/siril/ui/background_extraction_dialog.ui",
	"/org/siril/ui/bdeconv_dialog.ui",
	"/org/siril/ui/binxy_dialog.ui",
	"/org/siril/ui/canon_fixbanding_dialog.ui",
	"/org/siril/ui/ccm_dialog.ui",
	"/org/siril/ui/CLAHE_dialog.ui",
	"/org/siril/ui/color_calibration.ui",
	"/org/siril/ui/colorchooserdialog.ui",
	"/org/siril/ui/composition_dialog.ui",
	"/org/siril/ui/cosmetic_dialog.ui",
	"/org/siril/ui/crop_dialog.ui",
	"/org/siril/ui/curves_dialog.ui",
	"/org/siril/ui/cut_coords_dialog.ui",
	"/org/siril/ui/cut_dialog.ui",
	"/org/siril/ui/cut_spectroscopy_dialog.ui",
	"/org/siril/ui/data_dialog.ui",
	"/org/siril/ui/denoise_dialog.ui",
	"/org/siril/ui/dialog_FFT.ui",
	"/org/siril/ui/dialog_star_remix.ui",
	"/org/siril/ui/epf_dialog.ui",
	"/org/siril/ui/extract_channel_dialog.ui",
	"/org/siril/ui/extract_wavelet_layers_dialog.ui",
	"/org/siril/ui/file_information_dialog.ui",
	"/org/siril/ui/histogram_dialog.ui",
	"/org/siril/ui/keywords_dialog.ui",
	"/org/siril/ui/icc_dialog.ui",
	"/org/siril/ui/icc_gamut_dialog.ui",
	"/org/siril/ui/linearmatch_dialog.ui",
	"/org/siril/ui/livestacking_player.ui",
	"/org/siril/ui/mask_from_image.ui",
	"/org/siril/ui/mask_from_stars.ui",
	"/org/siril/ui/Median_dialog.ui",
	"/org/siril/ui/merge_cfa_dialog.ui",
	"/org/siril/ui/mouse_actions.ui",
	"/org/siril/ui/pixelmath.ui",
	"/org/siril/ui/resample_dialog.ui",
	"/org/siril/ui/rgradient_dialog.ui",
	"/org/siril/ui/rotation_dialog.ui",
	"/org/siril/ui/s_pcc_dialog.ui",
	"/org/siril/ui/satu_dialog.ui",
	"/org/siril/ui/savepopup.ui",
	"/org/siril/ui/SCNR_dialog.ui",
	"/org/siril/ui/seqlist_dialog.ui",
	"/org/siril/ui/settings_window.ui",
	"/org/siril/ui/siril.ui",
	"/org/siril/ui/spcc_details.ui",
	"/org/siril/ui/split_cfa_dialog.ui",
	"/org/siril/ui/starnet_dialog.ui",
	"/org/siril/ui/stars_list_window.ui",
	"/org/siril/ui/StatWindow.ui",
	"/org/siril/ui/unpurple_dialog.ui",
	"/org/siril/ui/wavelets_dialog.ui",
	// Must be terminated by a NULL string
	""
};

// UI files in this array will be processed differently and on MacOS
// GTK_CONTROL_MASK will be replaced with GDK_META_MASK prior to
// adding to the builder, in order to obtain the expected Cmd- shortcuts
const char* ui_files_with_primary_accelerator[] = {
	"/org/siril/ui/python_scratchpad.ui",
	// Must be terminated by a NULL string
	""
};

#endif
