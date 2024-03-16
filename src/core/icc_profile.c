/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include <glib.h>
#include "algos/lcms_acceleration/lcms2_fast_float.h"
#include "algos/lcms_acceleration/lcms2_threaded.h"
#include "core/siril.h"
#include "algos/colors.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/undo.h"
#include "icc_profile.h"
#include "icc_default_profiles.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/image_interactions.h"
#include "gui/siril-window.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "gui/siril_plot.h"
#include "gui/siril_preview.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/siril_plot.h"
#include "core/siril_log.h"
#include "core/siril_app_dirs.h"
#include "core/proto.h"

// For the log message about JPEG ICC profile support at startup
#ifdef HAVE_LIBJPEG
#include <jconfig.h>
#endif

static GMutex monitor_profile_mutex;
static GMutex soft_proof_profile_mutex;
static GMutex default_profiles_mutex;
static GMutex display_transform_mutex;

static cmsHPROFILE target = NULL; // Target profile for the GUI tool

////// Functions //////

static gchar *
	siril_color_profile_get_info (cmsHPROFILE profile,
								cmsInfoType info) {
	cmsUInt32Number  size;
	gchar           *text = NULL;

	size = cmsGetProfileInfoASCII (profile, info,
									"en", "US", NULL, 0);
	if (size > 0) {
		gchar *data = g_new (gchar, size + 1);

		size = cmsGetProfileInfoASCII (profile, info,
										"en", "US", data, size);
		if (size > 0)
			text = siril_any_to_utf8 (data, -1, NULL);

		g_free (data);
	}
  return text;
}

static gchar *
	siril_color_profile_get_copyright (cmsHPROFILE profile) {
	gchar *string = siril_color_profile_get_info (profile, cmsInfoCopyright);
	return string;
}

// This function is required outside this file
gchar *
	siril_color_profile_get_description (cmsHPROFILE profile) {
	gchar *string = siril_color_profile_get_info (profile, cmsInfoDescription);
	return string;
}

static gchar *
	siril_color_profile_get_manufacturer (cmsHPROFILE profile) {
	gchar *string = siril_color_profile_get_info (profile, cmsInfoManufacturer);
	return string;
}

static gchar *
	siril_color_profile_get_model (cmsHPROFILE profile) {
	gchar *string = siril_color_profile_get_info (profile, cmsInfoModel);
	return string;
}

void icc_profile_set_tag (cmsHPROFILE profile,
                                  cmsTagSignature sig,
                                  const gchar *tag) {
	cmsMLU *mlu;

	mlu = cmsMLUalloc (NULL, 1);
	cmsMLUsetASCII (mlu, "en", "US", tag);
	cmsWriteTag (profile, sig, mlu);
	cmsMLUfree (mlu);
}

cmsHPROFILE srgb_linear() {
	return cmsOpenProfileFromMem(sRGB_elle_V4_g10_icc, sRGB_elle_V4_g10_icc_len);
}
cmsHPROFILE srgb_trc() {
	return cmsOpenProfileFromMem(sRGB_elle_V4_srgbtrc_icc, sRGB_elle_V4_srgbtrc_icc_len);
}
static cmsHPROFILE srgb_trcv2() {
	return cmsOpenProfileFromMem(sRGB_elle_V2_srgbtrc_icc, sRGB_elle_V2_srgbtrc_icc_len);
}

cmsHPROFILE rec2020_linear() {
	return cmsOpenProfileFromMem(Rec2020_elle_V4_g10_icc, Rec2020_elle_V4_g10_icc_len);
}
cmsHPROFILE rec2020_trc() {
	return cmsOpenProfileFromMem(Rec2020_elle_V4_rec709_icc, Rec2020_elle_V4_rec709_icc_len);
}
static cmsHPROFILE rec2020_trcv2() {
	return cmsOpenProfileFromMem(Rec2020_elle_V2_rec709_icc, Rec2020_elle_V2_rec709_icc_len);
}

cmsHPROFILE gray_linear() {
	return cmsOpenProfileFromMem(Gray_elle_V4_g10_icc, Gray_elle_V4_g10_icc_len);
}
cmsHPROFILE gray_srgbtrc() {
	return cmsOpenProfileFromMem(Gray_elle_V4_srgbtrc_icc, Gray_elle_V4_srgbtrc_icc_len);
}
cmsHPROFILE gray_rec709trc() {
	return cmsOpenProfileFromMem(Gray_elle_V4_rec709_icc, Gray_elle_V4_rec709_icc_len);
}
static cmsHPROFILE gray_srgbtrcv2() {
	return cmsOpenProfileFromMem(Gray_elle_V2_srgbtrc_icc, Gray_elle_V2_srgbtrc_icc_len);
}
static cmsHPROFILE gray_rec709trcv2() {
	return cmsOpenProfileFromMem(Gray_elle_V2_rec709_icc, Gray_elle_V2_rec709_icc_len);
}

static cmsHPROFILE srgb_monitor_perceptual() {
	return cmsOpenProfileFromMem(sRGB_v4_ICC_preference_icc, sRGB_v4_ICC_preference_icc_len);
}

void color_manage(fits *fit, gboolean active) {
	fit->color_managed = active;
	if (fit == &gfit && !com.headless) {
		gchar *buffer = NULL, *monitor = NULL, *proof = NULL;
		gchar *name = g_build_filename("/org/siril/ui/", "pixmaps", active ? "color_management.svg" : "color_management_off.svg", NULL);
		gchar *tooltip = NULL;
		if (active) {
			if (fit->icc_profile) {
				buffer = siril_color_profile_get_description(fit->icc_profile);
				monitor = siril_color_profile_get_description(gui.icc.monitor);
			}
			if (gui.icc.soft_proof)
				proof = siril_color_profile_get_description(gui.icc.soft_proof);
			else
				proof = g_strdup(monitor);

			tooltip = g_strdup_printf(_("Image is color managed\nImage profile: %s\nMonitor profile: %s\nSoft proofing profile: %s"), buffer, monitor, proof);
			if (!tooltip)
				tooltip = g_strdup(_("Image is color managed\n\nLeft click: Color management dialog\nRight click: toggle ISO12646 color assessment mode"));
		} else {
			tooltip = g_strdup(_("Image is not color managed\n\nLeft click: Color management dialog\nRight click: toggle ISO12646 color assessment mode"));
		}
		GtkWidget *image = lookup_widget("color_managed_icon");
		GtkWidget *button = lookup_widget("icc_main_window_button");
		gtk_image_set_from_resource((GtkImage*) image, name);
		gtk_widget_set_tooltip_text(button, tooltip);
		g_free(name);
		g_free(buffer);
		g_free(monitor);
		g_free(proof);
		g_free(tooltip);
	}
}

static void export_profile(cmsHPROFILE profile, const char *provided_filename) {
	char *filename = NULL, *path = NULL;
	if (provided_filename != NULL && provided_filename[0] != '\0') {
		filename = strdup(provided_filename);
	} else {
		filename = siril_color_profile_get_description(profile);
		if (!g_str_has_suffix(filename, ".icc")) {
			gchar* temp = g_strdup_printf("%s.icc", filename);
			g_free(filename);
			filename = strdup(temp);
			g_free(temp);
		}
	}
	path = g_build_filename(com.wd, filename, NULL);
	free(filename);
	if (cmsSaveProfileToFile(profile, path)) {
		siril_log_color_message(_("Exported ICC profile to %s\n"), "green", path);
	} else {
		siril_log_color_message(_("Failed to export ICC profile to %s\n"), "red", path);
	}
	free(path);
}

static void export_elle_stone_profiles() {
	cmsHPROFILE profile;
	control_window_switch_to_tab(OUTPUT_LOGS);

	profile = srgb_linear();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = srgb_trc();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = srgb_trcv2();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = rec2020_linear();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = rec2020_trc();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = rec2020_trcv2();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_linear();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_srgbtrc();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_srgbtrcv2();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_rec709trc();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	profile = gray_rec709trcv2();
	export_profile(profile, NULL);
	cmsCloseProfile(profile);

	// Also we use this function to export the ICC sRGB perceptual display profile,
	// but use a confirmation dialog to confirm the ICC's terms of use.
	if (siril_confirm_dialog(_("Terms of Use"), _("To export a copy of the sRGB "
			"monitor profile with perceptual intent tables, please accept the ICC's terms "
			"of use:\n\nTo anyone who acknowledges that the file \"sRGB_v4_ICC_preference.icc\" "
			"is provided \"AS IS\" WITH NO EXPRESS OR IMPLIED WARRANTY, permission to use, "
			"copy and distribute this file for any purpose is hereby granted without fee, "
			"provided that the file is not changed including the ICC copyright notice tag, "
			"and that the name of ICC shall not be used in advertising or publicity "
			"pertaining to distribution of the software without specific, written prior "
			"permission. ICC makes no representations about the suitability of this software "
			"for any purpose."), _("Accept"))) {
		profile = srgb_monitor_perceptual();
		export_profile(profile, "sRGB_v4_ICC_preference.icc");
		cmsCloseProfile(profile);
	}
}

//Two functions to check if profiles are RGB or Gray
static gboolean siril_color_profile_is_rgb(cmsHPROFILE profile) {
	return (cmsGetColorSpace (profile) == cmsSigRgbData);
}

void lock_display_transform() {
	g_mutex_lock(&display_transform_mutex);
}
void unlock_display_transform() {
	g_mutex_unlock(&display_transform_mutex);
}

// This must be locked by the display_transform_mutex, but it is done from
// remap_all_vports() so the mutex lock covers all 3 calls to this function
void display_index_transform(BYTE* index, int vport) {
	BYTE buf[3 * (USHRT_MAX + 1)] = { 0 };
	BYTE* chan = &buf[0] + (vport * (USHRT_MAX + 1));
	memcpy(chan, index, USHRT_MAX + 1);
	cmsDoTransformLineStride(gui.icc.proofing_transform, &buf, &buf, USHRT_MAX + 1, 1, (USHRT_MAX + 1) * 3, (USHRT_MAX + 1) * 3, USHRT_MAX + 1, USHRT_MAX + 1);
	memcpy(index, chan, USHRT_MAX + 1);
}

static gboolean siril_color_profile_is_gray(cmsHPROFILE profile) {
	return (cmsGetColorSpace (profile) == cmsSigGrayData);
}

/* Check if a provided filename is non-null and points to a file that exists */
static gboolean validate_profile(gchar* filename) {
	if (!filename)
		return FALSE;
	if (filename[0] == '\0')
		return FALSE;
	if (!g_file_test(filename, G_FILE_TEST_EXISTS))
		return FALSE;
	return TRUE;
}

// Compares the colorant primaries of 2 profiles to see if they are the same
// Returns TRUE if they are the same or FALSE otherwise
// An optional third profile can be used to check soft proofing mode too
gboolean same_primaries(cmsHPROFILE a, cmsHPROFILE b, cmsHPROFILE c) {
	if (!(a && b))
		return FALSE;
	if (!(cmsIsTag(a, cmsSigRedColorantTag) && cmsIsTag(b, cmsSigRedColorantTag) &&
		cmsIsTag(a, cmsSigGreenColorantTag) && cmsIsTag(b, cmsSigGreenColorantTag) &&
		cmsIsTag(a, cmsSigBlueColorantTag) && cmsIsTag(b, cmsSigBlueColorantTag) &&
		((com.pref.icc.proofing_intent != INTENT_ABSOLUTE_COLORIMETRIC) | ((com.pref.icc.proofing_intent == INTENT_ABSOLUTE_COLORIMETRIC) && (cmsIsTag(a, cmsSigMediaWhitePointTag) && cmsIsTag(b, cmsSigMediaWhitePointTag))))))
		return FALSE; // Can't tell, so err on safe side and return FALSE
	cmsCIEXYZ *a_r = cmsReadTag(a, cmsSigRedColorantTag);
	cmsCIEXYZ *b_r = cmsReadTag(b, cmsSigRedColorantTag);
	cmsCIEXYZ *a_g = cmsReadTag(a, cmsSigGreenColorantTag);
	cmsCIEXYZ *b_g = cmsReadTag(b, cmsSigGreenColorantTag);
	cmsCIEXYZ *a_b = cmsReadTag(a, cmsSigBlueColorantTag);
	cmsCIEXYZ *b_b = cmsReadTag(b, cmsSigBlueColorantTag);
	cmsCIEXYZ *a_w = cmsReadTag(a, cmsSigMediaWhitePointTag);
	cmsCIEXYZ *b_w = cmsReadTag(b, cmsSigMediaWhitePointTag);
	if (memcmp(a_r, b_r, sizeof(cmsCIEXYZ)) ||
			memcmp(a_g, b_g, sizeof(cmsCIEXYZ)) ||
			memcmp(a_b, b_b, sizeof(cmsCIEXYZ)) ||
			((com.pref.icc.proofing_intent == INTENT_ABSOLUTE_COLORIMETRIC) && memcmp(a_w, b_w, sizeof(cmsCIEXYZ))))
		return FALSE;
	if (c) {
		if (!(cmsIsTag(c, cmsSigRedColorantTag) &&
				cmsIsTag(c, cmsSigGreenColorantTag) &&
				cmsIsTag(c, cmsSigBlueColorantTag)) &&
				((com.pref.icc.proofing_intent != INTENT_ABSOLUTE_COLORIMETRIC) | ((com.pref.icc.proofing_intent == INTENT_ABSOLUTE_COLORIMETRIC) && cmsIsTag(c, cmsSigMediaWhitePointTag))))
		return FALSE; // Can't tell, so err on safe side and return FALSE
		cmsCIEXYZ *c_r = cmsReadTag(c, cmsSigRedColorantTag);
		cmsCIEXYZ *c_g = cmsReadTag(c, cmsSigGreenColorantTag);
		cmsCIEXYZ *c_b = cmsReadTag(c, cmsSigBlueColorantTag);
		cmsCIEXYZ *c_w = cmsReadTag(c, cmsSigMediaWhitePointTag);
		if (memcmp(a_r, c_r, sizeof(cmsCIEXYZ)) ||
				memcmp(a_g, c_g, sizeof(cmsCIEXYZ)) ||
				memcmp(a_b, c_b, sizeof(cmsCIEXYZ)) ||
				((com.pref.icc.proofing_intent == INTENT_ABSOLUTE_COLORIMETRIC) && memcmp(a_w, c_w, sizeof(cmsCIEXYZ))))
			return FALSE;
	}
	siril_debug_print("Primaries are the same\n");
	return TRUE;
}

void reset_icc_transforms() {
	g_mutex_lock(&display_transform_mutex);
	if (gfit.color_managed) {
		if (gui.icc.proofing_transform) {
			cmsDeleteTransform(gui.icc.proofing_transform);
			gui.icc.proofing_transform = NULL;
		}
	}
	gui.icc.same_primaries = FALSE;
	gui.icc.profile_changed = TRUE;
	g_mutex_unlock(&display_transform_mutex);
}

void validate_custom_profiles() {
	if (com.pref.icc.icc_path_monitor && com.pref.icc.icc_path_monitor[0] != '\0' && com.pref.icc.custom_monitor_profile_active) {
		g_mutex_lock(&monitor_profile_mutex);
		if (validate_profile(com.pref.icc.icc_path_monitor)) {
			if (gui.icc.monitor)
				cmsCloseProfile(gui.icc.monitor);
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc.icc_path_monitor, "r");
			if (!gui.icc.monitor) {
				gui.icc.monitor = com.pref.icc.rendering_intent == INTENT_PERCEPTUAL ? srgb_monitor_perceptual() : srgb_trc();
				siril_log_color_message(_("Error opening custom monitor profile. "
								"Monitor profile set to sRGB.\n"), "red");
			}
		} else {
			if (gui.icc.monitor)
				cmsCloseProfile(gui.icc.monitor);
			gui.icc.monitor = srgb_trc();
			siril_log_message(_("Warning: custom monitor profile set but could not "
								"be loaded. Display will use a sRGB profile with "
								"the standard sRGB TRC.\n"));
		}
		g_mutex_unlock(&monitor_profile_mutex);
	} else {
		if (gui.icc.monitor)
			cmsCloseProfile(gui.icc.monitor);
		gui.icc.monitor = com.pref.icc.rendering_intent == INTENT_PERCEPTUAL ? srgb_monitor_perceptual() : srgb_trc();
	}

	if (com.pref.icc.icc_path_soft_proof && com.pref.icc.icc_path_soft_proof[0] != '\0') {
		g_mutex_lock(&soft_proof_profile_mutex);
		if (validate_profile(com.pref.icc.icc_path_soft_proof)) {
			if (gui.icc.soft_proof)
				cmsCloseProfile(gui.icc.soft_proof);
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc.icc_path_soft_proof, "r");
		} else {
			if (gui.icc.soft_proof)
				cmsCloseProfile(gui.icc.soft_proof);
			gui.icc.soft_proof = NULL;
			siril_log_message(_("Warning: soft proofing profile set but could not "
								"be loaded. Soft proofing will be unavailable.\n"));
		}
		g_mutex_unlock(&soft_proof_profile_mutex);
	}

	g_mutex_lock(&default_profiles_mutex);
	if (com.pref.icc.working_gamut == TYPE_SRGB) {
		if (com.icc.working_standard)
			cmsCloseProfile(com.icc.working_standard);
		com.icc.working_standard = srgb_trc();
		if (com.icc.working_out)
			cmsCloseProfile(com.icc.working_out);
		com.icc.working_out = srgb_trcv2();
		if (com.icc.mono_standard)
			cmsCloseProfile(com.icc.mono_standard);
		com.icc.mono_standard = gray_srgbtrc();
		if (com.icc.mono_out)
			cmsCloseProfile(com.icc.mono_out);
		com.icc.mono_out = gray_srgbtrcv2();
	} else if (com.pref.icc.working_gamut == TYPE_REC2020) {
		if (com.icc.working_standard)
			cmsCloseProfile(com.icc.working_standard);
		com.icc.working_standard = rec2020_trc();
		if (com.icc.working_out)
			cmsCloseProfile(com.icc.working_out);
		com.icc.working_out = rec2020_trcv2();
		if (com.icc.mono_standard)
			cmsCloseProfile(com.icc.mono_standard);
		com.icc.mono_standard = gray_rec709trc();
		if (com.icc.mono_out)
			cmsCloseProfile(com.icc.mono_out);
		com.icc.mono_out = gray_rec709trcv2();
	} else {
		if (validate_profile(com.pref.icc.custom_icc_trc)) {
			if (com.icc.working_standard)
				cmsCloseProfile(com.icc.working_standard);
			com.icc.working_standard = cmsOpenProfileFromFile(com.pref.icc.custom_icc_trc, "r");
			if (!com.icc.working_standard) {
				com.icc.working_standard = srgb_trc();
				siril_log_color_message(_("Error opening nonlinear working profile. Profile set to sRGB.\n"), "red");
			}
		} else {
			com.icc.working_standard = srgb_trc();
		}
		if (com.icc.working_out)
			cmsCloseProfile(com.icc.working_out);
		com.icc.working_out = copyICCProfile(com.icc.working_standard);

		if (validate_profile(com.pref.icc.custom_icc_gray)) {
			if (com.icc.mono_standard)
				cmsCloseProfile(com.icc.mono_standard);
			com.icc.mono_standard = cmsOpenProfileFromFile(com.pref.icc.custom_icc_gray, "r");
			if (!com.icc.mono_standard) {
				com.icc.mono_standard = gray_srgbtrc();
				siril_log_color_message(_("Error opening matched grayscale working profile. Profile set to Gray with srGB tone response curve.\n"), "red");
			}
		} else {
			com.icc.mono_standard = gray_srgbtrc();
		}
		com.icc.mono_out = copyICCProfile(com.icc.mono_standard);
	}
	g_mutex_unlock(&default_profiles_mutex);
}

void initialize_profiles_and_transforms() {
	// Enable the fast float plugin (as long as the OS / lcms2 version blacklist isn't triggered)
#ifndef EXCLUDE_FF
	com.icc.context_single = cmsCreateContext(cmsFastFloatExtensions(), NULL);
	com.icc.context_threaded = cmsCreateContext(cmsFastFloatExtensions(), NULL);
	cmsPluginTHR(com.icc.context_threaded, cmsThreadedExtensions(CMS_THREADED_GUESS_MAX_THREADS, 0));
#else
	com.icc.context_single = cmsCreateContext(NULL, NULL);
	com.icc.context_threaded = cmsCreateContext(cmsThreadedExtensions(CMS_THREADED_GUESS_MAX_THREADS, 0));
#endif

	// Initialize FITS sRGB hint
	com.icc.srgb_hint = FALSE;

	// Set alarm codes for soft proof out-of-gamut warning
	cmsUInt16Number alarmcodes[16] = { 65535, 0, 65535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	cmsSetAlarmCodesTHR(com.icc.context_single, alarmcodes); // Out of gamut colours will be shown in magenta

	com.icc.rendering_flags |= ((com.pref.icc.rendering_bpc * cmsFLAGS_BLACKPOINTCOMPENSATION) & !(com.pref.icc.rendering_intent == INTENT_ABSOLUTE_COLORIMETRIC));
	gui.icc.proofing_flags = ((com.pref.icc.rendering_bpc * cmsFLAGS_BLACKPOINTCOMPENSATION) & !(com.pref.icc.rendering_intent == INTENT_ABSOLUTE_COLORIMETRIC)) | cmsFLAGS_SOFTPROOFING;

	// Working profiles
	com.icc.mono_linear = gray_linear();
	com.icc.srgb_profile = srgb_trc();
	validate_custom_profiles();

	// Target profiles for embedding in saved files
	com.icc.srgb_out = srgb_trcv2();
	com.icc.mono_out = gray_srgbtrcv2();

	// ICC availability
	gboolean available = (com.icc.mono_linear && com.icc.working_standard && com.icc.mono_standard && com.icc.working_out && com.icc.mono_out);
	gboolean gui_available = available && gui.icc.monitor;
	if ((com.headless && !available) || (!com.headless && !gui_available)) {
		siril_log_message(_("Error: standard color management profiles "
							"could not be loaded. Cannot continue. "
							"Please report this error.\n"));
		exit(1);
	} else {
		siril_log_message(_("Color management active.\n"));
#ifdef HAVE_LIBJPEG
#if LIBJPEG_TURBO_VERSION_NUMBER >= 2000000
		siril_log_message(_("JPEG ICC profiles supported.\n"));
#else
		siril_log_message(_("JPEG ICC profiles unsupported: libjpeg-turbo 2.0.0 or higher is required.\n"));
#endif
#endif
	}
	if (gui_available) {
		gui.icc.same_primaries = FALSE;
		gui.icc.profile_changed = TRUE;
	}
}

cmsUInt32Number get_planar_formatter_type(cmsColorSpaceSignature tgt, data_type t, gboolean force_16) {
	if (force_16)
		t = DATA_USHORT;
	switch (tgt) {
		case cmsSigGrayData:
			return (t == DATA_FLOAT ? TYPE_GRAY_FLT : TYPE_GRAY_16);
		case cmsSigRgbData:
			return (t == DATA_FLOAT ? TYPE_RGB_FLT_PLANAR : TYPE_RGB_16_PLANAR);
		case cmsSigXYZData:
			return (t == DATA_FLOAT ? TYPE_XYZ_FLT_PLANAR : TYPE_XYZ_16_PLANAR);
		case cmsSigLabData:
			return (t == DATA_FLOAT ? TYPE_Lab_FLT_PLANAR : TYPE_Lab_16_PLANAR);
		case cmsSigLuvData:
			return (t == DATA_FLOAT ? TYPE_Luv_FLT_PLANAR : TYPE_Luv_16_PLANAR);
		case cmsSigYCbCrData:
			return (t == DATA_FLOAT ? TYPE_YCbCr_FLT_PLANAR : TYPE_YCbCr_16_PLANAR);
		case cmsSigYxyData:
			return (t == DATA_FLOAT ? TYPE_Yxy_FLT_PLANAR : TYPE_Yxy_16_PLANAR);
		case cmsSigHsvData:
			return (t == DATA_FLOAT ? TYPE_HSV_FLT_PLANAR : TYPE_HSV_16_PLANAR);
		case cmsSigHlsData:
			return (t == DATA_FLOAT ? TYPE_HLS_FLT_PLANAR : TYPE_HLS_16_PLANAR);
		case cmsSigCmyData:
			return (t == DATA_FLOAT ? TYPE_CMY_FLT_PLANAR : TYPE_CMY_16_PLANAR);
		default:
			return 0;
	}
}

/* Even for mono images with a Gray profile, the display datatype is always RGB;
 * lcms2 takes care of populating all 3 output channels for us.
 */
cmsHTRANSFORM initialize_display_transform() {
	g_assert(gui.icc.monitor);
	cmsHTRANSFORM transform = NULL;
	if (gfit.icc_profile == NULL || !gfit.color_managed) {
		siril_debug_print("NULL display transform\n");
		return NULL;
	}
	cmsUInt32Number gfit_signature = cmsGetColorSpace(gfit.icc_profile);
	cmsUInt32Number srctype = get_planar_formatter_type(gfit_signature, gfit.type, TRUE);
	g_mutex_lock(&monitor_profile_mutex);
	// The display transform is always single threaded as OpenMP is used within the remap function
	transform = cmsCreateTransformTHR(com.icc.context_single, gfit.icc_profile, srctype, gui.icc.monitor, TYPE_RGB_16_PLANAR, com.pref.icc.rendering_intent, com.icc.rendering_flags);
	g_mutex_unlock(&monitor_profile_mutex);
	if (transform == NULL)
		siril_log_message("Error: failed to create display_transform!\n");
	else
		siril_debug_print("Display transform created (gfit->icc_profile to gui.icc.monitor)\n");
	return transform;
}

cmsHTRANSFORM initialize_export8_transform(fits* fit, gboolean threaded) {
	g_assert(com.icc.working_standard);
	cmsHTRANSFORM transform = NULL;
	if (fit->icc_profile == NULL || !fit->color_managed)
		return NULL;
	cmsUInt32Number fit_signature = cmsGetColorSpace(fit->icc_profile);
	cmsUInt32Number srctype = get_planar_formatter_type(fit_signature, fit->type, TRUE);
	cmsUInt32Number desttype = (fit->naxes[2] == 1 ? TYPE_GRAY_16 : TYPE_RGB_16_PLANAR);
	transform = sirilCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, srctype, (fit->naxes[2] == 3 ? com.icc.working_standard : com.icc.mono_standard), desttype, com.pref.icc.rendering_intent, com.icc.rendering_flags);
	if (transform == NULL)
		siril_log_message("Error: failed to create export colorspace transform!\n");
	return transform;
}

cmsHTRANSFORM initialize_proofing_transform() {
	g_assert(gui.icc.monitor);
	if (gfit.icc_profile == NULL || gfit.color_managed == FALSE)
		return NULL;
	cmsUInt32Number flags = gui.icc.proofing_flags;
	if (fit_icc_is_linear(&gfit))
		flags |= cmsFLAGS_NOOPTIMIZE;
	gboolean gamutcheck = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkgamut")));
	if (gamutcheck) {
		flags |= cmsFLAGS_GAMUTCHECK;
	}
	cmsUInt32Number type = (gfit.naxes[2] == 1 ? TYPE_GRAY_8 : TYPE_RGB_8_PLANAR);
	g_mutex_lock(&soft_proof_profile_mutex);
	g_mutex_lock(&monitor_profile_mutex);
	cmsHPROFILE proofing_transform = cmsCreateProofingTransformTHR(
						com.icc.context_single,
						gfit.icc_profile,
						type,
						gui.icc.monitor,
						TYPE_RGB_8_PLANAR,
						(gui.icc.soft_proof && com.pref.icc.soft_proofing_profile_active) ? gui.icc.soft_proof : gui.icc.monitor,
						com.pref.icc.rendering_intent,
						com.pref.icc.proofing_intent,
						flags);
	g_mutex_unlock(&monitor_profile_mutex);
	g_mutex_unlock(&soft_proof_profile_mutex);
	return proofing_transform;
}

/* Refreshes the display and proofing transforms after a profile is changed. */
void refresh_icc_transforms() {
	if (!com.headless) {
		gui.icc.same_primaries = same_primaries(gfit.icc_profile, gui.icc.monitor, (gui.icc.soft_proof && com.pref.icc.soft_proofing_profile_active) ? gui.icc.soft_proof : NULL);
		g_mutex_lock(&display_transform_mutex);
		if (gui.icc.proofing_transform)
			cmsDeleteTransform(gui.icc.proofing_transform);
		gui.icc.proofing_transform = initialize_proofing_transform();
		g_mutex_unlock(&display_transform_mutex);
		gui.icc.profile_changed = TRUE;

	}
	if (is_preview_active())
		copy_gfit_icc_to_backup();
	check_gfit_profile_identical_to_monitor();
}

// Returns the full ICC profile data
unsigned char* get_icc_profile_data(cmsHPROFILE profile, guint32 *len) {
	unsigned char* block = NULL;
	cmsUInt32Number length;
	cmsBool ret = cmsSaveProfileToMem(profile, NULL, &length);
	if (length > 0) {
		block = malloc(length * sizeof(BYTE));
		ret = cmsSaveProfileToMem(profile, (void*) block, &length);
	}
	if (!ret) {
		siril_debug_print("Error preparing ICC profile for embedding...\n");
		return NULL;
	}
	*len = length;
	return block;
}

/* Intended for use if a fits has no profile, to decide what to assign.
 * This is not definitive, but it checks the FITS header HISTORY for signs of
 * GHT, Asinh or Histogram stretches (or autostretches) having been applied.
 */
static gboolean fit_appears_stretched(fits* fit) {
	GSList* entry = NULL;
	if (fit->history) {
		entry = fit->history;
		while (entry) {
			if (strstr(entry->data, "Histogram Transf."))
				return TRUE;
			if (strstr(entry->data, "Asinh"))
				return TRUE;
			if (strstr(entry->data, "Midtones"))
				return TRUE;
			if (strstr(entry->data, "Autostretch"))
				return TRUE;
			// this catches AutoGHS too
			if (strstr(entry->data, "GHS") && !strstr(entry->data, "LINEAR BP"))
				return TRUE;
			entry = entry->next;
		}
	}
	return FALSE;
}

/* Adapted from GIMP code */
cmsBool fit_icc_is_linear(fits *fit) {
	if (!(fit->color_managed) || fit->icc_profile == NULL)
		return FALSE;
	cmsToneCurve *curve;
	if (! cmsIsMatrixShaper (fit->icc_profile))
		return FALSE;

	if (cmsIsCLUT (fit->icc_profile, INTENT_PERCEPTUAL, LCMS_USED_AS_INPUT))
		return FALSE;

	if (cmsIsCLUT (fit->icc_profile, INTENT_PERCEPTUAL, LCMS_USED_AS_OUTPUT))
		return FALSE;

	if (siril_color_profile_is_rgb (fit->icc_profile))
		{
		curve = cmsReadTag(fit->icc_profile, cmsSigRedTRCTag);
		if (curve == NULL || ! cmsIsToneCurveLinear (curve))
			return FALSE;

		curve = cmsReadTag (fit->icc_profile, cmsSigGreenTRCTag);
		if (curve == NULL || ! cmsIsToneCurveLinear (curve))
			return FALSE;

		curve = cmsReadTag (fit->icc_profile, cmsSigBlueTRCTag);
		if (curve == NULL || ! cmsIsToneCurveLinear (curve))
			return FALSE;
		}
	else if (siril_color_profile_is_gray (fit->icc_profile))
		{
		curve = cmsReadTag(fit->icc_profile, cmsSigGrayTRCTag);
		if (curve == NULL || ! cmsIsToneCurveLinear (curve))
			return FALSE;
		}
	else
		{
		return FALSE;
		}

	return TRUE;
}

static gboolean siril_color_profile_get_rgb_matrix_colorants (cmsHPROFILE *profile, cmsCIEXYZTRIPLE *XYZtriple, cmsCIEXYZ *whitepoint) {
	cmsFloat64Number prior_adaptation_state = cmsSetAdaptationStateTHR(com.icc.context_single, 0);
	double redrgb[3] = { 1.0, 0.0, 0.0 };
	double greenrgb[3] = { 0.0, 1.0, 0.0 };
	double bluergb[3] = { 0.0, 0.0, 1.0 };
	double whitergb[3] = { 1.0, 1.0, 1.0 };
	cmsHPROFILE profileXYZ = cmsCreateXYZProfile();
	cmsHTRANSFORM transform = cmsCreateTransformTHR(com.icc.context_single, profile, TYPE_RGB_DBL, profileXYZ, TYPE_XYZ_DBL, INTENT_ABSOLUTE_COLORIMETRIC, cmsFLAGS_NOCACHE);
	cmsCloseProfile(profileXYZ);
	if (!transform) {
		cmsSetAdaptationStateTHR(com.icc.context_single, prior_adaptation_state);
		return FALSE;
	}
	cmsDoTransform(transform, &redrgb, &XYZtriple->Red, 1);
	cmsDoTransform(transform, &greenrgb, &XYZtriple->Green, 1);
	cmsDoTransform(transform, &bluergb, &XYZtriple->Blue, 1);
	cmsDoTransform(transform, &whitergb, whitepoint, 1);
	cmsDeleteTransform(transform);
	cmsSetAdaptationStateTHR(com.icc.context_single, prior_adaptation_state);
	return TRUE;
}

static void
	siril_color_profile_set_tag (cmsHPROFILE profile,
								cmsTagSignature sig,
								const gchar *tag) {
	cmsMLU *mlu;

	mlu = cmsMLUalloc (NULL, 1);
	cmsMLUsetASCII (mlu, "en", "US", tag);
	cmsWriteTag (profile, sig, mlu);
	cmsMLUfree (mlu);
}

static void
	siril_color_profile_make_tag (cmsHPROFILE profile,
								cmsTagSignature sig,
								const gchar *siril_tag,
								const gchar *siril_prefix,
								const gchar *siril_prefix_alt,
								const gchar *original_tag) {
	if (! original_tag || ! strlen (original_tag) ||
		! strcmp (original_tag, siril_tag)) {
		/* if there is no original tag (or it is the same as the new
		* tag), just use the new tag
		*/

		siril_color_profile_set_tag (profile, sig, siril_tag);
	} else {
		/* otherwise prefix the existing tag with a gimp prefix
		* indicating that the profile has been generated
		*/

		if (g_str_has_prefix (original_tag, siril_prefix)) {
			/* don't add multiple GIMP prefixes */
			siril_color_profile_set_tag (profile, sig, original_tag);
		}else if (siril_prefix_alt &&
				g_str_has_prefix (original_tag, siril_prefix_alt)) {
			/* replace GIMP prefix_alt by prefix */
			gchar *new_tag = g_strconcat (siril_prefix,
											original_tag + strlen (siril_prefix_alt),
											NULL);

			siril_color_profile_set_tag (profile, sig, new_tag);
			g_free (new_tag);
		} else {
			gchar *new_tag = g_strconcat (siril_prefix,
											original_tag,
											NULL);

			siril_color_profile_set_tag (profile, sig, new_tag);
			g_free (new_tag);
		}
	}
}

cmsHPROFILE siril_color_profile_linear_from_color_profile (cmsHPROFILE profile) {
	cmsHPROFILE target_profile;
	cmsCIEXYZTRIPLE XYZtriple;
	cmsCIEXYZ whitepoint;
	cmsToneCurve *curve;

	if (siril_color_profile_is_rgb (profile)) {
		if (! siril_color_profile_get_rgb_matrix_colorants (profile, &XYZtriple, &whitepoint))
			return NULL;
	} else if (! siril_color_profile_is_gray (profile)) {
		return NULL;
	}

	target_profile = cmsCreateProfilePlaceholder (0);

	cmsSetProfileVersion (target_profile, 4.3);
	cmsSetDeviceClass (target_profile, cmsSigDisplayClass);
	cmsSetPCS (target_profile, cmsSigXYZData);

	cmsWriteTag (target_profile, cmsSigMediaWhitePointTag, &whitepoint);

	curve = cmsBuildGamma (NULL, 1.00);

	siril_color_profile_make_tag (target_profile, cmsSigProfileDescriptionTag,
									"linear TRC from unnamed profile",
									"linear TRC from ",
									"sRGB TRC from ",
									siril_color_profile_get_description (profile));

	if (siril_color_profile_is_rgb (profile)) {

		cmsSetColorSpace (target_profile, cmsSigRgbData);

		cmsWriteTag (target_profile, cmsSigRedColorantTag,   &XYZtriple.Red);
		cmsWriteTag (target_profile, cmsSigGreenColorantTag, &XYZtriple.Green);
		cmsWriteTag (target_profile, cmsSigBlueColorantTag,  &XYZtriple.Blue);

		cmsWriteTag (target_profile, cmsSigRedTRCTag,   curve);
		cmsWriteTag (target_profile, cmsSigGreenTRCTag, curve);
		cmsWriteTag (target_profile, cmsSigBlueTRCTag,  curve);
	} else {
		cmsSetColorSpace (target_profile, cmsSigGrayData);

		cmsWriteTag (target_profile, cmsSigGrayTRCTag, curve);
	}

	cmsFreeToneCurve (curve);

	siril_color_profile_make_tag (target_profile, cmsSigDeviceMfgDescTag,
								"Siril",
								"Siril from ", NULL,
								siril_color_profile_get_manufacturer (profile));
	siril_color_profile_make_tag (target_profile, cmsSigDeviceModelDescTag,
								"Generated by GIMP",
								"Siril from ", NULL,
								siril_color_profile_get_model (profile));
	siril_color_profile_make_tag (target_profile, cmsSigCopyrightTag,
								"Public Domain",
								"Siril from ", NULL,
								siril_color_profile_get_copyright (profile));

	return target_profile;
}

/* Provides a sanity check of the ICC profile attached to a fits. This is used
 * because we shouldn't trust that an embedded profile in an imported file is
 * sensible. FIrst it checks that the channel count matches. Also, if there is
 * no profile attached to the fits it also looks for evidence that the image has
 * been stretched in the past, to decide whether to assign a linear or non-
 * linear profile.
 */
void check_profile_correct(fits* fit) {
	if (!fit->icc_profile) {
		if (com.icc.srgb_hint) {
			siril_log_message(_("FITS did not contain an ICC profile but is declared to be stretched. Assigning a sRGB color profile.\n"));
			// sRGB because this is the implicit assumption made in older versions
			fit->icc_profile = fit->naxes[2] == 1 ? gray_srgbtrc() : srgb_trc();
			// Clear the hint
			com.icc.srgb_hint = FALSE;
		} else if (fit_appears_stretched(fit)) {
			siril_log_message(_("FITS did not contain an ICC profile. It appears to have been stretched using an older version of Siril. Assigning a sRGB color profile.\n"));
			// sRGB because this is the implicit assumption made in older versions
			fit->icc_profile = fit->naxes[2] == 1 ? gray_srgbtrc() : srgb_trc();
			color_manage(fit, TRUE);
		} else {
			siril_debug_print("FITS did not contain an ICC profile and no hints were available in the HISTORY header.\n");
			fit->icc_profile = NULL;
			color_manage(fit, FALSE);
		}
	} else {
		cmsColorSpaceSignature sig = cmsGetColorSpace(fit->icc_profile);
		cmsUInt32Number chans = cmsChannelsOf(sig);
		if (chans != fit->naxes[2]) {
			cmsCloseProfile(fit->icc_profile);
			fit->icc_profile = NULL;
			color_manage(fit, FALSE);
			siril_log_color_message(_("Warning: embedded ICC profile channel count does not match image channel count. Color management is disabled for this image. To re-enable it, an ICC profile must be assigned using the Color Management menu item.\n"), "salmon");
		}
	}
	if (fit->color_managed && !fit->icc_profile) {
		color_manage(fit, FALSE);
		siril_debug_print("fit->color_managed inconsistent with missing profile");
	}
}

/* This function returns a separate copy of the cmsHPROFILE provided as the
 * argument. It is used so that a reference profile can be copied to many
 * fits images without worrying about freeing the original when the fits is
 * cleared.
 */
cmsHPROFILE copyICCProfile(cmsHPROFILE profile) {
	cmsUInt32Number length = 0;
	cmsUInt8Number* block = NULL;
	cmsBool ret = FALSE;
	cmsHPROFILE retval = NULL;
	if (!profile) {
		return NULL;
	} else {
		ret = cmsSaveProfileToMem(profile, NULL, &length);
	}
	if (length > 0) {
		block = malloc(length * sizeof(BYTE));
		ret = cmsSaveProfileToMem(profile, (void*) block, &length);
	}
	if (ret) {
		retval = cmsOpenProfileFromMem(block, length);
	}
	free(block);
	return retval;
}

/* This function is for initializing the profile during file
 * import from non-FITS formats. It assumes that if the image does not contain
 * a profile then it should be considered sRGB. This is an industry-standard
 * assumption.
 */
void fits_initialize_icc(fits *fit, cmsUInt8Number* EmbedBuffer, cmsUInt32Number EmbedLen) {
	if (EmbedBuffer) {
		// If there is an embedded profile we will use it
		fit->icc_profile = cmsOpenProfileFromMem(EmbedBuffer, EmbedLen);
		check_profile_correct(fit);
	} else {
		// If there is no embedded profile we assume the usual sRGB TRC
		fit->icc_profile = copyICCProfile((fit->naxes[2] == 1) ? com.icc.mono_standard : com.icc.srgb_profile);
	}
	color_manage(fit, TRUE);
}

cmsUInt8Number *siril_icc_profile_to_buffer(cmsHPROFILE profile, cmsUInt32Number *length) {
	if (!profile) {
		length = 0;
		return NULL;
	}
	cmsBool ret = cmsSaveProfileToMem(profile, NULL, length);
	if (!ret || length == 0)
		return NULL;
	void *buffer = malloc (*length * sizeof(cmsUInt8Number));
	if (!buffer) {
		PRINT_ALLOC_ERR;
		return NULL;
	} else {
		ret = cmsSaveProfileToMem(profile, buffer, length);
		if (!ret || !buffer)
			return NULL;
	}
	return buffer;
}

/* Compares two profiles. Returns TRUE if the profiles are identical or
 * FALSE if they are not. Note, this is a conservative check and requires
 * the profiles to be exactly identical in all respects. Even something
 * that makes no difference to the colorspace described by the profile, such
 * as the description tag, can trigger a FALSE return. This is not intended
 * as a rigorous check of the colorspaces described by two profiles, only as
 * an opportunistic means of avoiding unnecessary transforms if the profiles
 * are guaranteed to be the same.
 */
cmsBool profiles_identical(cmsHPROFILE a, cmsHPROFILE b) {
	cmsUInt8Number *block_a = NULL, *block_b = NULL;
	cmsUInt32Number length_a = 0, length_b = 0;
	cmsBool retval = FALSE;
	const gsize header_len = sizeof (cmsICCHeader);

	if ((!a && !b) || (!a || !b))
		goto ERROR_OR_FINISH;

	if (a) {
		block_a = siril_icc_profile_to_buffer(a, &length_a);
	}
	if (b) {
		block_b = siril_icc_profile_to_buffer(b, &length_b);
	}
	// If a profile can't be saved to a buffer or the lengths don't match
	// we can already return FALSE
	if (length_a != length_b)
		goto ERROR_OR_FINISH;

	retval = (memcmp(block_a + header_len, block_b + header_len, length_a - header_len) == 0) ? TRUE : FALSE;

ERROR_OR_FINISH:

	free(block_a);
	free(block_b);
	return retval;
}

static void set_source_information() {
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_current_profile_label");
	GtkLabel* mfr_label = (GtkLabel*) lookup_widget("icc_mfr_label");
	GtkLabel* copyright_label = (GtkLabel*) lookup_widget("icc_copyright_label");
	if (!gfit.color_managed) {
		siril_debug_print("Target image is not color managed\n");
		gtk_label_set_text(label, _("No ICC profile"));
		gtk_label_set_text(mfr_label, "");
		gtk_label_set_text(copyright_label, "");
		return;
	}
	if (!gfit.icc_profile) {
		siril_debug_print("Target profile is NULL\n");
		gtk_label_set_text(label, _("No ICC profile"));
		gtk_label_set_text(mfr_label, "");
		gtk_label_set_text(copyright_label, "");
		return;
	}
	// Set description
	gchar *buffer = siril_color_profile_get_description(gfit.icc_profile);
	if (buffer)
		gtk_label_set_text(label, buffer);
	free(buffer);

	// Set manufacturer
	buffer = siril_color_profile_get_manufacturer(gfit.icc_profile);
	if (buffer)
		gtk_label_set_text(mfr_label, buffer);
	free(buffer);

	// Set copyright
	buffer = siril_color_profile_get_copyright(gfit.icc_profile);
	if (buffer)
		gtk_label_set_text(copyright_label, buffer);
	free(buffer);
}

void set_icc_description_in_TIFF() {
	// Set description
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_save_label");
	gchar *buffer = NULL;
	if (gfit.icc_profile) {
		gtk_widget_set_tooltip_text((GtkWidget*) label, "");
		int length = cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoDescription, "en", "US", NULL, 0);
		if (length) {
			buffer = (char*) g_malloc(length * sizeof(char));
			cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoDescription, "en", "US", buffer, length);
		}
	} else {
			gtk_widget_set_tooltip_text((GtkWidget*) label, _("To write an ICC profile, assign a profile to the image using the Color Management dialog."));
			buffer = g_strdup(_("Image is not color managed. Will not write an ICC profile."));
	}
	gtk_label_set_text(label, buffer);
	g_free(buffer);
}

static void set_target_information() {
	if (!target) {
		return;
	}
	// Set description
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_target_profile_label");
	GtkLabel* mfr_label = (GtkLabel*) lookup_widget("icc_target_mfr_label");
	GtkLabel* copyright_label = (GtkLabel*) lookup_widget("icc_target_copyright_label");
	int length = cmsGetProfileInfoASCII(target, cmsInfoDescription, "en", "US", NULL, 0);
	char *buffer = NULL;
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(target, cmsInfoDescription, "en", "US", buffer, length);
		gtk_label_set_text(label, buffer);
		free(buffer);
	}

	// Set manufacturer
	length = cmsGetProfileInfoASCII(target, cmsInfoManufacturer, "en", "US", NULL, 0);
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(target, cmsInfoManufacturer, "en", "US", buffer, length);
		gtk_label_set_text(mfr_label, buffer);
		free(buffer);
	}

	// Set copyright
	length = cmsGetProfileInfoASCII(target, cmsInfoCopyright, "en", "US", NULL, 0);
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(target, cmsInfoCopyright, "en", "US", buffer, length);
		gtk_label_set_text(copyright_label, buffer);
		free(buffer);
	}
}

void siril_colorspace_transform(fits *fit, cmsHPROFILE profile) {

	// If profile is NULL, we remove the profile from fit to match it. This is an unusual
	// case but the behaviour is consistent.
	if (!profile) {
		if (fit->icc_profile) {
			cmsCloseProfile(fit->icc_profile);
			fit->history = g_slist_append(fit->history, g_strdup(_("ICC profile removed")));
		}
		fit->icc_profile = NULL;
		color_manage(fit, FALSE);
		return;
	}

	cmsUInt32Number target_colorspace = cmsGetColorSpace(profile);
	cmsUInt32Number target_colorspace_channels = cmsChannelsOf(target_colorspace);
	cmsUInt32Number fit_colorspace_channels;

	// If fit->color_managed is FALSE, we assign the profile rather than convert to it
	if (!fit->color_managed || !fit->icc_profile) {
		fit_colorspace_channels = fit->naxes[2];
		if (fit_colorspace_channels == target_colorspace_channels) {
			if (fit->icc_profile)
				cmsCloseProfile(fit->icc_profile);
			fit->icc_profile = copyICCProfile(profile);
			siril_debug_print("siril_colorspace_transform() assigned a profile\n");
			gchar *desc = siril_color_profile_get_description(profile);
			fit->history = g_slist_append(fit->history, g_strdup_printf(_("Assigned ICC profile: %s"), desc));
			g_free(desc);
			refresh_icc_transforms();
			color_manage(fit, TRUE);
			return;
		} else {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("Error"), _("Image number of channels does not match color profile number of channels. Cannot assign this profile to this image."));
			return;
		}
	}

	// fit is color managed, so we really have to do the transform
	cmsUInt32Number fit_colorspace = cmsGetColorSpace(fit->icc_profile);
	fit_colorspace_channels = cmsChannelsOf(fit_colorspace);

	if (target_colorspace != cmsSigGrayData && target_colorspace != cmsSigRgbData) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Siril only supports representing the image in Gray or RGB color spaces. You cannot assign or convert to non-RGB color profiles"));
		return;
	}
	void *data = NULL;
	cmsUInt32Number srctype, desttype;
	size_t npixels = fit->rx * fit->ry;
	// convert from fit->icc_profile to profile
	gboolean threaded = !get_thread_run();
	srctype = get_planar_formatter_type(fit_colorspace, fit->type, FALSE);
	desttype = get_planar_formatter_type(target_colorspace, fit->type, FALSE);
	cmsHTRANSFORM transform = cmsCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, srctype, profile, desttype, com.pref.icc.export_intent, com.icc.rendering_flags);
	if (transform) {
		if (fit_colorspace_channels < target_colorspace_channels)
			fits_change_depth(fit, target_colorspace_channels);
		data = (fit->type == DATA_FLOAT) ? (void *) fit->fdata : (void *) fit->data;
		cmsUInt32Number datasize = fit->type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
		cmsUInt32Number bytesperline = fit->rx * datasize;
		cmsUInt32Number bytesperplane = npixels * datasize;
		cmsDoTransformLineStride(transform, data, data, fit->rx, fit->ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
		cmsDeleteTransform(transform);
		cmsCloseProfile(fit->icc_profile);
		if (fit_colorspace_channels > target_colorspace_channels)
			fits_change_depth(fit, target_colorspace_channels);
		fit->icc_profile = copyICCProfile(profile);
		refresh_icc_transforms();
			gchar *desc = siril_color_profile_get_description(profile);
			fit->history = g_slist_append(fit->history, g_strdup_printf(_("Converted to ICC profile: %s"), desc));
			g_free(desc);
		color_manage(fit, TRUE);
		siril_debug_print("siril_colorspace_transform() converted a profile\n");
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Failed to create colorspace transform."));
	}
}

/* This function converts the current image to the working color space.
 * It honours the preferences com.pref.icc.autoassignment and
 * com.pref.icc.autoconversion.
 * It is used on loading a file and before stretching a file.
 */
void icc_auto_assign_or_convert(fits *fit, icc_assign_type occasion) {
	// Scripts are responsible for managing this themselves
	if (com.script)
		return;

	gboolean proceed = FALSE;

	// Handle images that have been extracted as channels from a 3-color image
	// It doesn't make sense to assign a color profile to these, they should
	// be treated as raw data until the user decides how to handle them.
	if (fit->naxes[2] == 1 && fit->history) {
		GSList *last_history = g_slist_last(fit->history);

		if (g_strrstr(last_history->data, "Extraction")) {
			siril_log_message(_("This image is a channel extracted from a 3-channel "
								"image. Cannot sensibly auto-assign an ICC profile. "
								"You must manually assign one at an appropriate time.\n"));
			return;
		}
	}

	// If there is no existing profile and the appropriate preference is set,
	// assign the working color profile
	if (!fit->color_managed || !fit->icc_profile) {
		if (com.pref.icc.autoassignment & occasion) {
			if (fit_appears_stretched(fit)) {
				// if the fit has previously been stretched, we assign the sRGB TRC:
				// this will then be converted to the working color space rather than
				// just assigned which would cause a color shift
				fit->icc_profile = fit->naxes[2] == 1 ? gray_srgbtrc() : srgb_trc();
				// color_manage() is called later from siril_colorspace_transform()
			} else if (com.pref.icc.pedantic_linear && !(occasion & ICC_ASSIGN_ON_STRETCH)) {
				fit->icc_profile = siril_color_profile_linear_from_color_profile (fit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.working_standard);
				// Unless we're about to stack, leave the image with a linear profile and return
				gchar *desc = siril_color_profile_get_description(fit->icc_profile);
				fit->history = g_slist_append(fit->history, g_strdup_printf(_("Assigned ICC profile: %s"), desc));
				g_free(desc);
				if (!(occasion & ICC_ASSIGN_ON_STACK)) {
					color_manage(fit, TRUE);
					return;
				}
			}
			proceed = TRUE;
		}
	} else {

		// If the preference is never to autoconvert, we have nothing to do
		if (com.pref.icc.autoconversion == ICC_NEVER_AUTOCONVERT)
			return;

		// If the image is already in the working color space, we have nothing to do
		if (fit->color_managed && profiles_identical(fit->icc_profile, fit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.working_standard))
			return;

		// If the preference is that we always convert, trigger the conversion
		if (com.pref.icc.autoconversion == ICC_ALWAYS_AUTOCONVERT)
			proceed = TRUE;

		else if (com.pref.icc.autoconversion == ICC_ASK_TO_CONVERT) {
			if (!fit->color_managed) {
				proceed = siril_confirm_dialog(_("Color Management"), _("The current image is not color managed. Do you want to assign the working color space?"), _("Assign"));
			} else {
				proceed = siril_confirm_dialog(_("Color Management"), _("Do you want to convert this image to the working color space?"), _("Convert"));
			}
		}
	}

	// If the conversion has been triggered automatically or via the dialog, do it.
	if (proceed) {
		set_cursor_waiting(TRUE);
		// siril_colorspace_transform takes care of hitherto non-color managed images, and assigns a profile instead of converting them
		siril_colorspace_transform(fit, (fit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.working_standard));
		if (fit == &gfit) {
			set_source_information();
			refresh_icc_transforms();
			notify_gfit_modified();
		}
		set_cursor_waiting(FALSE);
	}
}

/* This function automatically assigns the working color space or switches off
 * color management of the image, depending on the user preference. The
 * icc_assign_type parameter defines the occasion for which the function is
 * being called, it is checked against the preference com.pref.icc.autoassignment.
 * It is used when auto-assigning a profile to a new image created by stacking or
 * RGB composition or pixelmath composition. It takes no action if called from
 * a script (the icc_assign command must be used).
 */

void icc_auto_assign(fits *fit, icc_assign_type occasion) {
	// Check if the occasion matches the preference
	if (com.pref.icc.autoassignment & occasion) {
		siril_debug_print("Auto assigning working profile\n");
		set_cursor_waiting(TRUE);
		// siril_colorspace_transform takes care of hitherto non-color managed images, and assigns a profile instead of converting them
		fit->icc_profile = copyICCProfile((fit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.working_standard));
		color_manage(fit, TRUE);
	} else {
		if (fit->icc_profile)
			cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = NULL;
		color_manage(fit, FALSE);
	}
	if (fit == &gfit) {
		set_source_information();
		refresh_icc_transforms();
		notify_gfit_modified();
	}
	set_cursor_waiting(FALSE);
}

const char* default_system_icc_path() {
#ifdef _WIN32
	return "C:\\Windows\\System32\\spool\\drivers\\color";
#endif
#ifdef _MACOS
	return "/Library/ColorSync/Profiles";
#endif
	return "/usr/share/color/icc";
}

/* This function wraps cmsCreateTransform within a gMutex.
 * It should always be used when one of the profiles in com.icc
 * is used to create the transform, to prevent bad things happening
 * if the user changes the profile at the same time as Siril is trying
 * to create the transform. It is not required for transforms solely using
 * gui.icc.display or gui.icc.soft_proof, as they have their own
 * gMutices. It is also not required for creating transforms between
 * fit->icc_profile and another profile (such as the temporary one created
 * when carrying out autostretches), as there is no risk of a concurrency
 * clash.
 */

cmsHTRANSFORM sirilCreateTransformTHR(cmsContext Context, cmsHPROFILE Input, cmsUInt32Number InputFormat, cmsHPROFILE Output, cmsUInt32Number OutputFormat, cmsUInt32Number Intent, cmsUInt32Number dwFlags) {
	cmsHTRANSFORM transform = NULL;
	g_mutex_lock(&default_profiles_mutex);
	transform = cmsCreateTransformTHR(Context, Input, InputFormat, Output, OutputFormat, Intent, dwFlags);
	g_mutex_unlock(&default_profiles_mutex);
	return transform;
}

static void reset_working_profile_to_srgb() {
	if (com.icc.working_standard)
		cmsCloseProfile(com.icc.working_standard);
	com.icc.working_standard = srgb_trc();
	if (com.icc.working_out)
		cmsCloseProfile(com.icc.working_out);
	com.icc.working_out = srgb_trcv2();
	if (com.icc.mono_standard)
		cmsCloseProfile(com.icc.mono_standard);
	com.icc.mono_standard = gray_srgbtrc();
	if (com.icc.mono_out)
		cmsCloseProfile(com.icc.mono_out);
	com.icc.mono_out = gray_srgbtrcv2();
	com.pref.icc.working_gamut = TYPE_SRGB;
}

static void error_loading_profile() {
	siril_message_dialog(GTK_MESSAGE_ERROR, _("Error loading profile"),
						 _("The selected profile could not be loaded or did not contain a valid ICC profile. Defaulting to sRGB."));
	reset_working_profile_to_srgb();
}

static void reset_custom_to_srgb() {
	siril_log_color_message(_("Error: the preferred colorspace profiles are not all set, or some are not valid. "
							  "You need to set both a RGB and a Gray profile. Defaulting to sRGB.\n"), "red");
	reset_working_profile_to_srgb();
}

///// Preferences callbacks

// Being able to alter the monitor and soft_proof profiles and intents from the GTK thread means all operations
// that use these profiles need to go inside mutexes to prevent the profiles being ripped out from under them.
// This is not required for profiles in FITS structures as these are not subject to thread contention issues in
// the same way.

void update_profiles_after_gamut_change() {
	siril_log_message(_("Updating working profiles.\n"));
	working_gamut_type working_gamut = com.pref.icc.working_gamut;
	g_mutex_lock(&default_profiles_mutex);
	switch (working_gamut) {
		case TYPE_SRGB:
			if (com.icc.working_standard)
				cmsCloseProfile(com.icc.working_standard);
			com.icc.working_standard = srgb_trc();
			if (com.icc.working_out)
				cmsCloseProfile(com.icc.working_out);
			com.icc.working_out = srgb_trcv2();
			if (com.icc.mono_standard)
				cmsCloseProfile(com.icc.mono_standard);
			com.icc.mono_standard = gray_srgbtrc();
			if (com.icc.mono_out)
				cmsCloseProfile(com.icc.mono_out);
			com.icc.mono_out = gray_srgbtrcv2();
			break;
		case TYPE_REC2020:
			if (com.icc.working_standard)
				cmsCloseProfile(com.icc.working_standard);
			com.icc.working_standard = rec2020_trc();
			if (com.icc.working_out)
				cmsCloseProfile(com.icc.working_out);
			com.icc.working_out = rec2020_trcv2();
			if (com.icc.mono_standard)
				cmsCloseProfile(com.icc.mono_standard);
			com.icc.mono_standard = gray_rec709trc();
			if (com.icc.mono_out)
				cmsCloseProfile(com.icc.mono_out);
			com.icc.mono_out = gray_rec709trcv2();
			break;
		case TYPE_CUSTOM:
			if (!(com.pref.icc.custom_icc_trc && com.pref.icc.custom_icc_gray)) {
				reset_custom_to_srgb();
				break;
			}
			// Custom profiles will also be used for the output profile
			if (com.icc.working_standard)
				cmsCloseProfile(com.icc.working_standard);
			if (!(com.pref.icc.custom_icc_trc && (com.icc.working_standard =
						cmsOpenProfileFromFile(com.pref.icc.custom_icc_trc, "r")))) {
				error_loading_profile();
				break;
			} else {
				// copy to working_out as we don't require a separate profile for embedding in a custom profile set
				if (com.icc.working_out)
					cmsCloseProfile(com.icc.working_out);
				com.icc.working_out = copyICCProfile(com.icc.working_standard);
			}

			// Custom profiles will also be used for the output profile
			if (com.icc.mono_standard)
				cmsCloseProfile(com.icc.mono_standard);
			if (!(com.pref.icc.custom_icc_gray && (com.icc.mono_standard =
						cmsOpenProfileFromFile(com.pref.icc.custom_icc_gray, "r")))) {
				error_loading_profile();
			} else {
				if (com.icc.mono_out)
					cmsCloseProfile(com.icc.mono_out);
				com.icc.mono_out = copyICCProfile(com.icc.mono_standard);
			}
	}
	g_mutex_unlock(&default_profiles_mutex);
	refresh_icc_transforms();
}

void initialize_icc_preferences_widgets() {
	GtkToggleButton *monitortogglebutton = (GtkToggleButton*) lookup_widget("custom_monitor_profile_active");
	GtkFileChooser *monitorfilechooser = (GtkFileChooser*) lookup_widget("pref_custom_monitor_profile");

	GtkToggleButton *proofingtogglebutton = (GtkToggleButton*) lookup_widget("custom_proofing_profile_active");
	GtkFileChooser *proofingfilechooser = (GtkFileChooser*) lookup_widget("pref_soft_proofing_profile");

	if (!gtk_file_chooser_get_filename(monitorfilechooser)) {
		gtk_toggle_button_set_active(monitortogglebutton, FALSE);
		gtk_widget_set_sensitive((GtkWidget*) monitortogglebutton, FALSE);
	} else {
		gtk_toggle_button_set_active(monitortogglebutton, TRUE);
		gtk_widget_set_sensitive((GtkWidget*) monitortogglebutton, TRUE);
	}

	gtk_toggle_button_set_active(proofingtogglebutton, (gui.icc.soft_proof != NULL));
	if (!gtk_file_chooser_get_filename(proofingfilechooser)) {
		gtk_widget_set_sensitive((GtkWidget*) proofingtogglebutton, FALSE);
	} else {
		gtk_widget_set_sensitive((GtkWidget*) proofingtogglebutton, TRUE);
	}
}

void on_pref_custom_monitor_profile_file_set(GtkFileChooser* filechooser, gpointer user_data) {
	GtkToggleButton *togglebutton = (GtkToggleButton*) lookup_widget("custom_monitor_profile_active");
	gchar *filename = gtk_file_chooser_get_filename(filechooser);
	if (filename) {
		gtk_widget_set_sensitive((GtkWidget*) togglebutton, TRUE);
	}
}

void on_pref_soft_proofing_profile_file_set(GtkFileChooser* filechooser, gpointer user_data) {
	GtkToggleButton *togglebutton = (GtkToggleButton*) lookup_widget("custom_proofing_profile_active");
	gchar *filename = gtk_file_chooser_get_filename(filechooser);
	if (filename) {
		gtk_widget_set_sensitive((GtkWidget*) togglebutton, TRUE);
	}
}

void on_monitor_profile_clear_clicked(GtkButton* button, gpointer user_data) {
	GtkFileChooser *filechooser = (GtkFileChooser*) lookup_widget("pref_custom_monitor_profile");
	GtkToggleButton *togglebutton = (GtkToggleButton*) lookup_widget("custom_monitor_profile_active");
	gtk_file_chooser_unselect_all(filechooser);
	cmsHPROFILE old_monitor = copyICCProfile(gui.icc.monitor);
	g_mutex_lock(&monitor_profile_mutex);
	if (com.pref.icc.icc_path_monitor && com.pref.icc.icc_path_monitor[0] != '\0') {
		g_free(com.pref.icc.icc_path_monitor);
		com.pref.icc.icc_path_monitor = NULL;
		cmsCloseProfile(gui.icc.monitor);
		gui.icc.monitor = com.pref.icc.rendering_intent == INTENT_PERCEPTUAL ? srgb_monitor_perceptual() : srgb_trc();
		if (gui.icc.monitor) {
			siril_log_message(_("Monitor ICC profile set to sRGB\n"));
		} else {
			siril_log_color_message(_("Fatal error: standard sRGB ICC profile could not be loaded.\n"), "red");
			exit(1);
		}
	}
	g_mutex_unlock(&monitor_profile_mutex);
	gtk_toggle_button_set_active(togglebutton, FALSE);
	gtk_widget_set_sensitive((GtkWidget*) togglebutton, FALSE);
	if (!profiles_identical(old_monitor, gui.icc.monitor)) {
		refresh_icc_transforms();
		redraw(REMAP_ALL);
	}
	cmsCloseProfile(old_monitor);
}

void on_proofing_profile_clear_clicked(GtkButton* button, gpointer user_data) {
	GtkFileChooser *filechooser = (GtkFileChooser*) lookup_widget("pref_soft_proofing_profile");
	GtkToggleButton *togglebutton = (GtkToggleButton*) lookup_widget("custom_proofing_profile_active");
	g_mutex_lock(&soft_proof_profile_mutex);

	gtk_file_chooser_unselect_all(filechooser);
	if (com.pref.icc.icc_path_soft_proof && com.pref.icc.icc_path_soft_proof[0] != '\0') {
		g_free(com.pref.icc.icc_path_soft_proof);
		com.pref.icc.icc_path_soft_proof = NULL;
	}
	if (gui.icc.soft_proof)
		cmsCloseProfile(gui.icc.soft_proof);
	gui.icc.soft_proof = NULL;

	g_mutex_unlock(&soft_proof_profile_mutex);
	gtk_toggle_button_set_active(togglebutton, FALSE);
	gtk_widget_set_sensitive((GtkWidget*) togglebutton, FALSE);
	refresh_icc_transforms();
	redraw(REMAP_ALL);
	redraw_previews();
}

void on_custom_monitor_profile_active_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkFileChooser *filechooser = (GtkFileChooser*) lookup_widget("pref_custom_monitor_profile");
	gboolean no_file = FALSE;
	gboolean active = gtk_toggle_button_get_active(button);
	g_mutex_lock(&monitor_profile_mutex);
	if (gui.icc.monitor) {
		cmsCloseProfile(gui.icc.monitor);
	}
	if (active) {
		if (!com.pref.icc.icc_path_monitor || com.pref.icc.icc_path_monitor[0] == '\0') {
			com.pref.icc.icc_path_monitor = g_strdup(gtk_file_chooser_get_filename(filechooser));
		}
		if (!com.pref.icc.icc_path_monitor || com.pref.icc.icc_path_monitor[0] == '\0') {
			siril_log_color_message(_("Error: no filename specfied for custom monitor profile.\n"), "red");
			no_file = TRUE;
		} else {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc.icc_path_monitor, "r");
		}
		if (gui.icc.monitor) {
			siril_log_message(_("Monitor profile loaded from %s\n"), com.pref.icc.icc_path_monitor, "r");
			g_mutex_unlock(&monitor_profile_mutex);
			refresh_icc_transforms();
			return;
		} else {
			if (!no_file) {
				siril_log_color_message(_("Monitor profile could not be loaded from %s\n"), "red", com.pref.icc.icc_path_monitor);
			}
			gui.icc.monitor = srgb_monitor_perceptual();
			if (gui.icc.monitor) {
				siril_log_message(_("Monitor ICC profile set to sRGB (D65 whitepoint, gamma = 2.2)\n"));
			} else {
				siril_log_color_message(_("Fatal error: standard sRGB ICC profile could not be loaded.\n"), "red");
				exit(1);
			}
		}
	} else {
		gui.icc.monitor = srgb_monitor_perceptual();
		if (gui.icc.monitor) {
			siril_log_message(_("Monitor ICC profile set to sRGB (D65 whitepoint, gamma = 2.2)\n"));
		} else {
			siril_log_color_message(_("Fatal error: standard sRGB ICC profile could not be loaded.\n"), "red");
			exit(1);
		}
	}
	g_mutex_unlock(&monitor_profile_mutex);
	refresh_icc_transforms();
}

void on_custom_proofing_profile_active_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkFileChooser *filechooser = (GtkFileChooser*) lookup_widget("pref_soft_proofing_profile");
	gboolean no_file = FALSE;
	gboolean active = gtk_toggle_button_get_active(button);
	g_mutex_lock(&soft_proof_profile_mutex);
	if (gui.icc.soft_proof) {
		cmsCloseProfile(gui.icc.soft_proof);
		gui.icc.soft_proof = NULL;
	}
	if (active) {
		if (!com.pref.icc.icc_path_soft_proof || com.pref.icc.icc_path_soft_proof[0] == '\0') {
			com.pref.icc.icc_path_soft_proof = g_strdup(gtk_file_chooser_get_filename(filechooser));
		}
		if (!com.pref.icc.icc_path_soft_proof || com.pref.icc.icc_path_soft_proof[0] == '\0') {
			siril_log_color_message(_("Error: no filename specfied for output device proofing profile.\n"), "red");
			no_file = TRUE;
		} else {
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc.icc_path_soft_proof, "r");
		}
		if (gui.icc.soft_proof) {
			siril_log_message(_("Output device proofing profile loaded from %s\n"), com.pref.icc.icc_path_soft_proof);
			g_mutex_unlock(&soft_proof_profile_mutex);
			refresh_icc_transforms();
			queue_redraw(REMAP_ALL);
			return;
		} else {
			if (!no_file) {
				siril_log_color_message(_("Output device proofing profile could not be loaded from %s\n"), "red", com.pref.icc.icc_path_soft_proof);
			}
			siril_log_color_message(_("Soft proofing is not available while no soft proofing ICC profile is loaded.\n"), "salmon");
		}
	} else {
		siril_log_message(_("Output device proofing profile deactivated. Soft proofing will proof to the monitor profile.\n"));
	}
	g_mutex_unlock(&soft_proof_profile_mutex);
	refresh_icc_transforms();
	queue_redraw(REMAP_ALL);
}

void on_pref_icc_assign_never_toggled(GtkToggleButton *button, gpointer user_data);

void on_pref_icc_assign_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *never = (GtkToggleButton*) lookup_widget("pref_icc_assign_never");
	if (gtk_toggle_button_get_active(button)) {
		g_signal_handlers_block_by_func(never, on_pref_icc_assign_never_toggled, NULL);
		gtk_toggle_button_set_active(never, FALSE);
		g_signal_handlers_unblock_by_func(never, on_pref_icc_assign_never_toggled, NULL);
	}
}

void on_pref_icc_assign_never_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkToggleButton *load = (GtkToggleButton*) lookup_widget("pref_icc_assign_on_load");
	GtkToggleButton *stack = (GtkToggleButton*) lookup_widget("pref_icc_assign_on_stack");
	GtkToggleButton *stretch = (GtkToggleButton*) lookup_widget("pref_icc_assign_on_stretch");
	GtkToggleButton *composition = (GtkToggleButton*) lookup_widget("pref_icc_assign_on_composition");
	if (gtk_toggle_button_get_active(button)) {
		g_signal_handlers_block_by_func(load, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_block_by_func(stack, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_block_by_func(stretch, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_block_by_func(composition, on_pref_icc_assign_toggled, NULL);
		gtk_toggle_button_set_active(load, FALSE);
		gtk_toggle_button_set_active(stack, FALSE);
		gtk_toggle_button_set_active(stretch, FALSE);
		gtk_toggle_button_set_active(composition, FALSE);
		g_signal_handlers_unblock_by_func(load, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_unblock_by_func(stack, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_unblock_by_func(stretch, on_pref_icc_assign_toggled, NULL);
		g_signal_handlers_unblock_by_func(composition, on_pref_icc_assign_toggled, NULL);
	}
}
//////// GUI callbacks for the color management dialog

void on_icc_cancel_clicked(GtkButton* button, gpointer* user_data) {
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_target_profile_label");
	GtkLabel* label2 = (GtkLabel*) lookup_widget("icc_target_mfr_label");
	GtkLabel* label3 = (GtkLabel*) lookup_widget("icc_target_copyright_label");
	GtkLabel* label4 = (GtkLabel*) lookup_widget("icc_current_profile_label");
	GtkLabel* label5 = (GtkLabel*) lookup_widget("icc_mfr_label");
	GtkLabel* label6 = (GtkLabel*) lookup_widget("icc_copyright_label");
	GtkComboBox* target_combo = (GtkComboBox*) lookup_widget("icc_target_combo");
	GtkFileChooser* filechooser = (GtkFileChooser*) lookup_widget("icc_target_filechooser");
	gtk_combo_box_set_active(target_combo, 0);
	gtk_file_chooser_unselect_all(filechooser);
	if (target) {
		cmsCloseProfile(target);
		target = NULL;
	}
	gtk_label_set_text(label, "");
	gtk_label_set_text(label2, "");
	gtk_label_set_text(label3, "");
	gtk_label_set_text(label4, "");
	gtk_label_set_text(label5, "");
	gtk_label_set_text(label6, "");
siril_close_dialog("icc_dialog");
}

void on_icc_assign_clicked(GtkButton* button, gpointer* user_data) {
	if (!target) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No color profile chosen, nothing to assign."));
		return;
	}
	on_clear_roi();
	// We save the undo state as dealing with gfit
	undo_save_state(&gfit, _("Color profile assignment"));

	cmsUInt32Number gfit_colorspace_channels = gfit.naxes[2];
	cmsUInt32Number target_colorspace = cmsGetColorSpace(target);
	cmsUInt32Number target_colorspace_channels = cmsChannelsOf(target_colorspace);

	// Handle initial assignment of an ICC profile
	if (!gfit.color_managed || !gfit.icc_profile) {
		if (gfit.icc_profile) {
			cmsCloseProfile(gfit.icc_profile);
			gfit.icc_profile = NULL;
		}
		if (gfit_colorspace_channels != target_colorspace_channels) {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("Error"), _("Image number of channels does not match color profile number of channels. Cannot assign this profile to this image."));
			return;
		}
		goto FINISH;
	}

	cmsUInt32Number gfit_colorspace = cmsGetColorSpace(gfit.icc_profile);
	gfit_colorspace_channels = cmsChannelsOf(gfit_colorspace);

	if (target_colorspace != cmsSigGrayData && target_colorspace != cmsSigRgbData) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Color space not supported"), _("Siril only supports representing the image in Gray or RGB color spaces at present. You cannot assign or convert to non-RGB color profiles"));
		return;
	}
	if (gfit_colorspace_channels != target_colorspace_channels) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Transform not supported"), _("Image cannot be assigned a color profile with a different number of channels to its current color profile"));
		return;
	}
	if (gfit.icc_profile) {
		cmsCloseProfile(gfit.icc_profile);
		gfit.icc_profile = NULL;
	}
FINISH:
	gfit.icc_profile = copyICCProfile(target);
	if (gfit.icc_profile)
		color_manage(&gfit, TRUE);
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	gtk_widget_set_sensitive(lookup_widget("icc_remove"), gfit.color_managed);
	set_source_information();
	refresh_icc_transforms();
	notify_gfit_modified();
}

void on_icc_remove_clicked(GtkButton* button, gpointer* user_data) {
	on_clear_roi();
	// We save the undo state as dealing with gfit
	undo_save_state(&gfit, _("Color profile removal"));
	if (gfit.icc_profile) {
		cmsCloseProfile(gfit.icc_profile);
		gfit.icc_profile = NULL;
	}
	color_manage(&gfit, FALSE);
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	gtk_widget_set_sensitive(lookup_widget("icc_remove"), gfit.color_managed);
	set_source_information();
	refresh_icc_transforms();
	notify_gfit_modified();

}

void on_icc_convertto_clicked(GtkButton* button, gpointer* user_data) {
	if (!target) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("No color profile chosen, nothing to assign."));
		return;
	}
	on_clear_roi();
	if (!gfit.color_managed || !gfit.icc_profile) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No color profile set"), _("The current image has no color profile. You need to assign one first."));
		return;
	}

	cmsUInt32Number target_colorspace = cmsGetColorSpace(target);

	if (target_colorspace != cmsSigGrayData && target_colorspace != cmsSigRgbData) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Color space not supported"), _("Siril only supports representing the image in Gray or RGB color spaces at present. You cannot assign or convert to non-RGB color profiles"));
		return;
	}

	// Do the transform
	// We save the undo state if dealing with gfit
	undo_save_state(&gfit, _("Color profile conversion"));
	cmsUInt32Number temp_intent = com.pref.icc.processing_intent;
	com.pref.icc.processing_intent = com.pref.icc.export_intent;
	siril_colorspace_transform(&gfit, target);
	com.pref.icc.processing_intent = temp_intent;

	// Assign the new color space to gfit
	if (gfit.icc_profile)
		cmsCloseProfile(gfit.icc_profile);
	gfit.icc_profile = copyICCProfile(target);
	if (gfit.icc_profile)
		color_manage(&gfit, TRUE);
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	set_source_information();
	refresh_icc_transforms();
	close_tab();
	init_right_tab();
	notify_gfit_modified();
}

void on_icc_target_combo_changed(GtkComboBox* combo, gpointer* user_data) {
	GtkFileChooser* filechooser = (GtkFileChooser*) lookup_widget("icc_target_filechooser");
	internal_icc target_index = gtk_combo_box_get_active(combo);
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_target_profile_label");
	gtk_file_chooser_unselect_all(filechooser);
	if (target) {
		cmsCloseProfile(target);
		target = NULL;
	}
	GtkLabel *mfr_label = NULL, *copyright_label = NULL;
	switch (target_index) {
		case NONE:
			mfr_label = (GtkLabel*) lookup_widget("icc_target_mfr_label");
			copyright_label = (GtkLabel*) lookup_widget("icc_target_copyright_label");
			gtk_label_set_text(label, "");
			gtk_label_set_text(copyright_label, "");
			gtk_label_set_text(mfr_label, "");
			return;
		case SRGB_LINEAR:
			target = srgb_linear();
			break;
		case SRGB_TRC:
			target = srgb_trc();
			break;
		case REC2020_LINEAR:
			target = rec2020_linear();
			break;
		case REC2020_TRC:
			target = rec2020_trc();
			break;
		case GRAY_LINEAR:
			target = gray_linear();
			break;
		case GRAY_SRGBTRC:
			target = gray_srgbtrc();
			break;
		case GRAY_REC709TRC:
			target = gray_rec709trc();
			break;
	}
	set_target_information();
}

void on_icc_target_filechooser_file_set(GtkFileChooser* filechooser, gpointer* user_data) {
	if (target) {
		cmsCloseProfile(target);
		target = NULL;
	}
	GtkComboBox* target_combo = (GtkComboBox*) lookup_widget("icc_target_combo");
	g_signal_handlers_block_by_func(target_combo, on_icc_target_combo_changed, NULL);
	gtk_combo_box_set_active(target_combo, 0);
	g_signal_handlers_unblock_by_func(target_combo, on_icc_target_combo_changed, NULL);
	gchar *filename = siril_file_chooser_get_filename(filechooser);
	if (filename) {
		target = cmsOpenProfileFromFile(filename, "r");
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not load file"), _("Error getting filename from widget."));
		gtk_file_chooser_unselect_all(filechooser);
		return;
	}
	if (!target) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not load file"), _("Error loading selected file or it does not contain a valid ICC profile."));
		gtk_file_chooser_unselect_all(filechooser);
	} else {
		cmsColorSpaceSignature target_signature = cmsGetColorSpace(target);
		if (target_signature == cmsSigRgbData) {
			set_target_information();
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Profile error"), _("Error: profile does not describe a RGB color space. Using non-RGB color spaces (e.g. CIE La*b*, XYZ, HSL) as image working spaces is not supported, though these color spaces may be used internally for some operations."));
			cmsCloseProfile(target);
			target = NULL;
			gtk_file_chooser_unselect_all(filechooser);
		}
	}
	g_free(filename);
}

gboolean on_icc_main_window_button_clicked(GtkWidget *btn, GdkEventButton *event, gpointer userdata) {
	if (!(single_image_is_loaded() || sequence_is_loaded()))
		return FALSE;
	if (event->type == GDK_BUTTON_PRESS  &&  event->button == 3) {
		// Right mouse button press
        if (gui.icc.iso12646) {
			siril_debug_print("Disabling approximate ISO12646 viewing conditions\n");
			disable_iso12646_conditions(TRUE, TRUE, TRUE);
		} else {
			siril_debug_print("Enabling approximate ISO12646 viewing conditions\n");
			enable_iso12646_conditions();
		}
        return TRUE;
	}
	if (event->type == GDK_BUTTON_PRESS  &&  event->button == 1) {
		// Left mouse button press
		siril_open_dialog("icc_dialog");
		return TRUE;
	}
	return FALSE;
}

void on_icc_dialog_show(GtkWidget *dialog, gpointer user_data) {
	set_source_information();
	set_target_information();
	GtkFileChooser* fc = (GtkFileChooser*) lookup_widget("icc_target_filechooser");
	gtk_file_chooser_set_current_folder(fc, default_system_icc_path());
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	gtk_widget_set_sensitive(lookup_widget("icc_remove"), gfit.color_managed);
}

void on_icc_export_clicked(GtkButton *button, gpointer user_data) {
	control_window_switch_to_tab(OUTPUT_LOGS);
	if (!gfit.icc_profile) {
		siril_log_color_message(_("Error: no ICC profile associated with the current image. Cannot export.\n"), "red");
		return;
	}
	char *filename = NULL;
	if (com.uniq->filename && com.uniq->filename[0] != '\0')
		filename = remove_ext_from_filename(com.uniq->filename);
	else
		filename = strdup("image");
	gchar* temp = g_strdup_printf("%s.icc", filename);
	free(filename);
	filename = strdup(temp);
	g_free(temp);
	if (cmsSaveProfileToFile(gfit.icc_profile, filename))
		siril_log_color_message(_("Exported ICC profile to %s\n"), "green", filename);
	else
		siril_log_color_message(_("Failed to export ICC profile to %s\n"), "red", filename);
	free(filename);
}

void on_icc_export_builtin_clicked(GtkButton *button, gpointer user_data) {
	export_elle_stone_profiles();
}

void on_icc_plot_clicked(GtkButton *button, gpointer user_data) {
	if (gfit.icc_profile && siril_color_profile_is_rgb (gfit.icc_profile)) {
		siril_plot_colorspace(gfit.icc_profile, TRUE);
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Chromaticity plot only works with RGB color profiles"));
	}
}

static gboolean colorspace_comparison_image_set = FALSE;

void on_icc_gamut_visualisation_clicked() {
	GtkWidget *win = lookup_widget("icc_gamut_dialog");
	gtk_window_set_transient_for(GTK_WINDOW(win), GTK_WINDOW(lookup_widget("settings_window")));
	if (!colorspace_comparison_image_set) {
		GError *error = NULL;
		GdkPixbuf *pixbuf = gdk_pixbuf_new_from_resource("/org/siril/ui/pixmaps/CIE1931.svg", &error);
		if (error) {
			siril_debug_print("Error: %s\n", error->message);
			g_error_free(error);
		}
		gtk_image_set_from_pixbuf(GTK_IMAGE(lookup_widget("colorspace_comparison")), pixbuf);
		colorspace_comparison_image_set = TRUE;
	}
	/* Here this is wanted that we do not use siril_open_dialog */
	gtk_widget_show(win);
}

void on_icc_gamut_close_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *win = lookup_widget("icc_gamut_dialog");
	gtk_widget_hide(win);
}

gboolean iso_12646_draw_event(GtkWidget *widget, cairo_t *cr, gpointer data) {
	if (!(single_image_is_loaded() || (sequence_is_loaded())))
		return FALSE;
	GtkWidget *parent_widget = GTK_WIDGET(data);
	GtkWidget *red = lookup_widget("vbox_r");
	GtkWidget *green = lookup_widget("vbox_g");
	GtkWidget *blue = lookup_widget("vbox_b");
	GtkWidget *rgb = lookup_widget("vbox_rgb");
	GtkWidget *child_widget;
	if (parent_widget == red) {
		child_widget = lookup_widget("drawingarear");
	} else if (parent_widget == green) {
		child_widget = lookup_widget("drawingareag");
	} else if (parent_widget == blue) {
		child_widget = lookup_widget("drawingareab");
	} else if (parent_widget == rgb) {
		child_widget = lookup_widget("drawingareargb");
	} else {
		return FALSE;
	}

	// Get child widget's allocation
	GtkAllocation child_allocation;
	gtk_widget_get_allocation(child_widget, &child_allocation);

	// Translate child widget's coordinates to parent widget's drawing context
	gint white_border_width = 36;
	gdouble z = get_zoom_val();
	gdouble x0 = 0.0, y0 = 0.0;
	gint child_x0, child_y0, child_x1, child_y1;
	cairo_matrix_transform_point(&gui.display_matrix, &x0, &y0);
	gtk_widget_translate_coordinates(child_widget, parent_widget, (gint) x0, (gint) y0, &child_x0, &child_y0);
	gtk_widget_translate_coordinates(child_widget, parent_widget, (gint) ((z * gfit.rx) + white_border_width), (gint) ((z * gfit.ry) + white_border_width), &child_x1, &child_y1);

	// Get parent widget's allocation
	GtkAllocation parent_allocation;
	gtk_widget_get_allocation(parent_widget, &parent_allocation);

	// Fill parent widget's background with mid grey
	cairo_set_source_rgb(cr, 0.7647, 0.7647, 0.7647); // D50 Gray as recommended by ISO 12646
	cairo_rectangle(cr, 0, 0, parent_allocation.width, parent_allocation.height);
	cairo_fill(cr);

	// Draw white rectangle around the child widget with a margin of {white_border_width} pixels
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0); // White color
	cairo_rectangle(cr, child_x0 - white_border_width, child_y0 - white_border_width,
					child_x1 + white_border_width, child_y1 + white_border_width - 16);
	cairo_fill(cr);
	// Return FALSE to propagate the draw event further
	return FALSE;
}

gboolean on_iso12646_panel_hide_completed(GtkWidget *widget, gpointer data) {
	gboolean *remap = (gboolean *) data;
	// Set the zoom value to 80% of the zoom-to-fit value
	update_zoom_fit_button();
	gui.zoom_value = -1.0;
	double z = 0.8 * get_zoom_val();
	gui.zoom_value = z;
	// Set the display offset. This has to be calculated from the width of the control
	// window less that of the pane button because the width of the drawarea won't
	// have updated in time if we have to hide the side pane.
	int window_width = gtk_widget_get_allocated_width(gui.view[RED_VPORT].drawarea);
	int window_height = gtk_widget_get_allocated_height(gui.view[RED_VPORT].drawarea);
	gui.display_offset.x = (window_width / 2 - (gfit.rx / 2) * z);
	gui.display_offset.y = (window_height / 2) - (gfit.ry / 2 * z);
	adjust_vport_size_to_image();
	queue_redraw(remap ? REMAP_ALL : REDRAW_IMAGE); // Has to do the remap in case the sliders have been changed
	gtk_widget_queue_draw(lookup_widget("vbox_r"));
	gtk_widget_queue_draw(lookup_widget("vbox_g"));
	gtk_widget_queue_draw(lookup_widget("vbox_b"));
	gtk_widget_queue_draw(lookup_widget("vbox_rgb"));
	return FALSE;
}

static gboolean panel_state = TRUE;
static double prior_zoom = -1;
static point prior_offset = { 0.0 , 0.0 };
static sliders_mode prior_sliders = MINMAX;
static WORD prior_lo = 0, prior_hi = 65535;
static display_mode prior_rendering_mode = LINEAR_DISPLAY;
static gboolean mode_changed = FALSE;

// This function overrides any GTK theme to assign a white border and D50 Gray
// background to the 4 vports to approximate ISO 12646 viewing conditions.
// It is recommended in conjunction with the excellent Equilux GTK theme.
void enable_iso12646_conditions() {
	// Cache prior state
	prior_zoom = get_zoom_val();
	memcpy(&prior_offset, &gui.display_offset, sizeof(point));
	prior_sliders = gui.sliders;
	prior_lo = gui.lo;
	prior_hi = gui.hi;
	prior_rendering_mode = gui.rendering_mode;
	if (prior_rendering_mode != LINEAR_DISPLAY)
		mode_changed = TRUE;
	// Add draw callbacks
	GtkWidget *parent_widget = lookup_widget("vbox_rgb");
	gui.icc.sh_rgb = g_signal_connect(G_OBJECT(parent_widget), "draw", G_CALLBACK(iso_12646_draw_event), parent_widget);
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_b");
	gui.icc.sh_b = g_signal_connect(G_OBJECT(parent_widget), "draw", G_CALLBACK(iso_12646_draw_event), parent_widget);
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_g");
	gui.icc.sh_g = g_signal_connect(G_OBJECT(parent_widget), "draw", G_CALLBACK(iso_12646_draw_event), parent_widget);
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_r");
	gui.icc.sh_r = g_signal_connect(G_OBJECT(parent_widget), "draw", G_CALLBACK(iso_12646_draw_event), parent_widget);
	gtk_widget_queue_draw(parent_widget);
	// Set the text color of the labels within the vbox
	PangoAttrList *attrs = pango_attr_list_new();
	PangoAttribute *attr = pango_attr_foreground_new(15 * PANGO_SCALE, 15 * PANGO_SCALE, 15 * PANGO_SCALE);
	pango_attr_list_insert(attrs, attr);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_red")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_green")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_blue")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_rgb")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_red")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_green")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_blue")), attrs);
	gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_rgb")), attrs);
	// Get the panel out of the way
	GtkImage *image = GTK_IMAGE(gtk_bin_get_child(GTK_BIN(GTK_BUTTON(lookup_widget("button_paned")))));
	gtk_image_set_from_icon_name(image, "pan-start-symbolic", GTK_ICON_SIZE_BUTTON);
	GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
	GtkWidget *widget = gtk_paned_get_child2(paned);
	panel_state = gtk_widget_get_visible(widget);
	if (panel_state)
		gtk_widget_set_visible(widget, FALSE);
	// Set the sliders to min/max
	gboolean is_8bit = gfit.orig_bitpix == BYTE_IMG;
	gboolean remap = ((gui.lo == 0 && gui.hi == 65535) || (is_8bit && (gui.lo == 0 && gui.hi == 255))) || mode_changed;
	gui.sliders = USER;
	gui.lo = 0;
	gui.hi = is_8bit ? 255 : 65535;
	gui.rendering_mode = LINEAR_DISPLAY;
	set_display_mode();
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_user")), TRUE);
	set_cutoff_sliders_values(); // The redraw will happen in the idle
	gui.icc.iso12646 = TRUE;
	g_idle_add((GSourceFunc)on_iso12646_panel_hide_completed, &remap);
}

void disable_iso12646_conditions(gboolean revert_zoom, gboolean revert_panel, gboolean revert_rendering_mode) {
	GtkWidget *parent_widget = lookup_widget("vbox_rgb");
	if (gui.icc.sh_rgb)
		g_signal_handler_disconnect(G_OBJECT(parent_widget), gui.icc.sh_rgb);
	gui.icc.sh_rgb = 0;
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_r");
	if (gui.icc.sh_r)
		g_signal_handler_disconnect(G_OBJECT(parent_widget), gui.icc.sh_r);
	gui.icc.sh_r = 0;
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_g");
	if (gui.icc.sh_g)
		g_signal_handler_disconnect(G_OBJECT(parent_widget), gui.icc.sh_g);
	gui.icc.sh_g = 0;
	gtk_widget_queue_draw(parent_widget);
	parent_widget = lookup_widget("vbox_b");
	if (gui.icc.sh_b)
		g_signal_handler_disconnect(G_OBJECT(parent_widget), gui.icc.sh_b);
	gui.icc.sh_b = 0;
	// Revert the text color of the labels within the vbox
    PangoAttrList *attrs = NULL;
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_red")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_green")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_blue")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelfilename_rgb")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_red")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_green")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_blue")), attrs);
    gtk_label_set_attributes(GTK_LABEL(lookup_widget("labelwcs_rgb")), attrs);
	gui.icc.iso12646 = FALSE;
	if (revert_panel) {
		// Return the panel to its previous state
		GtkImage *image = GTK_IMAGE(gtk_bin_get_child(GTK_BIN(GTK_BUTTON(lookup_widget("button_paned")))));
		if (!panel_state)
			gtk_image_set_from_icon_name(image, "pan-end-symbolic", GTK_ICON_SIZE_BUTTON);
		GtkPaned *paned = GTK_PANED(lookup_widget("main_panel"));
		GtkWidget *widget = gtk_paned_get_child2(paned);
		if (panel_state)
			gtk_widget_set_visible(widget, TRUE);
	}
	if (revert_zoom)
		gui.zoom_value = prior_zoom;
	memcpy(&gui.display_offset, &prior_offset, sizeof(point));
	gui.sliders = prior_sliders;
	gui.lo = prior_lo;
	gui.hi = prior_hi;
	if (gui.sliders == MINMAX)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_minmax")), TRUE);
	else if (gui.sliders == MIPSLOHI)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_hilo")), TRUE);
	set_cutoff_sliders_values();
	if (revert_rendering_mode) {
		gui.rendering_mode = prior_rendering_mode;
		set_display_mode();
	}
	if (mode_changed)
		redraw(REMAP_ALL);
	gtk_widget_queue_draw(lookup_widget("control_window"));
}

void siril_plot_colorspace(cmsHPROFILE profile, gboolean compare_srgb) {
	cmsCIEXYZTRIPLE XYZtriple = { 0 };
	cmsCIEXYZ whitepoint = { 0 };
	cmsCIExyY redxyY, greenxyY, bluexyY, whitexyY;
	char *description = NULL;
	int length = cmsGetProfileInfoASCII(profile, cmsInfoDescription, "en", "US", NULL, 0);
	if (length) {
		description = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(profile, cmsInfoDescription, "en", "US", description, length);
	}

	if (!(siril_color_profile_is_rgb (profile))) {
		siril_log_message(_("This ICC profile is not RGB. Unable to plot the colorspace.\n"));
		return;
	}
	if (! siril_color_profile_get_rgb_matrix_colorants (profile, &XYZtriple, &whitepoint)) {
		siril_log_message(_("Error reading chromaticities\n"));
		return;
	}

	cmsXYZ2xyY(&redxyY, &XYZtriple.Red);
	cmsXYZ2xyY(&greenxyY, &XYZtriple.Green);
	cmsXYZ2xyY(&bluexyY, &XYZtriple.Blue);
	cmsXYZ2xyY(&whitexyY, &whitepoint);
	double white_x = whitexyY.x;
	double white_y = whitexyY.y;
	double *horseshoe_x = malloc(322 * sizeof(double)), *horseshoe_y = malloc(322 * sizeof(double));
	for (int i = 0 ; i < 321 ; i++) {
		double w = 380 + i;
		cmsCIEXYZ XYZ = { x1931(w), y1931(w), z1931(w)};
		cmsCIExyY xyY;
		cmsXYZ2xyY(&xyY, &XYZ);
		horseshoe_x[i] = xyY.x;
		horseshoe_y[i] = xyY.y;
	}

	horseshoe_x[321] = horseshoe_x[0];
	horseshoe_y[321] = horseshoe_y[0];
	double colorspace_x[4] = {redxyY.x, greenxyY.x, bluexyY.x, redxyY.x};
	double colorspace_y[4] = {redxyY.y, greenxyY.y, bluexyY.y, redxyY.y};
	double srgb_x[4] = {0.639998686, 0.300003784, 0.150002046, 0.639998686};
	double srgb_y[4] = {0.330010138, 0.600003357, 0.059997204, 0.330010138};
	siril_plot_data *spl_data = NULL;

	gchar *title1 = g_strdup_printf(_("Source Color Profile Chromaticity Diagram\n"
					"<span size=\"small\">"
					"%s"
					"</span>"), description);
	free(description);
	spl_data = malloc(sizeof(siril_plot_data));
	init_siril_plot_data(spl_data);
	siril_plot_set_xlabel(spl_data, _("CIE x"));
	siril_plot_set_savename(spl_data, "color_profile");
	siril_plot_set_title(spl_data, title1);
	siril_plot_set_ylabel(spl_data, _("CIE y"));
	int n = 1;
	siril_plot_add_xydata(spl_data, _("Color profile"), 4, colorspace_x, colorspace_y, NULL, NULL);
	siril_plot_set_nth_plot_type(spl_data, n, KPLOT_LINES);
	siril_plot_set_nth_color(spl_data, n, (double[3]) { 0.0, 0.5, 1.0 } );
	n++;
	siril_plot_add_xydata(spl_data, _("Color profile whitepoint"), 1, &white_x, &white_y, NULL, NULL);
	siril_plot_set_nth_plot_type(spl_data, n, KPLOT_POINTS);
	siril_plot_set_nth_color(spl_data, n, (double[3]) { 0.0, 0.5, 1.0 } );
	n++;
	siril_plot_add_xydata(spl_data, _("CIE 1931"), 322, horseshoe_x, horseshoe_y, NULL, NULL);
	siril_plot_set_nth_plot_type(spl_data, n, KPLOT_LINES);
	siril_plot_set_nth_color(spl_data, n, (double[3]) { 0.0, 0.0, 0.0 } );
	n++;
	if (!siril_plot_set_background(spl_data, "CIE1931xy.svg"))
		siril_log_color_message(_("Could not load background\n"), "red");
	if (compare_srgb) {
		siril_plot_add_xydata(spl_data, _("sRGB"), 4, srgb_x, srgb_y, NULL, NULL);
		siril_plot_set_nth_plot_type(spl_data, n, KPLOT_LINES);
	}
	spl_data->datamin = (point) { 0.0, 0.0 };
	spl_data->datamax = (point) { 0.8, 0.9 };
	spl_data->width = 600;
	spl_data->height = 600;
	spl_data->cfgdata.line.sz = 2;

	free(horseshoe_x);
	free(horseshoe_y);

	siril_add_idle(create_new_siril_plot_window, spl_data);
	siril_add_idle(end_generic, NULL);
}
