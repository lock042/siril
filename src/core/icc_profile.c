/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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

#include <glib.h>
#include "algos/lcms_acceleration/lcms2_fast_float.h"
#include "algos/lcms_acceleration/lcms2_threaded.h"
#include "core/siril.h"
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
#include "gui/utils.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "core/siril_log.h"
#include "core/siril_app_dirs.h"
#include "core/proto.h"

static GMutex monitor_profile_mutex;
static GMutex soft_proof_profile_mutex;
static GMutex default_profiles_mutex;

static cmsHPROFILE target = NULL; // Target profile for the GUI tool

////// Functions //////

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

void color_manage(fits *fit, gboolean active) {
	fit->color_managed = active;
	if (fit == &gfit && !com.headless) {
		gchar *name = g_build_filename(siril_get_system_data_dir(), "pixmaps", active ? "color_management.svg" : "color_management_off.svg", NULL);
		gchar *tooltip = NULL;
		if (active) {
			if (fit->icc_profile) {
				int length = cmsGetProfileInfoASCII(fit->icc_profile, cmsInfoDescription, "en", "US", NULL, 0);
				int length2 = cmsGetProfileInfoASCII(gui.icc.monitor, cmsInfoDescription, "en", "US", NULL, 0);
				char *monitor = malloc(length2 * sizeof(char));
				cmsGetProfileInfoASCII(gui.icc.monitor, cmsInfoDescription, "en", "US", monitor, length2);
				// TODO: update the tooltip to show the proofing profile if in soft proof mode
				if (gui.rendering_mode == SOFT_PROOF_DISPLAY) {
				}
				char *buffer = NULL;
				if (length) {
					buffer = malloc(length * sizeof(char));
					cmsGetProfileInfoASCII(fit->icc_profile, cmsInfoDescription, "en", "US", buffer, length);
					tooltip = g_strdup_printf(_("Image is color managed\nImage profile: %s\nMonitor profile: %s"), buffer, monitor);
					free(buffer);
					free(monitor);
				}
			}
			if (!tooltip)
				tooltip = g_strdup(_("Image is color managed\n"));
		} else {
			tooltip = g_strdup(_("Image is not color managed"));
		}
		GtkWidget *image = lookup_widget("color_managed_icon");
		GtkWidget *button = lookup_widget("icc_main_window_button");
		gtk_image_set_from_file((GtkImage*) image, name);
		gtk_widget_set_tooltip_text(button, tooltip);
		g_free(name);
		g_free(tooltip);
	}
}

static void export_profile(cmsHPROFILE profile) {
	char *filename = NULL;
	int length;
	length = cmsGetProfileInfoASCII(profile, cmsInfoDescription, "en", "US", NULL, 0);
	if (length) {
		filename = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(profile, cmsInfoDescription, "en", "US", filename, length);
	}
	if (!g_str_has_suffix(filename, ".icc")) {
		gchar* temp = g_strdup_printf("%s.icc", filename);
		free(filename);
		filename = strdup(temp);
		g_free(temp);
	}
	if (cmsSaveProfileToFile(profile, filename))
		siril_log_color_message(_("Exported ICC profile to %s\n"), "green", filename);
	else
		siril_log_color_message(_("Failed to export ICC profile to %s\n"), "red", filename);
	free(filename);
}

static void export_elle_stone_profiles() {
	cmsHPROFILE profile;
	control_window_switch_to_tab(OUTPUT_LOGS);

	profile = srgb_linear();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = srgb_trc();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = srgb_trcv2();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = rec2020_linear();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = rec2020_trc();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = rec2020_trcv2();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = gray_linear();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = gray_srgbtrc();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = gray_srgbtrcv2();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = gray_rec709trc();
	export_profile(profile);
	cmsCloseProfile(profile);

	profile = gray_rec709trcv2();
	export_profile(profile);
	cmsCloseProfile(profile);
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

void validate_custom_profiles() {
	if (com.pref.icc.icc_path_monitor && com.pref.icc.icc_path_monitor[0] != '\0') {
		g_mutex_lock(&monitor_profile_mutex);
		if (validate_profile(com.pref.icc.icc_path_monitor)) {
			if (gui.icc.monitor)
				cmsCloseProfile(gui.icc.monitor);
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc.icc_path_monitor, "r");
			if (!gui.icc.monitor) {
				gui.icc.monitor = srgb_trc();
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
		gui.icc.monitor = srgb_trc();
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
	cmsSetAlarmCodes(alarmcodes); // Out of gamut colours will be shown in magenta

	com.icc.rendering_flags |= ((com.pref.icc.rendering_bpc * cmsFLAGS_BLACKPOINTCOMPENSATION) & !(com.pref.icc.rendering_intent == INTENT_ABSOLUTE_COLORIMETRIC));
	gui.icc.proofing_flags = ((com.pref.icc.proofing_bpc * cmsFLAGS_BLACKPOINTCOMPENSATION) & !(com.pref.icc.proofing_intent == INTENT_ABSOLUTE_COLORIMETRIC)) | cmsFLAGS_SOFTPROOFING;

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
	if (gfit.icc_profile == NULL || !gfit.color_managed)
		return NULL;
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

/* This function creates a fallback display transform for when gfit is not
 * color managed but the monitor profile is different to sRGB */
cmsHTRANSFORM fallback_display_transform() {
	g_assert(gui.icc.monitor);
	cmsHTRANSFORM transform = NULL;
	cmsHPROFILE profile = srgb_trc();
	cmsUInt32Number fit_signature = cmsGetColorSpace(profile);
	cmsUInt32Number srctype = get_planar_formatter_type(fit_signature, gfit.type, TRUE);
	g_mutex_lock(&monitor_profile_mutex);
	// The display transform is always single threaded as OpenMP is used within the remap function
	transform = cmsCreateTransformTHR(com.icc.context_single, profile, srctype, gui.icc.monitor, TYPE_RGB_16_PLANAR, com.pref.icc.rendering_intent, com.icc.rendering_flags);
	g_mutex_unlock(&monitor_profile_mutex);
	if (transform == NULL)
		siril_log_message("Error: failed to create display_transform!\n");
	else
		siril_debug_print("Fallback display transform created (sRGB to monitor)\n");
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
	if (gfit.icc_profile == NULL || gfit.color_managed == FALSE || gui.icc.soft_proof == NULL)
		return NULL;
	cmsUInt32Number flags = gui.icc.proofing_flags;
	if (gui.cut_over) {
		flags |= cmsFLAGS_GAMUTCHECK;
	}
	cmsUInt32Number type = (gfit.naxes[2] == 1 ? TYPE_GRAY_16 : TYPE_RGB_16_PLANAR);
	g_mutex_lock(&soft_proof_profile_mutex);
	g_mutex_lock(&monitor_profile_mutex);
	cmsHPROFILE proofing_transform = cmsCreateProofingTransformTHR(
						com.icc.context_single,
						gfit.icc_profile,
						type,
						gui.icc.monitor,
						type,
						gui.icc.soft_proof,
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
		if (gui.icc.display_transform)
			cmsDeleteTransform(gui.icc.display_transform);
		gui.icc.display_transform = initialize_display_transform();
		if (gui.icc.proofing_transform)
			cmsDeleteTransform(gui.icc.proofing_transform);
		gui.icc.proofing_transform = initialize_proofing_transform();
	}
	check_gfit_profile_identical_to_monitor();
}

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

/* Wrapper for cmsIsToneCurveLinear() that reads the relevant tone curve from
 * fit->icc_profile and checks to see if it is linear
 */
cmsBool fit_icc_is_linear(fits *fit) {
	if (!(fit->color_managed) || fit->icc_profile == NULL)
		return FALSE;
	cmsToneCurve *tonecurve;
	if (fit->naxes[2] == 1) {
		tonecurve = cmsReadTag(fit->icc_profile, cmsSigGrayTRCTag);
	} else {
		tonecurve = cmsReadTag(fit->icc_profile, cmsSigRedTRCTag);
	}
	// If we fail to read a tonecurve then cmsIsToneCurveLinear will crash
	// Return FALSE as a conservative result - remapping will be done
	if (!tonecurve)
		return FALSE;
	return cmsIsToneCurveLinear(tonecurve);
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
		control_window_switch_to_tab(OUTPUT_LOGS);
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
			control_window_switch_to_tab(OUTPUT_LOGS);
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
	if (!a && !b)
		return TRUE;
	if (!a || !b)
		return FALSE;
	cmsUInt8Number *block_a = NULL, *block_b = NULL;
	cmsUInt32Number length_a, length_b;
	cmsBool ret_a, ret_b, retval;
	if (a) {
		ret_a = cmsSaveProfileToMem(a, NULL, &length_a);
	}
	if (b) {
		ret_b = cmsSaveProfileToMem(b, NULL, &length_b);
	}
	// If a profile can't be saved to a buffer or the lengths don't match
	// we can already return FALSE
	if (!ret_a || !ret_b || length_a != length_b)
		return FALSE;
	if (length_a > 0) {
		block_a = malloc(length_a * sizeof(BYTE));
		if (!block_a) {
			PRINT_ALLOC_ERR;
			return FALSE;
		}
		ret_a = cmsSaveProfileToMem(a, (void*) block_a, &length_a);
	}
	if (ret_a) {
		if (length_b > 0) {
			block_b = malloc(length_a * sizeof(BYTE));
			if (!block_b) {
				PRINT_ALLOC_ERR;
				free(block_a);
				return FALSE;
			}
			ret_b = cmsSaveProfileToMem(b, (void*) block_b, &length_b);
		}
		if (!ret_b) {
			free(block_a);
			free(block_b);
			return FALSE;
		}
	} else {
		free(block_a);
		return FALSE;
	}
	retval = (memcmp(block_a, block_b, length_a) == 0) ? TRUE : FALSE;
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
	int length = cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoDescription, "en", "US", NULL, 0);
	char *buffer = NULL;
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoDescription, "en", "US", buffer, length);
		gtk_label_set_text(label, buffer);
		free(buffer);
	}

	// Set manufacturer
	length = cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoManufacturer, "en", "US", NULL, 0);
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoManufacturer, "en", "US", buffer, length);
		gtk_label_set_text(mfr_label, buffer);
		free(buffer);
	}

	// Set copyright
	length = cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoCopyright, "en", "US", NULL, 0);
	if (length) {
		buffer = (char*) malloc(length * sizeof(char));
		cmsGetProfileInfoASCII(gfit.icc_profile, cmsInfoCopyright, "en", "US", buffer, length);
		gtk_label_set_text(copyright_label, buffer);
		free(buffer);
	}
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

	cmsUInt32Number target_colorspace = cmsGetColorSpace(profile);
	cmsUInt32Number target_colorspace_channels = cmsChannelsOf(target_colorspace);
	cmsUInt32Number fit_colorspace_channels;

	// If profile is NULL, we remove the profile from fit to match it. This is an unusual
	// case but the behaviour is consistent.
	if (!profile) {
		if (fit->icc_profile)
			cmsCloseProfile(fit->icc_profile);
		fit->icc_profile = NULL;
		color_manage(fit, FALSE);
		return;
	}

	// If fit->color_managed is FALSE, we assign the profile rather than convert to it
	if (!fit->color_managed || !fit->icc_profile) {
		fit_colorspace_channels = fit->naxes[2];
		if (fit_colorspace_channels == target_colorspace_channels) {
			if (fit->icc_profile)
				cmsCloseProfile(fit->icc_profile);
			fit->icc_profile = copyICCProfile(profile);
			siril_debug_print("siril_colorspace_transform() assigned a profile\n");
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
	cmsHTRANSFORM transform = cmsCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), fit->icc_profile, srctype, profile, desttype, com.pref.icc.processing_intent, com.icc.rendering_flags);
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
		color_manage(fit, TRUE);
		siril_debug_print("siril_colorspace_transform() converted a profile\n");
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Failed to create colorspace transform."));
	}
}

/* This function converts the current image to the working color space. The user
 * has the option to set a preference as to the behaviour to be adopted:
 * 0. Automatically convert the image to the working color space (if it isn't already);
 * 1. Always ask whether to convert to the working color space;
 * 2. Do nothing.
 * The function is intended to be triggered at the start of any GUI stretch function.
 * Scripts are responsible for specifying any color space assignment using icc_assign.
 */
void convert_with_approval(fits *fit) {
	// Scripts are responsible fo managing this themselves
	if (com.script)
		return;

	// If the preference is never to autoconvert, we have nothing to do
	if (com.pref.icc.autoconversion == ICC_NEVER_AUTOCONVERT)
		return;

	// If the image is already in the working color space, we have nothing to do
	if (fit->color_managed && profiles_identical(fit->icc_profile, fit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.working_standard))
		return;

	gboolean proceed = FALSE;

	// If the preference is that we always convert, trigger the conversion
	if (com.pref.icc.autoconversion == ICC_ALWAYS_AUTOCONVERT)
		proceed = TRUE;

	else if (com.pref.icc.autoconversion == ICC_ASK_TO_CONVERT) {
		if (!fit->color_managed) {
			proceed = siril_confirm_dialog(_("Recommend color space assignment"), _("The current image is not color managed. It looks like you're about to stretch the image: do you want to assign your nonlinear working color space now? (Recommended!)"), _("Assign"));
		} else {
			proceed = siril_confirm_dialog(_("Recommend color space conversion"), _("The current image has a linear ICC profile. It looks like you're about to stretch the image: do you want to convert it to your nonlinear working color space now? (Recommended!)"), _("Convert"));
		}
	}

	// If the conversion has been triggered automatically or via the dialog, do it.
	if (proceed) {
		set_cursor_waiting(TRUE);
		// siril_colorspace_transform takes care of hitherto non-color managed images, and assigns a profile instead of converting them
		siril_colorspace_transform(fit, (fit->naxes[2] == 1 ? com.icc.mono_standard : com.icc.working_standard));
		set_source_information();
		refresh_icc_transforms();
		notify_gfit_modified();
		set_cursor_waiting(FALSE);
	}
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
	siril_log_color_message(_("Error: the custom workspace profiles are not all set and / or valid. Defaulting to sRGB.\n"), "red");
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
	g_mutex_lock(&monitor_profile_mutex);
	if (com.pref.icc.icc_path_monitor && com.pref.icc.icc_path_monitor[0] != '\0') {
		g_free(com.pref.icc.icc_path_monitor);
		com.pref.icc.icc_path_monitor = NULL;
		cmsCloseProfile(gui.icc.monitor);
		gui.icc.monitor = srgb_trc();
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
	refresh_icc_transforms();
}

void on_proofing_profile_clear_clicked(GtkButton* button, gpointer user_data) {
	GtkFileChooser *filechooser = (GtkFileChooser*) lookup_widget("pref_soft_proofing_profile");
	GtkToggleButton *togglebutton = (GtkToggleButton*) lookup_widget("custom_proofing_profile_active");
	g_mutex_lock(&soft_proof_profile_mutex);

	gtk_file_chooser_unselect_all(filechooser);
	if (com.pref.icc.icc_path_soft_proof && com.pref.icc.icc_path_soft_proof[0] != '\0') {
		g_free(com.pref.icc.icc_path_soft_proof);
		com.pref.icc.icc_path_soft_proof = NULL;
		if (gui.icc.soft_proof)
			cmsCloseProfile(gui.icc.soft_proof);
		gui.icc.soft_proof = NULL;
	}
	g_mutex_unlock(&soft_proof_profile_mutex);
	gtk_toggle_button_set_active(togglebutton, FALSE);
	gtk_widget_set_sensitive((GtkWidget*) togglebutton, FALSE);
	refresh_icc_transforms();

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
			gui.icc.monitor = srgb_trc();
			if (gui.icc.monitor) {
				siril_log_message(_("Monitor ICC profile set to sRGB (D65 whitepoint, gamma = 2.2)\n"));
			} else {
				siril_log_color_message(_("Fatal error: standard sRGB ICC profile could not be loaded.\n"), "red");
				exit(1);
			}
		}
	} else {
		gui.icc.monitor = srgb_trc();
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
	}
	if (active) {
		if (!com.pref.icc.icc_path_soft_proof || com.pref.icc.icc_path_soft_proof[0] == '\0') {
			com.pref.icc.icc_path_soft_proof = g_strdup(gtk_file_chooser_get_filename(filechooser));
		}
		if (!com.pref.icc.icc_path_soft_proof || com.pref.icc.icc_path_soft_proof[0] == '\0') {
			siril_log_color_message(_("Error: no filename specfied for custom proofing profile.\n"), "red");
			no_file = TRUE;
		} else {
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc.icc_path_soft_proof, "r");
		}
		if (gui.icc.soft_proof) {
			siril_log_message(_("Soft proofing profile loaded from %s\n"), com.pref.icc.icc_path_soft_proof);
			g_mutex_unlock(&soft_proof_profile_mutex);
			refresh_icc_transforms();
			return;
		} else {
			if (!no_file) {
				siril_log_color_message(_("Soft proofing profile could not be loaded from %s\n"), "red", com.pref.icc.icc_path_soft_proof);
			}
			siril_log_color_message(_("Soft proofing is not available while no soft proofing ICC profile is loaded.\n"), "salmon");
		}
	} else {
		if (gui.icc.soft_proof) {
			cmsCloseProfile(gui.icc.soft_proof);
			gui.icc.soft_proof = NULL;
		}
		siril_log_color_message(_("Soft proofing ICC profile deactivated. Soft proofing is not available while no soft proofing ICC profile is loaded.\n"), "salmon");
	}
	g_mutex_unlock(&soft_proof_profile_mutex);
	refresh_icc_transforms();
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
	// Handle initial assignment of an ICC profile
	if (!gfit.color_managed || !gfit.icc_profile) {
		if (gfit.icc_profile) {
			cmsCloseProfile(gfit.icc_profile);
			gfit.icc_profile = NULL;
		}
		goto FINISH;
	}

	cmsUInt32Number gfit_colorspace = cmsGetColorSpace(gfit.icc_profile);
	cmsUInt32Number gfit_colorspace_channels = cmsChannelsOf(gfit_colorspace);
	cmsUInt32Number target_colorspace = cmsGetColorSpace(target);
	cmsUInt32Number target_colorspace_channels = cmsChannelsOf(target_colorspace);

	if (target_colorspace != cmsSigGrayData && target_colorspace != cmsSigRgbData) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Color space not supported"), _("Siril only supports representing the image in Gray or RGB color spaces at present. You cannot assign or convert to non-RGB color profiles"));
		return;
	}
	if (gfit_colorspace_channels != target_colorspace_channels) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Transform not supported"), _("Image cannot be assigned a color profile with a different number of channels to its current color profile"));
		return;
	}
	// We save the undo state if dealing with gfit
	undo_save_state(&gfit, _("Color profile assignment"));
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
	// We save the undo state if dealing with gfit

undo_save_state(&gfit, _("Color profile removal"));
	if (gfit.icc_profile)
		cmsCloseProfile(gfit.icc_profile);
	color_manage(&gfit, FALSE);
	gtk_widget_set_sensitive(lookup_widget("icc_convertto"), gfit.color_managed);
	gtk_widget_set_sensitive(lookup_widget("icc_remove"), gfit.color_managed);
	set_source_information();
	refresh_icc_transforms();
	notify_gfit_modified();

}

void on_icc_convertto_clicked(GtkButton* button, gpointer* user_data) {
	if (!gfit.color_managed || !gfit.icc_profile) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No color profile set"), _("The  current image has no color profile. You need to assign one first."));
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
	siril_colorspace_transform(&gfit, target);

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
	switch (target_index) {
		case NONE:
			if (target) {
				GtkLabel* mfr_label = (GtkLabel*) lookup_widget("icc_target_mfr_label");
				GtkLabel* copyright_label = (GtkLabel*) lookup_widget("icc_target_copyright_label");
				cmsCloseProfile(target);
				target = NULL;
				gtk_label_set_text(label, "");
				gtk_label_set_text(copyright_label, "");
				gtk_label_set_text(mfr_label, "");
			}
			return;
		case SRGB_LINEAR:
			if (target) {
				cmsCloseProfile(target);
			}
			target = srgb_linear();
			break;
		case SRGB_TRC:
			if (target) {
				cmsCloseProfile(target);
			}
			target = srgb_trc();
			break;
		case REC2020_LINEAR:
			if (target) {
				cmsCloseProfile(target);
			}
			target = rec2020_linear();
			break;
		case REC2020_TRC:
			if (target) {
				cmsCloseProfile(target);
			}
			target = rec2020_trc();
			break;
		case GRAY_LINEAR:
			if (target) {
				cmsCloseProfile(target);
			}
			target = gray_linear();
			break;
		case GRAY_SRGBTRC:
			if (target) {
				cmsCloseProfile(target);
			}
			target = gray_srgbtrc();
			break;
		case GRAY_REC709TRC:
			if (target) {
				cmsCloseProfile(target);
			}
			target = gray_rec709trc();
			break;
	}
	set_target_information();
}

void on_icc_target_filechooser_file_set(GtkFileChooser* filechooser, gpointer* user_data) {
	GtkComboBox* target_combo = (GtkComboBox*) lookup_widget("icc_target_combo");
	gtk_combo_box_set_active(target_combo, 0);
	gchar *filename = gtk_file_chooser_get_filename(filechooser);
	if (filename) {
		target = cmsOpenProfileFromFile(filename, "r");
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not load file"), _("Error: could not load selected file, or it does not contain a valid ICC profile."));
		g_free(filename);
		gtk_file_chooser_unselect_all(filechooser);
		return;
	}
	if (!target) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not load file"), _("Error: could not load selected file, or it does not contain a valid ICC profile."));
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

void on_icc_main_window_button_clicked(GtkButton *button, gpointer user_data) {
	siril_open_dialog("icc_dialog");
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

void on_icc_gamut_visualisation_clicked() {
	GtkWidget *win = lookup_widget("icc_gamut_dialog");
	gtk_window_set_transient_for(GTK_WINDOW(win), GTK_WINDOW(lookup_widget("settings_window")));
	/* Here this is wanted that we do not use siril_open_dialog */
	gtk_widget_show(win);
}

void on_icc_gamut_close_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *win = lookup_widget("icc_gamut_dialog");
	gtk_widget_hide(win);
}
