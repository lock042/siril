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

#include <glib.h>
#include "algos/lcms_acceleration/lcms2_fast_float.h"
#include "algos/lcms_acceleration/lcms2_threaded.h"
#include "core/siril.h"
#include "algos/colors.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "icc_default_profiles.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/icc_profile.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/registration_preview.h"
#include "gui/utils.h"
#include "gui/siril_plot.h"
#include "gui/siril_preview.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "io/siril_plot.h"
#include "core/siril_log.h"

// For the log message about JPEG ICC profile support at startup
#ifdef HAVE_LIBJPEG
#include <jconfig.h>
#endif

static GMutex monitor_profile_mutex;
static GMutex soft_proof_profile_mutex;
static GMutex default_profiles_mutex;
static GMutex display_transform_mutex;

static gboolean profile_check_verbose = TRUE;

////// Functions //////
struct cm_struct {
	fits *fit;
	gboolean active;
};

static gboolean cm_worker(gpointer user_data) {
	struct cm_struct *data = (struct cm_struct*) user_data;
	fits *fit = data->fit;
	gboolean active = data->active;
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
	return FALSE;
}

void color_manage(fits *fit, gboolean active) {
	fit->color_managed = active;
	struct cm_struct data = { fit, active };
	if (fit == &gfit && !com.script) {
		cm_worker(&data);
	}
}

static gchar *siril_color_profile_get_info (cmsHPROFILE profile, cmsInfoType info) {
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

gchar *
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

gchar *
	siril_color_profile_get_manufacturer (cmsHPROFILE profile) {
	gchar *string = siril_color_profile_get_info (profile, cmsInfoManufacturer);
	return string;
}

gchar *
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
cmsHPROFILE srgb_trcv2() {
	return cmsOpenProfileFromMem(sRGB_elle_V2_srgbtrc_icc, sRGB_elle_V2_srgbtrc_icc_len);
}

cmsHPROFILE rec2020_linear() {
	return cmsOpenProfileFromMem(Rec2020_elle_V4_g10_icc, Rec2020_elle_V4_g10_icc_len);
}
cmsHPROFILE rec2020_trc() {
	return cmsOpenProfileFromMem(Rec2020_elle_V4_rec709_icc, Rec2020_elle_V4_rec709_icc_len);
}
cmsHPROFILE rec2020_trcv2() {
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
cmsHPROFILE gray_srgbtrcv2() {
	return cmsOpenProfileFromMem(Gray_elle_V2_srgbtrc_icc, Gray_elle_V2_srgbtrc_icc_len);
}
cmsHPROFILE gray_rec709trcv2() {
	return cmsOpenProfileFromMem(Gray_elle_V2_rec709_icc, Gray_elle_V2_rec709_icc_len);
}

cmsHPROFILE srgb_monitor_perceptual() {
	return cmsOpenProfileFromMem(sRGB_v4_ICC_preference_icc, sRGB_v4_ICC_preference_icc_len);
}

void export_profile(cmsHPROFILE profile, const char *provided_filename) {
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

//Two functions to check if profiles are RGB or Gray
gboolean siril_color_profile_is_rgb(cmsHPROFILE profile) {
	return (cmsGetColorSpace (profile) == cmsSigRgbData);
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
//	if (gfit.color_managed) {
		if (gui.icc.proofing_transform) {
			cmsDeleteTransform(gui.icc.proofing_transform);
			gui.icc.proofing_transform = NULL;
		}
//	}
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
		g_mutex_lock(&monitor_profile_mutex);
		if (gui.icc.monitor)
			cmsCloseProfile(gui.icc.monitor);
		gui.icc.monitor = com.pref.icc.rendering_intent == INTENT_PERCEPTUAL ? srgb_monitor_perceptual() : srgb_trc();
		g_mutex_unlock(&monitor_profile_mutex);
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
				siril_log_color_message(_("Error opening matched grayscale working profile. Profile set to Gray with sRGB tone response curve.\n"), "red");
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

	// Target profiles for embedding in saved files
	com.icc.srgb_out = srgb_trcv2();
	com.icc.mono_out = gray_srgbtrcv2();

	validate_custom_profiles();

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

void cleanup_common_profiles() {
	if (com.icc.mono_linear)
		cmsCloseProfile(com.icc.mono_linear);
	if (com.icc.working_standard)
		cmsCloseProfile(com.icc.working_standard);
	if (com.icc.mono_standard)
		cmsCloseProfile(com.icc.mono_standard);
	if (com.icc.srgb_profile)
		cmsCloseProfile(com.icc.srgb_profile);
	if (com.icc.srgb_out)
		cmsCloseProfile(com.icc.srgb_out);
	if (com.icc.working_out)
		cmsCloseProfile(com.icc.working_out);
	if (com.icc.mono_out)
		cmsCloseProfile(com.icc.mono_out);
	if (gui.icc.monitor)
		cmsCloseProfile(gui.icc.monitor);
	if (gui.icc.soft_proof)
		cmsCloseProfile(gui.icc.soft_proof);
	if (gui.icc.proofing_transform)
		cmsDeleteTransform(gui.icc.proofing_transform);
	memset(&gui.icc, 0, sizeof(struct gui_icc));
	if (com.icc.context_single)
		cmsDeleteContext(com.icc.context_single);
	if (com.icc.context_threaded)
		cmsDeleteContext(com.icc.context_threaded);
	memset(&com.icc, 0, sizeof(struct common_icc));
	siril_debug_print("ICC profiles cleaned up\n");
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
	gui_function(redraw_previews, NULL);
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
			siril_log_color_message(_("Error: no filename specified for custom monitor profile.\n"), "red");
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
			siril_log_color_message(_("Error: no filename specified for output device proofing profile.\n"), "red");
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
	if (!profile)
		return NULL;
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

static gboolean siril_color_profile_get_rgb_matrix_colorants_d50 (cmsHPROFILE profile, cmsCIEXYZTRIPLE *XYZtriple, cmsCIEXYZ *whitepoint) {
	double redrgb[3] = { 1.0, 0.0, 0.0 };
	double greenrgb[3] = { 0.0, 1.0, 0.0 };
	double bluergb[3] = { 0.0, 0.0, 1.0 };
	double whitergb[3] = { 1.0, 1.0, 1.0 };

	// Create D50 XYZ profile (ICC PCS standard)
	cmsHPROFILE profileXYZ = cmsCreateXYZProfile();

	// Create transform WITH chromatic adaptation enabled (default behavior)
	// This ensures colorants are adapted to D50
	cmsHTRANSFORM transform = cmsCreateTransformTHR(com.icc.context_single,
	                                               profile, TYPE_RGB_DBL,
	                                               profileXYZ, TYPE_XYZ_DBL,
	                                               INTENT_RELATIVE_COLORIMETRIC,
	                                               cmsFLAGS_NOCACHE);

	cmsCloseProfile(profileXYZ);

	if (!transform) {
		return FALSE;
	}

	// Transform the primaries and white point
	// These will be automatically adapted to D50 by LCMS2
	cmsDoTransform(transform, &redrgb, &XYZtriple->Red, 1);
	cmsDoTransform(transform, &greenrgb, &XYZtriple->Green, 1);
	cmsDoTransform(transform, &bluergb, &XYZtriple->Blue, 1);
	cmsDoTransform(transform, &whitergb, whitepoint, 1);

	cmsDeleteTransform(transform);

	// Set the whitepoint to D50 since that's what the colorants are adapted to
	whitepoint->X = 0.9642;
	whitepoint->Y = 1.0000;
	whitepoint->Z = 0.8249;

	return TRUE;
}

static void siril_color_profile_set_tag (cmsHPROFILE profile,
								cmsTagSignature sig,
								const gchar *tag) {
	cmsMLU *mlu;

	mlu = cmsMLUalloc (NULL, 1);
	cmsMLUsetASCII (mlu, "en", "US", tag);
	cmsWriteTag (profile, sig, mlu);
	cmsMLUfree (mlu);
}

static void siril_color_profile_make_tag (cmsHPROFILE profile,
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
		/* otherwise prefix the existing tag with a siril prefix
		* indicating that the profile has been generated
		*/

		if (g_str_has_prefix (original_tag, siril_prefix)) {
			/* don't add multiple siril prefixes */
			siril_color_profile_set_tag (profile, sig, original_tag);
		}else if (siril_prefix_alt &&
				g_str_has_prefix (original_tag, siril_prefix_alt)) {
			/* replace siril prefix_alt by prefix */
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
		if (! siril_color_profile_get_rgb_matrix_colorants_d50 (profile, &XYZtriple, &whitepoint))
			return NULL;
	} else if (siril_color_profile_is_gray (profile)) {
		// For grayscale profiles, use D50 whitepoint
		whitepoint.X = 0.9642;
		whitepoint.Y = 1.0000;
		whitepoint.Z = 0.8249;
	} else {
		return NULL;
	}

	target_profile = cmsCreateProfilePlaceholder (0);
	if (!target_profile)
		return NULL;

	cmsSetProfileVersion (target_profile, 4.3);
	cmsSetDeviceClass (target_profile, cmsSigDisplayClass);
	cmsSetPCS (target_profile, cmsSigXYZData);

	// Use D50 whitepoint - this matches the adapted colorants
	cmsWriteTag (target_profile, cmsSigMediaWhitePointTag, &whitepoint);

	curve = cmsBuildGamma (NULL, 1.00);
	if (!curve) {
		cmsCloseProfile(target_profile);
		return NULL;
	}

	siril_color_profile_make_tag (target_profile, cmsSigProfileDescriptionTag,
									"linear TRC from unnamed profile",
									"linear TRC from ",
									"linear TRC from ",
									siril_color_profile_get_description (profile));

	if (siril_color_profile_is_rgb (profile)) {
		cmsSetColorSpace (target_profile, cmsSigRgbData);

		// Use the D50-adapted colorants
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
								"Generated by Siril",
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
			if (profile_check_verbose)
				siril_log_message(_("FITS did not contain an ICC profile but is declared to be stretched. Assigning a sRGB color profile.\n"));
			// sRGB because this is the implicit assumption made in older versions
			fit->icc_profile = fit->naxes[2] == 1 ? gray_srgbtrc() : srgb_trc();
			// Clear the hint
			com.icc.srgb_hint = FALSE;
		} else if (fit_appears_stretched(fit)) {
			if (profile_check_verbose)
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

void enable_profile_check_verbose() {
	profile_check_verbose = TRUE;
}

void disable_profile_check_verbose() {
	profile_check_verbose = FALSE;
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
		if (fit == &gfit && !com.headless) {
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

// This function overrides any GTK theme to assign a white border and D50 Gray
// background to the 4 vports to approximate ISO 12646 viewing conditions.
// It is recommended in conjunction with the excellent Equilux GTK theme.

void siril_plot_colorspace(cmsHPROFILE profile, gboolean compare_srgb) {
	siril_plot_data *spl_data = init_siril_plot_data();
	if (!spl_data) {
		return;
	}
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
		free_siril_plot_data(spl_data);
		return;
	}
	if (! siril_color_profile_get_rgb_matrix_colorants_d50 (profile, &XYZtriple, &whitepoint)) {
		siril_log_message(_("Error reading chromaticities\n"));
		free_siril_plot_data(spl_data);
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

	gchar *title1 = g_strdup_printf(_("Source Color Profile Chromaticity Diagram\n"
					"<span size=\"small\">"
					"%s"
					"</span>"), description);
	free(description);
	siril_plot_set_xlabel(spl_data, _("CIE x"));
	siril_plot_set_savename(spl_data, "color_profile");
	siril_plot_set_title(spl_data, title1);
	g_free(title1);
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

	siril_add_pythonsafe_idle(create_new_siril_plot_window, spl_data);
	siril_add_idle(end_generic, NULL);
}
