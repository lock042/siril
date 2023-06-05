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
#include "algos/lcms_fast_float/lcms2_fast_float.h"
#include "core/siril.h"
#include "core/OS_utils.h"
#include "icc_profile.h"
#include "icc_default_profiles.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "io/single_image.h"
#include "core/siril_log.h"
#include "core/siril_app_dirs.h"
#include "core/proto.h"

#define DEFAULT_PATH_GV2g22 "Gray-siril-V2-g22.icc"
#define DEFAULT_PATH_GV4g10 "Gray-siril-V4-g10.icc"
#define DEFAULT_PATH_GV4g22 "Gray-siril-V4-g22.icc"
#define DEFAULT_PATH_RGBV2g22 "sRGB-siril-V2-g22.icc"
#define DEFAULT_PATH_RGBV4g10 "sRGB-siril-V4-g10.icc"
#define DEFAULT_PATH_RGBV4g22 "sRGB-siril-V4-g22.icc"

#define INDEX_GV2G22 0
#define INDEX_GV4G10 1
#define INDEX_GV4G22 2
#define INDEX_SRGBV2G22 3
#define INDEX_SRGBV4G10 4
#define INDEX_SRGBV4G22 5

#define MAX_SYSTEM_ICC 5
// The following two indices refer to user-provided screen and proof profiles
// that may be stored in com.pref.icc but for which there are no defaults
#define INDEX_CUSTOM_MONITOR 6
#define INDEX_CUSTOM_PROOF 7
#define INDEX_CUSTOM_EXPORT 8

cmsHPROFILE copyICCProfile(cmsHPROFILE profile);
cmsHTRANSFORM sirilCreateTransform(cmsHPROFILE Input, cmsUInt32Number InputFormat, cmsHPROFILE Output, cmsUInt32Number OutputFormat, cmsUInt32Number Intent, cmsUInt32Number dwFlags);

const char* default_icc_paths[] = { DEFAULT_PATH_GV2g22, DEFAULT_PATH_GV4g10, DEFAULT_PATH_GV4g22, DEFAULT_PATH_RGBV2g22, DEFAULT_PATH_RGBV4g10, DEFAULT_PATH_RGBV4g22 };
static GMutex monitor_profile_mutex;
static GMutex soft_proof_profile_mutex;
static GMutex default_profiles_mutex;

static cmsHPROFILE target = NULL; // Target profile for the GUI tool

////// Functions //////

void initialize_icc_profiles_paths() {
	int nb_icc = sizeof(default_icc_paths) / sizeof(const char *);
	int maxpath = get_pathmax();
	for (int icc = 0; icc < nb_icc; icc++) {
		if (com.pref.icc.icc_paths[icc] &&
				com.pref.icc.icc_paths[icc][0] != '\0')
			continue;
		char path[maxpath];
		gchar *filename = g_build_filename(siril_get_system_data_dir(), "icc", default_icc_paths[icc], NULL);
		strncpy(path, filename, maxpath - 1);
		com.pref.icc.icc_paths[icc] = g_strdup(path);
		g_free(filename);
	}
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

void initialize_profiles_and_transforms() {
	// Enable the fast float plugin (as long as the OS / lcms2 version blacklist isn't triggered)
#ifndef EXCLUDE_FF
	cmsPlugin(cmsFastFloatExtensions());
#endif
	// Initialize paths to standard ICC profiles
	initialize_icc_profiles_paths();

	// Set alarm codes for soft proof out-of-gamut warning
	cmsUInt16Number alarmcodes[16] = { 65535, 0, 65535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	cmsSetAlarmCodes(alarmcodes); // Out of gamut colours will be shown in magenta
	int error = 0;

	// Intents
	com.icc.save_intent = com.pref.icc.export_intent;
	gui.icc.rendering_intent = com.pref.icc.rendering_intent;
	gui.icc.proofing_intent = com.pref.icc.proofing_intent;

	// Linear working profiles
	com.icc.working_linear = srgb_linear();
	com.icc.mono_linear = gray_linear();

	// Native TRC working profiles
	com.icc.working_standard = srgb_trc();
	com.icc.mono_standard = gray_srgbtrc();

	// Target profiles for embedding in saved files
	com.icc.srgb_out = srgb_trcv2();
	com.icc.working_out = srgb_trcv2();
	com.icc.mono_out = gray_srgbtrcv2();

	// ICC availability
	com.icc.available = (com.icc.working_linear && com.icc.mono_linear && com.icc.working_standard && com.icc.mono_standard && com.icc.working_out && com.icc.mono_out);
	gui.icc.available = (com.icc.available); // && gui.icc_profile_rgb && gui.icc_profile_mono);
	if ((com.headless && !com.icc.available) || (!com.headless && !gui.icc.available)) {
		error++;
		siril_log_message(_("Error: standard color management profiles could not be loaded. Siril will continue "
							"but colors in saved images will not be consistent when viewed in other applications "
							"or printed. Check your installation is correct.\n"));
	} else { // No point loading monitor and soft proof profiles if the standard ones are unavailable
		// Open the custom monitor and soft proofing profiles if there is a path set in preferences
		g_mutex_lock(&monitor_profile_mutex);
		if (com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR] && com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR][0] != '\0') {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR], "r");
			if (!gui.icc.monitor) {
				error++;
				siril_log_message(_("Warning: custom monitor profile set but could not be loaded. Display will use a "
									"sRGB profile with the standard sRGB TRC.\n"));
			} else {
				siril_log_message(_("Monitor ICC profile loaded from %s\n"), com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR]);
			}
		}
		if (!gui.icc.monitor) {
			gui.icc.monitor = srgb_trc();
			siril_log_message(_("Monitor ICC profile set to sRGB (standard sRGB TRC)\n"));
		}
		g_mutex_unlock(&monitor_profile_mutex);

		if (com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
			g_mutex_lock(&soft_proof_profile_mutex);
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF], "r");
			g_mutex_unlock(&soft_proof_profile_mutex);
		if (!gui.icc.soft_proof) {
				error++;
				siril_log_message(_("Warning: soft proofing profile set but could not be loaded. Soft proofing will be "
									"unavailable.\n"));
			} else {
				siril_log_message(_("Soft proofing ICC profile loaded from %s\n"), com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF]);
			}
		} else {
			siril_log_message(_("No soft proofing ICC profile set. Soft proofing is unavailable.\n"));
		}

	}
	if (!error) {
		siril_log_message(_("ICC profiles loaded correctly. Workflow will be color managed.\n"));
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

// Loads a custom monitor profile from a path in com.pref.icc.icc_paths
// The path must be set by the user in preferences, there is no default custom monitor profile
int load_monitor_icc_profile(const char* filename) {
	if (com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR] && com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR][0] != '\0') {
	g_mutex_lock(&monitor_profile_mutex);
		if (gui.icc.monitor)
			cmsCloseProfile(gui.icc.monitor);
		gui.icc.monitor = cmsOpenProfileFromFile(filename, "r");
	g_mutex_unlock(&soft_proof_profile_mutex);
		if (!gui.icc.monitor)
			return 1;
		else
			return 0;
	} else return 1;
}

// Loads a custom proof profile from a path in com.pref.icc.icc_paths
// The path must be set by the user in preferences, there is no default custom proof profile
int load_soft_proof_icc_profile(const char* filename) {
	if (com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
		g_mutex_lock(&soft_proof_profile_mutex);
		gui.icc.soft_proof = cmsOpenProfileFromFile(filename, "r");
		g_mutex_unlock(&soft_proof_profile_mutex);
		if (!gui.icc.soft_proof)
			return 1;
		else
			return 0;
	} else return 1;
}

void assign_linear_icc_profile(fits *fit) {
	if (fit->icc_profile) {
		cmsCloseProfile(fit->icc_profile);
	}
	fit->icc_profile = copyICCProfile(gfit.naxes[2] == 1 ? com.icc.mono_linear : com.icc.working_linear);
}

// Even for mono images with a Gray profile, the display datatype is always RGB;
// lcms2 takes care of populating all 3 output channels for us.
cmsHTRANSFORM initialize_display_transform() {
	g_assert(gui.icc.monitor);
	cmsHTRANSFORM transform = NULL;
	if (gfit.icc_profile == NULL)
		return NULL;
	cmsUInt32Number gfit_signature = cmsGetColorSpace(gfit.icc_profile);
	cmsUInt32Number srctype = get_planar_formatter_type(gfit_signature, gfit.type, TRUE);
	g_mutex_lock(&monitor_profile_mutex);
	transform = cmsCreateTransform(gfit.icc_profile, srctype, gui.icc.monitor, TYPE_RGB_16_PLANAR, gui.icc.rendering_intent, 0);
	g_mutex_unlock(&monitor_profile_mutex);
	if (transform == NULL)
		siril_log_message("Error: failed to create display_transform!\n");
	return transform;
}

cmsHTRANSFORM initialize_export8_transform(fits* fit) {
	g_assert(com.icc.working_standard);
	cmsHTRANSFORM transform = NULL;
	if (fit->icc_profile == NULL)
		return NULL;
	cmsUInt32Number fit_signature = cmsGetColorSpace(fit->icc_profile);
	cmsUInt32Number srctype = get_planar_formatter_type(fit_signature, fit->type, TRUE);
	cmsUInt32Number desttype = (fit->naxes[2] == 1 ? TYPE_GRAY_16 : TYPE_RGB_16_PLANAR);
	transform = sirilCreateTransform(fit->icc_profile, srctype, (fit->naxes[2] == 3 ? com.icc.working_standard : com.icc.mono_standard), desttype, gui.icc.rendering_intent, 0);
	if (transform == NULL)
		siril_log_message("Error: failed to create export colorspace transform!\n");
	return transform;
}

cmsHTRANSFORM initialize_proofing_transform() {
	g_assert(gui.icc.monitor);
	if (gfit.icc_profile == NULL || gui.icc.soft_proof == NULL)
		return NULL;
	cmsUInt32Number flags = cmsFLAGS_SOFTPROOFING;
	if (gui.cut_over) {
		flags |= cmsFLAGS_GAMUTCHECK;
	}
	cmsUInt32Number type = (gfit.naxes[2] == 1 ? TYPE_GRAY_16 : TYPE_RGB_16_PLANAR);
	g_mutex_lock(&soft_proof_profile_mutex);
	g_mutex_lock(&monitor_profile_mutex);
	cmsHPROFILE proofing_transform = cmsCreateProofingTransform(
						gfit.icc_profile,
						type,
						gui.icc.monitor,
						type,
						gui.icc.soft_proof,
						gui.icc.rendering_intent,
						gui.icc.proofing_intent,
						flags);
	g_mutex_unlock(&monitor_profile_mutex);
	g_mutex_unlock(&soft_proof_profile_mutex);
	return proofing_transform;
}

void refresh_icc_transforms() {
	if (gui.icc.available) {
		if (gui.icc.display_transform != NULL)
			cmsDeleteTransform(gui.icc.display_transform);
		gui.icc.display_transform = initialize_display_transform();

		if (gui.icc.proofing_transform != NULL)
			cmsDeleteTransform(gui.icc.proofing_transform);
		gui.icc.proofing_transform = initialize_proofing_transform();
	}
}

// Used where the max value of data in fits <=UCHAR_MAX
// (assume we are dealing with an 8 bit image)
BYTE uchar_pixel_icc_tx(BYTE in, int channel, int nchans, cmsHTRANSFORM transform) {
	g_assert(channel >= 0 && channel < 3);
	g_assert(channel < nchans);
	g_assert(transform);
	BYTE input[3] = { 0 };
	BYTE output[3] = { 0 };
	BYTE out;
	input[channel] = in;
	if (nchans > 2) {
		cmsDoTransform(transform, input, output, 1);
	} else {
		cmsDoTransform(transform, &in, &out, 1);
		output[0] = out;
		}
	return output[channel];
}
// Used for 16-bit images
WORD ushrt_pixel_icc_tx(WORD in, int channel, int nchans, cmsHTRANSFORM transform) {
	g_assert(channel >= 0 && channel < 3);
	g_assert(channel < nchans);
	g_assert(transform);
	WORD input[3] = { 0 };
	WORD output[3] = { 0 };
	WORD out;
	input[channel] = in;
	if (nchans > 2) {
		cmsDoTransform(transform, input, output, 1);
	} else {
		cmsDoTransform(transform, &in, &out, 1);
		output[0] = out;
		}
	return output[channel];
}

// Used for float images
float float_pixel_icc_tx(float in, int channel, int nchans, cmsHTRANSFORM transform) {
	g_assert(channel >= 0 && channel < 3);
	g_assert(channel < nchans);
	g_assert(transform);
	float input[3] = { 0 };
	float output[3] = { 0 };
	float out;
	input[channel] = in;
	if (nchans > 2) {
		cmsDoTransform(transform, input, output, 1);
	} else {
		cmsDoTransform(transform, &in, &out, 1);
		output[0] = out;
		}
	return output[channel];
}

// These functions populate the ICC profile buffers used in file save routines
unsigned char* get_sRGB_profile_data(guint32 *len, gboolean linear) {
	unsigned char* block = NULL;
	cmsUInt32Number length;
	cmsHPROFILE *profile = linear ? &com.icc.working_linear : &com.icc.working_out;
	cmsBool ret = cmsSaveProfileToMem(*profile, NULL, &length);
	if (length > 0) {
		block = malloc(length * sizeof(BYTE));
		ret = cmsSaveProfileToMem(*profile, (void*) block, &length);
	}
	if (!ret) {
		siril_debug_print("Error preparing ICC profile for embedding. This may happen if ICC profiles were not correctly loaded at startup.\n");
		return NULL;
	}
	*len = length;
	return block;
}

unsigned char* get_gray_profile_data(guint32 *len, gboolean linear) {
	unsigned char* block = NULL;
	cmsUInt32Number length;
	cmsHPROFILE *profile = linear ? &com.icc.mono_linear : &com.icc.mono_out;
	cmsBool ret = cmsSaveProfileToMem(*profile, NULL, &length);
	if (length > 0) {
		block = malloc(length * sizeof(BYTE));
		ret = cmsSaveProfileToMem(*profile, (void*) block, &length);
	}
	if (!ret) {
		siril_debug_print("Error preparing ICC profile for embedding. This may happen if ICC profiles were not correctly loaded at startup.\n");
		return NULL;
	}
	*len = length;
	return block;
}

unsigned char* get_icc_profile_data(cmsHPROFILE *profile, guint32 *len) {
	unsigned char* block = NULL;
	cmsUInt32Number length;
	cmsBool ret = cmsSaveProfileToMem(*profile, NULL, &length);
	if (length > 0) {
		block = malloc(length * sizeof(BYTE));
		ret = cmsSaveProfileToMem(*profile, (void*) block, &length);
	}
	if (!ret) {
		siril_debug_print("Error preparing ICC profile for embedding...\n");
		return NULL;
	}
	*len = length;
	return block;
}

cmsHPROFILE copyICCProfile(cmsHPROFILE profile) {
	cmsUInt32Number length = 0;
	cmsUInt8Number* block = NULL;
	cmsBool ret = FALSE;
	cmsHPROFILE retval = NULL;
	if (profile) {
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

// This function is for initializing the profile during file import, hence the fallback assumption of
// a sRGB profile
void fits_initialize_icc(fits *fit, cmsUInt8Number* EmbedBuffer, cmsUInt32Number EmbedLen) {
	if (com.icc.available) {
		if (EmbedBuffer) {
			// If there is an embedded profile we will use it
			fit->icc_profile = cmsOpenProfileFromMem(EmbedBuffer, EmbedLen);
		} else {
			// If there is no embedded profile we assume the usual sRGB D65 g22
			fit->icc_profile = copyICCProfile((fit->naxes[2] == 1) ? com.icc.mono_out : com.icc.srgb_out);
		}
	}
}

// This function is for checking that a fits has an associated ICC profile during image ops, hence the
// fallback assumption of a linear profile
void fits_check_icc(fits *fit) {
	if (!fit->icc_profile) {
		// If there is no embedded profile we assume the usual sRGB D65 g22
		fit->icc_profile = copyICCProfile((fit->naxes[2] == 1) ? com.icc.mono_linear : com.icc.working_linear);
	}
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
			ret_b = cmsSaveProfileToMem(a, (void*) block_b, &length_b);
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

///// Preferences callbacks

// Being able to alter the monitor and soft_proof profiles and intents from the GTK thread means all operations
// that use these profiles need to go inside mutexes to prevent the profiles being ripped out from under them.
// This is not required for profiles in FITS structures as these are not subject to thread contention issues in
// the same way.

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
	if (com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR] && com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR][0] != '\0') {
		g_free(com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR]);
		com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR] = NULL;
		cmsCloseProfile(gui.icc.monitor);
		if (com.icc.available) {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc.icc_paths[INDEX_SRGBV4G22], "r");
			if (gui.icc.monitor) {
				siril_log_message(_("Monitor ICC profile set to sRGB (D65 whitepoint, gamma = 2.2)\n"));
			} else {
				siril_log_color_message(_("Error: standard sRGB ICC profile could not be loaded.\n"), "red");
				com.icc.available = FALSE;
			}
		} else gui.icc.monitor = NULL;
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
	if (com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
		g_free(com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF]);
		com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF] = NULL;
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
		if (!com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR] || com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR][0] == '\0') {
			com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR] = g_strdup(gtk_file_chooser_get_filename(filechooser));
		}
		if (!com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR] || com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR][0] == '\0') {
			siril_log_color_message(_("Error: no filename specfied for custom monitor profile.\n"), "red");
			no_file = TRUE;
		} else {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR], "r");
		}
		if (gui.icc.monitor) {
			siril_log_message(_("Monitor profile loaded from %s\n"), com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR], "r");
			g_mutex_unlock(&monitor_profile_mutex);
			refresh_icc_transforms();
			return;
		} else {
			if (!no_file) {
				siril_log_color_message(_("Monitor profile could not be loaded from %s\n"), "red", com.pref.icc.icc_paths[INDEX_CUSTOM_MONITOR]);
			}
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc.icc_paths[INDEX_SRGBV4G22], "r");
			if (gui.icc.monitor) {
				siril_log_message(_("Monitor ICC profile set to sRGB (D65 whitepoint, gamma = 2.2)\n"));
			} else {
				siril_log_color_message(_("Error: standard sRGB ICC profile could not be loaded. Color management is disabled.\n"), "red");
				com.icc.available = FALSE;
			}
		}
	} else {
		gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc.icc_paths[INDEX_SRGBV4G22], "r");
		if (gui.icc.monitor) {
			siril_log_message(_("Monitor ICC profile set to sRGB (D65 whitepoint, gamma = 2.2)\n"));
		} else {
			siril_log_color_message(_("Error: standard sRGB ICC profile could not be loaded. Color management is disabled.\n"), "red");
			com.icc.available = FALSE;
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
		if (!com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF] || com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF][0] == '\0') {
			com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF] = g_strdup(gtk_file_chooser_get_filename(filechooser));
		}
		if (!com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF] || com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF][0] == '\0') {
			siril_log_color_message(_("Error: no filename specfied for custom proofing profile.\n"), "red");
			no_file = TRUE;
		} else {
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF], "r");
		}
		if (gui.icc.soft_proof) {
			siril_log_message(_("Soft proofing profile loaded from %s\n"), com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF]);
			g_mutex_unlock(&soft_proof_profile_mutex);
			refresh_icc_transforms();
			return;
		} else {
			if (!no_file) {
				siril_log_color_message(_("Soft proofing profile could not be loaded from %s\n"), "red", com.pref.icc.icc_paths[INDEX_CUSTOM_PROOF]);
			}
			siril_log_color_message(_("Soft proofing is not available while no soft proofing ICC profile is loaded.\n"), "salmon");
		}
	} else {
		gui.icc.soft_proof = NULL;
		siril_log_color_message(_("Soft proofing ICC profile closed. Soft proofing is not available while no soft proofing ICC profile is loaded.\n"), "salmon");
	}
	g_mutex_unlock(&soft_proof_profile_mutex);
	refresh_icc_transforms();
}

/* This function is used for compensating for differences
 * between themonitor profile and working profile XYZ primaries,
 * white point and colorants following an auto histogram
 * stretch. It is unlikely to be useful for anything else. */

cmsHPROFILE adjust_primaries (cmsHPROFILE working, cmsHPROFILE disp) {
	/* Copy the display profile to the primary adjustment temporary profile
	 * (this ensures the TRC is the same) */
	cmsHPROFILE temp = copyICCProfile(disp);

	/* Get and set the black and white points */
	cmsCIEXYZ  *white = cmsReadTag(working, cmsSigMediaWhitePointTag);
	if (white) {
		cmsCIEXYZ white_tag = *white;
		cmsWriteTag (temp, cmsSigMediaWhitePointTag, &white_tag);
	}
	cmsCIEXYZ  *black = cmsReadTag(working, cmsSigMediaBlackPointTag);
	if (black) {
		cmsCIEXYZ black_tag = *black;
		cmsWriteTag (temp, cmsSigMediaBlackPointTag, &black_tag);
	}

	/* Get and set the colorants */
	cmsCIEXYZ  *red = cmsReadTag(working, cmsSigRedColorantTag);
	if (red) {
		cmsCIEXYZ red_tag = *red;
		cmsWriteTag (temp, cmsSigRedColorantTag, &red_tag);
	}

	cmsCIEXYZ  *green = cmsReadTag(working, cmsSigGreenColorantTag);
	if (green) {
		cmsCIEXYZ green_tag = *green;
		cmsWriteTag (temp, cmsSigGreenColorantTag, &green_tag);
	}

	cmsCIEXYZ  *blue = cmsReadTag(working, cmsSigBlueColorantTag);
	if (blue) {
		cmsCIEXYZ blue_tag = *blue;
		cmsWriteTag (temp, cmsSigBlueColorantTag, &blue_tag);
	}

	return temp;
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

void error_loading_profile() {
	siril_message_dialog(GTK_MESSAGE_ERROR, _("Error loading profile"),
						 _("The selected profile could not be loaded or did not contain a valid ICC profile. Defaulting to sRGB."));
	if (com.icc.working_linear)
		cmsCloseProfile(com.icc.working_linear);
	com.icc.working_linear = srgb_linear();
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
}

void update_profiles_after_gamut_change() {
	working_gamut_type working_gamut = com.pref.icc.working_gamut;
	g_mutex_lock(&default_profiles_mutex);
	switch (working_gamut) {
		case TYPE_SRGB:
			if (com.icc.working_linear)
				cmsCloseProfile(com.icc.working_linear);
			com.icc.working_linear = srgb_linear();
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
			if (com.icc.working_linear)
				cmsCloseProfile(com.icc.working_linear);
			com.icc.working_linear = rec2020_linear();
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
			if (com.icc.working_linear)
				cmsCloseProfile(com.icc.working_linear);
			if (!(com.pref.icc.custom_icc_linear && (com.icc.working_linear = cmsOpenProfileFromFile(com.pref.icc.custom_icc_linear, "r")))) {
				error_loading_profile();
			}
			break;
			// Custom profiles will also be used for the output profile
			if (com.icc.working_standard)
				cmsCloseProfile(com.icc.working_standard);
			if (!(com.pref.icc.custom_icc_trc && (com.icc.working_standard =
						cmsOpenProfileFromFile(com.pref.icc.custom_icc_trc, "r")))) {
				error_loading_profile();
			} else {
				if (com.icc.working_out)
					cmsCloseProfile(com.icc.working_out);
				com.icc.working_out = copyICCProfile(com.icc.working_standard);
			}
			break;
			// Custom profiles will also be used for the output profile
			if (com.icc.mono_standard)
				cmsCloseProfile(com.icc.mono_standard);
			if (!(com.pref.icc.custom_icc_gray && (com.icc.working_standard =
						cmsOpenProfileFromFile(com.pref.icc.custom_icc_gray, "r")))) {
				error_loading_profile();
			} else {
				if (com.icc.mono_out)
					cmsCloseProfile(com.icc.mono_out);
				com.icc.mono_out = copyICCProfile(com.icc.mono_standard);
			}
			break;
	}
	g_mutex_unlock(&default_profiles_mutex);
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

cmsHTRANSFORM sirilCreateTransform(cmsHPROFILE Input, cmsUInt32Number InputFormat, cmsHPROFILE Output, cmsUInt32Number OutputFormat, cmsUInt32Number Intent, cmsUInt32Number dwFlags) {
	cmsHTRANSFORM transform = NULL;
	g_mutex_lock(&default_profiles_mutex);
	transform = cmsCreateTransform(Input, InputFormat, Output, OutputFormat, Intent, dwFlags);
	g_mutex_unlock(&default_profiles_mutex);
	return transform;
}

//////// GUI callbacks for the color management dialog
void set_source_information() {
	if (!gfit.icc_profile) {
		siril_debug_print("Target profile is NULL\n");
		return;
	}
	// Set description
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_current_profile_label");
	GtkLabel* mfr_label = (GtkLabel*) lookup_widget("icc_mfr_label");
	GtkLabel* copyright_label = (GtkLabel*) lookup_widget("icc_copyright_label");
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
void set_target_information() {
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

void on_icc_cancel_clicked(GtkButton* button, gpointer* user_data) {
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_target_profile_label");
	GtkLabel* label2 = (GtkLabel*) lookup_widget("icc_target_mfr_label");
	GtkLabel* label3 = (GtkLabel*) lookup_widget("icc_target_copyright_label");
	GtkLabel* label4 = (GtkLabel*) lookup_widget("icc_current_profile_label");
	GtkLabel* label5 = (GtkLabel*) lookup_widget("icc_mfr_label");
	GtkLabel* label6 = (GtkLabel*) lookup_widget("icc_copyright_label");
	GtkComboBox* target_combo = (GtkComboBox*) lookup_widget("icc_target_combo");
	GtkComboBox* operation = (GtkComboBox*) lookup_widget("icc_operation_combo");
	GtkFileChooser* filechooser = (GtkFileChooser*) lookup_widget("icc_target_filechooser");
	gtk_combo_box_set_active(target_combo, 0);
	gtk_combo_box_set_active(operation, 0);
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
	if (gfit.icc_profile) {
		cmsCloseProfile(gfit.icc_profile);
		gfit.icc_profile = NULL;
	}
	gfit.icc_profile = copyICCProfile(target);
	set_source_information();
	notify_gfit_modified();
}

void on_icc_convertto_clicked(GtkButton* button, gpointer* user_data) {
	cmsUInt32Number gfit_colorspace = cmsGetColorSpace(gfit.icc_profile);
	cmsUInt32Number gfit_colorspace_channels = cmsChannelsOf(gfit_colorspace);
	cmsUInt32Number target_colorspace = cmsGetColorSpace(target);
	cmsUInt32Number target_colorspace_channels = cmsChannelsOf(target_colorspace);

	if (target_colorspace != cmsSigGrayData && target_colorspace != cmsSigRgbData) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Color space not supported"), _("Siril only supports representing the image in Gray or RGB color spaces at present. You cannot assign or convert to non-RGB color profiles"));
		return;
	}
	void *data = NULL;
	cmsUInt32Number srctype, desttype;
	size_t npixels = gfit.rx * gfit.ry;
	// convert to profile
	data = (gfit.type == DATA_FLOAT) ? (void *) gfit.fdata : (void *) gfit.data;
	srctype = get_planar_formatter_type(gfit_colorspace, gfit.type, FALSE);
	desttype = get_planar_formatter_type(target_colorspace, gfit.type, FALSE);
	cmsHTRANSFORM transform = NULL;
	if (gfit_colorspace_channels == target_colorspace_channels) {
		transform = cmsCreateTransform(gfit.icc_profile, srctype, target, desttype, gui.icc.rendering_intent, 0);
	} else {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Transform not supported"), _("Transforms between color spaces with different numbers of channels not yet supported - this is coming soon..."));
		return;
	}
	if (transform) {
		cmsDoTransform(transform, data, data, npixels);
		cmsDeleteTransform(transform);
		gfit.icc_profile = copyICCProfile(target);
		set_source_information();
		notify_gfit_modified();
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Failed to create colorspace transform"));
	}

}

void icc_channels_mismatch() {
	siril_message_dialog(GTK_MESSAGE_ERROR, _("ICC Profile Mismatch"), _("The number of channels in the image and in the profile do not match."));
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
			if (gfit.naxes[2] != 3) {
				icc_channels_mismatch();
				break;
			}
			if (target) {
				cmsCloseProfile(target);
			}
			target = srgb_linear();
			break;
		case SRGB_TRC:
			if (gfit.naxes[2] != 3) {
				icc_channels_mismatch();
				break;
			}
			if (target) {
				cmsCloseProfile(target);
			}
			target = srgb_trc();
			break;
		case REC2020_LINEAR:
			if (gfit.naxes[2] != 3) {
				icc_channels_mismatch();
				break;
			}
			if (target) {
				cmsCloseProfile(target);
			}
			target = rec2020_linear();
			break;
		case REC2020_TRC:
			if (gfit.naxes[2] != 3) {
				icc_channels_mismatch();
				break;
			}
			if (target) {
				cmsCloseProfile(target);
			}
			target = rec2020_trc();
			break;
		case GRAY_LINEAR:
			if (gfit.naxes[2] != 1) {
				icc_channels_mismatch();
				break;
			}
			if (target) {
				cmsCloseProfile(target);
			}
			target = gray_linear();
			break;
		case GRAY_SRGBTRC:
			if (gfit.naxes[2] != 1) {
				icc_channels_mismatch();
				break;
			}
			if (target) {
				cmsCloseProfile(target);
			}
			target = gray_srgbtrc();
			break;
		case GRAY_REC709TRC:
			if (gfit.naxes[2] != 1) {
				icc_channels_mismatch();
				break;
			}
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

void on_enable_icc_toggled(GtkToggleButton* button, gpointer user_data) {
	gboolean status = gtk_toggle_button_get_active(button);
	if (!status) {
		if (siril_confirm_dialog(_("Are you sure?"), _("Disabling color management will result in an inconsistent appearance of your images when viewed in Siril and in other applications!"), _("Accept"))) {
		gui.icc.available = FALSE;
		com.icc.available = FALSE;
		siril_log_color_message(_("Warning: color management disabled.\n"), "salmon");
		} else {
			gtk_toggle_button_set_active(button, TRUE);
		}
	} else {
		com.icc.available = (com.icc.working_linear && com.icc.mono_linear && com.icc.working_out && com.icc.mono_out);
		gui.icc.available = (com.icc.available); // && gui.icc_profile_rgb && gui.icc_profile_mono);
		if (!com.icc.available)
			siril_log_color_message(_("Error: could not enable color management.\n"), "red");
		else if (!gui.icc.available)
			siril_log_color_message(_("Warning: could not enable display color management.\n"), "red");
		else
			siril_log_color_message(_("Color management enabled...\n"), "green");
	}
	notify_gfit_modified();
}

void update_icc_enabled_widget(GtkWidget *dialog, gpointer user_data) {
	// Set description
	set_source_information();
	GtkToggleButton* active_button = (GtkToggleButton*) lookup_widget("enable_icc");
	gtk_toggle_button_set_active(active_button, com.icc.available);
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

void on_working_gamut_changed(GtkComboBox *combo, gpointer user_data) {
	int choice = gtk_combo_box_get_active(combo);
	GtkWidget *lin = lookup_widget("custom_icc_linear_trc");
	GtkWidget *std = lookup_widget("custom_icc_standard_trc");
	GtkWidget *gray = lookup_widget("custom_gray_icc_matching_trc");
	gtk_widget_set_sensitive(lin, (choice == 2));
	gtk_widget_set_sensitive(std, (choice == 2));
	gtk_widget_set_sensitive(gray, (choice == 2));
}

void on_icc_dialog_show(GtkWidget *dialog, gpointer user_data) {
	set_source_information();
	set_target_information();
}
