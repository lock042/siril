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

// This file includes some code taken from Elle Stone's
// "make-elles-profiles.c" program (https://github.com/ellelstone/elles_icc_profiles/blob/master/code/make-elles-profiles.c)
// Particularly important is the definitions of the pre-quantized
// sRGB primaries which avoid issues with lcms
// This is licenced under the GPL (version 2 or any later version)
// Therefore it is included here and relicenced under the terms of the
// GPL (version 3 or any later version to align with the rest of Siril.

#include <glib.h>
#include "core/siril.h"
#include "core/OS_utils.h"
#include "icc_profile.h"
#include "gui/utils.h"
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
// that may be stored in com.pref but for which there are no defaults
#define INDEX_CUSTOM_MONITOR 6
#define INDEX_CUSTOM_PROOF 7

cmsHPROFILE copyICCProfile(cmsHPROFILE profile);

const char* default_icc_paths[] = { DEFAULT_PATH_GV2g22, DEFAULT_PATH_GV4g10, DEFAULT_PATH_GV4g22, DEFAULT_PATH_RGBV2g22, DEFAULT_PATH_RGBV4g10, DEFAULT_PATH_RGBV4g22 };

////// Functions //////

void initialize_profiles_and_transforms() {
//	cmsPlugin(cmsFastFloatExtensions());
	initialize_icc_profiles_paths();
	// Set alarm codes for out-of-gamut warning
	cmsUInt16Number alarmcodes[16] = { 65535, 0, 65535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	cmsSetAlarmCodes(alarmcodes); // Out of gamut colours will be shown in magenta
	int error = 0;
	// Intents
	com.icc.save_intent = INTENT_RELATIVE_COLORIMETRIC;
	gui.icc.rendering_intent = INTENT_PERCEPTUAL;
	// Linear working profiles
	com.icc.srgb_linear = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G10], "r");
	com.icc.mono_linear = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_GV4G10], "r");
	// Target profile for embedding in saved RGB and mono files
	com.icc.srgb_out = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV2G22], "r");
	com.icc.mono_out = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_GV2G22], "r");
	// ICC availability
	com.icc.available = (com.icc.srgb_linear && com.icc.mono_linear && com.icc.srgb_out && com.icc.mono_out);
	gui.icc.available = (com.icc.available); // && gui.icc_profile_rgb && gui.icc_profile_mono);
	if ((com.headless && !com.icc.available) || (!com.headless && !gui.icc.available)) {
		error++;
		siril_log_message(_("Error: standard color management profiles could not be loaded. Siril will continue "
							"but colors in saved images will not be consistent when viewed in other applications "
							"or printed. Check your installation is correct.\n"));
	} else { // No point loading monitor and soft proof profiles if the standard ones are unavailable
		// Open the custom monitor and soft proofing profiles if there is a path set in preferences
		if (com.pref.icc_paths[INDEX_CUSTOM_MONITOR] && com.pref.icc_paths[INDEX_CUSTOM_MONITOR][0] != '\0') {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_CUSTOM_MONITOR], "r");
			if (!gui.icc.monitor) {
				error++;
				siril_log_message(_("Warning: custom monitor profile set but could not be loaded. Display will use a "
									"standard sRGB profile with D65 white point and gamma = 2.2.\n"));
			} else {
				siril_log_message(_("Monitor ICC profile loaded from %s\n"), com.pref.icc_paths[INDEX_CUSTOM_MONITOR]);
			}
		}
		if (!gui.icc.monitor) {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G22], "r");
			siril_log_message(_("Monitor ICC profile set to sRGB (D65 whitepoint, gamma = 2.2)\n"));
		}
		if (com.pref.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_CUSTOM_PROOF], "r");
			if (!gui.icc.soft_proof) {
				error++;
				siril_log_message(_("Warning: soft proofing profile set but could not be loaded. Soft proofing will be "
									"unavailable.\n"));
			} else {
				siril_log_message(_("Soft proofing ICC profile loaded from %s\n"), com.pref.icc_paths[INDEX_CUSTOM_PROOF]);
			}
		} else {
			siril_log_message(_("No soft proofing ICC profile set. Soft proofing is unavailable.\n"));
			gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("soft_proof_item")), FALSE);
		}

	}
	if (!error) {
		siril_log_message(_("ICC profiles loaded correctly. Workflow will be color managed.\n"));
	}
}

// Loads a custom monitor profile from a path in com.pref.icc_paths
// The path must be set by the user in preferences, there is no default custom monitor profile
int load_monitor_icc_profile(const char* filename) {
	if (com.pref.icc_paths[INDEX_CUSTOM_MONITOR] && com.pref.icc_paths[INDEX_CUSTOM_MONITOR][0] != '\0') {
		if (gui.icc.monitor)
			cmsCloseProfile(gui.icc.monitor);
		gui.icc.monitor = cmsOpenProfileFromFile(filename, "r");
		if (!gui.icc.monitor)
			return 1;
		else
			return 0;
	} else return 1;
}

// Loads a custom proof profile from a path in com.pref.icc_paths
// The path must be set by the user in preferences, there is no default custom proof profile
int load_soft_proof_icc_profile(const char* filename) {
	if (com.pref.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
		gui.icc.soft_proof = cmsOpenProfileFromFile(filename, "r");
		if (!gui.icc.soft_proof)
			return 1;
		else
			return 0;
	} else return 1;
}

// Adapted from initialize_local_catalogues_paths()
void initialize_icc_profiles_paths() {
	int nb_icc = sizeof(default_icc_paths) / sizeof(const char *);
	int maxpath = get_pathmax();
	for (int icc = 0; icc < nb_icc; icc++) {
		if (com.pref.icc_paths[icc] &&
				com.pref.icc_paths[icc][0] != '\0')
			continue;
		char path[maxpath];
		gchar *filename = g_build_filename(siril_get_system_data_dir(), "icc", default_icc_paths[icc], NULL);
		strncpy(path, filename, maxpath - 1);
		com.pref.icc_paths[icc] = g_strdup(path);
		g_free(filename);
	}
}

void display_profile_transform(const void* src, void* dest, cmsUInt32Number pixels) {
	if (gfit.display_transform) {
		cmsDoTransform(gfit.display_transform, src, dest, pixels);
		return;
	}
	if (!gfit.icc_profile) {
		gfit.icc_profile = copyICCProfile(gfit.naxes[2] == 1 ? com.icc.mono_linear : com.icc.srgb_linear);
		gfit.display_transform = cmsCreateTransform(gfit.icc_profile, TYPE_BGRA_8, gui.icc.monitor, TYPE_BGRA_8, INTENT_PERCEPTUAL, 0);
		if (gfit.display_transform)
			cmsDoTransform(gfit.display_transform, src, dest, pixels);
	}
}

void assign_linear_icc_profile(fits *fit) {
	if (fit->icc_profile) {
		cmsCloseProfile(fit->icc_profile);
	}
	fit->icc_profile = copyICCProfile(gfit.naxes[2] == 1 ? com.icc.mono_linear : com.icc.srgb_linear);
}

cmsHTRANSFORM initialize_display_transform() {
	g_assert(gfit.icc_profile);
	g_assert(gui.icc.monitor);
	int norm = (int) get_normalized_value(&gfit);
	cmsHTRANSFORM transform;
	cmsUInt32Number type;
	switch (norm) {
		case 1:
			type = TYPE_RGB_FLT;
			break;
		case UCHAR_MAX:
			type = TYPE_RGB_8;
			break;
		default:
		// includes case USHRT_MAX:
			type = TYPE_RGB_16;
	}
	transform = cmsCreateTransform(gfit.icc_profile, type, gui.icc.monitor, type, gui.icc.rendering_intent, 0);
	if (transform == NULL)
		siril_log_message("Error: failed to create display_transform!\n");
	return transform;
}

cmsHTRANSFORM initialize_proofing_transform() {
	g_assert(gfit.icc_profile);
	g_assert(gui.icc.monitor);
	g_assert(gui.icc.soft_proof);
	cmsUInt32Number flags = cmsFLAGS_SOFTPROOFING;

	// Cut gamut check disabled, it tends to show everything as out of gamut and doesn't
	// work with anything over 16-bit precision
/*	if (gui.cut_over) {
		flags |= cmsFLAGS_GAMUTCHECK;
	}*/
	cmsUInt32Number type;
	int norm = (int) get_normalized_value(&gfit);
	switch (norm) {
		case 1:
			type = TYPE_RGB_FLT;
			break;
		case UCHAR_MAX:
			type = TYPE_RGB_8;
			break;
		default:
			// includes case USHRT_MAX:
			type = TYPE_RGB_16;
	}
	cmsHPROFILE proofing_transform = cmsCreateProofingTransform(
						gfit.icc_profile,
						type,
						gui.icc.monitor,
						type,
						gui.icc.soft_proof,
						gui.icc.rendering_intent,
						gui.icc.proofing_intent,
						flags);
	return proofing_transform;
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
	cmsHPROFILE *profile = linear ? &com.icc.srgb_linear : &com.icc.srgb_out;
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
	cmsBool ret;
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
