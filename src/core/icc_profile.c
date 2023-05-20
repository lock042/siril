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
/*
static cmsHPROFILE sRGB_g10;				// The working profile is always sRGB_g10
static cmsHPROFILE Gray_g10;				// or Gray_g10 for mono images
static cmsHPROFILE sRGB_g22_v4;				// This is the standard profile for display and for transforming RGB images before saving
static cmsHPROFILE Gray_g22_v4;				// This is the profile used for transforming mono images before saving
static cmsHPROFILE sRGB_g22_v2;				// This is the standard profile for embedding in saved RGB images
static cmsHPROFILE Gray_g22_v2;				// This is the standard profile for embedding in saved mono images
static cmsHPROFILE monitor_profile;			// This may store a user-provided monitor profile
static cmsBool monitor_profile_active; 		// This says whether or not to use the user-provided monitor profile
static cmsHPROFILE proof_profile;			// This may store a user-provided proof profile
static cmsBool proof_profile_active;		// This says whether or not the user-provided proof profile is active
static cmsHTRANSFORM displayTransform;		// This is the CMS transform used to remap the display
static cmsHTRANSFORM ucharTransform;		// This is the CMS transform used to remap the display if it holds an 8-bit image
static cmsHTRANSFORM ucharMonoTransform;	// This is the CMS transform used to remap the display if it holds an 8-bit image
static cmsHTRANSFORM ushrtTransform;		// This is the CMS transform used to remap the display
static cmsHTRANSFORM ushrtMonoTransform;		// This is the CMS transform used to remap the display
static cmsHTRANSFORM floatTransform;		// This is the CMS transform used to remap the display
static cmsHTRANSFORM floatMonoTransform;		// This is the CMS transform used to remap the display
*/
////// Functions //////

void initialize_profiles_and_transforms() {
//	cmsPlugin(cmsFastFloatExtensions());
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
	// Target profiles for display
//	gui.icc.sRGB = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G22], "r");
//	gui.icc.mono = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_GV4G22], "r");
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
			}
		}
		if (!gui.icc.monitor) {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G22], "r");
		}
		if (com.pref.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_CUSTOM_PROOF], "r");
			if (!gui.icc.soft_proof) {
				error++;
				siril_log_message(_("Warning: soft proofing profile set but could not be loaded. Soft proofing will be "
									"unavailable.\n"));
			}
		}
	}
	if (!error) {
		siril_log_message(_("ICC profiles loaded correctly. Workflow will be color managed.\n"));
	}
//	gui.icc.tx_display = cmsCreateTransform(gfit.icc_profile, TYPE_BGRA_8, gui.icc.srgb, TYPE_BGRA_8, INTENT_PERCEPTUAL, 0);
//	ushrtTransform = cmsCreateTransform(com.icc.mono_linear, TYPE_RGB_16, gui.icc.mono, TYPE_RGB_16, INTENT_PERCEPTUAL, 0);
//	ushrtMonoTransform = cmsCreateTransform(com.icc.mono_linear, TYPE_GRAY_16, Gray_g22_v4, TYPE_GRAY_16, INTENT_PERCEPTUAL, 0);
//	ucharTransform = cmsCreateTransform(sRGB_g10, TYPE_RGB_8, sRGB_g22_v4, TYPE_RGB_8, INTENT_PERCEPTUAL, 0);
//	ucharMonoTransform = cmsCreateTransform(com.icc.mono_linear, TYPE_GRAY_8, Gray_g22_v4, TYPE_GRAY_8, INTENT_PERCEPTUAL, 0);
//	floatTransform = cmsCreateTransform(sRGB_g10, TYPE_RGB_FLT, sRGB_g22_v4, TYPE_RGB_FLT, INTENT_PERCEPTUAL, 0);
//	floatMonoTransform = cmsCreateTransform(Gray_g10, TYPE_GRAY_FLT, Gray_g22_v4, TYPE_GRAY_FLT, INTENT_PERCEPTUAL, 0);
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
/*
BYTE uchar_pixel_icc_tx(BYTE in, int channel, int nchans) {
	g_assert(channel >= 0 && channel < 3);
	g_assert(channel < nchans);
	BYTE input[3] = { 0 };
	BYTE output[3] = { 0 };
	input[channel] = in;
	if (nchans == 1)
		cmsDoTransform(gfit.display_transform, input, output, 1);
	else
		cmsDoTransform(ucharTransform, input, output, 1);
	return output[channel];
}

WORD ushrt_pixel_icc_tx(WORD in, int channel, int nchans) {
	g_assert(channel >= 0 && channel < 3);
	g_assert(channel < nchans);
	WORD input[3] = { 0 };
	WORD output[3] = { 0 };
	input[channel] = in;
	if (nchans == 1)
		cmsDoTransform(ushrtMonoTransform, input, output, 1);
	else
		cmsDoTransform(ushrtTransform, input, output, 1);
	return output[channel];
}

float float_pixel_icc_tx(float in, int channel, int nchans) {
	g_assert(channel >= 0 && channel < 3);
	g_assert(channel < nchans);
	float input[3] = { 0 };
	float output[3] = { 0 };
	input[channel] = in;
	if (nchans == 2)
		cmsDoTransform(floatMonoTransform, input, output, 1);
	else
		cmsDoTransform(floatTransform, input, output, 1);
	return output[channel];
}
*/
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

// This function is used in converting buffers from fit.data or fit.fdata to a Gray or sRGB g22 profile prior to saving
/*
int transformBufferOnSave(void* src, void* dest, uint16_t src_bitspersample, uint16_t dest_bitspersample, uint16_t nsamples, size_t npixels, gboolean planar, gboolean linear) {
	cmsHPROFILE src_profile = nsamples == 1 ? Gray_g10 : sRGB_g10;
	cmsHPROFILE dest_profile;
	cmsHTRANSFORM transform;
	cmsUInt32Number src_type, dest_type;
	gboolean src_is_float = (src_bitspersample == 32);

	if (src_is_float) {
		// 32-bit
		if (nsamples == 1 || nsamples == 2) {
			// mono
			src_type = TYPE_GRAY_FLT;
			dest_profile = linear ? Gray_g10 : Gray_g22_v4;
		} else {
			// RGB
			src_type = planar ? TYPE_RGB_FLT_PLANAR : TYPE_RGB_FLT;
			dest_profile = linear ? sRGB_g10: sRGB_g22_v4;
		}
	} else {
		// 16-bit
		if (nsamples == 1 || nsamples == 2) {
			// mono
			src_type = TYPE_GRAY_16;
			dest_profile = linear ? Gray_g10 : Gray_g22_v4;
		} else {
			// RGB
			src_type = planar ? TYPE_RGB_16_PLANAR : TYPE_RGB_16;
			dest_profile = linear ? sRGB_g10 : sRGB_g22_v4;
		}
	}
	gboolean mono = (nsamples == 1 || nsamples == 2);
	switch (dest_bitspersample) {
		case 8:
			dest_type = mono ? TYPE_GRAY_8 : planar ? TYPE_RGB_8_PLANAR : TYPE_RGB_8;
			break;
		case 16:
			dest_type = mono ? TYPE_GRAY_16 : planar ? TYPE_RGB_16_PLANAR : TYPE_RGB_16;
			break;
		case 32:
			dest_type = mono ? TYPE_GRAY_FLT : planar ? TYPE_RGB_FLT_PLANAR : TYPE_RGB_FLT;
			break;
		default:
			return 1;
	}
	transform = cmsCreateTransform(src_profile, src_type, dest_profile, dest_type, INTENT_PERCEPTUAL, 0);
	cmsDoTransform(transform, src, dest, npixels);
	cmsDeleteTransform(transform);
	return 0;
}
*/
// This function is used in converting a loaded file into the working color profile
/*
void transformBufferOnLoad(void* buf, uint16_t bitdepth, cmsUInt8Number* EmbedBuffer, cmsUInt32Number EmbedLen, uint16_t nsamples, size_t npixels) {
	// Convert the buffer into a profile in memory
	cmsHPROFILE src_profile = cmsOpenProfileFromMem(EmbedBuffer, EmbedLen);
	cmsHPROFILE dest_profile;
	cmsHTRANSFORM transform;
	cmsUInt32Number src_type, dest_type;
	// Work out the src and dest data formats
	if (bitdepth != 32) {
		// 8- or 16-bit
		if (nsamples == 1 || nsamples == 2) {
			// mono
			src_type = dest_type = bitdepth == 16 ? TYPE_GRAY_16 : TYPE_GRAY_8;
			dest_profile = Gray_g10;
		} else {
			// RGB
			src_type = dest_type = bitdepth == 16 ? TYPE_RGB_16_PLANAR : TYPE_RGB_8_PLANAR;
			dest_profile = sRGB_g10;
		}
	} else {
		// 32-bit
		if (nsamples == 1 || nsamples == 2) {
			// mono
			src_type = dest_type = TYPE_GRAY_FLT;
			dest_profile = Gray_g10;
		} else {
			// RGB
			src_type = dest_type = TYPE_RGB_FLT_PLANAR;
			dest_profile = sRGB_g10;
		}
	}
	transform = cmsCreateTransform(src_profile, src_type, dest_profile, dest_type, INTENT_PERCEPTUAL, 0);
	cmsDoTransform(transform, buf, buf, npixels);
	cmsDeleteTransform(transform);
	cmsCloseProfile(src_profile);
}
*/
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
			fit->icc_profile = copyICCProfile((fit->naxes[2] == 1) ? com.icc.mono_standard : com.icc.srgb_standard);
		}
		if (gui.icc.available && fit->icc_profile) {
			fit->display_transform = cmsCreateTransform(fit->icc_profile, TYPE_BGRA_8, gui.icc.monitor, TYPE_BGRA_8, INTENT_PERCEPTUAL, 0);
		}
	}
}
