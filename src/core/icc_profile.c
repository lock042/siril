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
// that may be stored in com.pref but for which there are no defaults
#define INDEX_CUSTOM_MONITOR 6
#define INDEX_CUSTOM_PROOF 7
#define INDEX_CUSTOM_EXPORT 8

cmsHPROFILE copyICCProfile(cmsHPROFILE profile);

const char* default_icc_paths[] = { DEFAULT_PATH_GV2g22, DEFAULT_PATH_GV4g10, DEFAULT_PATH_GV4g22, DEFAULT_PATH_RGBV2g22, DEFAULT_PATH_RGBV4g10, DEFAULT_PATH_RGBV4g22 };
static GMutex monitor_profile_mutex;
static GMutex soft_proof_profile_mutex;

static cmsHPROFILE target = NULL; // Target profile for the GUI tool

////// Functions //////

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

void initialize_profiles_and_transforms() {
	// Enable the fast float plugin (as long as the OS / lcms2 version blacklist isn't triggered)
#ifndef EXCLUDE_FF
	cmsPlugin(cmsFastFloatExtensions());
#endif
	// Initialize paths to standard ICC profiles
	initialize_icc_profiles_paths();

	// Set alarm codes for out-of-gamut warning
	cmsUInt16Number alarmcodes[16] = { 65535, 0, 65535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	cmsSetAlarmCodes(alarmcodes); // Out of gamut colours will be shown in magenta
	int error = 0;

	// Intents
	com.icc.save_intent = com.pref.export_intent;
	gui.icc.rendering_intent = com.pref.rendering_intent;
	gui.icc.proofing_intent = com.pref.proofing_intent;

	// Linear working profiles
	com.icc.srgb_linear = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G10], "r");
	com.icc.mono_linear = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_GV4G10], "r");
	com.icc.mono_standard = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_GV4G22], "r");
	com.icc.srgb_standard = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G22], "r");

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
		g_mutex_lock(&monitor_profile_mutex);
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
		g_mutex_unlock(&monitor_profile_mutex);

		if (com.pref.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
			g_mutex_lock(&soft_proof_profile_mutex);
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_CUSTOM_PROOF], "r");
			g_mutex_unlock(&soft_proof_profile_mutex);
		if (!gui.icc.soft_proof) {
				error++;
				siril_log_message(_("Warning: soft proofing profile set but could not be loaded. Soft proofing will be "
									"unavailable.\n"));
			} else {
				siril_log_message(_("Soft proofing ICC profile loaded from %s\n"), com.pref.icc_paths[INDEX_CUSTOM_PROOF]);
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

// Loads a custom monitor profile from a path in com.pref.icc_paths
// The path must be set by the user in preferences, there is no default custom monitor profile
int load_monitor_icc_profile(const char* filename) {
	if (com.pref.icc_paths[INDEX_CUSTOM_MONITOR] && com.pref.icc_paths[INDEX_CUSTOM_MONITOR][0] != '\0') {
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

// Loads a custom proof profile from a path in com.pref.icc_paths
// The path must be set by the user in preferences, there is no default custom proof profile
int load_soft_proof_icc_profile(const char* filename) {
	if (com.pref.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
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
	fit->icc_profile = copyICCProfile(gfit.naxes[2] == 1 ? com.icc.mono_linear : com.icc.srgb_linear);
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
	g_assert(com.icc.srgb_standard);
	cmsHTRANSFORM transform = NULL;
	if (fit->icc_profile == NULL)
		return NULL;
	cmsUInt32Number fit_signature = cmsGetColorSpace(fit->icc_profile);
	cmsUInt32Number srctype = get_planar_formatter_type(fit_signature, fit->type, TRUE);
	cmsUInt32Number desttype = (fit->naxes[2] == 1 ? TYPE_GRAY_16 : TYPE_RGB_16_PLANAR);
	transform = cmsCreateTransform(fit->icc_profile, srctype, (fit->naxes[2] == 3 ? com.icc.srgb_standard : com.icc.mono_standard), desttype, gui.icc.rendering_intent, 0);
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
		fit->icc_profile = copyICCProfile((fit->naxes[2] == 1) ? com.icc.mono_linear : com.icc.srgb_linear);
	}
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
	if (com.pref.icc_paths[INDEX_CUSTOM_MONITOR] && com.pref.icc_paths[INDEX_CUSTOM_MONITOR][0] != '\0') {
		g_free(com.pref.icc_paths[INDEX_CUSTOM_MONITOR]);
		com.pref.icc_paths[INDEX_CUSTOM_MONITOR] = NULL;
		cmsCloseProfile(gui.icc.monitor);
		if (com.icc.available) {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G22], "r");
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
	if (com.pref.icc_paths[INDEX_CUSTOM_PROOF] && com.pref.icc_paths[INDEX_CUSTOM_PROOF][0] != '\0') {
		g_free(com.pref.icc_paths[INDEX_CUSTOM_PROOF]);
		com.pref.icc_paths[INDEX_CUSTOM_PROOF] = NULL;
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
		if (!com.pref.icc_paths[INDEX_CUSTOM_MONITOR] || com.pref.icc_paths[INDEX_CUSTOM_MONITOR][0] == '\0') {
			com.pref.icc_paths[INDEX_CUSTOM_MONITOR] = g_strdup(gtk_file_chooser_get_filename(filechooser));
		}
		if (!com.pref.icc_paths[INDEX_CUSTOM_MONITOR] || com.pref.icc_paths[INDEX_CUSTOM_MONITOR][0] == '\0') {
			siril_log_color_message(_("Error: no filename specfied for custom monitor profile.\n"), "red");
			no_file = TRUE;
		} else {
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_CUSTOM_MONITOR], "r");
		}
		if (gui.icc.monitor) {
			siril_log_message(_("Monitor profile loaded from %s\n"), com.pref.icc_paths[INDEX_CUSTOM_MONITOR], "r");
			g_mutex_unlock(&monitor_profile_mutex);
			refresh_icc_transforms();
			return;
		} else {
			if (!no_file) {
				siril_log_color_message(_("Monitor profile could not be loaded from %s\n"), "red", com.pref.icc_paths[INDEX_CUSTOM_MONITOR]);
			}
			gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G22], "r");
			if (gui.icc.monitor) {
				siril_log_message(_("Monitor ICC profile set to sRGB (D65 whitepoint, gamma = 2.2)\n"));
			} else {
				siril_log_color_message(_("Error: standard sRGB ICC profile could not be loaded. Color management is disabled.\n"), "red");
				com.icc.available = FALSE;
			}
		}
	} else {
		gui.icc.monitor = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_SRGBV4G22], "r");
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
	gboolean no_file;
	gboolean active = gtk_toggle_button_get_active(button);
	g_mutex_lock(&soft_proof_profile_mutex);
	if (gui.icc.soft_proof) {
		cmsCloseProfile(gui.icc.soft_proof);
	}
	if (active) {
		if (!com.pref.icc_paths[INDEX_CUSTOM_PROOF] || com.pref.icc_paths[INDEX_CUSTOM_PROOF][0] == '\0') {
			com.pref.icc_paths[INDEX_CUSTOM_PROOF] = g_strdup(gtk_file_chooser_get_filename(filechooser));
		}
		if (!com.pref.icc_paths[INDEX_CUSTOM_PROOF] || com.pref.icc_paths[INDEX_CUSTOM_PROOF][0] == '\0') {
			siril_log_color_message(_("Error: no filename specfied for custom proofing profile.\n"), "red");
			no_file = TRUE;
		} else {
			gui.icc.soft_proof = cmsOpenProfileFromFile(com.pref.icc_paths[INDEX_CUSTOM_PROOF], "r");
		}
		if (gui.icc.soft_proof) {
			siril_log_message(_("Soft proofing profile loaded from %s\n"), com.pref.icc_paths[INDEX_CUSTOM_PROOF]);
			g_mutex_unlock(&soft_proof_profile_mutex);
			refresh_icc_transforms();
			return;
		} else {
			if (!no_file) {
				siril_log_color_message(_("Soft proofing profile could not be loaded from %s\n"), "red", com.pref.icc_paths[INDEX_CUSTOM_PROOF]);
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

//////// GUI callbacks for the color management dialog
void set_source_information() {
	if (!gfit.icc_profile) {
		siril_debug_print("Error: target profile is NULL\n");
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
		siril_debug_print("Error: target profile is NULL\n");
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

void on_icc_apply_clicked(GtkButton* button, gpointer* user_data) {
	GtkComboBox* operation = (GtkComboBox*) lookup_widget("icc_operation_combo");
	int ui_operation = gtk_combo_box_get_active(operation);
	cmsUInt32Number gfit_colorspace = cmsGetColorSpace(gfit.icc_profile);
	cmsUInt32Number gfit_colorspace_channels = cmsChannelsOf(gfit_colorspace);
	cmsUInt32Number target_colorspace = cmsGetColorSpace(target);
	cmsUInt32Number target_colorspace_channels = cmsChannelsOf(target_colorspace);

/*	if (target_colorspace != cmsSigGrayData && target_colorspace != cmsSigRgbData) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Color space not supported"), _("Siril only supports representing the image in Gray or RGB color spaces at present. You cannot assign or convert to non-RGB color profiles"));
		return;
	}*/
	void *data = NULL;
	cmsUInt32Number srctype, desttype;
	size_t npixels = gfit.rx * gfit.ry;
	switch(ui_operation) {
		case 0:
			// assign profile
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
			break;
		case 1:
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
			break;
		default:
			siril_message_dialog(GTK_MESSAGE_WARNING, _("No operation selected"), _("Choose either \"Assign profile\" or \"Convert to profile\" from the dropdown."));
			return;
	}
}

void on_icc_target_combo_changed(GtkComboBox* combo, gpointer* user_data) {
	GtkFileChooser* filechooser = (GtkFileChooser*) lookup_widget("icc_target_filechooser");
	int target_index = gtk_combo_box_get_active(combo);
	GtkLabel* label = (GtkLabel*) lookup_widget("icc_target_profile_label");
	gtk_file_chooser_unselect_all(filechooser);
	switch (target_index) {
		case 0:
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
		case 1:
			if (target) {
				cmsCloseProfile(target);
			}
			if (gfit.naxes[2] == 1) {
				target = copyICCProfile(com.icc.mono_linear);
			} else {
				target = copyICCProfile(com.icc.srgb_linear);
			}
			break;
		case 2:
			if (target) {
				cmsCloseProfile(target);
			}
			if (gfit.naxes[2] == 1) {
				target = copyICCProfile(com.icc.mono_standard);
			} else {
				target = copyICCProfile(com.icc.srgb_standard);
			}
			break;
	}
	set_target_information();
}

void on_icc_operation_combo_changeded(GtkComboBox* combo, gpointer* user_data) {
}

void on_icc_target_filechooser_file_set(GtkFileChooser* filechooser, gpointer* user_data) {
	GtkComboBox* target_combo = (GtkComboBox*) lookup_widget("icc_target_combo");
	gtk_combo_box_set_active(target_combo, 0);
	gchar *filename = gtk_file_chooser_get_filename(filechooser);
	if (filename) {
		target = cmsOpenProfileFromFile(filename, "r");
	} else {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not load file"), _("Error:could not load selected file, or it does not contain a valid ICC profile."));
		g_free(filename);
		gtk_file_chooser_unselect_all(filechooser);
		return;
	}
	if (!target) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not load file"), _("Error:could not load selected file, or it does not contain a valid ICC profile."));
		gtk_file_chooser_unselect_all(filechooser);
	} else {
		set_target_information();
	}
	g_free(filename);
}

void on_icc_dialog_show(GtkWidget *dialog, gpointer user_data) {
	// Set description
	set_source_information();
}
