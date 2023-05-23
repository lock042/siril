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
#include "algos/lcms_fast_float/lcms2_fast_float.h"
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
#define INDEX_CUSTOM_EXPORT 8

cmsHPROFILE copyICCProfile(cmsHPROFILE profile);

const char* default_icc_paths[] = { DEFAULT_PATH_GV2g22, DEFAULT_PATH_GV4g10, DEFAULT_PATH_GV4g22, DEFAULT_PATH_RGBV2g22, DEFAULT_PATH_RGBV4g10, DEFAULT_PATH_RGBV4g22 };
static GMutex monitor_profile_mutex;
static GMutex soft_proof_profile_mutex;

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
	// Enablethe fast float plugin
	cmsPlugin(cmsFastFloatExtensions());

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
			gtk_widget_set_sensitive(GTK_WIDGET(lookup_widget("soft_proof_item")), FALSE);
		}

	}
	if (!error) {
		siril_log_message(_("ICC profiles loaded correctly. Workflow will be color managed.\n"));
	}
}

void refresh_icc_settings() {
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

cmsHTRANSFORM initialize_display_transform() {
	g_assert(gui.icc.monitor);
	cmsHTRANSFORM transform = NULL;
	if (gfit.icc_profile == NULL)
		return NULL;
	cmsUInt32Number type = (gfit.naxes[2] == 1 ? TYPE_GRAY_16 : TYPE_RGB_16_PLANAR);
	g_mutex_lock(&monitor_profile_mutex);
	transform = cmsCreateTransform(gfit.icc_profile, type, (gfit.naxes[2] == 3 ? gui.icc.monitor : com.icc.mono_out), type, gui.icc.rendering_intent, 0);
	g_mutex_unlock(&monitor_profile_mutex);
	if (transform == NULL)
		siril_log_message("Error: failed to create display_transform!\n");
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
// that use the profiles need to go inside mutexes to prevent the profiles being ripped out from under them.
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
