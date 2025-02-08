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
#ifndef SRC_CORE_ICC_PROFILE_H_
#define SRC_CORE_ICC_PROFILE_H_
#include <stdint.h>
#include <lcms2.h>

// Define some additional formatters that aren't defined in lcms2.h
#define TYPE_RGB_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_RGB)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))
#define TYPE_XYZ_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_XYZ)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))
#define TYPE_Lab_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_Lab)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))
#define TYPE_Luv_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_YUV)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))
#define TYPE_YCbCr_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_YCbCr)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))
#define TYPE_Yxy_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_Yxy)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))
#define TYPE_HSV_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_HSV)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))
#define TYPE_HLS_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_HLS)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))
#define TYPE_CMY_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_CMY)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))

#define TYPE_XYZ_16_PLANAR (COLORSPACE_SH(PT_XYZ)|CHANNELS_SH(3)|BYTES_SH(2)|PLANAR_SH(1))
#define TYPE_Lab_16_PLANAR (COLORSPACE_SH(PT_Lab)|CHANNELS_SH(3)|BYTES_SH(2)|PLANAR_SH(1))
#define TYPE_Luv_16_PLANAR (COLORSPACE_SH(PT_YUV)|CHANNELS_SH(3)|BYTES_SH(2)|PLANAR_SH(1))
#define TYPE_Yxy_16_PLANAR (COLORSPACE_SH(PT_Yxy)|CHANNELS_SH(3)|BYTES_SH(2)|PLANAR_SH(1))

typedef enum {
	NONE,
	SRGB_LINEAR,
	SRGB_TRC,
	REC2020_LINEAR,
	REC2020_TRC,
	GRAY_LINEAR,
	GRAY_SRGBTRC,
	GRAY_REC709TRC
} internal_icc;

typedef struct SirilMatrix3_d {
	gdouble coeff[3][3];
} SirilMatrix3_d;

void icc_profile_set_tag (cmsHPROFILE profile, cmsTagSignature sig, const gchar *tag);

cmsHPROFILE srgb_linear();
cmsHPROFILE gray_srgbtrc();
cmsHPROFILE srgb_trc();
cmsHPROFILE srgb_trcv2();
cmsHPROFILE gray_linear();
cmsHPROFILE rec2020_trc();
cmsHPROFILE rec2020_trcv2();
cmsHPROFILE rec2020_linear();
cmsHPROFILE gray_rec709trc();
cmsHPROFILE gray_srgbtrcv2();
cmsHPROFILE gray_rec709trcv2();
cmsHPROFILE srgb_monitor_perceptual();

void export_profile(cmsHPROFILE profile, const char *provided_filename);
void color_manage(fits *fit, gboolean active);
void lock_display_transform();
void unlock_display_transform();
void display_index_transform(BYTE* index, int vport);
cmsHTRANSFORM initialize_proofing_transform();
gboolean same_primaries(cmsHPROFILE a, cmsHPROFILE b, cmsHPROFILE c);
void reset_icc_transforms();
void validate_custom_profiles();
void initialize_profiles_and_transforms();
cmsUInt32Number get_planar_formatter_type(cmsColorSpaceSignature tgt, data_type t, gboolean force_16);
cmsHTRANSFORM initialize_export8_transform(fits* fit, gboolean threaded);
void refresh_icc_transforms();
char* siril_color_profile_get_copyright (cmsHPROFILE profile);
char* siril_color_profile_get_description (cmsHPROFILE profile);
char* siril_color_profile_get_manufacturer (cmsHPROFILE profile);
char* siril_color_profile_get_model (cmsHPROFILE profile);
gboolean siril_color_profile_is_rgb(cmsHPROFILE profile);
unsigned char* get_icc_profile_data(cmsHPROFILE profile, guint32 *len);
cmsBool fit_icc_is_linear(fits *fit);
cmsHPROFILE siril_color_profile_linear_from_color_profile (cmsHPROFILE profile);
void check_profile_correct(fits* fit);
void enable_profile_check_verbose();
void disable_profile_check_verbose();
cmsHPROFILE copyICCProfile(cmsHPROFILE profile);
void fits_initialize_icc(fits *fit, cmsUInt8Number* EmbedBuffer, cmsUInt32Number EmbedLen);
cmsBool profiles_identical(cmsHPROFILE a, cmsHPROFILE b);
void siril_colorspace_transform(fits *fit, cmsHPROFILE profile);
void icc_auto_assign_or_convert(fits *fit, icc_assign_type occasion);
void icc_auto_assign(fits *fit, icc_assign_type occasion);
const char* default_system_icc_path();
cmsHTRANSFORM sirilCreateTransformTHR(cmsContext Context, cmsHPROFILE Input, cmsUInt32Number InputFormat, cmsHPROFILE Output, cmsUInt32Number OutputFormat, cmsUInt32Number Intent, cmsUInt32Number dwFlags);
void update_profiles_after_gamut_change();
void siril_plot_colorspace(cmsHPROFILE profile, gboolean compare_srgb);
void cleanup_common_profiles();

#endif /* SRC_CORE_ICC_PROFILE_H_ */
