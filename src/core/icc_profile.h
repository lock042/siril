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

cmsHPROFILE srgb_linear();
cmsHPROFILE gray_srgbtrc();
cmsHPROFILE srgb_trc();
cmsHPROFILE gray_linear();

void validate_custom_profiles();
void initialize_profiles_and_transforms();
cmsUInt32Number get_planar_formatter_type(cmsColorSpaceSignature tgt, data_type t, gboolean force_16);
void assign_linear_icc_profile(fits *fit);
cmsHTRANSFORM initialize_display_transform();
cmsHTRANSFORM initialize_export8_transform(fits* fit, gboolean threaded);
cmsHTRANSFORM initialize_proofing_transform();
void refresh_icc_transforms();
unsigned char* get_icc_profile_data(cmsHPROFILE profile, guint32 *len);
cmsBool fit_icc_is_linear(fits *fit);
void check_profile_correct(fits* fit);
cmsHPROFILE copyICCProfile(cmsHPROFILE profile);
void fits_initialize_icc(fits *fit, cmsUInt8Number* EmbedBuffer, cmsUInt32Number EmbedLen);
cmsBool profiles_identical(cmsHPROFILE a, cmsHPROFILE b);
void convert_fit_colorspace(fits *fit, cmsHPROFILE *from, cmsHPROFILE *to);
void convert_fit_colorspace_to_reference_fit(fits* input, fits* reference);
void set_icc_description_in_TIFF();
cmsHPROFILE siril_construct_split_profile(fits *from, int channel);
void check_linear_and_convert_with_approval(fits *fit);
const char* default_system_icc_path();
cmsHTRANSFORM sirilCreateTransformTHR(cmsContext Context, cmsHPROFILE Input, cmsUInt32Number InputFormat, cmsHPROFILE Output, cmsUInt32Number OutputFormat, cmsUInt32Number Intent, cmsUInt32Number dwFlags);
void update_profiles_after_gamut_change();
void initialize_icc_preferences_widgets();

#endif /* SRC_CORE_ICC_PROFILE_H_ */
