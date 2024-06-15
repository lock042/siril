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
#ifndef SRC_CORE_ICC_PROFILE_H_
#define SRC_CORE_ICC_PROFILE_H_
#include <stdint.h>
#include <lcms2.h>

#define ICC_COPYRIGHT_TEXT "Copyright 2024, Team free-astro (http://siril.org/), CC-0 (https://creativecommons.org/publicdomain/zero/1.0/legalcode)."
#define SRGB_MANUFACTURER_TEXT "sRGB chromaticities from A Standard Default Color Space for the Internet - sRGB, http://www.w3.org/Graphics/Color/sRGB; also see http://www.color.org/specification/ICC1v43_2010-12.pdf"
#define P3_MANUFACTURER_TEXT "P3 chromaticities from https://en.wikipedia.org/wiki/DCI-P3"
#define ADOBE_MANUFACTURER_TEXT "Chromaticities as given in Adobe RGB (1998) Color Image Encoding, Version 2005-05, https://www.adobe.com/digitalimag/pdfs/AdobeRGB1998.pdf"
#define REC2020_MANUFACTURER_TEXT "Rec2020 chromaticities from https://www.itu.int/dms_pub/itu-r/opb/rep/R-REP-BT.2246-2-2012-PDF-E.pdf; https://www.itu.int/dms_pubrec/itu-r/rec/bt/R-REC-BT.2020-2-201510-I!!PDF-E.pdf"
#define ROMM_MANUFACTURER_TEXT "Chromaticities from Reference Input/Output Medium Metric RGB Color Encodings (RIMM/ROMM RGB), http://photo-lovers.org/pdf/color/romm.pdf"
#define ACESCG_MANUFACTURER_TEXT "ACEScg chromaticities from S-2014-004 v1.0.1, http://www.oscars.org/science-technology/aces/aces-documentation"
#define ACES_MANUFACTURER_TEXT "ACES chromaticities from TB-2014-004, http://www.oscars.org/science-technology/aces/aces-documentation"

#define SRGBPARAMS { 2.4, 1.0 / 1.055,  0.055 / 1.055, 1.0 / 12.92, 0.04045 }
#define REC709PARAMS { 1.0 / 0.45, 1.0 / 1.099,  0.099 / 1.099,  1.0 / 4.5, 0.081 }
#define LABLPARAMS { 3.0, 1.0 / 1.16,  0.16 / 1.16, 2700.0 / 24389.0, 0.08000 }

#define PRIMARIES_SRGB (cmsCIExyYTRIPLE){{0.639998686, 0.330010138, 1.0},{0.300003784, 0.600003357, 1.0},{0.150002046, 0.059997204, 1.0}}
#define PRIMARIES_P3 (cmsCIExyYTRIPLE){{0.680, 0.302, 1.0},{0.265, 0.690, 1.0},{0.150, 0.060, 1.0}}
#define PRIMARIES_ADOBE (cmsCIExyYTRIPLE){{0.639996511, 0.329996864, 1.0},{0.210005295, 0.710004866, 1.0},{0.149997606, 0.060003644, 1.0}}
#define PRIMARIES_REC2020 (cmsCIExyYTRIPLE){{0.708012540607, 0.291993664388, 1.0},{0.169991652439, 0.797007778423, 1.0},{0.130997824007, 0.045996550894, 1.0}}
#define PRIMARIES_ROMM (cmsCIExyYTRIPLE){{0.7347, 0.2653, 1.0},{0.1596, 0.8404, 1.0},{0.0366, 0.0001, 1.0}}
#define PRIMARIES_ACES_CG (cmsCIExyYTRIPLE){{0.713, 0.293,  1.0},{0.165, 0.830,  1.0},{0.128, 0.044,  1.0}}
#define PRIMARIES__ACES (cmsCIExyYTRIPLE){{0.734704192222, 0.265298276252,  1.0},{-0.000004945077, 0.999992850272,  1.0},{0.000099889199, -0.077007518685, 1.0}}

#define ROMMSPEC_WHITEPOINT (cmsCIExyY){0.3457, 0.3585, 1.0}
#define D50_ILLUMINANT_WHITEPOINT (cmsCIExyY){0.345702915, 0.358538597, 1.0}
#define D60_WHITEPOINT (cmsCIExyY){0.32168, 0.33767, 1.0}
#define D65_SRGB_WHITEPOINT (cmsCIExyY){0.3127, 0.3290, 1.0}

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

cmsHPROFILE sirilCreateRGBProfileV4(cmsCIExyY *whitepoint,
								  cmsCIExyYTRIPLE *primaries,
								  cmsToneCurve** curve,
								  const char* manufacturer_text,
								  const char* description_text);
cmsHPROFILE make_default_srgb_profile(gboolean is_linear);
cmsHPROFILE make_default_display_p3_profile(gboolean is_linear);
cmsHPROFILE make_default_adobergb_compat_profile(gboolean is_linear);
cmsHPROFILE make_default_rec2020_profile(gboolean is_linear);
cmsHPROFILE make_default_prophoto_compat_profile(gboolean is_linear);
cmsHPROFILE make_default_rec709_mono_profile(gboolean is_linear);
cmsHPROFILE make_default_srgb_mono_profile(gboolean is_linear);
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

#endif /* SRC_CORE_ICC_PROFILE_H_ */
