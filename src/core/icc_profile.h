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

//#include <lcms2_fast_float.h>

#define TYPE_RGB_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_RGB)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1))

unsigned char* get_sRGB_profile_data(guint32 *len, gboolean linear);
unsigned char* get_gray_profile_data(guint32 *len, gboolean linear);

//unsigned char* get_profile_buf(cmsHPROFILE* profile, uint32_t* profile_len);

void refresh_icc_settings();
void initialize_profiles_and_transforms();
void initialize_icc_preferences_widgets();
void assign_linear_icc_profile(fits *fit);
cmsHTRANSFORM initialize_display_transform();
cmsHTRANSFORM initialize_proofing_transform();
void refresh_icc_transforms();
int load_display_icc_profile(const char* filename);
int load_proof_icc_profile(const char* filename);
void initialize_icc_profiles_paths();
void transformBufferOnLoad(void* buf, uint16_t bitdepth, cmsUInt8Number* EmbedBuffer, cmsUInt32Number EmbedLen, uint16_t nsamples, size_t npixels);
int transformBufferOnSave(void* src, void* dest, uint16_t src_bitspersample, uint16_t dest_bitspersample, uint16_t nsamples, size_t npixels, gboolean planar, gboolean linear);
BYTE uchar_pixel_icc_tx(BYTE in, int channel, int nchans, cmsHTRANSFORM display_transform);
WORD ushrt_pixel_icc_tx(WORD in, int channel, int nchans, cmsHTRANSFORM display_transform);
float float_pixel_icc_tx(float in, int channel, int nchans, cmsHTRANSFORM display_transform);
void fits_initialize_icc(fits *fit, cmsUInt8Number* EmbedBuffer, cmsUInt32Number EmbedLen);
void fits_check_icc(fits *fit);
cmsHPROFILE copyICCProfile(cmsHPROFILE profile);
unsigned char* get_icc_profile_data(cmsHPROFILE *profile, guint32 *len);
#endif /* SRC_CORE_ICC_PROFILE_H_ */
