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

#define ICC_COPYRIGHT "Copyright (C) 2005-2011 Francois Meyer, (C) 2012-2023 team free-astro (website: https://free-astro.org/index.php.Siril/). This ICC profile is licensed under the GNU Public Licence, either version 3 of the License or (at your option) any later version (http://www.gnu.org/licenses/)"

#define GAMMA_LINEAR 1.0
#define GAMMA_DISPLAY 2.19921875
#define TYPE_RGB_FLT_PLANAR (FLOAT_SH(1)|COLORSPACE_SH(PT_RGB)|CHANNELS_SH(3)|BYTES_SH(4))|PLANAR_SH(1))

const unsigned char* get_sRGB_profile_data(guint32 *len);
const unsigned char* get_gray_profile_data(guint32 *len);

//unsigned char* get_profile_buf(cmsHPROFILE* profile, uint32_t* profile_len);

int load_display_icc_profile(const char* filename);
int load_proof_icc_profile(const char* filename);
void initialize_icc_profiles_paths();

#endif /* SRC_CORE_ICC_PROFILE_H_ */
