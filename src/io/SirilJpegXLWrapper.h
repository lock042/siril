/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2025 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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

#ifndef SRC_IO_SIRILJXLREADER_H_
#define SRC_IO_SIRILJXLREADER_H_

#ifdef __cplusplus
extern "C" {
#endif
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_LIBJXL

#include <glib.h>
#include <gdk-pixbuf/gdk-pixbuf.h>

int DecodeJpegXlOneShotWrapper(const uint8_t* jxl, size_t size,
                         float** pixels, size_t* xsize,
                         size_t* ysize, size_t* zsize, size_t* extra_channels, uint8_t* bitdepth,
                         uint8_t** icc_profile, size_t *icc_profile_length,
                         uint8_t** internal_icc_profile, size_t *internal_icc_profile_length);

int EncodeJpegXlOneshotWrapper(const void* pixels, const uint32_t xsize,
						const uint32_t ysize, const uint32_t zsize, const uint8_t bitdepth,
						uint8_t** compressed, size_t* compressed_length, uint32_t effort, const double quality,
						uint8_t *icc_profile, uint32_t icc_profile_length);

GdkPixbuf* get_thumbnail_from_jxl(uint8_t *jxl, gchar **descr, size_t size);

#ifdef __cplusplus
}
#endif

#endif /* HAVE_LIBJXL */

#endif /* SRC_IO_SIRILJXLREADER_H_ */
