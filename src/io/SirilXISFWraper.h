/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2020 team siril_free-astro (see more in AUTHORS file)
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
#ifndef SRC_IO_SIRILXISFWRAPER_H_
#define SRC_IO_SIRILXISFWRAPER_H_

#ifdef __cplusplus
extern "C" {
#endif
#include "config.h"

#include <gdk-pixbuf/gdk-pixbuf.h>
#include <glib.h>

#ifdef HAVE_LIBXISF

struct xisf_data {
	uint8_t *data;
    uint64_t width;
    uint64_t height;
    uint64_t channelCount;
    uint64_t sampleFormat;

    char *fitsHeader;
	uint8_t *icc_buffer;
	uint32_t icc_length;
};

int siril_get_xisf_buffer(const char *filename, struct xisf_data *xdata);
GdkPixbuf* get_thumbnail_from_xisf(char *filename, gchar **descr);

#ifdef __cplusplus
}
#endif

#endif /* HAVE_LIBXISF */

#endif /* SRC_IO_SIRILXISFWRAPER_H_ */
