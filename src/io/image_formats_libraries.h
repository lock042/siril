/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_IO_IMAGE_FORMATS_LIBRARIES_H_
#define SRC_IO_IMAGE_FORMATS_LIBRARIES_H_

#include "core/siril.h"

#ifdef HAVE_LIBTIFF
#include <tiffio.h>
#endif
#ifdef HAVE_LIBJPEG
#include <jpeglib.h>
#endif
#ifdef HAVE_LIBPNG
#include <png.h>
#include <setjmp.h>
#endif
#ifdef HAVE_LIBRAW
#include <libraw/libraw.h>
#include <libraw/libraw_version.h>
#endif
#ifdef HAVE_LIBHEIF
#include <libheif/heif.h>
#endif

/****************** image_formats_libraries.h ******************/
#ifdef HAVE_LIBTIFF
int readtif(const char *name, fits *fit, gboolean force_float);
int savetif(const char *name, fits *fit, uint16 bitspersample);
#endif

#ifdef HAVE_LIBJPEG
int readjpg(const char*, fits*);
int savejpg(const char*, fits*, int);
#endif

#ifdef HAVE_LIBPNG
int readpng(const char*, fits*);
int savepng(const char *filename, fits *fit, guint32 bytes_per_sample,
		gboolean is_colour);
#endif

#ifdef HAVE_LIBRAW
int open_raw_files(const char*, fits*, gboolean);
#endif

#ifdef HAVE_LIBHEIF
int readheif(const char *name, fits *fit, gboolean interactive);
#endif

#endif /* SRC_IO_IMAGE_FORMATS_LIBRARIES_H_ */
