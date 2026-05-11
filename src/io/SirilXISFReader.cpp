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

extern "C" {
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
}


#ifdef HAVE_LIBXISF

#include <iostream>
#include <iomanip>      // std::setw
#include <sstream>
#include <string>
#include <fitsio.h>	// fitsfile
#include <math.h>
#include <libintl.h>
#include <gdk-pixbuf/gdk-pixbuf.h>
#include <cstdint>

#include "io/SirilXISFWraper.h"
#include "io/SirilXISFReader.h"
#include "libxisf.h"


int siril_get_xisf_buffer(const char *filename, struct xisf_data *xdata) {
	try {
		LibXISF::XISFReader xisfReader;
		xisfReader.open(LibXISF::String(filename));

		if (xisfReader.imagesCount() == 0) {
			xisfReader.close();
			return -1;
		}

		const LibXISF::Image &image = xisfReader.getImage(0);

		// Retrieve the ICC profile, if there is one
		const LibXISF::ByteArray profile = image.iccProfile();
		xdata->icc_length = profile.size();
		if (xdata->icc_length > 0) {
			xdata->icc_buffer = (uint8_t*) malloc(xdata->icc_length * sizeof(uint8_t));
			if (!xdata->icc_buffer)
				return -1;
			memcpy(xdata->icc_buffer, profile.data(), xdata->icc_length);
		}
		else
			xdata->icc_buffer = nullptr;

		switch (image.sampleFormat()) {
		case LibXISF::Image::UInt8:
			xdata->sampleFormat = BYTE_IMG;
			break;
		case LibXISF::Image::UInt16:
			xdata->sampleFormat = USHORT_IMG;
			break;
		case LibXISF::Image::UInt32:
			xdata->sampleFormat = LONG_IMG;
			break;
		case LibXISF::Image::Float32:
			xdata->sampleFormat = FLOAT_IMG;
			break;
		case LibXISF::Image::Float64:
			xdata->sampleFormat = DOUBLE_IMG;
			break;
		default:
			xdata->sampleFormat = 0;
			xisfReader.close();
			return -1;
		}

		std::ostringstream fitsHeaderStream;
		xdata->fitsHeader = NULL;

		const auto& fitsKeywords = image.fitsKeywords();
		if (!image.fitsKeywords().empty()) {
			fitsHeaderStream << "SIMPLE  =                    T / file does conform to FITS standard" << std::endl;
			for (const auto& fitsKeyword : fitsKeywords) {
				if (fitsKeyword.name == "SIMPLE") continue;
				if (fitsKeyword.name == "COMMENT") continue;
				if (fitsKeyword.name == "END") continue;
				if (fitsKeyword.name.rfind("PV", 0) == 0) continue; // remove PV keywords as they mess up the WCS readout
				fitsHeaderStream << std::setw(8) << std::left << fitsKeyword.name;
				fitsHeaderStream << "= " << std::setw(20) << fitsKeyword.value;
				fitsHeaderStream << " / " << fitsKeyword.comment << std::endl;
			}
			fitsHeaderStream << "END";
			xdata->fitsHeader = strdup(fitsHeaderStream.str().c_str());
		}

		xdata->width = image.width();
		xdata->height = image.height();
		xdata->channelCount = image.channelCount();

		xdata->data = (uint8_t*) malloc(image.imageDataSize());
		if (!xdata->data) {
			xisfReader.close();
			return -1;
		}

		if (image.pixelStorage() == LibXISF::Image::Normal) {
			LibXISF::Image normalImage = image;
			normalImage.convertPixelStorageTo(LibXISF::Image::Planar);
			memcpy(xdata->data, normalImage.imageData(),
					normalImage.imageDataSize());

		} else {
			memcpy(xdata->data, image.imageData(), image.imageDataSize());
		}

		xisfReader.close();

	} catch (const LibXISF::Error &error) {
		std::cout << error.what() << std::endl;
		return -1;
	}
	return 0;
}

static GdkPixbufDestroyNotify free_preview_data(guchar *pixels, gpointer data) {
	free(pixels);
	free(data);
	return FALSE;
}

static int get_bit_depth(LibXISF::Image::SampleFormat depth) {
	switch (depth) {
	case LibXISF::Image::UInt8:
		return 8;
	case LibXISF::Image::UInt16:
		return 16;
	case LibXISF::Image::UInt32:
	case LibXISF::Image::Float32:
	case LibXISF::Image::Complex32:
		return 32;
	case LibXISF::Image::UInt64:
	case LibXISF::Image::Float64:
	case LibXISF::Image::Complex64:
		return 64;
	default:
		return -1;
	}
}

/* Core thumbnail extractor: returns a malloc'd RGB888 byte buffer plus
 * dimensions, or NULL on error.  Caller owns *data (free with free()) and
 * *descr (free with g_free()). */
extern "C" guchar *extract_thumbnail_from_xisf(const char *filename, gchar **descr,
                                                int *width_out, int *height_out) {
	gchar *description = NULL;
	guchar *pixbuf_data = NULL;
	int out_w = 0, out_h = 0;
	try {
		LibXISF::XISFReader xisfReader;
		xisfReader.open(LibXISF::String(filename));

		if (xisfReader.imagesCount() == 0) {
			xisfReader.close();
			return NULL;
		}

		/* Get thumbnail if available */
		const LibXISF::Image &thumbnail = xisfReader.getThumbnail();
		if (thumbnail.width() != 0 && thumbnail.height() != 0) {

			/* Only RGB is handled. Monochrome thumbnails are expanded to RGB. */
			size_t extra_size = 0;
			if (thumbnail.channelCount() == 1) {
				extra_size = 2;
			}
			pixbuf_data = (guchar*) malloc(thumbnail.imageDataSize() + extra_size * thumbnail.imageDataSize());
			if (!pixbuf_data) {
				xisfReader.close();
				return NULL;
			}

			LibXISF::Image planarThumbnail = thumbnail;
			planarThumbnail.convertPixelStorageTo(LibXISF::Image::Normal);
			if (thumbnail.channelCount() == 1) {
				uint8_t *buffer = (uint8_t *) planarThumbnail.imageData();
				for (size_t i = 0, j = 0; i < planarThumbnail.imageDataSize() * 3; i += 3, j++) {
					pixbuf_data[i + 0] = buffer[j];
					pixbuf_data[i + 1] = buffer[j];
					pixbuf_data[i + 2] = buffer[j];
				}
			} else {
				memcpy(pixbuf_data, planarThumbnail.imageData(), planarThumbnail.imageDataSize());
			}
			out_w = (int)thumbnail.width();
			out_h = (int)thumbnail.height();
		}

		const LibXISF::Image &image = xisfReader.getImage(0);

		description = g_strdup_printf("%" G_GUINT64_FORMAT " x %" G_GUINT64_FORMAT " %s\n%" G_GUINT64_FORMAT " %s (%d bits)",
						image.width(), image.height(), ngettext("pixel", "pixels", image.height()), image.channelCount(),
						ngettext("channel", "channels", image.channelCount()), get_bit_depth(image.sampleFormat()));
		xisfReader.close();
	} catch (const LibXISF::Error &error) {
		std::cout << error.what() << std::endl;
		if (pixbuf_data) free(pixbuf_data);
		return NULL;
	}
	*descr = description;
	if (width_out) *width_out = out_w;
	if (height_out) *height_out = out_h;
	return pixbuf_data;
}

/* GdkPixbuf shim around extract_thumbnail_from_xisf, kept for the GTK3
 * build. */
GdkPixbuf* get_thumbnail_from_xisf(char *filename, gchar **descr) {
	int w = 0, h = 0;
	guchar *data = extract_thumbnail_from_xisf(filename, descr, &w, &h);
	if (!data) return NULL;
	if (w <= 0 || h <= 0) { free(data); return NULL; }
	return gdk_pixbuf_new_from_data(data, GDK_COLORSPACE_RGB, FALSE, 8,
			w, h, w * 3,
			(GdkPixbufDestroyNotify) free_preview_data, NULL);
}

#endif
