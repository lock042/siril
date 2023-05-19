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
#include <string>
#include <fitsio.h>	// fitsfile
#include <gdk-pixbuf/gdk-pixbuf.h>
#include <glib.h>

#include "io/SirilXISFWraper.h"
#include "io/SirilXISFReader.h"
#include "libxisf.h"


SirilXISFReader::SirilXISFReader() {
    // Default constructor implementation
}

int siril_get_xisf_buffer(const char *filename, struct xisf_data *xdata) {
	try {
		LibXISF::XISFReader xisfReader;
		xisfReader.open(filename);

		if (xisfReader.imagesCount() == 0) {
			xisfReader.close();
			return false;
		}

		const LibXISF::Image &image = xisfReader.getImage(0);

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
		default:
			xdata->sampleFormat = 0;
			xisfReader.close();
			return false;
		}

		std::ostringstream fitsHeaderStream;

		fitsHeaderStream << "SIMPLE  =                    T / file does conform to FITS standard" << std::endl;

		auto &fitsKeywords = image.fitsKeywords();
		for (auto &fitsKeyword : fitsKeywords) {
			if (fitsKeyword.name == "SIMPLE") continue;
			if (fitsKeyword.name == "COMMENT") continue;
			if (fitsKeyword.name == "END") continue;
			fitsHeaderStream << std::setw(8) << std::left << fitsKeyword.name;
			fitsHeaderStream << "= " << std::setw(20) << fitsKeyword.value;
			fitsHeaderStream << " / " << fitsKeyword.comment << std::endl;
		}
		fitsHeaderStream << "END";

		xdata->fitsHeader = strdup(fitsHeaderStream.str().c_str());
		std::cout << xdata->fitsHeader << std::endl;

		xdata->width = image.width();
		xdata->height = image.height();
		xdata->channelCount = image.channelCount();

		xdata->data = (uint8_t*) malloc(image.imageDataSize());
		if (!xdata->data) {
			xisfReader.close();
			return 1;
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
		return 1;
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
		return 0;
	}
}

GdkPixbuf* get_thumbnail_from_xisf(char *filename, gchar **descr) {
	GdkPixbuf *pixbuf = NULL;
	gchar *description = NULL;
	try {
		LibXISF::XISFReader xisfReader;
		xisfReader.open(filename);

		if (xisfReader.imagesCount() == 0) {
			xisfReader.close();
			return NULL;
		}

		const LibXISF::Image &thumbnail = xisfReader.getThumbnail();
		if (thumbnail.width() == 0 || thumbnail.height() == 0) {
			xisfReader.close();
			return NULL;
		}

		// TODO: We need to convert to RGB pixbuf when monochrome
		uint8_t *pixbuf_data = (uint8_t*) malloc(thumbnail.imageDataSize());
		if (!pixbuf_data) {
			xisfReader.close();
			return NULL;
		}
		LibXISF::Image planarThumbnail = thumbnail;
		planarThumbnail.convertPixelStorageTo(LibXISF::Image::Normal);
		memcpy(pixbuf_data, planarThumbnail.imageData(), planarThumbnail.imageDataSize());

		pixbuf = gdk_pixbuf_new_from_data(pixbuf_data,	// guchar* data
				GDK_COLORSPACE_RGB,	// only this supported
				FALSE,				// no alpha
				8,				// number of bits
				thumbnail.width(), thumbnail.height(),				// size
				thumbnail.width() * 3,				// line length in bytes
				(GdkPixbufDestroyNotify) free_preview_data, // function (*GdkPixbufDestroyNotify) (guchar *pixels, gpointer data);
				NULL);

		const LibXISF::Image &image = xisfReader.getImage(0);

		description = g_strdup_printf("%ld x %ld %s\n%ld %s (%d bits)", thumbnail.width(),
				thumbnail.height(), ngettext("pixel", "pixels", thumbnail.height()), thumbnail.channelCount(),
				ngettext("channel", "channels", thumbnail.channelCount()), get_bit_depth(image.sampleFormat()));
		xisfReader.close();
	} catch (const LibXISF::Error &error) {
		std::cout << error.what() << std::endl;
		return NULL;
	}
	*descr = description;
	return pixbuf;
}

#endif
