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


#include <iostream>
#include <iomanip>      // std::setw
#include <string>
#include <fitsio.h>	// fitsfile

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
			return false;
		}

		std::ostringstream fitsHeaderStream;

		fitsHeaderStream << "SIMPLE  =                    T / file does conform to FITS standard" << std::endl;

        auto &fitsKeywords = image.fitsKeywords();
        for(auto &fitsKeyword : fitsKeywords) {
        	if (fitsKeyword.name == "SIMPLE") continue;
        	if (fitsKeyword.name == "COMMENT") continue;
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

		uint8_t *buffer = nullptr;

		if (image.pixelStorage() == LibXISF::Image::Normal) {
			LibXISF::Image normalImage = image;
			if (normalImage.pixelStorage() == LibXISF::Image::Normal) {
				normalImage.convertPixelStorageTo(LibXISF::Image::Planar);
			}
			buffer = (uint8_t *)normalImage.imageData();
		} else {
			buffer = (uint8_t *)image.imageData();
		}
		xdata->size = image.imageDataSize();

		xdata->data = (uint8_t *) malloc(xdata->size);
		if (!xdata->data) {
			return 1;
		}

		memcpy(xdata->data, buffer, xdata->size);

		xisfReader.close();

	} catch (const LibXISF::Error &error) {
		std::cout << error.what() << std::endl;
		return 1;
	}
	return 0;
}
