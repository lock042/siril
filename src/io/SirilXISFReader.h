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
#ifndef SRC_IO_SIRILXISFREADER_H_
#define SRC_IO_SIRILXISFREADER_H_

extern "C" {
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
}

#ifdef HAVE_LIBXISF

#include "libxisf.h"

using namespace LibXISF;


class SirilXISFReader {
public:
	SirilXISFReader();  // Default constructor
    explicit SirilXISFReader(const std::string& filename);  // Constructor with filename
    ~SirilXISFReader();  // Destructor


private:
    //std::string filename_;

    uint64_t m_width { 0 };
    uint64_t m_height { 0 };
    uint64_t m_samples_per_channel { 0 };
    uint64_t m_channels { 0 };
    uint8_t *m_ImageBuffer { nullptr };
    uint32_t m_ImageBufferSize { 0 };
};

#endif /* HAVE_LIBXISF */

#endif /* SRC_IO_SIRILXISFREADER_H_ */
