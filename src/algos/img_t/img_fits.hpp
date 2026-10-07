/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

/*
 * Adapters between Siril's `fits` struct and the shared img_t<T> type.
 *
 * Both Siril fits and img_t are PLANAR (channel-major), so the conversion is a
 * straight per-element copy. A genuine zero-copy "borrow" is intentionally NOT
 * offered: USHORT data needs ushort<->float conversion, img_t owns fftw-aligned
 * storage that the FFT ops depend on, and the fits may need its buffer
 * reallocated - so a copy is inherent.
 *
 * from_fits() is a normalising reader for the [0,1] float convention used by
 * the denoisers. to_fits() is a general type-aware writer (see below).
 */

#pragma once

#include <algorithm>
#include <cstdlib>
#include <type_traits>

#include "image.hpp"  // img_t + core/siril.h (fits, data_type, WORD, USHRT_MAX_SINGLE,
                      // BYTE_IMG/USHORT_IMG/FLOAT_IMG, PRINT_ALLOC_ERR)

// fit_replace_buffer is defined in image_format_fits.c (C linkage); it swaps in
// a new data buffer, frees the old one, NULLs the now-unused buffer, rewires the
// per-channel pdata/fpdata pointers and sets type/bitpix/orig_bitpix.
extern "C" void fit_replace_buffer(fits *fit, void *newbuf, data_type newtype);

namespace imgops {

//! Build an owning planar img_t<float> from a Siril fits. USHORT data is
//! normalised to [0,1]; with clip=true the result is clamped to [0,1] (matching
//! the denoise input sanitisation).
inline img_t<float> from_fits(const fits* fit, bool clip = true) {
    const int w = fit->naxes[0];
    const int h = fit->naxes[1];
    const int d = fit->naxes[2];
    img_t<float> out(w, h, d);
    const size_t n = static_cast<size_t>(w) * h * d;
    if (fit->type == DATA_FLOAT) {
        if (clip)
            for (size_t i = 0; i < n; i++)
                out.data[i] = std::max(std::min(1.f, fit->fdata[i]), 0.f);
        else
            for (size_t i = 0; i < n; i++)
                out.data[i] = fit->fdata[i];
    } else {
        const float invnorm = 1.f / USHRT_MAX_SINGLE;
        if (clip)
            for (size_t i = 0; i < n; i++)
                out.data[i] = std::max(std::min(1.f, static_cast<float>(fit->data[i]) * invnorm), 0.f);
        else
            for (size_t i = 0; i < n; i++)
                out.data[i] = static_cast<float>(fit->data[i]) * invnorm;
    }
    return out;
}

//! General-purpose img_t<T> -> fits writer. Updates an EXISTING fits IN PLACE
//! (it never allocates a new fits, so caller-side metadata - keywords, header,
//! WCS, ... - is preserved) to hold `im`, reallocating the pixel buffer and
//! fixing up rx/ry/naxes/naxis/bitpix/type when the size, channel count or bit
//! depth differs. The element type chooses the fits representation:
//!   - float                 -> DATA_FLOAT  (bitpix FLOAT_IMG)
//!   - double                -> DATA_FLOAT, narrowed to float
//!   - uint16_t / WORD        -> DATA_USHORT (bitpix USHORT_IMG)
//!   - uint8_t                -> DATA_USHORT (value cast to uint16), orig_bitpix
//!                              set to BYTE_IMG to record the 8-bit origin
//! This is a plain type conversion (no [0,1] range scaling); callers that work
//! in the normalised [0,1] convention quantise to the desired type first.
template <typename T>
void to_fits(const img_t<T>& im, fits* fit) {
    static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value ||
                  std::is_same<T, uint8_t>::value || std::is_same<T, uint16_t>::value,
                  "to_fits supports img_t of float, double, uint8_t or uint16_t");
    constexpr bool floaty = std::is_floating_point<T>::value;
    const size_t nbdata = static_cast<size_t>(im.w) * im.h * im.d;

    void* buf = malloc(nbdata * (floaty ? sizeof(float) : sizeof(WORD)));
    if (!buf) {
        PRINT_ALLOC_ERR;
        return;
    }
    if constexpr (floaty) {
        float* f = static_cast<float*>(buf);
        for (size_t i = 0; i < nbdata; i++)
            f[i] = static_cast<float>(im.data[i]);   // narrows double -> float
    } else {
        WORD* w = static_cast<WORD*>(buf);
        for (size_t i = 0; i < nbdata; i++)
            w[i] = static_cast<WORD>(im.data[i]);     // casts uint8 -> uint16
    }

    // fit_replace_buffer reads naxes for the per-channel pointer offsets, so set
    // the shape first.
    fit->rx = fit->naxes[0] = im.w;
    fit->ry = fit->naxes[1] = im.h;
    fit->naxes[2] = im.d;
    fit->naxis = (im.d == 1) ? 2 : 3;
    fit_replace_buffer(fit, buf, floaty ? DATA_FLOAT : DATA_USHORT);

    if constexpr (std::is_same<T, uint8_t>::value)
        fit->orig_bitpix = BYTE_IMG;  // record the 8-bit origin (data lives in a uint16 buffer)
}

}  // namespace imgops
