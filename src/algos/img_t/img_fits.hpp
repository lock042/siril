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
 * Adapters between Siril's `fits` struct and the shared img_t<float> type.
 *
 * Both Siril fits and img_t are PLANAR (channel-major), so the conversion is a
 * straight per-element copy. A genuine zero-copy "borrow" is intentionally NOT
 * offered: USHORT data needs ushort->float conversion, float data needs
 * clipping to [0,1] (which must not mutate the caller's fits), and img_t owns
 * fftw-aligned storage that the FFT-based operations depend on - so a copy is
 * inherent. These helpers replace the ad-hoc raw-buffer juggling that the
 * denoise/deconvolution entry points previously did by hand.
 */

#pragma once

#include <algorithm>

#include "image.hpp"  // img_t + core/siril.h (fits, DATA_FLOAT, USHRT_MAX_SINGLE)
#include "core/proto.h"  // roundf_to_WORD (self-guards with extern "C")

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

//! Write an img_t<float> result back into an existing fits, IN PLACE (it does
//! not allocate a new fits). The shape metadata (rx, ry, naxes, naxis, bitpix)
//! is updated to match `im`; the pixel data is then written, blended with the
//! fit's existing contents by `modulation` (fit = (1-modulation)*fit +
//! modulation*im), with USHORT de-normalised.
//!
//! Precondition: the fit's data buffers (fdata/data and the fpdata/pdata
//! per-channel pointers) must already be allocated for `im`'s dimensions - which
//! holds when `fit` is the same one `im` was built from via from_fits(). For a
//! differently-shaped target, (re)allocate the fit (e.g. with new_fit_image)
//! first. modulation < 1 blends over the fit's existing data, so it requires
//! that original data to be present.
inline void to_fits(const img_t<float>& im, fits* fit, float modulation = 1.f) {
    fit->rx = fit->naxes[0] = im.w;
    fit->ry = fit->naxes[1] = im.h;
    fit->naxes[2] = im.d;
    fit->naxis = (im.d == 1) ? 2 : 3;
    fit->bitpix = (fit->type == DATA_FLOAT) ? FLOAT_IMG : USHORT_IMG;

    const size_t n = static_cast<size_t>(im.w) * im.h * im.d;
    if (fit->type == DATA_FLOAT) {
        for (size_t i = 0; i < n; i++)
            fit->fdata[i] = (1.f - modulation) * fit->fdata[i] + modulation * im.data[i];
    } else {
        const float invnorm = 1.f / USHRT_MAX_SINGLE;
        for (size_t i = 0; i < n; i++)
            // roundf_to_WORD (float) is the type-appropriate conversion for our
            // float result, matching Siril's float_to_ushort_range; the old
            // denoise code used round_to_WORD (double), which needlessly
            // promoted the float and could differ by 1 LSB near saturation.
            fit->data[i] = roundf_to_WORD(USHRT_MAX_SINGLE * ((1.f - modulation) * (fit->data[i] * invnorm)
                                                             + modulation * im.data[i]));
    }
}

}  // namespace imgops
