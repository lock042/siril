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
#pragma once

#include <cstdio>
#include "image.hpp"
#include "vec2.hpp"
#include "optimization.hpp"
#include "fft.hpp"
#include "utils.hpp"
#include "chelperfuncs.h"
#include "rlstrings.h"

namespace wiener {
    template <typename T>
    void wiener_deconvolve(img_t<T>& x, const img_t<T>& f, const img_t<T>& K, T sigma) {

        // sigma is the noise std deviation. This assumes the noise is Gaussian i.e. its power spectrum is constant.
        // To extend this to other noise models, if needed, replace T sigma by const img_t<std::complex<T>> sigma
        // where the img_t contains the 2D noise power spectrum.

        assert(K.w % 2); // kernel must be odd in both dimensions
        assert(K.h % 2);
        x = f;

        // Initialize img_ts
        img_t<std::complex<T>> H(f.w, f.h, f.d);
        img_t<std::complex<T>> denom(f.w, f.h, f.d);
        img_t<std::complex<T>> G(f.w, f.h, f.d);

        // Generate OTF of kernel
        H.padcirc(K);
        H.map(H * std::complex<T>(K.d) / K.sum());
        H.fft(H);

        // Generate |H^2| = H * complex conjugate of H
        denom.map(std::conj(H) * H);
        denom.map(denom + sigma);
        denom.sanitize(); // Avoid NaNs and zeros in the denominator
        updateprogress(msg_wiener, 0.33);

        // Take the FFT of the image f, call this G
        G.map(f);
        G.fft(G);
        updateprogress(msg_wiener, 0.66);

        // Apply the Wiener filter
        G.map((G * std::conj(H)) / denom);
        G.ifft(G);
        updateprogress(msg_wiener, 1.0);

        x.map(std::real(G));
    }

}
