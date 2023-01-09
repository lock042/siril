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

        // Generate OTF of kernel
        img_t<std::complex<T>> H(f.w, f.h, f.d);
        H.padcirc(K);
        H.map(H * std::complex<T>(K.d) / K.sum());
        H.fft(H);

        // Generate |H^2| = H * complex conjugate of H
        img_t<std::complex<T>> H_abs_sq(H.w, H.h, H.d);
        H_abs_sq.map(std::conj(H) * H);
        img_t<std::complex<T>> denom(f.w, f.h, f.d);
        denom.map(H_abs_sq + sigma);
        denom.sanitize(); // Avoid NaNs and zeros in the denominator

        // Take the FFT of the image f, call this G
        img_t<std::complex<T>> G(f.w, f.h, f.d);
        G.map(f);
        G.fft(G);

        // Apply the Wiener filter
        G.map((G * std::conj(H)) / denom);
        G.ifft(G);
        x.map(std::real(G));
    }
}
