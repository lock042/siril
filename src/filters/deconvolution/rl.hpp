#pragma once

#include "image.hpp"
#include "vec2.hpp"
#include "optimization.hpp"
#include "fft.hpp"
#include "utils.hpp"
#include "chelperfuncs.h"

namespace richardsonlucy {

    template <typename T>
    void rl_deconvolve(img_t<T>& x, const img_t<T>& f, const img_t<T>& K,
                                T lambda, int maxiter) {

        x = f;
        // Generate OTF of kernel
        img_t<std::complex<T>> K_otf(f.w, f.h, f.d);
        K_otf.padcirc(K);
        K_otf.map(K_otf * std::complex<T>(K.d) / K.sum());
        K_otf.fft(K_otf);

        // Flip K and generate OTF
        img_t<T> Kf(K.w, K.h, K.d);
        Kf.flip(K);
        img_t<std::complex<T>> Kflip_otf(f.w, f.h, f.d);
        Kflip_otf.padcirc(Kf);
        Kflip_otf.map(Kflip_otf * std::complex<T>(Kf.d) / Kf.sum());
        Kflip_otf.fft(Kflip_otf);
        img_t<std::complex<T>> est(f.w, f.h, f.d);
        est.map(f);
        img_t<std::complex<T>> working(f.w, f.h, f.d);
        working.map(est);
        float regularization = 0.f;
        for (int iter = 0 ; iter < maxiter ; iter++) {
            printf("Iteration %d ", iter);
            working.fft(est);
            working.map(working * K_otf);
            working.ifft(working);
            working.map(f / working);
            working.fft(working);
            working.map(working * Kflip_otf);
            working.ifft(working);
            est.map(working * (est / (1 - lambda * regularization)));
            printf("complete...\n");
            updateprogress("Richardson-Lucy deconvolution", (iter / maxiter));
        }
        x.map(std::real(est)); // x needs to be real
    }
};

