#pragma once

#include "image.hpp"
#include "vec2.hpp"
#include "optimization.hpp"
#include "fft.hpp"
#include "utils.hpp"
#include "chelperfuncs.h"

namespace richardsonlucy {

    template <typename T>
    void rl_deconvolve(img_t<T>& x, const img_t<T>& f, const img_t<T>& K, T lambda, int maxiter) {

        assert(K.w % 2);
        assert(K.h % 2);
        x = f;
        optimization::operators::gradient<T> gradient(f);
//        using gradtype = typename decltype(gradient)::out_type;
        img_t<T> w(f.w, f.h, f.d);

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
        float reallambda = 1.f / lambda; // For consistency with other algorithms
        img_t<T> gx(f.w, f.h, f.d);
        img_t<T> gy(f.w, f.h, f.d);
        for (int iter = 0 ; iter < maxiter ; iter++) {
            if (is_thread_stopped())
                continue;
            // Calculate TV reglarization weighting
            w.map(std::real(est));
            gx.gradientx(w);
            gx.map(std::max(1.e-6f, gx)); // Avoid div/0
            for (int i = 0 ; i < gx.size; i++)
                gx[i] = std::max(1.e-6f, gx[i]);
            gy.gradienty(w);
            for (int i = 0 ; i < gy.size; i++)
                gy[i] = std::max(1.e-6f, gy[i]); // Avoid div/0
            w.map(std::hypot(gx, gy));
            gx.map(gx / w);
            gy.map(gy / w);
            w.divergence(gx, gy);

            // Richardson-Lucy iteration
            working.fft(est);
            working.map(working * K_otf);
            working.ifft(working);
            working.map(f / working);
            working.fft(working);
            working.map(working * Kflip_otf);
            working.ifft(working);

            est.map(working * est * (T(1) / (T(1) - reallambda * w)));

            // Stopping criterion?

            printf("complete...\n");
            updateprogress("Richardson-Lucy deconvolution", (static_cast<float>(iter) / static_cast<float>(maxiter)));
        }
        x.map(std::real(est)); // x needs to be real
    }
};

