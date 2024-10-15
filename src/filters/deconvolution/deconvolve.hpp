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
#include "deconvolution.h"
#include "fft.hpp"
#include "utils.hpp"
#include "core/siril.h"
#include "core/processing.h" // for get_thread_run()
#include "core/siril_log.h" // for siril_log_message()
#include "gui/progress_and_log.h" // for set_progress_bar_data()

namespace deconvolve {
    template <typename T>
    void wiener_deconvolve(img_t<T>& x, const img_t<T>& f, const img_t<T>& K, T sigma) {
        double sliceprogress = (double) f.slices_complete / f.total_slices;
        // sigma is the noise modelling parameter. This assumes the noise is Gaussian i.e. its power spectrum is
        // constant. To extend this to other noise models, if needed, replace T sigma by const
        // img_t<std::complex<T>> sigma where the img_t contains the 2D noise power spectrum.

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
        denom.map((img::conj(H) * H) + sigma);
        denom.sanitize(); // Avoid NaNs and zeros in the denominator
        set_progress_bar_data(_("Wiener deconvolution..."), std::min(1.0, (sliceprogress + 0.33 / f.total_slices)));
        // The std::min() call is used because on completion rounding errors can take the value marginally above 1.0 and cause an assert fail in set_progress_bar_data()
        // Similarly for other calls to this function

        // Take the FFT of the image f, call this G
        G.map(f);
        G.fft(G);
        set_progress_bar_data(_("Wiener deconvolution..."), std::min(1.0, (sliceprogress + 0.67 / f.total_slices)));

        // Apply the Wiener filter
        G.map((G * img::conj(H)) / denom);
        G.ifft(G);
        set_progress_bar_data(_("Wiener deconvolution..."), std::min(1.0, (sliceprogress + 1.0 / f.total_slices)));

        x.map(img::real(G));
    }

    template <typename T>
    void rl_deconvolve_fft(img_t<T>& x, const img_t<T>& f, img_t<T>& K, T lambda, int maxiter, T stopcriterion, regtype_t regtype, float stepsize, int stopcriterion_active) {
        assert(K.w % 2);
        assert(K.h % 2);
        double sliceprogress = (double) f.slices_complete / f.total_slices;
        x = f;

        // Generate OTF of kernel
        img_t<std::complex<T>> K_otf(f.w, f.h, f.d);
        K_otf.padcirc(K);
        K_otf.map(K_otf * std::complex<T>(K.d) / K.sum());
        K_otf.fft(K_otf);

        // Flip K and generate OTF
        img_t<std::complex<T>> Kflip_otf(f.w, f.h, f.d);
        {
            K.flip();
            Kflip_otf.padcirc(K);
            Kflip_otf.map(Kflip_otf * std::complex<T>(K.d) / K.sum());
            Kflip_otf.fft(Kflip_otf);
        }
        img_t<std::complex<T>> est(f.w, f.h, f.d);
        est.map(f);
        float reallambda = 1.f / lambda; // For consistency with other algorithms
        for (int iter = 0 ; iter < maxiter ; iter++) {
            if (!get_thread_run())
                continue;
            img_t<T> w(f.w, f.h, f.d);
            w.map(img::real(est));
            if (regtype == 0 || regtype == 3) {
                // Calculate TV weighting
                auto dx = to_img(gradientx(w));
                auto dy = to_img(gradienty(w));
                auto eps = std::numeric_limits<T>::epsilon();
                auto mag = img::hypot(dx, dy) + eps;
                w.map(divergence((dx / mag), (dy / mag))); // w now holds div(grad(est) / |grad(est)|)
            } else if (regtype == 1 || regtype == 4) {
                // Calculate Frobenius-Hessian weighting
                auto gxx = gradientxx(w);
                auto gxy = gradientxy(w);
                auto gyy = gradientyy(w);
                auto xxsq = img::pow(gxx, T(2)); // Avoid div/0
                auto yysq = img::pow(gyy, T(2));
                auto sumsq = to_img(xxsq + yysq); // We have to evaluate here or the compiler can't cope
                                                  // with the complex unevaluated img_expr_t
                auto sumgrads = img::fma(T(2), img::pow(gxy, T(2)), sumsq);
                auto weight = img::sqrt(sumgrads);
                w.map(weight);
                w.sanitize();
            }
            // Richardson-Lucy iteration
            img_t<std::complex<T>> ratio(f.w, f.h, f.d);
            ratio.fft(est);
            ratio.map(ratio * K_otf); // convolve
            ratio.ifft(ratio); // denominator
            ratio.sanitize();
            ratio.map(f / ratio); // divide
            ratio.fft(ratio);
            ratio.map(ratio * Kflip_otf); // correlate (convolve with flip)
            ratio.ifft(ratio);
            img_t<T> stopcrit;
            if(stopcriterion_active == 1) {
                stopcrit.resize(f.w, f.h, f.d);
                stopcrit = to_img(img::real(est));
            }
            T dt = T(stepsize);
            switch (regtype) {
                case REG_NONE_MULT:
                    est.map(ratio * est); // Basic multiplicative R-L: multiply old estimate by ratio and regularization factor to get new estimate
                    dt = T(1.);
                    break;
                case REG_FH_MULT:
                case REG_TV_MULT:
                    est.map(ratio * est * (T(1) / (T(1) - reallambda * w))); // Multiply old estimate by ratio and regularization factor to get new estimate
                    dt = T(1.);
                    break;
                case REG_NONE_GRAD:
                    est.map(est + dt * (T(-1) + ratio)); // Basic additive gradient-descent form, no regularization
                    break;
                case REG_FH_GRAD:
                case REG_TV_GRAD:
                    est.map(est + dt * (T(-1) + (reallambda * w) + ratio));
                    break;
                default:
                    break;
            }
            if (sequence_is_running == 0)
                set_progress_bar_data(_("Richardson-Lucy deconvolution..."), std::min(1.0, sliceprogress + ((static_cast<float>(iter + 1) / static_cast<float>(maxiter))/f.total_slices)));
            if (stopcriterion_active == 1) {
                // Stopping criterion?
                auto stopmeasure = (img::abs(img::real(est) - stopcrit) / img::abs(stopcrit));
                T stopping = (to_img(stopmeasure)).sum() / stopmeasure.size;
                if (stopping < stopcriterion) {
                    char msg[100];
                    sprintf(msg, "%s %d\n", _("Richardson-Lucy halted early by the stopping criterion after iteration"), iter+1);
                    siril_log_message(msg);
                    iter = maxiter;
                }
            }
        }
        x.map(img::real(est)); // x needs to be real
    }

    template <typename T>
    void rl_deconvolve_naive(img_t<T>& x, const img_t<T>& f, const img_t<T>& K, T lambda, int maxiter, T stopcriterion, regtype_t regtype, float stepsize, int stopcriterion_active) {
        assert(K.w % 2);
        assert(K.h % 2);
        double sliceprogress = (double) f.slices_complete / f.total_slices;
        x = f;
        img_t<T> w(f.w, f.h, f.d);
        float reallambda = 1.f / lambda; // For consistency with other algorithms

        // Flip K and generate OTF
        img_t<T> Kf(K.w, K.h, K.d);
        Kf.flip(K);
        img_t<T> ratio(f.w, f.h, f.d);
        img_t<T> gxx(f.w, f.h, f.d);
        img_t<T> gxy(f.w, f.h, f.d);
        img_t<T> gyy(f.w, f.h, f.d);
        for (int iter = 0 ; iter < maxiter ; iter++) {
            if (!get_thread_run())
                continue;
            // Regularization calcs
            w.map(x);
            if (regtype == 0 || regtype == 3) {
                // Calculate TV regularization weighting
                gxx.gradientx(w); // Use gxx, the name isn't quite appropriate but it
                                 // saves having to make img_ts within the loop
                gxx.sanitize(); // Avoid div/0
                gyy.gradienty(w);
                gyy.sanitize(); // Avoid div/0
                w.map(img::hypot(gxx, gyy)); // |grad(w)|
                gxx.map(gxx / w); // Together these 2 lines make gx, gy hold the
                gyy.map(gyy / w); // components of grad(est)
                w.divergence(gxx, gyy); // w now holds div(grad(est) / |grad(est)|)
            } else if (regtype == 1 || regtype == 4) {
                // Calculate Frobenius-Hessian weighting
                gxx.gradientxx(w);
                auto sanitizexxsq = img::pow(img::max(1.e-9f, gxx), T(2)); // Avoid div/0
                gxy.gradientxy(w);
                auto twosanitizexysq = T(2) * (img::pow(img::max(1.e-9f, gxy), T(2)));
                gyy.gradientyy(w);
                auto sanitizeyysq = img::pow(img::max(1.e-9f, gyy), T(2));
                auto sumgrads = (sanitizexxsq + (twosanitizexysq + sanitizeyysq));
                w.map(img::pow(sumgrads, T(0.5)));
            }

            // Richardson-Lucy iteration
            ratio.conv2(x, K); // convolve with kernel to get denominator
            ratio.map(f / ratio); // divide f by denominator
            ratio.map(img::max(T(1.e-9), ratio));
            ratio.conv2(ratio, Kf); // convolve by flipped kernel
            T dt = T(stepsize);
            switch (regtype) {
                case REG_NONE_MULT: // 5 and 4 are multiplicative RL with FH and TV reg
                    x.map(ratio * x); // Basic multiplicative R-L: multiply old estimate by ratio and regularization factor to get new estimate
                    break;
                case REG_FH_MULT:
                case REG_TV_MULT:
                    x.map(ratio * x * (T(1) / (T(1) - reallambda * w))); // Multiply old estimate by ratio and regularization factor to get new estimate
                    break;
                case REG_NONE_GRAD:
                    x.map(x + dt * (T(-1) + ratio)); // Basic additive gradient-descent form, no regularization
                    break;
                case REG_FH_GRAD: // FH, additive gradient descent
                case REG_TV_GRAD:
                    x.map(x + dt * (T(-1) + (reallambda * w) + ratio)); // TV, additive gradient-descent form
                    break;
                default:
                    break;
            }
            if (sequence_is_running == 0)
                set_progress_bar_data(_("Richardson-Lucy deconvolution..."), std::min(1.0, sliceprogress + ((static_cast<float>(iter + 1) / static_cast<float>(maxiter))/f.total_slices)));
            if (stopcriterion_active == 1) {
                gxy.map((img::abs(x - gxy)) / img::abs(gxy));
                T stopping = gxy.sum() / gxy.size;
                if (stopping < stopcriterion) {
                    char msg[100];
                    sprintf(msg, "%s %d\n", _("Richardson-Lucy halted early by the stopping criterion after iteration"), iter+1);
                    siril_log_message(msg);
                    iter = maxiter;
                }
            }
        }
    }

    template <typename T>
    void sb_deconvolve(img_t<T>& x, const img_t<T>& f, const img_t<T>& K,
                            T lambda, T b_0, T b_factor, T max_beta, const int iters) {
        assert(K.w % 2);
        assert(K.h % 2);
        double sliceprogress = (double) f.slices_complete / f.total_slices;
        optimization::operators::gradient<T> gradient(f);
        using gradtype = typename decltype(gradient)::out_type;
        img_t<gradtype> w(f.w, f.h, f.d);
        img_t<std::complex<T>> f_ft(f.w, f.h, f.d);
        f_ft.map(f);
        f_ft.fft(f_ft);

        if (!get_thread_run())
            return;

        img_t<std::complex<T>> dx_otf(f.w, f.h, f.d);
        {
            img_t<T> dx(3, 3, f.d);
            for (int c = 0; c < f.d; c++) {
                dx(0, 1, c) = 0; dx(1, 1, c) = -1; dx(2, 1, c) = 1;
            }
            dx_otf.padcirc(dx);
        }
        dx_otf.fft(dx_otf);

        if (!get_thread_run())
            return;

        img_t<std::complex<T>> dy_otf(f.w, f.h, f.d);
        {
            img_t<T> dy(3, 3, f.d);
            for (int c = 0; c < f.d; c++) {
                dy(1, 0, c) = 0; dy(1, 1, c) = -1; dy(1, 2, c) = 1;
            }
            dy_otf.padcirc(dy);
        }
        dy_otf.fft(dy_otf);

        if (!get_thread_run())
            return;

        img_t<std::complex<T>> K_otf(f.w, f.h, f.d);
        K_otf.padcirc(K);
        K_otf.map(K_otf * std::complex<T>(K.d) / K.sum());
        K_otf.fft(K_otf);

        if (!get_thread_run())
            return;

        auto Kf = img::conj(K_otf) * f_ft;

        img_t<T> ksq(f.w, f.h, f.d);
        ksq.map(img::pow(img::abs(K_otf), T(2)));

        auto dxdysq = (img::pow(img::abs(dx_otf), T(2)) + img::pow(img::abs(dy_otf), T(2)));

        T beta = b_0;
        x = f;
        img_t<std::complex<T>> w1_ft(f.w, f.h, f.d);
        img_t<std::complex<T>> x_ft(f.w, f.h, f.d);
        while (beta < max_beta) {
            if (!get_thread_run())
                break;
            T gamma = beta / lambda;
            auto denom = ksq + gamma * dxdysq;

            for (int inner = 0; inner < iters; inner++) {
                if (!get_thread_run())
                    break;
                auto grad = gradient.direct(x);
                w.map(grad * (T(1) - T(1) / (img::max(T(1), beta * img::hypot(grad)))));

                w1_ft.map(gradient.adjoint(w));
                w1_ft.fft(w1_ft);
                auto num = Kf + gamma * w1_ft;
                x_ft.map(num / denom);
                x_ft.ifft(x_ft);
                x.map(img::real(x_ft));
            }
            beta *= b_factor;
            if (sequence_is_running == 0)
                set_progress_bar_data("Split Bregman deconvolution...", std::min(1.0, sliceprogress + ((beta - b_0) / (max_beta - b_0))/f.total_slices));
        }
    }
}
