#pragma once

#include <cstdio>
#include "image.hpp"
#include "vec2.hpp"
#include "optimization.hpp"
#include "fft.hpp"
#include "utils.hpp"
#include "chelperfuncs.h"
#include "rlstrings.h"

namespace richardsonlucy {

    template <typename T>
    void rl_deconvolve(img_t<T>& x, const img_t<T>& f, const img_t<T>& K, T lambda, int maxiter, T stopcriterion, int regtype, float stepsize, int stopcriterion_active) {
        assert(K.w % 2);
        assert(K.h % 2);
        x = f;
        optimization::operators::gradient<T> gradient(f);
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
//        est.set_value(T(0.001));
        img_t<std::complex<T>> ratio(f.w, f.h, f.d);
        ratio.map(est);
        float reallambda = 1.f / lambda; // For consistency with other algorithms
        img_t<T> gxx(f.w, f.h, f.d);
        img_t<T> gxy(f.w, f.h, f.d);
        img_t<T> gyy(f.w, f.h, f.d);
        for (int iter = 0 ; iter < maxiter ; iter++) {
            if (is_thread_stopped())
                continue;
            w.map(std::real(est));
            if (regtype == 0 || regtype == 3) {
                // Calculate TV regularization weighting
                gxx.gradientx(w); // Use gxx, the name isn't quite appropriate but it
                                 // saves having to make img_ts within the loop
                gxx.map(std::max(1.e-9f, gxx)); // Avoid div/0
                for (int i = 0 ; i < gxx.size; i++)
                    gxx[i] = std::max(1.e-9f, gxx[i]);
                gyy.gradienty(w);
                for (int i = 0 ; i < gyy.size; i++)
                    gyy[i] = std::max(1.e-9f, gyy[i]); // Avoid div/0
                w.map(std::hypot(gxx, gyy)); // |grad(w)|
                gxx.map(gxx / w); // Together these 2 lines make gx, gy hold the
                gyy.map(gyy / w); // components of grad(est)
                w.divergence(gxx, gyy); // w now holds div(grad(est) / |grad(est)|)
            } else if (regtype == 1 || regtype == 4) {
                // Calculate Frobenius-Hessian weighting
                gxx.gradientxx(w);
                gxx.map(std::max(1.e-9f, gxx)); // Avoid div/0
                gxx.map(gxx * gxx);
                gxy.gradientxy(w);
                gxy.map(std::max(1.e-9f, gxy));
                gxy.map(gxy * gxy);
                gxy.map(gxy * T(2));
                gyy.gradientyy(w);
                gyy.map(std::max(1.e-9f, gyy));
                gyy.map(gyy * gyy);
                w.map(gxy + gyy);
                w.map(gxx + w);
                w.map(std::pow(w, T(0.5)));
            }

            // Richardson-Lucy iteration
            ratio.fft(est);
            ratio.map(ratio * K_otf); // convolve
            ratio.ifft(ratio); // denominator
            ratio.map(f / ratio); // divide
            ratio.fft(ratio);
            ratio.map(ratio * Kflip_otf); // correlate (convolve with flip)
            ratio.ifft(ratio);
            gxy.map(std::real(est)); // From here on, gxy is used for the stopping criterion
            T dt = T(stepsize);
            switch (regtype) {
                case 5: // 5 and 4 are multiplicative RL with FH and TV reg
                    est.map(ratio * est); // Basic multiplicative R-L: multiply old estimate by ratio and regularization factor to get new estimate
                    dt = T(1.);
                    break;
                case 4:
                case 3:
                    est.map(ratio * est * (T(1) / (T(1) - reallambda * w))); // Multiply old estimate by ratio and regularization factor to get new estimate
                    dt = T(1.);
                    break;
                case 2:
                    est.map(est + dt * (T(-1) + ratio)); // Basic additive gradient-descent form, no regularization
                    break;
                case 1: // FH, additive gradient descent
                case 0:
                    est.map(est + dt * (T(-1) + (reallambda * w) + ratio)); // TV, additive gradient-descent form
                    break;
                default:
                    break;
            }
            if (sequence_is_running == 0)
                updateprogress(msg_rl, (static_cast<float>(iter + 1) / static_cast<float>(maxiter)));
            if (stopcriterion_active == 1) {
                // Stopping criterion?
                gxy.map((std::abs(std::real(est) - gxy)) / std::abs(gxy));
                T stopping = gxy.sum() / gxy.size;
//                printf("stopping: %f\n", stopping);
                if (stopping < stopcriterion) {
                    char msg[100];
                    sprintf(msg, "%s %d\n", msg_earlystop, iter+1);
                    sirillog(msg);
                    iter = maxiter;
                }
            }
        }
        x.map(std::real(est)); // x needs to be real
    }

    template <typename T>
    void rl_deconvolve_naive(img_t<T>& x, const img_t<T>& f, const img_t<T>& K, T lambda, int maxiter, T stopcriterion, int regtype, float stepsize, int stopcriterion_active) {
        assert(K.w % 2);
        assert(K.h % 2);
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
            if (is_thread_stopped())
                continue;
            // Regularization calcs

            w.map(x);
            if (regtype == 0 || regtype == 3) {
                // Calculate TV regularization weighting
                gxx.gradientx(w); // Use gxx, the name isn't quite appropriate but it
                                 // saves having to make img_ts within the loop
                gxx.map(std::max(1.e-9f, gxx)); // Avoid div/0
                for (int i = 0 ; i < gxx.size; i++)
                    gxx[i] = std::max(1.e-9f, gxx[i]);
                gyy.gradienty(w);
                for (int i = 0 ; i < gyy.size; i++)
                    gyy[i] = std::max(1.e-9f, gyy[i]); // Avoid div/0
                w.map(std::hypot(gxx, gyy)); // |grad(w)|
                gxx.map(gxx / w); // Together these 2 lines make gx, gy hold the
                gyy.map(gyy / w); // components of grad(est)
                w.divergence(gxx, gyy); // w now holds div(grad(est) / |grad(est)|)
            } else if (regtype == 1 || regtype == 4) {
                // Calculate Frobenius-Hessian weighting
                gxx.gradientxx(w);
                gxx.map(std::max(1.e-9f, gxx)); // Avoid div/0
                gxx.map(gxx * gxx);
                gxy.gradientxy(w);
                gxy.map(std::max(1.e-9f, gxy));
                gxy.map(gxy * gxy);
                gxy.map(gxy * T(2));
                gyy.gradientyy(w);
                gyy.map(std::max(1.e-9f, gyy));
                gyy.map(gyy * gyy);
                w.map(gxy + gyy);
                w.map(gxx + w);
                w.map(std::pow(w, T(0.5)));
            }
            // Richardson-Lucy iteration
            ratio.conv2(x, K); // convolve with kernel to get denominator
            ratio.map(f / ratio); // divide f by denominator
            ratio.conv2(ratio, Kf); // convolve by flipped kernel
            x.map(ratio * x); // multiply up
            T dt = T(stepsize);
            switch (regtype) {
                case 5: // 5 and 4 are multiplicative RL with FH and TV reg
                    x.map(ratio * x); // Basic multiplicative R-L: multiply old estimate by ratio and regularization factor to get new estimate
                    break;
                case 4:
                case 3:
                    x.map(ratio * x * (T(1) / (T(1) - reallambda * w))); // Multiply old estimate by ratio and regularization factor to get new estimate
                    break;
                case 2:
                    x.map(x + dt * (T(-1) + ratio)); // Basic additive gradient-descent form, no regularization
                    break;
                case 1: // FH, additive gradient descent
                case 0:
                    x.map(x + dt * (T(-1) + (reallambda * w) + ratio)); // TV, additive gradient-descent form
                    break;
                default:
                    break;
            }
            if (sequence_is_running == 0)
                updateprogress(msg_rl, (static_cast<float>(iter + 1) / static_cast<float>(maxiter)));
            if (stopcriterion_active == 1) {
                // Stopping criterion?
                gxy.map((std::abs(x - gxy)) / std::abs(gxy));
                T stopping = gxy.sum() / gxy.size;
//                printf("stopping: %f\n", stopping);
                if (stopping < stopcriterion) {
                    char msg[100];
                    sprintf(msg, "%s %d\n", msg_earlystop, iter+1);
                    sirillog(msg);
                    iter = maxiter;
                }
            }
        }
    }

/*
        template <typename T>
    img_t<T> u_factor(img_t<T> lucy, img_t<T> raw_image, T multiplier=T(1)) {
        lucy.replacenan(T(0));

        double doublenoise;
        updatenoise(raw_image.data.data(), raw_image.w, raw_image.h, raw_image.d, &doublenoise);
        T noise_T = static_cast<T>(doublenoise);

        T first = T(-2) / (noise_T * noise_T);

        img_t<T> ratio(lucy.w, lucy.h, lucy.d);
        ratio.map(lucy / raw_image);
        printf("ratio mean: %f\n", ratio.mean());
        ratio.replacenan(T(0));
        ratio.map(std::log(ratio));
        ratio.replacenan(T(0));

        img_t<T> second(lucy.w, lucy.h, lucy.d);
        second.map(raw_image * ratio - lucy + raw_image);

        img_t<T> factor(lucy.w, lucy.h, lucy.d);
        factor.map(first * second);
        factor.replacenan(T(0));
        factor.map(std::min(factor, T(1)));

        return factor;
    }

    template <typename T>
    img_t<T> dampen(img_t<T> lucy, img_t<T> raw, int N=3, T noise_T=T(0), T multiplier=T(1)) {
        multiplier = T(1);
        img_t<T> first = u_factor(lucy, raw, multiplier);
        first.map(std::pow(first,T(N-1)));
        first.replacenan(T(0));
        printf("first.mean: %f\n", first.mean());

        img_t<T> second = u_factor(lucy, raw, multiplier);
        second.map(T(N) - T(N-1) * second);
        second.replacenan(T(0));
        printf ("second.mean: %f\n", second.mean());

        img_t<T> third(lucy.w, lucy.h, lucy.d);
        third.map((raw - lucy) / lucy);
        third.replacenan(T(0));
        printf ("third.mean: %f\n", third.mean());
        third.map(T(1) + first * second * third);
        return third;
    }

    template <typename T>
    void damped_rl_deconvolve(img_t<T>& x, const img_t<T>& f, const img_t<T>& K, T lambda, int maxiter, T stopcriterion) {

        assert(K.w % 2);
        assert(K.h % 2);
        int N = 10;
        T multiplier = T(1);
        T noise_T = T(0);
        img_t<T> f_init(f);
        img_t<T> lucy(f.w, f.h, f.d);
        lucy.set_value(f_init.mean()); // Start with grey image estimate
        x = f;
        optimization::operators::gradient<T> gradient(f);
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
        // Initial estimate is the observed image
        est.map(lucy);
        img_t<std::complex<T>> temp_lucy(f.w, f.h, f.d);
        img_t<std::complex<T>> ratio;
        for (int iter = 0 ; iter < maxiter ; iter++) {
            // Breakout if the stop button has been pressed
            if (is_thread_stopped())
                continue;

            // Richardson-Lucy iteration
            // -------------------------

            // Convolve the current estimate with the kernel
            temp_lucy.fft(est);
            temp_lucy.map(temp_lucy * K_otf);
            temp_lucy.ifft(temp_lucy);

            // Calculate damping ratio
            ratio = dampen(lucy, f, N, noise_T, multiplier);

            // Convolve damped ratio with the kernel
            img_t<std::complex<T>> top(f.w, f.h, f.d);
            {
                img_t<std::complex<T>> complexratio(f.w, f.h, f.d);
                complexratio.map(ratio);
                top.fft(complexratio);
            }
            top.map(top * K_otf);
            top.ifft(top);

            est.map(est * (top / temp_lucy));

            img_t<std::complex<T>> stop(f.w, f.h, f.d);
            stop.map(est);

            // New estimate is previous estimate * ratio
            est.map(est * ratio);

            // Stopping criterion
            stop.map((std::abs(std::real(est) - std::real(stop))) / std::abs(std::real(stop)));
            if (sequence_is_running == 0)
                updateprogress("Richardson-Lucy deconvolution...", (static_cast<float>(iter + 1) / static_cast<float>(maxiter)));
            T stopping = std::real(stop.sum()) / stop.size;
            if (stopping < stopcriterion) {
                sirillog("Richardson-Lucy halted early by the stopping criterion\n");
                iter = maxiter;
            }

        }


        x.map(std::real(est)); // x needs to be real
    }
*/

};
