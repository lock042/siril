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

#include <iostream>
#include <random>

#include "deconvolve.hpp"
#include "utils.hpp"
#include "edgetaper.hpp"

extern "C" int wienerdec(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float sigma, int max_threads) {
    img_t<float>::use_threading(max_threads);
    img_t<float> u;
    for (unsigned c = 0 ; c < nchans ; c++) {
        img_t<float> f(rx, ry, 1, fdata + c * rx * ry);
        unsigned kc = (c < kchans ? c : 0);
        img_t<float> K(kernelsize, kernelsize, 1, kernel + kc * kernelsize * kernelsize);
        K.map(K / K.sum());
        float max = f.max();
        if (max == 0.0f)
            return 1;
        if (max != 1.0f)
            f.map(f / max);
        f = utils::add_padding(f, K);
        edgetaper(f, f, K, 3);

        deconvolve::wiener_deconvolve(u, f, K, sigma);

        u = utils::remove_padding(u, K);
        if (max != 1.0f)
            u.map(u * max);
        for (unsigned i = 0; i < rx * ry; i++)
            fdata[c * rx * ry + i] = u.data[i];
    }
    if (sequence_is_running == 0)
        updateprogress("Ready.", 0.);
    return 0;
}

extern "C" int fft_richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int maxiter, float stopcriterion, int max_threads, int regtype, float stepsize, int stopcriterion_active) {
    img_t<float>::use_threading(max_threads);
    img_t<float> u;
    for (unsigned c = 0 ; c < nchans ; c++) {
        img_t<float> f(rx, ry, 1, fdata + c * rx * ry);
        unsigned kc = (c < kchans ? c : 0);
        img_t<float> K(kernelsize, kernelsize, 1, kernel + kc * kernelsize * kernelsize);
        K.map(K / K.sum());
        float max = f.max();
        if (max == 0.0f)
            return 1;
        if (max != 1.0f)
            f.map(f / max);
        f = utils::add_padding(f, K);
        edgetaper(f, f, K, 3);

        deconvolve::rl_deconvolve_fft(u, f, K, 2.f / lambda, maxiter, stopcriterion, regtype, stepsize, stopcriterion_active);

        u = utils::remove_padding(u, K);
        if (max != 1.0f)
            u.map(u * max);
        for (unsigned i = 0; i < rx * ry; i++) {
            fdata[i + rx * ry * c] = u.data[i];
        }
    }
    if (sequence_is_running == 0)
        updateprogress("Ready.", 0.);
    return 0;
}

extern "C" int naive_richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int maxiter, float stopcriterion, int max_threads, int regtype, float stepsize, int stopcriterion_active) {
    img_t<float>::use_threading(max_threads);
    img_t<float> u;
    for (unsigned c = 0 ; c < nchans ; c++) {
        img_t<float> f(rx, ry, 1, fdata + c * rx * ry);
        unsigned kc = (c < kchans ? c : 0);
        img_t<float> K(kernelsize, kernelsize, 1, kernel + kc * kernelsize * kernelsize);
        K.map(K / K.sum());
        float max = f.max();
        if (max == 0.0f)
            return 1;
        if (max != 1.0f)
            f.map(f / max);
        f = utils::add_padding(f, 2*K.w, 2*K.w);
        edgetaper(f, f, K, 3);

        deconvolve::rl_deconvolve_naive(u, f, K, 2.f / lambda, maxiter, stopcriterion, regtype, stepsize, stopcriterion_active);

        u = utils::remove_padding(u, 2*K.w, 2*K.w);
        if (max != 1.0f)
            u.map(u * max);
        for (unsigned i = 0; i < rx * ry; i++) {
            fdata[i + c * rx * ry] = u.data[i];
        }
    }
    if (sequence_is_running == 0)
        updateprogress("Ready.", 0.);
    return 0;
}

extern "C" int split_bregman(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int iters, int max_threads) {
    img_t<float>::use_threading(max_threads);
    img_t<float> u;
    for (unsigned c = 0 ; c < nchans ; c++) {
        unsigned kc = (c < kchans ? c : 0);
        img_t<float> f(rx, ry, 1, fdata + c * rx * ry);
        img_t<float> K(kernelsize, kernelsize, 1, kernel + kc * kernelsize * kernelsize);
        K.map(K / K.sum());
        float max = f.max();
        if (max == 0.0f)
            return 1;
        if (max != 1.0f)
            f.map(f / max);
        f = utils::add_padding(f, K);
        edgetaper(f, f, K, 3);

        deconvolve::sb_deconvolve(u, f, K, 2.f / lambda, 1.f, 2.f * std::sqrt(2.f), 256.f, iters);

        u = utils::remove_padding(u, K);
        if (max != 1.0f)
            u.map(u * max);
        for (unsigned i = 0; i < rx * ry; i++) {
            fdata[c * rx * ry + i] = u.data[i];
        }
    }
    if (sequence_is_running == 0)
        updateprogress("Ready.", 0.);
    return 0;
}
