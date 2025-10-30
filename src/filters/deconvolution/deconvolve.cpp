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
#include "core/OS_utils.h"

extern "C" int wienerdec(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float sigma, int max_threads) {
    const int num_copies_required = 8;
    img_t<float>::use_threading(max_threads);
    img_t<float> u;
    for (unsigned c = 0 ; c < nchans ; c++) {
        img_t<float> f(rx, ry, 1, fdata + c * rx * ry);
        unsigned kc = (c < kchans ? c : 0);
        img_t<float> K(kernelsize, kernelsize, 1, kernel + kc * kernelsize * kernelsize);
        K.map(K / K.sum());
        float max = f.max();
        if (max != 0.0f && max != 1.0f)
            f.map(f / max);
        f = utils::add_padding(f, K);
        f.process_in_slices(get_available_memory(), num_copies_required, u, kernelsize / 2, [&K, sigma](img_t<float>& slice) {
            edgetaper(slice, slice, K, 3);
            deconvolve::wiener_deconvolve(slice, slice, K, sigma);
        });
        u = utils::remove_padding(u, K);
        if (max != 1.0f)
            u.map(u * max);
        memcpy(fdata + c * rx * ry, u.data.data(), u.data.size() * sizeof(float));
    }
    if (sequence_is_running == 0)
        set_progress_bar_data("Ready.", 0.);
    return 0;
}

extern "C" int fft_richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int maxiter, float stopcriterion, int max_threads, regtype_t regtype, float stepsize, int stopcriterion_active) {
    const int num_copies_required = (regtype == 0 || regtype == 3) ? 12 : 10;
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
        f.process_in_slices(get_available_memory(), num_copies_required, u, kernelsize / 2, [&K, lambda, maxiter, stopcriterion, regtype, stepsize, stopcriterion_active](img_t<float>& slice) {
            edgetaper(slice, slice, K, 3);
            deconvolve::rl_deconvolve_fft(slice, slice, K, 2.f / lambda, maxiter, stopcriterion, regtype, stepsize, stopcriterion_active);
        });

        u = utils::remove_padding(u, K);
        if (max != 1.0f)
            u.map(u * max);
        memcpy(fdata + c * rx * ry, u.data.data(), u.data.size() * sizeof(float));
    }
    if (sequence_is_running == 0)
        set_progress_bar_data("Ready.", 0.);
    return 0;
}

extern "C" int naive_richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int maxiter, float stopcriterion, int max_threads, regtype_t regtype, float stepsize, int stopcriterion_active) {
    const int num_copies_required = 7;
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
        f.process_in_slices(get_available_memory(), num_copies_required, u, kernelsize / 2, [&K, lambda, maxiter, stopcriterion, regtype, stepsize, stopcriterion_active](img_t<float>& slice) {
            edgetaper(slice, slice, K, 3);
            deconvolve::rl_deconvolve_naive(slice, slice, K, 2.f / lambda, maxiter, stopcriterion, regtype, stepsize, stopcriterion_active);
        });

        u = utils::remove_padding(u, 2*K.w, 2*K.w);
        if (max != 1.0f)
            u.map(u * max);
        memcpy(fdata + c * rx * ry, u.data.data(), u.data.size() * sizeof(float));
    }
    if (sequence_is_running == 0)
        set_progress_bar_data("Ready.", 0.);
    return 0;
}

extern "C" int split_bregman(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int iters, int max_threads) {
    const int num_copies_required = 16;
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
        f.process_in_slices(get_available_memory(), num_copies_required, u, kernelsize / 2, [&K, lambda, iters](img_t<float>& slice) {
            edgetaper(slice, slice, K, 3);
            deconvolve::sb_deconvolve(slice, slice, K, 2.f / lambda, 1.f, 2.f * std::sqrt(2.f), 256.f, iters);
        });

        u = utils::remove_padding(u, K);
        if (max != 1.0f)
            u.map(u * max);
        memcpy(fdata + c * rx * ry, u.data.data(), u.data.size() * sizeof(float));
    }
    if (sequence_is_running == 0)
        set_progress_bar_data("Ready.", 0.);
    return 0;
}
