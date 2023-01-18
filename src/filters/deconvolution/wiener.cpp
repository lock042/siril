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

#include "wiener.hpp"
#include "utils.hpp"
#include "edgetaper.hpp"

extern "C" int wienerdec(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float sigma, int max_threads) {

    // These replace smart parameters from the original implementation
    // No need for smapa in Siril as we don't need to be able to overwrite options using
    // environment variables
    int NORMALIZE_INPUT = 1;
    int EDGETAPER = 1;

    img_t<float>::use_threading(max_threads);
    img_t<float> u;
    for (unsigned c = 0 ; c < nchans ; c++) {
        img_t<float> f(rx, ry, 1, fdata + c * rx * ry);
        unsigned kc = (c < kchans ? c : 0);
        img_t<float> K(kernelsize, kernelsize, 1, kernel + kc * kernelsize * kernelsize);
        K.map(K / K.sum());

        float max = f.max();
        if (NORMALIZE_INPUT)
            f.map(f / f.max());

        if (EDGETAPER) {
            f = utils::add_padding(f, K);
            edgetaper(f, f, K, 3);
        }

        wiener::wiener_deconvolve(u, f, K, sigma);

        if (EDGETAPER) {
            u = utils::remove_padding(u, K);
        }

        if (NORMALIZE_INPUT)
            u.map(u * max);

        // copy u.data.data back to image.fdata
        for (unsigned i = 0; i < rx * ry; i++)
            fdata[c * rx * ry + i] = u.data[i];
    }

    return 0;
}
