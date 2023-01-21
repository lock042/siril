/*
MIT License

Copyright (c) 2018 Jérémy Anger, Gabriele Facciolo, Mauricio Delbracio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <iostream>
#include <random>

#include "estimate_kernel.hpp"
#include "deconvolution.h"

#include <iostream>

extern "C" float *estimate_kernel(estk_data *args, int max_threads) {
    if (!cppfftwmultithreaded)
        max_threads = 1;
    img_t<float>::use_threading(max_threads);
    options opts;
    opts.ks = args->ks;
    opts.lambda = args->lambda;
    opts.lambda_ratio = args->lambda_ratio;
    opts.lambda_min = args->lambda_min;
    opts.gamma = args->gamma;
    opts.iterations = args->iterations;
    opts.multiscale = (bool) args->multiscale;
    opts.scalefactor = args->scalefactor;
    opts.kernel_threshold_max = args->kernel_threshold_max;
    opts.remove_isolated = (bool) args->remove_isolated;
    opts.verbose = true;
    opts.better_kernel = (bool) args->better_kernel;
    opts.warmg = false;
    opts.warmk = false;
    opts.upscaleblur = args->upscaleblur;
    opts.downscaleblur = args->downscaleblur;
    opts.k_l1 = args->k_l1;
    opts.use_filters = false;

    img_t<float> v(args->rx, args->ry, args->nchans, args->fdata);
    preprocess_image(v, v, opts);

    img_t<float> initu;
    initu = v;

    img_t<float> k;
    img_t<float> u;
    if (opts.multiscale) {
        multiscale_l0_kernel_estimation(k, u, v, opts);
    } else {
        l0_kernel_estimation(k, u, v, initu, opts);
    }
    opts.ks = k.w; // Update ks which was formerly an upper bound for the kernel size, now it contains the actual size of k

    float *kernel = (float*) calloc(k.w * k.h, sizeof(float)); // Kernel is passed as a NULL pointer, it is the responsibility of the calling function to free kernel
    for (int i = 0; i < k.w * k.h; i++)
        kernel[i] = k.data[i];

    return kernel;
}

