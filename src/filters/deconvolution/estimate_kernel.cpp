#include <iostream>
#include <random>

#include "estimate_kernel.hpp"
#include "deconvolution.h"

#include <iostream>

extern "C" float *estimate_kernel(estk_data *args) {
    img_t<float>::use_threading(1);
    options opts;
    opts.ks = args->ks;
    opts.input = ""; // Not used in siril
    opts.output = ""; // Not used in siril
    opts.lambda = args->lambda;
    opts.lambda_ratio = args->lambda_ratio;
    opts.lambda_min = args->lambda_min;
    opts.gamma = args->gamma;
    opts.iterations = args->iterations;
    opts.multiscale = (bool) args->multiscale;
    opts.scalefactor = args->scalefactor;
    opts.kernel_threshold_max = args->kernel_threshold_max;
    opts.remove_isolated = (bool) args->remove_isolated;
    opts.outputsharp = "";
    opts.verbose = true;
    opts.debug = "";
    opts.better_kernel = (bool) args->better_kernel;
    opts.warmg = false;
    opts.warmk = false;
    opts.upscaleblur = args->upscaleblur;
    opts.downscaleblur = args->downscaleblur;
    opts.initu = ""; // Not used in siril
    opts.k_l1 = args->k_l1;
    opts.use_filters = false;

    img_t<float> v(args->rx, args->ry, args->nchans, args->fdata);
    preprocess_image(v, v, opts);

    img_t<float> initu;
    initu = v; // initu allows for use of an initial sharpened image but this capability is not currently used in Siril

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

