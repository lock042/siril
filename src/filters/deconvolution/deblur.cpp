#include <iostream>
#include <random>

#include "deblur.hpp"
#include "utils.hpp"
#include "edgetaper.hpp"
#include "chelperfuncs.h"
#include "rlstrings.h"


extern "C" int split_bregman(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int iters, int max_threads) {

    // These replace smart parameters from the original implementation
    // No need for smapa in Siril as we don't need to be able to overwrite options using
    // environment variables
    int NORMALIZE_INPUT = 1;
    int EDGETAPER = 1;
    float CONTINUATION_BETA_INIT = 1.f;
    float CONTINUATION_BETA_RATE = 2.f * std::sqrt(2.f);
    float CONTINUATION_BETA_MAX = std::pow(2.f, 8.f);

    img_t<float>::use_threading(max_threads);
    img_t<float> u;

    for (unsigned c = 0 ; c < nchans ; c++) {
        unsigned kc = (c < kchans ? c : 0);
        img_t<float> f(rx, ry, 1, fdata + c * rx * ry);
        img_t<float> K(kernelsize, kernelsize, 1, kernel + kc * kernelsize * kernelsize);
        K.map(K / K.sum());

        float max = f.max();
        if (NORMALIZE_INPUT)
            f.map(f / f.max());

        if (EDGETAPER) {
            f = utils::add_padding(f, K);
            edgetaper(f, f, K, 3);
        }

        deblur::rof::split_continuation(u, f, K, 2.f / lambda, CONTINUATION_BETA_INIT,
                                        CONTINUATION_BETA_RATE, CONTINUATION_BETA_MAX, iters);

        if (EDGETAPER) {
            u = utils::remove_padding(u, K);
        }

        if (NORMALIZE_INPUT)
            u.map(u * max);
        // copy u.data.data back to image.fdata
        for (unsigned i = 0; i < rx * ry; i++)
            fdata[c * rx * ry + i] = u.data[i];
    }

    if (sequence_is_running == 0)
        updateprogress("Ready.", 0.);

    return 0;
}

