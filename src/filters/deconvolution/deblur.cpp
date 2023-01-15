#include <iostream>
#include <random>

#include "deblur.hpp"
#include "utils.hpp"
#include "edgetaper.hpp"
#include "chelperfuncs.h"
#include "rlstrings.h"


extern "C" int split_bregman(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda, int iters, int max_threads) {

    // These replace smart parameters from the original implementation
    // No need for smapa in Siril as we don't need to be able to overwrite options using
    // environment variables
    int NORMALIZE_INPUT = 1;
    int EDGETAPER = 1;
    float CONTINUATION_BETA_INIT = 1.f;
    float CONTINUATION_BETA_RATE = 2.f * std::sqrt(2.f);
    float CONTINUATION_BETA_MAX = std::pow(2.f, 8.f);

    img_t<float>::use_threading(max_threads);

    img_t<float> f(rx, ry, nchans, fdata);
    img_t<float> K(kernelsize, kernelsize, 1, kernel);
    K.map(K / K.sum());
    img_t<float> u;

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
    for (unsigned i = 0; i < rx * ry * nchans; i++)
        fdata[i] = u.data[i];
    if (sequence_is_running == 0)
        updateprogress("Ready.", 0.);

    return 0;
}

