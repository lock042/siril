#include <iostream>
#include <random>

#include "rl.hpp"
#include "utils.hpp"
#include "edgetaper.hpp"

extern "C" int richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda, int maxiter) {

    // These replace smart parameters from the original implementation
    // No need for smapa in Siril as we don't need to be able to overwrite options using
    // environment variables
    int NORMALIZE_INPUT = 1;
    int EDGETAPER = 1;

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

    richardsonlucy::rl_deconvolve(u, f, K, 2.f / lambda, maxiter);

    if (EDGETAPER) {
        u = utils::remove_padding(u, K);
    }

    if (NORMALIZE_INPUT)
        u.map(u * max);

    // copy u.data.data back to image.fdata
    for (unsigned i = 0; i < rx * ry * nchans; i++)
        fdata[i] = u.data[i];

    return 0;
}

