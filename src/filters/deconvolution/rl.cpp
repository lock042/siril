#include <iostream>
#include <random>

#include "rl.hpp"
#include "utils.hpp"
#include "edgetaper.hpp"

extern "C" int richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda, int maxiter, float stopcriterion, int max_threads, int regtype, float stepsize, int stopcriterion_active) {

    // These replace smart parameters from the original implementation
    // No need for smapa in Siril as we don't need to be able to overwrite options using
    // environment variables
    int NORMALIZE_INPUT = 1;
    int EDGETAPER = 1;

    printf("rx: %d, ry:%d\n", rx, ry);

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

    richardsonlucy::rl_deconvolve(u, f, K, 2.f / lambda, maxiter, stopcriterion, regtype, stepsize, stopcriterion_active);

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

extern "C" int naive_richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda, int maxiter, float stopcriterion, int max_threads, int regtype, float stepsize, int stopcriterion_active) {

    // These replace smart parameters from the original implementation
    // No need for smapa in Siril as we don't need to be able to overwrite options using
    // environment variables
    int NORMALIZE_INPUT = 1;
    int EDGETAPER = 1;

    img_t<float>::use_threading(max_threads);

    img_t<float> f(rx, ry, nchans, fdata);
    img_t<float> K(kernelsize, kernelsize, 1, kernel);
    K.map(K / K.sum());
    img_t<float> u;

    float max = f.max();
    if (NORMALIZE_INPUT)
        f.map(f / f.max());

    if (EDGETAPER) {
        f = utils::add_padding(f, 2*K.w, 2*K.w);
        edgetaper(f, f, K, 3);
    }

    richardsonlucy::rl_deconvolve_naive(u, f, K, 2.f / lambda, maxiter, stopcriterion, regtype, stepsize, stopcriterion_active);

    if (EDGETAPER) {
        u = utils::remove_padding(u, 2*K.w, 2*K.w);
    }

    if (NORMALIZE_INPUT)
        u.map(u * max);

    // copy u.data.data back to image.fdata
    for (unsigned i = 0; i < rx * ry * nchans; i++)
        fdata[i] = u.data[i];

    return 0;
}

extern "C" int richardson_lucy_damped(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, int maxiter) {

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
        f = utils::add_padding(f, 2*K.w, 2*K.w);
        edgetaper(f, f, K, 3);
    }
    u = f;
    richardsonlucy::rl_deconvolve_damped(u, f, K, maxiter);

    if (EDGETAPER) {
        u = utils::remove_padding(u, 2*K.w, 2*K.w);
    }

    if (NORMALIZE_INPUT)
        u.map(u * max);

    // copy u.data.data back to image.fdata
    for (unsigned i = 0; i < rx * ry * nchans; i++)
        fdata[i] = u.data[i];

    return 0;
}
