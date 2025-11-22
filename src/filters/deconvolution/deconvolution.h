/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
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

#pragma once

#include <glib.h>

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

typedef enum { DECONV_SB, DECONV_RL, DECONV_WIENER } nonblind_t;
typedef enum { PSF_BLIND, PSF_SELECTION, PSF_STARS, PSF_MANUAL, PSF_PREVIOUS } psftype_t;
typedef enum { BLIND_SI, BLIND_L0 } blind_t;
typedef enum { RL_MULT, RL_GD } rl_method_t;
typedef enum { PROFILE_GAUSSIAN, PROFILE_MOFFAT, PROFILE_DISK, PROFILE_AIRY } profile_t;
typedef enum { REG_TV_GRAD, REG_FH_GRAD, REG_NONE_GRAD, REG_TV_MULT, REG_FH_MULT, REG_NONE_MULT } regtype_t;
typedef enum { BOTTOM_UP, TOP_DOWN, UNDEFINED } orientation_t;

EXTERNC typedef struct estk_data {
    destructor destroy_fn;
    orientation_t kernelorientation;
    float* fdata; // Reference only, the struct does not own fdata and must not free it
    fits *fit;    // Reference to the FITS image being processed
    unsigned rx;
    unsigned ry;
    unsigned nchans;
    unsigned ndata;
    int max_threads;
    psftype_t psftype;
    psftype_t oldpsftype;
    int ks; // Kernel size
    int kchans; // Kernel channels
    blind_t blindtype;
    // Anger-Delbracio-Facciolo l0 Descent
    float lambda;// = 4e-3f; // Lambda
    float lambda_ratio;// = 1/1.1f; // Scaling ratio of lambda for multiscale
    float lambda_min;// = 1e-3f; // Min lambda
    float gamma;// = 20.f;
    float iterations;// = 2;
    gboolean multiscale;// = FALSE;
    float scalefactor;// = 0.5f;
    float kernel_threshold_max;// = 0.f;
    gboolean remove_isolated;// = FALSE;
    gboolean better_kernel;// = FALSE;
    float upscaleblur;// = 0.f;
    float downscaleblur;// = 1.6f;
    float k_l1;// = 0.5f;
    float psf_fwhm;
    float psf_beta;
    float psf_angle;
    float psf_ratio;
    gboolean symkern;
    profile_t profile;
    // Goldstein-Fattal Spectral Irregularity
    int ninner;
    int ntries;
    int nouter;
    float compensationfactor;
    int medianfilter;
    float finaldeconvolutionweight;
    float intermediatedeconvolutionweight;
    // Non-blind deconvoluton stage
    nonblind_t nonblindtype; // Type of final non-blind deconvolution to perform
    float alpha; // = 1/3000;
    int finaliters; // Iters for iterative final methods
    float stopcriterion; // Stopping distance for Richardson-Lucy
    int stopcriterion_active;
    rl_method_t rl_method; // 0 for multiplicative, 1 for gradient descent
    float stepsize; // Step size for gradient descent versions of R-L
    regtype_t regtype; // R-L regularization type
    float airy_diameter;
    float airy_fl;
    float airy_wl;
    float airy_pixelsize;
    float airy_pixelscale;
    float airy_obstruction;
    gboolean save_after; // for the makepsf -savepsf= option
    char* savepsf_filename; // for the makepsf -savepsf= option
    gboolean recalc_ks; // for the makepsf stars option
    gboolean stars_need_clearing; // for the makepsf stars option
    gboolean previewing; // Added field
} estk_data;

EXTERNC void free_estk_data(void *p);
EXTERNC estk_data *alloc_estk_data();
EXTERNC int deconvolve_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
EXTERNC int estimate_only_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

EXTERNC float *estimate_kernel(estk_data *args, int max_threads);
EXTERNC float *gf_estimate_kernel(estk_data *args, int max_threads);
EXTERNC int split_bregman(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int iters, int max_threads);
EXTERNC int fft_richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int maxiter, float stopcriterion, int max_threads, regtype_t regtype, float stepsize, int stopcriterion_active);
EXTERNC int naive_richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float lambda, int maxiter, float stopcriterion, int max_threads, regtype_t regtype, float stepsize, int stopcriterion_active);
EXTERNC int wienerdec(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, unsigned kchans, float sigma, int max_threads);
