#pragma once

#include <glib.h>

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif

EXTERNC typedef struct estk_data {
	float* fdata;
	unsigned rx;
	unsigned ry;
	unsigned nchans;
	int ks; // Kernel size
	// Anger-Delbracio-Facciolo
	float lambda;// = 4e-3f; // Lambda
	float lambda_ratio;// = 1/1.1f; // Scaling ratio of lambda for multiscale
	float lambda_min;// = 1e-2f; // Min lambda
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
	float alpha; // = 1/3000;
	float psf_fwhm;
	float psf_beta;
	float psf_angle;
	float psf_ratio;
	// Goldstein-Fattal
	int ninner;
	int ntries;
	int nouter;
	float compensationfactor;
	float medianfilter;
	float finaldeconvolutionweight;
	float intermediatedeconvolutionweight;
	int blindtype;
} estk_data;
#ifdef __cplusplus
}
#endif

EXTERNC float *estimate_kernel(estk_data *args);
#ifdef __cplusplus
}
#endif

EXTERNC float *gf_estimate_kernel(estk_data *args);
#ifdef __cplusplus
}
#endif

EXTERNC int split_bregman(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda);
#ifdef __cplusplus
}
#endif

EXTERNC int richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda, int maxiter);
#ifdef __cplusplus
}
#endif

EXTERNC int stochastic(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda);
#ifdef __cplusplus
}
#endif
