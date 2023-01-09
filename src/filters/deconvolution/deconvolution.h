#pragma once

#include <glib.h>

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif

typedef enum { DECONV_SB, DECONV_RL, DECONV_WIENER } nonblind_t;
typedef enum { PSF_BLIND, PSF_SELECTION, PSF_STARS, PSF_MANUAL, PSF_PREVIOUS } psftype_t;
typedef enum { BLIND_SI, BLIND_L0 } blind_t;
typedef enum { RL_MULT, RL_GD } rl_method_t;
typedef enum { PROFILE_GAUSSIAN, PROFILE_MOFFAT, PROFILE_DISK, PROFILE_AIRY } profile_t;
typedef enum { REG_TV_GRAD, REG_FH_GRAD, REG_NONE_GRAD, REG_TV_MULT, REG_FH_MULT, REG_NONE_MULT } regtype_t;


EXTERNC typedef struct estk_data {
	char* wisdom_file;
	float* fdata;
	unsigned rx;
	unsigned ry;
	unsigned nchans;
	unsigned ndata;
	int max_threads;
	psftype_t psftype;
	psftype_t oldpsftype;
	int ks; // Kernel size
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
	nonblind_t nonblindtype;
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
} estk_data;

#ifdef __cplusplus
}
#endif

EXTERNC float *estimate_kernel(estk_data *args, int max_threads);
#ifdef __cplusplus
}
#endif

EXTERNC float *gf_estimate_kernel(estk_data *args, int max_threads);
#ifdef __cplusplus
}
#endif

EXTERNC int split_bregman(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda, int iters, int max_threads);
#ifdef __cplusplus
}
#endif

EXTERNC int richardson_lucy(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda, int maxiter, float stopcriterion, int max_threads, int regtype, float stepsize, int stopcriterion_active);
#ifdef __cplusplus
}
#endif

EXTERNC int wienerdec(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float sigma, int max_threads);
#ifdef __cplusplus
}
#endif

EXTERNC int spectral_pre_adaption(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kdata, int kernelsize, float lambda, int max_threads, int deconv_algo);
#ifdef __cplusplus
}
#endif
