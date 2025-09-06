/*
 * The following copyright notice applies to drizzle code:

Copyright (C) 2011,2014 Association of Universities for Research in
Astronomy (AURA)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

    1. Redistributions of source code must retain the above
      copyright notice, this list of conditions and the following
      disclaimer.

    2. Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials
      provided with the distribution.

    3. The name of AURA and its representatives may not be used to
      endorse or promote products derived from this software without
      specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.

 * The code as integrated into Siril is modified from the original
 * AURA code by team free-astro. */

#ifndef CDRIZZLEUTIL_H
#define CDRIZZLEUTIL_H
#include "driz_portability.h"
#include "core/siril.h"
#include <assert.h>
#include <errno.h>
#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#if __STDC_VERSION__ >= 199901L
#include <stdint.h>
#endif
#include <stdlib.h>

#include "core/sequence_filtering.h" // for seq_image_filter
#include "registration/registration.h" // for framing_type

/*****************************************************************
 ERROR HANDLING
*/
#define MAX_DRIZ_ERROR_LEN 512

struct driz_error_t {
  char last_message[MAX_DRIZ_ERROR_LEN];
};

void driz_error_init(struct driz_error_t* error);
int driz_error_check(struct driz_error_t* error, const char* message, int test);
void driz_error_set_message(struct driz_error_t* error, const char* message);
void driz_error_format_message(struct driz_error_t* error, const char* format, ...);
const char* driz_error_get_message(struct driz_error_t* error);
int driz_error_is_set(struct driz_error_t* error);
void driz_error_unset(struct driz_error_t* error);

/*****************************************************************
 CONVENIENCE MACROS
*/
#if !defined(MIN)
  #define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif

#if !defined(MAX)
  #define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#define CLAMP_ABOVE(x, low)  (((x) < low) ? (low) : (x))
#define CLAMP_BELOW(x, high)  (((x) > high) ? (high) : (x))

#ifdef __GNUC__
#define UNUSED_PARAM __attribute__((unused))
#else
#define UNUSED_PARAM
#endif

#define MAX_SINGLE 3.4028234663852886e+38
#define MIN_SINGLE 1.1754943508222875e-38

#define MAX_COEFFS 128
#define COEFF_OFFSET 100

#undef TRUE
#define TRUE 1

#undef FALSE
#define FALSE 0

/*****************************************************************
 DATA TYPES
*/
typedef int integer_t;
#if __STDC_VERSION__ >= 199901L
typedef int_fast8_t bool_t;
#else
typedef unsigned char bool_t;
#endif

enum e_kernel_t {
  kernel_square,
  kernel_gaussian,
  kernel_point,
  kernel_turbo,
  kernel_lanczos2,
  kernel_lanczos3,
  kernel_LAST
};

enum e_unit_t {
  unit_counts,
  unit_cps
};

enum e_interp_t {
  interp_nearest,
  interp_bilinear,
  interp_poly3,
  interp_poly5,
  interp_spline3,
  interp_sinc,
  interp_lsinc,
  interp_lanczos3,
  interp_lanczos5,
  interp_LAST
};

/* Lanczos values */
struct lanczos_param_t {
  size_t nlut;
  float* lut;
  float sdp;
  integer_t nbox;
  float space;
  float misval;
};

typedef struct _imgmap_t {
	float* xmap;
	float* ymap;
	int rx;
	int ry;
} imgmap_t;

struct driz_args_t {
  bool_t is_bayer; /* Is this a Bayer drizzle? */
  /* Parameters to be copied into the driz_param_t for each frame */
  enum e_kernel_t kernel; /* Kernel shape and size */
  float scale;
  float weight_scale; /* Weight scale */
  fits *flat; /* Flat, for multiplying the pixel_count */
  gboolean use_flats; /* Whether to use master flat as weights */
  float pixel_fraction;
  BYTE cfa[36];
  size_t cfadim;
  float max_weight[3]; /* stores the max output weight per channel to renormalize the final output */
};

struct driz_param_t {
  struct driz_args_t *driz; /* sequence-wide drizzle args */
  regdata *current_regdata; /* Current reg data (to get Homography matrices, if required)*/
  int threads; /* threads allowed */

  /* Options */
  enum e_kernel_t kernel; /* Kernel shape and size */
  float          pixel_fraction; /* was: PIXFRAC */
  float           exposure_time; /* Exposure time was: EXPIN */
  // Siril doesn't use exposure_time because frame exposure weighting can already be
  // done in the stacking code, so this is always set to 1.f in the initializer.
  float           weight_scale; /* Weight scale was: WTSCL */
  enum e_unit_t   in_units; /* CPS / counts was: INCPS, either counts or CPS */
  enum e_unit_t   out_units; /* CPS / counts was: INCPS, either counts or CPS */
  // Siril typically works with images directly from a camera, with pixel values
  // representing counts. These are automatically initialised to counts in the
  // initializer but are left in the struct in case they are ever of use.

  /* Scaling */
  float scale;

  /* CFA filter pattern */
  BYTE cfa[36];
  size_t cfadim; /*1, 2 or 6*/

  /* Image subset */
  integer_t xmin;
  integer_t xmax;
  integer_t ymin;
  integer_t ymax;

  /* Input images */
  fits *data; // Per-image
  fits *weights; // Per-image
  imgmap_t *pixmap; // Per-image

  /* Output images */
  fits *output_data;
  fits *output_counts;  /* was: COU */

  /* Other output */
  integer_t nmiss;
  integer_t nskip;
  struct driz_error_t* error;
};

/**
Initialize all of the members of the drizzle_param_t to sane default
values, mostly zeroes.  Note, these are not *meaningful* values, just
ones that help with memory management etc.  It is up to users of the
struct, e.g. cdrizzle_, to fill the struct with valid parameters.
*/
#ifdef __cplusplus
extern "C" {
#endif
void
driz_param_init(struct driz_param_t* p);
#ifdef __cplusplus
}
#endif
void
driz_param_dump(struct driz_args_t* p);

/****************************************************************************/
/* ARRAY ACCESSORS */

/* Modified to use Siril struct fits-based accessors in place of numpy ones */

static inline_macro void
get_dimensions(fits *image, integer_t size[2]) {

  /* Put dimensions in xy order */
  size[0] = image->rx;
  size[1] = image->ry;

  return;
}

static inline_macro float
get_xmap(imgmap_t *p, integer_t xpix, integer_t ypix) {
  return *(float*) (p->xmap + ((xpix + ypix * p->rx)));
}

static inline_macro float
get_ymap(imgmap_t *p, integer_t xpix, integer_t ypix) {
  return *(float*) (p->ymap + ((xpix + ypix * p->rx)));
}

#if defined(LOGGING) && defined(CHECK_OOB)

static inline_macro int
oob_pixel(fits *image, integer_t xpix, integer_t ypix) {
    if ((xpix < 0 || xpix >= image->rx) || (ypix < 0 || ypix >= image->ry)) {
        siril_debug_print("Point [%d,%d] is outside of [%d, %d]",
                xpix, ypix, image->rx, image->ry);
        return 1;
    }
    return 0;
}

#else

#define oob_pixel(image, xpix, ypix) 0

#endif

static inline_macro float
get_pixel(fits *image, integer_t xpix, integer_t ypix, integer_t chan) {
  return *(float*) (image->fdata + xpix + ypix * image->rx + chan * image->rx * image->ry);
}

static inline_macro float
get_pixel_at_pos(fits *image, integer_t pos) {
  float *imptr;
  imptr = (float *) image->fdata;
  return imptr[pos];
}

static inline_macro void
set_pixel(fits *image, integer_t xpix, integer_t ypix, integer_t chan, float value) {
  *(float*) (image->fdata + xpix + ypix * image->rx + chan * image->rx * image->ry) = value;
  return;
}


/*****************************************************************
 STRING TO ENUMERATION CONVERSIONS
*/
int
kernel_str2enum(const char* s, enum e_kernel_t* result, struct driz_error_t* error);

int
unit_str2enum(const char* s, enum e_unit_t* result, struct driz_error_t* error);

int
interp_str2enum(const char* s, enum e_interp_t* result, struct driz_error_t* error);

const char*
kernel_enum2str(enum e_kernel_t value);

const char*
unit_enum2str(enum e_unit_t value);

const char*
interp_enum2str(enum e_interp_t value);

const char*
bool2str(bool_t value);

/*****************************************************************
 NUMERICAL UTILITIES
*/
/**
Fill up a look-up-table of Lanczos interpolation kernel values for
rapid weighting determination for kernel == kernel_lanczos.

@param kernel_order the order of the kernel.
@param npix the size of the lookup table
@param del the spacings of the sampling of the function
@param lanczos_lut 1d array of lookup values.  This is a single-sided Lanczos
   function with lanczos_lut[0] being the central value.

Note that no checking is done to see whether the values are sensible.

was: FILALU
*/
void
create_lanczos_lut(const int kernel_order, const size_t npix,
                   const float del, float* lanczos_lut);

/**
 Round to nearest integer in a way that mimics fortrans NINT
*/
static inline_macro integer_t
fortran_round(const float x) {
  return (x >= 0) ? (integer_t)floor(x + .5) : (integer_t)-floor(.5 - x);
}

static inline_macro float
min_floats(const float* a, const integer_t size) {
  const float* end = a + size;
  float value = MAX_SINGLE;
  for ( ; a != end; ++a)
    if (*a < value)
      value = *a;
  return value;
}

static inline_macro float
max_floats(const float* a, const integer_t size) {
  const float* end = a + size;
  float value = -MAX_SINGLE;
  for ( ; a != end; ++a)
    if (*a > value)
      value = *a;
  return value;
}

// High level function to apply drizzle
int apply_drizzle(struct driz_args_t *driz);

#endif /* CDRIZZLEUTIL_H */
