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
  kernel_tophat,
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
	float* pixmap;
	int rx;
	int ry;
} imgmap_t;

struct driz_args_t {
  /* Siril sequence data */
  sequence *seq; /* Sequence to operate on */
  int reference_image; /* reference image */
  bool_t is_bayer; /* Is this a Bayer drizzle? */
  bool_t use_wcs; /* Use WCS mapping? If not, Homography mapping will be used */
  regdata *ref_regdata; /* Reference reg data */
  struct wcsprm *refwcs; /* Reference WCS */
  gchar *prefix;
  gchar *new_seq_name;
  imgdata *imgparam;
  regdata *regparam;
  int new_total;
  gboolean load_new_sequence;
  float scale;
  BYTE* success;
};

struct driz_param_t {
  struct driz_args_t *driz; /* sequence-wide drizzle args */
  regdata *current_regdata; /* Current reg data */

  /* Options */
  enum e_kernel_t kernel; /* Kernel shape and size */
  float          pixel_fraction; /* was: PIXFRAC */
  float           exposure_time; /* Exposure time was: EXPIN */
  float           weight_scale; /* Weight scale was: WTSCL */
  float           fill_value; /* Filling was: FILVAL */
  bool_t          do_fill; /* was: FILL */
  enum e_unit_t   in_units; /* CPS / counts was: INCPS, either counts or CPS */
  enum e_unit_t   out_units; /* CPS / counts was: INCPS, either counts or CPS */
  integer_t       uuid; /* was: UNIQID */

  /* Scaling */
  float scale;

  /* Image subset */
  integer_t xmin;
  integer_t xmax;
  integer_t ymin;
  integer_t ymax;

  /* Blotting-specific parameters */
  enum e_interp_t interpolation; /* was INTERP */
  float ef;
  float misval;
  float sinscl;
  float kscale;

  /* Input images */
  fits *data; // Per-image
  fits *weights; // Per-image
  imgmap_t *pixmap; // Per-image

  /* Output images */
  fits *output_data;
  fits *output_counts;  /* was: COU */
  fits *output_context; /* was: CONTIM */

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
void
driz_param_init(struct driz_param_t* p);

void
driz_param_dump(struct driz_param_t* p);


/****************************************************************************/
/* LOGGING */

#ifdef LOGGING
extern FILE *driz_log_handle;


static inline_macro FILE*
driz_log_init(FILE *handle) {
    const char* dirs[] = {"TMPDIR", "TMP", "TEMP", "TEMPDIR"};
    int i;
    char *p;
    char buf[1024];
#ifdef _WIN32
    const char dirsep = '\\';
#else
    const char dirsep = '/';
#endif

/*    for (i = 0; i < 4; ++i) {
        p = getenv(dirs[i]);
        if (p) {
            sprintf(buf, "%s%cdrizzle.log", p, dirsep);
            break;
        }
    }
    if (!p) {
*/
        sprintf(buf, "drizzle.log");
//    }

    handle = fopen(buf, "a");
    if (handle) {
        setbuf(handle, 0);
    }
    return handle;
}


static inline_macro int
driz_log_close(FILE *handle) {
    if (handle) {
        return fclose(handle);
    }
}

static inline_macro int
driz_log_message(const char* message) {
    if (!driz_log_handle) {
        driz_log_handle = driz_log_init(driz_log_handle);
        if (!driz_log_handle) {
            return 1;
        }
    }

    fputs(message, driz_log_handle);
    fputs("\n", driz_log_handle);
    return 0;
}

#else
static inline_macro void *
driz_log_idem(void *ptr) {
    return ptr;
}

#define driz_log_init(handle) driz_log_idem(handle)
#define driz_log_close(handle) driz_log_idem(handle)
#define driz_log_message(message) driz_log_idem(message)

#endif

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

static inline_macro float*
get_pixmap(imgmap_t *p, integer_t xpix, integer_t ypix) {
  return (float*) p->pixmap + ((xpix + ypix * p->rx) * 2);
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
get_pixel(fits *image, integer_t xpix, integer_t ypix) {
  return *(float*) image->fdata + xpix + ypix * image->rx;
}

static inline_macro float
get_pixel_at_pos(fits *image, integer_t pos) {
  float *imptr;
  imptr = (float *) image->fdata;
  return imptr[pos];
}

static inline_macro void
set_pixel(fits *image, integer_t xpix, integer_t ypix, float value) {
  *(float*) (image->fdata + xpix + ypix * image->rx) = value;
  return;
}

/* For the context image we will blatantly misuse a FITS structure (because it's
 * convenient) by storing integer_t (aka uint8_t) data in the fdata array. This
 * applies to the next 3 functions.
 *
 * NOTE: the context image is optional, so for the moment we are simply not using it.
 *
 * TODO: this should really be made less ugly!! */

static inline_macro int
get_bit(fits *image, integer_t xpix, integer_t ypix, integer_t bitval) {
  integer_t value;
  value = *(integer_t*) (image->fdata + xpix + ypix * image->rx) & bitval;
  return value? 1 : 0;
}

static inline_macro void
set_bit(fits *image, integer_t xpix, integer_t ypix, integer_t bitval) {
  *(integer_t*) (image->fdata + xpix + ypix * image->rx) |= bitval;
  return;
}

static inline_macro void
unset_bit(fits *image, integer_t xpix, integer_t ypix) {
  *(integer_t*) (image->fpdata + ypix * image->rx +xpix) = 0;
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

void
put_fill(struct driz_param_t* p, const float fill_value);

/**
 Calculate the refractive index of MgF2 for a given C wavelength (in
 nm) using the formula given by Trauger (1995)
// Not required in Siril
float
mgf2(float lambda);
*/

/**
Weighted sum of 2 real vectors.

was: WSUMR
*/
static inline_macro void
weighted_sum_vectors(const integer_t npix,
                     const float* a /*[npix]*/, const float w1,
                     const float* b /*[npix]*/, const float w2,
                     /* Output arguments */
                     float* c /*[npix]*/) {
  float* c_end = c + npix;

  assert(a);
  assert(b);
  assert(c);

  while(c != c_end)
    *(c++) = *(a++) * w1 + *(b++) * w2;
}

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

/**
Evaluate a 3rd order radial geometric distortion in 2d
X version. Note that there is no zero order coefficient
as this is physically meaningless.

@param x The x coordinate

@param y The y coordinate

@param co An array of length 4 of coefficients

@param[out] xo The distorted x coordinate

@param[out] yo The distorted y coordinate
*/
static inline_macro void
rad3(const float x, const float y, const float* co,
     /* Output parameters */
     float* xo, float* yo) {
  float r, f;

  assert(co);
  assert(xo);
  assert(yo);

  r = sqrt(x*x + y*y);

  f = 1.0 + co[0] + co[1]*r + co[2]*r*r;
  *xo = f*x;
  *yo = f*y;
}

// High level function to apply drizzle
int apply_drizzle(struct driz_args_t *driz);

// FITS functions for saving 2-channel (x,y) mapping pixmaps as FITS bintables
int save_floats_to_fits_bintable(const char* filename, float* data, int rx, int ry, int numChannels);
int read_floats_from_fits_bintable(const char* filename, float** data, int* rx, int* ry, int* nchans);

#endif /* CDRIZZLEUTIL_H */
