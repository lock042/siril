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

#include "cdrizzlemap.h"
#include "cdrizzleutil.h"

#include "core/siril_log.h"

#include <assert.h>
#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*****************************************************************
 ERROR HANDLING
*/
void
driz_error_init(struct driz_error_t* error) {
  assert(error);

  error->last_message[0] = 0;
}

int
driz_error_check(struct driz_error_t* error, const char* message, int test) {
  if (! test) {
    driz_error_set_message(error, message);
    return 1;
  }

  return 0;
}

void
driz_error_set_message(struct driz_error_t* error, const char* message) {
  assert(error);
  assert(message);

  strncpy(error->last_message, message, MAX_DRIZ_ERROR_LEN - 1);
}

void
driz_error_format_message(struct driz_error_t* error, const char* format, ...) {
  /* See http://c-faq.com/varargs/vprintf.html
     for an explanation of how all this variable length argument list stuff
     works. */
  va_list argp;

  assert(error);
  assert(format);

  va_start(argp, format);
  (void)vsnprintf(error->last_message, MAX_DRIZ_ERROR_LEN - 1, format, argp);
  va_end(argp);
}

const char*
driz_error_get_message(struct driz_error_t* error) {
  assert(error);

  return error->last_message;
}

int
driz_error_is_set(struct driz_error_t* error) {
  assert(error);

  return error->last_message[0] != 0;
}

void
driz_error_unset(struct driz_error_t* error) {
  assert(error);

  driz_error_init(error);
}

/*****************************************************************
 DATA TYPES
*/
void
driz_param_dump(struct driz_args_t* p) {
  assert(p);

  siril_log_message(_("Drizzling parameters:\n"));
  siril_log_message(_("kernel:         %s\n"), kernel_enum2str(p->kernel));
  siril_log_message(_("pixel_fraction: %f\n"), p->pixel_fraction);
  siril_log_message(_("weight_scale:   %f\n"), p->weight_scale);
  siril_log_message(_("scale:          %f\n"), p->scale);
}

void
driz_param_init(struct driz_param_t* p) {
  assert(p);

  /* Kernel shape and size */
  p->kernel = kernel_square;
  p->pixel_fraction = 1.f;

  /* Exposure time */
  p->exposure_time = 1.f;

  /* Weight scale */
  p->weight_scale = 1.f;

  /* CPS / Counts */
  p->in_units = unit_counts;
  p->out_units = unit_counts;

  p->scale = 1.f;

  /* Input data */
  p->data = NULL;
  p->weights = NULL;
  p->pixmap = NULL;

  /* Output data */
  p->output_data = NULL;
  p->output_counts = NULL;

  p->nmiss = 0;
  p->nskip = 0;
  p->error = NULL;

  memset(p->cfa, 0, sizeof(p->cfa));
  p->cfadim = 1;
}

/*****************************************************************
 STRING TO ENUMERATION CONVERSIONS
*/
static const char* kernel_string_table[] = {
  "square",
  "gaussian",
  "point",
  "turbo",
  "lanczos2",
  "lanczos3",
  NULL
};

static const char* unit_string_table[] = {
  "counts",
  "cps",
  NULL
};

static const char* interp_string_table[] = {
  "nearest",
  "linear",
  "poly3",
  "poly5",
  "spline3",
  "sinc",
  "lsinc",
  "lan3",
  "lan5",
  NULL
};

static const char* bool_string_table[] = {
  "FALSE",
  "TRUE",
  NULL
};

static int
str2enum(const char* s, const char* table[], int* result, struct driz_error_t* error) {
  const char** it = table;

  assert(s);
  assert(table);
  assert(result);
  assert(error);

  while (*it != NULL) {
    if (strncmp(s, *it, 32) == 0) {
      *result = it - table;
      return 0;
    }
    ++it;
  }

  return 1;
}

int
kernel_str2enum(const char* s, enum e_kernel_t* result, struct driz_error_t* error) {
  if (str2enum(s, kernel_string_table, (int *)result, error)) {
    driz_error_format_message(error, "Unknown kernel type '%s'", s);
    return 1;
  }

  return 0;
}

int
unit_str2enum(const char* s, enum e_unit_t* result, struct driz_error_t* error) {
  if (str2enum(s, unit_string_table, (int *)result, error)) {
    driz_error_format_message(error, "Unknown unit type '%s'", s);
    return 1;
  }

  return 0;
}

int
interp_str2enum(const char* s, enum e_interp_t* result, struct driz_error_t* error) {
  if (str2enum(s, interp_string_table, (int *)result, error)) {
    driz_error_format_message(error, "Unknown interp type '%s'", s);
    return 1;
  }

  return 0;
}

const char*
kernel_enum2str(enum e_kernel_t value) {
  return kernel_string_table[value];
}

const char*
unit_enum2str(enum e_unit_t value) {
  return unit_string_table[value];
}

const char*
interp_enum2str(enum e_interp_t value) {
  return interp_string_table[value];
}

const char*
bool2str(bool_t value) {
  return bool_string_table[value ? 1 : 0];
}

/*****************************************************************
 NUMERICAL UTILITIES
*/
void
create_lanczos_lut(const int kernel_order, const size_t npix,
                   const float del, float* lanczos_lut) {
  size_t i;
  const float forder = (float)kernel_order;
  float poff;

  assert(lanczos_lut);
  assert(kernel_order < 6);

  /* Set the first value to avoid arithmetic problems */
  lanczos_lut[0] = 1.0;

  for (i = 1; i < npix; ++i) {
    poff = M_PI * (float)i * del;
    if (poff < M_PI * forder) {
      lanczos_lut[i] = sin(poff) / poff * sin(poff / forder) / (poff / forder);
    } else {
      lanczos_lut[i] = 0.0;
    }
  }
}
