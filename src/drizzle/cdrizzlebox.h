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
 * AURA code by team siril_free-astro. */

#ifndef CDRIZZLEBOX_H
#define CDRIZZLEBOX_H

#include "cdrizzleutil.h"

/**
dobox

This module does the actual mapping of input flux to output images
using "boxer", a code written by Bill Sparks for FOC geometric
distortion correction, rather than the "drizzling" approximation.

This works by calculating the positions of the four corners of a
quadrilateral on the output grid corresponding to the corners of the
input pixel and then working out exactly how much of each pixel in the
output is covered, or not.

In V1.6 this was simplified to use the DRIVAL routine and also to
include some limited multi-kernel support.
*/

integer_t
compute_bit_value(integer_t uuid);

#ifdef __cplusplus
extern "C" {
#endif
int
dobox(struct driz_param_t* p);
#ifdef __cplusplus
}
#endif

float
compute_area(float is, float js, const float x[4], const float y[4]);

float
boxer(float is, float js, const float x[4], const float y[4]);

typedef int (*kernel_handler_t)(struct driz_param_t*);

#endif /* CDRIZZLEBOX_H */
