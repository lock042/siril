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

#include "driz_portability.h"
#include "cdrizzlemap.h"
#include "cdrizzlebox.h"
#include "cdrizzleutil.h"
#include "algos/demosaicing.h"

#include <assert.h>
#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdlib.h>


//static char buf[1024];


/** ------------------------------------------------------------------------------
 * Update the flux and counts in the output image using a weighted average
 *
 * p:   structure containing options, input, and output
 * ii:  x coordinate in output images
 * jj:  y coordinate in output images
 * d:   new contribution to weighted flux
 * vc:  previous value of counts
 * dow: new contribution to weighted counts
 */

inline_macro static int
update_data(struct driz_param_t* p, const integer_t ii, const integer_t jj,
            const integer_t chan, const float d, const float vc, const float dow) {

  float vc_plus_dow;

  if (dow == 0.0f) return 0;

  vc_plus_dow = vc + dow;

  if (vc == 0.0f) {
    if (oob_pixel(p->output_data, ii, jj)) {
      driz_error_format_message(p->error, "OOB in output_data[%d,%d]", ii, jj);
      return 1;
    } else {
      set_pixel(p->output_data, ii, jj, chan, d);
    }

  } else {
    if (oob_pixel(p->output_data, ii, jj)) {
      driz_error_format_message(p->error, "OOB in output_data[%d,%d]", ii, jj);
      return 1;
    } else {
      float value;
      value = (get_pixel(p->output_data, ii, jj, chan) * vc + dow * d) / (vc_plus_dow);
      set_pixel(p->output_data, ii, jj, chan, value);
    }
  }

  if (oob_pixel(p->output_counts, ii, jj)) {
    driz_error_format_message(p->error, "OOB in output_counts[%d,%d]", ii, jj);
    return 1;
  } else {
    set_pixel(p->output_counts, ii, jj, chan, vc_plus_dow);
  }

  return 0;
}

/** --------------------------------------------------------------------------------------------------
 * Compute area of box overlap. Calculate the area common to input clockwise polygon x(n), y(n) with
 * square (is, js) to (is+1, js+1). This version is for a quadrilateral. Used by do_square_kernel.
 *
 * is: x coordinate of a pixel on the output image
 * js: y coordinate of a pixel on the output image
 * x:  x coordinates of endpoints of quadrilateral containing flux of input pixel
 * y:  y coordinates of endpoints of quadrilateral containing flux of input pixel
 */

float
compute_area(float is, float js, const float x[4], const float y[4]) {
  int ipoint, jpoint, idim, jdim, iside, outside, count;
  int positive[2];
  float area, width;
  float midpoint[2], delta[2];
  float border[2][2], segment[2][2];

  /* The area for a qadrilateral clipped to a square of unit length whose sides are
   * aligned with the axes. The area is computed by computing the area under each
   * line segment clipped to the boundary of three sides of the sqaure. Since the
   * computed width is positive for two of the sides and negative for the other two,
   * we subtract the area outside the quadrilateral without any extra code.
   */
  area = 0.0;

  border[0][0] = is - 0.5;
  border[0][1] = js - 0.5;
  border[1][0] = is + 0.5;
  border[1][1] = js + 0.5;

  for (ipoint = 0; ipoint < 4; ++ ipoint) {
    jpoint = (ipoint + 1) & 03; /* Next point in cyclical order */

    segment[0][0] = x[ipoint];
    segment[0][1] = y[ipoint];
    segment[1][0] = x[jpoint];
    segment[1][1] = y[jpoint];

    /* Compute the endpoints of the line segment that
     * lie inside the border (possibly the whole segment)
     */

    for (idim = 0, count = 3; idim < 2; ++ idim) {
      for (iside = 0; iside < 2; ++ iside, -- count) {

        delta[0] = segment[0][idim] - border[iside][idim];
        delta[1] = segment[1][idim] - border[iside][idim];

        positive[0] = delta[0] > 0.0;
        positive[1] = delta[1] > 0.0;

        /* If both deltas have the same signe there is no baundary crossing
         */
        if (positive[0] == positive[1]) {
          /* A diagram will convince that you decide a point is
           * inside or outside the boundary by the following test
           */
          if (positive[0] == iside) {
            /* Segment is entirely outside the boundary */
            if (count == 0) {
              /* Implicitly multiplied by 1.0, the square height */
              width = segment[1][0] - segment[0][0];
              area += width;
            } else {
              goto _nextsegment;
            }

          } else {
            /* Segment entirely within the boundary */
            if (count == 0) {
              /* Use the trapezoid formula to compute the area under the
               * segment. Delta is the distance to the top of the square
               * and is negative or zero for the segment inside the square
               */
              width = segment[1][0] - segment[0][0];
              area += 0.5 * width * ((1.0 + delta[0]) + (1.0 + delta[1]));
            }
          }

        } else {
          /* If ends of the line segment are on opposite sides of the
          * boundary, calculate midpoint, the point of intersection
          */
          outside = positive[iside];
          jdim = (idim + 1) & 01; /* the other dimension */

          midpoint[idim] = border[iside][idim];

          midpoint[jdim] =
            (delta[1] * segment[0][jdim] - delta[0] * segment[1][jdim]) /
            (delta[1] - delta[0]);

          if (count == 0) {
            /* If a segment cross the boundary the formula for its area
             * is a combination of the formulas for segments entirely
             * inside and outside.
             */
            if (outside == 0) {
              width = midpoint[0] - segment[0][0];
              area += width;
              width = segment[1][0] - midpoint[0];
              /* Delta[0] is at the crossing point and thus zero */
              area += 0.5 * width * (1.0 + (1.0 + delta[1]));
            } else {
              width = segment[1][0] - midpoint[0];
              area += width;
              width =  midpoint[0] - segment[0][0];
              /* Delta[1] is at the crossing point and thus zero */
              area += 0.5 * width * ((1.0 + delta[0]) + 1.0);
            }

          } else {
            /* Clip segment against each boundary except the last */
            segment[outside][0] = midpoint[0];
            segment[outside][1] = midpoint[1];
          }
        }
      }
    }

    _nextsegment: continue;
  }

  return fabs(area);
}

/** --------------------------------------------------------------------------------------------------
 * Calculate overlap between an arbitrary rectangle, aligned with the axes, and a pixel.
 * This is a simplified version of the compute_area, only valid if axes are nearly aligned.
 * Used by do_kernel_turbo.
 *
 * i:    the x coordinate of a pixel on the output image
 * j:    the y coordinate of a pixel on the output image
 * xmin: the x coordinate of the lower edge of rectangle containing flux of input pixel
 * xmax: the x coordinate of the upper edge of rectangle containing flux of input pixel
 * ymin: the y coordinate of the lower edge of rectangle containing flux of input pixel
 * ymax: the y coordinate of the upper edge of rectangle containing flux of input pixel
 */

static inline_macro float
over(const integer_t i, const integer_t j,
     const float xmin, const float xmax,
     const float ymin, const float ymax) {
    float dx, dy;

    assert(xmin <= xmax);
    assert(ymin <= ymax);

    dx = MIN(xmax, (float)(i) + 0.5) - MAX(xmin, (float)(i) - 0.5);
    dy = MIN(ymax, (float)(j) + 0.5) - MAX(ymin, (float)(j) - 0.5);

    if (dx > 0.0 && dy > 0.0)
        return dx*dy;

    return 0.0;
}

/** --------------------------------------------------------------------------------------------------
 * The kernel assumes all the flux in an input pixel is at the center
 *
 * p: structure containing options, input, and output
 */

static int
do_kernel_point(struct driz_param_t* p) {
    struct scanner s;
    integer_t i, j, ii, jj;
    integer_t /*ybounds[2],*/ osize[2];
    float scale2, vc[3], d, dow;
    int xmin, xmax, ymin, ymax, n;
	BYTE *cfa = p->cfa;
	size_t cfadim = p->cfadim;

    scale2 = p->scale * p->scale;

    if (init_image_scanner(p, &s, &ymin, &ymax)) return 1;

    p->nskip = (p->ymax - p->ymin) - (ymax - ymin);
    p->nmiss = p->nskip * (p->xmax - p->xmin);

    /* This is the outer loop over all the lines in the input image */
    get_dimensions(p->output_data, osize);
    for (j = ymin; j <= ymax; ++j) {
        /* Check the overlap with the output */
        n = get_scanline_limits(&s, j, &xmin, &xmax);
        if (n == 1) {
            // scan ended (y reached the top vertex/edge)
            p->nskip += (ymax + 1 - j);
            p->nmiss += (ymax + 1 - j) * (p->xmax - p->xmin);
            break;
        } else if (n == 2 || n == 3) {
            // pixel centered on y is outside of scanner's limits or image [0, height - 1]
            // OR: limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin);
            ++p->nskip;
            continue;
        } else {
            // limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin) - (xmax + 1 - xmin);
        }

        for (i = xmin; i <= xmax; ++i) {
            float ox, oy;
            int chan = FC_array(j, i, cfa, cfadim);
            if (map_pixel(p->pixmap, i, j, &ox, &oy)) {
              ++ p->nmiss;

            } else {
                ii = fortran_round(ox);
                jj = fortran_round(oy);

                /* Check it is on the output image */
                if (ii < 0 || ii >= osize[0] || jj < 0 || jj >= osize[1]) {
                    ++ p->nmiss;

                } else {
                    vc[chan] = get_pixel(p->output_counts, ii, jj, chan);

                    /* Allow for stretching because of scale change */
                    d = get_pixel(p->data, i, j, 0) * scale2;

                    /* Scale the weighting mask by the scale factor.  Note that we
                       DON'T scale by the Jacobian as it hasn't been calculated */
                    if (p->weights) {
                        dow = get_pixel(p->weights, i, j, 0) * p->weight_scale;
                    } else {
                        dow = 1.0;
                    }

                    if (update_data(p, ii, jj, chan, d, vc[chan], dow)) {
                        return 1;
                    }
                }
            }
        }
    }

    return 0;
}

/** --------------------------------------------------------------------------------------------------
 * This kernel assumes the flux is distributed acrass a gaussian around the center of an input pixel
 *
 * p: structure containing options, input, and output
 */

static int
do_kernel_gaussian(struct driz_param_t* p) {
    struct scanner s;
    integer_t i, j, ii, jj, nxi, nxa, nyi, nya, nhit;
    integer_t /*ybounds[2],*/ osize[2];
    float vc[3], d, dow;
    float gaussian_efac, gaussian_es;
    float pfo, ac,  scale2, xxi, xxa, yyi, yya, w, ddx, ddy, r2, dover;
    const float nsig = 2.5;
    int xmin, xmax, ymin, ymax, n;
    BYTE *cfa = p->cfa;
    size_t cfadim = p->cfadim;

    /* Added in V2.9 - make sure pfo doesn't get less than 1.2
       divided by the scale so that there are never holes in the
       output */

    pfo = nsig * p->pixel_fraction / 2.3548 / p->scale;
    pfo = CLAMP_ABOVE(pfo, 1.2 / p->scale);

    ac = 1.0 / (p->pixel_fraction * p->pixel_fraction);
    scale2 = p->scale * p->scale;

    gaussian_efac = (2.3548*2.3548) * scale2 * ac / 2.0;
    gaussian_es = gaussian_efac / M_PI;

    if (init_image_scanner(p, &s, &ymin, &ymax)) return 1;

    p->nskip = (p->ymax - p->ymin) - (ymax - ymin);
    p->nmiss = p->nskip * (p->xmax - p->xmin);

    /* This is the outer loop over all the lines in the input image */

    get_dimensions(p->output_data, osize);
    for (j = ymin; j <= ymax; ++j) {
        /* Check the overlap with the output */
        n = get_scanline_limits(&s, j, &xmin, &xmax);
        if (n == 1) {
            // scan ended (y reached the top vertex/edge)
            p->nskip += (ymax + 1 - j);
            p->nmiss += (ymax + 1 - j) * (p->xmax - p->xmin);
            break;
        } else if (n == 2 || n == 3) {
            // pixel centered on y is outside of scanner's limits or image [0, height - 1]
            // OR: limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin);
            ++p->nskip;
            continue;
        } else {
            // limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin) - (xmax + 1 - xmin);
        }

        for (i = xmin; i <= xmax; ++i) {
            float ox, oy;
            int chan = FC_array(j, i, cfa, cfadim);

            if (map_pixel(p->pixmap, i, j, &ox, &oy)) {
                nhit = 0;

            } else {
                /* Offset within the subset */
                xxi = ox - pfo;
                xxa = ox + pfo;
                yyi = oy - pfo;
                yya = oy + pfo;

                nxi = MAX(fortran_round(xxi), 0);
                nxa = MIN(fortran_round(xxa), osize[0]-1);
                nyi = MAX(fortran_round(yyi), 0);
                nya = MIN(fortran_round(yya), osize[1]-1);

                nhit = 0;

                /* Allow for stretching because of scale change */
                d = get_pixel(p->data, i, j, 0) * scale2;

                /* Scale the weighting mask by the scale factor and inversely by
                   the Jacobian to ensure conservation of weight in the output */
                if (p->weights) {
                    w = get_pixel(p->weights, i, j, 0) * p->weight_scale;
                } else {
                    w = 1.0;
                }

                /* Loop over output pixels which could be affected */
                for (jj = nyi; jj <= nya; ++jj) {
                    ddy = oy - (float)jj;
                    for (ii = nxi; ii <= nxa; ++ii) {
                        ddx = ox - (float)ii;
                        /* Radial distance */
                        r2 = ddx*ddx + ddy*ddy;

                        /* Weight is a scaled Gaussian function of radial
                           distance */
                        dover = gaussian_es * exp(-r2 * gaussian_efac);

                        /* Count the hits */
                        ++nhit;

                        vc[chan] = get_pixel(p->output_counts, ii, jj, chan);
                        dow = (float)dover * w;

                        if (update_data(p, ii, jj, chan, d, vc[chan], dow)) {
                            return 1;
                        }
                    }
                }
            }

            /* Count cases where the pixel is off the output image */
            if (nhit == 0) ++ p->nmiss;
        }
    }

    return 0;
}

/** --------------------------------------------------------------------------------------------------
 * This kernel assumes flux of input pixel is distributed according to lanczos function
 *
 * p: structure containing options, input, and output
 */

static int
do_kernel_lanczos(struct driz_param_t* p) {
    struct scanner s;
    integer_t i, j, ii, jj, nxi, nxa, nyi, nya, nhit, ix, iy;
    integer_t /*ybounds[2],*/ osize[2];
    float scale2, vc[3], d, dow;
    float pfo, xx, yy, xxi, xxa, yyi, yya, w, dx, dy, dover;
    int kernel_order;
    struct lanczos_param_t lanczos;
    const size_t nlut = 512;
    const float del = 0.01;
    int xmin, xmax, ymin, ymax, n;
    BYTE *cfa = p->cfa;
    size_t cfadim = p->cfadim;

    dx = 1.0;
    dy = 1.0;

    scale2 = p->scale * p->scale;
    kernel_order = (p->kernel == kernel_lanczos2) ? 2 : 3;
    pfo = (float)kernel_order * p->pixel_fraction / p->scale;

    if ((lanczos.lut = malloc(nlut * sizeof(float))) == NULL) {
        driz_error_set_message(p->error, _("Out of memory"));
        return driz_error_is_set(p->error);
    }

    /* Set up a look-up-table for Lanczos-style interpolation
       kernels */
    create_lanczos_lut(kernel_order, nlut, del, lanczos.lut);
    lanczos.sdp = p->scale / del / p->pixel_fraction;
    lanczos.nlut = nlut;

    if (init_image_scanner(p, &s, &ymin, &ymax)) return 1;

    p->nskip = (p->ymax - p->ymin) - (ymax - ymin);
    p->nmiss = p->nskip * (p->xmax - p->xmin);

    /* This is the outer loop over all the lines in the input image */

    get_dimensions(p->output_data, osize);
    for (j = ymin; j <= ymax; ++j) {
        /* Check the overlap with the output */
        n = get_scanline_limits(&s, j, &xmin, &xmax);
        if (n == 1) {
            // scan ended (y reached the top vertex/edge)
            p->nskip += (ymax + 1 - j);
            p->nmiss += (ymax + 1 - j) * (p->xmax - p->xmin);
            break;
        } else if (n == 2 || n == 3) {
            // pixel centered on y is outside of scanner's limits or image [0, height - 1]
            // OR: limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin);
            ++p->nskip;
            continue;
        } else {
            // limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin) - (xmax + 1 - xmin);
        }

        for (i = xmin; i <= xmax; ++i) {
            int chan = FC_array(j, i, cfa, cfadim);
            if (map_pixel(p->pixmap, i, j, &xx, &yy)) {
                nhit = 0;

            } else {
                xxi = xx - dx - pfo;
                xxa = xx - dx + pfo;
                yyi = yy - dy - pfo;
                yya = yy - dy + pfo;

                nxi = MAX(fortran_round(xxi), 0);
                nxa = MIN(fortran_round(xxa), osize[0]-1);
                nyi = MAX(fortran_round(yyi), 0);
                nya = MIN(fortran_round(yya), osize[1]-1);

                nhit = 0;

                /* Allow for stretching because of scale change */
                d = get_pixel(p->data, i, j, 0) * scale2;

                /* Scale the weighting mask by the scale factor and inversely by
                   the Jacobian to ensure conservation of weight in the output */
                if (p->weights) {
                    w = get_pixel(p->weights, i, j, 0) * p->weight_scale;
                } else {
                    w = 1.0;
                }

                /* Loop over output pixels which could be affected */
                for (jj = nyi; jj <= nya; ++jj) {
                    for (ii = nxi; ii <= nxa; ++ii) {
                        /* X and Y offsets */
                        ix = fortran_round(fabs(xx - (float)ii) * lanczos.sdp) + 1;
                        iy = fortran_round(fabs(yy - (float)jj) * lanczos.sdp) + 1;

                        /* Weight is product of Lanczos function values in X and Y */
                        dover = lanczos.lut[ix] * lanczos.lut[iy];

                        /* Count the hits */
                        ++nhit;

                        vc[chan] = get_pixel(p->output_counts, ii, jj, chan);
                        dow = (float)(dover * w);

                        if (update_data(p, ii, jj, chan, d, vc[chan], dow)) {
                            return 1;
                        }
                    }
                }
            }

            /* Count cases where the pixel is off the output image */
            if (nhit == 0) ++ p->nmiss;
        }
    }

    free(lanczos.lut);
    lanczos.lut = NULL;

    return 0;
}

/** --------------------------------------------------------------------------------------------------
 * This kernel assumes the input flux is evenly distributed over a rectangle whose sides are
 * aligned with the ouput pixel. Called turbo because it is fast, but approximate.
 *
 * p: structure containing options, input, and output
 */

static int
do_kernel_turbo(struct driz_param_t* p) {
    struct scanner s;
    integer_t i, j, ii, jj, nxi, nxa, nyi, nya, nhit, iis, iie, jjs, jje;
    integer_t osize[2];
    float vc[3], d, dow;
    float pfo, scale2, ac;
    float xxi, xxa, yyi, yya, w, dover;
    int xmin, xmax, ymin, ymax, n;
    BYTE *cfa = p->cfa;
    size_t cfadim = p->cfadim;

    siril_debug_print("starting do_kernel_turbo\n");
    ac = 1.0 / (p->pixel_fraction * p->pixel_fraction);
    pfo = p->pixel_fraction / p->scale / 2.0;
    scale2 = p->scale * p->scale;

    if (init_image_scanner(p, &s, &ymin, &ymax)) return 1;

    p->nskip = (p->ymax - p->ymin) - (ymax - ymin);
    p->nmiss = p->nskip * (p->xmax - p->xmin);

    /* This is the outer loop over all the lines in the input image */

    get_dimensions(p->output_data, osize);
    for (j = ymin; j <= ymax; ++j) {
        /* Check the overlap with the output */
        n = get_scanline_limits(&s, j, &xmin, &xmax);

        if (n == 1) {
            // scan ended (y reached the top vertex/edge)
            p->nskip += (ymax + 1 - j);
            p->nmiss += (ymax + 1 - j) * (p->xmax - p->xmin);
            break;
        } else if (n == 2 || n == 3) {
            // pixel centered on y is outside of scanner's limits or image [0, height - 1]
            // OR: limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin);
            ++p->nskip;
            continue;
        } else {
            // limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin) - (xmax + 1 - xmin);
        }

        for (i = xmin; i <= xmax; ++i) {
            float ox, oy;
            int chan = FC_array(j, i, cfa, cfadim);

            if (map_pixel(p->pixmap, i, j, &ox, &oy)) {
                nhit = 0;

            } else {
                /* Offset within the subset */
                xxi = ox - pfo;
                xxa = ox + pfo;
                yyi = oy - pfo;
                yya = oy + pfo;

                nxi = fortran_round(xxi);
                nxa = fortran_round(xxa);
                nyi = fortran_round(yyi);
                nya = fortran_round(yya);
                iis = MAX(nxi, 0);  /* Needed to be set to 0 to avoid edge effects */
                iie = MIN(nxa, osize[0]-1);
                jjs = MAX(nyi, 0);  /* Needed to be set to 0 to avoid edge effects */
                jje = MIN(nya, osize[1]-1);

                nhit = 0;

                /* Allow for stretching because of scale change */
                d = get_pixel(p->data, i, j, 0) * (float)scale2;

                /* Scale the weighting mask by the scale factor and inversely by
                   the Jacobian to ensure conservation of weight in the output. */
                if (p->weights) {
                    w = get_pixel(p->weights, i, j, 0) * p->weight_scale;
                } else {
                    w = 1.0;
                }

                /* Loop over the output pixels which could be affected */
                for (jj = jjs; jj <= jje; ++jj) {
                    for (ii = iis; ii <= iie; ++ii) {
                        /* Calculate the overlap using the simpler "aligned" box
                           routine */
                        dover = over(ii, jj, xxi, xxa, yyi, yya);

                        if (dover > 0.0) {
                            /* Correct for the pixfrac area factor */
                            dover *= scale2 * ac;

                            /* Count the hits */
                            ++nhit;

                            vc[chan] = get_pixel(p->output_counts, ii, jj, chan);
                            dow = (float)(dover * w);

                            if (update_data(p, ii, jj, chan, d, vc[chan], dow)) {
                                return 1;
                            }
                        }
                    }
                }
            }

            /* Count cases where the pixel is off the output image */
            if (nhit == 0) ++ p->nmiss;
        }
    }

    siril_debug_print("ending do_kernel_turbo\n");
    return 0;
}

/** --------------------------------------------------------------------------------------------------
 * This module does the actual mapping of input flux to output images. It works by calculating the
 * positions of the four corners of a quadrilateral on the output grid corresponding to the corners
 * of the input pixel and then working out exactly how much of each pixel in the output is covered, or not.
 *
 * p: structure containing options, input, and output
 */

int
do_kernel_square(struct driz_param_t* p) {
    integer_t i, j, ii, jj, min_ii, max_ii, min_jj, max_jj, nhit;
    integer_t osize[2];
    float scale2, vc[3], d, dow;
    float dh, jaco, tem, dover, w;
    float xin[4], yin[4], xout[4], yout[4];
    struct scanner s;
    int xmin, xmax, ymin, ymax, n;
    BYTE *cfa = p->cfa;
    size_t cfadim = p->cfadim;
	integer_t maxarea = 0, mnii = 0, mxii = 0, mnjj = 0, mxjj = 0;

    siril_debug_print("starting do_kernel_square\n");
    dh = 0.5 * p->pixel_fraction;
    scale2 = p->scale * p->scale;

    /* Next the "classic" drizzle square kernel...  this is different
       because we have to transform all four corners of the shrunken
       pixel */
    if (init_image_scanner(p, &s, &ymin, &ymax)) return 1;

    p->nskip = (p->ymax - p->ymin) - (ymax - ymin);
    p->nmiss = p->nskip * (p->xmax - p->xmin);

    /* This is the outer loop over all the lines in the input image */
    get_dimensions(p->output_data, osize);
    for (j = ymin; j <= ymax; ++j) {
        /* Check the overlap with the output */
        n = get_scanline_limits(&s, j, &xmin, &xmax);
        if (n == 1) {
            // scan ended (y reached the top vertex/edge)
            p->nskip += (ymax + 1 - j);
            p->nmiss += (ymax + 1 - j) * (p->xmax - p->xmin);
            break;
        } else if (n == 2 || n == 3) {
            // pixel centered on y is outside of scanner's limits or image [0, height - 1]
            // OR: limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin);
            ++p->nskip;
            continue;
        } else {
            // limits (x1, x2) are equal (line width is 0)
            p->nmiss += (p->xmax - p->xmin) - (xmax + 1 - xmin);
        }

        /* Set the input corner positions */

        yin[1] = yin[0] = (float) j + dh;
        yin[3] = yin[2] = (float) j - dh;

        for (i = xmin; i <= xmax; ++i) {
            nhit = 0;
            int chan = FC_array(j, i, cfa, cfadim);

            xin[3] = xin[0] = (float) i - dh;
            xin[2] = xin[1] = (float) i + dh;

            for (ii = 0; ii < 4; ++ii) {
                if (interpolate_point(p, xin[ii], yin[ii], xout + ii, yout + ii)) {
                    goto _miss;
                }
            }

            /* Work out the area of the quadrilateral on the output grid.
               Note that this expression expects the points to be in clockwise
               order */

            jaco = 0.5f * ((xout[1] - xout[3]) * (yout[0] - yout[2]) -
                           (xout[0] - xout[2]) * (yout[1] - yout[3]));

            if (jaco < 0.0) {
                jaco *= -1.0;
                /* Swap */
                tem = xout[1]; xout[1] = xout[3]; xout[3] = tem;
                tem = yout[1]; yout[1] = yout[3]; yout[3] = tem;
            }

            /* Allow for stretching because of scale change */
            d = get_pixel(p->data, i, j, 0) * scale2;

            /* Scale the weighting mask by the scale factor and inversely by
               the Jacobian to ensure conservation of weight in the output */
            if (p->weights) {
                w = get_pixel(p->weights, i, j,0) * p->weight_scale;
            } else {
                w = 1.0;
            }

            /* Loop over output pixels which could be affected */
            min_jj = MAX(fortran_round(min_floats(yout, 4)), 0);
            max_jj = MIN(fortran_round(max_floats(yout, 4)), osize[1]-1);
            min_ii = MAX(fortran_round(min_floats(xout, 4)), 0);
            max_ii = MIN(fortran_round(max_floats(xout, 4)), osize[0]-1);
			int area = (max_jj - min_jj) - (max_ii - min_ii);
			if (area > maxarea) {
				maxarea = area;
				mnjj = min_jj;
				mxjj = max_jj;
				mnii = min_ii;
				mxii = max_ii;
			}

            for (jj = min_jj; jj <= max_jj; ++jj) {
                for (ii = min_ii; ii <= max_ii; ++ii) {
                    /* Call compute_area to calculate overlap */
                    dover = compute_area((float)ii, (float)jj, xout, yout);

                    if (dover > 0.0) {
                        vc[chan] = get_pixel(p->output_counts, ii, jj, chan);

                        /* Re-normalise the area overlap using the Jacobian */
                        dover /= jaco;
                        dow = (float)(dover * w);

                        /* Count the hits */
                        ++nhit;

                        if (update_data(p, ii, jj, chan, d, vc[chan], dow)) {
                            return 1;
                        }
                    }
                }
            }

            /* Count cases where the pixel is off the output image */
            _miss:
            if (nhit == 0) {
                ++ p->nmiss;
            }
        }
    }

    printf("do_square max area: %d. (%d, %d) to (%d, %d)\n", maxarea, mnii, mnjj, mxii, mxjj);
    siril_debug_print("ending do_kernel_square\n");
    return 0;
}

/** --------------------------------------------------------------------------------------------------
 * The user selects a kernel to use for drizzling from a function in the following tables
 * The kernels differ in how the flux inside a single pixel is allocated: evenly spread
 * across the pixel, concentrated at the central point, or by some other function.
 */

static kernel_handler_t
kernel_handler_map[] = {
    do_kernel_square,
    do_kernel_gaussian,
    do_kernel_point,
    do_kernel_turbo,
    do_kernel_lanczos,
    do_kernel_lanczos
};

/** --------------------------------------------------------------------------------------------------
 * The executive function which calls the kernel which does the actual drizzling
 *
 * p: structure containing options, input, and output
 */

int
dobox(struct driz_param_t* p) {
    kernel_handler_t kernel_handler = NULL;
    siril_debug_print("starting dobox\n");

    /* Set up a function pointer to handle the appropriate kernel */
    if (p->kernel < kernel_LAST) {
        kernel_handler = kernel_handler_map[p->kernel];

        if (kernel_handler != NULL) {
            kernel_handler(p);
        }
    }

    if (kernel_handler == NULL) {
        driz_error_set_message(p->error, _("Invalid kernel type"));
    }

    siril_debug_print("ending dobox\n");
    return driz_error_is_set(p->error);
}
