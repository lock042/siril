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

#ifndef CDRIZZLEMAP_H
#define CDRIZZLEMAP_H

#include "driz_portability.h"
#include "cdrizzleutil.h"
#include <wcslib.h>

/* Line segment structure, used for computing overlap
 * The first index on line is the endpoint
 * The second index is {x, y) coordinate of the point
 * The valid flag is non-zero if it does not intersect the image
 */

struct segment {
    float  point[2][2];
    int     invalid;
};

// IMAGE_OUTLINE_NPTS - maximum number of vertices in the bounding polygon
// for input and resampled images
#define IMAGE_OUTLINE_NPTS 4

/** vertex structure.
 *
 *  This structure holds coordinates of a polygon vertex.
 *
 */
struct vertex {
    float x; /**< x-coordinate */
    float y; /**< y-coordinate */
};

/** polygon structure.
 *
 *  This structure holds information about polygon vertices. The maximum number
 *  of vertices that this structure can hold is determined by the constant
 *  IMAGE_OUTLINE_NPTS and it is float of IMAGE_OUTLINE_NPTS value.
 *
 *  NOTE: polygons must not be closed (that is last vertex != first vertex).
 *
 */
struct polygon {
    struct vertex v[2 * IMAGE_OUTLINE_NPTS];  /**< polygon vertices */
    int    npv;  /**< actual number of polygon vertices */
};

/** edge structure.
 *
 *  This structure holds information about vertices of the edge, edge position
 *  in the polygon (left or right of the line going through the top and bottom
 *  vertices), slope, and interceipts.
 *
 */
struct edge {
    struct vertex v1; /**< first vertex */
    struct vertex v2; /**< second vertex */
    float m; /**< edge's slope */
    float b; /**< edge's intercept */
    float c; /**< modified intercept */
    int p;  /**< edge's position: -1 for left-side edge and +1 for right-side edge */
};

/** scanner structure.
 *
 *  This structure holds information needed to "scan" or rasterize the
 *  intersection polygon formed by intersecting the bounding box of the output
 *  frame with the bounding box of the input frame in the input frame's
 *  coordinate system.
 *
 */
struct scanner {
    struct edge left_edges[2 * IMAGE_OUTLINE_NPTS]; /**< left edges */
    struct edge right_edges[2 * IMAGE_OUTLINE_NPTS]; /**< right edges */
    struct edge *left; /**< pointer to the current left edge; NULL when top polygon vertex was reached */
    struct edge *right; /**< pointer to the current right edge; NULL when top polygon vertex was reached */
    int nleft; /**< number of left edges */
    int nright; /**< number of right edges */
    float min_y; /**< minimum y-coordinate of all polygon vertices */
    float max_y; /**< maximum y-coordinate of all polygon vertices */
    int xmin; /**< min valid pixels' x-coord in pixmap (from bounding box carried over from driz_param_t) rounded to int */
    int xmax; /**< max valid pixels' x-coord in pixmap (from bounding box carried over from driz_param_t) rounded to int */
    int ymin; /**< min valid pixels' y-coord in pixmap (from bounding box carried over from driz_param_t) rounded to int */
    int ymax; /**< max valid pixels' y-coord in pixmap (from bounding box carried over from driz_param_t) rounded to int */
    // overlap_valid: 1 if polygon intersection and coord inversion worked;
    //                0 if computation of xmin, xmax, ymin, ymax has
    //                  failed in which case they are carried over from driz_param_t.
    int overlap_valid; /**< 1 if x/y min/max updated from polygon intersection and 0 if carried over from driz_param_t */
};

int
interpolate_point(struct driz_param_t *par, float xin, float yin,
                  float *xout, float *yout);

// int map_image_coordinates_wcs(int width, int height, struct wcsprm *wcs_i, struct wcsprm *wcs_o, imgmap_t *p, float scale);

int map_image_coordinates_h(fits *fit, Homography H, imgmap_t *p, int target_rx, int target_ry, float scale, disto_data *disto, int threads);

int
map_point(struct driz_param_t *par, float xin, float yin,
          float *xout, float *yout);

int
map_pixel(imgmap_t *pixmap, int i, int j, float *x, float *y);

int
shrink_image_section(fits *pixmap, int *xmin, int *xmax,
                     int *ymin, int *ymax);

int
invert_pixmap(struct driz_param_t* par, float xout, float yout,
              float *xin, float *yin);

int
clip_polygon_to_window(struct polygon *p, struct polygon *wnd,
                       struct polygon *cp);

int
init_scanner(struct polygon *p, struct driz_param_t* par, struct scanner *s);

int
get_scanline_limits(struct scanner *s, int y, int *x1, int *x2);

int
init_image_scanner(struct driz_param_t* par, struct scanner *s,
                   int *ymin, int *ymax);

#endif /* CDRIZZLEMAP_H */
