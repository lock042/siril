/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#ifndef SRC_REGISTRATION_MPP_DEROT_H_
#define SRC_REGISTRATION_MPP_DEROT_H_

#include <opencv2/core.hpp>

#include "registration/mpp/mpp_derot_sidecar.h"   /* mpp_derot_t */

/* Bridge between the .derot plan (mpp_derot_t: per-frame B/CM/P + disk fit, in
 * original full-frame pixel coordinates) and the geometry projection
 * (mpp_derot_geom). Builds the per-frame derotation remap into a target canvas
 * whose pixel (0,0) corresponds to original-frame coordinate (org_x, org_y)
 * and which is scaled by `scale` (the drizzle factor at stack time, 1.0 at
 * registration). Fills:
 *   mapx/mapy (CV_32F) — source position in the (MPP-aligned, same-scale) frame
 *                        canvas of the surface point seen at each output pixel;
 *   mu (CV_32F)        — cosine of the emission angle at the frame's time, 0
 *                        off-disk.
 * These are exactly the inputs the warp engine's DerotMapProvider expects (the
 * engine adds the global offset and AP-residual terms). */
void mpp_derot_frame_map(const mpp_derot_t *d, int frame,
                         int out_w, int out_h,
                         double org_x, double org_y, double scale,
                         cv::Mat &mapx, cv::Mat &mapy, cv::Mat &mu);

#endif /* SRC_REGISTRATION_MPP_DEROT_H_ */
