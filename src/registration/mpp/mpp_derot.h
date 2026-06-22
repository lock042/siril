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
 *   mapx/mapy (CV_32F) — source position to sample: the derotated globe point,
 *                        or the identity for pixels outside the fitted globe
 *                        (rings/sky pass through and stack normally);
 *   mu (CV_32F)        — 1 where usable, 0 only where a globe point rotated
 *                        behind the limb (masked).
 * These are exactly the inputs the warp engine's DerotMapProvider expects (the
 * engine adds the global offset and AP-residual terms). */
void mpp_derot_frame_map(const mpp_derot_t *d, int frame,
                         int out_w, int out_h,
                         double org_x, double org_y, double scale,
                         cv::Mat &mapx, cv::Mat &mapy, cv::Mat &mu);

/* Multi-source variant: the output canvas globe is the reference sequence's
 * epoch disk (out_d, placed by out_org/out_scale), while the source frame
 * belongs to a different sequence (src_d, frame `frame`, placed by src_org/
 * src_scale and oriented by its own ephemeris). One remap derotates the globe
 * to the epoch and relocates it from this sequence's frame into the common
 * reference canvas, so several sequences (e.g. R/G/B, or many of one filter)
 * stack as one. With out_d == src_d and equal placement this is exactly
 * mpp_derot_frame_map. */
void mpp_derot_frame_map_ms(const mpp_derot_t *out_d, const mpp_derot_t *src_d,
                            int frame, int out_w, int out_h,
                            double out_org_x, double out_org_y, double out_scale,
                            double src_org_x, double src_org_y, double src_scale,
                            cv::Mat &mapx, cv::Mat &mapy, cv::Mat &mu,
                            bool mask_outside = false);

#endif /* SRC_REGISTRATION_MPP_DEROT_H_ */
