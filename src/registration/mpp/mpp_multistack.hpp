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

/*
 * Multi-source derotation stack (Option B, C5+C6+C7). Combines several
 * sequences — each already derotated to a COMMON epoch via its own .derot — into
 * one stack in a designated reference sequence's epoch canvas. The standard
 * streamed MPP engines (global align, average reference, AP qualities, per-AP
 * shifts, warp stack) all consume a FrameProvider over one contiguous index
 * space, so the only multi-source-specific machinery is:
 *   - a global index space over all the sequences' frames (MultiSourceLayout);
 *   - analysis providers that read each source frame and RELOCATE+derotate it
 *     into the reference canvas via the C3 ms map, so the shared reference and
 *     the per-AP residual shifts are all measured consistently;
 *   - a warp-stack DerotMapProvider that maps each output pixel back to the
 *     right source sequence's native frame (again via the ms map).
 *
 * The core (multistack_channel) is provider-injected: callers supply per-source
 * frame readers, so it runs identically on real sequences (the CLI/GUI path)
 * and on synthetic in-memory frames (the unit test). No I/O or Siril sequence
 * types appear here.
 */

#ifndef SRC_REGISTRATION_MPP_MULTISTACK_HPP_
#define SRC_REGISTRATION_MPP_MULTISTACK_HPP_

#include <functional>
#include <vector>

#include <opencv2/core.hpp>

#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_derot_sidecar.h"   /* mpp_derot_t */

namespace mpp {

/* One source sequence in a multi-source combine. `derot` is this sequence's
 * .derot plan (all sources share one epoch). The two readers return a frame in
 * THIS sequence's own native pixel coordinates (the same coordinates its disk
 * fit in `derot` is expressed in):
 *   analysis_read(local) -> raw mono frame   (luminance; for align / AP / shift)
 *   full_read(local)     -> full frame        (1 or 3 layers; for the stack)
 * An empty Mat signals a read miss (the frame is then excluded). */
struct MsSource {
	const mpp_derot_t *derot = nullptr;
	int num_frames = 0;
	std::function<cv::Mat(int local)> analysis_read;
	std::function<cv::Mat(int local)> full_read;
};

struct MultiStackResult {
	cv::Mat image;            /* CV_16U / CV_16UC3; empty on cancel/error */
	cv::Vec4i intersection;   /* (y_low, y_high, x_low, x_high) in the reference
	                           * canvas — the image's placement (pre-drizzle) */
	bool oom = false;
	bool cancelled = false;
	bool error = false;
	int num_aps = 0;
	int num_frames = 0;       /* union total actually combined */
};

/* Combine the sources into one stack in the epoch canvas defined by `out_derot`
 * (its disk fit gives the output geometry, its frame_rows/cols the size). The
 * output canvas is decoupled from the source list so that several channels
 * (R/G/B) can each stack their OWN sequences into the SAME reference canvas and
 * stay mutually co-registered — `out_derot` is the user-designated reference
 * sequence's plan and need not be one of `srcs`. `num_layers` is the output
 * channel count (1 or 3 — full_read must match). Reuses the standard streamed
 * engines; every frame is relocated + derotated into the reference canvas. With
 * a single source whose derot is `out_derot` this reduces to the ordinary
 * single-sequence derotation stack. */
MultiStackResult multistack_channel(const std::vector<MsSource> &srcs,
                                    const mpp_derot_t *out_derot, int num_layers,
                                    const mpp_config_t &cfg,
                                    int max_threads = 1);

}  // namespace mpp

#endif /* SRC_REGISTRATION_MPP_MULTISTACK_HPP_ */
