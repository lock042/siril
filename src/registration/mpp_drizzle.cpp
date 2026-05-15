/*
 * mpp_drizzle — Phase 5b STScI drizzle backend for MPP stacking.
 *
 * Slice 5b.1 (this file): pixmap builder. Encodes the full MPP warp
 * (global frame shift + smooth interpolation of per-AP local shifts) as
 * a per-pixel output-coordinate map suitable for feeding to dobox().
 *
 * The actual dobox integration (5b.2) and Bayer-drizzle variant (5b.3)
 * land in follow-up slices and reuse the pixmap built here.
 */

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_drizzle.h"
#include "registration/mpp_shift.h"
}

#include "registration/mpp_stack_priv.hpp"   /* mpp::stack_one_dim_weight */

extern "C" mpp_status_t mpp_imgmap_alloc(imgmap_t *out, int rx, int ry) {
	if (!out || rx <= 0 || ry <= 0) return MPP_EINVAL;
	const size_t n = (size_t) rx * (size_t) ry;
	out->xmap = (float *) std::malloc(n * sizeof(float));
	if (!out->xmap) {
		out->ymap = nullptr;
		out->rx = 0;
		out->ry = 0;
		return MPP_ENOMEM;
	}
	out->ymap = (float *) std::malloc(n * sizeof(float));
	if (!out->ymap) {
		std::free(out->xmap);
		out->xmap = nullptr;
		out->rx = 0;
		out->ry = 0;
		return MPP_ENOMEM;
	}
	out->rx = rx;
	out->ry = ry;
	return MPP_OK;
}

extern "C" void mpp_imgmap_free(imgmap_t *m) {
	if (!m) return;
	std::free(m->xmap);
	std::free(m->ymap);
	m->xmap = nullptr;
	m->ymap = nullptr;
	m->rx = 0;
	m->ry = 0;
}

extern "C" mpp_status_t mpp_pixmap_build(const mpp_run_t *run,
                                         int frame_idx,
                                         double drizzle_scale,
                                         imgmap_t *out) {
	if (!run || !out || !out->xmap || !out->ymap) return MPP_EINVAL;
	if (!run->aps || !run->global_shifts || !run->shifts || !run->shifts->shifts)
		return MPP_EINVAL;
	if (frame_idx < 0 || frame_idx >= run->num_frames) return MPP_EINVAL;
	if (frame_idx >= run->shifts->num_frames) return MPP_EINVAL;
	if (out->rx != run->frame_cols || out->ry != run->frame_rows) return MPP_EINVAL;
	if (drizzle_scale <= 0.0) return MPP_EINVAL;

	const int rx = out->rx;
	const int ry = out->ry;
	const int M  = run->aps->count;

	const double gdy = (double) run->global_shifts[2 * frame_idx + 0];
	const double gdx = (double) run->global_shifts[2 * frame_idx + 1];

	/* Mean-frame origin in the original frame's coordinate system. */
	const double iy_lo = (double) run->intersection[0];
	const double ix_lo = (double) run->intersection[2];

	/* Per-AP 1D weight ramps, precomputed once for this frame and indexed
	 * directly per pixel inside the inner loop. The ramps depend only on
	 * patch bounds + AP centre + frame-border extend flags, all of which
	 * are constant across frames — but recomputing per frame is cheap
	 * (M * (patch_h + patch_w) float ops) and keeps the API surface to
	 * a single function call rather than a setup/build/teardown trio.
	 *
	 * Suppress empty-ramp pitfalls: stack_one_dim_weight handles
	 * zero-length cleanly so AP boxes degenerate to no influence rather
	 * than crashing.
	 */
	std::vector<std::vector<float>> wy(M), wx(M);
	for (int a = 0; a < M; ++a) {
		const auto &ap = run->aps->records[a];
		const bool extend_y_low  = (ap.patch_y_low  == 0);
		const bool extend_y_high = (ap.patch_y_high == ry);
		const bool extend_x_low  = (ap.patch_x_low  == 0);
		const bool extend_x_high = (ap.patch_x_high == rx);
		wy[a] = mpp::stack_one_dim_weight(ap.patch_y_low, ap.patch_y_high,
		                                  ap.y, extend_y_low, extend_y_high);
		wx[a] = mpp::stack_one_dim_weight(ap.patch_x_low, ap.patch_x_high,
		                                  ap.x, extend_x_low, extend_x_high);
	}

	/* Per-frame per-AP shifts: shifts[2 * (frame*M + ap) + {0=dy, 1=dx}]. */
	const size_t frame_base = (size_t) frame_idx * (size_t) M;
	const double *ap_shifts = &run->shifts->shifts[2 * frame_base];

	/* Per-pixel loop. collapse(2) parallelises over both axes so a single
	 * pixmap build saturates the available threads on small frames as
	 * well as large ones. Each (j, i) writes independently to xmap[j*rx+i]
	 * / ymap[j*rx+i] so there's no contention. */
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
	for (int j = 0; j < ry; ++j) {
		for (int i = 0; i < rx; ++i) {
			double wsum = 0.0, wsy = 0.0, wsx = 0.0;
			for (int a = 0; a < M; ++a) {
				const auto &ap = run->aps->records[a];
				if (j < ap.patch_y_low || j >= ap.patch_y_high) continue;
				if (i < ap.patch_x_low || i >= ap.patch_x_high) continue;
				const float wy_v = wy[a][j - ap.patch_y_low];
				const float wx_v = wx[a][i - ap.patch_x_low];
				const double w   = (double) std::min(wy_v, wx_v);
				if (w <= 0.0) continue;
				wsum += w;
				wsy  += w * ap_shifts[2 * a + 0];
				wsx  += w * ap_shifts[2 * a + 1];
			}
			double sy = gdy;
			double sx = gdx;
			if (wsum > 0.0) {
				sy += wsy / wsum;
				sx += wsx / wsum;
			}
			const size_t idx = (size_t) j * (size_t) rx + (size_t) i;
			out->xmap[idx] = (float) (drizzle_scale * ((double) i + sx - ix_lo));
			out->ymap[idx] = (float) (drizzle_scale * ((double) j + sy - iy_lo));
		}
	}
	return MPP_OK;
}
