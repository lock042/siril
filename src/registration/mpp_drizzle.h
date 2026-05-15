#ifndef SRC_REGISTRATION_MPP_DRIZZLE_H_
#define SRC_REGISTRATION_MPP_DRIZZLE_H_

#include "registration/mpp.h"
#include "drizzle/imgmap.h"   /* imgmap_t — minimal, C++-clean */

#ifdef __cplusplus
extern "C" {
#endif

/* -------- Pixmap (per-frame output-coordinate map) -------- */
/*
 * The Phase 5b STScI drizzle backend feeds dobox() per-frame pixmaps that
 * encode the full MPP warp: global frame shift + smooth interpolation of
 * the per-AP local shifts at each input pixel. dobox doesn't need to know
 * anything about APs — the pixmap is the only place that knows.
 *
 * One buffer is allocated for the whole stack run and overwritten in
 * place by every frame; see Phase 5b in pss_port_plan.md for the
 * lifetime / parallelism rationale.
 */

/* Allocate xmap/ymap as rx*ry float arrays and set out->rx, out->ry.
 * Pre-existing buffers in `out` are leaked — call mpp_imgmap_free first
 * if reusing the struct. Returns MPP_ENOMEM on allocation failure (in
 * which case both pointers are NULL and rx/ry are 0). */
mpp_status_t mpp_imgmap_alloc(imgmap_t *out, int rx, int ry);

/* Free xmap and ymap and zero the struct. Safe on an already-zeroed or
 * already-freed `m` (no-ops on NULL pointers via free()). */
void mpp_imgmap_free(imgmap_t *m);

/* Fill `out`'s xmap/ymap for input-frame `frame_idx`.
 *
 * Layout: row-major over the input frame. xmap[j*rx + i] is the column
 * (x) on the *drizzled output canvas* where input pixel (i, j) should
 * land; ymap[j*rx + i] is the row.
 *
 * Coordinate convention. Phase 2's `align_average_frame` defined the
 * intersection (mean-frame bounds) as a Vec4i (y_low, y_high, x_low,
 * x_high) in the original frame coordinate system, and aligned frame f
 * pixel (j_f, i_f) corresponds to mean-frame pixel
 *     (j_f + shift_y_f - intersection_y_low,
 *      i_f + shift_x_f - intersection_x_low).
 * The drizzled output canvas is the mean frame scaled by drizzle_scale.
 * So:
 *     xmap[j, i] = drizzle_scale * (i + sx - intersection_x_low)
 *     ymap[j, i] = drizzle_scale * (j + sy - intersection_y_low)
 * where sx, sy are the total per-pixel shift on this frame:
 *     sx = global_shifts[frame_idx].dx + ap_local_dx(i, j, frame_idx)
 *     sy = global_shifts[frame_idx].dy + ap_local_dy(i, j, frame_idx)
 * and ap_local_{dy, dx} is the weighted average of the per-AP shifts at
 * (i, j), with weights from PSS's `one_dim_weight` ramp on each axis,
 * combined with `min(wy, wx)` (the same convention Phase 5a uses for
 * per-AP buffer weighting). Pixels outside every AP's patch fall back
 * to the global shift only.
 *
 * `drizzle_scale` is a positive number, typically 1.0, 2.0, or 3.0.
 * Non-integer scales are accepted for future support but unused by the
 * current Phase 5b paths.
 *
 * Requires: run->aps, run->global_shifts, run->shifts (all with
 * matching dimensions); out pre-allocated to (run->frame_cols,
 * run->frame_rows). Returns MPP_EINVAL on inconsistent input. */
mpp_status_t mpp_pixmap_build(const mpp_run_t *run,
                              int frame_idx,
                              double drizzle_scale,
                              imgmap_t *out);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_DRIZZLE_H_ */
