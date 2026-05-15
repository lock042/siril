/*
 * Criterion tests for mpp_drizzle.
 *
 * Phase 5b slice 1: pixmap builder unit tests. Six small isolation tests
 * that prove the pixmap is what it claims before any dobox integration
 * (which lands in slice 5b.2 and will exercise the pixmap through a
 * real drizzle accumulation).
 *
 * The tests build minimal mpp_run_t fixtures by hand — no PSS oracle —
 * because the pixmap is a pure math function with no algorithm-faithful
 * reference. The point of these tests is "did we encode the coordinate
 * mapping correctly", not "do we match PSS".
 */

#include <criterion/criterion.h>
#include <cmath>
#include <cstdlib>
#include <cstring>

extern "C" {
#include "core/siril.h"
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_drizzle.h"
#include "registration/mpp_shift.h"
}

/* Slice 5b.2 pulls Stage-C orchestration into mpp_drizzle.cpp, which
 * transitively references `cominfo com` (gui_iface, threading prefs,
 * etc.). Tests linking against the new TU need a stub — same pattern
 * as imoper_test.c / stacking_blocks_test.c. */
extern "C" cominfo com;
extern "C" fits *gfit;
cominfo com;
fits   *gfit = nullptr;

namespace {

/* Minimal fixture: a run with frame dims (rx, ry), a single AP placed at
 * (ay, ax) with the given patch and box bounds, and explicit shifts.
 * The intersection is (0, ry, 0, rx) so identity maps without offset.
 *
 * Caller owns the returned run and must mpp_run_free it.
 */
mpp_run_t *make_fixture(int rx, int ry,
                        int num_frames,
                        int ay, int ax,
                        int patch_y_low, int patch_y_high,
                        int patch_x_low, int patch_x_high,
                        const int *global_dy_dx /* length 2*N or NULL */,
                        const double *ap_dy_dx  /* length 2*N or NULL */) {
	mpp_run_t *run = mpp_run_alloc();
	cr_assert_not_null(run);
	run->num_frames = num_frames;
	run->frame_rows = ry;
	run->frame_cols = rx;
	run->num_layers = 1;
	run->bitdepth   = 16;
	run->intersection[0] = 0;
	run->intersection[1] = ry;
	run->intersection[2] = 0;
	run->intersection[3] = rx;

	run->global_shifts = (int *) std::calloc(2 * num_frames, sizeof(int));
	if (global_dy_dx)
		std::memcpy(run->global_shifts, global_dy_dx,
		            2 * num_frames * sizeof(int));

	run->aps = (mpp_aps_t *) std::calloc(1, sizeof(mpp_aps_t));
	run->aps->count = 1;
	run->aps->records = (mpp_ap_record_t *) std::calloc(1, sizeof(mpp_ap_record_t));
	run->aps->records[0].y = ay;
	run->aps->records[0].x = ax;
	run->aps->records[0].patch_y_low  = patch_y_low;
	run->aps->records[0].patch_y_high = patch_y_high;
	run->aps->records[0].patch_x_low  = patch_x_low;
	run->aps->records[0].patch_x_high = patch_x_high;
	run->aps->records[0].box_y_low    = patch_y_low;
	run->aps->records[0].box_y_high   = patch_y_high;
	run->aps->records[0].box_x_low    = patch_x_low;
	run->aps->records[0].box_x_high   = patch_x_high;
	run->aps->records[0].structure    = 1.0;

	run->shifts = (mpp_shifts_t *) std::calloc(1, sizeof(mpp_shifts_t));
	run->shifts->num_frames = num_frames;
	run->shifts->num_aps    = 1;
	run->shifts->shifts = (double *) std::calloc(2 * num_frames * 1, sizeof(double));
	run->shifts->success = (uint8_t *) std::calloc(num_frames * 1, sizeof(uint8_t));
	for (int f = 0; f < num_frames; ++f) run->shifts->success[f] = 1;
	if (ap_dy_dx)
		std::memcpy(run->shifts->shifts, ap_dy_dx,
		            2 * num_frames * 1 * sizeof(double));
	return run;
}

/* Two-AP fixture for the between-APs interpolation test. */
mpp_run_t *make_two_ap_fixture(int rx, int ry,
                               int ay1, int ax1, int ay2, int ax2,
                               int patch_half_w, int patch_half_h,
                               double ap1_dy, double ap1_dx,
                               double ap2_dy, double ap2_dx) {
	mpp_run_t *run = mpp_run_alloc();
	cr_assert_not_null(run);
	run->num_frames = 1;
	run->frame_rows = ry;
	run->frame_cols = rx;
	run->num_layers = 1;
	run->bitdepth   = 16;
	run->intersection[0] = 0;
	run->intersection[1] = ry;
	run->intersection[2] = 0;
	run->intersection[3] = rx;

	run->global_shifts = (int *) std::calloc(2, sizeof(int));   /* zero */

	run->aps = (mpp_aps_t *) std::calloc(1, sizeof(mpp_aps_t));
	run->aps->count = 2;
	run->aps->records = (mpp_ap_record_t *) std::calloc(2, sizeof(mpp_ap_record_t));
	const int aps_y[2] = { ay1, ay2 };
	const int aps_x[2] = { ax1, ax2 };
	for (int a = 0; a < 2; ++a) {
		run->aps->records[a].y = aps_y[a];
		run->aps->records[a].x = aps_x[a];
		run->aps->records[a].patch_y_low  = std::max(0,  aps_y[a] - patch_half_h);
		run->aps->records[a].patch_y_high = std::min(ry, aps_y[a] + patch_half_h);
		run->aps->records[a].patch_x_low  = std::max(0,  aps_x[a] - patch_half_w);
		run->aps->records[a].patch_x_high = std::min(rx, aps_x[a] + patch_half_w);
		run->aps->records[a].box_y_low    = run->aps->records[a].patch_y_low;
		run->aps->records[a].box_y_high   = run->aps->records[a].patch_y_high;
		run->aps->records[a].box_x_low    = run->aps->records[a].patch_x_low;
		run->aps->records[a].box_x_high   = run->aps->records[a].patch_x_high;
		run->aps->records[a].structure    = 1.0;
	}

	run->shifts = (mpp_shifts_t *) std::calloc(1, sizeof(mpp_shifts_t));
	run->shifts->num_frames = 1;
	run->shifts->num_aps    = 2;
	run->shifts->shifts = (double *) std::calloc(2 * 1 * 2, sizeof(double));
	run->shifts->success = (uint8_t *) std::calloc(1 * 2, sizeof(uint8_t));
	run->shifts->success[0] = 1;
	run->shifts->success[1] = 1;
	/* frame 0, AP 0: dy, dx; frame 0, AP 1: dy, dx */
	run->shifts->shifts[0] = ap1_dy;
	run->shifts->shifts[1] = ap1_dx;
	run->shifts->shifts[2] = ap2_dy;
	run->shifts->shifts[3] = ap2_dx;
	return run;
}

}  // namespace

Test(mpp_pixmap, identity) {
	/* Zero global shifts, zero per-AP shifts, scale=1, intersection (0..ry,0..rx).
	 * xmap[j,i] should equal i; ymap[j,i] should equal j. */
	const int rx = 20, ry = 16;
	mpp_run_t *run = make_fixture(rx, ry, /*num_frames=*/1,
	                              /*ay=*/8, /*ax=*/10,
	                              /*patch=*/0, ry, 0, rx,
	                              /*global=*/nullptr, /*ap_shifts=*/nullptr);
	imgmap_t pm{};
	cr_assert_eq(mpp_imgmap_alloc(&pm, rx, ry), MPP_OK);
	cr_assert_eq(mpp_pixmap_build(run, 0, 1.0, &pm), MPP_OK);
	for (int j = 0; j < ry; ++j) {
		for (int i = 0; i < rx; ++i) {
			cr_assert_float_eq(pm.xmap[j * rx + i], (float) i, 1e-5f,
			                   "xmap(%d,%d) expected %d got %f", j, i, i,
			                   pm.xmap[j * rx + i]);
			cr_assert_float_eq(pm.ymap[j * rx + i], (float) j, 1e-5f,
			                   "ymap(%d,%d) expected %d got %f", j, i, j,
			                   pm.ymap[j * rx + i]);
		}
	}
	mpp_imgmap_free(&pm);
	mpp_run_free(run);
}

Test(mpp_pixmap, global_only) {
	/* Non-zero integer global shift (dy=3, dx=-2). xmap = i + dx, ymap = j + dy. */
	const int rx = 20, ry = 16;
	const int g[2] = { 3, -2 };
	mpp_run_t *run = make_fixture(rx, ry, 1, 8, 10, 0, ry, 0, rx,
	                              g, nullptr);
	imgmap_t pm{};
	cr_assert_eq(mpp_imgmap_alloc(&pm, rx, ry), MPP_OK);
	cr_assert_eq(mpp_pixmap_build(run, 0, 1.0, &pm), MPP_OK);
	for (int j = 0; j < ry; ++j) {
		for (int i = 0; i < rx; ++i) {
			cr_assert_float_eq(pm.xmap[j * rx + i], (float) (i - 2), 1e-5f);
			cr_assert_float_eq(pm.ymap[j * rx + i], (float) (j + 3), 1e-5f);
		}
	}
	mpp_imgmap_free(&pm);
	mpp_run_free(run);
}

Test(mpp_pixmap, per_ap_at_centre) {
	/* One AP at (8, 10), patch covers entire frame (extend-both on both axes),
	 * AP shift = (dy=2, dx=-1). The AP centre weight = 1 (extend_low/high
	 * push both ramps to 1.0 across the patch in the existing code, since
	 * the patch == frame and centre == ramp peak). Total per-pixel weight
	 * is min(1, 1) = 1 → effective per-AP shift = the AP shift everywhere
	 * the patch covers.
	 *
	 * Verify the AP centre gets exactly the AP shift. */
	const int rx = 20, ry = 16;
	const double ap[2] = { 2.0, -1.0 };
	mpp_run_t *run = make_fixture(rx, ry, 1, /*ay=*/8, /*ax=*/10,
	                              /*patch=*/0, ry, 0, rx,
	                              nullptr, ap);
	imgmap_t pm{};
	cr_assert_eq(mpp_imgmap_alloc(&pm, rx, ry), MPP_OK);
	cr_assert_eq(mpp_pixmap_build(run, 0, 1.0, &pm), MPP_OK);
	const int cy = 8, cx = 10;
	cr_assert_float_eq(pm.xmap[cy * rx + cx], (float) (cx + ap[1]), 1e-5f);
	cr_assert_float_eq(pm.ymap[cy * rx + cx], (float) (cy + ap[0]), 1e-5f);
	mpp_imgmap_free(&pm);
	mpp_run_free(run);
}

Test(mpp_pixmap, per_ap_between) {
	/* Two APs on the same row at x=10 and x=30, shifts (0, +2) and (0, -2).
	 * Patches overlap in the middle (half-width 10 each → x in [0, 20)
	 * and [20, 40) — back-to-back without overlap), so at x=20 only AP2
	 * is inside its patch. Test simpler invariants:
	 *  - At AP1 centre (y=10, x=10): xmap == 10 + 2 = 12.
	 *  - At AP2 centre (y=10, x=30): xmap == 30 - 2 = 28.
	 *  - At a point covered only by AP1 (y=10, x=5): xmap leans toward AP1.
	 *  - At a point covered only by AP2 (y=10, x=35): xmap leans toward AP2.
	 *  - At a point covered by neither (corner): xmap = i (global = 0). */
	const int rx = 40, ry = 20;
	mpp_run_t *run = make_two_ap_fixture(rx, ry,
	                                     /*ap1=*/10, 10, /*ap2=*/10, 30,
	                                     /*patch_half_w=*/10, /*patch_half_h=*/10,
	                                     /*ap1: dy, dx=*/0.0, 2.0,
	                                     /*ap2: dy, dx=*/0.0, -2.0);
	imgmap_t pm{};
	cr_assert_eq(mpp_imgmap_alloc(&pm, rx, ry), MPP_OK);
	cr_assert_eq(mpp_pixmap_build(run, 0, 1.0, &pm), MPP_OK);

	cr_assert_float_eq(pm.xmap[10 * rx + 10], 12.0f, 1e-4f, "AP1 centre");
	cr_assert_float_eq(pm.xmap[10 * rx + 30], 28.0f, 1e-4f, "AP2 centre");
	/* Point covered only by AP1 (x=5 is inside AP1's patch [0..20), outside AP2's [20..40)). */
	cr_assert_gt(pm.xmap[10 * rx +  5], 5.0f,  "AP1-only should add positive dx");
	/* Point covered only by AP2 (x=35 is outside AP1's, inside AP2's). */
	cr_assert_lt(pm.xmap[10 * rx + 35], 35.0f, "AP2-only should add negative dx");
	mpp_imgmap_free(&pm);
	mpp_run_free(run);
}

Test(mpp_pixmap, outside_aps_fallback) {
	/* One small AP in the centre; sample at the four corners which are
	 * outside its patch. Result should be identity (global = 0). No
	 * divide-by-zero or NaN anywhere. */
	const int rx = 40, ry = 30;
	mpp_run_t *run = make_fixture(rx, ry, 1,
	                              /*ay=*/15, /*ax=*/20,
	                              /*patch=*/13, 17, 18, 22,
	                              nullptr,
	                              (const double[]) { 0.5, -0.5 });
	imgmap_t pm{};
	cr_assert_eq(mpp_imgmap_alloc(&pm, rx, ry), MPP_OK);
	cr_assert_eq(mpp_pixmap_build(run, 0, 1.0, &pm), MPP_OK);
	/* No NaN anywhere. */
	for (int j = 0; j < ry; ++j)
		for (int i = 0; i < rx; ++i) {
			cr_assert(!std::isnan(pm.xmap[j * rx + i]));
			cr_assert(!std::isnan(pm.ymap[j * rx + i]));
		}
	/* Corners get identity. */
	cr_assert_float_eq(pm.xmap[0 * rx + 0],            0.0f, 1e-5f);
	cr_assert_float_eq(pm.ymap[0 * rx + 0],            0.0f, 1e-5f);
	cr_assert_float_eq(pm.xmap[(ry - 1) * rx + (rx - 1)], (float) (rx - 1), 1e-5f);
	cr_assert_float_eq(pm.ymap[(ry - 1) * rx + (rx - 1)], (float) (ry - 1), 1e-5f);
	mpp_imgmap_free(&pm);
	mpp_run_free(run);
}

Test(mpp_pixmap, drizzle_scale) {
	/* Zero shifts, scale=2 → xmap[j,i] = 2*i; ymap[j,i] = 2*j. */
	const int rx = 10, ry = 8;
	mpp_run_t *run = make_fixture(rx, ry, 1, 4, 5, 0, ry, 0, rx,
	                              nullptr, nullptr);
	imgmap_t pm{};
	cr_assert_eq(mpp_imgmap_alloc(&pm, rx, ry), MPP_OK);
	cr_assert_eq(mpp_pixmap_build(run, 0, 2.0, &pm), MPP_OK);
	for (int j = 0; j < ry; ++j) {
		for (int i = 0; i < rx; ++i) {
			cr_assert_float_eq(pm.xmap[j * rx + i], 2.0f * (float) i, 1e-5f);
			cr_assert_float_eq(pm.ymap[j * rx + i], 2.0f * (float) j, 1e-5f);
		}
	}
	mpp_imgmap_free(&pm);
	mpp_run_free(run);
}

Test(mpp_pixmap, invalid_inputs) {
	imgmap_t pm{};
	cr_assert_eq(mpp_imgmap_alloc(&pm, 10, 10), MPP_OK);

	cr_assert_eq(mpp_pixmap_build(nullptr, 0, 1.0, &pm), MPP_EINVAL);
	cr_assert_eq(mpp_pixmap_build(nullptr, 0, 1.0, nullptr), MPP_EINVAL);

	mpp_run_t *run = make_fixture(10, 10, 1, 4, 5, 0, 10, 0, 10, nullptr, nullptr);
	cr_assert_eq(mpp_pixmap_build(run, -1, 1.0, &pm), MPP_EINVAL);
	cr_assert_eq(mpp_pixmap_build(run,  1, 1.0, &pm), MPP_EINVAL);   /* out of range */
	cr_assert_eq(mpp_pixmap_build(run,  0, 0.0, &pm), MPP_EINVAL);
	cr_assert_eq(mpp_pixmap_build(run,  0, -1.0, &pm), MPP_EINVAL);
	mpp_run_free(run);

	/* Mismatched dimensions: alloc 10×10 but frame is 20×20. */
	mpp_run_t *run20 = make_fixture(20, 20, 1, 10, 10, 0, 20, 0, 20,
	                                nullptr, nullptr);
	cr_assert_eq(mpp_pixmap_build(run20, 0, 1.0, &pm), MPP_EINVAL);
	mpp_run_free(run20);

	mpp_imgmap_free(&pm);
}

Test(mpp_imgmap, alloc_free_zero) {
	imgmap_t pm{};
	cr_assert_eq(mpp_imgmap_alloc(&pm, 0, 10), MPP_EINVAL);
	cr_assert_eq(mpp_imgmap_alloc(&pm, 10, 0), MPP_EINVAL);
	cr_assert_eq(mpp_imgmap_alloc(nullptr, 10, 10), MPP_EINVAL);
	cr_assert_eq(mpp_imgmap_alloc(&pm, 8, 4), MPP_OK);
	cr_assert_eq(pm.rx, 8);
	cr_assert_eq(pm.ry, 4);
	cr_assert_not_null(pm.xmap);
	cr_assert_not_null(pm.ymap);
	mpp_imgmap_free(&pm);
	cr_assert_null(pm.xmap);
	cr_assert_null(pm.ymap);
	cr_assert_eq(pm.rx, 0);
	cr_assert_eq(pm.ry, 0);
	/* Safe to free twice. */
	mpp_imgmap_free(&pm);
	mpp_imgmap_free(nullptr);
}
