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
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include "registration/mpp_align_priv.hpp"
#include "registration/mpp_ap_priv.hpp"
#include "registration/mpp_drizzle_priv.hpp"
#include "registration/mpp_rank_priv.hpp"
#include "registration/mpp_shift_priv.hpp"
#include "registration/mpp_stack_priv.hpp"

extern "C" {
#include "core/siril.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_drizzle.h"
#include "registration/mpp_shift.h"
}

/* Slice 5b.2 pulls Stage-C orchestration into mpp_drizzle.cpp, which
 * transitively references `cominfo com` and `fits *gfit` (gui_iface,
 * threading prefs, etc.). Tests linking against the new TU need
 * definitions for both — same pattern as imoper_test.c /
 * stacking_blocks_test.c. The `extern` declarations live in
 * core/siril.h (C++ linkage). */
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

/* ===========================================================================
 * STScI bicubic-oracle sanity tests (Phase 5b plan, Group 2).
 *
 * Two real-SER tests that drive mpp::stack_apply_stsci on the same
 * `test-big.ser`-derived fixture as the existing Phase 5a end-to-end
 * oracle, then compare to the PSS oracle outputs. The comparison bars
 * are structural / 'rule out junk' rather than algorithmic — STScI ≠
 * bicubic, so an exact match isn't expected.
 *
 *   near-passthrough: drizzle_factor=1, pixfrac=1 vs PSS drizzle=Off
 *   real-2x:          drizzle_factor=2, pixfrac=0.6 vs PSS drizzle=2
 *
 * Skip cleanly if MPP_PSS_TEST_DATA_DIR isn't set or the fixture /
 * oracle binaries are missing.
 * ========================================================================= */

namespace {
namespace fs = std::filesystem;

mpp_config_t stsci_default_cfg() {
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	return cfg;
}

cv::Mat blurred_for_align(const cv::Mat &mono, const mpp_config_t &cfg) {
	return mpp::blur_mono_for_align(mono, cfg);
}

std::vector<int> read_ints_csv_local(const fs::path &p) {
	std::vector<int> v;
	std::ifstream fin(p);
	int x;
	while (fin >> x) v.push_back(x);
	return v;
}

/* Built-up run + the parallel raw-frame / brightness / included vectors
 * that mpp::stack_apply_stsci's signature needs. */
struct BuiltRun {
	mpp_run_t                *run = nullptr;
	std::vector<cv::Mat>     frames_raw;
	std::vector<int>         included;
	std::vector<double>      frame_brightness;
	~BuiltRun() { mpp_run_free(run); }
};

/* Run Phases 1–4 on a directory of mono PNG frames and assemble a run.
 * Mirrors the prep in mpp_stack_test.cpp::end_to_end_real_ser_oracle but
 * also computes per-AP per-frame shifts (Stage B) since mpp_pixmap_build
 * needs run->shifts to be populated. */
mpp_status_t build_real_run(const fs::path &frames_dir,
                            const mpp_config_t &cfg_in,
                            BuiltRun &out) {
	std::vector<fs::path> paths;
	for (const auto &e : fs::directory_iterator(frames_dir)) {
		const std::string n = e.path().filename().string();
		if (n.rfind("frame_", 0) == 0 && e.path().extension() == ".png")
			paths.push_back(e.path());
	}
	if (paths.empty()) return MPP_ENODATA;
	std::sort(paths.begin(), paths.end());

	mpp_config_t cfg = cfg_in;

	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q_rank, frame_brightness;
	for (const auto &p : paths) {
		cv::Mat m = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		if (m.empty()) return MPP_EIO;
		frames_raw.push_back(m);
		frames_blurred.push_back(blurred_for_align(m, cfg));
		q_rank.push_back(mpp::rank_score_normalized(m, cfg));
		frame_brightness.push_back(mpp::rank_average_brightness(m, cfg));
	}

	const auto align = mpp::align_global_from_frames(frames_blurred, q_rank, cfg);
	const auto avg   = mpp::align_average_frame(frames_raw, q_rank, align.shifts, cfg);
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	if (!aps || aps->count <= 0) { mpp_ap_free(aps); return MPP_ENODATA; }

	const auto ref_boxes = mpp::shift_prepare_ref_boxes(avg.mean_frame, *aps);
	const auto offsets   = mpp::shift_frame_offsets(align.shifts, avg.intersection);
	mpp_shifts_t *shifts = mpp::shift_compute_all(frames_blurred, *aps,
	                                              ref_boxes, offsets, cfg);
	if (!shifts) { mpp_ap_free(aps); return MPP_ENOMEM; }

	mpp_run_t *run = mpp_run_alloc();
	if (!run) {
		mpp_ap_free(aps);
		mpp_shift_free(shifts);
		return MPP_ENOMEM;
	}
	const int N = (int) frames_raw.size();
	run->cfg = (mpp_config_t *) std::malloc(sizeof(mpp_config_t));
	*run->cfg = cfg;
	run->num_frames = N;
	run->frame_rows = frames_raw[0].rows;
	run->frame_cols = frames_raw[0].cols;
	run->num_layers = frames_raw[0].channels();
	run->bitdepth   = cfg.bitdepth;
	run->aps        = aps;
	run->shifts     = shifts;
	run->global_shifts = (int *) std::calloc((size_t) 2 * N, sizeof(int));
	run->frame_brightness = (double *) std::calloc((size_t) N, sizeof(double));
	run->included   = (int *) std::calloc((size_t) N, sizeof(int));
	for (int i = 0; i < N; ++i) {
		run->global_shifts[2 * i + 0] = align.shifts[i][0];
		run->global_shifts[2 * i + 1] = align.shifts[i][1];
		run->frame_brightness[i] = frame_brightness[i];
		run->included[i] = 1;
	}
	for (int k = 0; k < 4; ++k) run->intersection[k] = avg.intersection[k];

	out.run = run;
	out.frames_raw = std::move(frames_raw);
	out.included.assign(run->included, run->included + N);
	out.frame_brightness = std::move(frame_brightness);
	return MPP_OK;
}

/* PSNR in dB between two uint16 images of the same size and channel
 * count. Uses cv::PSNR (computes via MSE; returns +inf for zero MSE). */
double psnr_u16(const cv::Mat &a, const cv::Mat &b) {
	return cv::PSNR(a, b, /*R=*/65535.0);
}

/* SSIM between two uint16 single-channel images. Standard Wang-2004
 * formulation with 11×11 Gaussian (σ=1.5) and dynamic range L=65535. */
double ssim_u16(const cv::Mat &a, const cv::Mat &b) {
	const double K1 = 0.01, K2 = 0.03, L = 65535.0;
	const double C1 = (K1 * L) * (K1 * L);
	const double C2 = (K2 * L) * (K2 * L);
	cv::Mat I1, I2;
	a.convertTo(I1, CV_32F);
	b.convertTo(I2, CV_32F);
	cv::Mat I1_2 = I1.mul(I1), I2_2 = I2.mul(I2), I1_I2 = I1.mul(I2);
	cv::Mat mu1, mu2;
	cv::GaussianBlur(I1, mu1, cv::Size(11, 11), 1.5);
	cv::GaussianBlur(I2, mu2, cv::Size(11, 11), 1.5);
	cv::Mat mu1_2 = mu1.mul(mu1), mu2_2 = mu2.mul(mu2), mu1_mu2 = mu1.mul(mu2);
	cv::Mat sigma1_2, sigma2_2, sigma12;
	cv::GaussianBlur(I1_2, sigma1_2, cv::Size(11, 11), 1.5);
	sigma1_2 -= mu1_2;
	cv::GaussianBlur(I2_2, sigma2_2, cv::Size(11, 11), 1.5);
	sigma2_2 -= mu2_2;
	cv::GaussianBlur(I1_I2, sigma12, cv::Size(11, 11), 1.5);
	sigma12 -= mu1_mu2;
	cv::Mat t1 = 2 * mu1_mu2 + C1;
	cv::Mat t2 = 2 * sigma12 + C2;
	cv::Mat numerator = t1.mul(t2);
	cv::Mat d1 = mu1_2 + mu2_2 + C1;
	cv::Mat d2 = sigma1_2 + sigma2_2 + C2;
	cv::Mat denominator = d1.mul(d2);
	cv::Mat ssim_map;
	cv::divide(numerator, denominator, ssim_map);
	return cv::mean(ssim_map)[0];
}

/* Helper: read PSS oracle stacked_u16.bin + dims.csv into a cv::Mat. */
bool read_oracle_u16(const fs::path &bin, const fs::path &dims,
                     cv::Mat &out) {
	auto v = read_ints_csv_local(dims);
	if (v.size() < 2) return false;
	const int H = v[0], W = v[1];
	out = cv::Mat(H, W, CV_16U);
	std::ifstream fin(bin, std::ios::binary);
	fin.read(reinterpret_cast<char *>(out.data),
	         (std::streamsize) H * W * sizeof(uint16_t));
	return fin.good();
}

}  // namespace

/* drizzle=1, pixfrac=1 — STScI should be very close to a no-drizzle
 * stack. Compared to PSS's bicubic stack at drizzle=Off
 * (oracle_out_real/stacked_u16.bin). Bar values are conservative
 * starting points; tighten or relax once we see real numbers. */
Test(mpp_stsci_oracle, near_passthrough) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	const fs::path data(data_root);
	const fs::path frames_dir = data / "test_data" / "real_jupiter";
	const fs::path oracle_bin = data / "oracle_out_real" / "stacked_u16.bin";
	const fs::path oracle_dim = data / "oracle_out_real" / "stacked_u16_dims.csv";
	if (!fs::exists(frames_dir) || !fs::exists(oracle_bin)) {
		cr_skip_test("real-SER fixture / oracle missing");
	}

	mpp_config_t cfg = stsci_default_cfg();
	cfg.bitdepth = 8;
	cfg.drizzle_mode    = MPP_DRIZZLE_STSCI;
	cfg.drizzle_factor  = 1;
	cfg.drizzle_pixfrac = 1.0;
	cfg.drizzle_kernel  = MPP_KERNEL_SQUARE;

	BuiltRun b;
	cr_assert_eq(build_real_run(frames_dir, cfg, b), MPP_OK);

	fits out{};
	const mpp_status_t rc = mpp::stack_apply_stsci(b.frames_raw, b.included,
	                                               b.frame_brightness,
	                                               b.run, &cfg, &out);
	cr_assert_eq(rc, MPP_OK);

	cv::Mat oracle;
	cr_assert(read_oracle_u16(oracle_bin, oracle_dim, oracle),
	          "failed to read PSS oracle");

	cv::Mat got(out.ry, out.rx, CV_16U, out.data);

	/* PSS's bicubic stack applies a 1-pixel border trim that the STScI
	 * path doesn't; STScI output may be a couple of pixels larger. Take
	 * the centre-aligned overlap of the two for comparison, then crop a
	 * further 10-px safety margin off each side — edge effects from
	 * dobox's partial-overlap weighting can be visibly different from
	 * the bicubic path's edge handling without indicating an actual
	 * regression. */
	const int common_h = std::min(got.rows, oracle.rows);
	const int common_w = std::min(got.cols, oracle.cols);
	const int crop = 10;
	cr_assert_gt(common_h, 2 * crop);
	cr_assert_gt(common_w, 2 * crop);
	cv::Rect roi_got(((got.cols    - common_w) / 2) + crop,
	                 ((got.rows    - common_h) / 2) + crop,
	                 common_w - 2 * crop, common_h - 2 * crop);
	cv::Rect roi_or (((oracle.cols - common_w) / 2) + crop,
	                 ((oracle.rows - common_h) / 2) + crop,
	                 common_w - 2 * crop, common_h - 2 * crop);
	cv::Mat got_c = got(roi_got), or_c = oracle(roi_or);

	const double psnr = psnr_u16(got_c, or_c);
	const double ssim = ssim_u16(got_c, or_c);
	cv::Mat absdiff;
	cv::absdiff(got_c, or_c, absdiff);
	cv::Mat absdiff_f;
	absdiff.convertTo(absdiff_f, CV_32F);
	const double mean_abs = cv::mean(absdiff_f)[0];
	double worst_d = 0;
	cv::minMaxLoc(absdiff_f, nullptr, &worst_d);
	const double mean_got = cv::mean(got_c)[0];
	const double mean_or  = cv::mean(or_c)[0];

	siril_log_color_message("[near_passthrough] dims got=%dx%d oracle=%dx%d "
	                        "means got=%.0f oracle=%.0f "
	                        "PSNR=%.2f dB SSIM=%.4f mean|Δ|=%.1f worst=%.0f\n",
	                        "green", got.cols, got.rows, oracle.cols, oracle.rows,
	                        mean_got, mean_or, psnr, ssim, mean_abs, worst_d);

	if (const char *dump = std::getenv("MPP_DUMP_RESULT_DIR")) {
		const fs::path out_dir(dump);
		cv::imwrite((out_dir / "stsci_passthrough_siril.png").string(), got);
		cv::imwrite((out_dir / "stsci_passthrough_oracle.png").string(), oracle);
	}

	cr_assert_geq(psnr, 35.0,
	              "PSNR %.2f dB below bar 35 — STScI passthrough drifting "
	              "from bicubic-Off oracle", psnr);
	cr_assert_geq(ssim, 0.95,
	              "SSIM %.4f below bar 0.95 — structural mismatch", ssim);
	cr_assert_leq(mean_abs, 200.0,
	              "mean|Δ| %.1f above bar 200 / 65535", mean_abs);
	cr_assert_leq(worst_d, 5000.0,
	              "worst|Δ| %.0f above bar 5000 / 65535", worst_d);

	clearfits(&out);
}

/* drizzle=2, pixfrac=1.0 — STScI 2× drizzle vs the PSS bicubic-2x mono
 * oracle.
 *
 * Earlier this test used `oracle_out_realser_drizzle2`, which was
 * generated from the Bayer SER directly and is 3-channel RGB framed
 * differently from the mono drizzle=Off oracle. That made the
 * comparison apples-to-oranges (offset planets, different sizes,
 * different colours). The new `oracle_out_real_drizzle2` is generated
 * from the same mono PNG fixture as oracle_out_real (drizzle=Off) but
 * with PSS's drizzle_factor=2 — so the inputs match and any remaining
 * difference is purely the algorithmic gap between STScI's boxer-
 * overlap accumulation and PSS's bicubic upsample-then-resample.
 *
 * Empirically that gap is small: max abs diff 1142 / 65535 (1.7%),
 * mean abs diff 186 / 65535 (0.28%). The diff is concentrated at the
 * planet's limb (high-gradient region — sub-pixel registration
 * differences show up there). */
Test(mpp_stsci_oracle, real_2x) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	const fs::path data(data_root);
	const fs::path frames_dir = data / "test_data" / "real_jupiter";
	/* oracle_out_real_drizzle2 — generated by run_pss.py on the same mono
	 * PNG dir as oracle_out_real (drizzle=Off) but with --drizzle=2x.
	 * oracle_out_realser_drizzle2 also exists but was generated from the
	 * Bayer SER directly (3-channel RGB output, different framing),
	 * which made the comparison apples-to-oranges. */
	const fs::path oracle_bin = data / "oracle_out_real_drizzle2" / "stacked_u16.bin";
	const fs::path oracle_dim = data / "oracle_out_real_drizzle2" / "stacked_u16_dims.csv";
	if (!fs::exists(frames_dir) || !fs::exists(oracle_bin)) {
		cr_skip_test("real-SER fixture / drizzle-2 oracle missing");
	}

	mpp_config_t cfg = stsci_default_cfg();
	cfg.bitdepth = 8;
	cfg.drizzle_mode    = MPP_DRIZZLE_STSCI;
	cfg.drizzle_factor  = 2;
	cfg.drizzle_pixfrac = 1.0;   /* full-pixel-fraction at 2x: cleanest comparison */
	cfg.drizzle_kernel  = MPP_KERNEL_SQUARE;

	BuiltRun b;
	cr_assert_eq(build_real_run(frames_dir, cfg, b), MPP_OK);

	fits out{};
	const mpp_status_t rc = mpp::stack_apply_stsci(b.frames_raw, b.included,
	                                               b.frame_brightness,
	                                               b.run, &cfg, &out);
	cr_assert_eq(rc, MPP_OK);

	cv::Mat oracle;
	cr_assert(read_oracle_u16(oracle_bin, oracle_dim, oracle),
	          "failed to read PSS oracle at drizzle=2");

	cv::Mat got(out.ry, out.rx, CV_16U, out.data);

	/* Centre-aligned overlap crop. The PSS 2x output's intersection
	 * dims may be off by 1 px vs ours due to PSS's dynamic border trim;
	 * crop a few px off each side beyond that to skip the edge ring
	 * (the gradient-heavy limb where most of the diff concentrates). */
	const int common_h = std::min(got.rows, oracle.rows);
	const int common_w = std::min(got.cols, oracle.cols);
	const int crop = 10;
	cr_assert_gt(common_h, 2 * crop);
	cr_assert_gt(common_w, 2 * crop);
	cv::Rect roi_got(((got.cols    - common_w) / 2) + crop,
	                 ((got.rows    - common_h) / 2) + crop,
	                 common_w - 2 * crop, common_h - 2 * crop);
	cv::Rect roi_or (((oracle.cols - common_w) / 2) + crop,
	                 ((oracle.rows - common_h) / 2) + crop,
	                 common_w - 2 * crop, common_h - 2 * crop);
	cv::Mat got_c = got(roi_got), or_c = oracle(roi_or);

	const double psnr = psnr_u16(got_c, or_c);
	const double ssim = ssim_u16(got_c, or_c);
	cv::Mat absdiff;
	cv::absdiff(got_c, or_c, absdiff);
	cv::Mat absdiff_f;
	absdiff.convertTo(absdiff_f, CV_32F);
	const double mean_abs = cv::mean(absdiff_f)[0];
	double worst_d = 0;
	cv::minMaxLoc(absdiff_f, nullptr, &worst_d);
	const double mean_got = cv::mean(got_c)[0];
	const double mean_or  = cv::mean(or_c)[0];

	siril_log_color_message("[real_2x] dims got=%dx%d oracle=%dx%d "
	                        "means got=%.0f oracle=%.0f "
	                        "PSNR=%.2f dB SSIM=%.4f mean|Δ|=%.1f worst=%.0f\n",
	                        "green",
	                        got.cols, got.rows, oracle.cols, oracle.rows,
	                        mean_got, mean_or, psnr, ssim, mean_abs, worst_d);

	if (const char *dump = std::getenv("MPP_DUMP_RESULT_DIR")) {
		const fs::path out_dir(dump);
		cv::imwrite((out_dir / "stsci_real_2x_siril.png").string(), got);
		cv::imwrite((out_dir / "stsci_real_2x_oracle.png").string(), oracle);
	}

	/* Bars derived from observed values; loosened ~3× for safety against
	 * minor algorithmic refinements. Observed: PSNR ~45 dB, SSIM ~0.99,
	 * mean|Δ| ~190 / 65535, worst|Δ| ~1200 / 65535. */
	cr_assert_geq(psnr, 30.0,
	              "PSNR %.2f dB below bar 30 — STScI 2x drifting from PSS 2x oracle",
	              psnr);
	cr_assert_geq(ssim, 0.95,
	              "SSIM %.4f below bar 0.95 — structural mismatch", ssim);
	cr_assert_leq(mean_abs, 600.0,
	              "mean|Δ| %.1f above bar 600 / 65535", mean_abs);
	cr_assert_leq(worst_d, 5000.0,
	              "worst|Δ| %.0f above bar 5000 / 65535", worst_d);

	clearfits(&out);
}

/* ===========================================================================
 * Bayer drizzle synthetic-mosaic smoke test (slice 5b.3 gate).
 *
 * Generates a synthetic RGGB Bayer mosaic from a known-constant RGB
 * ground truth and runs mpp::stack_apply_bayer at drizzle=1 pixfrac=1
 * with a minimal one-frame, zero-shift run. Verifies that dobox's CFA
 * routing distributes the input samples to the correct output channels
 * with our non-affine pixmap.
 *
 * This is the "does dobox's Bayer branch behave with the MPP pixmap"
 * gate test the plan calls for. The slanted-edge MTF test (5b.6) lands
 * later.
 * ========================================================================= */

namespace {
mpp_run_t *make_minimal_run(int rx, int ry, int num_layers, int bitdepth,
                            double brightness = 100.0) {
	mpp_run_t *run = mpp_run_alloc();
	if (!run) return nullptr;
	run->num_frames = 1;
	run->frame_rows = ry;
	run->frame_cols = rx;
	run->num_layers = num_layers;
	run->bitdepth   = bitdepth;
	run->intersection[0] = 0;
	run->intersection[1] = ry;
	run->intersection[2] = 0;
	run->intersection[3] = rx;
	run->cfg = (mpp_config_t *) std::calloc(1, sizeof(mpp_config_t));
	mpp_config_defaults(run->cfg);
	run->cfg->bitdepth = bitdepth;
	run->global_shifts = (int *) std::calloc(2, sizeof(int));
	run->frame_brightness = (double *) std::calloc(1, sizeof(double));
	run->frame_brightness[0] = brightness;
	run->included = (int *) std::calloc(1, sizeof(int));
	run->included[0] = 1;
	run->aps = (mpp_aps_t *) std::calloc(1, sizeof(mpp_aps_t));
	run->aps->count = 1;
	run->aps->records = (mpp_ap_record_t *) std::calloc(1, sizeof(mpp_ap_record_t));
	run->aps->records[0].y = ry / 2;
	run->aps->records[0].x = rx / 2;
	run->aps->records[0].patch_y_low  = 0;
	run->aps->records[0].patch_y_high = ry;
	run->aps->records[0].patch_x_low  = 0;
	run->aps->records[0].patch_x_high = rx;
	run->aps->records[0].box_y_low    = 0;
	run->aps->records[0].box_y_high   = ry;
	run->aps->records[0].box_x_low    = 0;
	run->aps->records[0].box_x_high   = rx;
	run->aps->records[0].structure    = 1.0;
	run->shifts = (mpp_shifts_t *) std::calloc(1, sizeof(mpp_shifts_t));
	run->shifts->num_frames = 1;
	run->shifts->num_aps = 1;
	run->shifts->shifts = (double *) std::calloc(2, sizeof(double));   /* zero */
	run->shifts->success = (uint8_t *) std::calloc(1, sizeof(uint8_t));
	run->shifts->success[0] = 1;
	return run;
}
}  // namespace

Test(mpp_bayer_drizzle, synthetic_mosaic_rggb) {
	/* Ground truth: 64×64 RGB constant — R=64, G=128, B=192 (in 0..255).
	 * Pick distinct, non-saturating values so cross-channel leakage in
	 * the output (if any) shows up clearly. */
	const int W = 64, H = 64;
	const uint8_t R_VAL = 64, G_VAL = 128, B_VAL = 192;

	/* Mosaic to RGGB single-channel:
	 *   row even, col even → R
	 *   row even, col odd  → G
	 *   row odd,  col even → G
	 *   row odd,  col odd  → B           */
	cv::Mat bayer(H, W, CV_8UC1);
	for (int y = 0; y < H; ++y) {
		uint8_t *row = bayer.ptr<uint8_t>(y);
		for (int x = 0; x < W; ++x) {
			const int kind = ((y & 1) << 1) | (x & 1);
			row[x] = (kind == 0) ? R_VAL
			       : (kind == 3) ? B_VAL
			       : G_VAL;
		}
	}

	mpp_config_t cfg = stsci_default_cfg();
	cfg.bitdepth = 8;
	cfg.drizzle_mode    = MPP_DRIZZLE_BAYER;
	cfg.drizzle_factor  = 1;
	cfg.drizzle_pixfrac = 1.0;
	cfg.drizzle_kernel  = MPP_KERNEL_SQUARE;

	mpp_run_t *run = make_minimal_run(W, H, /*num_layers=*/3, /*bitdepth=*/8,
	                                  /*brightness=*/(R_VAL + 2*G_VAL + B_VAL) / 4.0);
	cr_assert_not_null(run);

	const std::vector<cv::Mat> frames = { bayer };
	const std::vector<int> included = { 1 };
	const std::vector<double> brightness = { run->frame_brightness[0] };

	/* RGGB pattern (matches get_compiled_pattern output for BAYER_RGGB —
	 * demosaicing.c:333). RLAYER=0, GLAYER=1, BLAYER=2. */
	const unsigned char cfa[4] = { 0, 1, 1, 2 };

	fits out{};
	const mpp_status_t rc = mpp::stack_apply_bayer(frames, included, brightness,
	                                               run, &cfg, cfa, /*cfadim=*/2, &out);
	cr_assert_eq(rc, MPP_OK);
	cr_assert_eq(out.naxes[2], 3);
	cr_assert_eq(out.rx, W);
	cr_assert_eq(out.ry, H);

	/* Wrap each channel as a cv::Mat for stats. Output layout is planar
	 * (pdata[0..2] point into a single contiguous buffer). */
	const size_t plane = (size_t) W * (size_t) H;
	cv::Mat ch_r(H, W, CV_16U, out.data);
	cv::Mat ch_g(H, W, CV_16U, out.data + plane);
	cv::Mat ch_b(H, W, CV_16U, out.data + 2 * plane);

	/* Each channel receives input only at the corresponding CFA positions
	 * (R: 25 %, G: 50 %, B: 25 % of pixels). Where it doesn't receive,
	 * the cell stays at calloc'd zero — there's no interpolation in our
	 * drizzle path. So mean over the whole plane is:
	 *   mean_R = R_VAL × 0.25 × scale_8_to_16
	 *   mean_G = G_VAL × 0.50 × scale_8_to_16
	 *   mean_B = B_VAL × 0.25 × scale_8_to_16   */
	const double scale_8_to_16 = 256.0;
	const double exp_R = (double) R_VAL * 0.25 * scale_8_to_16;
	const double exp_G = (double) G_VAL * 0.50 * scale_8_to_16;
	const double exp_B = (double) B_VAL * 0.25 * scale_8_to_16;
	const double mean_R = cv::mean(ch_r)[0];
	const double mean_G = cv::mean(ch_g)[0];
	const double mean_B = cv::mean(ch_b)[0];

	siril_log_color_message(
	    "[bayer_rggb] means R=%.0f (exp %.0f)  G=%.0f (exp %.0f)  B=%.0f (exp %.0f)\n",
	    "green", mean_R, exp_R, mean_G, exp_G, mean_B, exp_B);

	const double tol = 1.0;   /* 1 / 65535 — should be near-bit-exact */
	cr_assert_lt(std::abs(mean_R - exp_R), tol,
	             "R channel mean %.2f drifts from expected %.2f", mean_R, exp_R);
	cr_assert_lt(std::abs(mean_G - exp_G), tol,
	             "G channel mean %.2f drifts from expected %.2f", mean_G, exp_G);
	cr_assert_lt(std::abs(mean_B - exp_B), tol,
	             "B channel mean %.2f drifts from expected %.2f", mean_B, exp_B);

	/* Cross-channel leakage: at a pixel where the Bayer position says
	 * "this is a R sample", the G and B channels should be zero (no
	 * input pixel mapped there for those channels). Check at (0, 0)
	 * which is an R position in RGGB. */
	cr_assert_eq(ch_r.at<uint16_t>(0, 0), (uint16_t)(R_VAL * scale_8_to_16),
	             "R(0,0) should be the R sample value");
	cr_assert_eq(ch_g.at<uint16_t>(0, 0), 0,
	             "G(0,0) should be zero (no G sample mapped here)");
	cr_assert_eq(ch_b.at<uint16_t>(0, 0), 0,
	             "B(0,0) should be zero (no B sample mapped here)");

	/* And at a G-position (0, 1) and B-position (1, 1): */
	cr_assert_eq(ch_r.at<uint16_t>(0, 1), 0,                                "R(0,1) zero");
	cr_assert_eq(ch_g.at<uint16_t>(0, 1), (uint16_t)(G_VAL * scale_8_to_16), "G(0,1) green");
	cr_assert_eq(ch_b.at<uint16_t>(0, 1), 0,                                "B(0,1) zero");
	cr_assert_eq(ch_r.at<uint16_t>(1, 1), 0,                                "R(1,1) zero");
	cr_assert_eq(ch_g.at<uint16_t>(1, 1), 0,                                "G(1,1) zero");
	cr_assert_eq(ch_b.at<uint16_t>(1, 1), (uint16_t)(B_VAL * scale_8_to_16), "B(1,1) blue");

	clearfits(&out);
	mpp_run_free(run);
}
