/*
 * Phase 2 — mpp_align Criterion tests.
 *
 * Covers PSS-faithful patch picker (compute_alignment_rect /
 * quality_measure_threshold_weighted) and two-phase MultiLevelCorrelation
 * via the private C++ entry points, plus the C-level mpp_align_global stub.
 */
#include <criterion/criterion.h>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "registration/mpp/mpp_align_priv.hpp"
#include "registration/mpp/mpp_rank_priv.hpp"
#include "registration/registration.h"   /* H_from_translation, for regdata seed test */

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp/mpp_align.h"
#include "registration/mpp/mpp_config.h"
}

cominfo com;
fits *gfit = nullptr;

namespace {

cv::Mat make_textured_frame(int rows = 200, int cols = 240, int seed = 1) {
	cv::Mat f(rows, cols, CV_16U, cv::Scalar(0));
	const double cy = rows / 2.0, cx = cols / 2.0;
	const double r0 = std::min(rows, cols) * 0.30;
	cv::RNG rng(seed);
	for (int y = 0; y < rows; ++y) {
		for (int x = 0; x < cols; ++x) {
			const double r = std::hypot(y - cy, x - cx);
			double v = 0.0;
			if (r < r0) {
				v = 0.55 + 0.25 * std::tanh((r0 - r) / 4.0);
				v += 0.03 * rng.gaussian(1.0);
			}
			v = std::clamp(v, 0.0, 1.0);
			f.at<uint16_t>(y, x) = static_cast<uint16_t>(v * 65535.0);
		}
	}
	rng.state = seed;
	for (int i = 0; i < 40; ++i) {
		const int sy = rng.uniform(int(cy - r0 * 0.6), int(cy + r0 * 0.6));
		const int sx = rng.uniform(int(cx - r0 * 0.6), int(cx + r0 * 0.6));
		cv::circle(f, {sx, sy}, 2, cv::Scalar(60000), -1);
	}
	return f;
}

cv::Mat blurred(const cv::Mat &mono, const mpp_config_t &cfg) {
	cv::Mat b;
	cv::GaussianBlur(mono, b, cv::Size(cfg.frames_gauss_width, cfg.frames_gauss_width), 0);
	return b;
}

cv::Mat shifted(const cv::Mat &mono, int dy, int dx) {
	const cv::Mat M = (cv::Mat_<double>(2, 3) << 1, 0, dx, 0, 1, dy);
	cv::Mat out;
	cv::warpAffine(mono, out, M, mono.size(), cv::INTER_NEAREST,
	               cv::BORDER_REPLICATE);
	return out;
}

mpp_config_t default_cfg() {
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	return cfg;
}

}  // namespace

/* ------------------------------------------------------------------------- */

Test(mpp_align, stub_global_returns_enotimpl) {
	int shifts[4] = {0};
	int patch[4] = {0};
	cr_assert_eq(mpp_align_global(nullptr, nullptr, nullptr, shifts, patch, nullptr),
	             MPP_ENOTIMPL);
}

Test(mpp_align, default_config_has_pss_values) {
	const auto cfg = default_cfg();
	cr_assert_eq(cfg.align_frames_search_width, 34);
	cr_assert_float_eq(cfg.align_frames_rectangle_scale_factor, 3.0, 1e-12);
	cr_assert_eq(cfg.align_frames_border_width, 10);
	cr_assert_eq(cfg.align_frames_rectangle_stride, 2);
	cr_assert_eq(cfg.align_frames_rectangle_black_threshold, 10240);
	cr_assert_float_eq(cfg.align_frames_rectangle_min_fraction, 0.7, 1e-12);
	cr_assert_eq(cfg.align_frames_average_frame_percent, 5);
	cr_assert_eq(cfg.align_frames_mode, MPP_ALIGN_SURFACE);
}

Test(mpp_align, patch_picker_returns_bounds_within_frame) {
	const auto cfg = default_cfg();
	const cv::Mat f = blurred(make_textured_frame(), cfg);
	const cv::Vec4i p = mpp::align_pick_patch(f, cfg);
	cr_assert_lt(p[0], p[1]);
	cr_assert_lt(p[2], p[3]);
	cr_assert_geq(p[0], cfg.align_frames_border_width + cfg.align_frames_search_width);
	cr_assert_leq(p[1], f.rows - (cfg.align_frames_border_width + cfg.align_frames_search_width));
	cr_assert_geq(p[2], cfg.align_frames_border_width + cfg.align_frames_search_width);
	cr_assert_leq(p[3], f.cols - (cfg.align_frames_border_width + cfg.align_frames_search_width));
}

/* MultiLevelCorrelation must recover a known integer shift. */
Test(mpp_align, multilevel_correlation_recovers_known_shift) {
	const auto cfg = default_cfg();
	const cv::Mat truth = blurred(make_textured_frame(), cfg);

	const cv::Vec4i patch = mpp::align_pick_patch(truth, cfg);
	cv::Mat ref_f32;
	truth.convertTo(ref_f32, CV_32F);
	const cv::Mat ref_window = ref_f32(cv::Range(patch[0], patch[1]),
	                                   cv::Range(patch[2], patch[3]));
	cv::Mat ref_window_first;
	{
		const int s = 2;
		const int new_h = (ref_window.rows + s - 1) / s;
		const int new_w = (ref_window.cols + s - 1) / s;
		ref_window_first.create(new_h, new_w, CV_32F);
		for (int y = 0; y < new_h; ++y)
			for (int x = 0; x < new_w; ++x)
				ref_window_first.at<float>(y, x) = ref_window.at<float>(y * s, x * s);
	}

	for (int dy : {-3, -1, 0, 2, 5}) {
		for (int dx : {-4, -2, 0, 1, 3}) {
			const cv::Mat moved = shifted(truth, dy, dx);
			const auto r = mpp::align_shift_one_frame(
			    ref_window, ref_window_first, moved,
			    patch[0], patch[1], patch[2], patch[3], cfg);
			cr_assert(r.success, "shift=(%d,%d) should have succeeded", dy, dx);
			/* PSS sign convention: shift_y reported = shift to apply to *moved*
			 * to get back to truth.  truth(y,x) ≈ moved(y+dy, x+dx).  Aligning
			 * moved to truth means shifting moved by (-dy, -dx).  Wait — let's
			 * verify experimentally: warpAffine with translation (dx,dy) means
			 * the output pixel at (y,x) = input pixel at (y-dy, x-dx).  To
			 * realign, we need shift = (-dy, -dx). Sub-pixel solver is on,
			 * so a parabolic-fit residual ≤0.5 px is expected; tolerance
			 * captures both that and any synth-floor noise. */
			cr_assert_float_eq(r.dy, (double)(-dy), 0.5,
			             "dy mismatch at injected (%d,%d): expected %d, got %g",
			             dy, dx, -dy, r.dy);
			cr_assert_float_eq(r.dx, (double)(-dx), 0.5,
			             "dx mismatch at injected (%d,%d): expected %d, got %g",
			             dy, dx, -dx, r.dx);
		}
	}
}

/* Multi-frame: each frame has known displacement, recovered shifts agree. */
Test(mpp_align, align_global_recovers_known_shifts_per_frame) {
	const auto cfg = default_cfg();
	const cv::Mat truth = blurred(make_textured_frame(), cfg);
	const std::vector<std::pair<int, int>> jit = {
	    {0, 0}, {-2, 1}, {3, -2}, {1, 2}, {-1, -1}, {0, 3}};

	std::vector<cv::Mat> frames;
	for (auto p : jit) frames.push_back(shifted(truth, p.first, p.second));

	/* Quality: synth makes them all "equal-ish"; bias so frame 0 wins (ref). */
	std::vector<double> q(frames.size(), 0.5);
	q[0] = 1.0;

	const auto r = mpp::align_global_from_frames(frames, q, cfg);
	cr_assert_eq(r.best_frame_idx, 0);
	for (size_t i = 0; i < jit.size(); ++i) {
		cr_assert_float_eq(r.shifts[i][0], (double)(-jit[i].first), 0.5,
		             "frame %zu dy: expected %d, got %g", i,
		             -jit[i].first, r.shifts[i][0]);
		cr_assert_float_eq(r.shifts[i][1], (double)(-jit[i].second), 0.5,
		             "frame %zu dx: expected %d, got %g", i,
		             -jit[i].second, r.shifts[i][1]);
	}
}

/* ---- Planet / center-of-gravity mode ----------------------------------- */

/* CoG of a symmetric disc is its geometric centre. */
Test(mpp_align, cog_recovers_disc_centre) {
	cv::Mat f(200, 240, CV_16U, cv::Scalar(0));
	const double cy = 99.0, cx = 119.0;   /* integer-grid centre of 200x240 */
	cv::circle(f, cv::Point((int) cx, (int) cy), 50, cv::Scalar(50000), -1);
	const cv::Vec2d cog = mpp::center_of_gravity(f);
	cr_assert_float_eq(cog[0], cy, 0.05, "cog y: expected %g, got %g", cy, cog[0]);
	cr_assert_float_eq(cog[1], cx, 0.05, "cog x: expected %g, got %g", cx, cog[1]);
}

/* The centroid must be sub-pixel — PSS rounds, the port does not. Two equal
 * bright pixels at x=10 and x=11 give an analytic centroid at x=10.5. */
Test(mpp_align, cog_is_subpixel_not_rounded) {
	cv::Mat f(20, 20, CV_16U, cv::Scalar(0));
	f.at<uint16_t>(10, 10) = 40000;
	f.at<uint16_t>(10, 11) = 40000;
	/* threshold = trunc((0+40000)/2) = 20000; weights = 20000 each.
	 * cog_x = (10+11)/2 = 10.5, cog_y = 10.0. */
	const cv::Vec2d cog = mpp::center_of_gravity(f);
	cr_assert_float_eq(cog[0], 10.0, 1e-9, "cog y: got %g", cog[0]);
	cr_assert_float_eq(cog[1], 10.5, 1e-9, "cog x: got %g (must be sub-pixel)", cog[1]);
}

/* Planet mode recovers known per-frame shifts via centroid displacement.
 * A rigidly translated image has a centroid that shifts by exactly the same
 * amount, so the internal texture is irrelevant. */
Test(mpp_align, planet_mode_recovers_known_shifts) {
	auto cfg = default_cfg();
	cfg.align_frames_mode = MPP_ALIGN_PLANET;
	const cv::Mat truth = blurred(make_textured_frame(), cfg);
	const std::vector<std::pair<int, int>> jit = {
	    {0, 0}, {-2, 1}, {3, -2}, {1, 2}, {-1, -1}, {0, 3}};

	std::vector<cv::Mat> frames;
	for (auto p : jit) frames.push_back(shifted(truth, p.first, p.second));

	std::vector<double> q(frames.size(), 0.5);
	q[0] = 1.0;   /* frame 0 is the reference */

	const auto r = mpp::align_global_from_frames(frames, q, cfg);
	cr_assert_eq(r.best_frame_idx, 0);
	cr_assert_eq(r.patch_yxyx[0], 0);   /* no patch in planet mode */
	cr_assert_eq(r.patch_yxyx[1], 0);
	for (size_t i = 0; i < jit.size(); ++i) {
		cr_assert_float_eq(r.shifts[i][0], (double)(-jit[i].first), 0.1,
		             "frame %zu dy: expected %d, got %g", i,
		             -jit[i].first, r.shifts[i][0]);
		cr_assert_float_eq(r.shifts[i][1], (double)(-jit[i].second), 0.1,
		             "frame %zu dx: expected %d, got %g", i,
		             -jit[i].second, r.shifts[i][1]);
	}
}

/* ---- Regdata gross-shift seed --------------------------------------------- */

/* A jump larger than align_frames_search_width cannot be recovered by the
 * bare MLC search, but a per-frame seed re-centres the window so only the
 * residual remains. The seed here is intentionally a few px off from the
 * true shift to prove the residual refinement runs on top of it. */
Test(mpp_align, seed_recovers_shift_beyond_search_width) {
	auto cfg = default_cfg();
	cfg.align_frames_search_width = 10;       /* shrink so a modest jump exceeds it */
	const cv::Mat truth = blurred(make_textured_frame(), cfg);
	const int big = 16;                        /* > 10 → out of bare MLC reach */
	std::vector<cv::Mat> frames = { truth, shifted(truth, big, 0) };
	std::vector<double> q = { 1.0, 0.5 };      /* best = frame 0 (reference) */

	/* Without a seed the jump is unreachable, so frame 1 is left near zero. */
	const auto bare = mpp::align_global_from_frames(frames, q, cfg);
	cr_assert_gt(std::abs(bare.shifts[1][0] - (double)(-big)), 5.0,
	             "without seed a >search_width jump should NOT be recovered (got dy=%g)",
	             bare.shifts[1][0]);

	/* Seed in the to-align convention (same as shifts), deliberately 3 px shy
	 * of the true -16 so the MLC has a non-trivial residual to solve. */
	std::vector<cv::Vec2d> seed = { cv::Vec2d(0.0, 0.0), cv::Vec2d(-(big - 3), 0.0) };
	const auto seeded = mpp::align_global_from_frames(frames, q, cfg, nullptr, nullptr, seed);
	cr_assert_float_eq(seeded.shifts[1][0], (double)(-big), 0.6,
	             "with seed the full shift should be recovered (got dy=%g)",
	             seeded.shifts[1][0]);
	cr_assert_float_eq(seeded.shifts[1][1], 0.0, 0.6,
	             "with seed dx should stay ~0 (got dx=%g)", seeded.shifts[1][1]);
}

/* An empty / wrong-sized seed must leave the cumulative path byte-for-byte
 * unchanged — the seed feature is strictly opt-in. */
Test(mpp_align, empty_seed_matches_unseeded) {
	const auto cfg = default_cfg();
	const cv::Mat truth = blurred(make_textured_frame(), cfg);
	std::vector<cv::Mat> frames = { truth, shifted(truth, 2, -3), shifted(truth, -1, 1) };
	std::vector<double> q = { 1.0, 0.5, 0.5 };

	const auto base = mpp::align_global_from_frames(frames, q, cfg);
	const std::vector<cv::Vec2d> empty_seed;            /* size 0 != N → ignored */
	const auto with_empty = mpp::align_global_from_frames(frames, q, cfg, nullptr,
	                                                      nullptr, empty_seed);
	for (size_t i = 0; i < frames.size(); ++i) {
		cr_assert_float_eq(with_empty.shifts[i][0], base.shifts[i][0], 1e-12);
		cr_assert_float_eq(with_empty.shifts[i][1], base.shifts[i][1], 1e-12);
	}
}

/* seed_from_regdata must reproduce the port's (dy, dx) to-align convention.
 * regdata is written with Siril's canonical H_from_translation; the seed must
 * invert it via translation_from_H (dx=h02, dy=-h12) and store (dy, dx). */
Test(mpp_align, seed_from_regdata_uses_port_convention) {
	sequence seq{};
	seq.number = 2;
	seq.nb_layers = 1;
	regdata *rd = (regdata *) calloc(2, sizeof(regdata));
	regdata *layers[1] = { rd };
	seq.regparam = layers;

	const double dx_ta = 8.0, dy_ta = -12.0;   /* an arbitrary to-align translation */
	rd[0].H = H_from_translation(0.0, 0.0);
	rd[1].H = H_from_translation(dx_ta, dy_ta);

	const std::vector<cv::Vec2d> seed = mpp::seed_from_regdata(&seq);
	cr_assert_eq(seed.size(), 2u);
	cr_assert_float_eq(seed[0][0], 0.0, 1e-9);
	cr_assert_float_eq(seed[0][1], 0.0, 1e-9);
	cr_assert_float_eq(seed[1][0], dy_ta, 1e-9, "seed dy: got %g", seed[1][0]);
	cr_assert_float_eq(seed[1][1], dx_ta, 1e-9, "seed dx: got %g", seed[1][1]);
	free(rd);
}

/* End-to-end, the careful one: a physical jump beyond search_width, regdata
 * written for that jump exactly as Siril's registration would, then the
 * regdata-derived seed must drive the real aligner to recover it. A wrong
 * sign would place the window ~2x the shift away and the assert would fail. */
Test(mpp_align, seed_from_regdata_drives_alignment) {
	auto cfg = default_cfg();
	cfg.align_frames_search_width = 10;
	const cv::Mat truth = blurred(make_textured_frame(), cfg);
	const int Dy = 12, Dx = -8;     /* content displacement; |Dy| > search_width */
	std::vector<cv::Mat> frames = { truth, shifted(truth, Dy, Dx) };
	std::vector<double> q = { 1.0, 0.5 };

	/* Siril stores the to-align translation, which for a content shift of
	 * (Dx, Dy) is (-Dx, -Dy). */
	sequence seq{};
	seq.number = 2;
	seq.nb_layers = 1;
	regdata *rd = (regdata *) calloc(2, sizeof(regdata));
	regdata *layers[1] = { rd };
	seq.regparam = layers;
	rd[0].H = H_from_translation(0.0, 0.0);
	rd[1].H = H_from_translation(-(double) Dx, -(double) Dy);
	const std::vector<cv::Vec2d> seed = mpp::seed_from_regdata(&seq);

	const auto bare = mpp::align_global_from_frames(frames, q, cfg);
	cr_assert_gt(std::abs(bare.shifts[1][0] - (double)(-Dy)), 4.0,
	             "unseeded should miss the >search_width jump (dy=%g)", bare.shifts[1][0]);

	const auto seeded = mpp::align_global_from_frames(frames, q, cfg, nullptr, nullptr, seed);
	cr_assert_float_eq(seeded.shifts[1][0], (double)(-Dy), 0.6,
	             "seeded dy: expected %d, got %g", -Dy, seeded.shifts[1][0]);
	cr_assert_float_eq(seeded.shifts[1][1], (double)(-Dx), 0.6,
	             "seeded dx: expected %d, got %g", -Dx, seeded.shifts[1][1]);
	free(rd);
}

/* find_best_frames must pick the highest-quality subset constrained to a
 * sliding window of size region_size. */
Test(mpp_align, find_best_frames_sliding_window) {
	/* Global top-2 (0.95 at idx 3, 0.9 at idx 1) are 2 apart, so they don't
	 * fit in a window of size 2. Sliding 2-windows have sums:
	 *   [0,1]:0.1+0.9=1.0   [1,2]:0.9+0.4=1.3
	 *   [2,3]:0.4+0.95=1.35 [3,4]:0.95+0.1=1.05
	 * Best is [2,3] → indices {2, 3}. */
	const std::vector<double> q = {0.1, 0.9, 0.4, 0.95, 0.1};
	auto best = mpp::align_find_best_frames(q, 2, 2);
	cr_assert_eq(best.size(), 2u);
	std::sort(best.begin(), best.end());
	cr_assert(best == (std::vector<int>{2, 3}),
	          "expected {2,3}, got {%d, %d}", best[0], best[1]);

	/* Window covering whole sequence should return the global top-2 ({1,3}). */
	auto best_full = mpp::align_find_best_frames(q, 2, (int) q.size());
	std::sort(best_full.begin(), best_full.end());
	cr_assert(best_full == (std::vector<int>{1, 3}));
}

Test(mpp_align, average_frame_basic) {
	const auto cfg = default_cfg();
	const cv::Mat truth = make_textured_frame();
	const std::vector<std::pair<int, int>> jit = {
	    {0, 0}, {-2, 1}, {3, -2}, {1, 2}, {-1, -1}, {0, 3}};
	std::vector<cv::Mat> frames;
	for (auto p : jit) frames.push_back(shifted(truth, p.first, p.second));
	const std::vector<double> q = {1.0, 0.7, 0.5, 0.85, 0.6, 0.55};
	std::vector<cv::Vec2d> shifts;
	for (auto p : jit) shifts.push_back(cv::Vec2d(-p.first, -p.second));

	auto cfg_no_fc = cfg;
	cfg_no_fc.align_frames_fast_changing_object = false;
	const auto r = mpp::align_average_frame(frames, q, shifts, cfg_no_fc);

	cr_assert_gt(r.mean_frame.rows, 0);
	cr_assert_gt(r.mean_frame.cols, 0);
	cr_assert_eq(r.mean_frame.type(), CV_32S);
	/* The intersection must cover the common region across all jittered
	 * frames; for {-2..3, -2..3} shifts and 200x240 frames, the intersection
	 * shape is (200 - 5) × (240 - 5) = 195 × 235. */
	cr_assert_eq(r.mean_frame.rows, 200 - 5);
	cr_assert_eq(r.mean_frame.cols, 240 - 5);
	const int n_target = std::max((int) std::ceil(6 * cfg.align_frames_average_frame_percent / 100.0), 1);
	cr_assert_eq((int) r.indices_used.size(), n_target);
}

