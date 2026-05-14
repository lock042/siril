/*
 * Phase 1 — mpp_rank Criterion tests.
 *
 * Covers PSS-faithful Laplace-σ ranking against the private C++ entry points,
 * plus the C-level mpp_rank_sequence stub (NB: sequence integration ships in
 * Phase 1.3; the stub still returns ENOTIMPL).
 */
#include <criterion/criterion.h>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "registration/mpp_rank_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_config.h"
#include "registration/mpp_rank.h"
}

namespace {

/* Synthesise a textured 16-bit mono frame with a bright disk + speckles, so
 * the Laplacian has plenty of structure to measure. Deterministic seed. */
cv::Mat make_sharp_frame(int rows = 200, int cols = 240, int seed = 1) {
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
	/* Sprinkle bright speckles to feed the Laplacian. */
	rng.state = seed;
	for (int i = 0; i < 40; ++i) {
		const int sy = rng.uniform(int(cy - r0 * 0.6), int(cy + r0 * 0.6));
		const int sx = rng.uniform(int(cx - r0 * 0.6), int(cx + r0 * 0.6));
		cv::circle(f, {sx, sy}, 2, cv::Scalar(60000), -1);
	}
	return f;
}

mpp_config_t default_cfg() {
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	return cfg;
}

}  // namespace

/* ------------------------------------------------------------------------- */

Test(mpp_rank, stub_sequence_returns_enotimpl) {
	double *q = nullptr;
	cr_assert_eq(mpp_rank_sequence(nullptr, nullptr, &q), MPP_ENOTIMPL);
	cr_assert_null(q);
}

Test(mpp_rank, default_config_has_pss_values) {
	const auto cfg = default_cfg();
	cr_assert_eq(cfg.frames_gauss_width, 7);
	cr_assert_eq(cfg.align_frames_sampling_stride, 2);
	cr_assert_float_eq(cfg.rank_laplacian_alpha, 1.0 / 256.0, 1e-12);
	cr_assert(cfg.frames_normalization);
	cr_assert_eq(cfg.frames_normalization_threshold, 15);
	cr_assert_eq(cfg.bitpix, 16);
}

Test(mpp_rank, score_is_positive_on_sharp_frame) {
	const auto cfg = default_cfg();
	const cv::Mat f = make_sharp_frame();
	const double s = mpp::rank_score_mat(f, cfg);
	cr_assert_gt(s, 0.0, "σ should be > 0 on a textured frame, got %.6g", s);
}

/* PSS sanity: blurring a frame more must reduce its Laplace-σ. */
Test(mpp_rank, blur_decreases_quality) {
	const auto cfg = default_cfg();
	const cv::Mat sharp = make_sharp_frame();
	cv::Mat soft;
	cv::GaussianBlur(sharp, soft, cv::Size(13, 13), 0);  /* extra blur beyond the 7×7 PSS uses */

	const double q_sharp = mpp::rank_score_mat(sharp, cfg);
	const double q_soft  = mpp::rank_score_mat(soft,  cfg);

	cr_assert_gt(q_sharp, q_soft,
	             "sharp σ (%.4f) should exceed softened σ (%.4f)",
	             q_sharp, q_soft);
}

/* PSS sanity: the normalized rank should be approximately invariant to a
 * uniform scaling of the input. The raw σ does scale with brightness, but
 * `rank_score_normalized` divides it back out. We test with `THRESH_TOZERO`-
 * style background to mimic PSS's brightness measure. */
Test(mpp_rank, brightness_normalization_is_scale_invariant) {
	auto cfg = default_cfg();
	cfg.frames_normalization = true;

	const cv::Mat f = make_sharp_frame();
	cv::Mat dim;
	f.convertTo(dim, CV_16U, 0.5);  /* halve every pixel */

	const double q_full = mpp::rank_score_normalized(f,   cfg);
	const double q_dim  = mpp::rank_score_normalized(dim, cfg);

	/* Quantisation through convertScaleAbs+uint8 means the equivalence isn't
	 * exact; we just want the dynamic range to be tamed compared to the raw σ. */
	const double q_full_raw = mpp::rank_score_mat(f,   cfg);
	const double q_dim_raw  = mpp::rank_score_mat(dim, cfg);

	const double normalized_rel = std::abs(q_full - q_dim) / q_full;
	const double raw_rel        = std::abs(q_full_raw - q_dim_raw) / q_full_raw;

	cr_assert_lt(normalized_rel, raw_rel,
	             "normalization should reduce the brightness sensitivity "
	             "(normalized Δrel=%.4f, raw Δrel=%.4f)",
	             normalized_rel, raw_rel);
}

/* Oracle test: rank the bundled synthetic dataset and compare to the
 * pre-computed PSS quality vector (tools/pss_reference/oracle_out/quality.csv).
 * Acceptance bar from pss_port_plan.md §1.4: per-frame relative error < 1e-4. */
Test(mpp_rank, oracle_equivalence_synthetic) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) {
		cr_skip_test("MPP_PSS_TEST_DATA_DIR unset — skipping oracle test "
		             "(meson sets this automatically; run via `meson test`)");
	}
	namespace fs = std::filesystem;
	const fs::path frames_dir = fs::path(data_root) / "test_data" / "synth_planet";
	const fs::path oracle_csv = fs::path(data_root) / "oracle_out" / "quality.csv";
	if (!fs::exists(frames_dir) || !fs::exists(oracle_csv)) {
		cr_skip_test("oracle artifacts not present at %s — regenerate via "
		             "tools/pss_reference/run_pss.py", data_root);
	}

	/* Collect frame_XXX.png paths in sorted order. */
	std::vector<fs::path> frame_paths;
	for (const auto &entry : fs::directory_iterator(frames_dir)) {
		const std::string n = entry.path().filename().string();
		if (n.rfind("frame_", 0) == 0 && entry.path().extension() == ".png")
			frame_paths.push_back(entry.path());
	}
	std::sort(frame_paths.begin(), frame_paths.end());
	cr_assert_gt(frame_paths.size(), 0u, "no frame_*.png found in %s",
	             frames_dir.c_str());

	/* Compute σ for each frame. PSS's ranking path matches our default cfg. */
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	std::vector<double> sigmas;
	sigmas.reserve(frame_paths.size());
	for (const auto &p : frame_paths) {
		cv::Mat f = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		cr_assert_eq(f.depth(), CV_16U,
		             "expected 16-bit synth PNGs, got depth %d for %s",
		             f.depth(), p.filename().c_str());
		/* PSS rank_score_normalized() applies brightness norm when
		 * frames_normalization is true — and the oracle was generated with
		 * the PSS default (true). */
		sigmas.push_back(mpp::rank_score_normalized(f, cfg));
	}

	/* PSS divides every value by the max, so reconstruct that here. */
	const double max_sigma = *std::max_element(sigmas.begin(), sigmas.end());
	cr_assert_gt(max_sigma, 0.0);
	std::vector<double> q_ours(sigmas.size());
	for (size_t i = 0; i < sigmas.size(); ++i)
		q_ours[i] = sigmas[i] / max_sigma;

	/* Read oracle quality.csv (one float per line). */
	std::vector<double> q_oracle;
	{
		std::ifstream fin(oracle_csv);
		cr_assert(fin.good(), "could not open %s", oracle_csv.c_str());
		double v;
		while (fin >> v) q_oracle.push_back(v);
	}
	cr_assert_eq(q_ours.size(), q_oracle.size(),
	             "frame count mismatch: ours=%zu oracle=%zu",
	             q_ours.size(), q_oracle.size());

	/* Compare per-frame. The acceptance bar is 1e-4 relative; we report the
	 * worst frame on failure for fast triage. */
	const double tol = 1e-4;
	double worst_rel = 0.0;
	int worst_i = -1;
	for (size_t i = 0; i < q_ours.size(); ++i) {
		const double rel = std::abs(q_ours[i] - q_oracle[i]) / q_oracle[i];
		if (rel > worst_rel) { worst_rel = rel; worst_i = (int) i; }
	}
	cr_assert_lt(worst_rel, tol,
	             "rank mismatch vs PSS at frame %d: ours=%.10g oracle=%.10g "
	             "(rel %.3e > %.3e)",
	             worst_i, q_ours[worst_i], q_oracle[worst_i], worst_rel, tol);
}

/* PSS uses 16-bit input by default; threshold scales by 256 in that case. */
Test(mpp_rank, brightness_threshold_scales_with_bitpix) {
	const cv::Mat f = make_sharp_frame();

	auto cfg16 = default_cfg();
	cfg16.bitpix = 16;
	cfg16.frames_normalization_threshold = 15;
	const double ab16 = mpp::rank_average_brightness(f, cfg16);

	auto cfg8 = default_cfg();
	cfg8.bitpix = 8;
	cfg8.frames_normalization_threshold = 15;
	const double ab8 = mpp::rank_average_brightness(f, cfg8);

	/* For 16-bit input, thr16 = 15 × 256 = 3840 in the 16-bit range; thr8 = 15
	 * (essentially no rejection at 16-bit scale). The 8-bit-threshold reading
	 * should therefore be ≥ the 16-bit-threshold reading. */
	cr_assert_geq(ab8, ab16,
	              "ab(thr=15)=%.2f should be ≥ ab(thr=3840)=%.2f", ab8, ab16);
}
