/*
 * Phase 5a — mpp_stack Criterion tests.
 *
 * Stacking is broken into separately-testable pieces; this file builds out
 * coverage in lock-step with the implementation so divergences from PSS
 * surface at the smallest unit they appear at.
 */
#include <criterion/criterion.h>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "registration/mpp_align_priv.hpp"
#include "registration/mpp_ap_priv.hpp"
#include "registration/mpp_rank_priv.hpp"
#include "registration/mpp_shift_priv.hpp"
#include "registration/mpp_stack_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_config.h"
#include "registration/mpp_stack.h"
}

namespace {

mpp_config_t default_cfg() {
	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	return cfg;
}

cv::Mat blurred(const cv::Mat &mono, const mpp_config_t &cfg) {
	cv::Mat b;
	cv::GaussianBlur(mono, b, cv::Size(cfg.frames_gauss_width, cfg.frames_gauss_width), 0);
	return b;
}

std::vector<double> read_doubles_csv(const std::filesystem::path &p) {
	std::vector<double> v;
	std::ifstream fin(p);
	double x;
	while (fin >> x) v.push_back(x);
	return v;
}

std::vector<int> read_ints_csv(const std::filesystem::path &p) {
	std::vector<int> v;
	std::ifstream fin(p);
	int x;
	while (fin >> x) v.push_back(x);
	return v;
}

int read_int_text(const std::filesystem::path &p) {
	std::ifstream fin(p);
	int x = -1;
	fin >> x;
	return x;
}

}  // namespace

/* ------------------------------------------------------------------------- */

Test(mpp_stack, stub_returns_enotimpl) {
	cr_assert_eq(mpp_stack(nullptr, nullptr, nullptr, nullptr, nullptr,
	                       MPP_RESAMPLE_BICUBIC, 1, nullptr),
	             MPP_ENOTIMPL);
}

Test(mpp_stack, default_config_phase5a) {
	const auto cfg = default_cfg();
	cr_assert_eq(cfg.alignment_points_frame_percent, 10);
	cr_assert_eq(cfg.alignment_points_frame_number, -1);
	cr_assert_eq(cfg.alignment_points_rank_pixel_stride, 2);
	cr_assert(cfg.alignment_points_de_warp);
	cr_assert_float_eq(cfg.alignment_points_penalty_factor, 0.00025, 1e-12);
	cr_assert_float_eq(cfg.stack_frames_background_fraction, 0.3, 1e-12);
	cr_assert_float_eq(cfg.stack_frames_background_blend_threshold, 0.2, 1e-12);
	cr_assert_eq(cfg.stack_frames_background_patch_size, 100);
	cr_assert_eq(cfg.drizzle_factor, 1);
}

/* PSS one_dim_weight reference values for a symmetric ramp:
 *   patch=[100, 172], box_center=136, no extend
 *   center_offset = 36; high_len = patch_size - center_offset = 36
 *   weights[0] = 1/37 ≈ 0.027027
 *   weights[35] = 36/37 ≈ 0.972973
 *   weights[36] = 1.0   (the centre — first index of high ramp)
 *   weights[37] = 35/36 ≈ 0.972222
 *   weights[71] = 1/36  ≈ 0.027778
 */
Test(mpp_stack, one_dim_weight_symmetric_ramp) {
	const auto w = mpp::stack_one_dim_weight(100, 172, 136, false, false);
	cr_assert_eq(w.size(), 72u);
	cr_assert_float_eq(w[0],  1.0f / 37.0f, 1e-6);
	cr_assert_float_eq(w[35], 36.0f / 37.0f, 1e-6);
	cr_assert_float_eq(w[36], 1.0f, 1e-6);
	cr_assert_float_eq(w[37], 35.0f / 36.0f, 1e-6);
	cr_assert_float_eq(w[71], 1.0f / 36.0f, 1e-6);
}

Test(mpp_stack, one_dim_weight_extend_low) {
	const auto w = mpp::stack_one_dim_weight(0, 50, 30, true, false);
	cr_assert_eq(w.size(), 50u);
	/* Lower ramp replaced by 1.0 from idx 0 to centre_offset (=30, exclusive). */
	for (int i = 0; i < 30; ++i)
		cr_assert_float_eq(w[i], 1.0f, 1e-9, "w[%d] should be 1.0", i);
	/* High ramp at idx 30 is 1.0; at idx 49 is 1/20. */
	cr_assert_float_eq(w[30], 1.0f, 1e-9);
	cr_assert_float_eq(w[49], 1.0f / 20.0f, 1e-6);
}

Test(mpp_stack, one_dim_weight_extend_high) {
	const auto w = mpp::stack_one_dim_weight(0, 50, 20, false, true);
	cr_assert_eq(w.size(), 50u);
	/* Low ramp: [1/21, 2/21, …, 20/21]. */
	cr_assert_float_eq(w[0], 1.0f / 21.0f, 1e-6);
	cr_assert_float_eq(w[19], 20.0f / 21.0f, 1e-6);
	/* High ramp replaced by 1.0 from idx 20 to end. */
	for (int i = 20; i < 50; ++i)
		cr_assert_float_eq(w[i], 1.0f, 1e-9, "w[%d] should be 1.0", i);
}

/* weight_matrix_first_phase: centre is exactly 1.0; the four corners are
 * 1 − penalty × 2 since (x/sw1 − 1)² + (y/sw1 − 1)² = 1 + 1 = 2 at
 * (0,0). For default search_width=14: sw1 = 5, extent = 11. */
Test(mpp_stack, weight_matrix_centre_and_corners) {
	const auto cfg = default_cfg();
	const cv::Mat m = mpp::stack_build_first_phase_weight_matrix(cfg);
	cr_assert_eq(m.rows, 11);
	cr_assert_eq(m.cols, 11);
	cr_assert_float_eq(m.at<float>(5, 5), 1.0f, 1e-9);

	const double pen = cfg.alignment_points_penalty_factor;
	const float expected_corner = (float) (1.0 - pen * 2.0);
	cr_assert_float_eq(m.at<float>(0, 0),  expected_corner, 1e-7);
	cr_assert_float_eq(m.at<float>(0, 10), expected_corner, 1e-7);
	cr_assert_float_eq(m.at<float>(10, 0), expected_corner, 1e-7);
	cr_assert_float_eq(m.at<float>(10, 10), expected_corner, 1e-7);

	/* (0, 5): x at centre (no x penalty), y at 0 (full y penalty of 1). */
	const float expected_edge = (float) (1.0 - pen * 1.0);
	cr_assert_float_eq(m.at<float>(0, 5), expected_edge, 1e-7);
}

/* remap_rigid: a 10×10 patch shifted (0, 0) lands exactly. */
Test(mpp_stack, remap_rigid_no_shift) {
	cv::Mat frame(50, 50, CV_32F, cv::Scalar(0));
	for (int y = 20; y < 30; ++y)
		for (int x = 20; x < 30; ++x)
			frame.at<float>(y, x) = (float) (y * 100 + x);
	cv::Mat buffer(10, 10, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	mpp::stack_remap_rigid(frame, buffer, 0, 0, 20, 30, 20, 30, b);

	for (int y = 0; y < 10; ++y)
		for (int x = 0; x < 10; ++x)
			cr_assert_float_eq(buffer.at<float>(y, x),
			                   (float) ((y + 20) * 100 + (x + 20)), 1e-6);
	cr_assert_eq(b.y_low,  0); cr_assert_eq(b.y_high, 0);
	cr_assert_eq(b.x_low,  0); cr_assert_eq(b.x_high, 0);
}

/* remap_rigid accumulates: calling twice doubles the contribution. */
Test(mpp_stack, remap_rigid_accumulates) {
	cv::Mat frame(20, 20, CV_32F, cv::Scalar(5.0f));
	cv::Mat buffer(8, 8, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	mpp::stack_remap_rigid(frame, buffer, 0, 0, 5, 13, 5, 13, b);
	mpp::stack_remap_rigid(frame, buffer, 0, 0, 5, 13, 5, 13, b);
	cr_assert_float_eq(buffer.at<float>(4, 4), 10.0f, 1e-6);
}

/* Oracle: per-AP per-frame quality matrix and stack_size match PSS on the
 * bundled synthetic dataset. */
Test(mpp_stack, ap_compute_frame_qualities_oracle) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	namespace fs = std::filesystem;
	const fs::path frames_dir = fs::path(data_root) / "test_data" / "synth_planet";
	const fs::path q_csv      = fs::path(data_root) / "oracle_out" / "ap_qualities.csv";
	const fs::path bfi_csv    = fs::path(data_root) / "oracle_out" / "best_frame_indices.csv";
	const fs::path bri_csv    = fs::path(data_root) / "oracle_out" / "frame_brightness.csv";
	const fs::path ss_txt     = fs::path(data_root) / "oracle_out" / "stack_size.txt";
	const fs::path ua_csv     = fs::path(data_root) / "oracle_out" / "used_alignment_points.csv";
	if (!fs::exists(q_csv) || !fs::exists(bfi_csv) || !fs::exists(bri_csv)
	 || !fs::exists(ss_txt) || !fs::exists(ua_csv))
		cr_skip_test("Phase 5a oracle artifacts missing — regenerate via run_pss.py");

	std::vector<fs::path> paths;
	for (const auto &e : fs::directory_iterator(frames_dir)) {
		const std::string n = e.path().filename().string();
		if (n.rfind("frame_", 0) == 0 && e.path().extension() == ".png")
			paths.push_back(e.path());
	}
	std::sort(paths.begin(), paths.end());

	const auto cfg = default_cfg();
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q_rank;
	for (const auto &p : paths) {
		cv::Mat m = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		frames_raw.push_back(m);
		frames_blurred.push_back(blurred(m, cfg));
		q_rank.push_back(mpp::rank_score_normalized(m, cfg));
	}
	const auto align = mpp::align_global_from_frames(frames_blurred, q_rank, cfg);
	const auto avg   = mpp::align_average_frame(frames_raw, q_rank, align.shifts, cfg);
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	cr_assert_not_null(aps);
	const auto offsets = mpp::shift_frame_offsets(align.shifts, avg.intersection);

	/* PSS frame_brightness comes from frames.frames_mono — our equivalent
	 * is rank_average_brightness on the raw mono. */
	std::vector<double> frame_brightness;
	for (const auto &m : frames_raw)
		frame_brightness.push_back(mpp::rank_average_brightness(m, cfg));

	const int H = frames_raw[0].rows, W = frames_raw[0].cols;
	const auto apq = mpp::ap_compute_frame_qualities(frames_raw, frame_brightness,
	                                                 *aps, offsets, H, W, cfg);

	/* Compare stack_size. */
	const int expected_stack_size = read_int_text(ss_txt);
	cr_assert_eq(apq.stack_size, expected_stack_size,
	             "stack_size mismatch: ours=%d oracle=%d",
	             apq.stack_size, expected_stack_size);

	/* Frame brightness — sanity vs oracle (PSS uses cv::mean+THRESH_TOZERO+1e-10,
	 * we use the same; should match within float precision). */
	const auto oracle_brightness = read_doubles_csv(bri_csv);
	cr_assert_eq(oracle_brightness.size(), frame_brightness.size());
	for (size_t i = 0; i < frame_brightness.size(); ++i) {
		const double rel = std::abs(frame_brightness[i] - oracle_brightness[i])
		                 / oracle_brightness[i];
		cr_assert_lt(rel, 1e-9,
		             "frame %zu brightness rel diff %.3e: ours=%.6f oracle=%.6f",
		             i, rel, frame_brightness[i], oracle_brightness[i]);
	}

	/* Per-AP per-frame qualities — matrix shape (M × N). */
	const auto q_flat = read_doubles_csv(q_csv);
	cr_assert_eq(q_flat.size(), (size_t) (aps->count * frames_raw.size()),
	             "quality matrix size mismatch: ours=%zu*%zu oracle=%zu",
	             (size_t) aps->count, frames_raw.size(), q_flat.size());
	int worst_a = -1, worst_f = -1;
	double worst_rel = 0.0;
	for (int a = 0; a < aps->count; ++a) {
		for (size_t f = 0; f < frames_raw.size(); ++f) {
			const double o = q_flat[(size_t) a * frames_raw.size() + f];
			const double rel = o == 0.0
			    ? std::abs(apq.qualities[a][f])
			    : std::abs(apq.qualities[a][f] - o) / std::abs(o);
			if (rel > worst_rel) { worst_rel = rel; worst_a = a; worst_f = (int) f; }
		}
	}
	cr_assert_lt(worst_rel, 1e-9,
	             "ap_qualities mismatch at AP %d frame %d: ours=%.10g oracle=%.10g (rel %.3e)",
	             worst_a, worst_f, apq.qualities[worst_a][worst_f],
	             q_flat[(size_t) worst_a * frames_raw.size() + worst_f], worst_rel);

	/* best_frame_indices: PSS sorts by quality desc, stable on ties. We use
	 * std::partial_sort which is not guaranteed stable, but for inputs with
	 * strictly distinct qualities (which our synthetic data delivers) the
	 * result is unique. Compare element-wise. */
	const auto bfi_flat = read_ints_csv(bfi_csv);
	cr_assert_eq(bfi_flat.size(), (size_t) (aps->count * apq.stack_size));
	for (int a = 0; a < aps->count; ++a) {
		for (int k = 0; k < apq.stack_size; ++k) {
			const int o = bfi_flat[(size_t) a * apq.stack_size + k];
			cr_assert_eq(apq.best_frame_indices[a][k], o,
			             "best_frame_indices mismatch at AP %d k=%d: ours=%d oracle=%d",
			             a, k, apq.best_frame_indices[a][k], o);
		}
	}

	/* used_alignment_points: jagged per-frame list (whitespace-separated). */
	std::ifstream ua_in(ua_csv);
	std::string line;
	int fidx = 0;
	while (std::getline(ua_in, line)) {
		std::vector<int> oracle_list;
		std::stringstream ss(line);
		int v;
		while (ss >> v) oracle_list.push_back(v);
		const auto &ours = apq.used_alignment_points[fidx];
		std::vector<int> ours_sorted = ours;
		std::sort(ours_sorted.begin(), ours_sorted.end());
		cr_assert_eq(ours_sorted.size(), oracle_list.size(),
		             "frame %d used_aps count: ours=%zu oracle=%zu",
		             fidx, ours_sorted.size(), oracle_list.size());
		for (size_t k = 0; k < oracle_list.size(); ++k)
			cr_assert_eq(ours_sorted[k], oracle_list[k],
			             "frame %d used_aps[%zu]: ours=%d oracle=%d",
			             fidx, k, ours_sorted[k], oracle_list[k]);
		++fidx;
	}
	mpp_ap_free(aps);
}

/* Oracle: prepare_for_stack_blending must produce the same
 * sum_single_frame_weights buffer as PSS, cell-for-cell. Also matches
 * number_stacking_holes. */
Test(mpp_stack, prepare_for_blending_oracle) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	namespace fs = std::filesystem;
	const fs::path frames_dir  = fs::path(data_root) / "test_data" / "synth_planet";
	const fs::path sw_bin      = fs::path(data_root) / "oracle_out" / "sum_single_frame_weights.bin";
	const fs::path sw_dims_csv = fs::path(data_root) / "oracle_out" / "sum_single_frame_weights_dims.csv";
	const fs::path holes_txt   = fs::path(data_root) / "oracle_out" / "number_stacking_holes.txt";
	const fs::path ss_txt      = fs::path(data_root) / "oracle_out" / "stack_size.txt";
	if (!fs::exists(sw_bin) || !fs::exists(sw_dims_csv)
	 || !fs::exists(holes_txt) || !fs::exists(ss_txt))
		cr_skip_test("oracle artifacts missing — regenerate via run_pss.py");

	std::vector<fs::path> paths;
	for (const auto &e : fs::directory_iterator(frames_dir)) {
		const std::string n = e.path().filename().string();
		if (n.rfind("frame_", 0) == 0 && e.path().extension() == ".png")
			paths.push_back(e.path());
	}
	std::sort(paths.begin(), paths.end());

	const auto cfg = default_cfg();
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q_rank;
	for (const auto &p : paths) {
		cv::Mat m = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		frames_raw.push_back(m);
		frames_blurred.push_back(blurred(m, cfg));
		q_rank.push_back(mpp::rank_score_normalized(m, cfg));
	}
	const auto align = mpp::align_global_from_frames(frames_blurred, q_rank, cfg);
	const auto avg   = mpp::align_average_frame(frames_raw, q_rank, align.shifts, cfg);
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	cr_assert_not_null(aps);
	const int stack_size = read_int_text(ss_txt);
	const auto state = mpp::stack_prepare_for_blending(*aps, avg.intersection,
	                                                   stack_size, cfg.drizzle_factor, cfg);

	const auto dims = read_ints_csv(sw_dims_csv);
	cr_assert_eq(dims.size(), 2u);
	const int H = dims[0], W = dims[1];
	cr_assert_eq(state.sum_single_frame_weights.rows, H);
	cr_assert_eq(state.sum_single_frame_weights.cols, W);

	cv::Mat oracle(H, W, CV_32F);
	std::ifstream fin(sw_bin, std::ios::binary);
	fin.read(reinterpret_cast<char *>(oracle.data),
	         (std::streamsize) H * W * sizeof(float));
	cr_assert(fin.good(), "failed to read sum_single_frame_weights.bin");

	float worst = 0.0f;
	int worst_y = -1, worst_x = -1;
	for (int y = 0; y < H; ++y) {
		const float *o = oracle.ptr<float>(y);
		const float *m = state.sum_single_frame_weights.ptr<float>(y);
		for (int x = 0; x < W; ++x) {
			const float d = std::abs(m[x] - o[x]);
			if (d > worst) { worst = d; worst_y = y; worst_x = x; }
		}
	}
	cr_assert_lt(worst, 1e-5f,
	             "sum_single_frame_weights diff at (%d, %d): ours=%.6f oracle=%.6f",
	             worst_y, worst_x,
	             state.sum_single_frame_weights.at<float>(worst_y, worst_x),
	             oracle.at<float>(worst_y, worst_x));

	const int oracle_holes = read_int_text(holes_txt);
	cr_assert_eq(state.number_stacking_holes, oracle_holes,
	             "number_stacking_holes mismatch: ours=%d oracle=%d",
	             state.number_stacking_holes, oracle_holes);
	mpp_ap_free(aps);
}

/* Oracle: full stack_frames_loop produces the same per-AP stacking_buffer,
 * averaged_background, and border counts as PSS on the synthetic dataset.
 *
 * This is the deepest oracle yet — it exercises everything from Phase 1
 * (quality), Phase 2 (global alignment + average frame), Phase 3 (APs),
 * Phase 4's correlation kernel (now with the penalty weight matrix), and
 * the new stacking loop. A divergence anywhere upstream will show up
 * here as a per-AP buffer mismatch. */
Test(mpp_stack, stack_frames_loop_oracle) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	namespace fs = std::filesystem;
	const fs::path frames_dir   = fs::path(data_root) / "test_data" / "synth_planet";
	const fs::path ap_bufs_bin  = fs::path(data_root) / "oracle_out" / "ap_stacking_buffers.bin";
	const fs::path ap_bufs_meta = fs::path(data_root) / "oracle_out" / "ap_stacking_buffers_meta.csv";
	const fs::path bg_patches_csv = fs::path(data_root) / "oracle_out" / "background_patches.csv";
	const fs::path borders_csv  = fs::path(data_root) / "oracle_out" / "stack_borders.csv";
	const fs::path avg_bg_bin   = fs::path(data_root) / "oracle_out" / "averaged_background.bin";
	const fs::path avg_bg_dims  = fs::path(data_root) / "oracle_out" / "averaged_background_dims.csv";
	if (!fs::exists(ap_bufs_bin) || !fs::exists(ap_bufs_meta)
	 || !fs::exists(borders_csv))
		cr_skip_test("oracle artifacts missing — regenerate via run_pss.py");

	std::vector<fs::path> paths;
	for (const auto &e : fs::directory_iterator(frames_dir)) {
		const std::string n = e.path().filename().string();
		if (n.rfind("frame_", 0) == 0 && e.path().extension() == ".png")
			paths.push_back(e.path());
	}
	std::sort(paths.begin(), paths.end());

	const auto cfg = default_cfg();
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q_rank;
	for (const auto &p : paths) {
		cv::Mat m = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		frames_raw.push_back(m);
		frames_blurred.push_back(blurred(m, cfg));
		q_rank.push_back(mpp::rank_score_normalized(m, cfg));
	}
	const auto align = mpp::align_global_from_frames(frames_blurred, q_rank, cfg);
	const auto avg   = mpp::align_average_frame(frames_raw, q_rank, align.shifts, cfg);
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	cr_assert_not_null(aps);
	const auto offsets = mpp::shift_frame_offsets(align.shifts, avg.intersection);

	std::vector<double> frame_brightness;
	for (const auto &m : frames_raw)
		frame_brightness.push_back(mpp::rank_average_brightness(m, cfg));

	const int H = frames_raw[0].rows, W = frames_raw[0].cols;
	const auto apq = mpp::ap_compute_frame_qualities(frames_raw, frame_brightness,
	                                                 *aps, offsets, H, W, cfg);

	/* Global quality sorted indices (descending). */
	std::vector<int> sorted_idx(frames_raw.size());
	std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
	std::sort(sorted_idx.begin(), sorted_idx.end(),
	          [&](int a, int b) { return q_rank[a] > q_rank[b]; });

	const auto loop = mpp::stack_frames_loop(frames_raw, frames_blurred,
	                                         avg.mean_frame, *aps, apq,
	                                         offsets, frame_brightness,
	                                         sorted_idx, avg.intersection, cfg);

	/* Background patches. */
	std::ifstream bg_in(bg_patches_csv);
	std::vector<mpp::StackBackgroundPatch> oracle_bg;
	{
		std::string line;
		while (std::getline(bg_in, line)) {
			if (line.empty()) continue;
			std::stringstream ss(line);
			mpp::StackBackgroundPatch p{};
			ss >> p.y_low >> p.y_high >> p.x_low >> p.x_high;
			oracle_bg.push_back(p);
		}
	}
	cr_assert_eq(loop.state.background_patches.size(), oracle_bg.size(),
	             "background_patches count: ours=%zu oracle=%zu",
	             loop.state.background_patches.size(), oracle_bg.size());
	for (size_t i = 0; i < oracle_bg.size(); ++i) {
		const auto &o = oracle_bg[i];
		const auto &m = loop.state.background_patches[i];
		cr_assert(m.y_low == o.y_low && m.y_high == o.y_high
		       && m.x_low == o.x_low && m.x_high == o.x_high,
		          "bg patch %zu: ours=(%d,%d,%d,%d) oracle=(%d,%d,%d,%d)",
		          i, m.y_low, m.y_high, m.x_low, m.x_high,
		          o.y_low, o.y_high, o.x_low, o.x_high);
	}

	/* Border counts. */
	{
		const auto b = read_ints_csv(borders_csv);
		cr_assert_eq(b.size(), 4u);
		cr_assert_eq(loop.border.y_low,  b[0]);
		cr_assert_eq(loop.border.y_high, b[1]);
		cr_assert_eq(loop.border.x_low,  b[2]);
		cr_assert_eq(loop.border.x_high, b[3]);
	}

	/* Per-AP stacking buffers. The oracle meta CSV has one row per AP
	 * with (rows, cols, byte_offset). Each AP's buffer is rows*cols*float32. */
	const auto meta_flat = read_ints_csv(ap_bufs_meta);
	cr_assert_eq(meta_flat.size(), (size_t) (aps->count * 3));
	std::ifstream bufs_in(ap_bufs_bin, std::ios::binary);
	cr_assert(bufs_in.good(), "cannot open ap_stacking_buffers.bin");
	double worst_rel = 0.0;
	int worst_a = -1;
	double worst_my = 0.0, worst_oracle = 0.0;
	for (int a = 0; a < aps->count; ++a) {
		const int rows = (int) meta_flat[3 * a + 0];
		const int cols = (int) meta_flat[3 * a + 1];
		cv::Mat oracle(rows, cols, CV_32F);
		bufs_in.read(reinterpret_cast<char *>(oracle.data),
		             (std::streamsize) rows * cols * sizeof(float));
		cr_assert(bufs_in.good(), "buffer read failed for AP %d", a);
		const auto &mine = loop.state.stacking_buffers[a];
		cr_assert_eq(mine.rows, rows, "AP %d row count: %d vs %d", a, mine.rows, rows);
		cr_assert_eq(mine.cols, cols, "AP %d col count: %d vs %d", a, mine.cols, cols);
		for (int y = 0; y < rows; ++y) {
			const float *o = oracle.ptr<float>(y);
			const float *m = mine.ptr<float>(y);
			for (int x = 0; x < cols; ++x) {
				const float ov = o[x];
				const float mv = m[x];
				const double rel = std::abs(ov) > 1e-3
				    ? std::abs(mv - ov) / std::abs(ov)
				    : std::abs(mv - ov);
				if (rel > worst_rel) {
					worst_rel = rel; worst_a = a;
					worst_my = mv; worst_oracle = ov;
				}
			}
		}
	}
	cr_assert_lt(worst_rel, 1e-3,
	             "per-AP buffer mismatch at AP %d: ours=%.6f oracle=%.6f (rel %.3e)",
	             worst_a, worst_my, worst_oracle, worst_rel);

	/* averaged_background. */
	if (fs::exists(avg_bg_bin) && fs::exists(avg_bg_dims)) {
		const auto dims = read_ints_csv(avg_bg_dims);
		cr_assert_eq(dims.size(), 2u);
		const int bH = dims[0], bW = dims[1];
		/* PSS dumps the *post-loop* averaged_background, i.e., after the
		 * `/= stack_size` step. Our state.averaged_background underwent the
		 * same step. */
		cv::Mat ob(bH, bW, CV_32F);
		std::ifstream fin(avg_bg_bin, std::ios::binary);
		fin.read(reinterpret_cast<char *>(ob.data),
		         (std::streamsize) bH * bW * sizeof(float));
		cr_assert(fin.good(), "averaged_background.bin read failed");
		cr_assert_eq(loop.state.averaged_background.rows, bH);
		cr_assert_eq(loop.state.averaged_background.cols, bW);
		double bg_worst_rel = 0.0;
		int by = -1, bx = -1;
		for (int y = 0; y < bH; ++y) {
			const float *o = ob.ptr<float>(y);
			const float *m = loop.state.averaged_background.ptr<float>(y);
			for (int x = 0; x < bW; ++x) {
				const double rel = std::abs(o[x]) > 1e-3
				    ? std::abs(m[x] - o[x]) / std::abs(o[x])
				    : std::abs(m[x] - o[x]);
				if (rel > bg_worst_rel) { bg_worst_rel = rel; by = y; bx = x; }
			}
		}
		cr_assert_lt(bg_worst_rel, 1e-3,
		             "avg_background diff at (%d,%d): ours=%.6f oracle=%.6f (rel %.3e)",
		             by, bx,
		             loop.state.averaged_background.at<float>(by, bx),
		             ob.at<float>(by, bx), bg_worst_rel);
	}

	mpp_ap_free(aps);
}

/* Phase 5a end-to-end oracle: the final uint16 stacked image must match
 * PSS's `stack.merge_alignment_point_buffers()` output cell-for-cell on the
 * synthetic dataset. Acceptance bar: every pixel within ≤2 levels (out of
 * 65535) of PSS, with at least 99% exact matches. */
Test(mpp_stack, end_to_end_stacked_oracle) {
	const char *data_root = std::getenv("MPP_PSS_TEST_DATA_DIR");
	if (!data_root) cr_skip_test("MPP_PSS_TEST_DATA_DIR unset");
	namespace fs = std::filesystem;
	const fs::path frames_dir = fs::path(data_root) / "test_data" / "synth_planet";
	const fs::path stacked_bin = fs::path(data_root) / "oracle_out" / "stacked_u16.bin";
	const fs::path stacked_dims = fs::path(data_root) / "oracle_out" / "stacked_u16_dims.csv";
	if (!fs::exists(stacked_bin) || !fs::exists(stacked_dims))
		cr_skip_test("stacked oracle missing — regenerate via run_pss.py");

	std::vector<fs::path> paths;
	for (const auto &e : fs::directory_iterator(frames_dir)) {
		const std::string n = e.path().filename().string();
		if (n.rfind("frame_", 0) == 0 && e.path().extension() == ".png")
			paths.push_back(e.path());
	}
	std::sort(paths.begin(), paths.end());

	const auto cfg = default_cfg();
	std::vector<cv::Mat> frames_raw, frames_blurred;
	std::vector<double> q_rank;
	for (const auto &p : paths) {
		cv::Mat m = cv::imread(p.string(), cv::IMREAD_UNCHANGED);
		frames_raw.push_back(m);
		frames_blurred.push_back(blurred(m, cfg));
		q_rank.push_back(mpp::rank_score_normalized(m, cfg));
	}
	const auto align = mpp::align_global_from_frames(frames_blurred, q_rank, cfg);
	const auto avg   = mpp::align_average_frame(frames_raw, q_rank, align.shifts, cfg);
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(avg.mean_frame, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	const auto offsets = mpp::shift_frame_offsets(align.shifts, avg.intersection);
	std::vector<double> frame_brightness;
	for (const auto &m : frames_raw)
		frame_brightness.push_back(mpp::rank_average_brightness(m, cfg));
	const int H = frames_raw[0].rows, W = frames_raw[0].cols;
	const auto apq = mpp::ap_compute_frame_qualities(frames_raw, frame_brightness,
	                                                 *aps, offsets, H, W, cfg);
	std::vector<int> sorted_idx(frames_raw.size());
	std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
	std::sort(sorted_idx.begin(), sorted_idx.end(),
	          [&](int a, int b) { return q_rank[a] > q_rank[b]; });
	const auto loop = mpp::stack_frames_loop(frames_raw, frames_blurred,
	                                         avg.mean_frame, *aps, apq,
	                                         offsets, frame_brightness,
	                                         sorted_idx, avg.intersection, cfg);
	const cv::Mat stacked = mpp::stack_merge_alignment_point_buffers(
	    loop.state, loop.border, *aps, cfg);

	const auto dims = read_ints_csv(stacked_dims);
	cr_assert_eq(dims.size(), 2u);
	const int eH = dims[0], eW = dims[1];
	cr_assert_eq(stacked.rows, eH, "height: ours=%d oracle=%d", stacked.rows, eH);
	cr_assert_eq(stacked.cols, eW, "width: ours=%d oracle=%d", stacked.cols, eW);

	cv::Mat oracle(eH, eW, CV_16U);
	std::ifstream fin(stacked_bin, std::ios::binary);
	fin.read(reinterpret_cast<char *>(oracle.data),
	         (std::streamsize) eH * eW * sizeof(uint16_t));
	cr_assert(fin.good(), "stacked_u16.bin read failed");

	int worst = 0, worst_y = -1, worst_x = -1;
	long long total = (long long) eH * eW;
	long long exact = 0;
	for (int y = 0; y < eH; ++y) {
		const uint16_t *o = oracle.ptr<uint16_t>(y);
		const uint16_t *m = stacked.ptr<uint16_t>(y);
		for (int x = 0; x < eW; ++x) {
			const int d = std::abs((int) o[x] - (int) m[x]);
			if (d == 0) ++exact;
			if (d > worst) { worst = d; worst_y = y; worst_x = x; }
		}
	}
	const double exact_pct = 100.0 * (double) exact / (double) total;
	cr_assert_leq(worst, 2,
	             "stacked u16 mismatch at (%d, %d): ours=%d oracle=%d (Δ=%d, exact=%.2f%%)",
	             worst_y, worst_x,
	             (int) stacked.at<uint16_t>(worst_y, worst_x),
	             (int) oracle.at<uint16_t>(worst_y, worst_x), worst, exact_pct);
	cr_assert_gt(exact_pct, 99.0,
	             "only %.2f%% of pixels exactly equal PSS (worst Δ=%d)",
	             exact_pct, worst);

	mpp_ap_free(aps);
}

/* remap_rigid clips when the shifted patch reaches off the frame and
 * records the maximum cropped extent. */
Test(mpp_stack, remap_rigid_clips_and_records_border) {
	cv::Mat frame(20, 20, CV_32F);
	for (int y = 0; y < 20; ++y)
		for (int x = 0; x < 20; ++x)
			frame.at<float>(y, x) = (float) (y * 100 + x);
	cv::Mat buffer(10, 10, CV_32F, cv::Scalar(0));
	mpp::RemapBorder b;
	/* Patch [0,10]×[0,10] with shift (-3, -3): src starts at (-3, -3), so
	 * the top-left 3 rows × 3 cols of the buffer remain zero (the clipped
	 * portion). Border counters record the clipped extents. */
	mpp::stack_remap_rigid(frame, buffer, -3, -3, 0, 10, 0, 10, b);
	cr_assert_eq(b.y_low, 3);
	cr_assert_eq(b.x_low, 3);
	cr_assert_eq(b.y_high, 0);
	cr_assert_eq(b.x_high, 0);
	for (int y = 0; y < 3; ++y)
		for (int x = 0; x < 10; ++x)
			cr_assert_float_eq(buffer.at<float>(y, x), 0.0f, 1e-9);
	/* The non-clipped region holds frame[y-3, x-3] for y ≥ 3, x ≥ 3. */
	cr_assert_float_eq(buffer.at<float>(3, 3), frame.at<float>(0, 0), 1e-6);
	cr_assert_float_eq(buffer.at<float>(9, 9), frame.at<float>(6, 6), 1e-6);
}
