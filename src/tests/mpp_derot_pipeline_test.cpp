/*
 * Derotation pipeline integration: synthesise a set of views of a planet
 * rotated by known central-meridian offsets, then stack them through the warp
 * engine driven by the real .derot -> geometry bridge. A correct pipeline
 * un-rotates every frame to the epoch and recovers the (sharp) base image; a
 * non-derotated stack of the same frames is visibly smeared. This exercises
 * C2 (geometry) + C4 (warp engine) + C5 (bridge) together.
 */
#include <criterion/criterion.h>
#include <opencv2/imgproc.hpp>

#include <cmath>
#include <numeric>
#include <vector>

#include "registration/mpp/mpp_derot_geom.h"
#include "registration/mpp/mpp_derot.h"   /* C++ (cv::Mat) — outside extern "C" */
#include "registration/mpp/mpp_ap_priv.hpp"
#include "registration/mpp/mpp_rank_priv.hpp"
#include "registration/mpp/mpp_shift_priv.hpp"
#include "registration/mpp/mpp_stack_priv.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp/mpp_ap.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_derot_sidecar.h"
}

cominfo com;
fits *gfit = nullptr;

#define DEG (M_PI / 180.0)

namespace {

cv::Mat make_planet(int H, int W, double cx, double cy, double R) {
	cv::Mat f(H, W, CV_16U, cv::Scalar(0));
	cv::RNG rng(7);
	for (int y = 0; y < H; ++y)
		for (int x = 0; x < W; ++x) {
			const double r = std::hypot(y - cy, x - cx);
			if (r < R) {
				double v = 0.45 + 0.2 * std::tanh((R - r) / 5.0);
				f.at<uint16_t>(y, x) = (uint16_t) std::clamp(v * 65535.0, 0.0, 65535.0);
			}
		}
	/* belts + spots so rotation is visible */
	for (int k = 0; k < 60; ++k) {
		const int sx = (int) rng.uniform(cx - R * 0.7, cx + R * 0.7);
		const int sy = (int) rng.uniform(cy - R * 0.7, cy + R * 0.7);
		if (std::hypot(sy - cy, sx - cx) < R * 0.85)
			cv::circle(f, {sx, sy}, 2, cv::Scalar(62000), -1);
	}
	return f;
}

double disk_mean_abs_diff(const cv::Mat &a, const cv::Mat &b,
                          double cx, double cy, double R) {
	double sum = 0.0; int n = 0;
	for (int y = 0; y < a.rows; ++y)
		for (int x = 0; x < a.cols; ++x)
			if (std::hypot(y - cy, x - cx) < R * 0.8) {
				sum += std::abs((double) a.at<uint16_t>(y, x) - (double) b.at<uint16_t>(y, x));
				++n;
			}
	return n ? sum / n : 0.0;
}

}  // namespace

Test(mpp_derot_pipeline, recovers_rotating_planet) {
	const int H = 200, W = 240;
	const double cx = 120, cy = 100, R = 80;
	const cv::Mat base = make_planet(H, W, cx, cy, R);

	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	cfg.frames_normalization = false;
	cfg.alignment_points_frame_percent = 100;

	/* geometry: epoch in the middle, frames rotated +/- in CM (System angle) */
	const int N = 5;
	const double B = 3.0, P = 0.0, CM0 = 40.0, dCM = 4.0;   /* +/-8 deg total span */
	derot_diskfit_t disk{cx, cy, R, 0.0, 1.0};
	derot_geom_t epoch{B * DEG, CM0 * DEG, 0.0};

	std::vector<cv::Mat> frames;
	for (int f = 0; f < N; ++f) {
		const double cm = CM0 + (f - (N - 1) / 2.0) * dCM;
		derot_geom_t fr{B * DEG, cm * DEG, 0.0};
		/* frame view = base resampled through the frame->epoch map, so that the
		 * epoch->frame map used by the stack inverts it back to base */
		cv::Mat mx, my, valid, mu;
		derot_build_map(W, H, disk, /*epoch=*/fr, /*frame=*/epoch, mx, my, valid, mu);
		cv::Mat view;
		cv::remap(base, view, mx, my, cv::INTER_LANCZOS4, cv::BORDER_REPLICATE);
		view.setTo(0, valid == 0);
		frames.push_back(view);
	}

	/* AP grid + per-AP qualities on the synthesised frames, zero global offset */
	const cv::Mat mean_blurred = mpp::blur_mean_frame_for_ap(base, cfg);
	mpp_aps_t *aps = mpp::ap_create_grid(mean_blurred, cfg);
	cr_assert_not_null(aps);
	cr_assert_gt(aps->count, 4);
	std::vector<mpp::FrameOffset> offsets(N);
	for (auto &o : offsets) { o.dy = 0.0; o.dx = 0.0; }
	std::vector<double> brightness;
	for (auto &f : frames) brightness.push_back(mpp::rank_average_brightness(f, cfg));
	const auto apq = mpp::ap_compute_frame_qualities(
	    frames, brightness, *aps, offsets, H, W, cfg);
	std::vector<int> sorted(N);
	std::iota(sorted.begin(), sorted.end(), 0);
	const cv::Vec4i isect(0, H, 0, W);
	mpp::FrameProvider provider = [&](int i) { return frames[i]; };

	/* .derot plan feeding the bridge provider */
	mpp_derot_t *d = mpp_derot_alloc(N);
	d->cx = cx; d->cy = cy; d->r_eq = R; d->flattening = 0.0; d->parity = 1.0;
	d->pole_angle_epoch = 0.0;
	d->epoch_sub_obs_lat = B; d->epoch_cm = CM0; d->epoch_pole_pa = P;
	for (int f = 0; f < N; ++f) {
		d->sub_obs_lat[f] = B;
		d->cm[f] = CM0 + (f - (N - 1) / 2.0) * dCM;
		d->pole_pa[f] = P;
	}
	mpp::DerotMapProvider bridge =
	    [d, &isect](int f, int DY, int DX, cv::Mat &mx, cv::Mat &my, cv::Mat &mu) {
		mpp_derot_frame_map(d, f, DX, DY, (double) isect[2], (double) isect[0],
		                    1.0, mx, my, mu);
		return true;
	};

	auto stack = [&](const mpp::DerotMapProvider *dp) {
		return mpp::stack_warp_apply_streamed(
		    provider, N, 1, *aps, apq, /*shifts=*/nullptr, offsets, brightness,
		    sorted, isect, cfg, nullptr, 1, false, 0, dp);
	};

	const auto derot = stack(&bridge);
	const auto plain = stack(nullptr);
	cr_assert(!derot.image.empty() && !plain.image.empty());

	const double err_derot = disk_mean_abs_diff(derot.image, base, cx, cy, R);
	const double err_plain = disk_mean_abs_diff(plain.image, base, cx, cy, R);

	/* The derotated stack must recover the base far better than the smeared,
	 * non-derotated stack of the same rotated frames. */
	cr_assert_lt(err_derot, err_plain * 0.5,
	             "derot err %.0f should beat plain err %.0f", err_derot, err_plain);
	cr_assert_lt(err_derot, 2500.0,
	             "derotated stack should recover the base (err %.0f LSB)", err_derot);

	mpp_derot_free(d);
	mpp_ap_free(aps);
}
