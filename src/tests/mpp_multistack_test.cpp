/*
 * Multi-source derotation stack, end to end (Option B, C5+C6+C7). Two synthetic
 * "sequences" view the same planet: each is rotated by its own central-meridian
 * offsets AND fitted at a DIFFERENT disk position/size in its own frame. Fed
 * through multistack_channel — which relocates + derotates every frame of both
 * sequences into the reference sequence's epoch canvas and stacks them as one —
 * the combined result recovers the sharp base image at the reference position.
 *
 * This exercises the real engine path (global align, average reference, AP
 * qualities, per-AP shifts, warp stack) over a multi-source FrameProvider with
 * the C3 split-disk derotation map, using provider-injected synthetic frames so
 * no sequence I/O is needed.
 */
#include <criterion/criterion.h>
#include <opencv2/imgproc.hpp>

#include <cmath>
#include <vector>

#include "registration/mpp/mpp_derot_geom.h"
#include "registration/mpp/mpp_multistack.hpp"

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp/mpp_config.h"
#include "registration/mpp/mpp_derot_sidecar.h"
}

cominfo com;
fits *gfit = nullptr;

#define DEG (M_PI / 180.0)

namespace {

/* Smooth limb-darkened disk with several SOFT (Gaussian) spots at assorted
 * longitudes — longitudinal structure makes rotation observable, while the
 * smoothness keeps the resampling floor low so misalignment, not ringing,
 * dominates the error metric. */
cv::Mat make_planet(int H, int W, double cx, double cy, double R) {
	const struct { double fx, fy, amp, sig; } spots[] = {
		{-0.45, -0.30, 0.22, 7.0}, { 0.40, -0.15, 0.20, 8.0},
		{-0.10,  0.35, 0.18, 6.0}, { 0.30,  0.40, 0.24, 7.0},
		{ 0.05, -0.50, 0.16, 6.0}, {-0.55,  0.10, 0.20, 7.0},
	};
	cv::Mat f(H, W, CV_16U, cv::Scalar(0));
	for (int y = 0; y < H; ++y)
		for (int x = 0; x < W; ++x) {
			const double r = std::hypot(y - cy, x - cx);
			if (r >= R) continue;
			double v = 0.45 + 0.18 * std::tanh((R - r) / 8.0);  /* limb darkening */
			for (const auto &s : spots) {
				const double dx = x - (cx + s.fx * R), dy = y - (cy + s.fy * R);
				v += s.amp * std::exp(-(dx * dx + dy * dy) / (2.0 * s.sig * s.sig));
			}
			f.at<uint16_t>(y, x) = (uint16_t) std::clamp(v * 65535.0, 0.0, 65535.0);
		}
	return f;
}

double disk_mean_abs_diff(const cv::Mat &a, const cv::Mat &b,
                          double cx, double cy, double R) {
	double sum = 0.0; int n = 0;
	for (int y = 0; y < a.rows && y < b.rows; ++y)
		for (int x = 0; x < a.cols && x < b.cols; ++x)
			if (std::hypot(y - cy, x - cx) < R * 0.75) {
				sum += std::abs((double) a.at<uint16_t>(y, x)
				              - (double) b.at<uint16_t>(y, x));
				++n;
			}
	return n ? sum / n : 0.0;
}

/* Render a source sequence's frame: place the base globe at this source's disk
 * and rotate it to `cm`. This is the INVERSE of the orchestrator's relocate
 * (epoch/frame swapped), so the orchestrator un-does it back to base. */
cv::Mat render_view(const cv::Mat &base, int W, int H,
                    const derot_diskfit_t &src_disk, const derot_diskfit_t &ref_disk,
                    double B, double cm, double CM0) {
	derot_geom_t epoch{B * DEG, cm * DEG, 0.0};   /* this frame's orientation */
	derot_geom_t frame{B * DEG, CM0 * DEG, 0.0};  /* the base (epoch) orientation */
	cv::Mat mx, my, valid, mu, view;
	derot_build_map_ms(W, H, src_disk, ref_disk, epoch, frame, mx, my, valid, mu);
	cv::remap(base, view, mx, my, cv::INTER_LANCZOS4, cv::BORDER_CONSTANT, cv::Scalar(0));
	view.setTo(0, valid == 0);
	return view;
}

mpp_derot_t *make_plan(int N, const derot_diskfit_t &disk, int W, int H,
                       double B, double CM0, double dCM) {
	mpp_derot_t *d = mpp_derot_alloc(N);
	d->num_frames = N; d->frame_rows = H; d->frame_cols = W;
	d->cx = disk.cx; d->cy = disk.cy; d->r_eq = disk.r_eq;
	d->flattening = 0.0; d->parity = 1.0; d->pole_angle_epoch = 0.0;
	d->epoch_sub_obs_lat = B; d->epoch_cm = CM0; d->epoch_pole_pa = 0.0;
	for (int f = 0; f < N; ++f) {
		d->sub_obs_lat[f] = B;
		d->cm[f] = CM0 + (f - (N - 1) / 2.0) * dCM;
		d->pole_pa[f] = 0.0;
	}
	return d;
}

}  // namespace

Test(mpp_multistack, combines_two_sources_into_reference_canvas) {
	const int W = 260, H = 220;
	const double B = 3.0, CM0 = 40.0, dCM = 5.0;
	const int N = 7;   /* central-meridian span +/-15 deg */

	/* reference disk (source 0) and a strongly displaced/smaller disk (src 1) */
	const derot_diskfit_t ref_disk{118.0, 100.0, 68.0, 0.0, 1.0};
	const derot_diskfit_t s1_disk { 156.0, 124.0, 60.0, 0.0, 1.0};

	const cv::Mat base = make_planet(H, W, ref_disk.cx, ref_disk.cy, ref_disk.r_eq);

	/* synthesise both sources' frames */
	std::vector<cv::Mat> f0(N), f1(N);
	for (int f = 0; f < N; ++f) {
		const double cm = CM0 + (f - (N - 1) / 2.0) * dCM;
		f0[f] = render_view(base, W, H, ref_disk, ref_disk, B, cm, CM0);
		f1[f] = render_view(base, W, H, s1_disk,  ref_disk, B, cm, CM0);
	}

	/* Naive un-derotated, un-relocated average of the whole union — the smeared
	 * baseline the derotation+relocation pipeline must beat. */
	cv::Mat acc(H, W, CV_32F, cv::Scalar(0));
	for (int f = 0; f < N; ++f) {
		cv::Mat a, b;
		f0[f].convertTo(a, CV_32F); f1[f].convertTo(b, CV_32F);
		acc += a; acc += b;
	}
	acc /= (2.0 * N);
	cv::Mat naive; acc.convertTo(naive, CV_16U);
	const double naive_err = disk_mean_abs_diff(naive, base,
	                                            ref_disk.cx, ref_disk.cy, ref_disk.r_eq);

	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	cfg.frames_normalization = false;
	cfg.alignment_points_frame_percent = 100;
	cfg.stack_method = MPP_STACK_WARP;
	cfg.drizzle_scale = 1.0;

	mpp_derot_t *d0 = make_plan(N, ref_disk, W, H, B, CM0, dCM);
	mpp_derot_t *d1 = make_plan(N, s1_disk,  W, H, B, CM0, dCM);

	std::vector<mpp::MsSource> srcs(2);
	srcs[0].derot = d0; srcs[0].num_frames = N;
	srcs[0].analysis_read = [&](int l) { return f0[l]; };
	srcs[0].full_read     = [&](int l) { return f0[l]; };
	srcs[1].derot = d1; srcs[1].num_frames = N;
	srcs[1].analysis_read = [&](int l) { return f1[l]; };
	srcs[1].full_read     = [&](int l) { return f1[l]; };

	const auto res = mpp::multistack_channel(srcs, /*out_derot=*/d0, /*num_layers=*/1, cfg, 1);
	cr_assert(!res.error, "multistack reported error");
	cr_assert(!res.oom && !res.cancelled, "multistack oom/cancelled");
	cr_assert(!res.image.empty(), "multistack produced no image");
	cr_assert_eq(res.num_frames, 2 * N);
	cr_assert_gt(res.num_aps, 4);

	/* the recovered stack, placed via the intersection origin, matches base */
	const double ox = res.intersection[2];   /* x_low */
	const double oy = res.intersection[0];   /* y_low */
	const double err = disk_mean_abs_diff(res.image, base,
	                                      ref_disk.cx - ox, ref_disk.cy - oy,
	                                      ref_disk.r_eq);
	cr_log_info("multistack: intersection=(%d,%d,%d,%d) err=%.0f naive_err=%.0f",
	            res.intersection[0], res.intersection[1], res.intersection[2],
	            res.intersection[3], err, naive_err);
	/* the derotated+relocated combine must clearly beat the smeared naive
	 * average of the same frames, and recover the base well in absolute terms */
	cr_assert_lt(err, naive_err * 0.5,
	             "combine (%.0f) should beat the naive smeared average (%.0f)",
	             err, naive_err);
	cr_assert_lt(err, 1100.0,
	             "combined multi-source stack should recover the base (err %.0f LSB)",
	             err);

	mpp_derot_free(d0);
	mpp_derot_free(d1);
}

/* Regression for the warp-stack source-map placement: when the frames drift
 * (so the global alignment crops the intersection to a non-zero origin), the
 * derotation map's source side must carry the same intersection origin/scale as
 * the output side. With the source placed at (0,0,1) instead, every frame is
 * sampled off by the intersection origin and the recovered planet smears. No
 * rotation here — pure per-frame translation — isolates the placement: the only
 * way recovery can succeed is if that origin cancels correctly. */
Test(mpp_multistack, drift_cropped_intersection_recovers) {
	const int W = 240, H = 220;
	const int N = 9;
	const double cx = 120, cy = 108, R = 64, B = 4.0, CM = 30.0;
	const cv::Mat base = make_planet(H, W, cx, cy, R);

	/* Frame 0 is the untranslated base (sharpest, so the aligner picks it as the
	 * reference); the rest are translated by NEGATIVE integer steps. Aligning
	 * them back to frame 0 needs positive shifts, so the intersection crops to a
	 * large positive origin — exactly the case that broke when the derotation
	 * map's source side dropped that origin. Integer steps keep the (integer)
	 * global/per-AP search exact. */
	std::vector<cv::Mat> frames(N);
	frames[0] = base.clone();
	for (int f = 1; f < N; ++f) {
		const int tx = -3 * f, ty = -2 * f;     /* down to (-24, -16) at f=8 */
		cv::Mat M = (cv::Mat_<double>(2, 3) << 1, 0, tx, 0, 1, ty);
		cv::warpAffine(base, frames[f], M, base.size(), cv::INTER_LANCZOS4,
		               cv::BORDER_CONSTANT, cv::Scalar(0));
	}

	mpp_config_t cfg{};
	mpp_config_defaults(&cfg);
	cfg.frames_normalization = false;
	cfg.alignment_points_frame_percent = 100;
	cfg.align_frames_mode = MPP_ALIGN_SURFACE;   /* tracks integer translation */
	cfg.stack_method = MPP_STACK_WARP;
	cfg.drizzle_scale = 1.0;

	mpp_derot_t *d = mpp_derot_alloc(N);
	d->num_frames = N; d->frame_rows = H; d->frame_cols = W;
	d->cx = cx; d->cy = cy; d->r_eq = R; d->flattening = 0.0; d->parity = 1.0;
	d->pole_angle_epoch = 0.0;
	d->epoch_sub_obs_lat = B; d->epoch_cm = CM; d->epoch_pole_pa = 0.0;
	for (int f = 0; f < N; ++f) {                 /* no rotation: all == epoch */
		d->sub_obs_lat[f] = B; d->cm[f] = CM; d->pole_pa[f] = 0.0;
	}

	std::vector<mpp::MsSource> srcs(1);
	srcs[0].derot = d; srcs[0].num_frames = N;
	srcs[0].analysis_read = [&](int l) { return frames[l]; };
	srcs[0].full_read     = [&](int l) { return frames[l]; };

	const auto res = mpp::multistack_channel(srcs, d, 1, cfg, 1);
	cr_assert(!res.error && !res.oom && !res.cancelled && !res.image.empty());

	const int ox = res.intersection[2], oy = res.intersection[0];
	cr_log_info("drift intersection=(%d,%d,%d,%d) img=%dx%d",
	            res.intersection[0], res.intersection[1], res.intersection[2],
	            res.intersection[3], res.image.cols, res.image.rows);
	/* the intersection must actually be cropped, else the test proves nothing */
	cr_assert(ox >= 1 && oy >= 1 && res.image.cols < W && res.image.rows < H,
	          "expected a cropped intersection (got x_low=%d y_low=%d)", ox, oy);
	/* crop-aware compare: res.image[y][x] corresponds to base[y+oy][x+ox] */
	double sum = 0.0; int n = 0;
	for (int y = 0; y < res.image.rows; ++y)
		for (int x = 0; x < res.image.cols; ++x) {
			const int by = y + oy, bx = x + ox;
			if (by < 0 || by >= base.rows || bx < 0 || bx >= base.cols) continue;
			if (std::hypot(y - (cy - oy), x - (cx - ox)) < R * 0.7) {
				sum += std::abs((double) res.image.at<uint16_t>(y, x)
				              - (double) base.at<uint16_t>(by, bx));
				++n;
			}
		}
	const double err = n ? sum / n : 1e9;
	cr_log_info("drift test: intersection=(%d,%d,%d,%d) img=%dx%d err=%.0f",
	            res.intersection[0], res.intersection[1],
	            res.intersection[2], res.intersection[3],
	            res.image.cols, res.image.rows, err);
	/* With the source placement + offset convention correct this recovers the
	 * base (residual is the Lanczos / integer-alignment floor of an aggressive
	 * +/-24 px synthetic drift). The pre-fix code (source at (0,0,1), raw shift
	 * offset) mis-samples by the ~21 px origin and roughly doubles the error. */
	cr_assert_lt(err, 2000.0,
	             "drifted frames should recover the base after global align "
	             "(err %.0f LSB)", err);

	mpp_derot_free(d);
}
