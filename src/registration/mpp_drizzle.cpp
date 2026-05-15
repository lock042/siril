/*
 * mpp_drizzle — Phase 5b STScI drizzle backend for MPP stacking.
 *
 * Slice 5b.1 (the pixmap builder): encodes the full MPP warp as a per-
 * pixel output-coordinate map suitable for handing to dobox().
 *
 * Slice 5b.2 (this file's STScI stack path): drives dobox() one frame
 * at a time, accumulating into a shared float output canvas + counts
 * buffer, then normalises and casts to uint16.
 *
 * Slice 5b.3 (Bayer-drizzle) lands in a follow-up.
 *
 * NB: include ordering. driz_portability.h does `#define private` for C
 * compat, which wrecks subsequent C++ standard-library includes (their
 * class definitions use `private:` as a keyword). All C++ standard
 * headers and C++ helpers MUST be included before the drizzle headers.
 * The drizzle block is included last and bookended with `#undef private`.
 */

/* === C++ std + Siril's C++-side headers FIRST === */
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_drizzle.h"
#include "registration/mpp_shift.h"

#include "core/proto.h"
#include "core/gui_iface.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "io/ser.h"
}

#include "registration/mpp_drizzle_priv.hpp"  /* mpp::stack_apply_stsci */
#include "registration/mpp_stack_priv.hpp"    /* mpp::stack_one_dim_weight */

/* read_full_frame is in mpp.cpp's mpp:: namespace (moved out of the
 * inner anonymous namespace in slice 5b.2). Forward-declared here so we
 * can read frames without exposing it through a public header. */
namespace mpp {
extern cv::Mat read_full_frame(sequence *seq, int idx);
}

/* === Drizzle headers LAST (they redefine `private`) === */
/* cdrizzlemap.h is deliberately NOT included — it pulls in wcslib.h
 * which has a conflicting `int64` typedef vs OpenCV's. cdrizzlebox.h
 * (with dobox + driz_param_t + e_kernel_t via cdrizzleutil.h) is all
 * we need for the stack path. */
extern "C" {
#include "drizzle/cdrizzleutil.h"
#include "drizzle/cdrizzlebox.h"
}
#undef private


/* ============================================================
 * Pixmap (slice 5b.1) — unchanged in this commit.
 * ============================================================ */

extern "C" mpp_status_t mpp_imgmap_alloc(imgmap_t *out, int rx, int ry) {
	if (!out || rx <= 0 || ry <= 0) return MPP_EINVAL;
	const size_t n = (size_t) rx * (size_t) ry;
	out->xmap = (float *) std::malloc(n * sizeof(float));
	if (!out->xmap) { out->ymap = nullptr; out->rx = 0; out->ry = 0; return MPP_ENOMEM; }
	out->ymap = (float *) std::malloc(n * sizeof(float));
	if (!out->ymap) {
		std::free(out->xmap);
		out->xmap = nullptr;
		out->rx = 0; out->ry = 0;
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

	const double iy_lo = (double) run->intersection[0];
	const double ix_lo = (double) run->intersection[2];

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

	const size_t frame_base = (size_t) frame_idx * (size_t) M;
	const double *ap_shifts = &run->shifts->shifts[2 * frame_base];

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


/* ============================================================
 * STScI stack path (slice 5b.2).
 * ============================================================ */

namespace {

enum e_kernel_t kernel_from_cfg(int mpp_kernel) {
	switch (mpp_kernel) {
		case MPP_KERNEL_GAUSSIAN: return kernel_gaussian;
		case MPP_KERNEL_POINT:    return kernel_point;
		case MPP_KERNEL_TURBO:    return kernel_turbo;
		case MPP_KERNEL_LANCZOS2: return kernel_lanczos2;
		case MPP_KERNEL_LANCZOS3: return kernel_lanczos3;
		case MPP_KERNEL_SQUARE:
		default:                  return kernel_square;
	}
}

/* Allocate a float-typed `fits` with the given shape, layout planar
 * (R then G then B for multi-channel) and pdata pointers set up for
 * dobox's get_pixel accessors. Returns false on allocation failure. */
bool alloc_float_fits(fits *out, int rx, int ry, int channels) {
	std::memset(out, 0, sizeof(fits));
	out->rx = rx;
	out->ry = ry;
	out->naxes[0] = rx;
	out->naxes[1] = ry;
	out->naxes[2] = channels;
	out->naxis = (channels == 1) ? 2 : 3;
	out->bitpix = FLOAT_IMG;
	out->orig_bitpix = FLOAT_IMG;
	out->type = DATA_FLOAT;
	const size_t plane = (size_t) rx * (size_t) ry;
	out->fdata = (float *) std::calloc(plane * (size_t) channels, sizeof(float));
	if (!out->fdata) return false;
	out->fpdata[0] = out->fdata;
	out->fpdata[1] = (channels >= 2) ? out->fdata + plane     : nullptr;
	out->fpdata[2] = (channels >= 3) ? out->fdata + 2 * plane : nullptr;
	return true;
}

/* Build a planar single-frame `fits` from a cv::Mat. The Mat may be
 * single-channel or 3-channel (interleaved BGR per OpenCV convention);
 * the fits is always planar. Applies the brightness-equalisation scale
 * during conversion. Caller owns `out->fdata` (clearfits to free). */
bool cv_to_planar_fits(const cv::Mat &raw, double scale, fits *out) {
	const int rx = raw.cols, ry = raw.rows, channels = raw.channels();
	if (!alloc_float_fits(out, rx, ry, channels)) return false;
	const size_t plane = (size_t) rx * (size_t) ry;
	if (channels == 1) {
		cv::Mat as_float;
		raw.convertTo(as_float, CV_32F, scale);
		std::memcpy(out->fdata, as_float.data, plane * sizeof(float));
	} else {
		/* Split into R/G/B planes. OpenCV stores BGR by default, but
		 * mpp::read_full_frame returns frames in Siril's R/G/B order
		 * (see mpp.cpp read_full_frame's cv::merge of pdata[0..2]). */
		cv::Mat as_float;
		raw.convertTo(as_float, CV_32FC3, scale);
		std::vector<cv::Mat> planes;
		cv::split(as_float, planes);
		for (int c = 0; c < channels; ++c)
			std::memcpy(out->fdata + (size_t) c * plane,
			            planes[c].data, plane * sizeof(float));
	}
	return true;
}

}  // namespace

/* Inner C++ entry point. Tests drive this directly with vector<cv::Mat>
 * so they can run on PNG-extracted frames without constructing a Siril
 * sequence. The extern "C" mpp_stack_apply_stsci wrapper just reads
 * frames via mpp::read_full_frame and forwards. */
mpp_status_t mpp::stack_apply_stsci(const std::vector<cv::Mat> &frames_raw,
                                    const std::vector<int> &included,
                                    const std::vector<double> &frame_brightness,
                                    const mpp_run_t *run,
                                    const mpp_config_t *cfg,
                                    fits *out) {
	if (!cfg || !run || !out || !run->aps || !run->shifts) return MPP_EINVAL;
	const int N = (int) frames_raw.size();
	if (N <= 0 || (int) included.size() != N
	 || (int) frame_brightness.size() != N) return MPP_EINVAL;
	if (N != run->num_frames) return MPP_EINVAL;

	const int frame_rx = run->frame_cols;
	const int frame_ry = run->frame_rows;
	if (frame_rx <= 0 || frame_ry <= 0) return MPP_EINVAL;

	const int intersection_rx = run->intersection[3] - run->intersection[2];
	const int intersection_ry = run->intersection[1] - run->intersection[0];
	if (intersection_rx <= 0 || intersection_ry <= 0) return MPP_EINVAL;

	const int factor   = cfg->drizzle_factor < 1 ? 1 : cfg->drizzle_factor;
	const int out_rx   = intersection_rx * factor;
	const int out_ry   = intersection_ry * factor;
	const int channels = run->num_layers > 0 ? run->num_layers : 1;

	fits output_data{};
	fits output_counts{};
	if (!alloc_float_fits(&output_data, out_rx, out_ry, channels)) return MPP_ENOMEM;
	if (!alloc_float_fits(&output_counts, out_rx, out_ry, channels)) {
		clearfits(&output_data);
		return MPP_ENOMEM;
	}

	imgmap_t pixmap{};
	if (mpp_imgmap_alloc(&pixmap, frame_rx, frame_ry) != MPP_OK) {
		clearfits(&output_data);
		clearfits(&output_counts);
		return MPP_ENOMEM;
	}

	std::vector<double> bright_inc;
	bright_inc.reserve(N);
	for (int i = 0; i < N; ++i)
		if (included[i]) bright_inc.push_back(frame_brightness[i]);
	if (bright_inc.empty()) {
		mpp_imgmap_free(&pixmap);
		clearfits(&output_data);
		clearfits(&output_counts);
		return MPP_ENODATA;
	}
	std::nth_element(bright_inc.begin(),
	                 bright_inc.begin() + bright_inc.size() / 2,
	                 bright_inc.end());
	const double median_brightness = bright_inc[bright_inc.size() / 2];

	struct driz_args_t driz_args;
	std::memset(&driz_args, 0, sizeof(driz_args));
	driz_args.kernel         = kernel_from_cfg(cfg->drizzle_kernel);
	driz_args.scale          = (float) factor;
	driz_args.pixel_fraction = (float) cfg->drizzle_pixfrac;
	driz_args.weight_scale   = 1.0f;
	driz_args.is_bayer       = FALSE;
	driz_args.cfadim         = 0;
	driz_args.use_flats      = FALSE;
	driz_args.flat           = nullptr;

	struct driz_error_t error;
	std::memset(&error, 0, sizeof(error));

	const int n_included = (int) bright_inc.size();
	int cur_nb = 0;
	gui_iface.set_progress(PROGRESS_RESET, _("Stack (STScI): drizzling frames"));

	for (int f = 0; f < N; ++f) {
		if (!included[f] || frames_raw[f].empty()) continue;

		const double scale = median_brightness
		                   / (frame_brightness[f] + 1e-7);

		fits frame_fit{};
		if (!cv_to_planar_fits(frames_raw[f], scale, &frame_fit)) {
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_ENOMEM;
		}

		if (mpp_pixmap_build(run, f, (double) factor, &pixmap) != MPP_OK) {
			clearfits(&frame_fit);
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EINVAL;
		}

		struct driz_param_t p;
		std::memset(&p, 0, sizeof(p));
		driz_param_init(&p);
		p.driz           = &driz_args;
		p.kernel         = driz_args.kernel;
		/* p.scale = 1.0 (not driz_args.scale). dobox uses p.scale only
		 * to compute a `scale²` flux-conservation multiplier on each
		 * input pixel value — at drizzle=2 every input value would be
		 * inflated 4× to preserve total energy across the spread of an
		 * input pixel's footprint over multiple output pixels. PSS's
		 * bicubic resampler doesn't do that (it just resamples each
		 * frame per-pixel), and the per-AP top-N weighted mean we feed
		 * back to dobox via the pixmap already keeps the values in
		 * input units. The 2x output canvas geometry is fully carried
		 * by the pixmap's drizzle_scale multiplier, so setting p.scale
		 * to 1.0 gives values that compare cleanly to PSS's drizzle=2
		 * reference. */
		p.scale          = 1.0f;
		p.pixel_fraction = driz_args.pixel_fraction;
		p.weight_scale   = 1.0f;
		p.in_units       = unit_counts;
		p.out_units      = unit_counts;
		p.exposure_time  = 1.0f;
		p.data           = &frame_fit;
		p.weights        = nullptr;
		p.pixmap         = &pixmap;
		p.output_data    = &output_data;
		p.output_counts  = &output_counts;
		p.xmin           = 0;
		p.ymin           = 0;
		p.xmax           = frame_rx - 1;
		p.ymax           = frame_ry - 1;
		p.threads        = 1;
		p.error          = &error;

		if (dobox(&p)) {
			siril_log_color_message(
			    _("Stack (STScI): dobox failed on frame %d: %s\n"),
			    "red", f, error.last_message[0] ? error.last_message : "(no detail)");
			clearfits(&frame_fit);
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EIO;
		}
		clearfits(&frame_fit);

		++cur_nb;
		gui_iface.set_progress((double) cur_nb / (double) n_included, NULL);
	}

	mpp_imgmap_free(&pixmap);

	/* update_data inside dobox accumulates a running weighted mean — by
	 * the time the per-frame loop exits, output_data already holds the
	 * weighted average of every contribution. Empty cells (no frame
	 * covered them) stay at the calloc'd zero, which is what we want.
	 * Do NOT divide by output_counts: that would mean / total-weight =
	 * mean × ~1/N which collapses everything to near-zero. */
	(void) output_counts;   /* kept around in case BAYER path needs per-channel counts */
	clearfits(&output_counts);

	clearfits(out);
	out->rx = out_rx;
	out->ry = out_ry;
	out->naxes[0] = out_rx;
	out->naxes[1] = out_ry;
	out->naxes[2] = channels;
	out->naxis = (channels == 1) ? 2 : 3;
	out->bitpix = USHORT_IMG;
	out->orig_bitpix = USHORT_IMG;
	out->type = DATA_USHORT;
	const size_t plane = (size_t) out_rx * (size_t) out_ry;
	out->data = (WORD *) std::calloc(plane * (size_t) channels, sizeof(WORD));
	if (!out->data) {
		clearfits(&output_data);
		return MPP_ENOMEM;
	}
	out->pdata[0] = out->data;
	out->pdata[1] = (channels >= 2) ? out->data + plane     : nullptr;
	out->pdata[2] = (channels >= 3) ? out->data + 2 * plane : nullptr;

	/* Scale to 16-bit uint range. For 8-bit input, output_data values
	 * are still in [0, 255] (the raw input range we fed dobox) so we
	 * multiply by 256 here — mirrors PSS's "mean_frame in 16-bit
	 * equivalent regardless of bitdepth" convention. For 16-bit input
	 * the scale is 1.0 and values pass through unchanged.
	 * NB: mpp_cfg_threshold_scale is the INVERSE direction (scales
	 * thresholds down to match input units), so we compute the output
	 * scale inline rather than reusing it. */
	const double cast_scale = (cfg->bitdepth == 8) ? 256.0 : 1.0;
	for (size_t k = 0; k < plane * (size_t) channels; ++k) {
		const double v = (double) output_data.fdata[k] * cast_scale;
		out->data[k] = (WORD) (v < 0.0 ? 0
		                     : (v > 65535.0 ? 65535
		                                    : (int) (v + 0.5)));
	}
	clearfits(&output_data);
	return MPP_OK;
}

/* ============================================================
 * Bayer drizzle (slice 5b.3).
 *
 * Same as the STScI path but:
 *   - input frames are single-channel raw Bayer (no debayer);
 *   - driz_args.is_bayer = TRUE with cfa[] + cfadim populated;
 *   - output is always 3-channel planar (R, G, B), even if input is
 *     single-channel.
 *
 * dobox's do_kernel_square (cdrizzlebox.c:998) reads the CFA channel
 * from the *input* pixel position, so the per-pixel routing is
 * intrinsic to (i, j) and works fine with our non-affine pixmap
 * (every input pixel still has a fixed colour identity regardless of
 * where the pixmap sends it).
 * ============================================================ */

mpp_status_t mpp::stack_apply_bayer(const std::vector<cv::Mat> &frames_raw_bayer,
                                    const std::vector<int> &included,
                                    const std::vector<double> &frame_brightness,
                                    const mpp_run_t *run,
                                    const mpp_config_t *cfg,
                                    const unsigned char *cfa,
                                    int cfadim,
                                    fits *out) {
	if (!cfg || !run || !out || !run->aps || !run->shifts || !cfa) return MPP_EINVAL;
	if (cfadim != 2) return MPP_EINVAL;     /* 2x2 Bayer only for now */
	const int N = (int) frames_raw_bayer.size();
	if (N <= 0 || (int) included.size() != N
	 || (int) frame_brightness.size() != N) return MPP_EINVAL;
	if (N != run->num_frames) return MPP_EINVAL;
	for (const auto &f : frames_raw_bayer) {
		if (!f.empty() && f.channels() != 1) return MPP_EINVAL;
	}

	const int frame_rx = run->frame_cols;
	const int frame_ry = run->frame_rows;
	if (frame_rx <= 0 || frame_ry <= 0) return MPP_EINVAL;

	const int intersection_rx = run->intersection[3] - run->intersection[2];
	const int intersection_ry = run->intersection[1] - run->intersection[0];
	if (intersection_rx <= 0 || intersection_ry <= 0) return MPP_EINVAL;

	const int factor   = cfg->drizzle_factor < 1 ? 1 : cfg->drizzle_factor;
	const int out_rx   = intersection_rx * factor;
	const int out_ry   = intersection_ry * factor;
	const int channels = 3;                  /* Bayer output is always RGB */

	fits output_data{};
	fits output_counts{};
	if (!alloc_float_fits(&output_data, out_rx, out_ry, channels)) return MPP_ENOMEM;
	if (!alloc_float_fits(&output_counts, out_rx, out_ry, channels)) {
		clearfits(&output_data);
		return MPP_ENOMEM;
	}

	imgmap_t pixmap{};
	if (mpp_imgmap_alloc(&pixmap, frame_rx, frame_ry) != MPP_OK) {
		clearfits(&output_data);
		clearfits(&output_counts);
		return MPP_ENOMEM;
	}

	std::vector<double> bright_inc;
	bright_inc.reserve(N);
	for (int i = 0; i < N; ++i)
		if (included[i]) bright_inc.push_back(frame_brightness[i]);
	if (bright_inc.empty()) {
		mpp_imgmap_free(&pixmap);
		clearfits(&output_data);
		clearfits(&output_counts);
		return MPP_ENODATA;
	}
	std::nth_element(bright_inc.begin(),
	                 bright_inc.begin() + bright_inc.size() / 2,
	                 bright_inc.end());
	const double median_brightness = bright_inc[bright_inc.size() / 2];

	struct driz_args_t driz_args;
	std::memset(&driz_args, 0, sizeof(driz_args));
	driz_args.kernel         = kernel_from_cfg(cfg->drizzle_kernel);
	driz_args.scale          = (float) factor;
	driz_args.pixel_fraction = (float) cfg->drizzle_pixfrac;
	driz_args.weight_scale   = 1.0f;
	driz_args.is_bayer       = TRUE;
	std::memcpy(driz_args.cfa, cfa, 4);
	driz_args.cfadim         = (size_t) cfadim;
	driz_args.use_flats      = FALSE;
	driz_args.flat           = nullptr;

	struct driz_error_t error;
	std::memset(&error, 0, sizeof(error));

	const int n_included = (int) bright_inc.size();
	int cur_nb = 0;
	gui_iface.set_progress(PROGRESS_RESET, _("Stack (Bayer drizzle): processing frames"));

	for (int f = 0; f < N; ++f) {
		if (!included[f] || frames_raw_bayer[f].empty()) continue;

		const double scale = median_brightness
		                   / (frame_brightness[f] + 1e-7);

		fits frame_fit{};
		if (!cv_to_planar_fits(frames_raw_bayer[f], scale, &frame_fit)) {
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_ENOMEM;
		}

		if (mpp_pixmap_build(run, f, (double) factor, &pixmap) != MPP_OK) {
			clearfits(&frame_fit);
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EINVAL;
		}

		struct driz_param_t p;
		std::memset(&p, 0, sizeof(p));
		driz_param_init(&p);
		p.driz           = &driz_args;
		p.kernel         = driz_args.kernel;
		p.scale          = 1.0f;            /* same as STScI path — see comment there */
		p.pixel_fraction = driz_args.pixel_fraction;
		p.weight_scale   = 1.0f;
		p.in_units       = unit_counts;
		p.out_units      = unit_counts;
		p.exposure_time  = 1.0f;
		p.data           = &frame_fit;
		p.weights        = nullptr;
		p.pixmap         = &pixmap;
		p.output_data    = &output_data;
		p.output_counts  = &output_counts;
		std::memcpy(p.cfa, cfa, 4);
		p.cfadim         = (size_t) cfadim;
		p.xmin           = 0;
		p.ymin           = 0;
		p.xmax           = frame_rx - 1;
		p.ymax           = frame_ry - 1;
		p.threads        = 1;
		p.error          = &error;

		if (dobox(&p)) {
			siril_log_color_message(
			    _("Stack (Bayer): dobox failed on frame %d: %s\n"),
			    "red", f, error.last_message[0] ? error.last_message : "(no detail)");
			clearfits(&frame_fit);
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EIO;
		}
		clearfits(&frame_fit);

		++cur_nb;
		gui_iface.set_progress((double) cur_nb / (double) n_included, NULL);
	}

	mpp_imgmap_free(&pixmap);
	clearfits(&output_counts);

	clearfits(out);
	out->rx = out_rx;
	out->ry = out_ry;
	out->naxes[0] = out_rx;
	out->naxes[1] = out_ry;
	out->naxes[2] = channels;
	out->naxis = 3;
	out->bitpix = USHORT_IMG;
	out->orig_bitpix = USHORT_IMG;
	out->type = DATA_USHORT;
	const size_t plane = (size_t) out_rx * (size_t) out_ry;
	out->data = (WORD *) std::calloc(plane * (size_t) channels, sizeof(WORD));
	if (!out->data) {
		clearfits(&output_data);
		return MPP_ENOMEM;
	}
	out->pdata[0] = out->data;
	out->pdata[1] = out->data + plane;
	out->pdata[2] = out->data + 2 * plane;

	const double cast_scale = (cfg->bitdepth == 8) ? 256.0 : 1.0;
	for (size_t k = 0; k < plane * (size_t) channels; ++k) {
		const double v = (double) output_data.fdata[k] * cast_scale;
		out->data[k] = (WORD) (v < 0.0 ? 0
		                     : (v > 65535.0 ? 65535
		                                    : (int) (v + 0.5)));
	}
	clearfits(&output_data);
	return MPP_OK;
}

/* Sequence-side wrapper: reads each included frame via Siril's
 * sequence I/O and forwards to the inner C++ function. */
extern "C" mpp_status_t mpp_stack_apply_stsci(sequence *seq,
                                              const mpp_config_t *cfg,
                                              const mpp_run_t *run,
                                              fits *out) {
	if (!seq || !cfg || !run || !run->included || !run->frame_brightness)
		return MPP_EINVAL;
	const int N = run->num_frames;
	std::vector<cv::Mat> frames_raw(N);
	std::vector<int> included(run->included, run->included + N);
	std::vector<double> brightness(run->frame_brightness,
	                                run->frame_brightness + N);
	int cur_nb = 0, n_inc = 0;
	for (int i = 0; i < N; ++i) if (included[i]) ++n_inc;
	gui_iface.set_progress(PROGRESS_RESET, _("Stack (STScI): reading frames"));
	for (int i = 0; i < N; ++i) {
		if (!included[i]) continue;
		cv::Mat m = mpp::read_full_frame(seq, i);
		if (m.empty()) return MPP_EIO;
		frames_raw[i] = m;
		++cur_nb;
		gui_iface.set_progress(0.5 * (double) cur_nb / (double) n_inc, NULL);
	}
	return mpp::stack_apply_stsci(frames_raw, included, brightness,
	                              run, cfg, out);
}

/* Sequence-side wrapper for Bayer drizzle. Reads each frame raw (no
 * debayer, bypassing com.pref.debayer.open_debayer) and forwards to
 * the inner mpp::stack_apply_bayer.
 *
 * Only SER sequences with a Bayer ColorID are supported. The CFA
 * pattern is derived from the SER header's color_id, matching the
 * 2x2 RLAYER/GLAYER/BLAYER convention Siril's get_compiled_pattern
 * produces (see demosaicing.c:333). */
namespace {
/* Returns FALSE if the sequence has no usable Bayer pattern. cfa
 * is filled with 4 bytes on success. */
bool bayer_cfa_from_ser(struct sequ *seq, unsigned char cfa[4]) {
	if (!seq || seq->type != SEQ_SER || !seq->ser_file) return false;
	const int R = 0, G = 1, B = 2;
	switch (seq->ser_file->color_id) {
		case SER_BAYER_RGGB: cfa[0]=R; cfa[1]=G; cfa[2]=G; cfa[3]=B; return true;
		case SER_BAYER_GRBG: cfa[0]=G; cfa[1]=R; cfa[2]=B; cfa[3]=G; return true;
		case SER_BAYER_GBRG: cfa[0]=G; cfa[1]=B; cfa[2]=R; cfa[3]=G; return true;
		case SER_BAYER_BGGR: cfa[0]=B; cfa[1]=G; cfa[2]=G; cfa[3]=R; return true;
		default: return false;
	}
}
}  // namespace

extern "C" mpp_status_t mpp_stack_apply_bayer(sequence *seq,
                                              const mpp_config_t *cfg,
                                              const mpp_run_t *run,
                                              fits *out) {
	if (!seq || !cfg || !run || !run->included || !run->frame_brightness)
		return MPP_EINVAL;
	if (seq->type != SEQ_SER || !seq->ser_file) {
		siril_log_color_message(_("Bayer drizzle: requires a SER sequence "
		                          "(other formats don't preserve raw Bayer pattern data).\n"),
		                        "red");
		return MPP_EINVAL;
	}
	unsigned char cfa[4];
	if (!bayer_cfa_from_ser(seq, cfa)) {
		siril_log_color_message(_("Bayer drizzle: SER ColorID is not RGGB / GRBG / "
		                          "GBRG / BGGR — cannot derive CFA pattern.\n"),
		                        "red");
		return MPP_EINVAL;
	}

	const int N = run->num_frames;
	std::vector<cv::Mat> frames_raw_bayer(N);
	std::vector<int> included(run->included, run->included + N);
	std::vector<double> brightness(run->frame_brightness,
	                                run->frame_brightness + N);
	int cur_nb = 0, n_inc = 0;
	for (int i = 0; i < N; ++i) if (included[i]) ++n_inc;

	gui_iface.set_progress(PROGRESS_RESET,
	                       _("Stack (Bayer drizzle): reading raw frames"));
	for (int i = 0; i < N; ++i) {
		if (!included[i]) continue;
		/* Read raw Bayer (open_debayer=FALSE) so we get the mosaic
		 * bytes regardless of the user's debayer-on-open preference.
		 * Stage A/B were computed on debayered frames upstream — that's
		 * fine; the run's pixmap is correctly scaled by Stage A/B's
		 * intersection geometry. */
		fits f{};
		if (ser_read_frame(seq->ser_file, i, &f, /*force_float=*/FALSE,
		                   /*open_debayer=*/FALSE)) {
			siril_log_color_message(_("Bayer drizzle: ser_read_frame failed "
			                          "on frame %d\n"), "red", i);
			clearfits(&f);
			return MPP_EIO;
		}
		/* Wrap the single-channel WORD buffer as cv::Mat (CV_16U). */
		cv::Mat m(f.ry, f.rx, CV_16U, f.data);
		frames_raw_bayer[i] = m.clone();   /* own the data; clearfits below frees f */
		clearfits(&f);
		++cur_nb;
		gui_iface.set_progress(0.5 * (double) cur_nb / (double) n_inc, NULL);
	}
	return mpp::stack_apply_bayer(frames_raw_bayer, included, brightness,
	                              run, cfg, cfa, /*cfadim=*/2, out);
}
