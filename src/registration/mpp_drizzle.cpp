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
}

#include "registration/mpp_stack_priv.hpp"   /* mpp::stack_one_dim_weight */

/* read_full_frame is declared in mpp.cpp's anonymous mpp:: namespace.
 * The pieces we need are accessed via a tiny extern declaration so we
 * don't have to expose the helper publicly. */
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

extern "C" mpp_status_t mpp_stack_apply_stsci(sequence *seq, const mpp_config_t *cfg,
                                              const mpp_run_t *run, fits *out) {
	if (!seq || !cfg || !run || !out || !run->aps || !run->shifts
	 || !run->included || !run->frame_brightness)
		return MPP_EINVAL;

	const int frame_rx = run->frame_cols;
	const int frame_ry = run->frame_rows;
	if (frame_rx <= 0 || frame_ry <= 0) return MPP_EINVAL;

	const int int_y_lo = run->intersection[0];
	const int int_y_hi = run->intersection[1];
	const int int_x_lo = run->intersection[2];
	const int int_x_hi = run->intersection[3];
	const int intersection_rx = int_x_hi - int_x_lo;
	const int intersection_ry = int_y_hi - int_y_lo;
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

	/* Median brightness across included frames — same per-frame
	 * normalisation that Phase 5a's bicubic path uses. */
	std::vector<double> bright_inc;
	bright_inc.reserve(run->num_frames);
	for (int i = 0; i < run->num_frames; ++i)
		if (run->included[i]) bright_inc.push_back(run->frame_brightness[i]);
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

	for (int f = 0; f < run->num_frames; ++f) {
		if (!run->included[f]) continue;

		cv::Mat raw = mpp::read_full_frame(seq, f);
		if (raw.empty()) {
			mpp_imgmap_free(&pixmap);
			clearfits(&output_data);
			clearfits(&output_counts);
			return MPP_EIO;
		}

		const double scale = median_brightness / (run->frame_brightness[f] + 1e-7);

		fits frame_fit{};
		if (!cv_to_planar_fits(raw, scale, &frame_fit)) {
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
		p.scale          = driz_args.scale;
		p.pixel_fraction = driz_args.pixel_fraction;
		p.weight_scale   = 1.0f;
		p.in_units       = unit_counts;
		p.out_units      = unit_counts;
		p.exposure_time  = 1.0f;
		p.data           = &frame_fit;
		p.weights        = nullptr;       /* unit per-pixel weight for now (5b.2 v1) */
		p.pixmap         = &pixmap;
		p.output_data    = &output_data;
		p.output_counts  = &output_counts;
		p.xmin           = 0;
		p.ymin           = 0;
		p.xmax           = frame_rx - 1;
		p.ymax           = frame_ry - 1;
		p.threads        = 1;             /* dobox runs single-threaded per frame */
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

	/* Per-channel normalise: data /= counts (with epsilon floor — empty
	 * cells stay zero). */
	{
		const size_t plane = (size_t) out_rx * (size_t) out_ry;
		for (int c = 0; c < channels; ++c) {
			float *od = output_data.fdata   + (size_t) c * plane;
			float *oc = output_counts.fdata + (size_t) c * plane;
			for (size_t k = 0; k < plane; ++k) {
				od[k] = (oc[k] > 1e-7f) ? od[k] / oc[k] : 0.0f;
			}
		}
	}
	clearfits(&output_counts);

	/* Pack the float-stacked result into the caller's uint16 fits.
	 * Layout: planar R/G/B with pdata[0..2] pointing into a single
	 * contiguous buffer — same shape Phase 5a's mpp_stack_apply
	 * produces. */
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

	for (size_t k = 0; k < plane * (size_t) channels; ++k) {
		const float v = output_data.fdata[k];
		out->data[k] = (WORD) (v < 0.f ? 0 : (v > 65535.f ? 65535 : (int) (v + 0.5f)));
	}
	clearfits(&output_data);
	return MPP_OK;
}
