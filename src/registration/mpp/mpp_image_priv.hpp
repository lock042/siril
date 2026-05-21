/*
 * Helpers bridging Siril's `fits` and `sequence` types to OpenCV.
 *
 * Siril stores image planes contiguously as `WORD*` (or `float*`) in
 * either `fit->data` / `fit->fdata` for mono or `fit->pdata[ch]` /
 * `fit->fpdata[ch]` for multi-layer (RGB) frames. The helpers here wrap
 * those buffers as cv::Mat headers without copying.
 *
 * Frame orientation: Siril SER frames are stored top-down (`fit->top_down
 * == true`), regular FITS bottom-up. The cv::Mat header doesn't carry that
 * flag — Phase 5a / PSS oracle equivalence was verified against PNG inputs
 * (top-down) and matches PSS exactly. The caller is responsible for
 * orientation correctness; mpp processing is rotation-invariant for
 * registration/stacking but the output preserves whichever orientation the
 * input had.
 */
#ifndef SRC_REGISTRATION_MPP_IMAGE_PRIV_HPP_
#define SRC_REGISTRATION_MPP_IMAGE_PRIV_HPP_

#include <opencv2/core.hpp>

#include "core/siril.h"
#include "io/sequence.h"
#include "registration/mpp/mpp_config.h"

namespace mpp {

/* Wrap a single layer of a `fits` image as a cv::Mat header (no copy).
 * The returned Mat's lifetime is bounded by the fits's pixel buffer. */
inline cv::Mat wrap_fits_layer(const fits *f, int layer) {
	if (!f) return cv::Mat();
	const int rows = (int) f->ry;
	const int cols = (int) f->rx;
	if (f->type == DATA_USHORT) {
		WORD *p = (layer == 0 && f->naxes[2] <= 1) ? f->data : f->pdata[layer];
		return cv::Mat(rows, cols, CV_16U, (void *) p);
	}
	if (f->type == DATA_FLOAT) {
		float *p = (layer == 0 && f->naxes[2] <= 1) ? f->fdata : f->fpdata[layer];
		return cv::Mat(rows, cols, CV_32F, (void *) p);
	}
	return cv::Mat();
}

/* Number of layers in the image (1 for mono, 3 for RGB). */
inline int fits_layer_count(const fits *f) {
	if (!f) return 0;
	return f->naxis >= 3 ? (int) f->naxes[2] : 1;
}

/* PSS bit-depth convention from fit->orig_bitpix (Siril stores both 8-bit
 * and 16-bit SER data in WORD; orig_bitpix preserves the source range). */
inline int fits_bitdepth(const fits *f) {
	if (!f) return 16;
	return mpp_bitdepth_from_fits_bitpix(f->orig_bitpix);
}

/* For colour input, the "analysis" frame used for ranking, alignment, AP
 * placement, and per-AP shift is the green channel by default (PSS's
 * default `color_index` for non-luminance pathways is 1 → green). For
 * mono, just layer 0. */
inline int analysis_layer_for(const fits *f) {
	return fits_layer_count(f) >= 3 ? 1 : 0;
}

}  // namespace mpp

#endif  /* SRC_REGISTRATION_MPP_IMAGE_PRIV_HPP_ */
