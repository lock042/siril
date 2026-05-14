#include "registration/mpp.h"
#include "registration/mpp_drizzle.h"

extern "C" mpp_status_t mpp_drizzle_resample(const struct fits *frame_in,
                                             int global_dy, int global_dx,
                                             int upscale, double pixfrac,
                                             int bayer_mode,
                                             struct fits *out_buf,
                                             struct fits *out_weights) {
	(void) frame_in; (void) global_dy; (void) global_dx;
	(void) upscale; (void) pixfrac; (void) bayer_mode;
	(void) out_buf; (void) out_weights;
	return MPP_ENOTIMPL;
}
