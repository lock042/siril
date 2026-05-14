#ifndef SRC_REGISTRATION_MPP_DRIZZLE_H_
#define SRC_REGISTRATION_MPP_DRIZZLE_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Wraps the STScI dobox() drizzle backend for both debayered and raw-Bayer
 * resampling paths. Real implementation lands in Phase 5b. */
mpp_status_t mpp_drizzle_resample(const fits *frame_in,
                                  int global_dy, int global_dx,
                                  int upscale, double pixfrac,
                                  int bayer_mode,                /* 0 = mono/RGB, 1 = raw Bayer */
                                  fits *out_buf,
                                  fits *out_weights);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_DRIZZLE_H_ */
