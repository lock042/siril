#ifndef SRC_REGISTRATION_MPP_STACK_H_
#define SRC_REGISTRATION_MPP_STACK_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sequence;
struct fits;

/* Pluggable per-frame resample backend; see mpp_drizzle.h for STScI variants. */
typedef enum {
	MPP_RESAMPLE_BICUBIC = 0,   /* Phase 5a — PSS-faithful (cv::resize INTER_LINEAR
	                              actually, despite the kind's name — PSS uses
	                              INTER_LINEAR for its drizzle resize, not cubic) */
	MPP_RESAMPLE_STSCI,         /* Phase 5b — debayered + dobox() */
	MPP_RESAMPLE_BAYER_STSCI    /* Phase 5b — raw Bayer + dobox() */
} mpp_resample_kind_t;

/* Top-N-per-AP weighted stacking. Phase 5. */
mpp_status_t mpp_stack(struct sequence *seq, const mpp_config_t *cfg,
                       const mpp_aps_t *aps, const mpp_shifts_t *shifts,
                       const int *global_shifts,
                       mpp_resample_kind_t backend, int upscale,
                       struct fits *stacked_out);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_STACK_H_ */
