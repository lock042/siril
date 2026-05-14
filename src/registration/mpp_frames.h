#ifndef SRC_REGISTRATION_MPP_FRAMES_H_
#define SRC_REGISTRATION_MPP_FRAMES_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sequence;
struct fits;

/* Load frame `idx` from the sequence as a mono float32 buffer suitable for
 * downstream analysis (rank, align, shift). Caller owns `out`. Bayer
 * sequences are debayered here. Full API contract finalised in Phase 1.3. */
mpp_status_t mpp_frames_load_mono(struct sequence *seq, int idx, const mpp_config_t *cfg,
                                  struct fits *out);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_FRAMES_H_ */
