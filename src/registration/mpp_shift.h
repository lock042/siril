#ifndef SRC_REGISTRATION_MPP_SHIFT_H_
#define SRC_REGISTRATION_MPP_SHIFT_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sequence;
struct fits;

/* For each AP in `aps`, compute per-frame local shifts via two-phase
 * multilevel correlation against the reference frame. Phase 4. */
mpp_status_t mpp_shift_compute(struct sequence *seq, const mpp_config_t *cfg,
                               const mpp_aps_t *aps, const struct fits *ref,
                               const int *global_shifts,
                               mpp_shifts_t **shifts_out);

void mpp_shift_free(mpp_shifts_t *shifts);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_SHIFT_H_ */
