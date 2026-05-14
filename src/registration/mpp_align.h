#ifndef SRC_REGISTRATION_MPP_ALIGN_H_
#define SRC_REGISTRATION_MPP_ALIGN_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sequence;
struct fits;

/* Global frame alignment. Picks the best-structured patch on the best frame,
 * cross-correlates every frame against it via two-phase multilevel matching,
 * and produces per-frame integer (dy, dx) plus an average reference frame.
 * Phase 2. */
mpp_status_t mpp_align_global(struct sequence *seq, const mpp_config_t *cfg,
                              const double *quality,
                              int *shifts_out,            /* 2 * seq->number */
                              int patch_yxyx_out[4],
                              struct fits *avg_ref_out);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_ALIGN_H_ */
