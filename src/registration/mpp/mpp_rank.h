#ifndef SRC_REGISTRATION_MPP_RANK_H_
#define SRC_REGISTRATION_MPP_RANK_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Compute per-frame Laplace-sigma quality for every frame in `seq`.
 * On success, *quality_out is a malloc'd array of length seq->number; caller
 * frees. Phase 1. */
mpp_status_t mpp_rank_sequence(sequence *seq, const mpp_config_t *cfg,
                               double **quality_out);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_RANK_H_ */
