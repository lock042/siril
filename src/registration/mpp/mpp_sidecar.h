#ifndef SRC_REGISTRATION_MPP_SIDECAR_H_
#define SRC_REGISTRATION_MPP_SIDECAR_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Persist a completed run (quality, shifts, AP coords, per-AP shifts) so it
 * can be inspected or resumed without re-processing the sequence. Optional
 * feature; deferable per the plan. */
mpp_status_t mpp_sidecar_write(const char *path, const mpp_run_t *run);
mpp_status_t mpp_sidecar_read(const char *path, mpp_run_t **run_out);
void mpp_sidecar_free(mpp_run_t *run);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_SIDECAR_H_ */
