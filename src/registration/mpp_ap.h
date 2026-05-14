#ifndef SRC_REGISTRATION_MPP_AP_H_
#define SRC_REGISTRATION_MPP_AP_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct fits;

/* Place a staggered AP grid over `ref` and filter by PSS's brightness/structure
 * thresholds. On success the caller owns the result; full API lands in Phase 3. */
mpp_status_t mpp_ap_place(const struct fits *ref, const mpp_config_t *cfg,
                          mpp_aps_t **aps_out);

void mpp_ap_free(mpp_aps_t *aps);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_AP_H_ */
