#ifndef SRC_REGISTRATION_MPP_CONFIG_H_
#define SRC_REGISTRATION_MPP_CONFIG_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Populate `cfg` with PSS defaults. Definition of mpp_config lives here once
 * Phase 1+ start consuming individual fields. */
mpp_status_t mpp_config_defaults(mpp_config_t *cfg);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_CONFIG_H_ */
