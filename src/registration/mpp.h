#ifndef SRC_REGISTRATION_MPP_H_
#define SRC_REGISTRATION_MPP_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct registration_args;

typedef enum {
	MPP_OK = 0,
	MPP_ENOTIMPL = -1,
	MPP_EINVAL = -2,
	MPP_ENOMEM = -3,
	MPP_EIO = -4,
	MPP_ENODATA = -5
} mpp_status_t;

/* Forward decls for the per-phase types. Definitions land with their owning
 * module as Phases 1+ flesh them out; consumers should treat them as opaque. */
typedef struct mpp_config mpp_config_t;
typedef struct mpp_run mpp_run_t;
typedef struct mpp_aps mpp_aps_t;
typedef struct mpp_shifts mpp_shifts_t;

/* Orchestrator: chains rank → global-align → AP grid → per-AP shifts → stack.
 * Wired up in Phase 6. */
int register_mpp(struct registration_args *regargs);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_H_ */
