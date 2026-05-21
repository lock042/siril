#ifndef SRC_REGISTRATION_MPP_SHIFT_H_
#define SRC_REGISTRATION_MPP_SHIFT_H_

#include <stdint.h>

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Per-AP per-frame shifts. Indexing convention (frame-major):
 *
 *   dy = shifts[2 * (frame * num_aps + ap) + 0];
 *   dx = shifts[2 * (frame * num_aps + ap) + 1];
 *   ok = success[frame * num_aps + ap];
 *
 * dy / dx are doubles to admit PSS's optional phase-2 sub-pixel parabolic
 * fit; in the default integer mode they're integer-valued doubles. The
 * sign convention matches Phase 2: positive (dy, dx) means the *frame*
 * needs that shift to align with the reference (post-global-alignment). */
struct mpp_shifts {
	int num_frames;
	int num_aps;
	double *shifts;       /* malloc'd, length 2 * num_frames * num_aps */
	uint8_t *success;     /* malloc'd, length     num_frames * num_aps */
	int failure_counter;  /* diagnostic: count of (frame, AP) shifts that failed */
};

/* For each AP in `aps`, compute per-frame local shifts via two-phase
 * multilevel correlation against the reference frame. Phase 4. */
mpp_status_t mpp_shift_compute(sequence *seq, const mpp_config_t *cfg,
                               const mpp_aps_t *aps, const fits *ref,
                               const int *global_shifts,
                               mpp_shifts_t **shifts_out);

void mpp_shift_free(mpp_shifts_t *shifts);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_SHIFT_H_ */
