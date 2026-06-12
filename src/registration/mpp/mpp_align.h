#ifndef SRC_REGISTRATION_MPP_ALIGN_H_
#define SRC_REGISTRATION_MPP_ALIGN_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Global frame alignment. Picks the best-structured patch on the best frame,
 * cross-correlates every frame against it via two-phase multilevel matching,
 * and produces per-frame integer (dy, dx) plus an average reference frame.
 * Phase 2. */
mpp_status_t mpp_align_global(sequence *seq, const mpp_config_t *cfg,
                              const double *quality,
                              int *shifts_out,            /* 2 * seq->number */
                              int patch_yxyx_out[4],
                              fits *avg_ref_out);

/* TRUE if `fit` shows a disc-like object: a bright body (round, or
 * elongated like Saturn) completely surrounded by dark sky, as opposed
 * to a lunar/solar surface that reaches the frame edges. Used by the
 * GUI to pick the default global alignment mode (Planet vs Surface). */
gboolean mpp_frame_has_disc(const fits *fit);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_ALIGN_H_ */
