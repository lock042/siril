#ifndef SRC_REGISTRATION_MPP_FRAMES_H_
#define SRC_REGISTRATION_MPP_FRAMES_H_

#include <stdbool.h>

#include "core/siril.h"            /* for `sequence` and `fits` typedefs */
#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Read frame `idx` from `seq` into `dest` (caller-allocated). Dispatches on
 * sequence type just like Siril's `seq_read_frame`, except for SER it
 * unconditionally passes `open_debayer=TRUE` so Bayer SERs come back as
 * 3-layer RGB (Phase 6 requirement). Pass `force_float=false` to keep
 * Siril's WORD storage convention. */
int mpp_seq_read_frame(sequence *seq, int idx, fits *dest,
                       bool force_float, int thread_id);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_FRAMES_H_ */
