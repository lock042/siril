#ifndef SRC_REGISTRATION_MPP_FRAMES_H_
#define SRC_REGISTRATION_MPP_FRAMES_H_

#include <stdbool.h>

#include "core/siril.h"            /* for `sequence` and `fits` typedefs */
#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Read frame `idx` from `seq` into `dest` (caller-allocated). Dispatches on
 * sequence type just like Siril's `seq_read_frame`, with explicit control
 * over SER debayering instead of the global debayer-on-open preference:
 *   ser_force_debayer = true  → a Bayer SER comes back as 3-layer RGB
 *     (RCD-debayered by ser_read_frame) regardless of the preference.
 *     Used by the full-frame stacking reads, where demosaic quality
 *     carries into the output.
 *   ser_force_debayer = false → the read follows ser_file->debayer_type_ser.
 *     Analysis callers force that to SER_MONO for the duration of a stage
 *     (see SerAnalysisRawGuard in mpp.cpp) so they get the raw mosaic and
 *     debayer it cheaply with cv::cvtColor.
 * Pass `force_float=false` to keep Siril's WORD storage convention. */
int mpp_seq_read_frame(sequence *seq, int idx, fits *dest,
                       bool force_float, int thread_id,
                       bool ser_force_debayer);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_FRAMES_H_ */
