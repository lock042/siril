#include "registration/mpp.h"
#include "registration/mpp_frames.h"

extern "C" mpp_status_t mpp_frames_load_mono(struct sequence *seq, int idx,
                                             const mpp_config_t *cfg,
                                             struct fits *out) {
	(void) seq; (void) idx; (void) cfg; (void) out;
	return MPP_ENOTIMPL;
}
