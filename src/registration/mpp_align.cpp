#include "registration/mpp.h"
#include "registration/mpp_align.h"

extern "C" mpp_status_t mpp_align_global(struct sequence *seq,
                                         const mpp_config_t *cfg,
                                         const double *quality,
                                         int *shifts_out,
                                         int patch_yxyx_out[4],
                                         struct fits *avg_ref_out) {
	(void) seq; (void) cfg; (void) quality;
	(void) shifts_out; (void) patch_yxyx_out; (void) avg_ref_out;
	return MPP_ENOTIMPL;
}
