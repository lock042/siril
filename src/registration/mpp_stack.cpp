#include "registration/mpp.h"
#include "registration/mpp_stack.h"

extern "C" mpp_status_t mpp_stack(struct sequence *seq,
                                  const mpp_config_t *cfg,
                                  const mpp_aps_t *aps,
                                  const mpp_shifts_t *shifts,
                                  const int *global_shifts,
                                  mpp_resample_kind_t backend, int upscale,
                                  struct fits *stacked_out) {
	(void) seq; (void) cfg; (void) aps; (void) shifts; (void) global_shifts;
	(void) backend; (void) upscale; (void) stacked_out;
	return MPP_ENOTIMPL;
}
