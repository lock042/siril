#include "registration/mpp.h"
#include "registration/mpp_shift.h"

extern "C" mpp_status_t mpp_shift_compute(struct sequence *seq,
                                          const mpp_config_t *cfg,
                                          const mpp_aps_t *aps,
                                          const struct fits *ref,
                                          const int *global_shifts,
                                          mpp_shifts_t **shifts_out) {
	(void) seq; (void) cfg; (void) aps; (void) ref; (void) global_shifts;
	if (shifts_out) *shifts_out = nullptr;
	return MPP_ENOTIMPL;
}

extern "C" void mpp_shift_free(mpp_shifts_t *shifts) {
	(void) shifts;
}
