#include "registration/mpp.h"
#include "registration/mpp_ap.h"

extern "C" mpp_status_t mpp_ap_place(const struct fits *ref,
                                     const mpp_config_t *cfg,
                                     mpp_aps_t **aps_out) {
	(void) ref; (void) cfg;
	if (aps_out) *aps_out = nullptr;
	return MPP_ENOTIMPL;
}

extern "C" void mpp_ap_free(mpp_aps_t *aps) {
	(void) aps;
}
