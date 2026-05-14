#include "registration/mpp.h"
#include "registration/mpp_rank.h"

extern "C" mpp_status_t mpp_rank_sequence(struct sequence *seq,
                                          const mpp_config_t *cfg,
                                          double **quality_out) {
	(void) seq; (void) cfg;
	if (quality_out) *quality_out = nullptr;
	return MPP_ENOTIMPL;
}
