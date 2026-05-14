#ifndef SRC_REGISTRATION_MPP_AP_H_
#define SRC_REGISTRATION_MPP_AP_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct fits;

/* One alignment point. Coordinates are pixel-space on the reference (mean)
 * frame. The patch bounds are extended to the frame border for edge APs
 * (PSS new_alignment_point extend_*_low/high logic). `structure` is the
 * post-normalisation Laplace/gradient measure used to filter weak APs. */
typedef struct mpp_ap_record {
	int y, x;
	int box_y_low, box_y_high, box_x_low, box_x_high;
	int patch_y_low, patch_y_high, patch_x_low, patch_x_high;
	double structure;        /* normalised to [0, 1] after grid construction */
} mpp_ap_record_t;

struct mpp_aps {
	int count;
	int dropped_dim;          /* APs rejected by brightness/contrast filter */
	int dropped_structure;    /* APs rejected by structure threshold */
	mpp_ap_record_t *records; /* malloc'd, length = count */
};

/* Place a staggered AP grid on the reference (mean) frame and apply PSS's
 * brightness, contrast, dim-fraction COM re-centring, and structure
 * filters. On success the caller owns the result and frees via mpp_ap_free. */
mpp_status_t mpp_ap_place(const struct fits *ref, const mpp_config_t *cfg,
                          mpp_aps_t **aps_out);

void mpp_ap_free(mpp_aps_t *aps);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_AP_H_ */
