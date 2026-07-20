#ifndef SRC_REGISTRATION_MPP_AP_H_
#define SRC_REGISTRATION_MPP_AP_H_

#include "registration/mpp.h"

#ifdef __cplusplus
extern "C" {
#endif


/* One alignment point. Coordinates are pixel-space on the reference (mean)
 * frame. The patch bounds are extended to the frame border for edge APs.
 * `structure` is the post-normalisation Laplace/gradient measure used to
 * filter weak APs. */
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
	/* Set by the AP editor's manual mutations (add/move/remove/resize). Once
	 * true, Register must not auto-regenerate the grid when the placement
	 * spinners differ from the analysis-time cfg — the manual grid (which may
	 * mix AP sizes) takes precedence. Cleared by a full auto-placement
	 * (mpp_ap_replace), which produces a uniform grid that matches cfg. */
	int user_edited;
	mpp_ap_record_t *records; /* malloc'd, length = count */
};

/* Place a staggered AP grid on the reference (mean) frame and apply the
 * brightness, contrast, dim-fraction COM re-centring, and structure
 * filters. On success the caller owns the result and frees via mpp_ap_free. */
mpp_status_t mpp_ap_place(const fits *ref, const mpp_config_t *cfg,
                          mpp_aps_t **aps_out);

void mpp_ap_free(mpp_aps_t *aps);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_AP_H_ */
