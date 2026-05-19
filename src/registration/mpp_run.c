/*
 * mpp_run_t alloc / free. Lives in its own TU (free of Siril runtime
 * dependencies like siril_log_message) so the sidecar I/O module can link
 * it without dragging mpp.cpp.o (which references com / siril_log_message
 * via the orchestrator stages) into test binaries.
 */

#include <stdlib.h>

#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_shift.h"

mpp_run_t *mpp_run_alloc(void) {
	return (mpp_run_t *) calloc(1, sizeof(mpp_run_t));
}

void mpp_run_free(mpp_run_t *run) {
	if (!run) return;
	free(run->cfg);
	free(run->quality);
	free(run->frame_brightness);
	free(run->included);
	free(run->global_shifts);
	free(run->mean_frame_data);
	free(run->best_frame_indices);
	mpp_ap_free(run->aps);
	mpp_shift_free(run->shifts);
	/* mpp_cache_free lives in mpp.cpp (needs the C++ cv::Mat dtor).
	 * Tests link with --unresolved-symbols=ignore-all and never
	 * populate the cache, so the guard skips the call before the
	 * unresolved symbol would matter. */
	if (run->cache) mpp_cache_free(run->cache);
	free(run);
}
