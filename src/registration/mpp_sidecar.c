#include "registration/mpp.h"
#include "registration/mpp_sidecar.h"

mpp_status_t mpp_sidecar_write(const char *path, const mpp_run_t *run) {
	(void) path; (void) run;
	return MPP_ENOTIMPL;
}

mpp_status_t mpp_sidecar_read(const char *path, mpp_run_t **run_out) {
	(void) path;
	if (run_out) *run_out = NULL;
	return MPP_ENOTIMPL;
}

void mpp_sidecar_free(mpp_run_t *run) {
	(void) run;
}
