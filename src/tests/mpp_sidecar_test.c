/*
 * Criterion tests for mpp_sidecar. Skeleton — populated alongside whichever
 * phase first persists run state to disk (deferable).
 */
#include <criterion/criterion.h>

#include "registration/mpp.h"
#include "registration/mpp_sidecar.h"

Test(mpp_sidecar, stub_returns_enotimpl) {
	mpp_run_t *r = NULL;
	cr_assert_eq(mpp_sidecar_read("/dev/null", &r), MPP_ENOTIMPL);
	cr_assert_null(r);
	mpp_sidecar_free(r);
}
