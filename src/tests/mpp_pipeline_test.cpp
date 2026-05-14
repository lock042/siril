/*
 * Integration test for the mpp orchestrator. Skeleton — Phase 6 will populate
 * with a synthetic-sequence end-to-end run against the PSS oracle.
 */
#include <criterion/criterion.h>

extern "C" {
#include "registration/mpp.h"
}

Test(mpp_pipeline, register_mpp_stub_returns_enotimpl) {
	cr_assert_eq(register_mpp(nullptr), MPP_ENOTIMPL);
}
