/*
 * Criterion tests for mpp_drizzle. Skeleton — Phase 5b will populate.
 */
#include <criterion/criterion.h>

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_drizzle.h"
}

Test(mpp_drizzle, stub_returns_enotimpl) {
	cr_assert_eq(mpp_drizzle_resample(nullptr, 0, 0, 2, 0.7, 0, nullptr, nullptr),
	             MPP_ENOTIMPL);
}
