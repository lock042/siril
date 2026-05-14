/*
 * Criterion tests for mpp_align. Skeleton — Phase 2 will populate.
 */
#include <criterion/criterion.h>

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_align.h"
}

Test(mpp_align, stub_returns_enotimpl) {
	int patch[4] = {0, 0, 0, 0};
	cr_assert_eq(mpp_align_global(nullptr, nullptr, nullptr, nullptr, patch, nullptr),
	             MPP_ENOTIMPL);
}
