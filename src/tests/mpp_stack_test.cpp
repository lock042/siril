/*
 * Criterion tests for mpp_stack. Skeleton — Phase 5a will populate.
 */
#include <criterion/criterion.h>

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_stack.h"
}

Test(mpp_stack, stub_returns_enotimpl) {
	cr_assert_eq(mpp_stack(nullptr, nullptr, nullptr, nullptr, nullptr,
	                       MPP_RESAMPLE_BICUBIC, 1, nullptr),
	             MPP_ENOTIMPL);
}
