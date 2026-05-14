/*
 * Criterion tests for mpp_shift. Skeleton — Phase 4 will populate.
 */
#include <criterion/criterion.h>

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_shift.h"
}

Test(mpp_shift, stub_returns_enotimpl) {
	mpp_shifts_t *s = nullptr;
	cr_assert_eq(mpp_shift_compute(nullptr, nullptr, nullptr, nullptr, nullptr, &s),
	             MPP_ENOTIMPL);
	cr_assert_null(s);
	mpp_shift_free(s);
}
