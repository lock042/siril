/*
 * Criterion tests for mpp_ap. Skeleton — Phase 3 will populate.
 */
#include <criterion/criterion.h>

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_ap.h"
}

Test(mpp_ap, stub_returns_enotimpl) {
	mpp_aps_t *aps = nullptr;
	cr_assert_eq(mpp_ap_place(nullptr, nullptr, &aps), MPP_ENOTIMPL);
	cr_assert_null(aps);
	mpp_ap_free(aps); /* no-op on null */
}
