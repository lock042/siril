/*
 * Criterion tests for mpp_rank. Skeleton — Phase 1 will populate.
 */
#include <criterion/criterion.h>

extern "C" {
#include "registration/mpp.h"
#include "registration/mpp_rank.h"
}

Test(mpp_rank, stub_returns_enotimpl) {
	double *q = nullptr;
	cr_assert_eq(mpp_rank_sequence(nullptr, nullptr, &q), MPP_ENOTIMPL);
	cr_assert_null(q);
}
