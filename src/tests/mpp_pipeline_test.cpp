/*
 * Phase 6's orchestrator pulls in Siril runtime globals (`com` etc.) via
 * siril_log_message that the unit-test binary can't resolve. Integration
 * testing of register_mpp / mpp_analyze / mpp_compute_shifts /
 * mpp_stack_apply happens via the real Siril CLI on the test-big.ser
 * fixture (Phase 6 commit 5).
 *
 * Per-phase algorithmic equivalence is locked in by mpp_rank_test,
 * mpp_align_test, mpp_ap_test, mpp_shift_test, and mpp_stack_test —
 * those don't go through register_mpp and so don't pull mpp.cpp.o into
 * their link footprint.
 */
#include <criterion/criterion.h>

Test(mpp_pipeline, placeholder_real_tests_live_in_phase6_commit5) {
	cr_assert(true);
}
