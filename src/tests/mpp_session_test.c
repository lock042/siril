/*
 * mpp_session: the multi-sequence derotation session data model — adding,
 * channel tagging, reference designation, the shared-epoch computation and
 * channel grouping. Pure logic, no I/O.
 */
#include <criterion/criterion.h>

#include "core/siril.h"
#include "registration/mpp/mpp_session.h"

/* Linked siril objects reference these globals. */
cominfo com;
fits *gfit = NULL;

Test(mpp_session, add_tag_reference_epoch) {
	mpp_session_t *s = mpp_session_new(PLANET_JUPITER, 2);
	cr_assert_not_null(s);
	cr_assert_eq(s->reference, -1);

	cr_assert_eq(mpp_session_add(s, "R1", MPP_CHAN_R, 100.0, 100.1, 1000), 0);
	cr_assert_eq(s->reference, 0);                      /* first added = reference */
	cr_assert_eq(mpp_session_add(s, "R2", MPP_CHAN_R, 100.2, 100.3, 1000), 1);
	cr_assert_eq(mpp_session_add(s, "G1", MPP_CHAN_G, 100.4, 100.5, 1000), 2);
	cr_assert_eq(mpp_session_add(s, "R1", MPP_CHAN_B, 0, 0, 0), -1);  /* duplicate */
	cr_assert_eq(s->count, 3);

	cr_assert(mpp_session_set_reference(s, 2));
	cr_assert_eq(s->reference, 2);

	mpp_channel_t ch[MPP_CHAN_COUNT];
	cr_assert_eq(mpp_session_channels(s, ch), 2);
	cr_assert_eq(ch[0], MPP_CHAN_R);                    /* R before G */
	cr_assert_eq(ch[1], MPP_CHAN_G);

	int idx[8];
	cr_assert_eq(mpp_session_channel_seqs(s, MPP_CHAN_R, idx), 2);
	cr_assert_eq(idx[0], 0);
	cr_assert_eq(idx[1], 1);
	cr_assert_eq(mpp_session_channel_seqs(s, MPP_CHAN_G, idx), 1);
	cr_assert_eq(idx[0], 2);

	/* union span [100.0, 100.5] -> midpoint 100.25 */
	cr_assert_float_eq(mpp_session_compute_epoch(s), 100.25, 1e-9);
	cr_assert_float_eq(s->epoch_jd, 100.25, 1e-9);

	/* disk fit is recorded and flags the entry */
	cr_assert_not(s->seqs[1].has_fit);
	cr_assert(mpp_session_set_fit(s, 1, 50.0, 60.0, 30.0, 28.0, 12.0, -1.0,
	                              480, 640, /*bayer=*/-1, 25.0));
	cr_assert(s->seqs[1].has_fit);
	cr_assert_float_eq(s->seqs[1].cx, 50.0, 1e-9);
	cr_assert_float_eq(s->seqs[1].r_eq, 30.0, 1e-9);
	cr_assert_float_eq(s->seqs[1].parity, -1.0, 1e-9);
	cr_assert_eq(s->seqs[1].frame_cols, 640);

	mpp_session_free(s);
}

Test(mpp_session, remove_keeps_reference_consistent) {
	mpp_session_t *s = mpp_session_new(PLANET_SATURN, 3);
	mpp_session_add(s, "A", MPP_CHAN_MONO, 0, 1, 10);
	mpp_session_add(s, "B", MPP_CHAN_MONO, 0, 1, 10);
	mpp_session_add(s, "C", MPP_CHAN_MONO, 0, 1, 10);
	mpp_session_set_reference(s, 2);                    /* C */

	cr_assert(mpp_session_remove(s, 0));                /* drop A */
	cr_assert_eq(s->count, 2);
	cr_assert_eq(s->reference, 1);                      /* C shifted down to 1 */

	cr_assert(mpp_session_remove(s, 1));                /* drop the reference (C) */
	cr_assert_eq(s->reference, 0);                      /* falls back to first */

	mpp_session_free(s);
}
