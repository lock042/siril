/*
 * mpp_multisource: the global-index layout that concatenates several
 * sequences' frames into one FrameProvider index space (Option B, C4). Pure
 * index arithmetic + provider dispatch; no I/O.
 */
#include <criterion/criterion.h>

#include "registration/mpp/mpp_multisource.hpp"

using namespace mpp;

/* layout of three sequences [3, 2, 4] -> 9 frames, contiguous */
Test(mpp_multisource, layout_offsets) {
	MultiSourceLayout L;
	cr_assert_eq(L.add(3), 0);
	cr_assert_eq(L.add(2), 1);
	cr_assert_eq(L.add(4), 2);

	cr_assert_eq(L.count(), 3);
	cr_assert_eq(L.total(), 9);
	cr_assert_eq(L.frames_of(0), 3);
	cr_assert_eq(L.frames_of(1), 2);
	cr_assert_eq(L.frames_of(2), 4);
	cr_assert_eq(L.base_of(0), 0);
	cr_assert_eq(L.base_of(1), 3);
	cr_assert_eq(L.base_of(2), 5);
}

/* every global index maps to the right (seq, local); boundaries included */
Test(mpp_multisource, locate_all) {
	MultiSourceLayout L;
	L.add(3); L.add(2); L.add(4);

	const int exp_seq[9]   = { 0, 0, 0, 1, 1, 2, 2, 2, 2 };
	const int exp_local[9] = { 0, 1, 2, 0, 1, 0, 1, 2, 3 };
	for (int g = 0; g < 9; ++g) {
		int s = -1, l = -1;
		cr_assert(L.locate(g, &s, &l), "global %d", g);
		cr_assert_eq(s, exp_seq[g], "global %d seq", g);
		cr_assert_eq(l, exp_local[g], "global %d local", g);
	}
	/* out of range leaves outputs untouched and returns false */
	int s = 42, l = 42;
	cr_assert_not(L.locate(-1, &s, &l));
	cr_assert_not(L.locate(9, &s, &l));
	cr_assert_eq(s, 42);
	cr_assert_eq(l, 42);
}

/* the composed provider dispatches each global index to read(seq, local);
 * out-of-range -> empty Mat (engine treats it as an excluded frame) */
Test(mpp_multisource, provider_dispatch) {
	MultiSourceLayout L;
	L.add(3); L.add(2); L.add(4);

	auto read = [](int seq, int local) -> cv::Mat {
		cv::Mat m(1, 1, CV_32S);
		m.at<int>(0, 0) = seq * 100 + local;   /* encode the routing */
		return m;
	};
	auto provider = make_multisource_provider(L, read);

	const int exp[9] = { 0, 1, 2, 100, 101, 200, 201, 202, 203 };
	for (int g = 0; g < 9; ++g) {
		cv::Mat m = provider(g);
		cr_assert(!m.empty(), "global %d", g);
		cr_assert_eq(m.at<int>(0, 0), exp[g], "global %d", g);
	}
	cr_assert(provider(9).empty());
	cr_assert(provider(-1).empty());
}

/* a single sequence reduces to a 1:1 pass-through provider */
Test(mpp_multisource, single_sequence_passthrough) {
	MultiSourceLayout L;
	cr_assert_eq(L.add(5), 0);
	cr_assert_eq(L.total(), 5);
	for (int g = 0; g < 5; ++g) {
		int s = -1, l = -1;
		cr_assert(L.locate(g, &s, &l));
		cr_assert_eq(s, 0);
		cr_assert_eq(l, g);
	}
}
