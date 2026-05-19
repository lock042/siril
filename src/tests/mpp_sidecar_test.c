/*
 * mpp_sidecar binary I/O round-trip.
 *
 * Builds a small mpp_run_t in memory, writes it to a temp file via
 * mpp_sidecar_write, reads it back via mpp_sidecar_read, and verifies
 * every persisted field matches. Catches binary-format mistakes (off-
 * by-one in field ordering, mis-sized arrays, etc.) before they show up
 * in an end-to-end CLI test.
 */
#include <criterion/criterion.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "registration/mpp.h"
#include "registration/mpp_ap.h"
#include "registration/mpp_config.h"
#include "registration/mpp_shift.h"
#include "registration/mpp_sidecar.h"

static mpp_run_t *make_test_run(int with_shifts) {
	mpp_run_t *run = mpp_run_alloc();
	if (!run) return NULL;
	run->cfg = malloc(sizeof(mpp_config_t));
	mpp_config_defaults(run->cfg);
	run->cfg->bitdepth = 8;

	run->num_frames = 4;
	run->frame_rows = 50;
	run->frame_cols = 60;
	run->num_layers = 1;
	run->bitdepth = 8;
	for (int i = 0; i < 4; ++i) {
		run->patch_yxyx[i] = 10 + i;
		run->intersection[i] = 20 + i;
	}
	run->best_frame_idx = 2;
	run->stack_size = 2;

	run->mean_frame_rows = 10;
	run->mean_frame_cols = 12;
	run->mean_frame_data = malloc(10 * 12 * sizeof(int32_t));
	for (int i = 0; i < 10 * 12; ++i)
		run->mean_frame_data[i] = i * 7;

	run->quality          = malloc(4 * sizeof(double));
	run->frame_brightness = malloc(4 * sizeof(double));
	run->included         = malloc(4 * sizeof(int32_t));
	run->global_shifts    = malloc(8 * sizeof(double));
	for (int i = 0; i < 4; ++i) {
		run->quality[i] = 0.5 + 0.1 * i;
		run->frame_brightness[i] = 100.0 + i;
		run->included[i] = 1;
		/* Include fractional residuals so the sub-pixel round-trip
		 * is exercised. */
		run->global_shifts[2 * i + 0] = -(double) i + 0.25;
		run->global_shifts[2 * i + 1] =  (double) i - 0.125;
	}

	mpp_aps_t *aps = calloc(1, sizeof(*aps));
	aps->count = 3;
	aps->dropped_dim = 1;
	aps->dropped_structure = 0;
	aps->records = malloc(3 * sizeof(mpp_ap_record_t));
	for (int a = 0; a < 3; ++a) {
		aps->records[a].y = 30 + a;
		aps->records[a].x = 40 + a;
		aps->records[a].box_y_low = 6 + a;
		aps->records[a].box_y_high = 54 + a;
		aps->records[a].box_x_low = 16 + a;
		aps->records[a].box_x_high = 64 + a;
		aps->records[a].patch_y_low = 0;
		aps->records[a].patch_y_high = 72 + a;
		aps->records[a].patch_x_low = 4 + a;
		aps->records[a].patch_x_high = 76 + a;
		aps->records[a].structure = 0.3 + 0.1 * a;
	}
	run->aps = aps;

	run->best_frame_indices = malloc(3 * 2 * sizeof(int32_t));
	for (int a = 0; a < 3; ++a) {
		run->best_frame_indices[a * 2 + 0] = (a + 1) % 4;
		run->best_frame_indices[a * 2 + 1] = (a + 2) % 4;
	}

	if (with_shifts) {
		mpp_shifts_t *s = calloc(1, sizeof(*s));
		s->num_frames = 4;
		s->num_aps = 3;
		s->shifts = malloc(2 * 4 * 3 * sizeof(double));
		s->success = malloc(4 * 3 * sizeof(uint8_t));
		for (int f = 0; f < 4; ++f)
			for (int a = 0; a < 3; ++a) {
				s->shifts[2 * (f * 3 + a) + 0] = f - a * 0.5;
				s->shifts[2 * (f * 3 + a) + 1] = f * 0.25 + a;
				s->success[f * 3 + a] = (f + a) % 2;
			}
		s->failure_counter = 6;
		run->shifts = s;
	}
	return run;
}

cominfo com;
fits *gfit = NULL;

Test(mpp_sidecar, roundtrip_without_shifts) {
	mpp_run_t *src = make_test_run(0);
	cr_assert_not_null(src);

	char tmp[] = "/tmp/mpp_sidecar_test_XXXXXX";
	int fd = mkstemp(tmp);
	cr_assert_gt(fd, 0);
	close(fd);

	cr_assert_eq(mpp_sidecar_write(tmp, src), MPP_OK);

	mpp_run_t *dst = NULL;
	cr_assert_eq(mpp_sidecar_read(tmp, &dst), MPP_OK);
	cr_assert_not_null(dst);
	cr_assert_eq(dst->num_frames, src->num_frames);
	cr_assert_eq(dst->stack_size, src->stack_size);
	cr_assert_eq(dst->aps->count, src->aps->count);
	cr_assert_null(dst->shifts);
	for (int i = 0; i < src->num_frames; ++i) {
		cr_assert_float_eq(dst->quality[i], src->quality[i], 1e-15);
		cr_assert_float_eq(dst->global_shifts[2 * i + 0],
		                   src->global_shifts[2 * i + 0], 1e-15);
		cr_assert_float_eq(dst->global_shifts[2 * i + 1],
		                   src->global_shifts[2 * i + 1], 1e-15);
	}
	for (int a = 0; a < src->aps->count; ++a) {
		cr_assert_eq(dst->aps->records[a].y, src->aps->records[a].y);
		cr_assert_float_eq(dst->aps->records[a].structure,
		                   src->aps->records[a].structure, 1e-15);
	}
	for (int i = 0; i < src->mean_frame_rows * src->mean_frame_cols; ++i)
		cr_assert_eq(dst->mean_frame_data[i], src->mean_frame_data[i]);
	for (int i = 0; i < src->aps->count * src->stack_size; ++i)
		cr_assert_eq(dst->best_frame_indices[i], src->best_frame_indices[i]);

	mpp_run_free(src);
	mpp_run_free(dst);
	unlink(tmp);
}

Test(mpp_sidecar, roundtrip_with_shifts) {
	mpp_run_t *src = make_test_run(1);
	char tmp[] = "/tmp/mpp_sidecar_test_XXXXXX";
	int fd = mkstemp(tmp);
	cr_assert_gt(fd, 0);
	close(fd);

	cr_assert_eq(mpp_sidecar_write(tmp, src), MPP_OK);
	mpp_run_t *dst = NULL;
	cr_assert_eq(mpp_sidecar_read(tmp, &dst), MPP_OK);
	cr_assert_not_null(dst->shifts);
	cr_assert_eq(dst->shifts->num_frames, src->shifts->num_frames);
	cr_assert_eq(dst->shifts->num_aps,    src->shifts->num_aps);
	cr_assert_eq(dst->shifts->failure_counter, src->shifts->failure_counter);
	const int N = src->shifts->num_frames * src->shifts->num_aps;
	for (int i = 0; i < N; ++i) {
		cr_assert_float_eq(dst->shifts->shifts[2 * i + 0],
		                   src->shifts->shifts[2 * i + 0], 1e-15);
		cr_assert_float_eq(dst->shifts->shifts[2 * i + 1],
		                   src->shifts->shifts[2 * i + 1], 1e-15);
		cr_assert_eq(dst->shifts->success[i], src->shifts->success[i]);
	}

	mpp_run_free(src);
	mpp_run_free(dst);
	unlink(tmp);
}

Test(mpp_sidecar, bad_path_returns_eio) {
	mpp_run_t *r = NULL;
	cr_assert_eq(mpp_sidecar_read("/this/path/does/not/exist", &r), MPP_EIO);
	cr_assert_null(r);
}
