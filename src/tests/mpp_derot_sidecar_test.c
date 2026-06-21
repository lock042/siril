/*
 * mpp_derot_sidecar binary I/O round-trip.
 *
 * Builds a synthetic mpp_derot_t, writes it, reads it back, and checks every
 * persisted field. Also verifies the strict magic/version check and that a
 * truncated file is rejected.
 */
#include <criterion/criterion.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "registration/mpp/mpp_derot_sidecar.h"

static mpp_derot_t *make(int n) {
	mpp_derot_t *d = mpp_derot_alloc(n);
	d->body = 0; d->rot_system = 1; d->ephem_version = 3; d->num_frames = n;
	d->frame_rows = 480; d->frame_cols = 640;
	d->epoch_jd = 2461212.5;
	d->obs_lat = 48.85; d->obs_lon = 2.35; d->obs_elev = 35.0;
	d->cx = 511.5; d->cy = 383.25; d->r_eq = 201.7; d->flattening = 0.06487;
	d->parity = -1.0; d->pole_angle_epoch = 0.1234;
	d->epoch_sub_obs_lat = 1.0; d->epoch_cm = 164.9; d->epoch_pole_pa = 13.25;
	for (int i = 0; i < n; ++i) {
		d->jd[i] = 2461212.5 + i / 86400.0;
		d->sub_obs_lat[i] = 1.0 + i * 1e-5;
		d->cm[i] = fmod(164.9 + i * 0.006, 360.0);
		d->pole_pa[i] = 13.25 - i * 1e-4;
	}
	return d;
}

Test(mpp_derot_sidecar, roundtrip) {
	const int n = 2000;
	mpp_derot_t *src = make(n);
	char tmp[] = "/tmp/derot_test_XXXXXX";
	int fd = mkstemp(tmp); cr_assert_geq(fd, 0); close(fd);

	cr_assert_eq(mpp_derot_write(tmp, src), MPP_OK);
	mpp_derot_t *dst = NULL;
	cr_assert_eq(mpp_derot_read(tmp, &dst), MPP_OK);
	cr_assert_not_null(dst);

	cr_assert_eq(dst->body, src->body);
	cr_assert_eq(dst->rot_system, src->rot_system);
	cr_assert_eq(dst->ephem_version, src->ephem_version);
	cr_assert_eq(dst->num_frames, src->num_frames);
	cr_assert_eq(dst->frame_rows, src->frame_rows);
	cr_assert_eq(dst->frame_cols, src->frame_cols);
	cr_assert_float_eq(dst->epoch_jd, src->epoch_jd, 1e-9);
	cr_assert_float_eq(dst->obs_lat, src->obs_lat, 1e-12);
	cr_assert_float_eq(dst->obs_lon, src->obs_lon, 1e-12);
	cr_assert_float_eq(dst->obs_elev, src->obs_elev, 1e-12);
	cr_assert_float_eq(dst->cx, src->cx, 1e-12);
	cr_assert_float_eq(dst->cy, src->cy, 1e-12);
	cr_assert_float_eq(dst->r_eq, src->r_eq, 1e-12);
	cr_assert_float_eq(dst->flattening, src->flattening, 1e-12);
	cr_assert_float_eq(dst->parity, src->parity, 1e-12);
	cr_assert_float_eq(dst->pole_angle_epoch, src->pole_angle_epoch, 1e-12);
	cr_assert_float_eq(dst->epoch_sub_obs_lat, src->epoch_sub_obs_lat, 1e-12);
	cr_assert_float_eq(dst->epoch_cm, src->epoch_cm, 1e-12);
	cr_assert_float_eq(dst->epoch_pole_pa, src->epoch_pole_pa, 1e-12);
	for (int i = 0; i < n; ++i) {
		cr_assert_float_eq(dst->jd[i], src->jd[i], 1e-9);
		cr_assert_float_eq(dst->sub_obs_lat[i], src->sub_obs_lat[i], 1e-12);
		cr_assert_float_eq(dst->cm[i], src->cm[i], 1e-12);
		cr_assert_float_eq(dst->pole_pa[i], src->pole_pa[i], 1e-12);
	}
	mpp_derot_free(src); mpp_derot_free(dst);
	unlink(tmp);
}

/* The fingerprint must change for ANY change to the plan — this is what lets
 * Stage C reject a .derot that was edited (disk fit, epoch, per-frame geometry)
 * without re-registering. */
Test(mpp_derot_sidecar, fingerprint_detects_changes) {
	mpp_derot_t *a = make(100);
	const uint64_t fa = mpp_derot_fingerprint(a);
	cr_assert_neq(fa, 0u);
	cr_assert_eq(mpp_derot_fingerprint(a), fa, "fingerprint must be stable");

	mpp_derot_t *disk = make(100); disk->cx += 0.5;
	cr_assert_neq(mpp_derot_fingerprint(disk), fa, "disk-fit change undetected");

	mpp_derot_t *epoch = make(100); epoch->epoch_jd += 1e-5;
	cr_assert_neq(mpp_derot_fingerprint(epoch), fa, "epoch change undetected");

	mpp_derot_t *frame = make(100); frame->cm[50] += 1e-6;
	cr_assert_neq(mpp_derot_fingerprint(frame), fa, "per-frame change undetected");

	mpp_derot_t *sys = make(100); sys->rot_system = 2;
	cr_assert_neq(mpp_derot_fingerprint(sys), fa, "system change undetected");

	mpp_derot_free(a); mpp_derot_free(disk); mpp_derot_free(epoch);
	mpp_derot_free(frame); mpp_derot_free(sys);
}

Test(mpp_derot_sidecar, zero_frames) {
	mpp_derot_t *src = make(0);
	char tmp[] = "/tmp/derot_test0_XXXXXX";
	int fd = mkstemp(tmp); cr_assert_geq(fd, 0); close(fd);
	cr_assert_eq(mpp_derot_write(tmp, src), MPP_OK);
	mpp_derot_t *dst = NULL;
	cr_assert_eq(mpp_derot_read(tmp, &dst), MPP_OK);
	cr_assert_eq(dst->num_frames, 0);
	mpp_derot_free(src); mpp_derot_free(dst);
	unlink(tmp);
}

Test(mpp_derot_sidecar, bad_magic) {
	char tmp[] = "/tmp/derot_bad_XXXXXX";
	int fd = mkstemp(tmp); cr_assert_geq(fd, 0);
	FILE *f = fdopen(fd, "wb");
	fwrite("NOTADEROT....", 1, 13, f);
	fclose(f);
	mpp_derot_t *dst = NULL;
	cr_assert_eq(mpp_derot_read(tmp, &dst), MPP_EIO);
	cr_assert_null(dst);
	unlink(tmp);
}

Test(mpp_derot_sidecar, truncated) {
	mpp_derot_t *src = make(500);
	char tmp[] = "/tmp/derot_trunc_XXXXXX";
	int fd = mkstemp(tmp); cr_assert_geq(fd, 0); close(fd);
	cr_assert_eq(mpp_derot_write(tmp, src), MPP_OK);
	/* chop the file in half */
	FILE *f = fopen(tmp, "rb"); fseek(f, 0, SEEK_END); long sz = ftell(f); fclose(f);
	cr_assert_eq(truncate(tmp, sz / 2), 0);
	mpp_derot_t *dst = NULL;
	cr_assert_eq(mpp_derot_read(tmp, &dst), MPP_EIO);
	cr_assert_null(dst);
	mpp_derot_free(src);
	unlink(tmp);
}
