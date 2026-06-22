/*
 * mpp_derot_autodetect_disk: the disk seed must lock onto the planet globe,
 * not its surroundings. The hard case is Saturn, whose bright rings extend far
 * past the globe along one axis; a naive bounding box sizes to the rings. The
 * detector takes the minor axis of the bright region (perpendicular to the
 * rings) and recovers the equatorial radius from the body flattening, so it
 * should return ~the globe radius regardless of the ring span.
 */
#include <criterion/criterion.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "core/siril.h"
#include "registration/mpp/mpp_derot_build.h"

/* The linked siril objects reference these globals; the detector itself does
 * not use them, but they must resolve at load time. */
cominfo com;
fits *gfit = NULL;

/* Minimal in-memory USHORT image; the detector only reads rx/ry/type/data. */
static fits make_fits(int w, int h, WORD *buf) {
	fits f;
	memset(&f, 0, sizeof(f));
	f.rx = (unsigned) w; f.ry = (unsigned) h;
	f.naxes[0] = w; f.naxes[1] = h; f.naxes[2] = 1;
	f.type = DATA_USHORT;
	f.data = buf;
	return f;
}

/* A filled oblate globe plus thin bright rings running along x out to 2.3 Req. */
Test(mpp_derot_build, autodetect_globe_not_rings) {
	const int W = 400, H = 300;
	const double cx = 200.0, cy = 150.0, req = 40.0, rpol = 36.0;
	const double sat_flat = 0.09796;
	WORD *buf = calloc((size_t) W * H, sizeof(WORD));
	cr_assert_not_null(buf);
	for (int y = 0; y < H; ++y) {
		for (int x = 0; x < W; ++x) {
			const double nx = (x - cx) / req, ny = (y - cy) / rpol;
			if (nx * nx + ny * ny <= 1.0) {
				buf[y * W + x] = 50000;                      /* globe */
			} else {
				const double dxr = fabs(x - cx), dyr = fabs(y - cy);
				if (dxr >= 1.2 * req && dxr <= 2.3 * req && dyr <= 5.0)
					buf[y * W + x] = 45000;                  /* rings (low tilt) */
			}
		}
	}
	fits f = make_fits(W, H, buf);

	double dcx, dcy, r;
	cr_assert(mpp_derot_autodetect_disk(&f, sat_flat, TRUE, &dcx, &dcy, &r, NULL));
	cr_assert_float_eq(dcx, cx, 8.0, "centre x %.1f", dcx);
	cr_assert_float_eq(dcy, cy, 8.0, "centre y %.1f", dcy);
	/* Must be near the globe equatorial radius, and well below the ring span
	 * (~2.3*req = 92) that a bounding-box detector would return. */
	cr_assert(r > 0.7 * req && r < 1.5 * req,
	          "radius %.1f should be ~globe %.0f, not the ring span ~%.0f",
	          r, req, 2.3 * req);
	free(buf);
}

/* A bare round disk (no rings) is recovered with f = 0. */
Test(mpp_derot_build, autodetect_plain_disk) {
	const int W = 256, H = 256;
	const double cx = 128.0, cy = 120.0, R = 50.0;
	WORD *buf = calloc((size_t) W * H, sizeof(WORD));
	cr_assert_not_null(buf);
	for (int y = 0; y < H; ++y)
		for (int x = 0; x < W; ++x)
			if (hypot(x - cx, y - cy) <= R) buf[y * W + x] = 60000;
	fits f = make_fits(W, H, buf);

	double dcx, dcy, r;
	cr_assert(mpp_derot_autodetect_disk(&f, 0.0, FALSE, &dcx, &dcy, &r, NULL));
	cr_assert_float_eq(dcx, cx, 4.0);
	cr_assert_float_eq(dcy, cy, 4.0);
	cr_assert_float_eq(r, R, 0.2 * R, "radius %.1f vs %.0f", r, R);
	free(buf);
}

/* The shared (multi-sequence) reference epoch is the midpoint of the union of
 * the spans, regardless of order or overlap. */
Test(mpp_derot_build, union_epoch_midpoint) {
	/* three R/G/B spans, captured back to back */
	const double first[3] = { 100.0, 100.2, 100.4 };
	const double last[3]  = { 100.1, 100.3, 100.5 };
	/* union = [100.0, 100.5] -> midpoint 100.25 */
	cr_assert_float_eq(mpp_derot_union_epoch(first, last, 3), 100.25, 1e-12);
	/* order-independent */
	const double f2[3] = { 100.4, 100.0, 100.2 };
	const double l2[3] = { 100.5, 100.1, 100.3 };
	cr_assert_float_eq(mpp_derot_union_epoch(f2, l2, 3), 100.25, 1e-12);
	/* single sequence collapses to its own midpoint */
	const double f1 = 50.0, l1 = 51.0;
	cr_assert_float_eq(mpp_derot_union_epoch(&f1, &l1, 1), 50.5, 1e-12);
	cr_assert_float_eq(mpp_derot_union_epoch(&f1, &l1, 0), 0.0, 1e-12);
}
