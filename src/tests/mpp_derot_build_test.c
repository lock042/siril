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

/* Alt-az field rotation folded into the plan: relative to a plan built with
 * field rotation off, pole_pa[i] must differ by exactly −(q_i − q_epoch),
 * with q evaluated from the same site/time/target the builder used. */
Test(mpp_derot_build, altaz_field_rotation_fold) {
	const int N = 5;
	double jd[5];
	const double jd0 = 2461212.5;              /* 2026-06-21T00:00:00 UTC */
	for (int i = 0; i < N; i++)
		jd[i] = jd0 + i * (2.0 / 1440.0);      /* 2-minute cadence */
	const double epoch = 0.5 * (jd[0] + jd[N - 1]);
	const double lat = 52.0, lon = -2.0;

	mpp_derot_t *base = mpp_derot_build(PLANET_JUPITER, 2, epoch, jd, N,
	                                    600, 600, 300, 300, 180, 0.0, 1.0,
	                                    lat, lon, NAN, MPP_FIELD_ROT_NONE);
	mpp_derot_t *rot  = mpp_derot_build(PLANET_JUPITER, 2, epoch, jd, N,
	                                    600, 600, 300, 300, 180, 0.0, 1.0,
	                                    lat, lon, NAN, MPP_FIELD_ROT_ALTAZ);
	cr_assert_not_null(base);
	cr_assert_not_null(rot);

	planet_geom_t ge;
	cr_assert_eq(planet_ephemeris(PLANET_JUPITER, epoch, lat, lon, NAN, &ge), 0);
	const double q0 = planet_parallactic_angle(epoch, ge.ra, ge.dec, lat, lon);

	double span = 0.0;
	for (int i = 0; i < N; i++) {
		planet_geom_t g;
		cr_assert_eq(planet_ephemeris(PLANET_JUPITER, jd[i], lat, lon, NAN, &g), 0);
		const double qi = planet_parallactic_angle(jd[i], g.ra, g.dec, lat, lon);
		const double expect = base->pole_pa[i] - (qi - q0);
		cr_assert_float_eq(rot->pole_pa[i], expect, 1e-9,
		                   "frame %d: pole_pa %.9f expected %.9f",
		                   i, rot->pole_pa[i], expect);
		span = qi - q0;
	}
	/* And the fold is non-trivial: 8 minutes of hour angle moves q. */
	cr_assert(fabs(span) > 1e-3, "q span %.6f unexpectedly tiny", span);

	/* Requiring the site: ALTAZ without lat/lon must fail cleanly. */
	mpp_derot_t *bad = mpp_derot_build(PLANET_JUPITER, 2, epoch, jd, N,
	                                   600, 600, 300, 300, 180, 0.0, 1.0,
	                                   NAN, NAN, NAN, MPP_FIELD_ROT_ALTAZ);
	cr_assert_null(bad);

	mpp_derot_free(base);
	mpp_derot_free(rot);
}

/* Rotation-only (field-rotation) plan: body sentinel, zero planetary
 * geometry, and per-frame pole angles equal to -(q_i - q_epoch) computed
 * directly for the same target/site/times. */
Test(mpp_derot_build, field_rotation_only_plan) {
	const int N = 4;
	double jd[4];
	const double jd0 = 2461212.5;
	for (int i = 0; i < N; i++)
		jd[i] = jd0 + i * (3.0 / 1440.0);
	const double lat = -37.81, lon = 144.96;

	mpp_derot_t *d = mpp_derot_build_field_rotation(EPHEM_TARGET_SUN, jd, N,
	                                                480, 640, lat, lon);
	cr_assert_not_null(d);
	cr_assert_eq(d->body, MPP_DEROT_BODY_FIELDROT);
	cr_assert_eq(d->rot_system, EPHEM_TARGET_SUN);
	cr_assert_float_eq(d->flattening, 0.0, 1e-12);
	cr_assert_float_eq(d->parity, 1.0, 1e-12);

	const double epoch = 0.5 * (jd[0] + jd[N - 1]);
	double ra0, dec0;
	cr_assert_eq(ephem_target_radec(EPHEM_TARGET_SUN, epoch, &ra0, &dec0), 0);
	const double q0 = planet_parallactic_angle(epoch, ra0, dec0, lat, lon);
	for (int i = 0; i < N; i++) {
		double ra, dec;
		cr_assert_eq(ephem_target_radec(EPHEM_TARGET_SUN, jd[i], &ra, &dec), 0);
		const double qi = planet_parallactic_angle(jd[i], ra, dec, lat, lon);
		cr_assert_float_eq(d->pole_pa[i], -(qi - q0), 1e-9,
		                   "frame %d: %.9f vs %.9f", i, d->pole_pa[i], -(qi - q0));
		cr_assert_float_eq(d->sub_obs_lat[i], 0.0, 1e-12);
		cr_assert_float_eq(d->cm[i], 0.0, 1e-12);
	}
	mpp_derot_free(d);

	/* Missing site fails cleanly. */
	cr_assert_null(mpp_derot_build_field_rotation(EPHEM_TARGET_SUN, jd, N,
	                                              480, 640, NAN, lon));
}

/* Registration-side sync against an existing PLANETARY plan: enabling folds
 * the alt-az drift into the plan's pole angles (derotation + field rotation
 * compose), re-enabling with the same site is a no-op, and disabling
 * restores the original angles exactly. */
Test(mpp_derot_build, sync_folds_into_planetary_plan) {
	const int N = 4;
	double jd[4];
	for (int i = 0; i < N; i++)
		jd[i] = 2461212.5 + i * (3.0 / 1440.0);
	const double epoch = 0.5 * (jd[0] + jd[N - 1]);
	const double lat = -37.81, lon = 144.96;

	mpp_derot_t *plan = mpp_derot_build(PLANET_JUPITER, 2, epoch, jd, N,
	                                    600, 600, 300, 300, 180, 10.0, 1.0,
	                                    NAN, NAN, NAN, MPP_FIELD_ROT_NONE);
	cr_assert_not_null(plan);
	double orig[4];
	for (int i = 0; i < N; i++) orig[i] = plan->pole_pa[i];

	gchar *base = g_build_filename(g_get_tmp_dir(), "siril_fold_test", NULL);
	gchar *path = g_strdup_printf("%s.derot", base);
	cr_assert_eq(mpp_derot_write(path, plan), MPP_OK);
	mpp_derot_free(plan);

	sequence seq;
	memset(&seq, 0, sizeof(seq));
	seq.seqname = base;
	seq.number = N;

	/* enable: fold */
	cr_assert_eq(mpp_derot_sync_field_rotation_plan(&seq, TRUE,
	             EPHEM_TARGET_SUN /* ignored for planetary plans */,
	             lat, lon), MPP_OK);
	mpp_derot_t *folded = NULL;
	cr_assert_eq(mpp_derot_read(path, &folded), MPP_OK);
	cr_assert(folded->flags & MPP_DEROT_FLAG_FIELDROT);
	cr_assert_eq(folded->body, PLANET_JUPITER);
	cr_assert_float_eq(folded->obs_lat, lat, 1e-12);
	double ra, dec;
	cr_assert_eq(ephem_target_radec(EPHEM_TARGET_JUPITER, epoch, &ra, &dec), 0);
	const double q0 = planet_parallactic_angle(epoch, ra, dec, lat, lon);
	for (int i = 0; i < N; i++) {
		cr_assert_eq(ephem_target_radec(EPHEM_TARGET_JUPITER, jd[i], &ra, &dec), 0);
		const double qi = planet_parallactic_angle(jd[i], ra, dec, lat, lon);
		cr_assert_float_eq(folded->pole_pa[i], orig[i] - (qi - q0), 1e-9,
		                   "frame %d fold mismatch", i);
	}
	mpp_derot_free(folded);

	/* enable again, same site: idempotent */
	cr_assert_eq(mpp_derot_sync_field_rotation_plan(&seq, TRUE,
	             EPHEM_TARGET_SUN, lat, lon), MPP_OK);
	mpp_derot_t *again = NULL;
	cr_assert_eq(mpp_derot_read(path, &again), MPP_OK);
	cr_assert(again->flags & MPP_DEROT_FLAG_FIELDROT);
	mpp_derot_free(again);

	/* disable: unfold restores the original angles, keeps the plan */
	cr_assert_eq(mpp_derot_sync_field_rotation_plan(&seq, FALSE, 0,
	             NAN, NAN), MPP_OK);
	mpp_derot_t *restored = NULL;
	cr_assert_eq(mpp_derot_read(path, &restored), MPP_OK);
	cr_assert_not(restored->flags & MPP_DEROT_FLAG_FIELDROT);
	for (int i = 0; i < N; i++)
		cr_assert_float_eq(restored->pole_pa[i], orig[i], 1e-9,
		                   "frame %d unfold mismatch", i);
	mpp_derot_free(restored);

	g_unlink(path);
	g_free(path);
	g_free(base);
}

/* Timestamp-order classifier: a permuted-but-coherent capture timeline
 * (quality-sorted output) is benign; garbage trailers are not. */
Test(mpp_derot_build, timestamp_order_classifier) {
	const int N = 1000;
	double jd[1000];
	const double jd0 = 2461212.5, dt = 0.02 / 86400.0;   /* 50 fps */
	/* Deterministic large-scale permutation of a perfect cadence (stride
	 * walk over the whole span — like a quality sort). */
	for (int i = 0; i < N; i++)
		jd[i] = jd0 + ((i * 617) % N) * dt;
	cr_assert(mpp_derot_report_timestamp_order(jd, N, N / 2),
	          "coherent permuted timeline must classify as benign");

	/* Broken: half the stamps collapsed to one value (dup run). */
	for (int i = 0; i < N; i++)
		jd[i] = (i % 2) ? jd0 : jd0 + ((i * 617) % N) * dt;
	cr_assert_not(mpp_derot_report_timestamp_order(jd, N, N / 2),
	              "duplicate-run trailer must classify as broken");

	/* Broken: one stamp light-years away dominates the span. */
	for (int i = 0; i < N; i++)
		jd[i] = jd0 + ((i * 617) % N) * dt;
	jd[500] = jd0 + 300.0;   /* +300 days */
	cr_assert_not(mpp_derot_report_timestamp_order(jd, N, N / 2),
	              "dominating gap must classify as broken");
}
