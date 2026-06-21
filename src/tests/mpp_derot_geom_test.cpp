/*
 * mpp_derot_geom: oblate-spheroid projection + derotation map invariants.
 *
 * Pure geometry, no fixtures: round-trip consistency of forward/inverse
 * projection, identity at zero rotation, the r_eq * dCM displacement of a
 * central feature, emission-angle behaviour, and limb/off-disk masking.
 */
#include <criterion/criterion.h>
#include <cmath>

#include "registration/mpp/mpp_derot_geom.h"
#include "registration/mpp/mpp_derot.h"
#include "registration/mpp/mpp_derot_sidecar.h"

#define DEG (M_PI / 180.0)

/* forward then inverse recovers (lat,lon) on the visible hemisphere */
Test(mpp_derot_geom, project_unproject_roundtrip) {
	const double B = 12.0 * DEG, CM = 47.0 * DEG, f = 0.06487;
	int tested = 0;
	for (double lat = -80; lat <= 80; lat += 20)
		for (double dlon = -80; dlon <= 80; dlon += 20) {
			double lon = CM + dlon * DEG;
			double u, v, mu;
			if (!derot_project(lat * DEG, lon, B, CM, f, &u, &v, &mu))
				continue;
			double lat2, lon2;
			cr_assert(derot_unproject(u, v, B, CM, f, &lat2, &lon2));
			cr_assert_float_eq(lat2, lat * DEG, 1e-9, "lat %g", lat);
			double dl = lon2 - lon;
			while (dl >  M_PI) dl -= 2 * M_PI;
			while (dl < -M_PI) dl += 2 * M_PI;
			cr_assert_float_eq(dl, 0.0, 1e-9, "lon at lat=%g dlon=%g", lat, dlon);
			tested++;
		}
	cr_assert_gt(tested, 20);
}

/* inverse then forward recovers (u,v) inside the disk */
Test(mpp_derot_geom, unproject_project_roundtrip) {
	const double B = -5.0 * DEG, CM = 130.0 * DEG, f = 0.05;
	for (double u = -0.9; u <= 0.9; u += 0.3)
		for (double v = -0.9; v <= 0.9; v += 0.3) {
			if (u * u + v * v > 0.85) continue;
			double lat, lon;
			if (!derot_unproject(u, v, B, CM, f, &lat, &lon)) continue;
			double u2, v2, mu;
			cr_assert(derot_project(lat, lon, B, CM, f, &u2, &v2, &mu));
			cr_assert_float_eq(u2, u, 1e-9);
			cr_assert_float_eq(v2, v, 1e-9);
		}
}

/* a sphere's sub-observer point projects to the centre with mu=1 */
Test(mpp_derot_geom, subobserver_centre) {
	double u, v, mu;
	cr_assert(derot_project(20 * DEG, 75 * DEG, 20 * DEG, 75 * DEG, 0.0, &u, &v, &mu));
	cr_assert_float_eq(u, 0.0, 1e-12);
	cr_assert_float_eq(v, 0.0, 1e-12);
	cr_assert_float_eq(mu, 1.0, 1e-12);
}

/* zero rotation -> identity map; disk pixels valid, off-disk masked */
Test(mpp_derot_geom, identity_map) {
	derot_diskfit_t disk = { 100.0, 100.0, 80.0, 0.06487, +1.0 };
	derot_geom_t g = { 8 * DEG, 33 * DEG, 0.0 };
	cv::Mat mx, my, valid, mu;
	derot_build_map(200, 200, disk, g, g, mx, my, valid, mu);
	/* epoch == frame: the globe maps to itself and everything outside passes
	 * through, so the whole canvas is the identity and fully usable. */
	for (int y = 0; y < 200; ++y)
		for (int x = 0; x < 200; ++x) {
			cr_assert_eq(valid.at<uint8_t>(y, x), 1);
			cr_assert_float_eq(mx.at<float>(y, x), (float) x, 1e-3f);
			cr_assert_float_eq(my.at<float>(y, x), (float) y, 1e-3f);
			cr_assert_float_eq(mu.at<float>(y, x), 1.0f, 1e-6f);
		}
}

/* rotating CM by a small angle shifts the central feature by r_eq * dCM */
Test(mpp_derot_geom, central_displacement) {
	derot_diskfit_t disk = { 256.0, 256.0, 200.0, 0.0, +1.0 };
	double dCM = 0.5 * DEG;
	derot_geom_t e = { 0.0, 0.0,   0.0 };
	derot_geom_t fr = { 0.0, dCM,   0.0 };
	cv::Mat mx, my, valid, mu;
	derot_build_map(512, 512, disk, e, fr, mx, my, valid, mu);
	cr_assert(valid.at<uint8_t>(256, 256));
	/* epoch sub-observer point (centre) is at body lon 0; at frame CM=dCM it
	 * sits at u = sin(-dCM) -> source x = cx - r_eq*sin(dCM) */
	double expect_x = 256.0 - 200.0 * sin(dCM);
	cr_assert_float_eq(mx.at<float>(256, 256), (float) expect_x, 0.05f,
	                   "mapx=%.4f expect=%.4f", mx.at<float>(256, 256), expect_x);
	cr_assert_float_eq(my.at<float>(256, 256), 256.0f, 0.05f);
}

/* emission cosine: ~1 at centre, small near the limb */
Test(mpp_derot_geom, emission_angle) {
	double u, v, mu_c, mu_l;
	derot_project(0, 0, 0, 0, 0.0, &u, &v, &mu_c);          /* centre */
	cr_assert_float_eq(mu_c, 1.0, 1e-12);
	cr_assert(derot_project(0, 85 * DEG, 0, 0, 0.0, &u, &v, &mu_l)); /* near limb */
	cr_assert_lt(mu_l, 0.10f);
	cr_assert_gt(mu_l, 0.0f);
}

/* The .derot bridge: when a frame's geometry equals the epoch geometry, the
 * generated map is the identity over the disk and masks the background. Also
 * checks the canvas-origin/scale transform of the disk fit. */
Test(mpp_derot_geom, bridge_identity_at_epoch) {
	mpp_derot_t *d = mpp_derot_alloc(1);
	cr_assert_not_null(d);
	d->cx = 100.0; d->cy = 100.0; d->r_eq = 70.0;
	d->flattening = 0.06487; d->parity = 1.0; d->pole_angle_epoch = 0.0;
	d->epoch_sub_obs_lat = 5.0; d->epoch_cm = 40.0; d->epoch_pole_pa = 12.0;
	d->sub_obs_lat[0] = 5.0; d->cm[0] = 40.0; d->pole_pa[0] = 12.0;  /* == epoch */

	cv::Mat mx, my, mu;
	mpp_derot_frame_map(d, 0, 200, 200, /*org*/0.0, 0.0, /*scale*/1.0, mx, my, mu);
	int nvalid = 0;
	for (int y = 0; y < 200; ++y)
		for (int x = 0; x < 200; ++x)
			if (mu.at<float>(y, x) > 0.0f) {
				cr_assert_float_eq(mx.at<float>(y, x), (float) x, 1e-3f);
				cr_assert_float_eq(my.at<float>(y, x), (float) y, 1e-3f);
				nvalid++;
			}
	cr_assert_gt(nvalid, 5000);
	cr_assert_eq(mu.at<float>(2, 2), 1.0f);   /* corner is off-globe: pass-through */

	/* canvas origin + scale: a 2x canvas offset by (10,10) puts the disk
	 * centre at ((100-10)*2, (100-10)*2) = (180,180). */
	mpp_derot_frame_map(d, 0, 400, 400, 10.0, 10.0, 2.0, mx, my, mu);
	cr_assert_gt(mu.at<float>(180, 180), 0.9f);   /* near disk centre, mu~1 */
	mpp_derot_free(d);
}

/* The map only derotates the globe: pixels outside the fitted disk pass
 * through as identity (rings/sky survive), and a small CM rotation displaces a
 * central globe feature while masking only the sliver that rotated off-limb. */
Test(mpp_derot_geom, off_globe_passthrough_and_limb_mask) {
	derot_diskfit_t disk = { 256.0, 256.0, 120.0, 0.0648, +1.0 };
	derot_geom_t e = { 0.0, 0.0,        0.0 };
	derot_geom_t fr = { 0.0, 6.0 * DEG, 0.0 };   /* 6 deg rotation */
	cv::Mat mx, my, valid, mu;
	derot_build_map(512, 512, disk, e, fr, mx, my, valid, mu);

	/* corner is far outside the globe -> identity pass-through, usable */
	cr_assert_eq(valid.at<uint8_t>(8, 8), 1);
	cr_assert_float_eq(mx.at<float>(8, 8), 8.0f, 1e-3f);
	cr_assert_float_eq(my.at<float>(8, 8), 8.0f, 1e-3f);
	cr_assert_float_eq(mu.at<float>(8, 8), 1.0f, 1e-6f);

	/* centre globe feature is displaced (not identity) and masked nowhere */
	cr_assert(valid.at<uint8_t>(256, 256));
	cr_assert_neq(mx.at<float>(256, 256), 256.0f);

	/* somewhere along the leading limb a globe point has rotated out of view */
	int masked = 0;
	for (int x = 0; x < 512; ++x)
		for (int y = 200; y < 312; ++y)
			if (valid.at<uint8_t>(y, x) == 0) masked++;
	cr_assert_gt(masked, 0, "expected a masked off-limb sliver");
}

/* a feature that rotates onto the far side becomes invalid in the map */
Test(mpp_derot_geom, rotate_off_limb) {
	/* epoch: feature visible near the trailing limb (lon = CM + 80 deg) */
	double u, v, mu;
	cr_assert(derot_project(0, 80 * DEG, 0, 0, 0.0, &u, &v, &mu));
	double lat, lon;
	cr_assert(derot_unproject(u, v, 0, 0, 0.0, &lat, &lon));
	/* move CM to -40 deg -> feature separation lon-CM = 120 deg -> far side */
	double u2, v2, mu2;
	cr_assert_not(derot_project(lat, lon, 0, -40 * DEG, 0.0, &u2, &v2, &mu2));
}
