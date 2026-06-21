/*
 * planet_ephem: built-in analytic ephemeris validation.
 *
 * Offline assertions use reference values captured from the IMCCE Miriade web
 * service (ephemcc + ephemph, geocentre, default IAU 2015 physical model) for
 * 2026-06-21T00:00:00 UTC (JD 2461212.5). The network-gated test re-fetches
 * Miriade live and compares; enable it with SIRIL_TEST_EPHEM_NETWORK=1.
 *
 * Conventions: RA/Dec/distance are astrometric J2000 and match Miriade to
 * < 0.001 deg. Sub-observer latitude B matches to < 0.01 deg. The pole position
 * angle P and the absolute central-meridian longitude carry a constant frame
 * offset (J2000 vs of-date, ~0.1 deg) and time-scale offset respectively; both
 * cancel in the per-frame *differences* derotation actually uses, so we assert
 * the CM rotation *rate* tightly and the absolute P loosely.
 */
#include <criterion/criterion.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "algos/planet_ephem.h"

#define JD_REF 2461212.5   /* 2026-06-21T00:00:00 UTC */

Test(planet_ephem, jupiter_geometry) {
	planet_geom_t g;
	cr_assert_eq(planet_ephemeris(PLANET_JUPITER, JD_REF, NAN, NAN, NAN, &g), 0);
	cr_assert_float_eq(g.ra,          119.85123, 0.002, "RA=%.5f", g.ra);
	cr_assert_float_eq(g.dec,          21.04249, 0.002, "Dec=%.5f", g.dec);
	cr_assert_float_eq(g.dist_au,    6.14454624, 5e-5, "dist=%.7f", g.dist_au);
	cr_assert_float_eq(g.sub_obs_lat,      1.00, 0.05, "B=%.3f", g.sub_obs_lat);
	cr_assert_float_eq(g.pole_pa,         13.25, 0.30, "P=%.3f", g.pole_pa);
	cr_assert_float_eq(g.ang_diam_eq, 2 * 16.0423, 0.02, "diam=%.4f", g.ang_diam_eq);
	cr_assert_float_eq(g.flattening,    0.06487, 1e-4);
	for (int s = 0; s < 3; ++s)
		cr_assert(g.cm[s] >= 0.0 && g.cm[s] < 360.0, "cm[%d]=%.3f", s, g.cm[s]);
}

Test(planet_ephem, mars_geometry) {
	planet_geom_t g;
	cr_assert_eq(planet_ephemeris(PLANET_MARS, JD_REF, NAN, NAN, NAN, &g), 0);
	cr_assert_float_eq(g.ra,           51.76147, 0.002, "RA=%.5f", g.ra);
	cr_assert_float_eq(g.dec,          18.37350, 0.002, "Dec=%.5f", g.dec);
	cr_assert_float_eq(g.dist_au,    2.13493751, 5e-5, "dist=%.7f", g.dist_au);
	cr_assert_float_eq(g.sub_obs_lat,    -12.14, 0.05, "B=%.3f", g.sub_obs_lat);
	cr_assert_float_eq(g.pole_pa,        322.11, 0.30, "P=%.3f", g.pole_pa);
}

Test(planet_ephem, saturn_position) {
	planet_geom_t g;
	cr_assert_eq(planet_ephemeris(PLANET_SATURN, JD_REF, NAN, NAN, NAN, &g), 0);
	cr_assert_float_eq(g.ra,           13.15678, 0.002, "RA=%.5f", g.ra);
	cr_assert_float_eq(g.dec,           3.09739, 0.002, "Dec=%.5f", g.dec);
	cr_assert_float_eq(g.dist_au,    9.66018417, 5e-5, "dist=%.7f", g.dist_au);
}

/* The rotation amount over an interval is what derotation applies. Over 1 h it
 * must equal the system's sidereal rate (minus the tiny synodic correction). */
Test(planet_ephem, cm_rotation_rate) {
	planet_geom_t a, b;
	const double dt_h = 1.0 / 24.0;
	planet_ephemeris(PLANET_JUPITER, JD_REF,        NAN, NAN, NAN, &a);
	planet_ephemeris(PLANET_JUPITER, JD_REF + dt_h, NAN, NAN, NAN, &b);
	const double rate[3] = { 877.900 / 24.0, 870.270 / 24.0, 870.536 / 24.0 };
	for (int s = 0; s < 3; ++s) {
		double d = b.cm[s] - a.cm[s];
		while (d <= -180.0) d += 360.0;
		while (d >   180.0) d -= 360.0;
		/* tolerance covers the ~0.01 deg/h synodic correction */
		cr_assert_float_eq(fabs(d), rate[s], 0.05,
		                   "system %d: |dCM|=%.4f deg/h, expected %.4f", s, fabs(d), rate[s]);
	}
}

Test(planet_ephem, delta_t_contemporary) {
	double dt = planet_ephem_delta_t(JD_REF);
	cr_assert(dt > 67.0 && dt < 71.0, "delta-T(2026)=%.2f s", dt);
}

/* Live cross-check against Miriade. Gated: shells out to curl. */
Test(planet_ephem, miriade_live) {
	if (!getenv("SIRIL_TEST_EPHEM_NETWORK")) {
		cr_log_info("set SIRIL_TEST_EPHEM_NETWORK=1 to enable the live Miriade check");
		return;
	}
	const char *url =
	    "https://ssp.imcce.fr/webservices/miriade/api/ephemcc.php?"
	    "-name=p:jupiter&-mime=text/csv&-ep=2026-06-21T00:00:00&-observer=@500&-tcoor=5";
	char cmd[1024];
	snprintf(cmd, sizeof(cmd), "curl -sS --max-time 25 '%s' 2>/dev/null", url);
	FILE *p = popen(cmd, "r");
	cr_assert_not_null(p, "popen(curl) failed");
	char line[512]; double ra = NAN, dec = NAN, dist = NAN;
	while (fgets(line, sizeof(line), p)) {
		if (strncmp(line, "2026", 4) == 0) {
			char *t = strtok(line, ",");          /* date */
			char *sra = strtok(NULL, ",");
			char *sdec = strtok(NULL, ",");
			char *sd = strtok(NULL, ",");
			if (sra && sdec && sd) { ra = atof(sra); dec = atof(sdec); dist = atof(sd); }
			(void) t;
		}
	}
	pclose(p);
	cr_assert(!isnan(ra), "no Miriade response (network?)");
	planet_geom_t g;
	planet_ephemeris(PLANET_JUPITER, JD_REF, NAN, NAN, NAN, &g);
	cr_assert_float_eq(g.ra,   ra,   0.005, "RA mine=%.5f miriade=%.5f", g.ra, ra);
	cr_assert_float_eq(g.dec,  dec,  0.005, "Dec mine=%.5f miriade=%.5f", g.dec, dec);
	cr_assert_float_eq(g.dist_au, dist, 1e-4, "dist mine=%.7f miriade=%.7f", g.dist_au, dist);
}
