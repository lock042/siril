/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Unit tests for the colour-space transforms in algos/colors.h / algos/colors.c.
 *
 * Strategy: the *branchless* float cores that were moved to colors.h as static
 * inline are compared against independent, locally embedded "oracle" copies of
 * the original branchy double implementations. This guards two things at once:
 *   - numerical agreement with the trusted reference to float precision, and
 *   - the exclusive-mask hue dispatch (so a two-channel tie for the maximum,
 *     e.g. a clipped (1,1,0) highlight, no longer reports a hue 60 degrees off).
 *
 * New transform families append their own oracle + Test() block below.
 */

#include <criterion/criterion.h>
#include <math.h>
#include <stdlib.h>

#include "core/siril.h"
#include "algos/colors.h"

/* ============================ test utilities ============================= */

#define N_FUZZ 500000

static unsigned long s_rng = 0x2545F4914F6CDD1DUL;
static void rng_reset(void) { s_rng = 0x2545F4914F6CDD1DUL; }
static float frand(void) {              /* deterministic xorshift in [0,1] */
	s_rng ^= s_rng << 13; s_rng ^= s_rng >> 7; s_rng ^= s_rng << 17;
	return (float)((s_rng >> 11) * (1.0 / 9007199254740992.0));
}
static float fquant(void) {             /* snap to force exact ties / clipping */
	static const float q[] = {0.f, 0.25f, 0.5f, 0.75f, 1.f};
	s_rng = s_rng * 6364136223846793005UL + 1;
	return q[(s_rng >> 60) % 5];
}
static double hue_err(double a, double b) {  /* circular distance on [0,1) */
	double d = fabs(a - b);
	return (d > 0.5) ? 1.0 - d : d;
}

/* ============================ HSL <-> RGB =============================== */
/* Oracle: original branchy double implementations (max/min as ternaries). */

static void rgb_to_hsl_oracle(double r, double g, double b,
                              double *h, double *s, double *l) {
	double v = (r > g) ? r : g; v = (v > b) ? v : b;
	double m = (r < g) ? r : g; m = (m < b) ? m : b;
	double vm, r2, g2, b2;
	*h = 0.0; *s = 0.0;
	if ((*l = (m + v) / 2.0) <= 0.0) { *l = 0.0; return; }
	if ((*s = vm = v - m) > 0.0)
		*s /= (*l <= 0.5) ? (v + m) : (2.0 - v - m);
	else return;
	r2 = (v - r) / vm; g2 = (v - g) / vm; b2 = (v - b) / vm;
	if (r == v)      *h = (g == m ? 5.0 + b2 : 1.0 - g2);
	else if (g == v) *h = (b == m ? 1.0 + r2 : 3.0 - b2);
	else             *h = (r == m ? 3.0 + g2 : 5.0 - r2);
	*h /= 6;
}

static void hsl_to_rgb_oracle(double h, double sl, double l,
                              double *r, double *g, double *b) {
	double v;
	if (h >= 1.0) h -= 1.0;
	v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);
	if (v <= 0) { *r = *g = *b = 0.0; return; }
	double m = l + l - v;
	double sv = (v - m) / v;
	h *= 6.0;
	int sextant = (int)h;
	double fract = h - sextant;
	double vsf = v * sv * fract;
	double mid1 = m + vsf, mid2 = v - vsf;
	switch (sextant) {
		case 0: *r=v;    *g=mid1; *b=m;    break;
		case 1: *r=mid2; *g=v;    *b=m;    break;
		case 2: *r=m;    *g=v;    *b=mid1; break;
		case 3: *r=m;    *g=mid2; *b=v;    break;
		case 4: *r=mid1; *g=m;    *b=v;    break;
		default:*r=v;    *g=m;    *b=mid2; break;
	}
}

Test(color_hsl, forward_matches_double_oracle) {
	rng_reset();
	double maxh = 0, maxs = 0, maxl = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r, g, b;
		if (i & 1) { r=fquant(); g=fquant(); b=fquant(); }
		else       { r=frand();  g=frand();  b=frand();  }
		double hr, sr, lr; rgb_to_hsl_oracle(r, g, b, &hr, &sr, &lr);
		float hn, sn, ln;  rgb_to_hslf(r, g, b, &hn, &sn, &ln);
		double eh = hue_err(hr, hn);
		if (eh > maxh) maxh = eh;
		if (fabs(sr-sn) > maxs) maxs = fabs(sr-sn);
		if (fabs(lr-ln) > maxl) maxl = fabs(lr-ln);
	}
	cr_expect(maxh < 1e-4, "max hue err %.3e too large", maxh);
	cr_expect(maxs < 1e-4, "max sat err %.3e too large", maxs);
	cr_expect(maxl < 1e-5, "max lum err %.3e too large", maxl);
}

Test(color_hsl, inverse_matches_double_oracle) {
	rng_reset();
	double maxr = 0, maxg = 0, maxb = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float h=frand(), s=frand(), l=frand();
		double rr, gr, br; hsl_to_rgb_oracle(h, s, l, &rr, &gr, &br);
		float rn, gn, bn;  hsl_to_rgbf(h, s, l, &rn, &gn, &bn);
		if (fabs(rr-rn) > maxr) maxr = fabs(rr-rn);
		if (fabs(gr-gn) > maxg) maxg = fabs(gr-gn);
		if (fabs(br-bn) > maxb) maxb = fabs(br-bn);
	}
	cr_expect(maxr < 1e-4 && maxg < 1e-4 && maxb < 1e-4,
	          "inverse err R %.3e G %.3e B %.3e", maxr, maxg, maxb);
}

Test(color_hsl, tie_for_max_hue_is_exact) {
	/* Two-channel ties at the clip boundary: the old independent-mask version
	 * was 60 deg off here; the exclusive-mask version must match the oracle. */
	struct { float r, g, b; } ties[] = {
		{1,1,0}, {0,1,1}, {1,0,1}, {1,1,0.5f}, {0.5f,1,1}, {0.8f,0.8f,0.2f},
	};
	for (unsigned c = 0; c < sizeof(ties)/sizeof(ties[0]); c++) {
		double hr, sr, lr; rgb_to_hsl_oracle(ties[c].r, ties[c].g, ties[c].b, &hr, &sr, &lr);
		float hn, sn, ln;  rgb_to_hslf(ties[c].r, ties[c].g, ties[c].b, &hn, &sn, &ln);
		cr_expect(hue_err(hr, hn) < 1e-5,
		          "tie (%.2f,%.2f,%.2f): oracle h=%.5f core h=%.5f",
		          ties[c].r, ties[c].g, ties[c].b, hr, hn);
	}
}

Test(color_hsl, round_trip_is_stable) {
	rng_reset();
	double rtmax = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r=frand(), g=frand(), b=frand();
		float h, s, l;    rgb_to_hslf(r, g, b, &h, &s, &l);
		float r2, g2, b2; hsl_to_rgbf(h, s, l, &r2, &g2, &b2);
		double e = fmax(fmax(fabs(r-r2), fabs(g-g2)), fabs(b-b2));
		if (e > rtmax) rtmax = e;
	}
	cr_expect(rtmax < 1e-5, "round-trip max err %.3e", rtmax);
}

/* ===================== double HSL shipped vs oracle ==================== */
/* The shipped rgb_to_hsl/hsl_to_rgb are now branchless too; verify they still
 * match the original branchy double implementation (the oracle). */

Test(color_hsl_double, shipped_matches_oracle) {
	rng_reset();
	double maxh = 0, maxs = 0, maxl = 0, rt = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r, g, b;
		if (i & 1) { r=fquant(); g=fquant(); b=fquant(); }
		else       { r=frand();  g=frand();  b=frand();  }
		double ho, so, lo; rgb_to_hsl_oracle(r, g, b, &ho, &so, &lo);
		double hn, sn, ln; rgb_to_hsl((double)r, (double)g, (double)b, &hn, &sn, &ln);
		if (hue_err(ho, hn) > maxh) maxh = hue_err(ho, hn);
		if (fabs(so-sn) > maxs) maxs = fabs(so-sn);
		if (fabs(lo-ln) > maxl) maxl = fabs(lo-ln);
		double rr, gg, bb; hsl_to_rgb(hn, sn, ln, &rr, &gg, &bb);
		double e = fmax(fmax(fabs(r-rr), fabs(g-gg)), fabs(b-bb));
		/* round trip only meaningful for chromatic, non-clipped pixels */
		if (ln > 0 && sn > 0 && e > rt) rt = e;
	}
	cr_expect(maxh < 1e-12 && maxs < 1e-12 && maxl < 1e-12,
	          "double HSL vs oracle: h %.3e s %.3e l %.3e", maxh, maxs, maxl);
	cr_expect(rt < 1e-9, "double HSL round-trip err %.3e", rt);
}

/* ===================== float_sat HSL (h on [0,6]) ====================== */

static void rgb_to_hsl_float_sat_oracle(float r, float g, float b, float low,
                                        float *h, float *s, float *l) {
	float v = (r > g) ? r : g; v = (v > b) ? v : b;
	float m = (r < g) ? r : g; m = (m < b) ? m : b;
	float vm;
	if (m + v < low + low) { *l = 0.f; return; }
	*l = (m + v) / 2.f; *h = 0.f; *s = 0.f;
	if ((*s = vm = v - m) > 0.f) *s /= (*l <= 0.5f) ? (v + m) : (2.f - v - m);
	else return;
	if (r == v)      { float g2=(v-g)/vm, b2=(v-b)/vm; *h = (g==m ? 5.f+b2 : 1.f-g2); }
	else if (g == v) { float r2=(v-r)/vm, b2=(v-b)/vm; *h = (b==m ? 1.f+r2 : 3.f-b2); }
	else             { float r2=(v-r)/vm, g2=(v-g)/vm; *h = (r==m ? 3.f+g2 : 5.f-r2); }
}

static double hue_err6(double a, double b) {  /* circular distance on [0,6) */
	double d = fabs(a - b);
	return (d > 3.0) ? 6.0 - d : d;
}

Test(color_hsl_float_sat, forward_matches_oracle_no_floor) {
	rng_reset();
	double maxh = 0, maxs = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r, g, b;
		if (i & 1) { r=fquant(); g=fquant(); b=fquant(); }
		else       { r=frand();  g=frand();  b=frand();  }
		float ho=0, so=0, lo=0; rgb_to_hsl_float_sat_oracle(r, g, b, 0.f, &ho, &so, &lo);
		float hn=0, sn=0, ln=0; rgb_to_hsl_float_sat(r, g, b, 0.f, &hn, &sn, &ln);
		float vm = ((r>g?r:g)>b?(r>g?r:g):b) - ((r<g?r:g)<b?(r<g?r:g):b);
		if (vm > 0.f) {                   /* h,s only defined for chromatic pixels */
			if (hue_err6(ho, hn) > maxh) maxh = hue_err6(ho, hn);
			if (fabs(so-sn) > maxs) maxs = fabs(so-sn);
		}
	}
	cr_expect(maxh < 1e-4, "float_sat hue (0..6) err %.3e", maxh);
	cr_expect(maxs < 1e-4, "float_sat sat err %.3e", maxs);
}

Test(color_hsl_float_sat, floor_zeroes_luminance) {
	/* Below the 2*low floor, l must be 0 so the saturation caller skips. */
	float h, s, l;
	rgb_to_hsl_float_sat(0.10f, 0.12f, 0.08f, 0.5f, &h, &s, &l);
	cr_expect(l == 0.f, "expected l==0 below floor, got %.5f", l);
	rgb_to_hsl_float_sat(0.80f, 0.82f, 0.78f, 0.5f, &h, &s, &l);
	cr_expect(l > 0.f, "expected l>0 above floor, got %.5f", l);
}

Test(color_hsl_float_sat, round_trip_is_stable) {
	rng_reset();
	double rtmax = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r=frand(), g=frand(), b=frand();
		float h, s, l;    rgb_to_hsl_float_sat(r, g, b, 0.f, &h, &s, &l);
		float r2, g2, b2; hsl_to_rgb_float_sat(h, s, l, &r2, &g2, &b2);
		double e = fmax(fmax(fabs(r-r2), fabs(g-g2)), fabs(b-b2));
		if (e > rtmax) rtmax = e;
	}
	cr_expect(rtmax < 1e-5, "float_sat round-trip max err %.3e", rtmax);
}

/* ============================ HSV <-> RGB ============================== */

static void rgb_to_hsv_oracle(double r, double g, double b,
                              double *h, double *s, double *v) {
	double cmax = (r > g) ? r : g; cmax = (cmax > b) ? cmax : b;
	double cmin = (r < g) ? r : g; cmin = (cmin < b) ? cmin : b;
	double delta = cmax - cmin; *v = cmax;
	if (delta == 0.0) { *s = 0.0; *h = 0.0; return; }
	*s = delta / cmax;
	if (cmax == r)      *h = (((g - b) / delta)) / 6.0;
	else if (cmax == g) *h = (((b - r) / delta) + 2.0) / 6.0;
	else                *h = (((r - g) / delta) + 4.0) / 6.0;
	if (*h < 0.0) *h += 1.0;
}

static void hsv_to_rgb_oracle(double h, double s, double v,
                              double *r, double *g, double *b) {
	if (h >= 1.0) h -= 1.0;
	h *= 6.0; int i = (int)h; double f = h - (double)i;
	double p = v*(1.0-s), q = v*(1.0-(s*f)), t = v*(1.0-(s*(1.0-f)));
	switch (i) {
		case 0: *r=v; *g=t; *b=p; break;
		case 1: *r=q; *g=v; *b=p; break;
		case 2: *r=p; *g=v; *b=t; break;
		case 3: *r=p; *g=q; *b=v; break;
		case 4: *r=t; *g=p; *b=v; break;
		default:*r=v; *g=p; *b=q; break;
	}
}

Test(color_hsv, forward_matches_oracle_and_ties) {
	rng_reset();
	double maxh = 0, maxs = 0, maxv = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r, g, b;
		if (i & 1) { r=fquant(); g=fquant(); b=fquant(); }
		else       { r=frand();  g=frand();  b=frand();  }
		double ho, so, vo; rgb_to_hsv_oracle(r, g, b, &ho, &so, &vo);
		float hn, sn, vn;  rgb_to_hsvf(r, g, b, &hn, &sn, &vn);
		if (hue_err(ho, hn) > maxh) maxh = hue_err(ho, hn);
		if (fabs(so-sn) > maxs) maxs = fabs(so-sn);
		if (fabs(vo-vn) > maxv) maxv = fabs(vo-vn);
	}
	cr_expect(maxh < 1e-4, "HSVf hue err %.3e", maxh);
	cr_expect(maxs < 1e-4, "HSVf sat err %.3e", maxs);
	cr_expect(maxv < 1e-6, "HSVf val err %.3e", maxv);
}

Test(color_hsv, inverse_matches_oracle) {
	rng_reset();
	double maxr = 0, maxg = 0, maxb = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		double h=frand(), s=frand(), v=frand();
		double rr, gr, br; hsv_to_rgb_oracle(h, s, v, &rr, &gr, &br);
		double rn, gn, bn; hsv_to_rgb(h, s, v, &rn, &gn, &bn);
		if (fabs(rr-rn) > maxr) maxr = fabs(rr-rn);
		if (fabs(gr-gn) > maxg) maxg = fabs(gr-gn);
		if (fabs(br-bn) > maxb) maxb = fabs(br-bn);
	}
	cr_expect(maxr < 1e-12 && maxg < 1e-12 && maxb < 1e-12,
	          "HSV->RGB vs oracle: R %.3e G %.3e B %.3e", maxr, maxg, maxb);
}

Test(color_hsv, float_round_trip_is_stable) {
	rng_reset();
	double rtmax = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r=frand(), g=frand(), b=frand();
		float h, s, v;    rgb_to_hsvf(r, g, b, &h, &s, &v);
		float r2, g2, b2; hsv_to_rgbf(h, s, v, &r2, &g2, &b2);
		double e = fmax(fmax(fabs(r-r2), fabs(g-g2)), fabs(b-b2));
		if (e > rtmax) rtmax = e;
	}
	cr_expect(rtmax < 1e-5, "HSV float round-trip err %.3e", rtmax);
}

/* ============================ YUV <-> RGB =============================== */

Test(color_yuv, round_trip_is_near_exact) {
	rng_reset();
	double rtmax = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r=frand(), g=frand(), b=frand();
		float y, u, v;    rgb_to_yuvf(r, g, b, &y, &u, &v);
		float r2, g2, b2; yuv_to_rgbf(y, u, v, &r2, &g2, &b2);
		double e = fmax(fmax(fabs(r-r2), fabs(g-g2)), fabs(b-b2));
		if (e > rtmax) rtmax = e;
	}
	cr_expect(rtmax < 1e-5, "YUV round-trip max err %.3e", rtmax);
}

/* ============================ RGB <-> XYZ <-> LAB ====================== */

Test(color_xyz, rgb_round_trip_float) {
	rng_reset();
	double rtmax = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		float r=frand(), g=frand(), b=frand();
		float x, y, z, r2, g2, b2;
		rgb_to_xyzf(r, g, b, &x, &y, &z);
		xyz_to_rgbf(x, y, z, &r2, &g2, &b2);
		double e = fmax(fmax(fabs(r-r2), fabs(g-g2)), fabs(b-b2));
		if (e > rtmax) rtmax = e;
	}
	cr_expect(rtmax < 1e-3, "RGB->XYZ->RGB (float) round-trip err %.3e", rtmax);
}

Test(color_lab, xyz_round_trip_float) {
	rng_reset();
	double rtmaxf = 0;
	for (int i = 0; i < N_FUZZ; i++) {
		/* feed physically meaningful XYZ by routing through rgb_to_xyzf */
		float r=frand(), g=frand(), b=frand();
		float x, y, z; rgb_to_xyzf(r, g, b, &x, &y, &z);
		float L, A, B, x2, y2, z2;
		xyz_to_LABf(x, y, z, &L, &A, &B);
		LAB_to_xyzf(L, A, B, &x2, &y2, &z2);
		double ef = fmax(fmax(fabs(x-x2), fabs(y-y2)), fabs(z-z2));
		if (ef > rtmaxf) rtmaxf = ef;
	}
	cr_expect(rtmaxf < 1e-2, "XYZ->LAB->XYZ (float) round-trip err %.3e", rtmaxf);
}

/* ============================ linRGB <-> XYZ =========================== */

Test(color_linrgb, xyz_round_trip_both_scales) {
	rng_reset();
	for (int sc = 0; sc <= 1; sc++) {
		gboolean scale = sc ? TRUE : FALSE;
		double rtmaxf = 0;
		for (int i = 0; i < N_FUZZ / 2; i++) {
			float r=frand(), g=frand(), b=frand();
			float xf, yf, zf, rf, gf, bf;
			linrgb_to_xyzf(r, g, b, &xf, &yf, &zf, scale);
			xyz_to_linrgbf(xf, yf, zf, &rf, &gf, &bf, scale);
			double ef = fmax(fmax(fabs(r-rf), fabs(g-gf)), fabs(b-bf));
			if (ef > rtmaxf) rtmaxf = ef;
		}
		cr_expect(rtmaxf < 1e-5, "linRGB<->XYZ (float, scale=%d) err %.3e", sc, rtmaxf);
	}
}
