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

#include <math.h>

#include "core/siril.h"
#include "core/siril_date.h"
#include "io/ser.h"

#include "registration/mpp/mpp_derot_build.h"

gboolean mpp_derot_frame_times(sequence *seq, double fps, double start_jd,
                               double *jd_out) {
	if (!seq || !jd_out || seq->number <= 0) return FALSE;
	const int N = seq->number;

	if (fps > 0.0) {
		double jd0 = start_jd;
		if (isnan(jd0)) {
			if (seq->type == SEQ_SER && seq->ser_file) {
				GDateTime *t0 = ser_read_frame_date(seq->ser_file, 0);
				if (!t0) return FALSE;
				jd0 = date_time_to_Julian(t0);
				g_date_time_unref(t0);
			} else {
				return FALSE;
			}
		}
		for (int i = 0; i < N; i++)
			jd_out[i] = jd0 + (double) i / (fps * 86400.0);
		return TRUE;
	}

	if (seq->type == SEQ_SER && seq->ser_file) {
		for (int i = 0; i < N; i++) {
			GDateTime *t = ser_read_frame_date(seq->ser_file, i);
			if (!t) return FALSE;
			jd_out[i] = date_time_to_Julian(t);
			g_date_time_unref(t);
		}
		return TRUE;
	}
	return FALSE;
}

double mpp_derot_midpoint_epoch(const double *jd, int n) {
	if (!jd || n <= 0) return 0.0;
	return 0.5 * (jd[0] + jd[n - 1]);
}

gboolean mpp_derot_sequence_span(sequence *seq, double fps, double start_jd,
                                 double *first_jd, double *last_jd) {
	if (!seq || seq->number <= 0 || !first_jd || !last_jd) return FALSE;
	const int N = seq->number;

	if (fps > 0.0) {
		double jd0 = start_jd;
		if (isnan(jd0)) {
			if (seq->type == SEQ_SER && seq->ser_file) {
				GDateTime *t0 = ser_read_frame_date(seq->ser_file, 0);
				if (!t0) return FALSE;
				jd0 = date_time_to_Julian(t0);
				g_date_time_unref(t0);
			} else {
				return FALSE;
			}
		}
		*first_jd = jd0;
		*last_jd  = jd0 + (double) (N - 1) / (fps * 86400.0);
		return TRUE;
	}

	if (seq->type == SEQ_SER && seq->ser_file) {
		GDateTime *t0 = ser_read_frame_date(seq->ser_file, 0);
		GDateTime *tn = ser_read_frame_date(seq->ser_file, N - 1);
		if (!t0 || !tn) {
			if (t0) g_date_time_unref(t0);
			if (tn) g_date_time_unref(tn);
			return FALSE;
		}
		*first_jd = date_time_to_Julian(t0);
		*last_jd  = date_time_to_Julian(tn);
		g_date_time_unref(t0);
		g_date_time_unref(tn);
		return TRUE;
	}
	return FALSE;
}

double mpp_derot_union_epoch(const double *first_jd, const double *last_jd, int n) {
	if (!first_jd || !last_jd || n <= 0) return 0.0;
	double gmin = first_jd[0], gmax = last_jd[0];
	for (int i = 1; i < n; i++) {
		if (first_jd[i] < gmin) gmin = first_jd[i];
		if (last_jd[i]  > gmax) gmax = last_jd[i];
	}
	return 0.5 * (gmin + gmax);
}

mpp_derot_t *mpp_derot_build(planet_body_t body, int system, double epoch_jd,
                             const double *jd, int num_frames,
                             int frame_rows, int frame_cols,
                             double cx, double cy, double radius,
                             double pa_deg, double parity,
                             double obs_lat, double obs_lon, double obs_elev) {
	if (!jd || num_frames <= 0 || system < 1 || system > 3) return NULL;

	planet_geom_t ge;
	if (planet_ephemeris(body, epoch_jd, obs_lat, obs_lon, obs_elev, &ge) != 0)
		return NULL;

	mpp_derot_t *d = mpp_derot_alloc(num_frames);
	if (!d) return NULL;
	d->body = body;
	d->rot_system = system - 1;
	d->ephem_version = 1;
	d->frame_rows = frame_rows;
	d->frame_cols = frame_cols;
	d->epoch_jd = epoch_jd;
	d->obs_lat = obs_lat; d->obs_lon = obs_lon; d->obs_elev = obs_elev;
	d->cx = cx; d->cy = cy; d->r_eq = radius; d->parity = parity;
	d->pole_angle_epoch = pa_deg * M_PI / 180.0;
	d->flattening = ge.flattening;
	d->epoch_sub_obs_lat = ge.sub_obs_lat;
	d->epoch_cm = ge.cm[system - 1];
	d->epoch_pole_pa = ge.pole_pa;

	for (int i = 0; i < num_frames; i++) {
		planet_geom_t g;
		planet_ephemeris(body, jd[i], obs_lat, obs_lon, obs_elev, &g);
		d->jd[i] = jd[i];
		d->sub_obs_lat[i] = g.sub_obs_lat;
		d->cm[i] = g.cm[system - 1];
		d->pole_pa[i] = g.pole_pa;
	}
	return d;
}

gboolean mpp_derot_autodetect_disk(const fits *fit, double flattening,
                                   gboolean has_rings,
                                   double *cx, double *cy, double *radius,
                                   double *major_axis_deg) {
	if (!fit || !cx || !cy || !radius || fit->rx <= 0 || fit->ry <= 0)
		return FALSE;
	const int rx = (int) fit->rx, ry = (int) fit->ry;
	const gboolean is_float = (fit->type == DATA_FLOAT);
	const float *fp = is_float ? fit->fdata : NULL;
	const WORD  *wp = is_float ? NULL : fit->data;
	if ((is_float && !fp) || (!is_float && !wp)) return FALSE;
	const size_t npix = (size_t) rx * ry;
#define PIXV(i) (is_float ? (double) fp[(i)] : (double) wp[(i)])

	/* Peak, and a sky-background estimate from the frame border (the planet is
	 * roughly centred, so the edge pixels are sky). Works for 8/16-bit USHORT
	 * and 32-bit float buffers. */
	double maxv = 0.0;
	for (size_t i = 0; i < npix; i++) { const double v = PIXV(i); if (v > maxv) maxv = v; }
	if (maxv <= 0.0) return FALSE;
	double bsum = 0.0; long bn = 0;
	for (int x = 0; x < rx; x++) {
		bsum += PIXV((size_t) x) + PIXV((size_t)(ry - 1) * rx + x); bn += 2;
	}
	for (int y = 0; y < ry; y++) {
		bsum += PIXV((size_t) y * rx) + PIXV((size_t) y * rx + (rx - 1)); bn += 2;
	}
	const double bg = bn > 0 ? bsum / (double) bn : 0.0;
	if (maxv <= bg) return FALSE;
	/* Threshold as a fraction of the peak above sky. The fraction is low on
	 * purpose: the contour must reach the actual limb, which on a real (limb-
	 * darkened, seeing-blurred) disk sits well below a fraction of the bright
	 * disk peak — and the minor axis used below runs through the dimmest part
	 * of the globe (the poles), so a high fraction clips it and underestimates
	 * the radius on both Jupiter and Saturn. 0.10 was tuned against Jupiter and
	 * Saturn frames to land on the visible limb (the disc the user sees on a log
	 * stretch); it stays above the ~0.07 regime where faint glow starts to
	 * inflate the moments. */
	const double thr = bg + 0.10 * (maxv - bg);

	/* Centroid and second moments of the thresholded region (geometric, binary
	 * mask: the second moment along a semi-axis r of a uniform region is
	 * r^2/4). */
	double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
	long count = 0;
	for (int y = 0; y < ry; y++) {
		for (int x = 0; x < rx; x++) {
			if (PIXV((size_t) y * rx + x) >= thr) {
				sx += x; sy += y;
				sxx += (double) x * x; syy += (double) y * y; sxy += (double) x * y;
				count++;
			}
		}
	}
	if (count < 16) return FALSE;
	const double n = (double) count;
	const double mxc = sx / n, myc = sy / n;
	const double cxx = sxx / n - mxc * mxc;
	const double cyy = syy / n - myc * myc;
	const double cxy = sxy / n - mxc * myc;
	*cx = mxc;
	*cy = myc;
	/* Major-axis orientation (Saturn's ring plane / Jupiter's equator), image
	 * coordinates with y up; the caller turns it into the pole angle. */
	const double theta = 0.5 * atan2(2.0 * cxy, cxx - cyy);
	if (major_axis_deg) *major_axis_deg = theta * 180.0 / M_PI;

	const double f = (flattening > 0.0 && flattening < 0.9) ? flattening : 0.0;

	if (has_rings) {
		/* Saturn: the minor axis is unreliable here — the rings pile mass near
		 * the ring plane and shrink the perpendicular variance, under-sizing the
		 * globe. Instead measure the bright rings (their tips define the major
		 * axis) and scale by the precisely known A-ring-outer / equatorial-radius
		 * ratio. The tip is the 99th percentile of the on-major-axis distance
		 * over a generous-threshold mask — robust to a stray moon (which the max
		 * is not). */
		const double thr_ring = bg + 0.03 * (maxv - bg);
		const double ct = cos(theta), st = sin(theta);
		const int maxd = rx + ry;
		int *hist = (int *) calloc((size_t) maxd + 1, sizeof(int));
		long rc = 0;
		if (hist) {
			for (int y = 0; y < ry; y++) {
				for (int x = 0; x < rx; x++) {
					if (PIXV((size_t) y * rx + x) >= thr_ring) {
						double pr = fabs((x - mxc) * ct + (y - myc) * st);
						int bin = (int) (pr + 0.5);
						if (bin < 0) bin = 0; else if (bin > maxd) bin = maxd;
						hist[bin]++; rc++;
					}
				}
			}
			if (rc > 0) {
				const long target = (long) (0.99 * (double) rc);
				long cum = 0; int tip = 0;
				for (int bbin = 0; bbin <= maxd; bbin++) {
					cum += hist[bbin];
					if (cum >= target) { tip = bbin; break; }
				}
				/* A-ring outer edge / equatorial radius = 136775 km / 60268 km. */
				if (tip > 1) *radius = (double) tip / 2.269;
			}
			free(hist);
		}
		return (*radius > 1.0);
	}
#undef PIXV

	/* No rings: smaller eigenvalue of the covariance = the minor (≈ polar) axis;
	 * recover the equatorial radius from the body flattening. */
	const double tr = cxx + cyy;
	const double disc = sqrt(fmax(0.0, 0.25 * (cxx - cyy) * (cxx - cyy) + cxy * cxy));
	const double lam_minor = 0.5 * tr - disc;
	if (lam_minor <= 0.0) return FALSE;
	const double b = 2.0 * sqrt(lam_minor);
	*radius = b / (1.0 - f);
	return (*radius > 1.0);
}
