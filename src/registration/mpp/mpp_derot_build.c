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

gboolean mpp_derot_autodetect_disk(const fits *fit, double *cx, double *cy,
                                   double *radius) {
	if (!fit || !cx || !cy || !radius || fit->rx <= 0 || fit->ry <= 0)
		return FALSE;
	const int rx = (int) fit->rx, ry = (int) fit->ry;

	/* Threshold at a fraction of the per-image maximum (use the first channel;
	 * a planet disk is bright against a dark sky). Works for 8/16-bit USHORT
	 * and 32-bit float buffers. */
	double maxv = 0.0;
	const gboolean is_float = (fit->type == DATA_FLOAT);
	const size_t npix = (size_t) rx * ry;
	if (is_float) {
		const float *p = fit->fdata;
		if (!p) return FALSE;
		for (size_t i = 0; i < npix; i++) if (p[i] > maxv) maxv = p[i];
	} else {
		const WORD *p = fit->data;
		if (!p) return FALSE;
		for (size_t i = 0; i < npix; i++) if (p[i] > maxv) maxv = p[i];
	}
	if (maxv <= 0.0) return FALSE;
	const double thr = 0.25 * maxv;

	long sx = 0, sy = 0, count = 0;
	int minx = rx, maxx = -1, miny = ry, maxy = -1;
	for (int y = 0; y < ry; y++) {
		for (int x = 0; x < rx; x++) {
			double v = is_float ? (double) fit->fdata[(size_t) y * rx + x]
			                    : (double) fit->data[(size_t) y * rx + x];
			if (v >= thr) {
				sx += x; sy += y; count++;
				if (x < minx) minx = x; if (x > maxx) maxx = x;
				if (y < miny) miny = y; if (y > maxy) maxy = y;
			}
		}
	}
	if (count < 16 || maxx < minx || maxy < miny) return FALSE;
	*cx = (double) sx / (double) count;
	*cy = (double) sy / (double) count;
	*radius = 0.25 * ((maxx - minx) + (maxy - miny));   /* mean semi-extent */
	return (*radius > 1.0);
}
