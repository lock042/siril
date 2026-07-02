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

#include <glib/gstdio.h>

#include "core/siril.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
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
		int disorder = 0;
		for (int i = 0; i < N; i++) {
			GDateTime *t = ser_read_frame_date(seq->ser_file, i);
			if (!t) return FALSE;
			jd_out[i] = date_time_to_Julian(t);
			g_date_time_unref(t);
			if (i > 0 && jd_out[i] < jd_out[i - 1])
				disorder++;
		}
		/* Derotation itself is fine with unordered stamps (each frame's own
		 * JD drives its map), but they usually mean the capture software
		 * wrote a broken trailer, and anything downstream that assumes
		 * chronological order (span estimates, the midpoint epoch) is
		 * quietly degraded. Tell the user how to repair it. */
		if (disorder > 0)
			siril_log_message(_("derotate: %d of %d frame timestamps are out "
			                    "of order — derotation uses each frame's own "
			                    "timestamp, but the capture trailer looks "
			                    "unreliable; `ser_fix_timestamps` can rewrite "
			                    "a uniform cadence\n"), disorder, N);
		return TRUE;
	}
	return FALSE;
}

double mpp_derot_midpoint_epoch(const double *jd, int n) {
	if (!jd || n <= 0) return 0.0;
	/* min/max rather than first/last: real SER trailers are sometimes out
	 * of chronological order, and the epoch should sit mid-span anyway. */
	double lo = jd[0], hi = jd[0];
	for (int i = 1; i < n; i++) {
		if (jd[i] < lo) lo = jd[i];
		if (jd[i] > hi) hi = jd[i];
	}
	return 0.5 * (lo + hi);
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
                             double obs_lat, double obs_lon, double obs_elev,
                             mpp_field_rot_t field_rot) {
	if (!jd || num_frames <= 0 || system < 1 || system > 3) return NULL;
	if (field_rot == MPP_FIELD_ROT_ALTAZ && (isnan(obs_lat) || isnan(obs_lon))) {
		siril_log_error(_("derotate: alt-az field rotation needs the observer "
		                  "site — supply -obs-lat and -obs-lon\n"));
		return NULL;
	}

	if (field_rot == MPP_FIELD_ROT_NONE
	    && (!isnan(obs_lat) || !isnan(obs_lon) || !isnan(obs_elev)))
		/* Without field rotation the observer coordinates only ride along
		 * for provenance; the ephemeris itself is geocentric (the
		 * sub-observer geometry differs by well under 0.001 deg). */
		siril_log_message(_("derotate: observer coordinates are recorded for "
		                    "provenance only — the ephemeris is geocentric "
		                    "(difference < 0.001° for these bodies)\n"));

	planet_geom_t ge;
	if (planet_ephemeris(body, epoch_jd, obs_lat, obs_lon, obs_elev, &ge) != 0)
		return NULL;
	const double q0 = (field_rot == MPP_FIELD_ROT_ALTAZ)
	    ? planet_parallactic_angle(epoch_jd, ge.ra, ge.dec, obs_lat, obs_lon)
	    : 0.0;

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
	d->pole_angle_epoch = pa_deg * G_PI / 180.0;
	d->flattening = ge.flattening;
	d->epoch_sub_obs_lat = ge.sub_obs_lat;
	d->epoch_cm = ge.cm[system - 1];
	d->epoch_pole_pa = ge.pole_pa;

	for (int i = 0; i < num_frames; i++) {
		planet_geom_t g;
		if (planet_ephemeris(body, jd[i], obs_lat, obs_lon, obs_elev, &g) != 0) {
			/* A frame with a nonsense timestamp (unset / corrupt) must not
			 * bake garbage geometry into the plan. */
			siril_log_error(_("derotate: frame %d has an unusable timestamp "
			                  "(JD %.6f) — fix the timestamps or supply "
			                  "-fps/-start\n"), i, jd[i]);
			mpp_derot_free(d);
			return NULL;
		}
		d->jd[i] = jd[i];
		d->sub_obs_lat[i] = g.sub_obs_lat;
		d->cm[i] = g.cm[system - 1];
		d->pole_pa[i] = g.pole_pa;
		if (field_rot == MPP_FIELD_ROT_ALTAZ) {
			/* Alt-az field rotation, folded into the per-frame pole angle.
			 * The map builder turns the ON-SKY drift (pole_pa[i] −
			 * epoch_pole_pa) into the in-image angle (with parity), so any
			 * additional on-sky-frame rotation of the IMAGE belongs in the
			 * same slot. With the camera fixed to the horizon frame, the
			 * in-image direction of the zenith is constant while sky north
			 * sits at −q from it (q = PA of the zenith, N→E positive):
			 * the in-image orientation of everything on the sky drifts by
			 * −Δq. Hence subtract (q_i − q_epoch). Per-frame RA/Dec keeps
			 * long -epoch-from spans honest. */
			const double qi = planet_parallactic_angle(jd[i], g.ra, g.dec,
			                                           obs_lat, obs_lon);
			d->pole_pa[i] -= qi - q0;
		}
	}
	if (field_rot == MPP_FIELD_ROT_ALTAZ) {
		/* Mark the fold (and the site it used, already stored in obs_lat /
		 * obs_lon) so the registration-side option recognises it, keeps it
		 * idempotent, and can unfold or re-base it. */
		d->flags |= MPP_DEROT_FLAG_FIELDROT;
		siril_log_message(_("derotate: alt-az field rotation folded in — "
		                    "%.2f° across the span (parallactic angle %.2f° "
		                    "at the epoch)\n"),
		                  fabs(planet_parallactic_angle(jd[num_frames - 1],
		                           ge.ra, ge.dec, obs_lat, obs_lon)
		                       - planet_parallactic_angle(jd[0], ge.ra, ge.dec,
		                           obs_lat, obs_lon)),
		                  q0);
	}
	return d;
}

/* The ephemeris target whose sky position drives a plan's parallactic
 * angle: the planet itself for a planetary plan (planet body ids match
 * ephem_target_t), the stored target for a rotation-only plan. */
static int derot_q_target(const mpp_derot_t *d) {
	return d->body == MPP_DEROT_BODY_FIELDROT ? d->rot_system : d->body;
}

/* Fold (direction +1) or unfold (direction -1) the alt-az parallactic
 * drift -(q_i - q_epoch) into the plan's per-frame pole angles, with q
 * evaluated for the plan's own target and timestamps at the given site.
 * Because q is recomputable exactly from what the plan stores, a fold done
 * with the site recorded in obs_lat/obs_lon is exactly reversible. */
static mpp_status_t derot_apply_field_rotation(mpp_derot_t *d, int direction,
                                               double lat, double lon) {
	if (!d || isnan(lat) || isnan(lon)) return MPP_EINVAL;
	double ra, dec;
	const int target = derot_q_target(d);
	if (ephem_target_radec((ephem_target_t) target, d->epoch_jd, &ra, &dec) != 0)
		return MPP_EINVAL;
	const double q0 = planet_parallactic_angle(d->epoch_jd, ra, dec, lat, lon);
	for (int i = 0; i < d->num_frames; i++) {
		if (ephem_target_radec((ephem_target_t) target, d->jd[i], &ra, &dec) != 0)
			return MPP_EINVAL;
		const double qi = planet_parallactic_angle(d->jd[i], ra, dec, lat, lon);
		d->pole_pa[i] -= direction * (qi - q0);
	}
	return MPP_OK;
}

mpp_derot_t *mpp_derot_build_field_rotation(int target, const double *jd,
                                            int num_frames,
                                            int frame_rows, int frame_cols,
                                            double obs_lat, double obs_lon) {
	if (!jd || num_frames <= 0) return NULL;
	if (isnan(obs_lat) || isnan(obs_lon)) {
		siril_log_error(_("field rotation: the observer site is required — "
		                  "supply latitude and east longitude\n"));
		return NULL;
	}

	const double epoch = mpp_derot_midpoint_epoch(jd, num_frames);
	double ra0, dec0;
	if (ephem_target_radec((ephem_target_t) target, epoch, &ra0, &dec0) != 0) {
		siril_log_error(_("field rotation: cannot compute the target position "
		                  "(bad timestamps?)\n"));
		return NULL;
	}
	const double q0 = planet_parallactic_angle(epoch, ra0, dec0,
	                                           obs_lat, obs_lon);

	mpp_derot_t *d = mpp_derot_alloc(num_frames);
	if (!d) return NULL;
	d->body = MPP_DEROT_BODY_FIELDROT;
	d->rot_system = target;          /* provenance: which target drove q */
	d->ephem_version = 1;
	d->frame_rows = frame_rows;
	d->frame_cols = frame_cols;
	d->epoch_jd = epoch;
	d->obs_lat = obs_lat; d->obs_lon = obs_lon; d->obs_elev = NAN;
	/* Rotation-only geometry: with B/CM equal on both sides of the map and
	 * flattening 0, every pixel — on or off the nominal disk — maps through
	 * the same rigid similarity, so the disk parameters only pick the
	 * rotation centre. */
	d->cx = frame_cols / 2.0;
	d->cy = frame_rows / 2.0;
	d->r_eq = (frame_rows < frame_cols ? frame_rows : frame_cols) / 4.0;
	d->flattening = 0.0;
	d->parity = 1.0;
	d->pole_angle_epoch = 0.0;
	d->epoch_sub_obs_lat = 0.0;
	d->epoch_cm = 0.0;
	d->epoch_pole_pa = 0.0;

	for (int i = 0; i < num_frames; i++) {
		d->jd[i] = jd[i];
		d->sub_obs_lat[i] = 0.0;
		d->cm[i] = 0.0;
		d->pole_pa[i] = 0.0;
	}
	/* Same sign convention as the planetary fold: in-image sky orientation
	 * drifts by −Δq on an alt-az mount. */
	if (derot_apply_field_rotation(d, +1, obs_lat, obs_lon) != MPP_OK) {
		siril_log_error(_("field rotation: unusable timestamps in the "
		                  "sequence\n"));
		mpp_derot_free(d);
		return NULL;
	}
	d->flags |= MPP_DEROT_FLAG_FIELDROT;
	siril_log_message(_("field rotation: alt-az plan built — %.2f° across the "
	                    "span (parallactic angle %.2f° at the epoch)\n"),
	                  fabs(d->pole_pa[num_frames - 1] - d->pole_pa[0]), q0);
	return d;
}

mpp_status_t mpp_derot_sync_field_rotation_plan(sequence *seq, gboolean enable,
                                                int target,
                                                double obs_lat, double obs_lon) {
	if (!seq || !seq->seqname) return MPP_EINVAL;
	gchar *path = g_strdup_printf("%s.derot", seq->seqname);

	mpp_derot_t *existing = NULL;
	const gboolean have = (mpp_derot_read(path, &existing) == MPP_OK && existing);
	const gboolean is_synth = have
	    && existing->body == MPP_DEROT_BODY_FIELDROT;
	const gboolean is_folded = have
	    && (existing->flags & MPP_DEROT_FLAG_FIELDROT);
	mpp_status_t rc = MPP_OK;

	if (!enable) {
		if (is_synth) {
			/* Remove what this option created. */
			g_unlink(path);
			siril_log_message(_("field rotation: removed the synthesised "
			                    "%s\n"), path);
		} else if (is_folded) {
			/* Planetary plan with a folded correction: unfold with the
			 * site it was folded for (recorded in the plan) and rewrite. */
			rc = derot_apply_field_rotation(existing, -1,
			                                existing->obs_lat,
			                                existing->obs_lon);
			if (rc == MPP_OK) {
				existing->flags &= ~MPP_DEROT_FLAG_FIELDROT;
				rc = mpp_derot_write(path, existing);
				if (rc == MPP_OK)
					siril_log_message(_("field rotation: removed the alt-az "
					                    "correction from %s\n"), path);
			}
		}
		goto out;
	}

	if (isnan(obs_lat) || isnan(obs_lon)) {
		siril_log_error(_("field rotation: the observer site is required — "
		                  "supply latitude and east longitude\n"));
		rc = MPP_EINVAL;
		goto out;
	}

	if (have && !is_synth) {
		/* A real planetary derotation plan: fold the alt-az correction into
		 * it — derotation and field rotation compose (the near-zenith
		 * planet case). The plan's own body drives the parallactic angle;
		 * the requested target is ignored. Re-basing to a new site unfolds
		 * the recorded one first. */
		if (is_folded) {
			if (fabs(existing->obs_lat - obs_lat) < 1e-9
			    && fabs(existing->obs_lon - obs_lon) < 1e-9)
				goto out;   /* already folded for this site */
			rc = derot_apply_field_rotation(existing, -1,
			                                existing->obs_lat,
			                                existing->obs_lon);
			if (rc != MPP_OK) goto out;
		}
		rc = derot_apply_field_rotation(existing, +1, obs_lat, obs_lon);
		if (rc != MPP_OK) goto out;
		existing->obs_lat = obs_lat;
		existing->obs_lon = obs_lon;
		existing->flags |= MPP_DEROT_FLAG_FIELDROT;
		rc = mpp_derot_write(path, existing);
		if (rc == MPP_OK) {
			siril_log_message(_("field rotation: alt-az correction folded "
			                    "into the existing derotation plan %s\n"),
			                  path);
			if (target != existing->body)
				siril_log_message(_("field rotation: the plan's body drives "
				                    "the parallactic angle; the Target "
				                    "selection is ignored\n"));
		}
		goto out;
	}

	if (is_synth && existing->num_frames == seq->number) {
		/* Rebuild the rotation-only plan from its own stored timestamps
		 * (target or site may have changed) — no sequence re-read. */
		mpp_derot_t *d = mpp_derot_build_field_rotation(target, existing->jd,
		                                                existing->num_frames,
		                                                existing->frame_rows,
		                                                existing->frame_cols,
		                                                obs_lat, obs_lon);
		rc = d ? mpp_derot_write(path, d) : MPP_EINVAL;
		mpp_derot_free(d);
		goto out;
	}

	{
		double *jd = malloc((size_t) seq->number * sizeof(double));
		if (!jd) { rc = MPP_ENOMEM; goto out; }
		if (!mpp_derot_frame_times(seq, 0.0, NAN, jd)) {
			siril_log_error(_("field rotation: no usable per-frame timestamps "
			                  "— write them first (`ser_fix_timestamps ... "
			                  "-fps=F`)\n"));
			free(jd);
			rc = MPP_EINVAL;
			goto out;
		}
		mpp_derot_t *d = mpp_derot_build_field_rotation(target, jd, seq->number,
		                                                (int) seq->ry,
		                                                (int) seq->rx,
		                                                obs_lat, obs_lon);
		free(jd);
		rc = d ? mpp_derot_write(path, d) : MPP_EINVAL;
		mpp_derot_free(d);
	}

out:
	if (rc != MPP_OK && rc != MPP_EINVAL)
		siril_log_error(_("field rotation: failed to update %s\n"), path);
	if (existing) mpp_derot_free(existing);
	g_free(path);
	return rc;
}

gboolean mpp_derot_autodetect_disk(const fits *fit, double flattening,
                                   gboolean has_rings,
                                   double *cx, double *cy, double *radius,
                                   double *major_axis_deg) {
	if (!fit || !cx || !cy || !radius || fit->rx <= 0 || fit->ry <= 0)
		return FALSE;
	/* Every failure exit below must leave *radius deterministic — the ring
	 * path returns (*radius > 1.0) and previously read the caller's
	 * uninitialised variable when the tip estimate never landed. */
	*radius = 0.0;
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
	if (major_axis_deg) *major_axis_deg = theta * 180.0 / G_PI;

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
