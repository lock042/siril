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

#ifndef SRC_REGISTRATION_MPP_DEROT_BUILD_H_
#define SRC_REGISTRATION_MPP_DEROT_BUILD_H_

#include "core/siril.h"                          /* sequence */
#include "algos/planet_ephem.h"                  /* planet_body_t */
#include "registration/mpp/mpp_derot_sidecar.h"  /* mpp_derot_t */

#ifdef __cplusplus
extern "C" {
#endif

/* Shared builder for the .derot plan — used by the `derotate` command and the
 * derotation GUI window so both produce identical sidecars from the same
 * ephemeris path. */

/* Fill jd_out[seq->number] with each frame's UTC Julian date. When fps > 0 a
 * uniform cadence is used, anchored at start_jd (or, if start_jd is NaN, the
 * SER frame-0 timestamp). Otherwise per-frame SER timestamps are read. Returns
 * FALSE when no timing source is available (caller should ask for -fps). */
gboolean mpp_derot_frame_times(sequence *seq, double fps, double start_jd,
                               double *jd_out);

/* Midpoint of a per-frame Julian-date array (the default reference epoch). */
double mpp_derot_midpoint_epoch(const double *jd, int n);

/* First/last frame UTC Julian dates of a sequence, without reading every frame
 * (the endpoints are all the span needs). Same timing sources as
 * mpp_derot_frame_times. Returns FALSE if no timing source is available. */
gboolean mpp_derot_sequence_span(sequence *seq, double fps, double start_jd,
                                 double *first_jd, double *last_jd);

/* Common reference epoch for a set of sequences whose spans are given by
 * first_jd[]/last_jd[]: the midpoint of their union (0.5·(min first, max last)).
 * This is the epoch that minimises the largest rotation any one frame across
 * all the sequences is derotated through. Returns 0 for n <= 0. */
double mpp_derot_union_epoch(const double *first_jd, const double *last_jd, int n);

/* Build a fully-populated plan (per-frame + epoch geometry) from the built-in
 * ephemeris. `system` is 1..3. Disk fit is in original full-frame pixels;
 * pa_deg is the in-image pole position angle, parity ±1. obs_* may be NaN for
 * geocentric. Returns a heap plan (free with mpp_derot_free) or NULL. */
mpp_derot_t *mpp_derot_build(planet_body_t body, int system, double epoch_jd,
                             const double *jd, int num_frames,
                             int frame_rows, int frame_cols,
                             double cx, double cy, double radius,
                             double pa_deg, double parity,
                             double obs_lat, double obs_lon, double obs_elev);

/* Rough planet-disk detector. Thresholds the bright region (half-max above a
 * border-estimated sky background) and takes its flux centroid and the minor
 * axis of the brightness distribution; the equatorial radius is recovered from
 * that minor (≈ polar) axis using the body's known `flattening`. Using the
 * minor axis makes it robust to Saturn's rings, which inflate the major axis
 * but not the perpendicular extent. Fills centre/equatorial-radius (full-frame
 * pixels) and returns TRUE when a plausible disk was found. Pass flattening = 0
 * for a spherical estimate. When has_rings is TRUE (Saturn) the equatorial
 * radius is taken from the bright rings — their tips give the major axis, scaled
 * by the known A-ring-outer / equatorial-radius ratio — because the rings
 * corrupt the minor-axis estimate. If major_axis_deg is non-NULL it receives the
 * orientation of the bright region's major axis (Saturn's ring plane /
 * Jupiter's equator) in degrees, image coordinates with y up — the caller turns
 * this into the pole position angle (perpendicular). Used to seed the GUI disk
 * fit. */
gboolean mpp_derot_autodetect_disk(const fits *fit, double flattening,
                                   gboolean has_rings,
                                   double *cx, double *cy, double *radius,
                                   double *major_axis_deg);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_DEROT_BUILD_H_ */
