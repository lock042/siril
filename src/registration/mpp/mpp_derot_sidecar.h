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

#ifndef SRC_REGISTRATION_MPP_DEROT_SIDECAR_H_
#define SRC_REGISTRATION_MPP_DEROT_SIDECAR_H_

#include <stdint.h>
#include "registration/mpp.h"   /* mpp_status_t */

#ifdef __cplusplus
extern "C" {
#endif

/* Per-sequence derotation plan, written as <seqname>.derot by the derotation
 * tool (or the `derotate` command) and consumed by the stack stage. It is
 * deliberately independent of the .mpp sidecar so re-running registration
 * (which rewrites .mpp) never clobbers it.
 *
 * The per-frame disk geometry (B, CM, P) is *resolved* at compute time from the
 * built-in ephemeris and stored, so stacking is reproducible regardless of
 * later ephemeris-model changes or timestamp re-reads. Angles are degrees.
 *
 * Pole-angle handling: `pole_angle_epoch` is the fitted position angle of the
 * projected north pole *in the image* at the epoch (the disk-fit orientation,
 * radians). A frame's in-image pole angle is pole_angle_epoch + (pole_pa[i] -
 * epoch_pole_pa) — the ephemeris pole PA drifts only slowly across a derotation
 * span. */
typedef struct mpp_derot {
	int32_t body;            /* planet_body_t */
	int32_t rot_system;      /* rot_system_t: which CM column drives derotation */
	int32_t ephem_version;   /* ephemeris model version, for provenance */
	int32_t num_frames;
	int32_t frame_rows;      /* sequence frame height the disk fit is in */
	int32_t frame_cols;      /* sequence frame width  the disk fit is in */

	double epoch_jd;         /* reference epoch, UTC Julian date */
	double obs_lat, obs_lon, obs_elev;   /* observer; NAN => geocentric */

	/* disk fit (reference/epoch canvas pixels) */
	double cx, cy, r_eq, flattening, parity;
	double pole_angle_epoch; /* theta0, radians (in-image pole angle) */

	/* epoch geometry (degrees) */
	double epoch_sub_obs_lat, epoch_cm, epoch_pole_pa;

	/* per-frame resolved geometry, each length num_frames (degrees / JD) */
	double *jd;
	double *sub_obs_lat;
	double *cm;
	double *pole_pa;
} mpp_derot_t;

/* Allocate a record with the four per-frame arrays sized to num_frames
 * (zero-filled), parity defaulted to +1. NULL on allocation failure. */
mpp_derot_t *mpp_derot_alloc(int num_frames);
void mpp_derot_free(mpp_derot_t *d);

/* Binary I/O (magic "SIRILDRT"). Round-trips every field; strict magic/version
 * check, short reads -> MPP_EIO. mpp_derot_read allocates *out (caller frees
 * with mpp_derot_free). */
mpp_status_t mpp_derot_write(const char *path, const mpp_derot_t *d);
mpp_status_t mpp_derot_read(const char *path, mpp_derot_t **out);

/* Stable 64-bit fingerprint over everything in the plan that affects the
 * derotation maps (body, system, epoch, disk fit, observer, and the per-frame
 * geometry). Registration stamps this into the .mpp so Stage C can detect a
 * .derot that changed in ANY way since the AP shifts were measured — not just a
 * different epoch. Never returns 0 (reserved as the "registered without
 * derotation" sentinel). Returns 0 only for a NULL plan. */
uint64_t mpp_derot_fingerprint(const mpp_derot_t *d);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_DEROT_SIDECAR_H_ */
