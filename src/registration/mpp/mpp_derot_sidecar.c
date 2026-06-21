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
 * Binary sidecar (`.derot`) read/write for an mpp_derot_t.
 *
 * Layout (little-endian, explicit per-field serialisation):
 *
 *   8 bytes : magic "SIRILDRT"
 *   uint32  : version (current = 1)
 *   uint32  : flags (reserved; 0)
 *   int32[6]: body, rot_system, ephem_version, num_frames, frame_rows, frame_cols
 *   double[14]: epoch_jd, obs_lat, obs_lon, obs_elev,
 *               cx, cy, r_eq, flattening, parity, pole_angle_epoch,
 *               epoch_sub_obs_lat, epoch_cm, epoch_pole_pa, (reserved 0)
 *   double[num_frames] : jd
 *   double[num_frames] : sub_obs_lat
 *   double[num_frames] : cm
 *   double[num_frames] : pole_pa
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpp_derot_sidecar.h"

#define DEROT_MAGIC   "SIRILDRT"
#define DEROT_VERSION 2u   /* v2 added frame_rows/frame_cols */

static int write_all(FILE *f, const void *p, size_t n) {
	return fwrite(p, 1, n, f) == n ? 0 : -1;
}
static int read_all(FILE *f, void *p, size_t n) {
	return fread(p, 1, n, f) == n ? 0 : -1;
}

/* FNV-1a 64-bit over a byte range, folded into the running hash. */
static void fnv1a(uint64_t *h, const void *p, size_t n) {
	const unsigned char *b = (const unsigned char *) p;
	for (size_t i = 0; i < n; i++) {
		*h ^= b[i];
		*h *= 1099511628211ULL;
	}
}

uint64_t mpp_derot_fingerprint(const mpp_derot_t *d) {
	if (!d) return 0;
	uint64_t h = 1469598103934665603ULL;  /* FNV offset basis */
	const int32_t ints[6] = { d->body, d->rot_system, d->ephem_version,
	                          d->num_frames, d->frame_rows, d->frame_cols };
	fnv1a(&h, ints, sizeof(ints));
	const double scal[13] = {
		d->epoch_jd, d->obs_lat, d->obs_lon, d->obs_elev,
		d->cx, d->cy, d->r_eq, d->flattening, d->parity, d->pole_angle_epoch,
		d->epoch_sub_obs_lat, d->epoch_cm, d->epoch_pole_pa
	};
	fnv1a(&h, scal, sizeof(scal));
	if (d->num_frames > 0) {
		const size_t nb = (size_t) d->num_frames * sizeof(double);
		fnv1a(&h, d->jd, nb);
		fnv1a(&h, d->sub_obs_lat, nb);
		fnv1a(&h, d->cm, nb);
		fnv1a(&h, d->pole_pa, nb);
	}
	return h ? h : 1;   /* 0 is the "no derotation" sentinel */
}

mpp_derot_t *mpp_derot_alloc(int num_frames) {
	if (num_frames < 0) return NULL;
	mpp_derot_t *d = calloc(1, sizeof(*d));
	if (!d) return NULL;
	d->num_frames = num_frames;
	d->parity = 1.0;
	if (num_frames > 0) {
		d->jd          = calloc((size_t) num_frames, sizeof(double));
		d->sub_obs_lat = calloc((size_t) num_frames, sizeof(double));
		d->cm          = calloc((size_t) num_frames, sizeof(double));
		d->pole_pa     = calloc((size_t) num_frames, sizeof(double));
		if (!d->jd || !d->sub_obs_lat || !d->cm || !d->pole_pa) {
			mpp_derot_free(d);
			return NULL;
		}
	}
	return d;
}

void mpp_derot_free(mpp_derot_t *d) {
	if (!d) return;
	free(d->jd);
	free(d->sub_obs_lat);
	free(d->cm);
	free(d->pole_pa);
	free(d);
}

mpp_status_t mpp_derot_write(const char *path, const mpp_derot_t *d) {
	if (!path || !d || d->num_frames < 0) return MPP_EINVAL;
	if (d->num_frames > 0 && (!d->jd || !d->sub_obs_lat || !d->cm || !d->pole_pa))
		return MPP_EINVAL;

	FILE *f = fopen(path, "wb");
	if (!f) return MPP_EIO;

	const uint32_t ver = DEROT_VERSION, flags = 0;
	const int32_t hdr[6] = { d->body, d->rot_system, d->ephem_version,
	                         d->num_frames, d->frame_rows, d->frame_cols };
	const double dbl[14] = {
		d->epoch_jd, d->obs_lat, d->obs_lon, d->obs_elev,
		d->cx, d->cy, d->r_eq, d->flattening, d->parity, d->pole_angle_epoch,
		d->epoch_sub_obs_lat, d->epoch_cm, d->epoch_pole_pa, 0.0
	};
	const size_t nb = (size_t) d->num_frames * sizeof(double);

	if (write_all(f, DEROT_MAGIC, 8) != 0) goto err;
	if (write_all(f, &ver, 4) != 0) goto err;
	if (write_all(f, &flags, 4) != 0) goto err;
	if (write_all(f, hdr, sizeof(hdr)) != 0) goto err;
	if (write_all(f, dbl, sizeof(dbl)) != 0) goto err;
	if (d->num_frames > 0) {
		if (write_all(f, d->jd, nb) != 0) goto err;
		if (write_all(f, d->sub_obs_lat, nb) != 0) goto err;
		if (write_all(f, d->cm, nb) != 0) goto err;
		if (write_all(f, d->pole_pa, nb) != 0) goto err;
	}
	if (fclose(f) != 0) return MPP_EIO;
	return MPP_OK;
err:
	fclose(f);
	return MPP_EIO;
}

mpp_status_t mpp_derot_read(const char *path, mpp_derot_t **out) {
	if (!path || !out) return MPP_EINVAL;
	*out = NULL;

	FILE *f = fopen(path, "rb");
	if (!f) return MPP_EIO;

	char magic[8];
	uint32_t ver = 0, flags = 0;
	int32_t hdr[6];
	double dbl[14];
	if (read_all(f, magic, 8) != 0 || memcmp(magic, DEROT_MAGIC, 8) != 0) goto eio;
	if (read_all(f, &ver, 4) != 0 || ver != DEROT_VERSION) goto eio;
	if (read_all(f, &flags, 4) != 0) goto eio;
	if (read_all(f, hdr, sizeof(hdr)) != 0) goto eio;
	if (read_all(f, dbl, sizeof(dbl)) != 0) goto eio;
	if (hdr[3] < 0) goto eio;   /* num_frames */

	mpp_derot_t *d = mpp_derot_alloc(hdr[3]);
	if (!d) { fclose(f); return MPP_ENOMEM; }
	d->body = hdr[0]; d->rot_system = hdr[1]; d->ephem_version = hdr[2];
	d->frame_rows = hdr[4]; d->frame_cols = hdr[5];
	d->epoch_jd = dbl[0]; d->obs_lat = dbl[1]; d->obs_lon = dbl[2]; d->obs_elev = dbl[3];
	d->cx = dbl[4]; d->cy = dbl[5]; d->r_eq = dbl[6]; d->flattening = dbl[7];
	d->parity = dbl[8]; d->pole_angle_epoch = dbl[9];
	d->epoch_sub_obs_lat = dbl[10]; d->epoch_cm = dbl[11]; d->epoch_pole_pa = dbl[12];

	const size_t nb = (size_t) d->num_frames * sizeof(double);
	if (d->num_frames > 0) {
		if (read_all(f, d->jd, nb) != 0 ||
		    read_all(f, d->sub_obs_lat, nb) != 0 ||
		    read_all(f, d->cm, nb) != 0 ||
		    read_all(f, d->pole_pa, nb) != 0) {
			mpp_derot_free(d);
			fclose(f);
			return MPP_EIO;
		}
	}
	fclose(f);
	*out = d;
	return MPP_OK;
eio:
	fclose(f);
	return MPP_EIO;
}
