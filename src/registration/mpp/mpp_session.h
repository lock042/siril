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
 * Multi-sequence derotation session (Option B). A pure data model: the set of
 * sequences to combine, each tagged with a colour channel and its capture span,
 * one designated reference sequence shared across all channels, the shared body
 * and the reference epoch. The engine and GUI build on this; it does no I/O.
 */

#ifndef SRC_REGISTRATION_MPP_MPP_SESSION_H_
#define SRC_REGISTRATION_MPP_MPP_SESSION_H_

#include <glib.h>
#include "algos/planet_ephem.h"   /* planet_body_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
	MPP_CHAN_MONO = 0,   /* untagged: a single-channel output */
	MPP_CHAN_R,
	MPP_CHAN_G,
	MPP_CHAN_B,
	MPP_CHAN_LUM,
	MPP_CHAN_OSC,        /* one-shot colour (CFA): stacks to a 3-layer result */
	MPP_CHAN_COUNT
} mpp_channel_t;

typedef struct {
	char *seqname;            /* sequence base name (owned) */
	mpp_channel_t channel;
	double first_jd, last_jd; /* capture span, UTC Julian dates (caller-filled) */
	int num_frames;

	/* Disk fit in this sequence's own frame pixels (filled via
	 * mpp_session_set_fit once the user has fitted the planet). r_pol <= 0
	 * means "use the body's default flattening"; fps > 0 supplies a uniform
	 * cadence for sequences without per-frame timestamps. */
	gboolean has_fit;
	double cx, cy, r_eq, r_pol, pa_deg, parity;
	int frame_rows, frame_cols;
	int bayer;                /* CFA signature for compatibility (SER color_id; -1 = none) */
	double fps;
} mpp_session_seq_t;

typedef struct {
	mpp_session_seq_t *seqs;
	int count, capacity;
	int reference;            /* index of the reference sequence, -1 if none */
	planet_body_t body;
	int system;               /* rotation system 1..3 */
	double epoch_jd;          /* shared reference epoch (0 = not computed yet) */
} mpp_session_t;

mpp_session_t *mpp_session_new(planet_body_t body, int system);
void mpp_session_free(mpp_session_t *s);

/* Add a sequence with its channel and capture span (first/last frame JD, from
 * mpp_derot_sequence_span). Returns the new index, or -1 on error or duplicate
 * name. The first sequence added becomes the reference by default. */
int mpp_session_add(mpp_session_t *s, const char *seqname, mpp_channel_t channel,
                    double first_jd, double last_jd, int num_frames);

gboolean mpp_session_remove(mpp_session_t *s, int index);
gboolean mpp_session_set_reference(mpp_session_t *s, int index);
gboolean mpp_session_set_channel(mpp_session_t *s, int index, mpp_channel_t channel);

/* Record the disk fit (and frame dims / cadence) for a sequence already added.
 * Marks the entry as fitted. */
gboolean mpp_session_set_fit(mpp_session_t *s, int index,
                             double cx, double cy, double r_eq, double r_pol,
                             double pa_deg, double parity,
                             int frame_rows, int frame_cols, int bayer, double fps);

/* Shared epoch = midpoint of the union of every sequence's span. Stores it on
 * the session and returns it (0 when the session is empty). */
double mpp_session_compute_epoch(mpp_session_t *s);

/* Distinct channels present, in R, G, B, LUM, MONO order; fills out (size
 * MPP_CHAN_COUNT) and returns the count. */
int mpp_session_channels(const mpp_session_t *s, mpp_channel_t *out);

/* Indices of the sequences in a given channel; fills out (size >= s->count) and
 * returns how many. */
int mpp_session_channel_seqs(const mpp_session_t *s, mpp_channel_t channel, int *out);

#ifdef __cplusplus
}
#endif

#endif /* SRC_REGISTRATION_MPP_MPP_SESSION_H_ */
