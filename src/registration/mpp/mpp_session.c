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

#include <stdlib.h>
#include <string.h>

#include "registration/mpp/mpp_session.h"
#include "registration/mpp/mpp_derot_build.h"   /* mpp_derot_union_epoch */

mpp_session_t *mpp_session_new(planet_body_t body, int system) {
	mpp_session_t *s = calloc(1, sizeof(*s));
	if (!s) return NULL;
	s->reference = -1;
	s->body = body;
	s->system = (system >= 1 && system <= 3) ? system : 1;
	return s;
}

void mpp_session_free(mpp_session_t *s) {
	if (!s) return;
	for (int i = 0; i < s->count; i++)
		free(s->seqs[i].seqname);
	free(s->seqs);
	free(s);
}

static int session_find(const mpp_session_t *s, const char *seqname) {
	for (int i = 0; i < s->count; i++)
		if (!g_strcmp0(s->seqs[i].seqname, seqname)) return i;
	return -1;
}

int mpp_session_add(mpp_session_t *s, const char *seqname, mpp_channel_t channel,
                    double first_jd, double last_jd, int num_frames) {
	if (!s || !seqname || !*seqname) return -1;
	if (channel < 0 || channel >= MPP_CHAN_COUNT) return -1;
	if (session_find(s, seqname) >= 0) return -1;   /* no duplicates */

	if (s->count >= s->capacity) {
		const int ncap = s->capacity ? s->capacity * 2 : 4;
		mpp_session_seq_t *ns = realloc(s->seqs, (size_t) ncap * sizeof(*ns));
		if (!ns) return -1;
		s->seqs = ns;
		s->capacity = ncap;
	}
	mpp_session_seq_t *e = &s->seqs[s->count];
	e->seqname = g_strdup(seqname);
	e->channel = channel;
	e->first_jd = first_jd;
	e->last_jd = last_jd;
	e->num_frames = num_frames;
	const int idx = s->count++;
	if (s->reference < 0) s->reference = idx;   /* first added is the reference */
	return idx;
}

gboolean mpp_session_remove(mpp_session_t *s, int index) {
	if (!s || index < 0 || index >= s->count) return FALSE;
	free(s->seqs[index].seqname);
	for (int i = index; i < s->count - 1; i++)
		s->seqs[i] = s->seqs[i + 1];
	s->count--;
	if (s->reference == index)
		s->reference = (s->count > 0) ? 0 : -1;
	else if (s->reference > index)
		s->reference--;
	return TRUE;
}

gboolean mpp_session_set_reference(mpp_session_t *s, int index) {
	if (!s || index < 0 || index >= s->count) return FALSE;
	s->reference = index;
	return TRUE;
}

gboolean mpp_session_set_channel(mpp_session_t *s, int index, mpp_channel_t channel) {
	if (!s || index < 0 || index >= s->count
	    || channel < 0 || channel >= MPP_CHAN_COUNT)
		return FALSE;
	s->seqs[index].channel = channel;
	return TRUE;
}

double mpp_session_compute_epoch(mpp_session_t *s) {
	if (!s || s->count <= 0) { if (s) s->epoch_jd = 0.0; return 0.0; }
	double *firsts = malloc((size_t) s->count * sizeof(double));
	double *lasts = malloc((size_t) s->count * sizeof(double));
	if (!firsts || !lasts) { free(firsts); free(lasts); return s->epoch_jd; }
	for (int i = 0; i < s->count; i++) {
		firsts[i] = s->seqs[i].first_jd;
		lasts[i] = s->seqs[i].last_jd;
	}
	s->epoch_jd = mpp_derot_union_epoch(firsts, lasts, s->count);
	free(firsts);
	free(lasts);
	return s->epoch_jd;
}

int mpp_session_channels(const mpp_session_t *s, mpp_channel_t *out) {
	if (!s || !out) return 0;
	/* report in a stable, display-friendly order */
	static const mpp_channel_t order[] = {
		MPP_CHAN_R, MPP_CHAN_G, MPP_CHAN_B, MPP_CHAN_LUM, MPP_CHAN_MONO
	};
	int n = 0;
	for (size_t k = 0; k < sizeof(order) / sizeof(order[0]); k++) {
		for (int i = 0; i < s->count; i++) {
			if (s->seqs[i].channel == order[k]) { out[n++] = order[k]; break; }
		}
	}
	return n;
}

int mpp_session_channel_seqs(const mpp_session_t *s, mpp_channel_t channel, int *out) {
	if (!s || !out) return 0;
	int n = 0;
	for (int i = 0; i < s->count; i++)
		if (s->seqs[i].channel == channel) out[n++] = i;
	return n;
}
