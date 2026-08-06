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
#ifndef _NDE_CAT_H_
#define _NDE_CAT_H_

/**
 * \file nde_cat.h
 * \brief Embedded star catalogues for the recompute-at-replay calibrations.
 *
 * A PCC/SPCC record's INTENT can only be recomputed offline if the star data
 * the original analysis used is still available: B/V magnitudes (or Teff for
 * the Gaia catalogues) for PCC, the full sampled flux spectra for SPCC.  At
 * analysis time all of it sits in memory in the pipeline's ref_stars — so it
 * is stashed right there, adopted by the NDE capture under the new record's
 * id, persisted into the .flis as a FLIS_NDE_CAT binary-table HDU (the file
 * is then fully self-contained and portable — user decision 2026-08-05), and
 * re-registered on load.  Replay re-projects the stored sky positions with
 * the replayed image's WCS and runs the ordinary pipeline, no network.
 *
 * The registry is the in-session half: record_id → deep-copied catalogue.
 * The HDU I/O lives with the rest of the FLIS format (image_format_flis.c).
 *
 * Threading: leaf mutex; all entry points callable from any thread.
 */

#include <glib.h>
#include "io/siril_catalogues.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Deep-copy the fields replay needs (identity, photometry, spectra; not the
 *  projection, which is re-derived).  Caller siril_catalog_free()s. */
siril_catalogue *nde_cat_copy(const siril_catalogue *cat);

/**
 * Stash a copy of @cat as the PENDING catalogue of the capture about to
 * happen.  Called by the photometric pipeline right after a successful
 * analysis (never during replay).  The very next NDE capture of a
 * photometric op adopts it; any other capture discards it as stale.
 */
void nde_cat_stash_pending(const siril_catalogue *cat);

/** Adopt or discard the pending catalogue — the capture side (nde_history.c).
 *  Adopts under @record_id when the captured op is one of the photometric
 *  ones (@op_id decides); discards otherwise. */
void nde_cat_adopt_pending(gint64 record_id, const char *op_id);

/** Register @cat (ownership transferred) under @record_id — the load path. */
void nde_cat_register(gint64 record_id, siril_catalogue *cat);

/** A deep copy of @record_id's catalogue, or NULL.  Caller frees. */
siril_catalogue *nde_cat_get_copy(gint64 record_id);

/** Borrowed pointer for the save path (single-threaded there), or NULL. */
siril_catalogue *nde_cat_peek(gint64 record_id);

gboolean nde_cat_has(gint64 record_id);

/** Drop one entry / everything (document change). */
void nde_cat_drop(gint64 record_id);
void nde_cat_purge(void);

#ifdef __cplusplus
}
#endif

#endif /* _NDE_CAT_H_ */
