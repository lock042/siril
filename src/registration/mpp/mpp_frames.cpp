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
 *
 * This implementation of multipoint registration & stacking is based on
 * PlanetarySystemStacker by Rolf Hempel:
 *     https://github.com/Rolf-Hempel/PlanetarySystemStacker
 */

/*
 * Bridge between Siril's sequence/fits machinery and mpp's cv::Mat
 * pipeline. Phase 6.1.
 *
 * The headline reason this exists rather than calling seq_read_frame
 * directly is the debayer-on-open behaviour. Siril's seq_read_frame
 * routes SER debayering through the global `com.pref.debayer.open_debayer`
 * preference which we don't want to mutate (it's user-visible and not
 * thread-safe under our parallel reads). We dispatch the same way
 * seq_read_frame does but pass our own debayer flag for SER reads:
 * stacking reads force RCD debayer on; analysis reads leave the decision
 * to ser_file->debayer_type_ser (which the mpp stages pin to raw via
 * SerAnalysisRawGuard and demosaic with cv::cvtColor instead — far
 * cheaper in transient memory and CPU for a luminance-only consumer).
 */

#include <opencv2/core.hpp>

/* Siril headers don't have `extern "C"` guards; wrap them so the C symbols
 * link from a C++ TU. core/siril.h goes first (outside the wrap) because
 * it pulls in <omp.h>, whose libgomp 15 headers declare C++ templates
 * that can't sit inside an extern "C" block. The header guard then
 * suppresses the re-include from inside the wrap. */
#include "core/siril.h"

extern "C" {
#include "io/sequence.h"
#include "io/ser.h"
#include "io/fits_sequence.h"
#include "io/image_format_fits.h"
#ifdef HAVE_FFMS2
#include "io/films.h"
#endif

#include "registration/mpp.h"
#include "registration/mpp/mpp_frames.h"
}

extern "C" int mpp_seq_read_frame(sequence *seq, int idx, fits *dest,
                                  bool force_float, int thread_id,
                                  bool ser_force_debayer) {
	if (!seq || !dest || idx < 0 || idx >= seq->number)
		return MPP_EINVAL;

	const gboolean ff = force_float ? TRUE : FALSE;
	switch (seq->type) {
		case SEQ_REGULAR: {
			char filename[256];
			fit_sequence_get_image_filename_checkext(seq, idx, filename);
			if (readfits(filename, dest, NULL, ff))
				return MPP_EIO;
			break;
		}
		case SEQ_SER:
			if (!seq->ser_file) return MPP_EINVAL;
			/* TRUE forces debayer-on-read regardless of the user's global
			 * pref (ser_read_frame honours the flag); FALSE follows
			 * ser_file->debayer_type_ser — see the header for the analysis
			 * raw-read contract. */
			if (ser_read_frame(seq->ser_file, idx, dest, ff,
			                   ser_force_debayer ? TRUE : FALSE))
				return MPP_EIO;
			break;
		case SEQ_FITSEQ:
			if (!seq->fitseq_file) return MPP_EINVAL;
			if (fitseq_read_frame(seq->fitseq_file, idx, dest, ff, thread_id))
				return MPP_EIO;
			break;
#ifdef HAVE_FFMS2
		case SEQ_AVI:
			if (!seq->film_file) return MPP_EINVAL;
			if (film_read_frame(seq->film_file, idx, dest))
				return MPP_EIO;
			break;
#endif
		case SEQ_INTERNAL:
			/* In-memory sequence — used by tests. Defer to Siril's helper. */
			return seq_read_frame(seq, idx, dest, ff, thread_id) == 0
			    ? MPP_OK : MPP_EIO;
		default:
			return MPP_EINVAL;
	}
	return MPP_OK;
}
