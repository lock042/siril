#ifndef _LINE_DETECTION_H
#define _LINE_DETECTION_H

#include <glib.h>
#include "core/siril.h"
#include "opencv/tracks.h"

struct streak_detector {
	int min_allowed_length;		// negative for unset
	float min_allowed_ksigma;	// negative for unset
	int initial_length;		// expected segment length (<= 0 for default)
	gboolean bright_streak;		// use linear pixel thresholds and larger box

	int max_allowed_segments;	// should always be set, an amount that would be unusual
	int enough_segments;		// should always be set, stop if we find that many
	gboolean can_recurse;		// true (default) enables candidate confirmation
	gboolean compute_flux;		// true for streak photometry, done on confirmation
	float fwhm;			// for the width of the streak, for photometry
};

struct results {
	struct track *tracks;
	int nb_tracks;
	float ksigma;		// the ksigma at which detection for this result was made
};

/* the initial sigma coefficient will be this number times the minimum value (min_allowed_ksigma) */
#define KSIGMA_INITIAL_FACTOR 3.1f
#define KSIGMA_MINIMAL_VALUE 2.4f

struct results *detect_streaks(fits *fit, int layer, struct streak_detector *conf, int nb_threads);
void clear_results(struct results *result);

#endif
