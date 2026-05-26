#include <math.h>
#include "line_detection.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "opencv/tracks.h"
#include "algos/sorting.h"
#include "algos/statistics_float.h"
#include "algos/statistics.h"
#include "io/image_format_fits.h"
#include "algos/photometry.h"	// for getMagnitude()

/* simple line detection interface, uses the opencv Hough lines algorithm from opencv/tracks.cpp in a
 * loop with decreasing detection thresholds
 */

static int compute_linear_threshold(fits *fit, int layer, int threads, int iteration, float min_allowed_ksigma, float *ksigma, float *threshold);
static int compute_knoise_threshold(fits *fit, int layer, int threads, float knoise, float *threshold);
static void confirm_streaks(struct track *tracks, int *nbtracks, fits *fit, int layer, struct streak_detector *conf, float ksigma);

/* this function provides an exponential decrease in the ksigma value to try, which limits the number of
 * tries compared to a linear decrease while concentrating the efforts on the low end of the range, where
 * small changes have large effects on noise. It is defined between 0 and 10.
 */
#define MAX_KSIGMA_IDX 10
static float ksigma_f(int x, float initial, float minimal) {
	g_assert(x >= 0);
	if (x >= MAX_KSIGMA_IDX)
		return minimal;
        return minimal + (initial - minimal) * expf(-0.50f*(float)x);
}

// no prior line detection
struct results *detect_streaks(fits *fit, int layer, struct streak_detector *conf, int nb_threads) {
	/* current strategy: start with a high k.sigma and length; decrease k.sigma down
	 * to some noise limit; if still not ok, lower the length and restart with high k.
	 * The bright streak mode decreases k with a linear law starting with saturated
	 * pixels, the normal mode decreases with a logarithmic law and starts closer to
	 * the median.
	 */
	if (conf->initial_length > 0 && conf->min_allowed_length > 0 &&
			conf->initial_length < conf->min_allowed_length) {
		siril_log_message("Line detection configuration invalid: min length too high\n");
		return NULL;
	}
	int min_length = conf->initial_length;
	if (min_length <= 0)
		min_length = 500;
	int ks_idx = 0;

	int min_allowed_length = conf->min_allowed_length;
	float min_allowed_ksigma = conf->min_allowed_ksigma;

	if (min_allowed_length <= 0)
		min_allowed_length = max(25, round_to_int(fit->rx * 0.005));
	if (min_allowed_ksigma <= 0.0f)
		min_allowed_ksigma = conf->bright_streak ? 3.0f : KSIGMA_MINIMAL_VALUE;	// that's the real limit
	float initial_ksigma = KSIGMA_INITIAL_FACTOR * min_allowed_ksigma;

	siril_log_debug("min allowed length: %d, min allowed_ksigma: %f (%s mode)\n",
			min_allowed_length, min_allowed_ksigma, conf->bright_streak ? "bright" : "normal");
	struct results *result = NULL;
	int loop = 1;
	WORD previous_threshold = 0;

	do {
		if (ks_idx > MAX_KSIGMA_IDX) {
			// we use a polynomial function to decrease the length, it could also be an
			// exponential
			// defined for lengths <= 500:
			// 500	0.5
			// 250	0.7
			// 100	0.9
			//  50	0.94
			//  40	0.95
			double factor;
			if (min_length > 500)
				factor = 0.4;
			else if (min_length < 40)
				factor = 0.90;
			else factor = 0.9628319 - min_length * 1.504368e-4 - (min_length * min_length) * 5.654642e-6 + (min_length * min_length * min_length) * 8.208338e-9;
			siril_log_debug("segment length multiplier: %f\n", factor);
			int new_min_length = round_to_int(factor * min_length);
			if (new_min_length < min_allowed_length && min_length != min_allowed_length)
				min_length = min_allowed_length;
			else min_length = new_min_length;
			ks_idx = 0;
		}
		if (min_length < min_allowed_length) {
			if (!result) {
				if (conf->can_recurse)
					siril_log_message("No detected line segment within length limits (%d)\n", min_allowed_length);
				else siril_log_debug("Could not confirm a streak with minimal length %d\n", min_allowed_length);
			}
			else siril_log_debug("end of the easy detection loop\n");
			break;
		}

		float threshold, ksigma;
		int retval = 0;
		if (conf->bright_streak) {
			// ksigma is used for the logs and return values here too
			compute_linear_threshold(fit, layer, nb_threads, ks_idx++,
					min_allowed_ksigma, &ksigma, &threshold);
		} else {
			ksigma = ksigma_f(ks_idx++, initial_ksigma, min_allowed_ksigma);
			retval = (compute_knoise_threshold(fit, layer, nb_threads, ksigma, &threshold));
		}
		if (retval < 0)	// stats error or threshold lower than median
			break;
		if (retval > 0)	// threshold higher than max value
			continue;
		if (fit->type == DATA_USHORT) {
			// avoid running the same detection if integer threshold doesn't change
			WORD current_threshold = (WORD)threshold;
			if (previous_threshold == current_threshold) {
				siril_log_debug("skipping detection for this threshold\n");
				continue;
			}
			previous_threshold = current_threshold;
		}
		siril_log_debug("iteration %d, min length: %d, ksigma %.2f (%.4f)\n", loop, min_length, ksigma, threshold);

		struct track *cur_tracks;
		int cur_nbtracks = cvHoughLines(fit, layer, threshold, min_length, &cur_tracks);
		siril_log_debug("iteration %d: %d line segments detected\n", loop, cur_nbtracks);
		if (cur_nbtracks < 0) {
			// too many pixels above the threshold
			min_allowed_ksigma = ksigma;
			if (ks_idx == 1) {
				siril_log_message("The number of pixels above the provided threshold is too high for this image, cannot process\n");
				break;
			}
			continue;
		}
		if (cur_nbtracks > conf->max_allowed_segments) {
			// revert to best
			if (!result)
				siril_log_message("Too many segments detected (%d)\n", cur_nbtracks);
			free(cur_tracks);
			break;
		}
		if (cur_nbtracks > 0) {
			// check that detected streaks are actual streaks, remove them if not
			confirm_streaks(cur_tracks, &cur_nbtracks, fit, layer, conf, ksigma);
			if (cur_nbtracks <= 0)
				break;

			// store current and continue, we may get a better solution
			if (cur_nbtracks > 0 && (!result || cur_nbtracks > result->nb_tracks)) {
				if (!result)
					result = calloc(1, sizeof(struct results));
				else free(result->tracks);
				result->nb_tracks = cur_nbtracks;
				result->tracks = cur_tracks;
				result->ksigma = ksigma;
				siril_log_debug("stored a solution for %d tracks\n", cur_nbtracks);
			}
		}
		if (cur_nbtracks >= conf->enough_segments) {
			siril_log_debug("Target reached, %d tracks (%d loops)\n", cur_nbtracks, loop);
			if (cur_nbtracks < 5 && conf->can_recurse) {
				for (int i = 0; i < cur_nbtracks; i++) {
					double dx = cur_tracks[i].start.x - cur_tracks[i].end.x;
					double dy = cur_tracks[i].start.y - cur_tracks[i].end.y;
					double len = sqrt(dx * dx + dy * dy);
					siril_log_message("Track found at (%d,%d -> %d,%d), length %d, angle %.1f°\n",
							cur_tracks[i].start.x, cur_tracks[i].start.y,
							cur_tracks[i].end.x, cur_tracks[i].end.y,
							(int)len, cur_tracks[i].angle);
				}
			}
			break;
		}
		loop++;
	} while (1);

	return result;
}

void clear_results(struct results *result) {
	if (!result) return;
	free(result->tracks);
	free(result);
}

// for bright targets
static int compute_linear_threshold(fits *fit, int layer, int threads, int iteration, float min_allowed_ksigma, float *ksigma, float *threshold) {
	// we use stats to check that we're not going too low if the image is bright
	imstats* stat = statistics(NULL, -1, fit, layer, NULL, STATS_BASIC, threads);
	if (!stat) {
		siril_log_error(_("Error: statistics computation failed.\n"));
		*threshold = 0.0f; *ksigma = 0.0f;
		return -1;
	}
	float factor = powf(0.66f, (float)iteration);
	float max = fit->type == DATA_USHORT ? (fit->bitpix == BYTE_IMG ? 254.0f : 65000.0f) : 0.99f;
	*threshold = max * factor;
	*ksigma = (*threshold - stat->median) / stat->bgnoise;
	int retval = 0;
	// sanity checks
	if (*threshold > stat->max) {
		siril_log_debug("Detection threshold %f is larger than max value %f.\n",
				*threshold, stat->max);
		retval = 1;
	}
	if (*threshold < stat->median + min_allowed_ksigma * stat->bgnoise) {
		siril_log_error(_("Detection threshold %f is lower than %.2f ksigma.\n"),
				*threshold, min_allowed_ksigma);
		retval = -1;
	}
	free_stats(stat);
	return retval;
}

// for regular targets (not `bright')
static int compute_knoise_threshold(fits *fit, int layer, int threads, float knoise, float *threshold) {
	imstats* stat = statistics(NULL, -1, fit, layer, NULL, STATS_BASIC, threads);
	if (!stat) {
		siril_log_error(_("Error: statistics computation failed.\n"));
		return -1;
	}
	*threshold = stat->median + knoise * stat->bgnoise;
	int retval = 0;
	// sanity checks
	if (*threshold > stat->max) {
		siril_log_debug("Detection threshold %f is larger than max value %f.\n",
				*threshold, stat->max);
		retval = 1;
	}
	if (*threshold < stat->median) {
		siril_log_error(_("Detection threshold %f is lower than median value %f.\n"),
				*threshold, stat->median);
		retval = -1;
	}
	free_stats(stat);
	return retval;
}

static void force_area_in_image(fits *fit, rectangle *a) {
	if (a->x < 0) a->x = 0;
	if (a->y < 0) a->y = 0;
	if (a->w > fit->rx) a->w = fit->rx;
	if (a->h > fit->ry) a->h = fit->ry;
	if (a->x + a->w > fit->rx) a->x = fit->rx - a->w;
	if (a->y + a->h > fit->ry) a->y = fit->ry - a->h;
}

static struct streak_detector *clone_ldconf(struct streak_detector *conf) {
	struct streak_detector *ret = malloc(sizeof(struct streak_detector));
	if (ret)
		memcpy(ret, conf, sizeof(struct streak_detector));
	return ret;
}

#define DEBUG_MASKS 0
#if DEBUG_MASKS
static void show_masks_on_image(fits *fit, BYTE *streak_mask, BYTE *bkg_mask, float streak_level) {
	int i = 0;
	if (streak_level < 0.0f)
		streak_level = 1.0f;
	else streak_level *= 1.5;
	if (fit->type == DATA_USHORT) {
		WORD level = roundf_to_WORD(streak_level * 65535.0f);
		siril_log_debug("mask annotation level: %d\n", level);
		for (int y = 0; y < fit->ry; y++) {
			for (int x = 0; x < fit->rx; x++) {
				if (streak_mask[i])
					fit->data[i] = fit->data[i] < level ? 0 : fit->data[i] - level;
				else if (bkg_mask[i])
					fit->data[i] = level >= 65535 - fit->data[i] ? 65535 : fit->data[i] + level;
				i++;
			}
		}
	}
	else if (fit->type == DATA_FLOAT) {
		siril_log_debug("mask annotation level: %.4f\n", streak_level);
		for (int y = 0; y < fit->ry; y++) {
			for (int x = 0; x < fit->rx; x++) {
				if (streak_mask[i])
					fit->fdata[i] -= streak_level;
				else if (bkg_mask[i])
					fit->fdata[i] += streak_level;
				i++;
			}
		}
	}
	savefits("annotated_streak.fit", fit);
}
#endif

static void streak_photometry(fits *fit, struct track *track, struct track *result, float fwhm, struct phot_config *pset) {
#if DEBUG_MASKS
	savefits("photometry.fit", fit);
#endif
	BYTE *streak_mask = NULL, *bkg_mask = NULL;
	if (create_streak_masks(fit->rx, fit->ry, track->start.x, track->start.y,
				track->end.x, track->end.y, fwhm, &streak_mask, &bkg_mask)) {
		siril_log_message("Error in flux computation (mask creation failed)\n");
		return;
	}

	int i = 0;
	double flux = 0.0;
	double *bkg = malloc(fit->rx * fit->ry * sizeof(double));
	int bkg_count = 0, ap_count = 0;
	gboolean saturated = FALSE;
	if (fit->type == DATA_USHORT) {
		double invnorm = fit->bitpix == BYTE_IMG ? INV_UCHAR_MAX_DOUBLE : INV_USHRT_MAX_DOUBLE;
		WORD sat = fit->bitpix == BYTE_IMG ?
			round_to_WORD(pset->maxval / INV_UCHAR_MAX_DOUBLE) :
			round_to_WORD(pset->maxval);
		for (int y = 0; y < fit->ry; y++) {
			for (int x = 0; x < fit->rx; x++) {
				if (streak_mask[i]) {
					WORD pixel = fit->data[i];
					flux += pixel * invnorm;
					ap_count++;
					if (pixel >= sat)
						saturated = TRUE;
				}
				else if (bkg_mask[i])
					bkg[bkg_count++] = fit->data[i] * invnorm;
				i++;
			}
		}
	}
	else if (fit->type == DATA_FLOAT) {
		for (int y = 0; y < fit->ry; y++) {
			for (int x = 0; x < fit->rx; x++) {
				if (streak_mask[i]) {
					float pixel = fit->fdata[i];
					flux += pixel;
					ap_count++;
					if (pixel >= pset->maxval)
						saturated = TRUE;
				}
				else if (bkg_mask[i])
					bkg[bkg_count++] = fit->fdata[i];
				// should we check for minval in the background?
				i++;
			}
		}
	}

	double background = 0.0, stdev = 0.0;
	if (bkg_count < 5 /* = MIN_SKY */|| robustmean(bkg_count, bkg, &background, &stdev)) {
		siril_log_message("failed to compute flux for streak (%d background samples)\n", bkg_count);
	} else {
		siril_log_debug("photometry background computed to %f (flux sum: %f)\n", background, flux);
		flux -= background * ap_count;
		if (flux < 0.0)
			siril_log_message("Error in flux computation (negative)\n");
		else {
			if (fit->type == DATA_USHORT)
				flux *= fit->bitpix == BYTE_IMG ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
			double mag = getMagnitude(flux);
			double snr;
			double mag_err = getMagErr(flux, ap_count, bkg_count, stdev, pset->gain, &snr);
			siril_log_message("streak mag%s: %.3f (%.3f corrected) +/- %.3f, (flux %f from %d pixels)\n",
					saturated? " (saturated)" : "",
					mag, mag+com.magOffset, mag_err, flux, ap_count);
			result->mag = mag+com.magOffset;
			result->mag_err = mag_err;
			result->mag_is_absolute = com.magOffset != 0.0;
			result->mag_is_accurate = !saturated;
			result->snr = snr;
		}
	}
	// see also getPhotometryData() for star flux computation

#if DEBUG_MASKS
	// produce an image showing the masks
	show_masks_on_image(fit, streak_mask, bkg_mask, flux / (float)ap_count);
#endif
	free(streak_mask);
	free(bkg_mask);
	free(bkg);
}

/* look closely at detected coordinates to verify there is actually a streak, also measure flux */
static gboolean confirm_streak(fits *fit, int layer, struct track *track, struct streak_detector *conf, struct phot_config *pset) {
	int sx = track->start.x, sy = track->start.y;
	int ex = track->end.x, ey = track->end.y;
	int dx = ex - sx, dy = ey - sy;
	double middle_x = (sx + ex) * 0.5;
	double middle_y = (sy + ey) * 0.5;
	double length = sqrt(dx * dx + dy * dy);
	double box_half_sz = length * 1.5;
	if (conf->bright_streak || (fit->rx <= 2.0 * box_half_sz && fit->ry <= 2.0 * box_half_sz)) {
		siril_log_debug("Streak doesn't need to be confirmed\n");
		if (conf->compute_flux)
			streak_photometry(fit, track, track,
					conf->bright_streak ? 5.0f * conf->fwhm : conf->fwhm, pset);
		return TRUE;
	}
	rectangle surroundings = {
		.x = middle_x - box_half_sz, .y = middle_y - box_half_sz,
		.w = 2.0 * box_half_sz, .h = 2.0 * box_half_sz
	};
	force_area_in_image(fit, &surroundings);
	fits newfit = { 0 };
	extract_region_from_fits(fit, layer, &newfit, &surroundings);

	conf->initial_length = length;
	conf->min_allowed_length = max(25, length * 0.65);
	struct results *newstreaks = detect_streaks(&newfit, layer, conf, SINGLE_THREADED);
	gboolean retval = FALSE;
	if (newstreaks) {
		if (newstreaks->nb_tracks > 0) {
			siril_log_debug("Confirmed streak at (%d,%d -> %d,%d), angle %.1f°\n",
					newstreaks->tracks[0].start.x + surroundings.x,
					newstreaks->tracks[0].start.y + surroundings.y,
					newstreaks->tracks[0].end.x + surroundings.x,
					newstreaks->tracks[0].end.y + surroundings.y,
					newstreaks->tracks[0].angle);
			// TODO: compare results (with region offset)
			retval = TRUE;
		}

		/* streak photometry: computing the flux only */
		if (retval && conf->compute_flux) {
			streak_photometry(&newfit, &(newstreaks->tracks[0]), track, conf->fwhm, pset);
		}
		clear_results(newstreaks);
	}
	clearfits(&newfit);
	return retval;
}

static void confirm_streaks(struct track *tracks, int *nbtracks, fits *fit, int layer, struct streak_detector *conf, float ksigma) {
	if (!conf->can_recurse) return;
	int i = 0, j = 0;
	struct streak_detector *conf2 = clone_ldconf(conf);
	conf2->can_recurse = FALSE;
	conf2->enough_segments = 1;
	conf2->min_allowed_ksigma = ksigma * 0.48f;
	struct phot_config *pset = phot_set_adjusted_for_image(fit);
	siril_log_debug("confirming streaks\n");
	while (i < *nbtracks) {
		struct track *track = &tracks[i++];
		if (!confirm_streak(fit, layer, track, conf2, pset))
			continue;
		if (i != j)
			memcpy(&tracks[j], track, sizeof(struct track));
		j++;
	}
	*nbtracks = j;
	free(pset);
	free(conf2);
}

