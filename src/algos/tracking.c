#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/siril_date.h"
#include "opencv/tracks.h"
#include "algos/astrometry_solver.h"
#include "core/siril_world_cs.h"
#include "gui/image_display.h"
#include "algos/statistics.h"
#include "algos/statistics_float.h"
#include "algos/siril_wcs.h"
#include "algos/tracking.h"
#include "algos/sorting.h"

/* 1. find stars, if there are 1 or a few, continue
 * 2. find lines, if there are 10 at least, continue
 * 3. sort lines by angle to keep only stars
 * 4. plate solve for line start
 * 5. get object coordinates for line start
 * 6. plate solve for line end
 * 7. get object coordinates for line end
 */

static int plate_solve_for_stars(fits *fit, psf_star **stars) {
	struct astrometry_data *args = NULL;
	args = calloc(1, sizeof(struct astrometry_data));
	args->fit = fit;
	args->pixel_size = max(fit->pixel_size_x, fit->pixel_size_y);
	args->focal_length = 1500.0;	// FIXME
	args->use_local_cat = TRUE;
	args->onlineCatalog = NOMAD;
	args->for_photometry_cc = FALSE;
	args->already_in_a_thread = TRUE;
	args->cat_center = get_eqs_from_header(fit);
	args->downsample = FALSE;
	args->autocrop = FALSE;
	args->flip_image = FALSE;
	args->manual = TRUE;
	args->stars = stars;
	args->auto_magnitude = TRUE;
	args->pcc = NULL;
	process_plate_solver_input(args);
	return GPOINTER_TO_INT(match_catalog(args));
}

struct tracked_object {
	psf_star *fixed_target;
	struct track *moving_track;
	SirilWorldCS *start, *end;
	double max_motion;
};

struct tracked_object *get_moving_targets(struct track *tracks, int nbtracks, float mean, int *nb) {
	*nb = 0;
	struct tracked_object *tracked = calloc(nbtracks, sizeof(struct tracked_object));
	for (int i = 0; i < nbtracks; i++) {
		/* excluding tracks too close from the bunch of trails, and
		 * vertical ones which are often just bright pixel columns */
		if (fabsf(mean - tracks[i].angle) > 20.0 && tracks[i].angle != -90.0) {
			tracked[*nb].moving_track = &tracks[i];
			(*nb)++;
		}
	}
	if (*nb > 0)
		tracked = realloc(tracked, *nb * sizeof(struct tracked_object));
	else {
		free(tracked);
		tracked = NULL;
	}
	siril_log_message(_("Considering %d moving targets in the image\n"), *nb);
	return tracked;
}

psf_star **extract_stars(struct track *tracks, int nbtracks, gboolean start, float mean, psf_star **stars) {
	/* we have a list of tracks, some may not be stars, only keep the tracks that have the
	 * majority of similar orientation */
	int j = 0;
	if (!stars && start) {
		// assuming stars is NULL
		stars = malloc((nbtracks + 1) * sizeof(psf_star *));
		for (int i = 0; i < nbtracks; i++) {
			if (fabsf(mean - tracks[i].angle) > 10.0) {
				/* consider the track as a rogue object to track */

				continue;
			}
			stars[j] = new_psf_star();
			stars[j]->xpos = tracks[i].start.x;
			stars[j]->ypos = tracks[i].start.y;
			stars[j]->fwhmx = 5.0;
			stars[j]->has_saturated = FALSE;
			j++;
		}
		stars[j] = NULL;
	} else if (stars && !start) {
		// assuming stars is non-NULL
		for (int i = 0; i < nbtracks; i++) {
			if (fabsf(mean - tracks[i].angle) > 10.0)
				continue;
			stars[j]->xpos = tracks[i].end.x;
			stars[j]->ypos = tracks[i].end.y;
			j++;
		}
	}
	else {
		siril_log_message("ERROR\n");
		return NULL;
	}
	return stars;
}

/* call free_star_list() before calling that */
gpointer tracking_worker(gpointer ptr) {
	struct linetrack_conf *arg = (struct linetrack_conf *)ptr;

	int retval = 0;
	struct track *tracks;
	int nblines = cvHoughLines(arg->fit, arg->layer, arg->threshold, arg->minlen, &tracks);
	if (nblines > 0 && nblines < 200) {
		siril_log_message(_("Found %d trail(s) in current frame, displaying start points\n"), nblines);
		if (nblines > 1000) nblines = 1000;
		if (arg->display_lines) {
			com.stars = malloc((2 * nblines + 1) * sizeof(psf_star *));
			for (int i = 0; i < nblines; i++) {
				com.stars[2*i] = new_psf_star();
				com.stars[2*i]->xpos = tracks[i].start.x;
				com.stars[2*i]->ypos = tracks[i].start.y;
				com.stars[2*i]->fwhmx = 5.0;
				com.stars[2*i]->has_saturated = FALSE;
				com.stars[2*i+1] = new_psf_star();
				com.stars[2*i+1]->xpos = tracks[i].end.x;
				com.stars[2*i+1]->ypos = tracks[i].end.y;
				com.stars[2*i+1]->fwhmx = 5.0;
				com.stars[2*i+1]->has_saturated = TRUE;
			}
			com.stars[2*nblines] = NULL;
			redraw(REDRAW_OVERLAY);
		}

		// compute main angle
		float *angles = malloc(nblines * sizeof(float));
		for (int i = 0; i < nblines; i++)
			angles[i] = tracks[i].angle;
		quicksort_f(angles, nblines);
		float mean = siril_stats_robust_mean(angles, 1, nblines, NULL);
		free(angles);
		siril_log_message(_("Considering trails with angle close to %f for stars\n"), mean);

		// extract moving targets from the list of tracks
		int nb_moving_targets;
		struct tracked_object *moving_targets = get_moving_targets(tracks, nblines, mean, &nb_moving_targets);

		int nb_targets = arg->nb_fixed_targets + nb_moving_targets;
		struct tracked_object *targets = calloc(nb_targets, sizeof(struct tracked_object));
		for (int i = 0; i < arg->nb_fixed_targets; i++)
			targets[i].fixed_target = arg->fixed_targets[i];
		for (int i = 0; i < nb_moving_targets; i++)
			memcpy(targets + (i + arg->nb_fixed_targets), &moving_targets[i], sizeof(struct tracked_object));

		gboolean start_solve_failed = FALSE, end_solve_failed = FALSE;
		psf_star **stars = extract_stars(tracks, nblines, TRUE, mean, NULL);
		if (plate_solve_for_stars(arg->fit, stars)) {
			siril_log_message(_("Plate solving using the detected trails as stars failed, adjust parameters or make sure the file metadata is correct\n"));
			start_solve_failed = TRUE;
		} else {
			for (int i = 0; i < nb_targets; i++) {
				double ra, dec;
				if (targets[i].fixed_target)
					pix2wcs(arg->fit, targets[i].fixed_target->xpos, targets[i].fixed_target->ypos, &ra, &dec);
				else {
					// using start of segment here, arbitrarily
					pix2wcs(arg->fit, targets[i].moving_track->start.x, targets[i].moving_track->start.y, &ra, &dec);
				}
				targets[i].start = siril_world_cs_new_from_a_d(ra, dec);
			}
			siril_log_message(_("Plate solving for start of movement succeeded.\n"));
		}

		stars = extract_stars(tracks, nblines, FALSE, mean, stars);
		if (plate_solve_for_stars(arg->fit, stars)) {
			siril_log_message(_("Plate solving using the detected trails as stars failed, adjust parameters or make sure the file metadata is correct\n"));
			end_solve_failed = TRUE;
		} else {
			for (int i = 0; i < nb_targets; i++) {
				double ra, dec;
				if (targets[i].fixed_target) {
					pix2wcs(arg->fit, targets[i].fixed_target->xpos, targets[i].fixed_target->ypos, &ra, &dec);
				} else {
					// using end of segment
					pix2wcs(arg->fit, targets[i].moving_track->end.x, targets[i].moving_track->end.y, &ra, &dec);
				}
				targets[i].end = siril_world_cs_new_from_a_d(ra, dec);
			}
		}
		free_fitted_stars(stars);

		/* TODO: estimate the width of tracks, which would give some error margin for
		 * the position of the targets, in addition to their FWHM and elongation */
		if (!start_solve_failed && !end_solve_failed) {
			for (int i = 0; i < nb_targets; i++) {
				gchar *start_ra = siril_world_cs_alpha_format(targets[i].start, " %02dh%02dm%02ds");
				gchar *start_dec = siril_world_cs_delta_format(targets[i].start, "%c%02d°%02d\'%02d\"");
				gchar *end_ra = siril_world_cs_alpha_format(targets[i].end, " %02dh%02dm%02ds");
				gchar *end_dec = siril_world_cs_delta_format(targets[i].end, "%c%02d°%02d\'%02d\"");
				gchar *start_time = date_time_to_FITS_date(arg->fit->date_obs);
				GDateTime *end_date = g_date_time_add_seconds(arg->fit->date_obs, arg->fit->exposure);
				gchar *end_time = date_time_to_FITS_date(end_date);
				g_date_time_unref(end_date);
				if (targets[i].fixed_target) {
					siril_log_message(_("Target %d (fixed): moved from %s, %s to %s, %s between %s and %s\n"),
							i, start_ra, start_dec, end_ra, end_dec, start_time, end_time);
				} else {
					siril_log_message(_("Target %d (moving, angle %f): moved from %s, %s to %s, %s between %s and %s\n"),
							i, targets[i].moving_track->angle, start_ra, start_dec, end_ra, end_dec, start_time, end_time);
				}
				g_free(start_ra); g_free(start_dec); g_free(end_ra);
				g_free(end_dec); g_free(start_time); g_free(end_time);
			}
		}
		else if ((start_solve_failed && !end_solve_failed) || (!start_solve_failed && end_solve_failed)) {
			for (int i = 0; i < nb_targets; i++) {
				gchar *ra = siril_world_cs_alpha_format(end_solve_failed ? targets[i].start : targets[i].end, " %02dh%02dm%02ds");
				gchar *dec = siril_world_cs_delta_format(end_solve_failed ? targets[i].start: targets[i].end, "%c%02d°%02d\'%02d\"");
				gchar *start_time = date_time_to_FITS_date(arg->fit->date_obs);
				GDateTime *end_date = g_date_time_add_seconds(arg->fit->date_obs, arg->fit->exposure);
				gchar *end_time = date_time_to_FITS_date(end_date);
				if (targets[i].fixed_target) {
					siril_log_message(_("Target %d (fixed): recorded only one position at %s, %s between %s and %s\n"),
							i, ra, dec, start_time, end_time);
				} else {
					siril_log_message(_("Target %d (moving, angle %f): recorded only one position at %s, %s between %s and %s\n"),
							i, targets[i].moving_track->angle, ra, dec, start_time, end_time);
				}
				g_free(ra); g_free(dec); g_free(start_time); g_free(end_time);
			}
		}
		else {
			retval = 1;
		}

		free(tracks);
		free(moving_targets);
	} else {
		siril_log_message(_("%d trail(s) found\n"), nblines);
	}

	free_fitted_stars(arg->fixed_targets);
	if (arg->use_idle)
		siril_add_idle(end_generic, NULL);
	free(arg);
	return GINT_TO_POINTER(retval);
}
