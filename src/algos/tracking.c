#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
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
	psf_star *target;
	SirilWorldCS *start, *end;
	double max_motion;
};

psf_star **extract_stars(struct track *tracks, int nbtracks, gboolean start, psf_star **stars) {
	/* we have a list of tracks, some may not be stars, only keep the tracks that have the
	 * majority of similar orientation */
	float *angles = malloc(nbtracks * sizeof(float));
	for (int i = 0; i < nbtracks; i++)
		angles[i] = tracks[i].angle;
	quicksort_f(angles, nbtracks);
	float mean = siril_stats_robust_mean(angles, 1, nbtracks, NULL);
	free(angles);
	siril_log_message(_("considering trails with angle close to %f for stars\n"), mean);

	int j = 0;
	if (!stars && start) {
		// assuming stars is NULL
		stars = malloc((nbtracks + 1) * sizeof(psf_star *));
		for (int i = 0; i < nbtracks; i++) {
			if (fabsf(mean - tracks[i].angle) > 10.0)
				continue;
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

	struct tracked_object *targets = calloc(arg->nb_targets, sizeof(struct tracked_object));
	for (int i = 0; i < arg->nb_targets; i++)
		targets[i].target = arg->targets[i];

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

		psf_star **stars = extract_stars(tracks, nblines, TRUE, NULL);
		if (plate_solve_for_stars(arg->fit, stars)) {
			siril_log_message(_("Plate solving using the detected trails as stars failed, adjust parameters or make sure the file metadata is correct\n"));
			free(tracks);
			free_fitted_stars(arg->targets);
			free_fitted_stars(stars);
			if (arg->use_idle)
				siril_add_idle(end_generic, NULL);
			return GINT_TO_POINTER(1);
		}
		for (int i = 0; i < arg->nb_targets; i++) {
			double ra, dec;
			pix2wcs(arg->fit, targets[i].target->xpos, targets[i].target->ypos, &ra, &dec);
			targets[i].start = siril_world_cs_new_from_a_d(ra, dec);
		}
		siril_log_message(_("Plate solving for start of movement succeeded.\n"));

		stars = extract_stars(tracks, nblines, FALSE, stars);
		if (plate_solve_for_stars(arg->fit, stars)) {
			siril_log_message(_("Plate solving using the detected trails as stars failed, adjust parameters or make sure the file metadata is correct\n"));
			free(tracks);
			free_fitted_stars(arg->targets);
			free_fitted_stars(stars);
			if (arg->use_idle)
				siril_add_idle(end_generic, NULL);
			return GINT_TO_POINTER(1);
		}
		free_fitted_stars(stars);
		for (int i = 0; i < arg->nb_targets; i++) {
			double ra, dec;
			pix2wcs(arg->fit, targets[i].target->xpos, targets[i].target->ypos, &ra, &dec);
			targets[i].end = siril_world_cs_new_from_a_d(ra, dec);
		}

		/* TODO: estimate the width of tracks, which would give some error margin for
		 * the position of the targets, in addition to their FWHM and elongation */
		for (int i = 0; i < arg->nb_targets; i++) {
			gchar *start_ra = siril_world_cs_alpha_format(targets[i].start, " %02dh%02dm%02ds");
			gchar *start_dec = siril_world_cs_delta_format(targets[i].start, "%c%02d°%02d\'%02d\"");
			gchar *end_ra = siril_world_cs_alpha_format(targets[i].end, " %02dh%02dm%02ds");
			gchar *end_dec = siril_world_cs_delta_format(targets[i].end, "%c%02d°%02d\'%02d\"");
			siril_log_message(_("Target %d: moved from %s, %s to %s, %s\n"), i,
					start_ra, start_dec, end_ra, end_dec);
			g_free(start_ra); g_free(start_dec); g_free(end_ra); g_free(end_dec);
		}
		free(tracks);
	} else {
		siril_log_message(_("%d trail(s) found\n"), nblines);
	}

	free_fitted_stars(arg->targets);
	if (arg->use_idle)
		siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(0);
}

