#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/siril_date.h"
#include "core/proto.h"
#include "opencv/tracks.h"
#include "algos/astrometry_solver.h"
#include "core/siril_world_cs.h"
#include "core/gui_iface.h"
#include "algos/statistics.h"
#include "algos/statistics_float.h"
#include "algos/siril_wcs.h"
#include "algos/streaks.h"
#include "algos/sorting.h"
#include "gui/PSF_list.h"
#include "filters/deconvolution/deconvolution.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "algos/line_detection.h"
//#include "io/gps_parser.h"

// side of the square in which sampling will be done for magnitde estimation
#define MAG_SAMPLE_SIZE 7

// move somewhere else
void ssr_internal(fits *fit, int layer, double median, double bgnoise, double pixvalue, psf_star **stars, int nb_stars) {
	double alpha = bgnoise;// * 0.1;
	size_t pixels_written = 0;
	if (fit->type == DATA_FLOAT) {
		float **rows = malloc(fit->ry * sizeof(float *));
		for (int i = 0; i < fit->ry; i++)
			rows[fit->ry - i - 1] = fit->fpdata[layer] + i * fit->rx;

		for (int i = 0; i < nb_stars; i++) {
			double radius = psf_get_star_radius(stars[i], alpha, median + bgnoise);
			//siril_log_debug("radius = %f\n", radius);
			int ymin = max(0, round_to_int(stars[i]->ypos - radius - 1.0));
			int ymax = min(fit->ry-1, round_to_int(stars[i]->ypos + radius));
			int xmin = max(0, round_to_int(stars[i]->xpos - radius - 1.0));
			int xmax = min(fit->rx-1, round_to_int(stars[i]->xpos + radius));
			for (int y = ymin; y <= ymax; y++) {
				double dy = stars[i]->ypos - y - 0.5;
				for (int x = xmin; x <= xmax; x++) {
					double dx = stars[i]->xpos - x - 0.5;
					double dist = sqrt(dx * dx + dy * dy);
					if (dist > radius)
						continue;
					rows[y][x] = pixvalue;
					pixels_written++;
				}
			}
		}
		free(rows);
	} else {
		WORD pixval = round_to_WORD(pixvalue);
		WORD **rows = malloc(fit->ry * sizeof(WORD *));
		for (int i = 0; i < fit->ry; i++)
			rows[fit->ry - i - 1] = fit->pdata[layer] + i * fit->rx;

		for (int i = 0; i < nb_stars; i++) {
			double radius = psf_get_star_radius(stars[i], alpha, median + bgnoise);
			//siril_log_debug("radius = %f\n", radius);
			int ymin = max(0, round_to_int(stars[i]->ypos - radius - 1.0));
			int ymax = min(fit->ry-1, round_to_int(stars[i]->ypos + radius));
			int xmin = max(0, round_to_int(stars[i]->xpos - radius - 1.0));
			int xmax = min(fit->rx-1, round_to_int(stars[i]->xpos + radius));
			for (int y = ymin; y <= ymax; y++) {
				double dy = stars[i]->ypos - y - 0.5;
				for (int x = xmin; x <= xmax; x++) {
					double dx = stars[i]->xpos - x - 0.5;
					double dist = sqrt(dx * dx + dy * dy);
					if (dist > radius)
						continue;
					rows[y][x] = pixval;
					pixels_written++;
				}
			}
		}
		free(rows);
	}
	siril_log_message("Erased %d stars (%zd pixels written).\n", nb_stars, pixels_written);
}

// size is the number of frames in the sequence. Indices between 0 and size-1 are assumed.
struct result_set *alloc_results(int size) {
	struct result_set *set = malloc(sizeof(struct result_set));
	if (!set) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	set->data = calloc(size, sizeof(GSList *));
	if (!set->data) {
		free(set);
		PRINT_ALLOC_ERR;
		return NULL;
	}
	set->size = size;
	set->has_data = FALSE;
	return set;
}

void free_results(struct result_set *set) {
	for (int i = 0; i < set->size; i++) {
		g_slist_free_full(set->data[i], free);
	}
	free(set->data);
	free(set);
}

/* do we want to store the ends of the streaks or the middle?
 * in most cases the middle is more accurate and does not need to know the direction of travel.
 * But when the dataset is very small, it can be useful to also have the ends.
 */
static void add_result_middle_and_ends(struct result_set *set, int image_idx, int filenum,
		int target_idx, float angle, double start_ra, double start_dec, GDateTime *start_date,
		double middle_ra, double middle_dec, GDateTime *middle_date, double end_ra, double
		end_dec, GDateTime *end_date, double image_center_ra, double image_center_dec,
		int streak_1x, int streak_1y, int streak_2x, int streak_2y, float mag,
		float mag_err, gboolean mag_reliable, gboolean mag_absolute, float snr, gboolean
		middle_date_is_gps) {
	if (set->size != 1 && image_idx >= set->size) {
		siril_log_error(_("Cannot store result, bad index %d, size is %d\n"),
				image_idx, set->size);
		return;
	}
	struct streak_result *cur = malloc(sizeof(struct streak_result));
	cur->filenum = filenum;
	cur->target_idx = target_idx;
	cur->angle = angle;
	cur->start_ra = start_ra;
	cur->start_dec = start_dec;
	cur->start_date = start_date;
	cur->middle_ra = middle_ra;
	cur->middle_dec = middle_dec;
	cur->middle_date = middle_date;
	cur->middle_date_is_gps = middle_date_is_gps;
	cur->end_ra = end_ra;
	cur->end_dec = end_dec;
	cur->end_date = end_date;
	cur->image_center_ra = image_center_ra;
	cur->image_center_dec = image_center_dec;
	cur->streak_start.x = streak_1x;
	cur->streak_start.y = streak_1y;
	cur->streak_end.x = streak_2x;
	cur->streak_end.y = streak_2y;
	cur->magnitude = mag;
	cur->magnitude_error = mag_err;
	cur->magnitude_reliable = mag_reliable;
	cur->magnitude_absolute = mag_absolute;
	cur->snr = snr;
	if (set->size == 1)
		set->data[0] = g_slist_prepend(set->data[0], cur);
	else set->data[image_idx] = g_slist_prepend(set->data[image_idx], cur);
	set->has_data = TRUE;
}

/*static void add_result_middle(struct result_set *set, int image_idx, int filenum,
		int target_idx, float angle, double middle_ra,
		double middle_dec, GDateTime *middle_date,
		double image_center_ra, double image_center_dec,
		int streak_1x, int streak_1y, int streak_2x, int streak_2y,
		gboolean date_is_gps) {
	add_result_middle_and_ends(set, image_idx, filenum, target_idx, angle,
			0.0, 0.0, NULL, middle_ra, middle_dec, middle_date, 0.0, 0.0, NULL,
			image_center_ra, image_center_dec,
			streak_1x, streak_1y, streak_2x, streak_2y,
			99.9f, 0.0f, FALSE, FALSE, 0.0f, date_is_gps);
}*/

// also frees the result_set
void dump_results(struct result_set *set, const char *filename) {
	FILE *fd = fopen(filename, "w");
	if (!fd) {
		siril_log_error(_("Could not save the result CSV file %s\n"), filename);
		fd = stdout;
	}
	else siril_log_message("Saving results to %s\n", filename);

	fprintf(fd, "# T,image,target,angle (deg),start_ra (deg),start_dec (deg),start_date,mid_ra (deg),mid_dec (deg),mid_date,end_ra (deg),end_dec (deg),end_date,image center ra (deg),image center dec (deg),streak start x (px),streak start y (px),streak end x (px),streak end y (px),magnitude,mag_error,mag_is_reliable,mag_is_absolute,snr,mid_time is GPS\n");
	for (int i = 0; i < set->size; i++) {
		set->data[i] = g_slist_reverse(set->data[i]);
		GSList *cur = set->data[i];
		while (cur) {
			struct streak_result *res = cur->data;
			if (res->start_date && res->middle_date && res->end_date) {
				gchar *date1 = date_time_to_FITS_date(res->start_date);
				gchar *date2 = date_time_to_FITS_date(res->middle_date);
				gchar *date3 = date_time_to_FITS_date(res->end_date);
				fprintf(fd, "T,%d,%d,%f,%f,%f,%sZ,%f,%f,%sZ,%f,%f,%sZ,%f,%f,%d,%d,%d,%d,%f,%f,%d,%d,%f,%d\n",
						res->filenum, res->target_idx, res->angle, res->start_ra,
						res->start_dec, date1, res->middle_ra, res->middle_dec,
						date2, res->end_ra, res->end_dec, date3,
						res->image_center_ra, res->image_center_dec,
						res->streak_start.x, res->streak_start.y,
						res->streak_end.x, res->streak_end.y,
						res->magnitude, res->magnitude_error, res->magnitude_reliable,
						res->magnitude_absolute, res->snr, res->middle_date_is_gps);
				g_free(date1);
				g_free(date2);
				g_free(date3);
				g_date_time_unref(res->start_date);
				g_date_time_unref(res->middle_date);
				g_date_time_unref(res->end_date);
			}
			else if (res->start_date && res->end_date) {
				gchar *date1 = date_time_to_FITS_date(res->start_date);
				gchar *date3 = date_time_to_FITS_date(res->end_date);
				fprintf(fd, "T,%d,%d,%f,%f,%f,%sZ,0,0,0,%f,%f,%sZ,%f,%f,%d,%d,%d,%d,%f,%f,%d,%d,%f,0\n",
						res->filenum, res->target_idx, res->angle, res->start_ra,
						res->start_dec, date1, res->end_ra, res->end_dec, date3,
						res->image_center_ra, res->image_center_dec,
						res->streak_start.x, res->streak_start.y,
						res->streak_end.x, res->streak_end.y,
						res->magnitude, res->magnitude_error, res->magnitude_reliable,
						res->magnitude_absolute, res->snr);
				g_free(date1);
				g_free(date3);
				g_date_time_unref(res->start_date);
				g_date_time_unref(res->end_date);

			} else {	// only middle
				gchar *date = date_time_to_FITS_date(res->middle_date);
				fprintf(fd, "T,%d,%d,%f,0,0,0,%f,%f,%sZ,0,0,0,%f,%f,%d,%d,%d,%d,%f,%f,%d,%d,%f,%d\n",
						res->filenum, res->target_idx, res->angle, res->middle_ra,
						res->middle_dec, date, res->image_center_ra,
						res->image_center_dec, res->streak_start.x,
						res->streak_start.y, res->streak_end.x, res->streak_end.y,
						res->magnitude, res->magnitude_error, res->magnitude_reliable,
						res->magnitude_absolute, res->snr, res->middle_date_is_gps);
				g_free(date);
				g_date_time_unref(res->middle_date);
			}
			cur = cur->next;
		}
		g_slist_free_full(set->data[i], free);
	}

	free(set->data);
	free(set);
	fclose(fd);
}

gpointer streak_detection_worker(gpointer ptr) {
	struct streak_detection_conf *arg = (struct streak_detection_conf *)ptr;
	struct results *lines;
	int initial_length = arg->initial_segment_length > 0 ? arg->initial_segment_length : 60;
	struct streak_detector conf = {
		.min_allowed_length = arg->minimum_segment_length,
		.min_allowed_ksigma = -1,
		.initial_length = initial_length,
		.bright_streak = arg->bright_target,
		.max_allowed_segments = 5,
		.enough_segments = 1,
		.can_recurse = TRUE,
		.compute_flux = TRUE,
		.fwhm = arg->fwhm,
	};

	lines = detect_streaks(arg->fit, arg->layer, &conf, arg->nb_threads);
	if (!lines || lines->nb_tracks == 0) {
		siril_log_message("No satellite detected for image %d\n", arg->im_idx);
		if (arg->use_idle)
			siril_add_idle(end_generic, NULL);
		if (arg->free_fit) {
			clearfits(arg->fit);
			free(arg->fit);
		}
		free(arg);
		return GINT_TO_POINTER(1);
	}
	siril_log_message(_("Found %d trail(s) in frame %2d\n"), lines->nb_tracks, arg->im_idx);

	int nblines = lines->nb_tracks;
	if (arg->display_streaks) {
		psf_star **stars = malloc((2 * nblines + 1) * sizeof(psf_star *));
		if (stars) {
			for (int i = 0; i < nblines; i++) {
				stars[2*i] = new_psf_star();
				stars[2*i]->xpos = lines->tracks[i].start.x;
				stars[2*i]->ypos = lines->tracks[i].start.y;
				stars[2*i]->fwhmx = 5.0;
				stars[2*i]->fwhmy = 5.0;
				stars[2*i]->has_saturated = FALSE;
				// we could also add an angle and uneven fwhms
				stars[2*i+1] = new_psf_star();
				stars[2*i+1]->xpos = lines->tracks[i].end.x;
				stars[2*i+1]->ypos = lines->tracks[i].end.y;
				stars[2*i+1]->fwhmx = 5.0;
				stars[2*i+1]->fwhmy = 5.0;
				stars[2*i+1]->has_saturated = TRUE;
			}
			stars[2*nblines] = NULL;
		}
		gui_iface.update_star_list(stars, FALSE, FALSE);
	}

	//TODO: check before calling the function that has_wcs(arg->fit)
	double image_center_ra = arg->fit->keywords.wcsdata.ra, image_center_dec = arg->fit->keywords.wcsdata.dec;
	for (int i = 0; i < nblines; i++) {
		double start_ra, start_dec, end_ra, end_dec, middle_ra, middle_dec;
		double fx, fy;
		display_to_siril(lines->tracks[i].start.x, lines->tracks[i].start.y, &fx, &fy, arg->fit->ry);
		pix2wcs(arg->fit, fx, fy, &start_ra, &start_dec);

		display_to_siril(lines->tracks[i].end.x, lines->tracks[i].end.y, &fx, &fy, arg->fit->ry);
		pix2wcs(arg->fit, fx, fy, &end_ra, &end_dec);

		double middle_x = (lines->tracks[i].start.x + lines->tracks[i].end.x) * 0.5;
		double middle_y = (lines->tracks[i].start.y + lines->tracks[i].end.y) * 0.5;
		display_to_siril(middle_x, middle_y, &fx, &fy, arg->fit->ry);
		pix2wcs(arg->fit, fx, fy, &middle_ra, &middle_dec);
		if ((start_ra == 0.0 && start_dec == 0.0) || (middle_ra == 0.0 && middle_dec == 0.0) ||
				(end_ra == 0.0 && end_dec == 0.0))
			continue;

		GDateTime *center_date = NULL;
		gboolean date_is_gps = FALSE;
		/*if (arg->fit->date_and_exp_from_gps && !arg->fit->gps_data) {
			date_is_gps = TRUE;
		} else {
			// check for rolling shutter QHY GPS camera data
			int mid_y = round_to_int(middle_y);
			siril_log_debug("getting GPS timestamp for middle exposure of row %d\n", mid_y);
			center_date = get_timestamp_for_pixel(arg->fit->gps_data, EXP_MIDDLE, 0, mid_y);
			date_is_gps = center_date != NULL;
		}*/
		/* fallback: use DATE-OBS + EXPOSURE / 2 */
		double exposure;
		if (arg->fit->keywords.expstart > 0.0 && arg->fit->keywords.expend > 0.0) {
			// warning: this is precise to 0.1 seconds at best
			exposure = julian_date_difference_sec(arg->fit->keywords.expend, arg->fit->keywords.expstart);
			siril_log_debug("using EXPSTART and EXPEND %f as exposure\n", exposure);
		}
		else exposure = arg->fit->keywords.exposure;
		if (!center_date)
			center_date = g_date_time_add_seconds(arg->fit->keywords.date_obs,
					exposure * 0.5);

		GDateTime *start_date = g_date_time_ref(arg->fit->keywords.date_obs);
		GDateTime *end_date = g_date_time_add_seconds(arg->fit->keywords.date_obs, exposure);

		add_result_middle_and_ends(arg->results, arg->im_idx, arg->filenum, i,
				lines->tracks[i].angle, start_ra, start_dec, start_date,
				middle_ra, middle_dec, center_date, end_ra, end_dec, end_date,
				image_center_ra, image_center_dec,
				lines->tracks[i].start.x, lines->tracks[i].start.y,
				lines->tracks[i].end.x, lines->tracks[i].end.y,
				lines->tracks[i].mag, lines->tracks[i].mag_err,
				lines->tracks[i].mag_is_accurate, lines->tracks[i].mag_is_absolute,
				lines->tracks[i].snr, date_is_gps);
	} // end of streak loop

	clear_results(lines);
	if (arg->use_idle) {
		/* from using the single-file command in the GUI */
		if (arg->results->has_data) {
			gchar *filename = replace_ext(arg->filename, ".streaks");
			dump_results(arg->results, filename);
			g_free(filename);
		}
		siril_add_idle(end_generic, NULL);
	}
	if (arg->free_fit) {
		clearfits(arg->fit);
		free(arg->fit);
	}
	free(arg);
	return GINT_TO_POINTER(0);
}

void display_tracks(struct track *tracks, int nblines) {
	psf_star **stars = malloc((2 * nblines + 1) * sizeof(psf_star *));
	if (stars) {
		for (int i = 0; i < nblines; i++) {
			stars[2*i] = new_psf_star();
			stars[2*i]->xpos = tracks[i].start.x;
			stars[2*i]->ypos = tracks[i].start.y;
			stars[2*i]->fwhmx = 5.0;
			stars[2*i]->fwhmy = 5.0;
			stars[2*i]->has_saturated = FALSE;
			// we could also add an angle and uneven fwhms
			stars[2*i+1] = new_psf_star();
			stars[2*i+1]->xpos = tracks[i].end.x;
			stars[2*i+1]->ypos = tracks[i].end.y;
			stars[2*i+1]->fwhmx = 5.0;
			stars[2*i+1]->fwhmy = 5.0;
			stars[2*i+1]->has_saturated = TRUE; // color change
		}
		stars[2*nblines] = NULL;
	}
	gui_iface.update_star_list(stars, FALSE, FALSE);
}

gboolean has_streaks(fits *fit, int layer, int nb_threads) {
	struct streak_detector conf = {
		.min_allowed_length = -1,
		.min_allowed_ksigma = -1,
		.initial_length = 500,
		.max_allowed_segments = 20,
		.enough_segments = 1,
		.can_recurse = TRUE
	};

	struct results *lines = detect_streaks(fit, layer, &conf, nb_threads);
	gboolean retval = lines != NULL && lines->nb_tracks > 0;
	if (!com.script && retval)
		display_tracks(lines->tracks, lines->nb_tracks);
	clear_results(lines);
	return retval;
}
