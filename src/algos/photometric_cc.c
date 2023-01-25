/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

#include <gsl/gsl_statistics.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "core/siril_log.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "algos/statistics_float.h"
#include "algos/photometry.h"
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
#include "algos/star_finder.h"
#include "algos/siril_wcs.h"
#include "io/single_image.h"
#include "io/catalogues.h"
#include "gui/progress_and_log.h"
#include "registration/matching/misc.h" // for catalogue parsing helpers
#include "photometric_cc.h"

enum {
	RED, GREEN, BLUE
};

static void bv2rgb(float *r, float *g, float *b, float bv) { // RGB <0,1> <- BV <-0.4,+2.0> [-]
	float t;
	*r = 0.f;
	*g = 0.f;
	*b = 0.f;
	if (bv < -0.4f)
		bv = -0.4f;
	if (bv > 2.f)
		bv = 2.f;
	if ((bv >= -0.4f) && (bv < 0.0f)) {
		t = (bv + 0.4f) / (0.f + 0.4f);
		*r = 0.61f + (0.11f * t) + (0.1f * t * t);
	} else if ((bv >= 0.f) && (bv < 0.4f)) {
		t = (bv - 0.0f) / (0.4f - 0.f);
		*r = 0.83f + (0.17f * t);
	} else if ((bv >= 0.4f) && (bv < 2.1f)) {
		*r = 1.f;
	}
	if ((bv >= -0.4f) && (bv < 0.f)) {
		t = (bv + 0.4f) / (0.f + 0.4f);
		*g = 0.7f + (0.07f * t) + (0.1f * t * t);
	} else if ((bv >= 0.f) && (bv < 0.4f)) {
		t = (bv - 0.f) / (0.4f - 0.f);
		*g = 0.87f + (0.11f * t);
	} else if ((bv >= 0.4f) && (bv < 1.6f)) {
		t = (bv - 0.4f) / (1.6f - 0.4f);
		*g = 0.98f - (0.16f * t);
	} else if ((bv >= 1.6f) && (bv < 2.f)) {
		t = (bv - 1.6f) / (2.f - 1.6f);
		*g = 0.82f - (0.5f * t * t);
	}
	if ((bv >= -0.4f) && (bv < 0.4f)) {
		*b = 1.f;
	} else if ((bv >= 0.4f) && (bv < 1.5f)) {
		t = (bv - 0.4f) / (1.5f - 0.4f);
		*b = 1.f - (0.47f * t) + (0.1f * t * t);
	} else if ((bv >= 1.5f) && (bv < 1.94f)) {
		t = (bv - 1.5f) / (1.94f - 1.5f);
		*b = 0.63f - (0.6f * t * t);
	}
}

static int make_selection_around_a_star(pcc_star star, rectangle *area, fits *fit) {
	/* make a selection around the star */
	double outer = com.pref.phot_set.outer;
	area->x = round_to_int(star.x - outer);
	area->y = round_to_int(star.y - outer);
	area->w = area->h = (int)ceil(outer * 2.0);

	/* Don't want stars too close to the edge */
	if (area->x <= 0 || area->x + area->w >= fit->rx - 1)
		return 1;
	if (area->y <= 0 || area->y + area->h >= fit->ry - 1)
		return 1;

	return 0;
}

static int get_white_balance_coeff(pcc_star *stars, int nb_stars, fits *fit, float *kw, int norm_channel) {
	float *data[3];
	data[RED] = malloc(sizeof(float) * nb_stars);
	data[GREEN] = malloc(sizeof(float) * nb_stars);
	data[BLUE] = malloc(sizeof(float) * nb_stars);
	if (!data[RED] || !data[GREEN] || !data[BLUE]) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	for (int k = 0; k < nb_stars; k++) {
		data[RED][k] = FLT_MAX;
		data[GREEN][k] = FLT_MAX;
		data[BLUE][k] = FLT_MAX;
	}

	gchar *str = ngettext("Applying aperture photometry to %d star.\n", "Applying aperture photometry to %d stars.\n", nb_stars);
	str = g_strdup_printf(str, nb_stars);
	siril_log_message(str);
	g_free(str);

	struct phot_config *ps = phot_set_adjusted_for_image(fit);
	siril_debug_print("aperture: %2.1f%s\tinner: %2.1f\touter: %2.1f\n", ps->aperture, ps->force_radius?"":" (auto)", ps->inner, ps->outer);
	gint ngood = 0, progress = 0;
	gint errors[PSF_ERR_MAX_VALUE] = { 0 };

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) shared(progress, ngood)
#endif
	for (int i = 0; i < nb_stars; i++) {
		if (!get_thread_run())
			continue;
		rectangle area = { 0 };
		float flux[3] = { 0.f, 0.f, 0.f };
		float r, g, b, bv;
		if (!(g_atomic_int_get(&progress) % 16))	// every 16 iterations
			set_progress_bar_data(NULL, (double) progress / (double) nb_stars);
		g_atomic_int_inc(&progress);

		if (make_selection_around_a_star(stars[i], &area, fit)) {
			siril_debug_print("star %d is outside image or too close to border\n", i);
			g_atomic_int_inc(errors+PSF_ERR_OUT_OF_WINDOW);
			continue;
		}

		gboolean no_phot = FALSE;
		psf_error error = PSF_NO_ERR;
		for (int chan = 0; chan < 3 && !no_phot; chan ++) {
			psf_star *photometry = psf_get_minimisation(fit, chan, &area, TRUE, ps, FALSE, com.pref.starfinder_conf.profile, &error);
			g_atomic_int_inc(errors+error);
			if (!photometry || !photometry->phot_is_valid || error != PSF_NO_ERR)
				no_phot = TRUE;
			else flux[chan] = powf(10.f, -0.4f * (float) photometry->mag);
			if (photometry)
				free_psf(photometry);
		}
		if (no_phot) {
			siril_debug_print("photometry failed for star %d, error %d\n", i, error);
			continue;
		}
		/* get r g b coefficient from bv color index */
		bv = stars[i].BV;
		bv2rgb(&r, &g, &b, bv);

		/* get Color calibration factors for current star */
		data[RED][i] = (flux[norm_channel] / flux[RED]) * r;
		data[GREEN][i] = (flux[norm_channel] / flux[GREEN]) * g;
		data[BLUE][i] = (flux[norm_channel] / flux[BLUE]) * b;

		if (xisnanf(data[RED][i]) || xisnanf(data[GREEN][i]) || xisnanf(data[BLUE][i])) {
			data[RED][i] = FLT_MAX;
			data[GREEN][i] = FLT_MAX;
			data[BLUE][i] = FLT_MAX;
			continue;
		}
		g_atomic_int_inc(&ngood);
	}

	free(ps);
	if (!get_thread_run()) {
		free(data[RED]);
		free(data[GREEN]);
		free(data[BLUE]);
		return 1;
	}
	int excl = nb_stars - ngood;
	str = ngettext("%d star excluded from the calculation\n", "%d stars excluded from the calculation\n", excl);
	str = g_strdup_printf(str, excl);
	siril_log_message(str);
	g_free(str);
	if (excl > 0)
		print_psf_error_summary(errors);

	if (ngood == 0) {
		siril_log_message(_("No valid stars found.\n"));
		free(data[RED]);
		free(data[GREEN]);
		free(data[BLUE]);
		return 1;
	}
	/* sort in ascending order before using siril_stats_mean_from_linearFit
	 * Hence, DBL_MAX are at the end of the tab */
	quicksort_f(data[RED], nb_stars);
	quicksort_f(data[GREEN], nb_stars);
	quicksort_f(data[BLUE], nb_stars);

	double deviation[3];
	/* we do not take into account FLT_MAX values */
	kw[RED] = siril_stats_robust_mean(data[RED], 1, ngood, &(deviation[RED]));
	kw[GREEN] = siril_stats_robust_mean(data[GREEN], 1, ngood, &(deviation[GREEN]));
	kw[BLUE] = siril_stats_robust_mean(data[BLUE], 1, ngood, &(deviation[BLUE]));
	if (kw[RED] < 0.f || kw[GREEN] < 0.f || kw[BLUE] < 0.f) {
		free(data[RED]);
		free(data[GREEN]);
		free(data[BLUE]);
		return 1;
	}

	/* normalize factors */
	kw[RED] /= (kw[norm_channel]);
	kw[GREEN] /= (kw[norm_channel]);
	kw[BLUE] /= (kw[norm_channel]);
	siril_log_message(_("Found a solution for color calibration using %d stars. Factors:\n"), ngood);
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3lf\t(deviation: %.3f)\n", chan, kw[chan], deviation[chan]);
	}

	if (ngood < 20)
		siril_log_color_message(_("The photometric color correction has found a solution which may not be perfect because it did not rely on many stars\n"), ngood < 5 ? "red" : "salmon");
	else if (deviation[RED] > 0.1 || deviation[GREEN] > 0.1 || deviation[BLUE] > 0.1)
		siril_log_message(_("The photometric color correction seems to have found an imprecise solution, consider correcting the image gradient first\n"));
	free(data[RED]);
	free(data[GREEN]);
	free(data[BLUE]);
	return 0;
}

static int cmp_coeff(const void *a, const void *b) {
	coeff *a1 = (coeff *) a;
	coeff *a2 = (coeff *) b;
	if (a1->value > a2->value)
		return 1;
	if (a1->value < a2->value)
		return -1;
	return 0;
}

/*
Gets bg, min and max values per channel and sets the chennel with middle bg value
*/
static int get_stats_coefficients(fits *fit, rectangle *area, coeff *bg, float *mins, float *maxs, int *norm_channel) {
	// we cannot use compute_all_channels_statistics_single_image because of the area
	for (int chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, area, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg[chan].value = stat->median;
		bg[chan].channel = chan;
		if (!area) {
			mins[chan] = stat->min;
			maxs[chan] = stat->max;
		}
		free_stats(stat);
	}
	coeff tmp[3];
	memcpy(tmp, bg, 3 * sizeof(coeff));
	/* ascending order */
	qsort(tmp, 3, sizeof(tmp[0]), cmp_coeff);
	//selecting middle channel for norm
	*norm_channel = tmp[1].channel;

	// if no selection for background we have all the stats required
	if (!area) return 0;

	// otherwise, we compute image min/max
	for (int chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, NULL, STATS_MINMAX, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		mins[chan] = stat->min;
		maxs[chan] = stat->max;
	}
	return 0;
}

static int apply_photometric_color_correction(fits *fit, const float *kw, const coeff *bg, const float *mins, const float *maxs, int norm_channel) {
	float maximum = -FLT_MAX;
	float minimum = FLT_MAX;
	float scale[3];
	float offset[3];
	float invrange;

	for (int chan = 0; chan < 3; chan++) {
		maximum = max(maximum, kw[chan] * (maxs[chan] - bg[chan].value) + bg[norm_channel].value);
		minimum = min(minimum, kw[chan] * (mins[chan] - bg[chan].value) + bg[norm_channel].value);
	}
	invrange = ((fit->type == DATA_USHORT) ? USHRT_MAX_SINGLE : 1.f) / (maximum - minimum);

	for (int chan = 0; chan < 3; chan++) {
		scale[chan] = kw[chan] * invrange;
		offset[chan] = (-bg[chan].value * kw[chan] + bg[norm_channel].value  - minimum) * invrange;
	}

	siril_log_message("After renormalization, the following coefficients are applied\n");
	siril_log_color_message(_("White balance factors:\n"), "green");
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3f\n", chan, scale[chan]);
	}
	siril_log_color_message(_("Background reference:\n"), "green");
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("B%d: %+.5e\n", chan, offset[chan]);
	}

	for (int chan = 0; chan < 3; chan++) {
		size_t n = fit->naxes[0] * fit->naxes[1];
		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[chan];
			for (size_t i = 0; i < n; ++i) {
				buf[i] = roundf_to_WORD((float)buf[i] * scale[chan] + offset[chan]);
			}
		}
		else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[chan];
			for (size_t i = 0; i < n; ++i) {
				buf[i] = buf[i] * scale[chan] + offset[chan];
			}
		}
		else return 1;
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/* run the PCC using the existing star list of the image from the provided file */
int photometric_cc(struct photometric_cc_data *args) {
	float kw[3];
	coeff bg[3];
	float mins[3];
	float maxs[3];
	int norm_channel;

	if (!isrgb(args->fit)) {
		siril_log_message(_("Photometric color correction will do nothing for monochrome images\n"));
		return 0;
	}

	for (int i = 0; i < args->nb_stars && i < 40; i++)
		siril_debug_print("star %d: %.2f, %.2f\tmag %.3f, BV %3f\n", i,
				args->stars[i].x, args->stars[i].y,
				args->stars[i].mag, args->stars[i].BV);

	rectangle *bkg_sel = NULL;
	if (!args->bg_auto)
		bkg_sel = &(args->bg_area);

	/* we use the median of each channel to sort them by level and select
	 * the reference channel expressed in terms of order of middle median value */
	if (get_stats_coefficients(args->fit, bkg_sel, bg, mins, maxs, &norm_channel)) {
		siril_log_message(_("failed to compute statistics on image, aborting\n"));
		free(args);
		return 1;
	}
	siril_log_message(_("Normalizing on %s channel.\n"), (norm_channel == 0) ? _("red") : ((norm_channel == 1) ? _("green") : _("blue")));

	/* set photometry parameters to values adapted to the image */
	struct phot_config backup = com.pref.phot_set;
	com.pref.phot_set.force_radius = FALSE;
	com.pref.phot_set.inner = max(7.0, 3.0 * args->fwhm);
	com.pref.phot_set.outer = com.pref.phot_set.inner + 10;
	siril_log_message(_("Photometry radii set to %.1f for inner and %.1f for outer\n"),
			com.pref.phot_set.inner, com.pref.phot_set.outer);

	set_progress_bar_data(_("Photometry color calibration in progress..."), PROGRESS_RESET);
	int ret = get_white_balance_coeff(args->stars, args->nb_stars, args->fit, kw, norm_channel);

	if (!ret) {
		apply_photometric_color_correction(args->fit, kw, bg, mins, maxs, norm_channel);
	} else {
		set_progress_bar_data(_("Photometric Color Calibration failed"), PROGRESS_DONE);
	}

	com.pref.phot_set = backup;
	free(args);
	return ret;
}

float measure_image_FWHM(fits *fit) {
	float fwhm[3];
	/*fits downsampled = { 0 };
	copyfits(fit, &downsampled, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	cvResizeGaussian(&downsampled, DOWNSAMPLE_FACTOR * args->fit->rx, DOWNSAMPLE_FACTOR * args->fit->ry, OPENCV_AREA, FALSE, 0.0);*/
	image im = { .fit = fit, .from_seq = NULL, .index_in_seq = -1 };
	gboolean failed = FALSE;
#ifdef _OPENMP
	int *threads = compute_thread_distribution(3, com.max_thread);
#pragma omp parallel for num_threads(com.max_thread)
#endif
	for (int chan = 0; chan < 3; chan++) {
		int nb_stars;
		int nb_subthreads;
#ifdef _OPENMP
		nb_subthreads = threads[chan];
#else
		nb_subthreads = com.max_thread;
#endif
		psf_star **stars = peaker(&im, chan, &com.pref.starfinder_conf, &nb_stars, NULL, FALSE, TRUE, 200, com.pref.starfinder_conf.profile, nb_subthreads);
		if (stars) {
			fwhm[chan] = filtered_FWHM_average(stars, nb_stars);
			siril_debug_print("FWHM for channel %d: %.3f\n", chan, fwhm[chan]);

			for (int i = 0; i < nb_stars; i++)
				free_psf(stars[i]);
			free(stars);
		}
		else failed = TRUE;
	}
	// clearfits(&downsampled);
	if (failed)
		return 0.0f;
	return max(fwhm[0], max(fwhm[1], fwhm[2]));
}

/* photometric_cc is the entry point for the PCC following the call to the
 * plate solver which gives the star list. We can also run the PCC on a
 * plate-solved image without running plate solving again, this is what this
 * function does.
 */
gpointer photometric_cc_standalone(gpointer p) {
	struct photometric_cc_data *args = (struct photometric_cc_data *)p;
	if (!has_wcs(args->fit)) {
		siril_log_color_message(_("Cannot run the standalone photometric color correction on this image because it has no WCS data or it is not supported\n"), "red");
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	/* run peaker to measure FWHM of the image to adjust photometry settings */
	args->fwhm = measure_image_FWHM(args->fit);
	if (args->fwhm <= 0.0f) {
		siril_log_message(_("Error computing FWHM for photometry settings adjustment\n"));
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	/* get stars from a photometric catalog */
	double ra, dec;
	center2wcs(args->fit, &ra, &dec);
	double resolution = get_wcs_image_resolution(args->fit);
	if ((ra == -1.0 && dec == -1.0) || resolution <= 0.0) {
		siril_log_color_message(_("Cannot run the standalone photometric color correction on this image because it has no WCS data or it is not supported\n"), "red");
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	pcc_star *stars = NULL;
	int nb_stars = 0;
	gboolean image_is_gfit = args->fit == &gfit;

	uint64_t sqr_radius = (gfit.rx * gfit.rx + gfit.ry * gfit.ry) / 4;
	double radius = resolution * sqrt((double)sqr_radius);	// in degrees
	double mag = args->mag_mode == LIMIT_MAG_ABSOLUTE ?
		args->magnitude_arg : compute_mag_limit_from_fov(radius * 2.0);
	if (args->mag_mode == LIMIT_MAG_AUTO_WITH_OFFSET)
		mag += args->magnitude_arg;

	int retval = 0;
	if (args->use_local_cat) {
		siril_log_message(_("Getting stars from local catalogues for PCC, with a radius of %.2f degrees and limit magnitude %.2f\n"), radius * 2.0,  mag);
		if (get_photo_stars_from_local_catalogues(ra, dec, radius, args->fit, mag, &stars, &nb_stars)) {
			siril_log_color_message(_("Failed to get data from the local catalogue, is it installed?\n"), "red");
			retval = 1;
		}
	} else {
		const gchar *cat = NULL;
		switch (args->catalog) {
			case CAT_APASS:
				cat = "APASS";
				mag = min(mag, 17.0);	// in APASS, B is available for V < 17
				break;
			case CAT_NOMAD:
				cat = "NOMAD";
				mag = min(mag, 18.0);	// in NOMAD, B is available for V < 18
				break;
			default:
				siril_log_color_message(_("No valid catalog found.\n"), "red");
				return GINT_TO_POINTER(1);
		}
		siril_log_message(_("Image has a field of view of %.2f degrees, using a limit magnitude of %.2f\n"), radius * 2.0, mag);

		SirilWorldCS *center = siril_world_cs_new_from_a_d(ra, dec);
		GFile *catalog_file = download_catalog(args->catalog, center, radius * 60.0, mag);
		siril_world_cs_unref(center);
		if (!catalog_file) {
			siril_log_message(_("Could not download the online star catalog.\n"));
			siril_add_idle(end_generic, NULL);
			return GINT_TO_POINTER(1);
		}
		siril_log_message(_("The %s catalog has been successfully downloaded.\n"), cat);

		/* project using WCS */
		retval = project_catalog_with_WCS(catalog_file, args->fit, &stars, &nb_stars);
	}

	if (!retval) {
		if (!com.script) {
			undo_save_state(args->fit, _("Photometric CC"));
		}
		args->stars = stars;
		args->nb_stars = nb_stars;
		retval = photometric_cc(args);	// args is freed from here
	}
	free(stars);
	args = NULL;

	if (!retval && image_is_gfit) {
		set_progress_bar_data(_("Photometric Color Calibration succeeded"), PROGRESS_DONE);
		notify_gfit_modified();
	}
	else siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

int project_catalog_with_WCS(GFile *catalog_file, fits *fit, pcc_star **ret_stars, int *ret_nb_stars) {
	GError *error = NULL;
	GInputStream *input_stream = NULL;
	/* catalog format should be 5 columns: distance from centre, RA, Dec, V, B */
	if (!(input_stream = (GInputStream*) g_file_read(catalog_file, NULL, &error))) {
		if (error) {
			siril_log_message(_("Can't open catalog file %s for PCC: %s\n"), g_file_peek_path(catalog_file), error->message);
			g_clear_error(&error);
		}
		*ret_stars = NULL;
		*ret_nb_stars = 0;
		return 1;
	}

	int nb_alloc = 1200, nb_stars = 0;
	pcc_star *stars = malloc(nb_alloc * sizeof(pcc_star));

	/* see also proc_star_file() or read_NOMAD_catalog() */
	// TODO: merge the three codes? they are not quite the same
	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	gchar *line;
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		if (line[0] == COMMENT_CHAR || is_blank(line) || g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		double r = 0.0, ra = 0.0, dec = 0.0, Vmag = 0.0, Bmag = 0.0;
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &ra, &dec, &Vmag, &Bmag);
		g_free(line);
		if (n == 5 && Bmag < 30.0) {	// 30 sometimes means not available in NOMAD
			if (nb_stars >= nb_alloc) {
				nb_alloc *= 2;
				pcc_star *new_array = realloc(stars, nb_alloc * sizeof(pcc_star));
				if (!new_array) {
					PRINT_ALLOC_ERR;
					g_object_unref(data_input);
					free(stars);
					*ret_stars = NULL;
					*ret_nb_stars = 0;
					return 1;
				}
				stars = new_array;
			}

			double x, y;
			if (!wcs2pix(fit, ra, dec, &x, &y)) {
				stars[nb_stars].x = x;
				stars[nb_stars].y = fit->ry - y - 1;
				stars[nb_stars].mag = Vmag;
				stars[nb_stars].BV = Bmag - Vmag;
				nb_stars++;
			}
		}
	}

	g_object_unref(data_input);
	siril_debug_print("projected %d stars from the provided catalogue\n", nb_stars);
	*ret_stars = stars;
	*ret_nb_stars = nb_stars;
	return 0;
}


