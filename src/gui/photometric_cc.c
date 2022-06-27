/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/sleef.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "algos/photometry.h"
#include "algos/PSF.h"
#include "algos/astrometry_solver.h"
#include "algos/star_finder.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/registration_preview.h"
#include "io/single_image.h"
#include "photometric_cc.h"
#include "registration/matching/misc.h" // for catalogue parsing helpers
#include "algos/siril_wcs.h"

enum {
	RED, GREEN, BLUE
};

static rectangle get_bkg_selection();

static void start_photometric_cc() {
	GtkComboBox *norm_box = GTK_COMBO_BOX(lookup_widget("combo_box_cc_norm"));
	GtkToggleButton *auto_bkg = GTK_TOGGLE_BUTTON(lookup_widget("button_cc_bkg_auto"));
	GtkToggleButton *use_wcs_ips_button = GTK_TOGGLE_BUTTON(lookup_widget("use_wcs_ips_button"));
	gboolean plate_solve;

	plate_solve = !gtk_toggle_button_get_active(use_wcs_ips_button);

	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("There is no valid WCS information in the header. Let's make a plate solving.\n"), "salmon");
		plate_solve = TRUE;
	}

	struct astrometry_data *args = NULL;
	struct photometric_cc_data *pcc_args = calloc(1, sizeof(struct photometric_cc_data));
	if (plate_solve) {
		args = calloc(1, sizeof(struct astrometry_data));
		args->fit = &gfit;

		args->for_photometry_cc = TRUE;

		args->pcc = malloc(sizeof(struct photometric_cc_data));
		args->pcc->fit = &gfit;
		args->pcc->bg_auto = gtk_toggle_button_get_active(auto_bkg);
		args->pcc->bg_area = get_bkg_selection();
		args->pcc->n_channel = (normalization_channel) gtk_combo_box_get_active(norm_box);
		args->pcc = pcc_args;
	}

	pcc_args->fit = &gfit;
	pcc_args->bg_auto = gtk_toggle_button_get_active(auto_bkg);
	pcc_args->n_channel = (normalization_channel) gtk_combo_box_get_active(norm_box);
	pcc_args->catalog = get_photometry_catalog();

	set_cursor_waiting(TRUE);

	if (plate_solve) {
		if (!fill_plate_solver_structure_from_GUI(args)) {
			start_in_new_thread(match_catalog, args);
		}
	} else {
		control_window_switch_to_tab(OUTPUT_LOGS);
		start_in_new_thread(photometric_cc_standalone, pcc_args);
	}
}

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

static float Qn0(const float sorted_data[], const size_t stride, const size_t n) {
	const size_t wsize = n * (n - 1) / 2;
	const size_t n_2 = n / 2;
	const size_t k = ((n_2 + 1) * n_2) / 2;
	size_t idx = 0;

	if (n < 2)
		return (0.0);

	float *work = malloc(wsize * sizeof(float));
	if (!work) {
		PRINT_ALLOC_ERR;
		return -1.0f;
	}

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j)
			work[idx++] = fabsf(sorted_data[i] - sorted_data[j]);
	}

	quicksort_f(work, idx);
	float Qn = work[k - 1];

	free(work);
	return Qn;
}

static float siril_stats_robust_mean(const float sorted_data[],
		const size_t stride, const size_t size, double *deviation) {
	float mx = (float) gsl_stats_float_median_from_sorted_data(sorted_data, stride, size);
	float qn0 = Qn0(sorted_data, 1, size);
	if (qn0 < 0)
		return -1.0f;
	float sx = 2.2219f * qn0;
	float *x, mean;
	int i, j;

	x = malloc(size * sizeof(float));
	if (!x) {
		PRINT_ALLOC_ERR;
		return -1.0f;
	}

	for (i = 0, j = 0; i < size; ++i) {
		if (fabsf(sorted_data[i] - mx) <= 3.f * sx) {
			x[j++] = sorted_data[i];
		}
	}
	siril_debug_print("keeping %d samples on %zu for the robust mean (mx: %f, sx: %f)\n", j, size, mx, sx);
	/* not enough stars, try something anyway */
	if (j < 5) {
		mean = siril_stats_trmean_from_sorted_data(0.3f, sorted_data, stride, size);
	} else {
		mean = (float) gsl_stats_float_mean(x, stride, j);
	}
	free(x);

	/* compute the deviation of the mean against the values */
	if (deviation) {
		double dev = 0.0;
		int inliers = 0;
		for (i = 0, j = 0; i < size; ++i) {
			if (fabsf(sorted_data[i] - mx) <= 3.f * sx) {
				dev += fabs(sorted_data[i] - mean);
				inliers++;
			}
		}
		if (inliers)
			*deviation = dev / inliers;
		else *deviation = 1.0;
	}

	return mean;
}

static int get_white_balance_coeff(pcc_star *stars, int nb_stars, fits *fit, float kw[], int n_channel) {
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
		if (!(i % 16))	// every 16 iterations
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
			psf_star *photometry = psf_get_minimisation(fit, chan, &area, TRUE, com.pref.phot_set.force_radius, FALSE, &error);
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
		data[RED][i] = (flux[n_channel] / flux[RED]) * r;
		data[GREEN][i] = (flux[n_channel] / flux[GREEN]) * g;
		data[BLUE][i] = (flux[n_channel] / flux[BLUE]) * b;

		if (xisnanf(data[RED][i]) || xisnanf(data[GREEN][i]) || xisnanf(data[BLUE][i])) {
			data[RED][i] = FLT_MAX;
			data[GREEN][i] = FLT_MAX;
			data[BLUE][i] = FLT_MAX;
			continue;
		}
		g_atomic_int_inc(&ngood);
	}
	if (!get_thread_run())
		return 1;
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
	kw[RED] /= (kw[n_channel]);
	kw[GREEN] /= (kw[n_channel]);
	kw[BLUE] /= (kw[n_channel]);
	siril_log_message(_("Found a solution for color calibration using %d stars. Factors:\n"), ngood);
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3lf\t(deviation: %.3f)\n", chan, kw[chan], deviation[chan]);
	}

	if (ngood < 20)
		siril_log_color_message(_("The photometric color correction has found a solution which may not be perfect because it did not rely on many stars\n"), "salmon");
	else if (deviation[RED] > 0.1 || deviation[GREEN] > 0.1 || deviation[BLUE] > 0.1)
		siril_log_message(_("The photometric color correction seem to have found an imprecise solution, consider correcting the image gradient first\n"));
	free(data[RED]);
	free(data[GREEN]);
	free(data[BLUE]);

	return 0;
}

static int get_background_coefficients(fits *fit, rectangle *area, coeff bg[], gboolean verbose) {
	if (verbose)
		siril_log_message(_("Background reference:\n"));
	// we cannot use compute_all_channels_statistics_single_image because of the area
	for (int chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, area, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		bg[chan].value = stat->median / stat->normValue;
		bg[chan].channel = chan;
		if (verbose)
			siril_log_message("B%d: %.5e\n", chan, bg[chan].value);
		free_stats(stat);
	}
	return 0;
}

static int apply_white_balance(fits *fit, const float kw[]) {
	for (int chan = 0; chan < 3; chan++) {
		float scale = kw[chan];
		if (scale == 1.0) continue;

		size_t i, n = fit->naxes[0] * fit->naxes[1];
		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[chan];
			for (i = 0; i < n; ++i) {
				buf[i] = roundf_to_WORD((float)buf[i] * scale);
			}
		}
		else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[chan];
			for (i = 0; i < n; ++i) {
				buf[i] = buf[i] * scale;
			}
		}
		else return 1;
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/* This function equalize the background by giving equal value for all layers */
static void background_neutralize(fits* fit, coeff bg[], int n_channel, double norm) {
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->type == DATA_USHORT) {
		for (int chan = 0; chan < 3; chan++) {
			float offset = (bg[chan].value - bg[n_channel].value) * norm;
			siril_debug_print("offset: %d, %f\n", chan, offset);
			WORD *buf = fit->pdata[chan];
			for (i = 0; i < n; ++i) {
				buf[i] = roundf_to_WORD((float)buf[i] - offset);
			}
		}
	}
	else if (fit->type == DATA_FLOAT) {
		for (int chan = 0; chan < 3; chan++) {
			float offset = bg[chan].value - bg[n_channel].value;
			siril_debug_print("offset: %d, %f\n", chan, offset);
			float *buf = fit->fpdata[chan];
			for (i = 0; i < n; ++i) {
				buf[i] = buf[i] - offset;
			}
		}
	}
	invalidate_stats_from_fit(fit);
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

static int determine_chan_for_norm(coeff bg[], normalization_channel n_channel) {
	/* make a copy of bg coefficients because we don't
	 * want to sort original data */
	coeff tmp[3];
	memcpy(tmp, bg, 3 * sizeof(coeff));
	/* ascending order */
	qsort(tmp, 3, sizeof(tmp[0]), cmp_coeff);

	if (n_channel == CHANNEL_HIGHEST)
		return tmp[2].channel;
	if (n_channel == CHANNEL_MIDDLE)
		return tmp[1].channel;
	/* on lowest */
	return tmp[0].channel;
}

/* run the PCC using the existing star list of the image from the provided file */
int photometric_cc(struct photometric_cc_data *args) {
	float kw[3];
	coeff bg[3];

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
	 * the reference channel expressed in terms of order of median value */
	if (get_background_coefficients(args->fit, bkg_sel, bg, FALSE)) {
		siril_log_message(_("failed to compute statistics on image, aborting\n"));
		free(args);
		return 1;
	}
	int chan = determine_chan_for_norm(bg, args->n_channel);
	siril_log_message(_("Normalizing on %s channel.\n"), (chan == 0) ? _("red") : ((chan == 1) ? _("green") : _("blue")));

	/* set photometry parameters to values adapted to the image */
	struct phot_config backup = com.pref.phot_set;
	com.pref.phot_set.force_radius = FALSE;
	com.pref.phot_set.inner = max(7.0, 3.0 * args->fwhm);
	com.pref.phot_set.outer = com.pref.phot_set.inner + 10;
	siril_debug_print("set photometry inner radius to %.2f\n", com.pref.phot_set.inner);

	set_progress_bar_data(_("Photometry color calibration in progress..."), PROGRESS_RESET);
	int ret = get_white_balance_coeff(args->stars, args->nb_stars, args->fit, kw, chan);
	if (!ret) {
		double norm = get_normalized_value(args->fit);
		apply_white_balance(args->fit, kw);

		if ((ret = get_background_coefficients(args->fit, bkg_sel, bg, TRUE))) {
			siril_log_message(_("failed to compute statistics on image, aborting\n"));
			set_progress_bar_data(_("Photometric Color Calibration failed"), PROGRESS_DONE);
		} else {
			background_neutralize(args->fit, bg, chan, norm);
			set_progress_bar_data(_("Photometric Color Calibration applied"), PROGRESS_DONE);
		}
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
	cvResizeGaussian(&downsampled, DOWNSAMPLE_FACTOR * args->fit->rx, DOWNSAMPLE_FACTOR * args->fit->ry, OPENCV_AREA);*/
	image im = { .fit = fit, .from_seq = NULL, .index_in_seq = -1 };
	gboolean failed = FALSE;
#ifdef _OPENMP
	int *threads = compute_thread_distribution(3, com.max_thread);
#pragma omp parallel for num_threads(com.max_thread)
#endif
	for (int chan = 0; chan < 3; chan++) {
		int nb_stars;
		int nb_subthreads = com.max_thread;
#ifdef _OPENMP
		nb_subthreads = threads[chan];
#endif
		psf_star **stars = peaker(&im, chan, &com.starfinder_conf, &nb_stars, NULL, FALSE, TRUE, 200, nb_subthreads);
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
	if ((ra == -1.0 && dec == -1.0) || resolution == -1.0) {
		siril_log_color_message(_("Cannot run the standalone photometric color correction on this image because it has no WCS data or it is not supported\n"), "red");
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	/* for now, only from online sources */
	const gchar *cat = NULL;
	SirilWorldCS *center = siril_world_cs_new_from_a_d(ra, dec);
	double fov = resolution  * args->fit->rx;	// fov in degrees
	double mag = compute_mag_limit_from_fov(fov);

	switch(args->catalog) {
	case APASS:
		cat = "APASS";
		mag = min(mag, 17.0);	// in APASS, B is available for V < 17
		break;
	case NOMAD:
		cat = "NOMAD";
		mag = min(mag, 18.0);	// in NOMAD, B is available for V < 18
		break;
	default:
		siril_log_color_message(_("No valid catalog found.\n"), "red");
		return GINT_TO_POINTER(1);
	}
	siril_log_message(_("Image has a field of view of %.2f degrees, using a limit magnitude of %.2f\n"), fov, mag);

	GFile *catalog_file = download_catalog(args->catalog, center, fov * 60.0, mag);
	siril_world_cs_unref(center);
	if (!catalog_file) {
		siril_log_message(_("Could not download the online star catalog."));
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}
	siril_log_message(_("The %s catalog has been successfully downloaded.\n"), cat);

	/* project using WCS */
	pcc_star *stars;
	int nb_stars;
	gboolean image_is_gfit = args->fit == &gfit;
	int retval = project_catalog_with_WCS(catalog_file, args->fit, &stars, &nb_stars);
	if (!retval) {
		args->stars = stars;
		args->nb_stars = nb_stars;
		retval = photometric_cc(args);	// args is freed from here
		free(stars);
		args = NULL;
	}
	if (!retval && image_is_gfit)
		notify_gfit_modified();
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
			wcs2pix(fit, ra, dec, &x, &y);
			if (x >= 0.0 && y >= 0.0 && x < fit->rx && y < fit->ry) {
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

static rectangle get_bkg_selection() {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}

	rectangle black_selection;
	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);
	return black_selection;
}

static gboolean is_selection_ok() {
	rectangle selection = get_bkg_selection();
	return selection.w != 0 && selection.h != 0;
}

/****
 * PUBLIC FUNCTIONS
 */

void initialize_photometric_cc_dialog() {
	GtkWidget *button_ips_ok, *button_cc_ok, *catalog_label, *catalog_box_ips,
			*catalog_box_pcc, *catalog_auto, *frame_cc_bkg, *frame_cc_norm,
			*catalog_label_pcc;
	GtkWindow *parent;
	GtkAdjustment *selection_cc_black_adjustment[4];

	button_ips_ok = lookup_widget("buttonIPS_ok");
	button_cc_ok = lookup_widget("button_cc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_ips = lookup_widget("ComboBoxIPSCatalog");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_auto = lookup_widget("GtkCheckButton_OnlineCat");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	frame_cc_norm = lookup_widget("frame_cc_norm");

	parent = GTK_WINDOW(lookup_widget("ImagePlateSolver_Dial"));

	selection_cc_black_adjustment[0] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_x"));
	selection_cc_black_adjustment[1] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_y"));
	selection_cc_black_adjustment[2] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_w"));
	selection_cc_black_adjustment[3] = GTK_ADJUSTMENT(gtk_builder_get_object(gui.builder, "adjustment_cc_bkg_h"));

	gtk_widget_set_visible(button_ips_ok, FALSE);
	gtk_widget_set_visible(button_cc_ok, TRUE);
	gtk_widget_set_visible(catalog_label, FALSE);
	gtk_widget_set_visible(catalog_label_pcc, TRUE);
	gtk_widget_set_visible(catalog_box_ips, FALSE);
	gtk_widget_set_visible(catalog_box_pcc, TRUE);
	gtk_widget_set_visible(catalog_auto, FALSE);
	gtk_widget_set_visible(frame_cc_bkg, TRUE);
	gtk_widget_set_visible(frame_cc_norm, TRUE);

	gtk_window_set_title(parent, _("Photometric Color Calibration"));

	gtk_adjustment_set_upper(selection_cc_black_adjustment[0], gfit.rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[1], gfit.ry);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[2], gfit.rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[3], gfit.ry);
	gtk_adjustment_set_value(selection_cc_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[3], 0);

}

int get_photometry_catalog() {
	GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("ComboBoxPCCCatalog"));
	if (gtk_combo_box_get_active(box) == 1)
		return APASS;
	return NOMAD;
}

/*****
 * CALLBACKS FUNCTIONS
 */

void on_button_cc_ok_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_bkg;

	auto_bkg = GTK_TOGGLE_BUTTON(lookup_widget("button_cc_bkg_auto"));

	if ((!gtk_toggle_button_get_active(auto_bkg)) && (!is_selection_ok())) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
	}
	else start_photometric_cc();
}

void on_button_cc_bkg_auto_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *box_cc_manual_bkg;

	box_cc_manual_bkg = lookup_widget("box_cc_manual_bkg");
	gtk_widget_set_sensitive(box_cc_manual_bkg, !gtk_toggle_button_get_active(button));
}

void on_button_cc_bkg_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_cc_bkg_value[4] = { NULL, NULL, NULL, NULL };

	if ((!com.selection.h) || (!com.selection.w)) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}

	if (!selection_cc_bkg_value[0]) {
		selection_cc_bkg_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_cc_bkg_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_cc_bkg_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_cc_bkg_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}

	gtk_spin_button_set_value(selection_cc_bkg_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_cc_bkg_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_cc_bkg_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_cc_bkg_value[3], com.selection.h);

	delete_selected_area(); // needed because we don't want the selection being used for astrometry
}
