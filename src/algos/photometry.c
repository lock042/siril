/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_statistics_double.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "algos/sorting.h"
#include "algos/PSF.h"
#include "algos/photometry.h"
#include "algos/astrometry_solver.h"
#include "algos/statistics_float.h"
#include "algos/siril_wcs.h"
#include "algos/search_objects.h"
#include "algos/comparison_stars.h"
#include "gui/PSF_list.h"
#include "gui/plot.h"
#include "gui/image_display.h"
#include "gui/siril_plot.h"
#include "gui/occultation.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/siril_plot.h"
#include "opencv/opencv.h"

#define MIN_SKY    5	// min number of backgroun pixels for valid photometry

static double getMagnitude(double intensity) {
	return -2.5 * log10(intensity);
}

static double getMagErr(double intensity, double area, int nsky, double skysig, double cvf, double *SNR) {
	double skyvar = skysig * skysig;/* variance of the sky brightness */
	double sigsq = skyvar / nsky;	/* square of the standard error of the mean sky brightness */
	double err1 = area * skyvar;
	double err2 = intensity / cvf;
	double err3 = sigsq * area * area;
	double noise = sqrt(err1 + err2 + err3);

	*SNR = 10.0 * log10(intensity / noise);

	return fmin(9.999, 1.0857 * noise / intensity);
}

struct phot_config *phot_set_adjusted_for_image(const fits *fit) {
	struct phot_config *retval = malloc(sizeof(struct phot_config));
	if (!retval) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	memcpy(retval, &com.pref.phot_set, sizeof(struct phot_config));
	if (fit->cvf > 0.0) {
		// always overwrite with value from the image if available
		retval->gain = fit->cvf;
	}
	if (fit->type == DATA_FLOAT) {
		retval->gain *= USHRT_MAX_DOUBLE;
		retval->minval /= USHRT_MAX_DOUBLE;
		retval->maxval /= USHRT_MAX_DOUBLE;
	}
	//siril_debug_print("phot_set min=%f, max=%f\n", retval->minval, retval->maxval);
	return retval;
}

/* Function that compute all photometric data. The result must be freed */
photometry *getPhotometryData(gsl_matrix* z, const psf_star *psf,
		struct phot_config *phot_set, gboolean verbose, psf_error *error) {
	int width = z->size2;
	int height = z->size1;
	int n_sky = 0, ret;
	int x, y, x1, y1, x2, y2;
	double r1, r2, r, rmin_sq, appRadius;
	double apmag = 0.0, mean = 0.0, stdev = 0.0, area = 0.0;
	gboolean valid = TRUE;
	photometry *phot;

	if (!phot_set) {
		fprintf(stderr, "invalid pointer for phot_set\n");
		if (error) *error = PSF_ERR_ALLOC;
		return NULL;
	}

	double xc = psf->x0;
	double yc = psf->y0;

	if (xc <= 0.0 || yc <= 0.0 || xc >= width || yc >= height) {
		if (error) *error = PSF_ERR_OUT_OF_WINDOW;
		return NULL;
	}

	r1 = phot_set->inner;
	r2 = phot_set->outer;
	appRadius = phot_set->force_radius ? phot_set->aperture : psf->fwhmx * 2.0;
	if (appRadius >= r1 && !phot_set->force_radius) {
		if (verbose) {
			/* Translator note: radii is plural for radius */
			siril_log_message(_("Inner and outer radii are too small (%d required for inner). Please update values in preferences or with setphot.\n"), round_to_int(appRadius));
		}
		if (error) *error = PSF_ERR_INNER_TOO_SMALL;
		return NULL;
	}

	/* compute the bounding box of the outer radius around the star */
	x1 = xc - r2;
	if (x1 < 1)
		x1 = 1;
	x2 = xc + r2;
	if (x2 > width - 1)
		x2 = width - 1;
	y1 = yc - r2;
	if (y1 < 1)
		y1 = 1;
	y2 = yc + r2;
	if (y2 > height - 1)
		y2 = height - 1;

	int ndata = (y2 - y1) * (x2 - x1);
	if (ndata <= 0) {
		siril_log_color_message(_("An error occurred in your selection. Please make another selection.\n"), "red");
		if (error) *error = PSF_ERR_OUT_OF_WINDOW;
		return NULL;
	}
	double *data = calloc(ndata, sizeof(double));
	if (!data) {
		PRINT_ALLOC_ERR;
		if (error) *error = PSF_ERR_ALLOC;
		return NULL;
	}

	r1 *= r1;	// we square the radii to avoid doing
	r2 *= r2;	// sqrts for radius checks in the loop
	rmin_sq = (appRadius - 0.5) * (appRadius - 0.5);
	double lo = phot_set->minval, hi = phot_set->maxval;

	/* from the matrix containing pixel data, we extract pixels within
	 * limits of pixel value and of distance to the star centre for
	 * background evaluation */
	for (y = y1; y <= y2; ++y) {
		int yp = (y - yc) * (y - yc);
		for (x = x1; x <= x2; ++x) {
			r = yp + (x - xc) * (x - xc);
			double pixel = gsl_matrix_get(z, y, x);
			if (pixel > lo && pixel < hi) {
				double f = (r < rmin_sq ? 1 : appRadius - sqrt(r) + 0.5);
				if (f >= 0) {
					area += f;
					apmag += pixel * f;
				}
				/* annulus */
				if (r < r2 && r > r1) {
					data[n_sky] = pixel;
					n_sky++;
				}
			} else {
				valid = FALSE;
				if (error) *error = PSF_ERR_INVALID_PIX_VALUE;
			}
		}
	}
	if (area < 1.0) {
		siril_debug_print("area is < 1: not enough pixels of star data, too small aperture?\n");
		free(data);
		if (error) *error = PSF_ERR_APERTURE_TOO_SMALL;
		return NULL;
	}
	if (n_sky < MIN_SKY) {
		if (verbose)
			siril_log_message(_("Warning: There aren't enough pixels"
						" in the sky annulus. You need to make a larger selection.\n"));
		if (error) *error = PSF_ERR_TOO_FEW_BG_PIX;
		free(data);
		return NULL;
	}

	ret = robustmean(n_sky, data, &mean, &stdev);
	free(data);
	if (ret > 0) {
		if (error) *error = PSF_ERR_MEAN_FAILED;
		return NULL;
	}

	phot = calloc(1, sizeof(photometry));
	if (!phot) {
		if (error) *error = PSF_ERR_ALLOC;
		PRINT_ALLOC_ERR;
	}
	else {
		double SNR = 0.0;
		double signalIntensity = apmag - (area * mean);

		phot->mag = getMagnitude(signalIntensity);
		phot->s_mag = getMagErr(signalIntensity, area, n_sky, stdev, phot_set->gain, &SNR);
		if (phot->s_mag < 9.999) {
			phot->SNR = SNR;
			if (valid && error) *error = PSF_NO_ERR;
		} else {
			phot->SNR = 0.0;
			valid = FALSE;
			if (error && *error == PSF_NO_ERR)
				*error = PSF_ERR_INVALID_STD_ERROR;
		}
		phot->valid = valid;
	}

	return phot;
}

/* used only by unit tests */
void initialize_photometric_param() {
	com.pref.phot_set.inner = 20;
	com.pref.phot_set.outer = 30;
	com.pref.phot_set.aperture = 10;
	com.pref.phot_set.force_radius = FALSE;
	com.pref.phot_set.auto_inner_factor = 4.2;
	com.pref.phot_set.auto_outer_factor = 6.3;
	com.pref.phot_set.gain = 2.3;
	com.pref.phot_set.minval = -1000;
	com.pref.phot_set.maxval = 60000;
}

const char *psf_error_to_string(psf_error err) {
	switch (err) {
		case PSF_NO_ERR:
			return _("no error");
		case PSF_ERR_ALLOC:
			return _("memory allocation");
		case PSF_ERR_UNSUPPORTED:
			return _("unsupported image type");
		case PSF_ERR_DIVERGED:
			return _("PSF fit failed");
		case PSF_ERR_OUT_OF_WINDOW:
			return _("not in area");
		case PSF_ERR_INNER_TOO_SMALL:
			return _("inner radius too small");
		case PSF_ERR_APERTURE_TOO_SMALL:
			return _("aperture too small");
		case PSF_ERR_TOO_FEW_BG_PIX:
			return _("not enough background");
		case PSF_ERR_MEAN_FAILED:
			return _("statistics failed");
		case PSF_ERR_INVALID_STD_ERROR:
			return _("invalid measurement error");
		case PSF_ERR_INVALID_PIX_VALUE:
			return _("pixel out of range");
		case PSF_ERR_WINDOW_TOO_SMALL:
			return _("area too small");
		case PSF_ERR_INVALID_IMAGE:
			return _("image is invalid");
		case PSF_ERR_FLUX_RATIO:
			return _("flux ratio failed");
		default:
			return _("unknown error");
	}
}

void print_psf_error_summary(gint *code_sums) {
	GString *msg = g_string_new("Distribution of errors: ");
	gboolean first = TRUE;
	for (int i = 0; i < PSF_ERR_MAX_VALUE; i++) {
		if (code_sums[i] > 0) {
			if (!first)
				msg = g_string_append(msg, ", ");
			g_string_append_printf(msg, "%d %s", code_sums[i], psf_error_to_string(i));
			first = FALSE;
		}
	}

	gchar *str = g_string_free(msg, FALSE);
	siril_log_message("%s\n", str);
	g_free(str);
}

// save the light curve to a dat file
// with formatting compatible with ETD
static gboolean siril_plot_save_ETD_light_curve(siril_plot_data *spl_data, const char *datfilename, gboolean add_title) {
	GString *header = NULL;
	FILE* fileout = NULL;
	gboolean retval = TRUE;
	double *data = NULL;
	if (g_list_length(spl_data->plots) > 1) {
		siril_debug_print("Light curve should not hold more than one data series, aborting\n");
		return FALSE;
	}

	int nbpoints = 0, nbcols = 3;
	if (add_title && spl_data->title) {
		// spl_data->title is assumed to have the # signs at each line start as necessary
		// and to finish by a \n character
		header = g_string_new(spl_data->title);
		g_string_append_printf(header, "# JD_UT V-C err");
	} else
		header = g_string_new("# JD_UT V-C err");

	// xy points with y error bars
	splxyerrdata *plots = (splxyerrdata *)spl_data->plots->data;
	nbpoints = plots->nb;

	// gathering all the data
	data = malloc(nbpoints * nbcols * sizeof(double));
	if (!data) {
		PRINT_ALLOC_ERR;
		retval = FALSE;
		goto clean_and_exit;
	}

	// writing JD
	int index = 0;
	for (int i = 0; i < nbpoints; i++) {
		data[index] = plots->plots[0]->data[i].x + plots->plots[0]->x_offset; // adding JD
		index += nbcols;
	}

	// writing V-C and error
	index = 1;
	for (int i = 0; i < nbpoints; i++) {
		for (int k = 0; k < 2; k++)
			data[index + k] = plots->plots[k]->data[i].y;
		index += nbcols;
	}

	fileout = g_fopen(datfilename, "w");
	if (fileout == NULL) {
		siril_log_message(_("Could not create %s, aborting\n"));
		retval = FALSE;
		goto clean_and_exit;
	}
	fprintf(fileout, "%s", header->str);
	index = 0;
	for (int r = 0 ; r < nbpoints ; r++) {
		fprintf(fileout, "\n%f", data[index++]); // print newline and x
		for (int c = 1 ; c < nbcols ; c++)
			fprintf(fileout, " %g", data[index++]);
	}
	fclose(fileout);
	siril_log_message(_("%s has been saved.\n"), datfilename);

clean_and_exit:
	g_string_free(header, TRUE);
	free(data);
	return retval;
}

/****************** making a light curve from sequence-stored data ****************/
/* It uses data stored in the sequence, in seq->photometry, which is
 * populated by successive calls to seqpsf on the opened sequence;
 */
int new_light_curve(const char *filename, struct light_curve_args *lcargs) {
	int i, j;
	siril_plot_data *spl_data = NULL;
	sequence *seq = lcargs->seq;

	if (!seq->photometry[0]) {
		siril_log_color_message(_("No photometry data found, error\n"), "red");
		return -1;
	}

	/* get number of valid frames for each star */
	int ref_valid_count[MAX_SEQPSF] = { 0 };
	int nbImages = 0;
	gboolean ref_valid[MAX_SEQPSF] = { FALSE };
	for (i = 0; i < seq->number; i++) {
		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid)
			continue;
		++nbImages;
		for (int ref = 1; ref < MAX_SEQPSF && seq->photometry[ref]; ref++) {
			if (seq->photometry[ref][i] && seq->photometry[ref][i]->phot_is_valid)
				ref_valid_count[ref]++;
		}
	}
	siril_debug_print("we have %d images with a valid photometry for the variable star\n", nbImages);
	if (nbImages < 1) {
		siril_log_color_message(_("There are not enough valid stars to make a photometric analysis.\n"), "red");
		return -1;
	}

	int nb_ref_stars = 0;
	// select reference stars that are only available at least 4/5 of the time
	for (int ref = 1; ref < MAX_SEQPSF && seq->photometry[ref]; ref++) {
		ref_valid[ref] = ref_valid_count[ref] >= round_to_int(nbImages * 4.0 / 5.0);
		siril_debug_print("reference star %d has %d/%d valid measures, %s\n", ref, ref_valid_count[ref], nbImages, ref_valid[ref] ? "including" : "discarding");
		if (ref_valid[ref])
			nb_ref_stars++;
	}

	if (nb_ref_stars == 0) {
		siril_log_color_message(_("The reference stars are not good enough, probably out of the configured valid pixel range, cannot calibrate the light curve\n"), "red");
		return -1;
	}
	if (nb_ref_stars == 1)
		siril_log_color_message(_("Only one reference star was validated, this will not result in an accurate light curve. Try to add more reference stars or check the configured valid pixel range\n"), "salmon");
	else siril_log_message(_("Using %d stars to calibrate the light curve\n"), nb_ref_stars);


	// arrays containing the graph data: X, Y and Y error bars
	double *date = calloc(nbImages, sizeof(double));	// X is the julian date
	double *vmag = calloc(nbImages, sizeof(double));	// Y is the calibrated magnitude
	double *err = calloc(nbImages, sizeof(double));		// Y error bar
	if (!date || !vmag || !err) {
		PRINT_ALLOC_ERR;
		free(date); free(vmag); free(err);
		return -1;
	}
	double min_date = DBL_MAX;
	// i is index in dataset, j is index in output
	for (i = 0, j = 0; i < seq->number; i++) {
		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid)
			continue;

		// X value: the date
		if (seq->imgparam[i].date_obs) {
			double julian;
			GDateTime *tsi = g_date_time_ref(seq->imgparam[i].date_obs);
			if (seq->exposure > 0.0) {
				GDateTime *new_dt = g_date_time_add_seconds(tsi, seq->exposure * 0.5);
				julian = date_time_to_Julian(new_dt);
				g_date_time_unref(new_dt);
			} else {
				julian = date_time_to_Julian(tsi);
			}
			g_date_time_unref(tsi);
			date[j] = julian;
			if (julian < min_date)
				min_date = julian;
		} else {
			date[j] = (double) i + 1; // should not happen.
		}

		// Y value: the magnitude and error and their calibration
		double target_mag = seq->photometry[0][i]->mag;
		double target_err = seq->photometry[0][i]->s_mag;

		double cmag = 0.0, cerr = 0.0;
		int nb_ref = 0;
		/* First data plotted are variable data, others are references
		 * Variable is done above, now we compute references */
		for (int ref = 1; ref < MAX_SEQPSF && seq->photometry[ref]; ref++) {
			if (ref_valid[ref] && seq->photometry[ref][i] && seq->photometry[ref][i]->phot_is_valid) {
				/* inversion of Pogson's law to get back to the flux
				 * Flux = 10^(-0.4 * mag)
				 */
				cmag += pow(10, -0.4 * seq->photometry[ref][i]->mag);
				cerr += seq->photometry[ref][i]->s_mag;
				++nb_ref;
			}
		}
		/* Converting back to magnitude */
		if (nb_ref == nb_ref_stars) {
			/* we consider an image to be invalid if all references are not valid,
			 * because it changes the mean otherwise and makes nonsense data */
			cmag = -2.5 * log10(cmag / nb_ref);
			cerr = (cerr / nb_ref) / sqrt((double) nb_ref);

			vmag[j] = target_mag - cmag;
			err[j] = fmin(9.999, sqrt(target_err * target_err + cerr * cerr));
			j++;
		}
	}
	int nb_valid_images = j;
	int julian0 = 0;
	if (min_date != DBL_MAX)
		julian0 = (int)min_date;

	siril_log_message(_("Calibrated data for %d points of the light curve, %d excluded because of invalid photometry\n"), nb_valid_images, seq->selnum - nb_valid_images);

	gchar *subtitleimg = generate_lc_subtitle(lcargs->metadata, TRUE);
	gchar *titleimg;
	if (!lcargs->target_descr) {
		titleimg = g_strdup_printf("%s %s", _("Light curve of star"), subtitleimg);
	} else {
		titleimg = g_strdup_printf("%s %s%s", _("Light curve of star"), lcargs->target_descr, subtitleimg);
	}
	gchar *subtitledat = generate_lc_subtitle(lcargs->metadata, FALSE);
	gchar *titledat = g_strdup_printf("%s#JD_UT (+ %d)\n", subtitledat, julian0);
	gchar *xlabel = g_strdup_printf("JD_UT (+ %d)", julian0);

	spl_data = malloc(sizeof(siril_plot_data));
	init_siril_plot_data(spl_data);
	siril_plot_set_title(spl_data, titledat);
	siril_plot_set_xlabel(spl_data, xlabel);
	spl_data->revertY = TRUE;
	siril_plot_set_savename(spl_data, "light_curve");
	spl_data->forsequence = TRUE;
	double *date0 = malloc(nb_valid_images * sizeof(double));
	for (int k = 0; k < nb_valid_images; k++)
		date0[k] = date[k] - julian0;
	siril_plot_add_xydata(spl_data, _("V-C"), nb_valid_images, date0, vmag, err, NULL);
	splxyerrdata *lc = (splxyerrdata *)spl_data->plots->data;
	lc->plots[0]->x_offset = (double)julian0;
	free(date0);

	// now we sort to have all dates ascending
	siril_plot_sort_x(spl_data);
	// saving dat
	int ret = 0;
	if (!siril_plot_save_ETD_light_curve(spl_data, filename, TRUE)) {
		ret = 1;
		free_siril_plot_data(spl_data);
		spl_data = NULL; // just in case we try to use it later on
	} else {
		// now saving the plot if required
		siril_plot_set_title(spl_data, titleimg);
		if (!lcargs->display_graph) { // if not used for display we can free spl_data now
			gchar *image_name = replace_ext(filename, ".png");
			siril_plot_save_png(spl_data, image_name, 0, 0);
			free_siril_plot_data(spl_data);
			spl_data = NULL; // just in case we try to use it later on
			g_free(image_name);
		}
	}
	g_free(xlabel);
	g_free(subtitleimg);
	g_free(titleimg);
	g_free(subtitledat);
	g_free(titledat);

	if (lcargs->display_graph && spl_data)
		lcargs->spl_data = spl_data;

	free(date);
	free(vmag);
	free(err);
	return ret;
}

static gboolean end_light_curve_worker(gpointer p) {
	if (sequence_is_loaded()) {
		drawPlot();
		notify_new_photometry();	// switch to and update plot tab
		redraw(REDRAW_OVERLAY);
	}
	return end_generic(NULL);
}

void free_light_curve_args(struct light_curve_args *args) {
	if (args->seq && !check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);
	free(args->areas);
	g_free(args->target_descr);
	if (args->metadata)
		free(args->metadata);
	free(args);
	return;
}

gpointer light_curve_worker(gpointer arg) {
	int retval = 0;
	struct light_curve_args *args = (struct light_curve_args *)arg;

	framing_mode framing = REGISTERED_FRAME;
	if (framing == REGISTERED_FRAME && !args->seq->regparam[args->layer])
		framing = FOLLOW_STAR_FRAME;

	// someday we should move the area in the seqpsf args, not needed for now
	for (int star_index = 0; star_index < args->nb; star_index++) {
		com.selection = args->areas[star_index];

		if (seqpsf(args->seq, args->layer, FALSE, FALSE, framing, FALSE, TRUE)) {
			if (star_index == 0) {
				siril_log_message(_("Failed to analyse the variable star photometry\n"));
				retval = 1;
				break;
			}
			else siril_log_message(_("Failed to analyse the photometry of reference star %d\n"),
					star_index);
		}

		if (args->seq == &com.seq)
			queue_redraw(REDRAW_OVERLAY);
	}
	memset(&com.selection, 0, sizeof(rectangle));

	/* analyse data and create the light curve */
	if (!retval)
		retval = new_light_curve("light_curve.dat", args);
	if (!retval && args->display_graph && args->spl_data) {
		siril_add_idle(create_new_siril_plot_window, args->spl_data);
	}
	free_light_curve_args(args); // this will not free args->spl_data which is free by siril_plot window upon closing
	siril_add_idle(end_light_curve_worker, NULL);
	return GINT_TO_POINTER(retval);
}

////////////////////////////////////////
// All the stuff for occultation
////////////////////////////////////////

static int pulse_detection (struct occultation_args *pulse, double *vmag, int *orig_ind, struct occ_res *timing) {
	gboolean verbose = FALSE;
	int k = 0;
	double sum = 0.0;
	gboolean start_pulse = FALSE;
	gboolean stop_pulse = FALSE;

	/* Loop over the valid images */
	for (int i = 1, j = 0; i + 1 < timing->valid_images; i++) {
		if (verbose) siril_log_color_message(_("1- Valid_images= %i, Vmag= %lf\n"), "red", orig_ind[i + 1], vmag[i]);
		if (!start_pulse && (vmag[i-1] < 0.5 * vmag[i]) && (vmag[i-1] < 0.5 * vmag[i + 1]) && (vmag[i + 1] > 0.85 * timing->hi_val)) {		// Identify raising edges
			start_pulse = TRUE;
			stop_pulse = FALSE;
			pulse[k].start_ind_inseq = orig_ind[i + 1];
			pulse[k].start_ind = i;
			if (verbose) siril_log_color_message(_("Raise!!\n"), "salmon");
		}

		if (start_pulse && (vmag[i-1] > 1.06 * vmag[i]) && (vmag[i-1] > 20.0 * vmag[i + 1]) && (vmag[i + 1] < 0.05 * timing->hi_val)) {		// Identify falling edges
			stop_pulse = TRUE;
			if (verbose) siril_log_color_message(_("Fall!!\n"), "salmon");
		}

		if (start_pulse) {
			if (!stop_pulse) {
				sum += vmag[i];
				j++;
			} else {
				sum += vmag[i];
				start_pulse = FALSE;
				stop_pulse = FALSE;
				j++;	// Number of pictures in a PPS pulse
				pulse[k].pls_nbr = j;
				pulse[k].sum_flux = sum;
				if (verbose) siril_log_color_message(_("Stoppic= %i, Vmag(%lf), Startpic= %i, nbr_pic= %i, sum= %lf\n"), "salmon", orig_ind[i + 1], vmag[orig_ind[i +1]], orig_ind[i +1] - pulse[k].pls_nbr +1, pulse[k].pls_nbr, pulse[k].sum_flux);
				k++;	// Index of the pulse
				j = 0;
				sum = 0.0;
			}
		}
	}	
	return (k - 1);		// number of detected PPS in the sequence
}

static GDateTime *round_date (GDateTime *date){
		/* Check date and time */
	gint year, month, day;

	GTimeZone *tz = g_date_time_get_timezone(date);
	g_date_time_get_ymd(date, &year, &month, &day);
	gint hour = g_date_time_get_hour(date);
	gint minute = g_date_time_get_minute(date);
	double seconds = g_date_time_get_seconds(date);
	double seconds_prime = round(seconds);
	if (seconds_prime == 60.0) {
		minute ++;
		seconds_prime = 0.0;
	}

	GDateTime *ret = g_date_time_new (tz,year, month, day, hour, minute, seconds_prime);
	g_time_zone_unref(tz);
	return ret;
}

static int time_offset (struct occultation_args *pulse, double *vmag, sequence *seq, struct occ_res *timing) {
	gboolean verbose = FALSE;
	double *unity_flux = malloc(timing->th_pls_nbr * sizeof(double));
	double *first_p = malloc(timing->th_pls_nbr * sizeof(double));
	double exposure = timing->exposure;

	for (int i = 0; i < timing->th_pls_nbr; i++) {
		if (i >= timing->det_pulses) continue;	// Avoid an overflow, even if it should not occur here...
		unity_flux[i] = pulse[i].sum_flux / 100.0;	// in ADU/ms, assuming the pulse duration = 100ms
		first_p[i] = vmag[pulse[i].start_ind] / unity_flux[i];	// Lenght of the first pulse (ms)
		if (verbose) siril_log_color_message(_("unit_flux[%i] %lf, sum_flux(%i) = %lf\n"), "red", i, unity_flux[i], i, pulse[i].sum_flux);
		if (verbose) siril_log_color_message(_("Vmag(%i)= %lf, first_p[%i] %lf (ms)\n"), "salmon", pulse[i].start_ind_inseq, vmag[pulse[i].start_ind], i, first_p[i]);

		GDateTime *begin_frame = g_date_time_ref(seq->imgparam[pulse[i].start_ind_inseq - 1].date_obs);	// Timestamp at the beginning of the frame
		GDateTime *end_frame = g_date_time_add_seconds (begin_frame, exposure);	// (computed) Timestamp at the end of the frame
		GDateTime *end_frame2 = g_date_time_add_seconds (end_frame, -1.0 * first_p[i] / 1000.0);	// (Computed) Timestamp of the PPS
		GDateTime *pps_time = round_date (end_frame2);	// (Observed) Timestamp of the precise PPS time the pulse refers to


		gchar *str = date_time_to_FITS_date(begin_frame);
		if (verbose) siril_log_color_message(_("begin_frame %s, exposure= %lf (ms)\n"), "salmon", str, 1000.0 * exposure);
		str = date_time_to_FITS_date(end_frame);
		if (verbose) siril_log_color_message(_("end_frame %s\n"), "salmon", str);
		str = date_time_to_FITS_date(end_frame2);
		if (verbose) siril_log_color_message(_("end_frame2 (supposed PPS) %s, with first_p[%i]= %lf (ms)\n"), "salmon", str, i, first_p[i]);
		str = date_time_to_FITS_date(pps_time);
		if (verbose) siril_log_color_message(_("Real pps_time %s\n"), "salmon", str);
		free(str);

		pulse[i].delay_comp = 1000.0 * timediff_in_s(end_frame2, pps_time);	// delay after the real PPS timestamp
		if (verbose) siril_log_color_message(_("diff time (ms) %lf\n\n"), "salmon", pulse[i].delay_comp);

		g_date_time_unref(begin_frame);
		g_date_time_unref(end_frame);
		g_date_time_unref(end_frame2);
		g_date_time_unref(pps_time);

	}

	/*Compute the median and sigma of the delay to be used later*/
	double *tmp_dbl, lo_val, hi_val;
	tmp_dbl = malloc(timing->det_pulses * sizeof(double));
	for (int i = 0; i < timing->det_pulses; i++)
		tmp_dbl[i] = pulse[i].delay_comp;
	gsl_sort (tmp_dbl, 1, timing->det_pulses);
	timing->median_seq = gsl_stats_median_from_sorted_data (tmp_dbl, 1, timing->det_pulses);
	timing->sig_seq = gsl_stats_variance(tmp_dbl, 1, timing->det_pulses);
	gsl_stats_minmax (&lo_val, &hi_val, tmp_dbl, 1, timing->det_pulses);
	free(tmp_dbl);

	siril_log_color_message(_("Time Offset: %0.3lf(ms), sigma= %0.3lf(ms)\n"), "green", timing->median_seq, timing->sig_seq);
	if (verbose) siril_log_color_message(_("mini diff time %0.3lf(ms), maxi diff time %0.3lf(ms)\n"), "green", lo_val, hi_val);

	free(unity_flux);
	free(first_p);
	return 0;
}




/****************** making a light curve from sequence-stored data ****************/
/* It uses data stored in the sequence, in seq->photometry, which is
 * populated by successive calls to seqpsf on the opened sequence;
 */
void free_occultation_args(struct occultation_args *args) {
	free(args);
	return;
}

void free_occ_res_args(struct occ_res *args) {
	free(args);
	return;
}

int occult_curve(struct light_curve_args *lcargs) {
	int i, j;
	int ret = 0;	// return value

	sequence *seq = lcargs->seq;

	if (!seq->photometry[0]) {
		siril_log_color_message(_("No photometry data found, error\n"), "red");
		return -1;
	}

	// Get number of valid frames, nbImages
	int nbImages = 0;
	for (i = 0; i < seq->number; i++) {
		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid){
			continue;
		}
		++nbImages;
	}
	siril_debug_print("we have %d images with a valid photometry for the variable star\n", nbImages);
	if (nbImages < 1) {
		siril_log_color_message(_("There are not enough valid stars to make a photometric analysis.\n"), "red");
		return -1;
	}

	// Arrays containing the graph data: X, Y and the original index
	double *date = malloc(nbImages * sizeof(double));	// X is the julian date
	double *vmag = malloc(nbImages * sizeof(double));	// Y is the calibrated magnitude
	int *orig_ind = malloc(nbImages * sizeof(int));		// Original index in the sequence
	if (!date || !vmag) {
		PRINT_ALLOC_ERR;
		free(date); free(vmag);
		return -1;
	}

	double min_date = DBL_MAX;
	// i is index in dataset, j is index in output
	for (i = 0, j = 0; i < seq->number; i++) {
		if (!seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid)	// Discard unusable images
			continue;

		// X value: the date
		if (seq->imgparam[i].date_obs) {
			double julian;
			GDateTime *tsi = g_date_time_ref(seq->imgparam[i].date_obs);
			if (seq->exposure > 0.0) {
				GDateTime *new_dt = g_date_time_add_seconds(tsi, seq->exposure * 0.5);
				julian = date_time_to_Julian(new_dt);
				g_date_time_unref(new_dt);
			} else {		// this is the case for ser files. So the timestamp is the beginning of the frame
				julian = date_time_to_Julian(tsi);
			}
			g_date_time_unref(tsi);
			date[j] = julian;
			if (julian < min_date)
				min_date = julian;
		} else {
			date[j] = (double) i + 1; // should not happen.
		}

		// Y value: the the PSF amplitude
		double target_amp = seq->photometry[0][i]->A;
		double target_bck = seq->photometry[0][i]->B;
		vmag[j] = target_amp - target_bck;
		orig_ind[j] = i;	// Keep trace of the original index 
		j++;	// Count for usable images
	
	}

	
	// Retrieve the expected number of pulses in the sequence
	double expos = timediff_in_s(seq->imgparam[0].date_obs, seq->imgparam[seq->number - 1].date_obs);	// Get the time difference in seconds
	int exp_pls_nbr = (int)expos;	// Floor it to get the total expected number of pulses
	expos = 0.01 * expos / (double)exp_pls_nbr;	// Compute the (average) exposition time of the frames

	struct occ_res *timing = NULL;	// Structure embedding the final results
	timing = malloc(sizeof(struct occ_res));
	timing->exposure = expos;
	timing->th_pls_nbr = exp_pls_nbr;
	timing->valid_images = j;

//	siril_log_message(_("Number of PPS we should find in the sequence: %d, exposure(s)= %lf \n"), timing->th_pls_nbr, timing->exposure);
	
	// Retrieve the amplitude range
	gsl_stats_minmax (&timing->lo_val, &timing->hi_val, vmag, 1, timing->valid_images);
//	siril_log_color_message(_("min= %lf, max= %lf, nb_valid_images= %i, nbr pulse= %i\n"), "red", smallest_val, largest_val, nb_valid_images, th_pls_nbr);

	struct occultation_args *pulse = NULL;	// Structure embedding the results for each pulse
	pulse = malloc(timing->th_pls_nbr * sizeof(struct occultation_args));

	// Seeking the pulses in the sequence
	timing->det_pulses = pulse_detection (pulse, vmag, orig_ind, timing);

	// Does it match the forseen value
	siril_log_color_message(_("Detected pulses= %i vs Expected pulses= %i\n"), "salmon", timing->det_pulses, timing->th_pls_nbr - 1);
	if (timing->det_pulses < (timing->th_pls_nbr / 2)) {
		siril_log_color_message(_("Not enought pulses detected. Enlarge the selection.\n"), "red");
		return 0;
	}

	// Computes offset
	ret = time_offset(pulse, vmag, seq, timing);
	lcargs->JD_offset = timing->median_seq;

	// Frees what needs to be
	free_occultation_args(pulse);
	free_occ_res_args(timing);
	free(date);
	free(vmag);

	return ret;
}

/*
void free_occultation_args(struct light_curve_args *args) {
	if (args->seq && !check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);
	free(args->areas);
	if (args->target_descr)
		g_free(args->target_descr);
	if (args->metadata)
		free(args->metadata);
	free(args);
	return;
}
*/
gpointer occultation_worker(gpointer arg) {
	int retval = 0;
	struct light_curve_args *args = (struct light_curve_args *)arg;

	siril_log_message(_("Entering occulation_worker\n"));
	framing_mode framing = REGISTERED_FRAME;

	// Set predifined aperture data to be independant o fthe lasting user parameters 
	com.pref.phot_set.inner = 20;
	com.pref.phot_set.outer = 30;
	com.pref.phot_set.aperture = 10;
	com.pref.phot_set.force_radius = TRUE;

	if (seqpsf(args->seq, args->layer, FALSE, FALSE, framing, FALSE, TRUE)) {
		siril_log_color_message(_("Something went wrong with the PSF analisys. Try to enlarge the selection or check the blinking star is not too faint\n"), "red");
		return GINT_TO_POINTER(retval);
	}

	if (args->seq == &com.seq)
		queue_redraw(REDRAW_OVERLAY);

	memset(&com.selection, 0, sizeof(rectangle));

	/* analyse data and compute the offset */
	if (!retval)
		retval = occult_curve(args);

	if (sequence_is_loaded()) {
		drawPlot();
		notify_new_photometry();	// switch to and update plot tab
		redraw(REDRAW_OVERLAY);
	}

	siril_add_idle(end_occultation_worker, args);
	return GINT_TO_POINTER(retval);
}