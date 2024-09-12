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
#include "gui/utils.h"
#include "io/sequence.h"
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
	if (fit->keywords.cvf > 0.0) {
		// always overwrite with value from the image if available
		retval->gain = fit->keywords.cvf;
	}
	if (fit->type == DATA_FLOAT) {
		retval->gain *= USHRT_MAX_DOUBLE;
		retval->minval /= USHRT_MAX_DOUBLE;
		retval->maxval /= USHRT_MAX_DOUBLE;
	}
	//siril_debug_print("phot_set min=%f, max=%f\n", retval->minval, retval->maxval);
	return retval;
}

void fluxCut_factors (const psf_star *psf, double beta_ref, double fwhm_ref, double* in_rad, double* out_rad, double* ap_rad){
	double threshold = 0.01 * com.pref.phot_set.flux_cut_factor;
	if (psf->profile == PSF_GAUSSIAN) {
		*ap_rad = fwhm_ref *INV_2_SQRT_2_LOG2 * sqrt(-2.0 * log(threshold));
	}
	if (psf->profile == PSF_MOFFAT_BFREE || psf->profile == PSF_MOFFAT_BFIXED) {
		double inv_beta = 1. / beta_ref;
		double const1 = (pow(threshold, -inv_beta) - 1.);
		double const2 = (pow(2., inv_beta) - 1.);
		*ap_rad = 0.5 * fwhm_ref * sqrt(const1 / const2);
	}
	*in_rad = *ap_rad * com.pref.phot_set.flux_inner_factor;
	*out_rad = *ap_rad * com.pref.phot_set.flux_outer_factor;
}

struct radii_set *radii_strat (struct phot_config *phot_set, const psf_star *psf) {
	struct radii_set *retval = malloc(sizeof(struct radii_set));
	if (!retval) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	///*******************************************
///	The goal would be to set fwhm_ref to a particular value:
///	This is the value defines for each star in the currently loaded picture.
///	 It can be the reference image after registration or another one.
	double r11 = 0.0, r21 = 0.0, appRadius1 = 0.0;	// Dump values
	double fwhm_ref = com.pref.phot_set.dump_fwhmx;
	double beta_ref = com.pref.phot_set.dump_beta;
	double in_rad = 0.0, out_rad = 0.0, ap_rad = 0.0;
///*******************************************
///	According to the choosen startegy, computation of the radii:
	switch (phot_set->ape_strat){
		case FIXED_AP:
			r11 = phot_set->inner;
			r21 = phot_set->outer;
			appRadius1 = phot_set->aperture;
			break;
		case FWHM_VAR:	
			r11 = 0.5 * fwhm_ref * phot_set->auto_inner_factor;
			r21 = 0.5 * fwhm_ref * phot_set->auto_outer_factor;
			appRadius1 = 0.5 * fwhm_ref * phot_set->auto_aperture_factor;
			break;
		case FLUX_CUT:
			fluxCut_factors (psf, beta_ref, fwhm_ref, &in_rad, &out_rad, &ap_rad);
			r11 = in_rad;
			r21 = out_rad;
			appRadius1 = ap_rad;
			break;
	}

	if (com.pref.phot_set.dump_fwhmx == 0.0) {	// Workaround for the photometry test in the CI and the first pass of the loop (to set initial values)
		r11 = phot_set->inner;
		r21 = phot_set->outer;
		appRadius1 = phot_set->force_radius ? phot_set->aperture : 0.5 * psf->fwhmx * phot_set->auto_aperture_factor;
	}

	retval->in_Radius = r11;
	retval->out_Radius = r21;
	retval->ape_Radius = appRadius1;

//	siril_debug_print("fwhm_ref: %lf, phot_set->ape_strat: %i \n", fwhm_ref, phot_set->ape_strat);
//	siril_debug_print("phot_set->auto_inner_factor: %lf, phot_set->auto_outer_factor: %lf, phot_set->auto_aperture_factor: %lf \n", phot_set->auto_inner_factor, phot_set->auto_outer_factor, phot_set->auto_aperture_factor);
//	siril_debug_print("Aperture: %lf, Inner: %lf Outer: %lf\n", appRadius1, r11, r21);
//	siril_log_message(_("Inner: %.2lf, Outer: %.2lf, Aperture: %.2lf\n"), r11, r21, appRadius1);

	return retval;
}

photometry *getPhotometryData(gsl_matrix* z, const psf_star *psf,
		struct phot_config *phot_set, gboolean verbose, psf_error *error) {
	int width = z->size2;
	int height = z->size1;
	int n_sky = 0, ret;
	int x, y, x1, y1, x2, y2;
	double r1 = 0.0, r2 = 0.0, r = 0.0, rmin_sq = 0.0, appRadius = 0.0;
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

	struct radii_set *r_set = radii_strat (phot_set, psf);
	r1 = r_set->in_Radius;
	r2 = r_set->out_Radius;
	appRadius = r_set->ape_Radius;

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
	com.pref.phot_set.ape_strat = FWHM_VAR;
	com.pref.phot_set.force_radius = FALSE;
	com.pref.phot_set.auto_inner_factor = 4.2;
	com.pref.phot_set.auto_outer_factor = 6.3;
	com.pref.phot_set.auto_aperture_factor = 4.0;
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
	double *snr_opt = calloc(nbImages, sizeof(double));	// SNR
	if (!date || !vmag || !err || !snr_opt) {
		PRINT_ALLOC_ERR;
		free(date); free(vmag); free(err); free(snr_opt);
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
			snr_opt[j] = seq->photometry[0][i]->SNR;
			j++;
		}
	}
	int nb_valid_images = j;
	int julian0 = 0;

	// Additionnal information on the error bars and variable SNR
	// 
//	if (!lcargs->force_rad) {
	double median_err, largest_err, smallest_err;
	gsl_stats_minmax (&smallest_err, &largest_err, err, 1, nb_valid_images);
	median_err = quickmedian_double(err, nb_valid_images);

	double median_snr, largest_snr, smallest_snr;
	gsl_stats_minmax (&smallest_snr, &largest_snr, snr_opt, 1, nb_valid_images);
	median_snr = quickmedian_double(snr_opt, nb_valid_images);

	siril_log_color_message(_("Error bars-- (%d images) median: %.2lfmmag, max: %.2lfmmag, min: %.2lfmmag, delta: %.2lf\n"), "blue",
		nb_valid_images,
		1000.0 * median_err,
		1000.0 * largest_err,
		1000.0 * smallest_err,
		1000.0 * (largest_err - smallest_err));
	siril_log_color_message(_("Variable star SNR-- (%d images) median: %.2lfdB, max: %.2lfdB, min: %.2lfdB\n"), "blue",
		nb_valid_images,
		median_snr,
		largest_snr,
		smallest_snr);
//	}

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


// Gets the fwhmx of the selected star in the currently loaded image only
int one_psf(int star_index) {
	if (!check_ok_if_cfa())
		return 1;
	psf_star *result = NULL;
	int layer = select_vport(gui.cvport);

	if (layer == -1)
		return 1;
	if (!(com.selection.h && com.selection.w))
		return 1;

	struct phot_config *ps = phot_set_adjusted_for_image(&gfit);
	int ape_strat_bkp = com.pref.phot_set.ape_strat;
	ps->ape_strat = FIXED_AP;
	result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE, ps, TRUE, com.pref.starfinder_conf.profile, NULL);
	ps->fwhm_ref[star_index] = max(result->fwhmx, result->fwhmy);
	ps->beta_ref[star_index] = result->beta;
	com.pref.phot_set.fwhm_ref[star_index] = max(result->fwhmx, result->fwhmy);
	com.pref.phot_set.beta_ref[star_index] = result->beta;
	siril_debug_print("ICI-result->fwhmx: %lf\n", result->fwhmx);
	ps->ape_strat = ape_strat_bkp;
	free(ps);
	if (!result) 
		return 0;

//	popup_psf_result(result, &com.selection, &gfit);
	free_psf(result);
	return 1;
}

gpointer light_curve_worker(gpointer arg) {
	int retval = 0;
	struct light_curve_args *args = (struct light_curve_args *)arg;
	framing_mode framing = REGISTERED_FRAME;
	if (framing == REGISTERED_FRAME && !args->seq->regparam[args->layer])
		framing = FOLLOW_STAR_FRAME;
	// someday we should move the area in the seqpsf args, not needed for now


//	Retrieves in an array, the fwhm of all the stars in the list. O beeing the target star
	for (int star_index = 0; star_index < args->nb; star_index++) {
		//Performed over the current image
		com.selection = args->areas[star_index];
		if(!one_psf(star_index)) {
			siril_log_message(_("Failed to analyse the photometry of reference star %d\n"),
				star_index);
		}
		siril_debug_print("LA-com.pref.phot_set.fwhm_ref[%d]: %lf\n", star_index, com.pref.phot_set.fwhm_ref[star_index]);
		if (args->seq == &com.seq)
			queue_redraw(REDRAW_OVERLAY);
	}


	for (int star_index = 0; star_index < args->nb; star_index++) {
		com.selection = args->areas[star_index];
		com.pref.phot_set.dump_fwhmx = com.pref.phot_set.fwhm_ref[star_index];	// Sets the fwhm_ref for the seqpsf process for this particular star
		com.pref.phot_set.dump_beta = com.pref.phot_set.beta_ref[star_index];	// Sets the fwhm_ref for the seqpsf process for this particular star
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
	args->force_rad = com.pref.phot_set.force_radius;	// Retrieve the Aperture state (fixed/dynamic)
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

