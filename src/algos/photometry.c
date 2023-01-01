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

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
#include "algos/PSF.h"
#include "algos/photometry.h"
#include "algos/statistics_float.h"
#include "algos/siril_wcs.h"
#include "io/sequence.h"
#include "opencv/opencv.h"
#include "gui/PSF_list.h"
#include "gui/plot.h"
#include "gui/image_display.h"
#include "io/gnuplot_i.h"

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

struct phot_config *phot_set_adjusted_for_image(fits *fit) {
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

// dynamic radius: aperture = 2*fwhm, inner = 3*fwhm, outer = 4*fwhm
rectangle compute_dynamic_area_for_psf(psf_star *psf, struct phot_config *original, struct phot_config *phot_set, Homography H, Homography Href) {
	phot_set->gain = original->gain;
	phot_set->aperture = psf->fwhmx * 2.0;
	phot_set->inner = psf->fwhmx * 3.0;
	phot_set->outer = psf->fwhmx * 4.0;
	phot_set->force_radius = TRUE;
	phot_set->minval = original->minval;
	phot_set->maxval = original->maxval;

	double start = 1.5 * phot_set->outer;
	double size = 3 * phot_set->outer;
	double x = psf->xpos, y = psf->ypos;
	cvTransfPoint(&x, &y, Href, H);

	rectangle area = {
		.x = x - start,
		.y = y - start,
		.w = size,
		.h = size
	};
	return area;
}

/* Function that compute all photometric data. The result must be freed */
photometry *getPhotometryData(gsl_matrix* z, psf_star *psf,
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

/* used only be the old (libconfig) settings file reading and unit tests */
void initialize_photometric_param() {
	com.pref.phot_set.inner = 20;
	com.pref.phot_set.outer = 30;
	com.pref.phot_set.aperture = 10;
	com.pref.phot_set.force_radius = FALSE;
	com.pref.phot_set.gain = 2.3;
	com.pref.phot_set.minval = 0;
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
		case PSF_ERR_OUT_OF_IMAGE:
			return _("not in image");
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

/****************** making a light curve from sequence-stored data ****************/
/* contrary to light_curve() in gui/plot.c, this one does not use preprocessed data and the
 * kplot data structure. It only uses data stored in the sequence, in seq->photometry, which is
 * populated by successive calls to seqpsf on the opened sequence;
 * TODO: replace light_curve() by this?
 */
int new_light_curve(sequence *seq, const char *filename, const char *target_descr, gboolean display_graph) {
	int i, j;
	gboolean use_gnuplot = gnuplot_is_available();
	if (!use_gnuplot) {
		siril_log_message(_("Gnuplot was not found, the light curve data will be produced in %s but no image will be created.\n"), filename);
	}
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
	if (nbImages < 1)
		return -1;

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
	if (nb_ref_stars < 1)
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

	/* Exporting data in a dat file */
	int ret = gnuplot_write_xyyerr_dat(filename, date, vmag, err, nb_valid_images, "JD_UT V-C err");
	if (ret) {
		if (com.script)
			siril_log_color_message(_("Failed to create the light curve data file %s\n"), "red", filename);
	} else {
		siril_log_message(_("%s has been saved.\n"), filename);

		/*  data are computed, now plot the graph. */
		if (use_gnuplot) {
			gnuplot_ctrl *gplot = gnuplot_init();
			if (gplot) {
				/* Plotting light curve */
				gchar *title = g_strdup_printf("Light curve of star %s", target_descr);
				gnuplot_set_title(gplot, title);
				gchar *xlabel = g_strdup_printf("Julian date (+ %d)", julian0);
				gnuplot_set_xlabel(gplot, xlabel);
				gnuplot_reverse_yaxis(gplot);
				gnuplot_setstyle(gplot, "errorbars");
				if (display_graph) {
					gnuplot_plot_xyyerr_from_datfile(gplot, filename, "relative magnitude", julian0);
				} else {
					gchar *image_name = replace_ext(filename, ".png");
					gnuplot_plot_datfile_to_png(gplot, filename, "relative magnitude", julian0, image_name);
					siril_log_message(_("%s has been generated.\n"), image_name);
					g_free(image_name);
				}
				g_free(title);
				g_free(xlabel);
			}
			else siril_log_message(_("Communicating with gnuplot failed, still creating the data file\n"));
		}
	}

	free(date);
	free(vmag);
	free(err);
	return ret;
}

static int get_photo_area_from_ra_dec(fits *fit, double ra, double dec, rectangle *ret_area) {
	double x, y;
	if (wcs2pix(fit, ra, dec, &x, &y)) {
		siril_debug_print("star is outside image\n");
		return 1;
	}
	y = fit->ry - y - 1;
	double start = 1.5 * com.pref.phot_set.outer;
	double size = 3 * com.pref.phot_set.outer;
	rectangle area;
	area.x = x - start;
	area.y = y - start;
	area.w = size;
	area.h = size;
	if (area.x < 0 || area.y < 0 ||
			area.h <= 0 || area.w <= 0 ||
			area.x + area.w >= fit->rx ||
			area.y + area.h >= fit->ry) {
		siril_debug_print("star is outside image\n");
		return 1;
	}
	*ret_area = area;
	siril_debug_print("Pixel coordinates of a star: %.1f, %.1f\n", x, y);
	return 0;
}

// area is not recovering another at least half their size (assumed square and identical)
static int area_is_unique(rectangle *area, rectangle *areas, int nb_areas) {
	int half_size = area->w / 2;
	for (int i = 0; i < nb_areas; i++) {
		if (abs(area->x - areas[i].x) < half_size && abs(area->y - areas[i].y) < half_size)
			return 0;
	}
	return 1;
}

int parse_nina_stars_file_using_WCS(struct light_curve_args *args, const char *file_path, gboolean use_comp1, gboolean use_comp2, fits *first) {
	/* The file is a CSV with these fields:
	 * Type,Name,HFR,xPos,yPos,AvgBright,MaxBright,Background,Ra,Dec
	 *
	 * Type can be 'Target' for the variable star to analyse, 'Var' for variable stars to
	 * absolutely exclude as calibration reference, 'Comp1' are reference stars obtained
	 * from SIMBAD based on color, 'Comp2' are stars obtained from the AAVSO site.
	 * We just need this and Ra,Dec. xPos and yPos are in pixels, but we'll use wcs2pix
	 * to get them with our plate solve and our star fitting.
	 */
	FILE *fd = g_fopen(file_path, "r");
	if (!fd) {
		siril_log_message(_("Could not open file %s: %s\n"), file_path, strerror(errno));
		return 1;
	}
	char buf[512];
	rectangle *areas = malloc(MAX_REF_STARS * sizeof(rectangle));
	areas[0].x = 0; areas[0].y = 0;
	int ra_index = -1, dec_index = -1, name_index = 2;
	int stars_count = 0;
	gboolean ready_to_parse = FALSE, target_acquired = FALSE;
	while (fgets(buf, 512, fd)) {
		if (buf[0] == '\0' || buf[0] == '\r' || buf[0] == '\n' || buf[0] == '#')
			continue;
		remove_trailing_eol(buf);
		gchar **tokens = g_strsplit(buf, ",", -1);
		int length = g_strv_length(tokens);
		if (!tokens[0] || length <= ra_index || length <= dec_index) {
			siril_debug_print("malformed line: %s\n", buf);
			g_strfreev(tokens);
			continue;
		}
		gchar *type = tokens[0];
		if (!ready_to_parse) {
			if (!strcasecmp(type, "type")) {
				siril_debug_print("header from the NINA file: %s\n", buf);
				for (int i = 1; tokens[i]; i++) {
					if (!strcasecmp(tokens[i], "ra"))
						ra_index = i;
					else if (!strcasecmp(tokens[i], "dec"))
						dec_index = i;
					else if (!strcasecmp(tokens[i], "name"))
						name_index = i;
				}
				g_strfreev(tokens);

				if (ra_index < 1 || dec_index < 1) {
					siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
					fclose(fd);
					return 1;
				}
				siril_debug_print("Found RA and Dec indices in file: %d and %d\n", ra_index, dec_index);
				ready_to_parse = TRUE;
				continue;
			}
			else {
				siril_debug_print("malformed line: %s\n", buf);
				siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
				g_strfreev(tokens);
				fclose(fd);
				return 1;
			}
		}

		if (!strcasecmp(type, "target")) {
			gchar *end1, *end2;
			double ra = g_ascii_strtod(tokens[ra_index], &end1);
			double dec = g_ascii_strtod(tokens[dec_index], &end2);
			if (end1 == tokens[ra_index] || end2 == tokens[dec_index]) {
				siril_debug_print("malformed line: %s\n", buf);
				siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
				g_strfreev(tokens);
				fclose(fd);
				return 1;
			}

			args->target_descr = g_strdup(tokens[name_index]);
			if (!get_photo_area_from_ra_dec(first, ra, dec, &areas[0])) {
				target_acquired = TRUE;
				stars_count++;
				siril_log_message(_("Target star identified: %s\n"), tokens[name_index]);
			} else {
				siril_log_message(_("There was a problem finding the target star in the image, cannot continue with the light curve\n"));
			}
		}
		else if (!strcasecmp(type, "var")) {
			// we don't use them for this, but we could add them in the
			// user catalogue for annotations, or a local database
		}
		else if (!strcasecmp(type, "comp1")) {
			if (!use_comp1) {
				siril_debug_print("ignoring comp1 star\n");
				g_strfreev(tokens);
				continue;
			}
			gchar *end1, *end2;
			double ra = g_ascii_strtod(tokens[ra_index], &end1);
			double dec = g_ascii_strtod(tokens[dec_index], &end2);
			if (end1 == tokens[ra_index] || end2 == tokens[dec_index]) {
				siril_debug_print("malformed line: %s\n", buf);
				siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
				g_strfreev(tokens);
				fclose(fd);
				return 1;
			}
			int index = target_acquired ? stars_count : stars_count + 1;
			if (!get_photo_area_from_ra_dec(first, ra, dec, &areas[index])) {
				if (area_is_unique(&areas[index], areas, index)) {
					stars_count++;
					siril_log_message(_("Star %s added as a reference star\n"), tokens[name_index]);
				}
				else siril_log_message(_("Star %s ignored because it was too close to another\n"), tokens[name_index]);
			}
			else siril_log_message(_("Star %s could not be used because it's on the borders or outside\n"), tokens[name_index]);
		}
		else if (!strcasecmp(type, "comp2")) {
			if (!use_comp2) {
				siril_debug_print("ignoring comp2 star\n");
				g_strfreev(tokens);
				continue;
			}
			gchar *end1, *end2;
			double ra = g_ascii_strtod(tokens[ra_index], &end1);
			double dec = g_ascii_strtod(tokens[dec_index], &end2);
			if (end1 == tokens[ra_index] || end2 == tokens[dec_index]) {
				siril_debug_print("malformed line: %s\n", buf);
				siril_log_message(_("The NINA star information file did not contain all expected data (RA and Dec)\n"));
				g_strfreev(tokens);
				fclose(fd);
				return 1;
			}
			int index = target_acquired ? stars_count : stars_count + 1;
			if (!get_photo_area_from_ra_dec(first, ra, dec, &areas[index])) {
				if (area_is_unique(&areas[index], areas, index)) {
					stars_count++;
					siril_log_message(_("Star %s added as a reference star\n"), tokens[name_index]);
				}
				else siril_log_message(_("Star %s ignored because it was too close to another\n"), tokens[name_index]);
			}
			else siril_log_message(_("Star %s could not be used because it's on the borders or outside\n"), tokens[name_index]);
		}
		else {
			siril_debug_print("malformed line: %s\n", buf);
		}
		g_strfreev(tokens);
		if (stars_count >= MAX_REF_STARS)
			break;
	}
	if (target_acquired) {
		args->areas = areas;
		args->nb = stars_count;
	}
	fclose(fd);
	return !target_acquired;
}

static gboolean end_light_curve_worker(gpointer p) {
	drawPlot();
	notify_new_photometry();	// switch to and update plot tab
	redraw(REDRAW_OVERLAY);
	return end_generic(NULL);
}

gpointer light_curve_worker(gpointer arg) {
	int retval = 0;
	struct light_curve_args *args = (struct light_curve_args *)arg;

	framing_mode framing = REGISTERED_FRAME;
	if (framing == REGISTERED_FRAME && !args->seq->regparam[args->layer])
		framing = FOLLOW_STAR_FRAME;

	/* for now, we use seqpsf as many times as needed and the GUI way of
	 * generating the light curve. Maybe someday it would be wise to move to
	 * all_stars_psf instead, depending on the number of reference stars */
	for (int star_index = 0; star_index < args->nb; star_index++) {
		com.selection = args->areas[star_index];

		if (seqpsf(args->seq, args->layer, FALSE, FALSE, framing, FALSE, TRUE)) {
			if (star_index == 0) {
				siril_log_message(_("Failed to analyse the variable star photometry\n"));
				retval = 1;
				break;
			}
			else siril_log_message(_("Failed to analyse the photometry of reference star %d\n"), star_index);
		}

		if (args->seq == &com.seq)
			queue_redraw(REDRAW_OVERLAY);
	}

	/* analyse data and create the light curve */
	if (!retval)
		retval = new_light_curve(args->seq, "light_curve.dat", args->target_descr, args->display_graph);

	if (args->seq != &com.seq)
		free_sequence(args->seq, TRUE);
	free(args);
	siril_add_idle(end_light_curve_worker, NULL);
	return GINT_TO_POINTER(retval);
}

