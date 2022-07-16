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

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "algos/PSF.h"
#include "algos/photometry.h"
#include "algos/statistics_float.h"
#include "io/sequence.h"
#include "opencv/opencv.h"
#include "gui/PSF_list.h"
#include "gui/gnuplot_i/gnuplot_i.h"

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
	cvTransfPoint(&x, &y, H, Href);

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

	double xc = psf->x0 - 1;
	double yc = psf->y0 - 1;

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
			return _("Gaussian fit failed");
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

/********************** PSF analysis for all stars in the image ****************/

struct all_stars_struct {
	int ref_image;		// index of the image used to find stars below
	regdata *regparam;
	psf_star **stars_in_ref;// non-saturated stars in first image, provider by peaker
	int nbstars;		// number of stars in stars_in_ref
	//siril_wcs **coords;	// their equatorial J2000 coords		TODO

	gint *valid_counts;	// allocated to nbstars size, number of images a star is valid in
	int validity_thresh;	// the threshold in number of images where a star is valid to then use it

	psf_star ***photo;	// measurements: array of stars per image, only valid or NULL

	int min_number_of_ref;	// minimum allowed number of reference stars to be considered calibrated
	double mag_delta_lim;	// the value around a mag to create a valid range for reference stars
	int distance_thresh;	// upper threshold for pixel distance between a star and its reference
	int **reference_stars;	// indices of reference stars found for each star
	BYTE *ref_stars_count;	// number of reference stars for each star, size of *reference_stars

	//int julian0;		// the day of the start of the sequence		TODO
	gboolean use_gnuplot;	// if FALSE, just create data files
};

int all_stars_psf_image_hook(struct generic_seq_args *args, int out_index, int index, fits *fit, rectangle *area, int threads) {
	struct all_stars_struct *asargs = (struct all_stars_struct *)args->user;
	asargs->photo[out_index] = calloc(asargs->nbstars, sizeof(psf_star *));
	if (!asargs->photo[out_index]) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	Homography Href = asargs->regparam[asargs->ref_image].H;
	struct phot_config *ps = phot_set_adjusted_for_image(fit);

	// save the date for X axis (seqpsf does that too)
	if (!args->seq->imgparam[index].date_obs && fit->date_obs)
		args->seq->imgparam[index].date_obs = g_date_time_ref(fit->date_obs);

	int nb_valid = 0;
	psf_error errors[PSF_ERR_MAX_VALUE] = { 0 };
	for (int i = 0; i < asargs->nbstars; i++) {
		psf_star *star = asargs->stars_in_ref[i];
		struct phot_config phot_set;
		Homography H = asargs->regparam[index].H;
		rectangle area = compute_dynamic_area_for_psf(star, ps, &phot_set, H, Href);
		if (enforce_area_in_image(&area, args->seq, index)) {
			errors[PSF_ERR_OUT_OF_WINDOW]++;
			continue;
		}

		psf_error error;
		psf_star *phot_psf = psf_get_minimisation(fit, 0, &area, FALSE, TRUE, &phot_set, FALSE, &error);
		if (phot_psf) {
			errors[error]++;
			if (phot_psf->phot_is_valid && phot_psf->SNR > 5.0) {
				g_atomic_int_inc(&asargs->valid_counts[i]);
				asargs->photo[out_index][i] = phot_psf;
				nb_valid++;
			}
			else free_psf(phot_psf);
		}
	}
	if (nb_valid < (asargs->nbstars * 2 / 3))
		print_psf_error_summary((gint*)errors);
	siril_debug_print("Got %d stars with valid photometry in image %d\n", nb_valid, index);
	free(ps);
	return 0;
}

int all_stars_psf_finalize_hook(struct generic_seq_args *args) {
	struct all_stars_struct *asargs = (struct all_stars_struct *)args->user;
	int nb_useful_stars = 0;
	for (int i = 0; i < asargs->nbstars; i++) {
		siril_debug_print("star %d was valid in %d images (%.1f %%, rel mag: %.1f)\n",
				i, asargs->valid_counts[i],
				asargs->valid_counts[i] / (double)args->seq->selnum * 100.0,
				asargs->stars_in_ref[i]->mag);
		if (asargs->valid_counts[i] > asargs->validity_thresh)
			nb_useful_stars++;
	}
	siril_log_message(_("%d stars had a valid photometry throughout the sequence (on at least %d images)\n"), nb_useful_stars, asargs->validity_thresh);
	return 0;
}

int all_stars_psf(sequence *seq, int layer, struct all_stars_struct *asargs) {
	struct generic_seq_args *args = create_default_seqargs(seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = seq->selnum;
	args->image_hook = all_stars_psf_image_hook;
	args->finalize_hook = all_stars_psf_finalize_hook;
	args->description = "all stars PSF";
	args->already_in_a_thread = TRUE;
	asargs->valid_counts = calloc(asargs->nbstars, sizeof(gint));
	if (!asargs->valid_counts) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	asargs->photo = malloc(seq->selnum * sizeof(psf_star **));
	if (!asargs->photo) {
		PRINT_ALLOC_ERR;
		free(asargs->valid_counts);
		asargs->valid_counts = NULL;
		return -1;
	}
	asargs->reference_stars = calloc(asargs->nbstars, sizeof(int *));
	if (!asargs->reference_stars) {
		PRINT_ALLOC_ERR;
		free(asargs->valid_counts);
		asargs->valid_counts = NULL;
		free(asargs->photo);
		asargs->photo = NULL;
		return -1;
	}
	asargs->ref_stars_count = malloc(asargs->nbstars * sizeof(BYTE));
	if (!asargs->ref_stars_count) {
		PRINT_ALLOC_ERR;
		free(asargs->valid_counts);
		asargs->valid_counts = NULL;
		free(asargs->photo);
		asargs->photo = NULL;
		free(asargs->reference_stars);
		asargs->reference_stars = NULL;
		return -1;
	}
	args->user = asargs;

	return GPOINTER_TO_INT(generic_sequence_worker(args));
}

struct ref_star {
	int distance;	// in pixels
	int index;	// in asargs->stars
};

int star_distance_cmp(const void *a, const void *b) {
	struct ref_star *s1 = (struct ref_star *)a;
	struct ref_star *s2 = (struct ref_star *)b;
	if (s1->distance > s2->distance)
		return 1;
	if (s1->distance < s2->distance)
		return -1;
	return 0;
}

static int pixel_distance(psf_star *s1, psf_star *s2) {
	double dx = s1->xpos - s2->xpos;
	double dy = s1->ypos - s2->ypos;
	return round_to_int(sqrt(dx * dx + dy * dy));
}

#define MAX_REF_STARS 30	// uchar
int find_reference_stars(struct all_stars_struct *asargs) {
	struct ref_star *ref_candidates = malloc(asargs->nbstars * sizeof(struct ref_star));
	if (!ref_candidates) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	for (int i = 0; i < asargs->nbstars; i++) {
		/* for each star, find a list of stars of similar magnitude and not too
		 * far in the image if possible. With WCS and a catalogue/request, we
		 * could filter out the variable stars too */
		if (asargs->valid_counts[i] < asargs->validity_thresh) {
			asargs->reference_stars[i] = NULL;
			asargs->ref_stars_count[i] = 0;
			continue;
		}
		psf_star *star = asargs->stars_in_ref[i];
		int nb_refs = 0;
		for (int j = 0; j < asargs->nbstars; j++) {
			if (i == j) continue;
			if (asargs->valid_counts[j] < asargs->validity_thresh)
				continue;
			psf_star *ref = asargs->stars_in_ref[j];
			// mag here is the estimate, not a photometric measure, but always available
			if (fabs(ref->mag - star->mag) < asargs->mag_delta_lim) {
				ref_candidates[nb_refs].index = j;
				ref_candidates[nb_refs].distance = pixel_distance(ref, star);
				nb_refs++;
			}
		}

		if (nb_refs < asargs->min_number_of_ref) {
			asargs->reference_stars[i] = NULL;
			asargs->ref_stars_count[i] = 0;
			siril_debug_print("star %d didn't have enough valid reference stars\n", i);
			continue;
		}
		int kept_refs = min(MAX_REF_STARS, nb_refs);
		asargs->reference_stars[i] = malloc(kept_refs * sizeof(int));
		if (!asargs->reference_stars[i]) {
			asargs->ref_stars_count[i] = 0;
			PRINT_ALLOC_ERR;
			break;
		}
		asargs->ref_stars_count[i] = kept_refs;

		// sort by distance and keep the MAX_REF_STARS best
		qsort(ref_candidates, nb_refs, sizeof(struct ref_star), star_distance_cmp);
		for (int j = 0; j < kept_refs; j++) {
			if (ref_candidates[j].distance <= asargs->distance_thresh || j < asargs->min_number_of_ref)
				asargs->reference_stars[i][j] = ref_candidates[j].index;
		}
		siril_debug_print("found %d reference stars for star %d\n", kept_refs, i);
	}
	free(ref_candidates);
	return 0;
}

static int create_light_curve(sequence *seq, struct all_stars_struct *asargs, int star_to_plot);

gpointer crazy_photo_worker(gpointer arg) {
	int reglayer = GPOINTER_TO_INT(arg);
	int layer;
	if (com.headless)
		layer = 0;
	else layer = gui.cvport == RGB_VPORT ? GLAYER : gui.cvport;

	// 1. detect stars
	image im = { .fit = &gfit, .from_seq = &com.seq, .index_in_seq = com.seq.current };

	int nbstars = 0;
	rectangle *selection = NULL;
	if (com.selection.w != 0 && com.selection.h != 0)
		selection = &com.selection;
	psf_star **stars = peaker(&im, layer, &com.pref.starfinder_conf, &nbstars, selection, FALSE, FALSE, -1, com.max_thread);
	if (!stars) {
		siril_log_message(_("Cannot run photometry analysis if stars are not found\n"));
		return GINT_TO_POINTER(1);
	}

	// 2. run seqpsf on all non-saturated
	psf_star **photo_stars = malloc((nbstars + 1) * sizeof(psf_star*));
	if (!photo_stars) {
		PRINT_ALLOC_ERR;
		return GINT_TO_POINTER(1);
	}
	// filter them a little first, only not saturated
	struct phot_config *ps = phot_set_adjusted_for_image(&gfit);
	if (!ps) return GINT_TO_POINTER(1);
	int valid_candidates = 0;
	for (int i = 0; i < nbstars; i++) {
		if (stars[i]->A > ps->maxval) {
			free_psf(stars[i]);
			continue;
		}
		// magnitude check? SNR check? we don't have reliable photometric info yet
		photo_stars[valid_candidates++] = stars[i];
	}
	photo_stars[valid_candidates] = NULL;	// for display as com.stars
	free(stars);
	photo_stars = realloc(photo_stars, (valid_candidates + 1) * sizeof(psf_star*));
	siril_log_message(_("Found %d photometry candidates in %s, channel #%d\n"), valid_candidates,
			selection ? _("selection") : _("image"), layer);
	if (valid_candidates > 5) {
		clear_stars_list(FALSE);
		com.stars = photo_stars;
	} else {
		siril_log_message(_("Not enough stars found in the reference image to do a photometry analysis\n"));
		return GINT_TO_POINTER(1);
	}

	struct all_stars_struct asargs = {
		.ref_image = com.seq.current,
		.regparam = com.seq.regparam[reglayer],
		.stars_in_ref = photo_stars,
		.nbstars = valid_candidates,
		.valid_counts = NULL,
		.validity_thresh = com.seq.selnum * 4 / 5,
		.photo = NULL,
		.min_number_of_ref = 5,
		.mag_delta_lim = 0.6,
		.distance_thresh = com.seq.rx / 2,
		.reference_stars = NULL,
		.ref_stars_count = NULL,
	};
	int retval = all_stars_psf(&com.seq, layer, &asargs);
	if (retval)
		return GINT_TO_POINTER(retval);

	// 3. magic happens
	find_reference_stars(&asargs);
	asargs.use_gnuplot = FALSE; //gnuplot_is_available();

	for (int i = 0; i < asargs.nbstars; i++) {
		// filter down references based on common availability along the sequence?

		// make a light curve
		create_light_curve(&com.seq, &asargs, i);
	}

	free(asargs.valid_counts);
	for (int i = 0; i < com.seq.selnum; i++) {
		if (asargs.photo[i]) {
			for (int j = 0; j < asargs.nbstars; j++) {
				if (asargs.photo[i][j])
					free_psf(asargs.photo[i][j]);
			}
			free(asargs.photo[i]);
		}
	}
	free(asargs.photo);
	for (int i = 0; i < asargs.nbstars; i++)
		free(asargs.reference_stars[i]);
	free(asargs.reference_stars);
	free(asargs.ref_stars_count);
	return NULL;
}

/************************** headless light curves ***************************
 * light curves are historically made in gui/plot.c by filling the kplot data
 * structure with X and Y data then calling a function that will identify the
 * data in all plots, calibrate the data using the reference stars and produce
 * a .dat file and a plot calling gnuplot.
 * Here we pass directly the sequence and psf_star data, do the same kind of
 * calibration and generate the .dat file.
 ****************************************************************************/
struct photo_data_point {
	double julian_date;
	double mag;
	double mag_err;
};

static int write_photometry_data_file(const char *filename, struct photo_data_point *points, int nb_points, gboolean use_gnuplot);

static int create_light_curve(sequence *seq, struct all_stars_struct *asargs, int star_to_plot) {
	if (!asargs->reference_stars[star_to_plot])
		return 1;
	struct photo_data_point *points = malloc(seq->selnum * sizeof(struct photo_data_point));
	if (!points) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	int j = 0;	// index in points
	for (int i = 0; i < seq->selnum; i++) {
		if (!seq->imgparam[i].date_obs)
			continue;
		double julian;
		GDateTime *tsi = g_date_time_ref(seq->imgparam[i].date_obs);
		if (seq->exposure) {
			GDateTime *new_dt = g_date_time_add_seconds(tsi, seq->exposure / 2.0);
			julian = date_time_to_Julian(new_dt);
			g_date_time_unref(new_dt);
		} else {
			julian = date_time_to_Julian(tsi);
		}

		// X axis: value is 'julian'
		//plot->julian[j] = julian - (double)julian0;

		if (!asargs->photo[i][star_to_plot])
			continue;
		double mag = asargs->photo[i][star_to_plot]->mag;
		double err = asargs->photo[i][star_to_plot]->s_mag;

		double cmag = 0.0, cerr = 0.0;
		int ref, nb_ref = asargs->ref_stars_count[star_to_plot];

		for (ref = 0; ref < nb_ref; ref++) {
			int ref_index = asargs->reference_stars[star_to_plot][ref];
			if (!asargs->photo[i][ref_index])
				break;
			if (!asargs->photo[i][ref_index]->phot_is_valid) {
				siril_debug_print("Invalid photometry in phot\n");
				return 1;
			}

			double rmag = asargs->photo[i][ref_index]->mag;
			double rerr = asargs->photo[i][ref_index]->s_mag;
			/* inversion of Pogson's law: Flux = 10^(-0.4 * mag) */
			cmag += pow(10, -0.4 * rmag);
			cerr += rerr;
		}
		/* we consider an image to be invalid if all references are not valid,
		 * because it changes the mean otherwise and makes nonsense data */
		if (ref != nb_ref)
			continue;

		/* Converting back to magnitude */
		cmag = -2.5 * log10(cmag / nb_ref);
		cerr = (cerr / nb_ref) / sqrt((double) nb_ref);

		mag -= cmag;
		err = fmin(9.999, sqrt(err * err + cerr * cerr));

		// Y axis: value is mag, error is err
		points[j].julian_date = julian;
		points[j].mag = mag;
		points[j].mag_err = err;
		j++;

		g_date_time_unref(tsi);
	}

	gchar *plot_file = g_strdup_printf("plot_%05d.dat", star_to_plot);
	int retval = write_photometry_data_file(plot_file, points, j, asargs->use_gnuplot);
	if (retval > 0)
		asargs->use_gnuplot = FALSE;
	g_free(plot_file);
	free(points);
	return retval;
}

int write_photometry_data_file(const char *filename, struct photo_data_point *points, int nb_points, gboolean use_gnuplot) {
	if (nb_points < 1)
		return -2;

	FILE *fd = g_fopen(filename, "w");
	if (!fd) {
		perror("creating data file");
		return -1;
	}

	// TODO: write star coordinates or even name some day
	fprintf(fd, "# julian_date magnitude error_bar_of_mag\n");


	for (int i = 0; i < nb_points; i++) {
		fprintf(fd, "%f\t%f\t%f\n", points[i].julian_date, points[i].mag, points[i].mag_err);
		//fprintf(fileHandle, "%14.6f %8.6f %8.6f\n", x[i], y[i], yerr[i]);
	}

	fclose(fd);

	if (use_gnuplot) {
		gnuplot_ctrl *gplot = NULL;
		/* First, close the graph if already exists */
		if (gplot != NULL) {
			gnuplot_close(gplot);
		}

		if ((gplot = gnuplot_init())) {
			/* Plotting light curve */
			gnuplot_set_title(gplot, _("Light Curve"));
			gnuplot_set_xlabel(gplot, "JD");
			gnuplot_reverse_yaxis(gplot);
			gnuplot_setstyle(gplot, "errorbars");
			gnuplot_plot_xyyerr_from_datfile(gplot, filename);
		}
		else return 1;
	}
	return 0;
}
