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
//#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
//#include <gsl/gsl_statistics_double.h>
#include <string.h>

//#include "core/siril.h"
//#include "core/proto.h"
#include "core/siril_date.h"
#include "core/siril_log.h"
//#include "algos/sorting.h"
//#include "algos/PSF.h"
#include "algos/photometry.h"
#include "algos/occult_time.h"
//#include "algos/astrometry_solver.h"
//#include "algos/statistics_float.h"
//#include "algos/siril_wcs.h"
//#include "algos/search_objects.h"
//#include "algos/comparison_stars.h"
//#include "gui/PSF_list.h"
//#include "gui/plot.h"
#include "gui/image_display.h"
//#include "gui/siril_plot.h"
#include "gui/occultation.h"
#include "io/sequence.h"
#include "io/ser.h"
//#include "io/siril_plot.h"
//#include "opencv/opencv.h"


#define PULSE_WIDTH	100.0	// Pulse width (ms)

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
		if (!start_pulse && (vmag[i-1] < 0.5 * vmag[i]) && (vmag[i-1] < 0.5 * vmag[i + 1]) && (vmag[i + 1] > 0.55 * timing->hi_val)) {		// Identify raising edges
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
				if (sum == 0.0) return 0;
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

	GTimeZone *tz = g_time_zone_new_utc();
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

static int time_offset_calc (struct occultation_args *pulse, double *vmag, sequence *seq, struct occ_res *timing) {
	gboolean verbose = FALSE;
	double *unity_flux = malloc(timing->th_pls_nbr * sizeof(double));
	double *first_p = malloc(timing->th_pls_nbr * sizeof(double));
	double exposure = timing->exposure;

	for (int i = 0; i < timing->th_pls_nbr; i++) {
		if (i >= timing->det_pulses) continue;	// Avoid an overflow, even if it should not occur here...
		unity_flux[i] = pulse[i].sum_flux / PULSE_WIDTH;	// in ADU/ms, assuming the pulse duration = 100ms
		first_p[i] = vmag[pulse[i].start_ind] / unity_flux[i];	// Lenght of the first pulse (ms)

		GDateTime *begin_frame = g_date_time_ref(seq->imgparam[pulse[i].start_ind_inseq - 1].date_obs);	// Timestamp at the beginning of the frame
		GDateTime *end_frame = g_date_time_add_seconds (begin_frame, exposure / 1000.0);	// (computed) Timestamp at the end of the frame
		GDateTime *end_frame2 = g_date_time_add_seconds (end_frame, -1.0 * first_p[i] / 1000.0);	// (Computed) Timestamp of the PPS
//		GDateTime *end_frame2 = g_date_time_add_seconds (begin_frame, -1.0 * first_p[i] / 1000.0);	// (Computed) Timestamp of the PPS
		GDateTime *pps_time = round_date (end_frame2);	// (Observed) Timestamp of the precise PPS time the pulse refers to

		if (verbose) {
			siril_log_color_message(_("unit_flux[%i] %lf, sum_flux(%i) = %lf\n"), "red", i, unity_flux[i], i, pulse[i].sum_flux);
			siril_log_color_message(_("Vmag(%i)= %lf, first_p[%i] %lf (ms)\n"), "salmon", pulse[i].start_ind_inseq, vmag[pulse[i].start_ind], i, first_p[i]);

			gchar *str = date_time_to_FITS_date(begin_frame);
			siril_log_color_message(_("begin_frame %s, exposure= %lf (ms)\n"), "salmon", str, 1000.0 * exposure);
			str = date_time_to_FITS_date(end_frame);
			siril_log_color_message(_("end_frame %s\n"), "salmon", str);
			str = date_time_to_FITS_date(end_frame2);
			siril_log_color_message(_("end_frame2 (supposed PPS) %s, with first_p[%i]= %lf (ms)\n"), "salmon", str, i, first_p[i]);
			str = date_time_to_FITS_date(pps_time);
			siril_log_color_message(_("Real pps_time %s\n"), "salmon", str);
			free(str);
		}

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
	timing->std_seq = gsl_stats_sd(tmp_dbl, 1, timing->det_pulses);
	gsl_stats_minmax (&lo_val, &hi_val, tmp_dbl, 1, timing->det_pulses);
	free(tmp_dbl);

	siril_log_color_message(_("Time Offset: %0.3lf(ms), Std= %0.3lf(ms)\n"), "green", timing->median_seq, timing->std_seq);
	if (verbose) siril_log_color_message(_("mini diff time %0.3lf(ms), maxi diff time %0.3lf(ms)\n"), "green", lo_val, hi_val);

	free(unity_flux);
	free(first_p);
	return 0;
}




/****************** making a light curve from sequence-stored data ****************/
/* It uses data stored in the sequence, in seq->photometry, which is
 * populated by successive calls to seqpsf on the opened sequence;
 */
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
//		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid){
		if (!seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid){
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
		free(date); free(vmag); free(orig_ind);
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
//				GDateTime *new_dt = g_date_time_add_seconds(tsi, seq->exposure * 0.5);
				GDateTime *new_dt = g_date_time_add_seconds(tsi, com.pref.phot_set.t_delayed ? seq->exposure * 0.5 + com.pref.phot_set.time_offset / 1000.0 : seq->exposure * 0.5);
				julian = date_time_to_Julian(new_dt);
				g_date_time_unref(new_dt);
			} else {		// this is the case for ser files. So the timestamp is the beginning of the frame
				GDateTime *new_dt = g_date_time_add_seconds(tsi, com.pref.phot_set.t_delayed ? com.pref.phot_set.time_offset / 1000.0 : 0.0);
//				julian = date_time_to_Julian(tsi);
				julian = date_time_to_Julian(new_dt);
				g_date_time_unref(new_dt);
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

	// Retrieve the fps, the expopsition time and expected number of pulses in the sequence
	// Depends on the sequnece type, from .ser or from fits sequence
	double fps, expos;	// fps, exposition time
	int exp_pls_nbr;	// expected pulses number

	if (com.seq.ser_file != NULL) {
		fps = com.seq.ser_file->fps;
		expos = 1000.0 / fps;
		exp_pls_nbr = (int)(expos * 1e-3 * (double)seq->number);	// Floors it to get the total expected number of pulses
	} else {
		if (!seq->imgparam[1].date_obs || !seq->imgparam[seq->number - 1].date_obs) {
			siril_log_color_message(_("Something went wrong in the fits sequence. Should use the .ser file as a sequence.\n"), "red");
			return 0;
		}
		GTimeSpan timespan = g_date_time_difference(seq->imgparam[seq->number - 1].date_obs, seq->imgparam[1].date_obs);	// Given in µs. Note: smthg strange happens with [0]. [1] is a workaround.
		expos = 1e-3 * timespan / (double)seq->number;	// Expressed in ms
		expos = expos * (1.0 + 2.0 / (double)seq->number);	// Correction because of the [0] issue. Still in ms
		exp_pls_nbr = (int)(expos * 1e-3 * (double)seq->number);	// Floors it to get the total expected number of pulses
		fps = 1000.0 / expos;
	}

	siril_log_color_message(_("Exposure= %0.3lf(ms), fps= %0.3lf\n"), "salmon", expos, fps);
	siril_log_color_message(_("Number of frames = %i, Number of valid images= %i\n"), "salmon", seq->number, j);
	
	struct occ_res *timing = NULL;	// Structure embedding the final results
	timing = malloc(sizeof(struct occ_res));
	timing->exposure = expos;
	timing->th_pls_nbr = exp_pls_nbr;
	timing->valid_images = j;

	// Retrieve the amplitude range
	gsl_stats_minmax (&timing->lo_val, &timing->hi_val, vmag, 1, timing->valid_images);
//	siril_log_color_message(_("min= %lf, max= %lf, nb_valid_images= %i, nbr pulse= %i, exposure(s)= %lf\n"), "red", timing->lo_val, timing->hi_val, timing->valid_images, timing->th_pls_nbr, timing->exposure);
	if ((int)timing->hi_val > com.pref.phot_set.maxval) {
		siril_log_color_message(_("The blinking star looks saturated.\n"), "red");
		return 0;
	}

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
	ret = time_offset_calc (pulse, vmag, seq, timing);
	com.pref.phot_set.time_offset = timing->median_seq;

	// Frees what needs to be
	free(date);
	free(vmag);
	free(orig_ind);

	return ret;
}

gpointer occultation_worker(gpointer arg) {
	int retval = 0;
	struct light_curve_args *args = (struct light_curve_args *)arg;
	siril_log_color_message(_("Entering the Occulation Time precedure.\n"), "green");
	framing_mode framing = REGISTERED_FRAME;

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

	siril_add_idle(end_occultation_worker, args);
	return GINT_TO_POINTER(retval);
}