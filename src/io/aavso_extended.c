/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include <glib/gprintf.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "algos/PSF.h"
#include "core/siril_log.h"
#include "io/siril_plot.h"
#include "gui/siril_plot.h"

#include "aavso_extended.h"

#define MAX_HEADER_LENGTH 1024 // Maximum length for the data string
#define _NA_ "na"

// We want to implement this format: https://www.aavso.org/aavso-extended-file-format

//   #TYPE=EXTENDED
//   #OBSCODE=TST01
//   #SOFTWARE=IRAF 12.4
//   #DELIM=,
//   #DATE=JD
//   #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
//   SS CYG,2450702.1234,11.135,0.003,V,NO,STD,105,na,na,na,na,na,X16382L,This is a test


static void generate_aavso_header(const aavso_parameters *data, char *result) {
    snprintf(result, MAX_HEADER_LENGTH,
             "#TYPE=%s\n#OBSCODE=%s\n#SOFTWARE=%s\n#DELIM=%c\n#DATE=%s\n#OBSTYPE=%s\n%s",
             data->type, data->obscode, data->software, data->delim, data->date, data->obstype, data->aavso_data_header);
}

// save the data to a txt file
static gboolean siril_plot_save_aavso(siril_plot_data *spl_data,
		const char *datfilename, const char *head_data, aavso_data *adata) {
	GString *header = g_string_new(head_data);
	gboolean retval = TRUE;

	if (g_list_length(spl_data->plots) != 1 && g_list_length(spl_data->plot) != 3) {
		siril_debug_print("AAVSO data are not as expected, aborting\n");
		g_string_free(header, TRUE);
		return FALSE;
	}

	int nbcols = 15; // aavso has 15 columns

	if (spl_data->plots != NULL && spl_data->plot) {
		splxyerrdata *plots = (splxyerrdata*) spl_data->plots->data;
		int nbpoints = plots->nb;

		GList *list = spl_data->plot;

		splxydata *cplot = (splxydata*) list->data;
		splxydata *kplot = (splxydata*) list->next->data;
		splxydata *airmass = (splxydata*) list->next->next->data;

		if (nbpoints != cplot->nb || nbpoints != kplot->nb
				|| nbpoints != airmass->nb) {
			g_string_free(header, TRUE);
			return FALSE;
		}

		gchar **data = g_new(gchar*, nbpoints * nbcols);
		if (!data) {
			PRINT_ALLOC_ERR;
			g_string_free(header, TRUE);
			return FALSE;
		}

		int index = 0;
		for (int i = 0; i < nbpoints; i++) {
			data[index++] = g_strdup(adata->starid); // NAME
			data[index++] = g_strdup_printf("%12.6lf", plots->plots[0]->data[i].x); // DATE
			for (int k = 0; k < 2; k++) // MAG, MAGERR
				data[index++] = g_strdup_printf("%.3lf", plots->plots[k]->data[i].y);
			data[index++] = g_strdup(adata->filter); // FILT
			data[index++] = g_strdup("NO"); // TRANS
			data[index++] = g_strdup("STD"); // MTYPE
			data[index++] = g_strdup(adata->cname); // CNAME
			data[index++] = g_strdup_printf("%.3lf", cplot->data[i].y); // CMAG
			data[index++] = g_strdup(adata->kname); // KNAME
			data[index++] = kplot->data[i].y == 0 ? g_strdup(_NA_) : g_strdup_printf("%.3lf", kplot->data[i].y); // KMAG
			data[index++] = airmass->data[i].y == 0 ? g_strdup(_NA_) :  g_strdup_printf("%.3lf", airmass->data[i].y); // AMASS
			data[index++] = g_strdup(_NA_); // GROUP
			data[index++] = g_strdup(adata->chart); // CHART
			data[index++] = g_strdup(adata->notes); // NOTES
		}

		FILE *fileout = g_fopen(datfilename, "w");
		if (fileout == NULL) {
			siril_log_message(_("Could not create %s, aborting\n"), datfilename);
			retval = FALSE;
		} else {
			fprintf(fileout, "%s", header->str);
			index = 0;
			for (int r = 0; r < nbpoints; r++) {
				fprintf(fileout, "\n%s", data[index++]); // print newline and x
				for (int c = 1; c < nbcols; c++) {
					if (index < nbpoints * nbcols) {
						fprintf(fileout, ",%s", data[index++]);
					}
				}
			}

			fclose(fileout);
			siril_log_message(_("%s has been saved.\n"), datfilename);
		}

		// Free memory
		for (int i = 0; i < nbpoints * nbcols; i++) {
			g_free(data[i]);
		}
		g_free(data);
	}
	g_string_free(header, TRUE);
	return retval;
}


static gboolean export_to_aavso_extended(siril_plot_data *data, aavso_dlg *aavso_ptr, const char *datfilename) {
	aavso_parameters header;
	aavso_data adata;
	char aavso_param[MAX_HEADER_LENGTH];
	const char *data_header = "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES";

	strcpy(header.type, "EXTENDED");
	gboolean truncated = FALSE;
	if (g_strlcpy(header.obscode, aavso_ptr->obscode, 12) >= 12) truncated = TRUE;
	if (g_strlcpy(header.software, PACKAGE_STRING, 256) >= 256) truncated = TRUE;
	header.delim = ',';
	strcpy(header.date, "JD");
	if (g_strlcpy(header.obstype, aavso_ptr->obstype, 5) >= 5) truncated = TRUE;
	if (g_strlcpy(header.aavso_data_header, data_header, 512) >= 512) truncated = TRUE;

	if (truncated) {
		siril_log_color_message(_("Warning: at least one AAVSO header element was truncated.\n"), "salmon");
		truncated = FALSE;
	}
	generate_aavso_header(&header, aavso_param);

	if (g_strlcpy(adata.starid, aavso_ptr->starid, 31) >= 31) truncated = TRUE;
	if (g_strlcpy(adata.filter, aavso_ptr->filter, 6) >= 6) truncated = TRUE;
	if (g_strlcpy(adata.cname, aavso_ptr->cname, 21) >= 21) truncated = TRUE;
	if (g_strlcpy(adata.kname, aavso_ptr->kname, 21) >= 21) truncated = TRUE;
	if (g_strlcpy(adata.chart, aavso_ptr->chart, 21) >= 21) truncated = TRUE;
	gsize len = strlen(aavso_ptr->notes) + 1;
	if (g_strlcpy(adata.notes, aavso_ptr->notes, len) >= len) truncated = TRUE;

	if (truncated) {
		siril_log_color_message(_("Warning: at least one AAVSO data element was truncated.\n"), "salmon");
		truncated = FALSE;
	}
	siril_plot_save_aavso(data, datfilename, aavso_param, &adata);

	return TRUE;
}

int export_AAVSO(pldata *plot, sequence *seq, gchar *filename, gchar **error, void *ptr) {
	int i, j, nbImages = 0, c_idx, k_idx;
	aavso_dlg *aavso_ptr = NULL;
	double c_std = 0.0;
	double *vmag = NULL, *cstar = NULL, *kstar = NULL, *err = NULL, *x = NULL,
			*airmass = NULL;
	siril_plot_data *spl_data = NULL;
	if (!seq->photometry[0]) {
		gchar *msg = siril_log_color_message(_("No photometry data found, error\n"), "red");
		*error = g_strdup(msg);
		return -1;
	}

	aavso_ptr = (aavso_dlg *)ptr;
	/* define c_std */
	c_std = aavso_ptr->c_std;
	c_idx = aavso_ptr->c_idx + 1;
	k_idx = aavso_ptr->k_idx + 1;

	/* get number of valid frames for c_star and k_star */
	int ref_valid_count[MAX_SEQPSF] = { 0 };
	gboolean ref_valid[MAX_SEQPSF] = { FALSE };
	for (i = 0; i < seq->number; i++) {
		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid)
			continue;
		++nbImages;
		if (seq->photometry[c_idx][i] && seq->photometry[c_idx][i]->phot_is_valid) {
			ref_valid_count[c_idx]++;
		}
		if (seq->photometry[k_idx][i] && seq->photometry[k_idx][i]->phot_is_valid) {
			ref_valid_count[k_idx]++;
		}

	}

	siril_debug_print("we have %d images with a valid photometry for the variable star\n", nbImages);
	if (nbImages < 1) {
		gchar *msg = siril_log_color_message(_("Number of images is not enough\n"), "red");
		*error = g_strdup(msg);
		return -1;
	}
	int nb_ref_stars = 0;

	// select reference c_star and k_star that are only available at least 3/4 of the time
	ref_valid[c_idx] = ref_valid_count[c_idx] >= nbImages * 3 / 4;
	siril_debug_print("reference star %d has %d/%d valid measures, %s\n", c_idx, ref_valid_count[c_idx], nbImages, ref_valid[c_idx] ? "including" : "discarding");
	if (ref_valid[c_idx])
	    nb_ref_stars++;
	ref_valid[k_idx] = ref_valid_count[k_idx] >= nbImages * 3 / 4;
	siril_debug_print("reference star %d has %d/%d valid measures, %s\n", k_idx, ref_valid_count[k_idx], nbImages, ref_valid[k_idx] ? "including" : "discarding");
	if (ref_valid[k_idx])
	    nb_ref_stars++;

	if (nb_ref_stars < 2) { // we want both c_idx and k_idx valid
		gchar *msg = siril_log_color_message(_("The reference stars are not good enough, probably out of the configured valid pixel range, cannot calibrate the data\n"), "red");
		*error = g_strdup(msg);
		return -1;
	}
	else siril_log_message(_("Using %d stars to calibrate the data\n"), nb_ref_stars);

	vmag = calloc(nbImages, sizeof(double));
	err = calloc(nbImages, sizeof(double));
	x = calloc(nbImages, sizeof(double));
	cstar = calloc(nbImages, sizeof(double));
	kstar = calloc(nbImages, sizeof(double));
	airmass = calloc(nbImages, sizeof(double));
	if (!vmag || !err || !x || !cstar || !kstar || !airmass) {
		PRINT_ALLOC_ERR;
		gchar *msg = siril_log_color_message(_("Out of memory\n"), "red");
		*error = g_strdup(msg);
		free(vmag);
		free(err);
		free(x);
		free(airmass);
		free(cstar);
		free(kstar);
		return -1;
	}

	double min_date = DBL_MAX;
	// i is index in dataset, j is index in output
	for (i = 0, j = 0; i < seq->number; i++) {
		if (!seq->imgparam[i].incl)
			continue;
		// if target or compstar are invalid, we won't log the point
		if (!seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid ||
			!seq->photometry[c_idx][i] || !seq->photometry[c_idx][i]->phot_is_valid) {
				siril_debug_print("Image %d discarded because of invalid calibration data\n", i + 1);
				continue;
		}
		double cmag = 0.0, cerr = 0.0, kmag = 0.0;

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
			x[j] = julian;
			if (julian < min_date)
				min_date = julian;
		} else {
			x[j] = (double) i + 1; // should not happen.
		}


		// Y value: the magnitude and error and their calibration
		double target_mag = seq->photometry[0][i]->mag;
		double target_err = seq->photometry[0][i]->s_mag;
		airmass[j] = seq->imgparam[i].airmass;

		/* First data plotted are variable data, c_idx is comparison star
		 * and k_idx stands for the check one */
		cmag = pow(10, -0.4 * seq->photometry[c_idx][i]->mag);
		cerr = seq->photometry[c_idx][i]->s_mag;

		if (seq->photometry[k_idx][i] && seq->photometry[k_idx][i]->phot_is_valid) {
			/* if non-valid, kmag is equal to 0. It will be "na" in the file */
				kmag = pow(10, -0.4 * seq->photometry[k_idx][i]->mag);
		}

		/* Converting back to magnitude */
		cmag = -2.5 * log10(cmag);
		vmag[j] = target_mag - cmag + c_std;
		err[j] = fmin(9.999, sqrt(target_err * target_err + cerr * cerr));

		if (kmag != 0.0) {
			kmag = -2.5 * log10(kmag);
		}

		cstar[j] = cmag;
		kstar[j] = kmag;
		j++;
	}
	int nb_valid_images = j;

	siril_log_message(_("Calibrated data for %d points of the output data, %d excluded because of invalid calibration\n"), nb_valid_images, plot->nb - nb_valid_images);

	/*  data are computed, now plot the graph. */

	spl_data = init_siril_plot_data();
	int ret;
	if (!spl_data) {
		PRINT_ALLOC_ERR;
		ret = -1;
	} else {
		spl_data->revertY = TRUE;
		siril_plot_set_xlabel(spl_data, "Julian date");
		siril_plot_add_xydata(spl_data, "V-C", nb_valid_images, x, vmag, err, NULL);

		siril_plot_add_xydata(spl_data, "cmag", nb_valid_images, x, cstar, NULL, NULL);
		siril_plot_add_xydata(spl_data, "kmag", nb_valid_images, x, kstar, NULL, NULL);
		siril_plot_add_xydata(spl_data, "airmass", nb_valid_images, x, airmass, NULL, NULL);

		// now we sort to have all dates ascending
		siril_plot_sort_x(spl_data);

		ret = export_to_aavso_extended(spl_data, aavso_ptr, filename);

//		siril_plot_set_title(spl_data, "AAVSO data");
//		create_new_siril_plot_window(spl_data);

		free_siril_plot_data(spl_data);
	}

	free(vmag);
	free(err);
	free(x);
	free(airmass);
	free(cstar);
	free(kstar);
	return ret;
}
