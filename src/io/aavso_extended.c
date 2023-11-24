/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include <glib/gprintf.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "algos/PSF.h"
#include "core/siril_log.h"
#include "io/siril_plot.h"

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

	if (g_list_length(spl_data->plots) > 1) {
		siril_debug_print("AAVSO should not hold more than one data series, aborting\n");
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
			data[index++] = g_strdup_printf("%12.6lf",
					plots->plots[0]->data[i].x); // DATE
			for (int k = 0; k < 2; k++) // MAG, MAGERR
				data[index++] = g_strdup_printf("%8.6lf", plots->plots[k]->data[i].y);
			data[index++] = g_strdup(adata->filter); // FILT
			data[index++] = g_strdup("NO"); // TRANS
			data[index++] = g_strdup("STD"); // MTYPE
			data[index++] = g_strdup(adata->cname); // CNAME
			data[index++] = g_strdup_printf("%8.6lf", cplot->data[i].y); // CMAG
			data[index++] = g_strdup(adata->kname); // KNAME
			data[index++] = g_strdup_printf("%8.6lf", kplot->data[i].y); // KMAG
			data[index++] = g_strdup_printf("%.3lf", airmass->data[i].y); // AMASS
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
		g_string_free(header, TRUE);
		for (int i = 0; i < nbpoints * nbcols; i++) {
			g_free(data[i]);
		}
		g_free(data);
	}

	return retval;
}


static gboolean export_to_aavso_extended(siril_plot_data *data, aavso_dlg *aavso_ptr, const char *datfilename) {
	aavso_parameters header;
	aavso_data adata;
	char aavso_param[MAX_HEADER_LENGTH];
    const char *data_header = "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES";

    strcpy(header.type, "EXTENDED");
	g_strlcpy(header.obscode, aavso_ptr->obscode, 12);
	g_strlcpy(header.software, PACKAGE_STRING, 256);
	header.delim = ',';
	strcpy(header.date, "JD");
	g_strlcpy(header.obstype, aavso_ptr->obstype, 5);
	g_strlcpy(header.aavso_data_header, data_header, 512);

    generate_aavso_header(&header, aavso_param);

    g_strlcpy(adata.starid, aavso_ptr->starid, 31);
    g_strlcpy(adata.filter,aavso_ptr->filter, 6);
    g_strlcpy(adata.cname, aavso_ptr->cname, 21);
    g_strlcpy(adata.kname, aavso_ptr->kname, 21);
    g_strlcpy(adata.chart, aavso_ptr->chart, 21);
    g_strlcpy(adata.notes, aavso_ptr->notes, strlen(aavso_ptr->notes));

    siril_plot_save_aavso(data, datfilename, aavso_param, &adata);

	return TRUE;
}

int export_AAVSO(pldata *plot, sequence *seq, gchar *filename, void *ptr) {
	int i, j, nbImages = 0, c_idx, k_idx;
	aavso_dlg *aavso_ptr = NULL;
	double c_std = 0.0;
	double *vmag = NULL, *cstar = NULL, *kstar = NULL, *err = NULL, *x = NULL,
			*airmass = NULL;
	siril_plot_data *spl_data = NULL;
	if (!seq->photometry[0]) {
		siril_log_color_message(_("No photometry data found, error\n"), "red");
		return -1;
	}

	aavso_ptr = (aavso_dlg *)ptr;
	/* define c_std */
	c_std = aavso_ptr->c_std;
	c_idx = aavso_ptr->c_idx + 1;
	k_idx = aavso_ptr->k_idx + 1;

	/* get number of valid frames for each star */
	int ref_valid_count[MAX_SEQPSF] = { 0 };
	gboolean ref_valid[MAX_SEQPSF] = { FALSE };
	for (i = 0; i < plot->nb; i++) {
		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid)
			continue;
		++nbImages;
		for (int r = 1; r < MAX_SEQPSF && seq->photometry[r]; r++) {
			if (seq->photometry[r][i] && seq->photometry[r][i]->phot_is_valid)
				ref_valid_count[r]++;
		}
	}

	siril_debug_print("we have %d images with a valid photometry for the variable star\n", nbImages);
	if (nbImages < 1)
		return -1;

	int nb_ref_stars = 0;
	// select reference stars that are only available at least 3/4 of the time
	for (int r = 1; r < MAX_SEQPSF && seq->photometry[r]; r++) {
		ref_valid[r] = ref_valid_count[r] >= nbImages * 3 / 4;
		siril_debug_print("reference star %d has %d/%d valid measures, %s\n", r, ref_valid_count[r], nbImages, ref_valid[r] ? "including" : "discarding");
		if (ref_valid[r])
			nb_ref_stars++;
	}

	if (nb_ref_stars == 0) {
		siril_log_color_message(_("The reference stars are not good enough, probably out of the configured valid pixel range, cannot calibrate the light curve\n"), "red");
		return -1;
	}
	if (nb_ref_stars == 1)
		siril_log_color_message(_("Only one reference star was validated, this will not result in an accurate light curve. Try to add more reference stars or check the configured valid pixel range\n"), "salmon");
	else siril_log_message(_("Using %d stars to calibrate the light curve\n"), nb_ref_stars);

	vmag = calloc(nbImages, sizeof(double));
	err = calloc(nbImages, sizeof(double));
	x = calloc(nbImages, sizeof(double));
	cstar = calloc(nbImages, sizeof(double));
	kstar = calloc(nbImages, sizeof(double));
	airmass = calloc(nbImages, sizeof(double));
	if (!vmag || !err || !x || !cstar || !kstar || !airmass) {
		PRINT_ALLOC_ERR;
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
	for (i = 0, j = 0; i < plot->nb; i++) {
		if (!seq->imgparam[i].incl || !seq->photometry[0][i] || !seq->photometry[0][i]->phot_is_valid)
			continue;
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

		/* First data plotted are variable data, others are references
		 * Variable is done above, now we compute references */
		for (int r = 1; r < MAX_SEQPSF && seq->photometry[r]; r++) {
			if (ref_valid[r] && seq->photometry[r][i] && seq->photometry[r][i]->phot_is_valid) {
				/* inversion of Pogson's law to get back to the flux
				 * Flux = 10^(-0.4 * mag)
				 */
				if (r == c_idx) {
					cmag = pow(10, -0.4 * seq->photometry[r][i]->mag);
					cerr = seq->photometry[r][i]->s_mag;
				} else if (r == k_idx) {
					kmag = pow(10, -0.4 * seq->photometry[r][i]->mag);
				}
			}
		}

		/* Converting back to magnitude */

		if (cmag != 0.0) {
			cmag = -2.5 * log10(cmag);
			vmag[j] = target_mag - cmag + c_std;
			err[j] = fmin(9.999, sqrt(target_err * target_err + cerr * cerr));
		}

		if (kmag != 0.0) {
			kmag = -2.5 * log10(kmag);
		}

		cstar[j] = cmag;
		kstar[j] = kmag;
		j++;
	}
	int nb_valid_images = j;

	siril_log_message(_("Calibrated data for %d points of the light curve, %d excluded because of invalid calibration\n"), nb_valid_images, plot->nb - nb_valid_images);

	/*  data are computed, now plot the graph. */

	spl_data = malloc(sizeof(siril_plot_data));
	init_siril_plot_data(spl_data);

	spl_data->revertY = TRUE;
	siril_plot_set_xlabel(spl_data, "Julian date");
	siril_plot_add_xydata(spl_data, "V-C", nb_valid_images, x, vmag, err, NULL);

	siril_plot_add_xydata(spl_data, "cmag", nb_valid_images, x, cstar, NULL, NULL);
	siril_plot_add_xydata(spl_data, "kmag", nb_valid_images, x, kstar, NULL, NULL);
	siril_plot_add_xydata(spl_data, "airmass", nb_valid_images, x, airmass, NULL, NULL);

	int ret = export_to_aavso_extended(spl_data, aavso_ptr, filename);

	free_siril_plot_data(spl_data);

	free(vmag);
	free(err);
	free(x);
	free(airmass);
	free(cstar);
	free(kstar);
	return ret;
}
