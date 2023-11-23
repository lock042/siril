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

#include "core/siril.h"
#include "core/proto.h"
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
static gboolean siril_plot_save_aavso(siril_plot_data *spl_data, const char *datfilename, const char *head_data, aavso_data *adata) {
    GString *header = NULL;
    FILE *fileout = NULL;
    gboolean retval = TRUE;

    if (g_list_length(spl_data->plots) > 1) {
        siril_debug_print("Light curve should not hold more than one data series, aborting\n");
        return FALSE;
    }

    int nbpoints = 0, nbcols = 15; // aavso has 15 columns

    header = g_string_new(head_data);

	// Check if the data series is available
	if (spl_data->plots != NULL && spl_data->plot) {
		splxyerrdata *plots = (splxyerrdata*) spl_data->plots->data;
		nbpoints = plots->nb;

		GList *list = spl_data->plot;

		/* Get cstar data */
		splxydata *cplot = (splxydata*) list->data;
		if (nbpoints != cplot->nb) {
			retval = FALSE;
			goto clean_and_exit;
		}
		/* move to kstar data */
		list = list->next;
		splxydata *kplot = (splxydata*) list->data;
		if (nbpoints != kplot->nb) {
			retval = FALSE;
			goto clean_and_exit;
		}

		/* move to airmass data */
		list = list->next;
		splxydata *airmass = (splxydata*) list->data;
		if (nbpoints != airmass->nb) {
			retval = FALSE;
			goto clean_and_exit;
		}

		// Allocate memory for data array
		gchar **data = g_new(gchar*, nbpoints * nbcols);
		if (!data) {
			PRINT_ALLOC_ERR;
			retval = FALSE;
            goto clean_and_exit;
        }

        // Fill the data array
        int index = 0;
        for (int i = 0; i < nbpoints; i++) {
            data[index++] = g_strdup(adata->starid); // NAME
            data[index++] = g_strdup_printf("%12.6lf", plots->plots[0]->data[i].x + plots->plots[0]->x_offset); // DATE
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
            data[index++] = g_strdup(_NA_); // CHART
            data[index++] = g_strdup(_NA_); // NOTE
        }

        // Write data to the file
        fileout = g_fopen(datfilename, "w");
        if (fileout == NULL) {
            siril_log_message(_("Could not create %s, aborting\n"), datfilename);
            retval = FALSE;
            goto clean_and_exit;
        }

        fprintf(fileout, "%s", header->str);
        index = 0;
        for (int r = 0; r < nbpoints; r++) {
            fprintf(fileout, "\n%s", data[index++]); // print newline and x
            for (int c = 1; c < nbcols; c++)
                fprintf(fileout, ",%s", data[index++]);
        }

        fclose(fileout);
        siril_log_message(_("%s has been saved.\n"), datfilename);

clean_and_exit:
        // Free memory
        g_string_free(header, TRUE);
        if (data != NULL) {
            for (int i = 0; i < nbpoints * nbcols; i++) {
                g_free(data[i]);
            }
            g_free(data);
        }
    }

    return retval;
}

gboolean export_to_aavso_extended(siril_plot_data *data, aavso_dlg *aavso_ptr, const char *datfilename) {
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

    siril_plot_save_aavso(data, datfilename, aavso_param, &adata);

	return TRUE;
}
