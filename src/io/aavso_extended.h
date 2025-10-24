/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at siril_free.fr)
 * Copyright (C) 2012-2020 team siril_free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is siril_free software: you can redistribute it and/or modify
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
#ifndef SRC_IO_AAVSO_EXTENDED_H_
#define SRC_IO_AAVSO_EXTENDED_H_

#include "io/siril_plot.h"

typedef struct aavso_parameters_struct {
	char type[9]; // TYPE: Should always say Extended for this format.
	char obscode[12]; // The official AAVSO Observer Code for the observer which was previously assigned by the AAVSO. // FIXME: what is the size?
	char software[256];
	char delim; // ','
	char date[3]; // JD
	char obstype[5]; // CCD, DSLR, PEP
	char aavso_data_header[512]; // header of aavso_data
} aavso_parameters;

typedef struct aavso_data_struct {
	char starid[31];       // STARID: Limit: 30 characters
	char date[17];         // DATE: Limit: 16 characters
	char magnitude[9];     // MAGNITUDE: Limit: 8 characters
	char magerr[7];        // MAGERR: Limit: 6 characters
	char filter[6];        // FILTER: list is in filter_list
	char trans[4];         // TRANS: Limit: 3 characters
	char mtype[4];         // MTYPE: Limit: 3 characters
	char cname[21];        // CNAME: Limit: 20 characters
	char cmag[9];          // CMAG: Limit: 8 characters
	char kname[21];        // KNAME: Limit: 20 characters
	char kmag[9];          // KMAG: Limit: 8 characters
	char airmass[8];       // AIRMASS: Limit: 7 characters
	char group[6];         // GROUP: Limit: 5 characters
	char chart[21];        // CHART: Limit: 20 characters
	char notes[4000];      // NOTES: Limit: Several thousand characters
} aavso_data;

typedef struct aavso_dlg_struct {
	const gchar *obscode;
	const gchar *obstype;
	const gchar *filter;
	const gchar *starid;
	const gchar *cname;
	int c_idx;
	const gchar *kname;
	int k_idx;
	double c_std;
	const gchar *chart;
	const gchar *notes;
} aavso_dlg;

int export_AAVSO(pldata *plot, sequence *seq, gchar *filename, gchar ** error, void *ptr);

#endif /* SRC_IO_AAVSO_EXTENDED_H_ */
