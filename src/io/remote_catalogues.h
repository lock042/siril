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
#ifndef _REMOTE_CATALOGUES_H
#define _REMOTE_CATALOGUES_H

#include "io/siril_catalogues.h"

// new queries
#define VIZIER_TAP_QUERY "http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=SELECT+"
#define EXOPLANETARCHIVE_TAP_QUERY "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=select+"
#define SIMBAD_TAP_QUERY "https://simbad.u-strasbg.fr/simbad/sim-tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=SELECT+"
#define IMCCE_QUERY "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?&-mime=text&-output=basic&-filter=0&-objFilter=111&-refsys=EQJ2000&-from=Siril"
#define AAVSOCHART_QUERY "https://app.aavso.org/vsp/api/chart/?format=json"

// only the first 9 columns are valid TAP queries
// fields after this are used in other catalogues
#define MAX_TAP_QUERY_COLUMNS 9

// The structure used to declare the columns to be queried from the tables
// for TAP queries only!
typedef struct {
	gchar *catcode;
	gchar *tap_columns[MAX_TAP_QUERY_COLUMNS];
	gchar *tap_server;
} cat_tap_query_fields;

gchar *download_catalog(siril_catalogue *siril_cat);

gchar *fetch_url(const gchar *url);
void free_fetch_result(gchar *result);

gpointer search_in_online_conesearch(gpointer p);
gpointer catsearch_worker(gpointer p);

int siril_catalog_get_stars_from_online_catalogues(siril_catalogue *siril_cat);

#endif
