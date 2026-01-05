/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#ifndef _SEARCH_OBJECTS_H
#define _SEARCH_OBJECTS_H

#include <glib.h>
#include "core/siril_world_cs.h"
#include "io/annotation_catalogues.h"
#include "algos/PSF.h"

#define CDSSESAME "http://cds.unistra.fr/cgi-bin/nph-sesame"
#define VIZIERSESAME "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame"
#define SIMBADSESAME "http://simbad.cds.unistra.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=TSV&query=SELECT basic.OID, ra, dec, main_id FROM basic JOIN ident ON ident.oidref = oid WHERE id ='"

#define SIMBAD "http://simbad.cds.unistra.fr/simbad/sim-id?output.format=ASCII&Ident="
#define EPHEMCC "https://ssp.imcce.fr/webservices/miriade/api/ephemcc.php?-tcoor=5&-mime=text/csv&-output=--jd&-from=Siril"

typedef enum {
	QUERY_SERVER_UNSET= -1,
	QUERY_SERVER_CDS,
	QUERY_SERVER_VIZIER,
	QUERY_SERVER_SIMBAD,
	QUERY_SERVER_EPHEMCC,
	QUERY_SERVER_SIMBAD_PHOTO
} query_server;

typedef enum {
	RESOLVER_UNSET = -1,
	RESOLVER_NED = 0,
	RESOLVER_SIMBAD,
	RESOLVER_VIZIER,
	RESOLVER_LOCAL,
	RESOLVER_NUMBER,
} resolver_t;

struct sky_object {
	gchar *name;
	double radius;
	int maxRecords;
	SirilWorldCS *world_cs;
	point imageCenter;
	gboolean south;
};

int parse_catalog_buffer(const gchar *buffer, sky_object_query_args *args);
int cached_object_lookup(sky_object_query_args *args);
char *search_in_online_catalogs(sky_object_query_args *args);
void search_object(GtkEntry *entry);

void add_plated_from_annotations(const cat_item *obj);
void free_Platedobject();
gboolean has_nonzero_coords();
int parse_resolver_buffer(const gchar *buffer, struct sky_object *obj);

gpointer catsearch_worker(gpointer p);

#endif
