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

#ifndef _SEARCH_OBJECTS_H
#define _SEARCH_OBJECTS_H

#include <glib.h>
#include "core/siril.h"
#include "core/siril_world_cs.h"
#include "algos/annotate.h"
#include "algos/PSF.h"

#define CDSSESAME "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame"
#define VIZIERSESAME "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame"
#define SIMBADSESAME "http://simbad.cds.unistra.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=TSV&query=SELECT basic.OID, ra, dec, main_id FROM basic JOIN ident ON ident.oidref = oid WHERE id ='"

#define SIMBAD "http://simbad.u-strasbg.fr/simbad/sim-id?output.format=ASCII&Ident="
#define EPHEMCC "https://ssp.imcce.fr/webservices/miriade/api/ephemcc.php?"

typedef enum {
	QUERY_SERVER_CDS,
	QUERY_SERVER_VIZIER,
	QUERY_SERVER_SIMBAD,
	QUERY_SERVER_EPHEMCC,
	QUERY_SERVER_SIMBAD_PHOTO,
	QUERY_SERVER_SKYBOT, // In case of adding other items, leave this one at the end of the list
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

int parse_catalog_buffer(const gchar *buffer, psf_star **result);
int parse_conesearch_buffer(const gchar *buffer, double lim_mag);
int cached_object_lookup(const gchar *name, psf_star **opt_result);
gchar *search_in_online_catalogs(const gchar *object, query_server server);

void add_plated_from_annotations(const CatalogObjects *obj);
void free_Platedobject();
gboolean has_nonzero_coords();
int parse_resolver_buffer(const gchar *buffer, struct sky_object *obj);

#endif
