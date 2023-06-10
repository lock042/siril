#ifndef _REMOTE_CATALOGUES_H
#define _REMOTE_CATALOGUES_H

#include <glib.h>
#include <gio/gio.h>
#include "core/siril_world_cs.h"
#include "algos/PSF.h"

//#define SIMBADPHOTO "http://simbad.cds.unistra.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=TSV&query=SELECT B, V, R ,I , J from allfluxes JOIN ident USING(oidref) WHERE id ='"
#define SKYBOT "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?"
#define VIZIER_QUERY "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source="


#define AAVSO_QUERY "https://app.aavso.org/vsp/api/chart/?"

typedef enum {
	CAT_TYCHO2,
	CAT_NOMAD,
	CAT_GAIADR3,
	CAT_PPMXL,
	CAT_BRIGHT_STARS,
	CAT_APASS,
	CAT_AAVSO,
	CAT_AUTO = 98,
	CAT_LOCAL = 99,		// siril local (KStars Tycho-2 and NOMAD)
	CAT_ASNET = 100,	// solve-field local (astrometry.net)
} online_catalog;	// TODO: rename?

const char *catalog_to_str(online_catalog cat);

GFile *download_catalog(online_catalog onlineCatalog, SirilWorldCS *catalog_center, double radius, double mag);
gchar *get_catalog_url(SirilWorldCS *center, double mag_limit, double dfov, int type);

gchar *fetch_url(const gchar *url);
void free_fetch_result(gchar *result);

int read_projected_catalog(GInputStream *stream, psf_star **cstars, online_catalog cat);

gpointer search_in_online_conesearch(gpointer p);

// temp
//struct compstars_arg;
int load_catalog(GFile *catalog_file, gboolean phot, psf_star **ret_stars, int *ret_nb_stars);
//int read_photo_catalog_buffer(const char *buffer, struct compstars_arg *args);
//int read_photo_aavso_buffer(const char *buffer, struct compstars_arg *args);

#endif
