#ifndef _CATALOGUES_H
#define _CATALOGUES_H

#include <glib.h>
#include "core/siril.h"
#include "core/siril_world_cs.h"
#include "algos/photometry.h"
#include "io/siril_catalogues.h"

int get_stars_from_local_catalogues(double ra, double dec, double radius, fits *fit, float max_mag, pcc_star **stars, int *nb_stars, gboolean photo_only);

void initialize_local_catalogues_paths();
gboolean local_catalogues_available();

gchar *get_and_project_local_catalog(SirilWorldCS *catalog_center, double radius, double max_mag, gboolean for_photometry);
siril_catalogue *siril_catalog_get_raw_stars_from_local_catalogues(SirilWorldCS *catalog_center, double radius, double max_mag,
		gboolean for_photometry);

#endif
