#ifndef _CATALOGUES_H
#define _CATALOGUES_H

#include <glib.h>
#include "core/siril.h"
#include "core/siril_world_cs.h"
#include "algos/photometry.h"
#include "io/siril_catalogues.h"

void initialize_local_catalogues_paths();
gboolean local_catalogues_available();

int siril_catalog_get_stars_from_local_catalogues(siril_catalogue *siril_cat);
gpointer write_trixels(gpointer p);
gpointer list_trixels(gpointer p);

#endif
