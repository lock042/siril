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
#ifndef SRC_ALGOS_ANNOTATE_H_
#define SRC_ALGOS_ANNOTATE_H_

#include "core/siril_world_cs.h"

typedef struct _CatalogObjects CatalogObjects;

#define USER_DSO_CAT_INDEX 6
#define USER_SSO_CAT_INDEX 7
#define USER_TEMP_CAT_INDEX 8

GSList *find_objects(fits *fit);
void add_object_in_catalogue(gchar *code, SirilWorldCS *wcs, gboolean is_solar_system, gboolean is_in_field);
gchar *get_catalogue_object_code(CatalogObjects *object);
guint get_catalogue_object_cat(CatalogObjects *object);
gchar *get_catalogue_object_name(CatalogObjects *object);
gdouble get_catalogue_object_ra(CatalogObjects *object);
gchar *retrieve_site_coord (fits *fit);
gdouble get_catalogue_object_dec(CatalogObjects *object);
gdouble get_catalogue_object_radius(CatalogObjects *object);
void refresh_found_objects();
void free_catalogue_object(CatalogObjects *object);
void purge_temp_user_catalogue();
gboolean is_inside(fits *fit, double ra, double dec);
void load_csv_targets_to_temp(const gchar *filename);

#endif /* SRC_ALGOS_ANNOTATE_H_ */
