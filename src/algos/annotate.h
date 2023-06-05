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

typedef enum {
	ANCAT_NONE = -1,
	ANCAT_MESSIER = 0,
	ANCAT_NGC = 1,
	ANCAT_IC = 2,
	ANCAT_LDN = 3,
	ANCAT_SH2 = 4,
	ANCAT_STARS = 5,
	USER_DSO_CAT_INDEX = 6,
	USER_SSO_CAT_INDEX = 7,
	USER_TEMP_CAT_INDEX = 8
} annotations_cat;

GSList *find_objects_in_field(fits *fit);
void add_object_in_catalogue(gchar *code, SirilWorldCS *wcs, gboolean check_duplicates, annotations_cat cat);
const CatalogObjects *search_in_annotations_by_name(const char *target);
const char *cat_index_to_name(annotations_cat index);

gchar *get_catalogue_object_code(const CatalogObjects *object);
gchar *get_catalogue_object_code_pretty(CatalogObjects *object);
annotations_cat get_catalogue_object_cat(const CatalogObjects *object);
gchar *get_catalogue_object_name(const CatalogObjects *object);
gdouble get_catalogue_object_ra(const CatalogObjects *object);
gdouble get_catalogue_object_dec(const CatalogObjects *object);
gdouble get_catalogue_object_radius(const CatalogObjects *object);

void refresh_found_objects();
void free_catalogue_object(CatalogObjects *object);
void purge_temp_user_catalogue();
gboolean is_inside(fits *fit, double ra, double dec);
gboolean is_inside2(fits *fit, double ra, double dec, double *x, double *y);
int load_csv_targets_to_temp(const gchar *filename);

#endif /* SRC_ALGOS_ANNOTATE_H_ */
