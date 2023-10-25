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
#ifndef SRC_IO_ANNOTATION_CATALOGUES_H_
#define SRC_IO_ANNOTATION_CATALOGUES_H_

#include "core/siril_world_cs.h"
#include "io/siril_catalogues.h"

typedef struct _CatalogObjects CatalogObjects;

#define NB_AN_CATALOGUES 8

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

typedef struct annotations_catalogue {
	siril_catalogue *cat;
	gboolean show;
} annotations_catalogue_t;

GSList *find_objects_in_field(fits *fit);
void cleanup_annotation_catalogues();
gchar *get_annotation_catalog_filename(annotations_cat cat_index, gboolean for_reading);
void add_item_in_catalogue(cat_item *item, annotations_cat cat_index, gboolean check_duplicates);
cat_item *search_in_annotations_by_name(const char *input, object_catalog *cattype);
cat_item *search_in_solar_annotations(sky_object_query_args *args);
const char *cat_index_to_name(annotations_cat index);
void refresh_annotation_visibility();
void refresh_annotation_to_temp();

gchar *get_catalogue_object_code(const CatalogObjects *object);
gchar *get_catalogue_object_code_pretty(CatalogObjects *object);
annotations_cat get_catalogue_object_cat(const CatalogObjects *object);
gdouble get_catalogue_object_x(const CatalogObjects *object);
gdouble get_catalogue_object_y(const CatalogObjects *object);
gdouble get_catalogue_object_radius(const CatalogObjects *object);

void refresh_found_objects();
void purge_user_catalogue(object_catalog cattype);
gboolean is_inside(fits *fit, double ra, double dec);
int load_csv_targets_to_temp(const gchar *filename);
int load_siril_cat_to_temp(siril_catalogue *siril_cat);

#endif /* SRC_IO_ANNOTATION_CATALOGUES_H_ */
