/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_IO_ANNOTATION_CATALOGUES_H_
#define SRC_IO_ANNOTATION_CATALOGUES_H_

#include "core/siril_world_cs.h"
#include "io/siril_catalogues.h"

typedef struct _CatalogObjects CatalogObjects;

typedef struct annotations_catalogue {
	siril_catalogue *cat;
	gboolean show;
} annotations_catalogue_t;

GSList *find_objects_in_field(fits *fit);
gchar *get_annotation_catalog_filename(siril_cat_index cat_index, gboolean for_reading);
void add_item_in_catalogue(cat_item *item, siril_cat_index cat_index, gboolean check_duplicates);
cat_item *search_in_annotations_by_name(const char *input, siril_cat_index *cat_index);
cat_item *search_in_solar_annotations(sky_object_query_args *args);
void set_annotation_visibility(siril_cat_index cat_index, gboolean visible);
void refresh_annotation_visibility();
void refresh_annotation_to_temp();
void cleanup_annotation_catalogues(gboolean purge_temp);
void refresh_annotations(gboolean purge_temp);

gchar *get_catalogue_object_code(const CatalogObjects *object);
gchar *get_catalogue_object_code_pretty(CatalogObjects *object);
siril_cat_index get_catalogue_object_cat(const CatalogObjects *object);
gdouble get_catalogue_object_x(const CatalogObjects *object);
gdouble get_catalogue_object_y(const CatalogObjects *object);
gdouble get_catalogue_object_x1(const CatalogObjects *object);
gdouble get_catalogue_object_y1(const CatalogObjects *object);
gdouble get_catalogue_object_radius(const CatalogObjects *object);

void refresh_found_objects();
void purge_user_catalogue(siril_cat_index cat_index);
gboolean is_inside(fits *fit, double ra, double dec);
int load_csv_targets_to_temp(const gchar *filename);
int load_siril_cat_to_temp(siril_catalogue *siril_cat);

#endif /* SRC_IO_ANNOTATION_CATALOGUES_H_ */
