/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include <math.h>

#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "algos/siril_wcs.h"

#include "annotate.h"

static const gchar *cat[] = { "messier.txt"/*, "ngc.txt", "ic.txt"*/ };

struct _CatalogObjects {
	gchar *code;
	gdouble ra;
	gdouble dec;
	gdouble radius;
	gchar *name;
};

static CatalogObjects *new_catalog_object(gchar *code, gdouble ra, gdouble dec, gdouble radius, gchar *name) {
	CatalogObjects *object = g_new(CatalogObjects, 1);
	object->code = g_strdup(code);
	object->ra = ra;
	object->dec = dec;
	object->radius = radius;
	if (name) {
		object->name = g_strdup(name);
	}
	return object;
}

/* A utility function to calculate area of
 triangle formed by (x1, y1), (x2, y2) and
 (x3, y3) */
static float area(float x1, float y1, float x2, float y2, float x3, float y3) {
	return fabs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0);
}
/* A function to check whether point P(x, y)
 lies inside the rectangle formed by A(x1, y1),
 B(x2, y2), C(x3, y3) and D(x4, y4) */
static gboolean is_in_picture(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4,
		float x, float y) {
	/* Calculate area of rectangle ABCD */
	float A = area(x1, y1, x2, y2, x3, y3) + area(x1, y1, x4, y4, x3, y3);

	/* Calculate area of triangle PAB */
	float A1 = area(x, y, x1, y1, x2, y2);

	/* Calculate area of triangle PBC */
	float A2 = area(x, y, x2, y2, x3, y3);

	/* Calculate area of triangle PCD */
	float A3 = area(x, y, x3, y3, x4, y4);

	/* Calculate area of triangle PAD */
	float A4 = area(x, y, x1, y1, x4, y4);

	/* Check if sum of A1, A2, A3 and A4
	 is same as A */
	return (A >= A1 + A2 + A3 + A4);
}

static GSList *load_catalog(const gchar *catalogue) {
	GFile *file;
	gchar *line;
	GSList *list = NULL;
	GError *error = NULL;

	file = g_file_new_build_filename(siril_get_system_data_dir(), "catalogue", catalogue, NULL);
	GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, &error);

	if (input_stream == NULL) {
		if (error != NULL) {
			g_clear_error(&error);
			siril_log_message(_("File [%s] does not exist\n"), g_file_peek_path(file));
		}
		g_object_unref(file);
		return NULL;
	}

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL,
				NULL, NULL))) {
		if (g_str_has_prefix (line, "Code")) {
			continue;
		}
		CatalogObjects *object = g_new(CatalogObjects, 1);
		gchar **token = g_strsplit(line, "\t", -1);
		gint nargs = g_strv_length(token);
		object->code = g_strdup(token[0]);
		object->ra = g_ascii_strtod(token[1], NULL) * 15.0;
		object->dec = g_strcmp0(token[2], "-") ? g_ascii_strtod(token[3], NULL) : g_ascii_strtod(token[3], NULL) * -1.0;
		object->radius = g_ascii_strtod(token[4], NULL) * 0.5;
		if (nargs > 6) {
			object->name = g_strdup(token[6]);
		} else {
			object->name = NULL;
		}

		list = g_slist_prepend(list, (gpointer) object);

		g_strfreev(token);
		g_free(line);
	}
	list = g_slist_reverse(list);

	g_object_unref(file);
	return list;
}

GSList *find_objects(fits *fit) {
	if (!has_wcs()) return NULL;
	GSList *targets = NULL;

	for (int i = 0; i < G_N_ELEMENTS(cat); i++) {
		GSList *list = load_catalog(cat[i]);

		for (GSList *l = list; l; l = l->next) {
			gdouble x1, y1, x2, y2, x3, y3, x4, y4;

			pix2wcs(0, fit->ry, &x1, &y1);
			pix2wcs(fit->rx, fit->ry, &x2, &y2);
			pix2wcs(0, 0, &x3, &y3);
			pix2wcs(fit->rx, 0, &x4, &y4);

			CatalogObjects *cur = (CatalogObjects *)l->data;

			printf("%s: %lf, %lf, %lf, %lf\t%lf %lf\n", cur->code, x1, y1, x4, y4, cur->ra, cur->dec);
			if (is_in_picture(x1, y1, x2, y2, x3, y3, x4, y4, cur->ra, cur->dec)) {
				CatalogObjects *new_object = new_catalog_object(cur->code, cur->ra, cur->dec, cur->radius, cur->name);
				targets = g_slist_prepend(targets, new_object);
			}
		}

		g_slist_free_full(list, (GDestroyNotify) free_object);
	}

	if (targets) {
		targets = g_slist_reverse(targets);
	}
	return targets;
}

gchar *get_catalogue_object_code(CatalogObjects *object) {
	return object->code;
}

gchar *get_catalogue_object_name(CatalogObjects *object) {
	return object->name;
}

gdouble get_catalogue_object_ra(CatalogObjects *object) {
	return object->ra;
}

gdouble get_catalogue_object_dec(CatalogObjects *object) {
	return object->dec;
}

gdouble get_catalogue_object_radius(CatalogObjects *object) {
	return object->radius;
}

void free_object(CatalogObjects *object) {
	g_free(object->code);
	g_free(object->name);
	g_free(object);
}
