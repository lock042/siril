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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_world_cs.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "gui/message_dialog.h"
#include "gui/utils.h"
#include "algos/siril_wcs.h"
#include "gui/image_display.h"

#include "annotate.h"

#define CATALOG_DIST_EPSILON (1/3600.0)	// 1 arcsec

static GSList *siril_catalogue_list = NULL; // loaded data from all annotation catalogues

static gboolean show_catalog(int catalog);

static const gchar *cat[] = {
	"messier.txt",
	"ngc.txt",
	"ic.txt",
	"ldn.txt",
	"sh2.txt",
	"stars.txt",
	/* below this line, user catalogues are in the user directory, not the
	 * siril install dir. See also USER_DSO_CAT_INDEX which gives the index
	 * of this separation */
	"user-DSO-catalogue.txt",
	"user-SSO-catalogue.txt"
};
// update USER_DSO_CAT_INDEX and USER_SSO_CAT_INDEX from the .h in case of change

struct _CatalogObjects {
	gchar *code;	// displayed name
	gdouble ra;	// in degrees, but in the files it's in hours...
	gdouble dec;	// degrees
	gdouble radius;	// in degrees but in the files it's the diameter. 0 for point-like,
			// negative for no accurate size
	gchar *name;
	gchar *alias;
	gint catalogue;	// index from the list of catalogues above
};

static CatalogObjects* new_catalog_object(const gchar *code, gdouble ra,
		gdouble dec, gdouble radius, const gchar *name, const gchar *alias,
		gint catalogue) {
	CatalogObjects *object = g_new(CatalogObjects, 1);

	object->code = g_strdup(code);
	object->ra = ra;
	object->dec = dec;
	object->radius = radius;
	object->name = g_strdup(name);
	object->alias = g_strdup(alias);
	object->catalogue = catalogue;
	return object;
}

static const char *cat_index_to_name(int index) {
	switch (index) {
		case 0:
			return "Messier";
		case 1:
			return "NGC";
		case 2:
			return "IC";
		case 3:
			return "LDN";
		case 4:
			return "Sh2";
		case 5:
			return "stars";
		case 6:
			return "user-DSO";
		case 7:
			return "user-SSO";
		case 8:
			return "user-temp";
		default:
			return "(undefined)";
	}
}

gboolean is_inside(fits *fit, double ra, double dec) {
	return wcs2pix(fit, ra, dec, NULL, NULL) == 0;
}

/* compare two objects, looking for duplicates based on the alias names againts the code of the object to search */
static gint object_compare(gconstpointer *a, gconstpointer *b) {
	const CatalogObjects *s1 = (const CatalogObjects *) a;
	const CatalogObjects *s2 = (const CatalogObjects *) b;

	if (!s1->alias) return 1;

	if (!strchr(s1->alias, '/'))
		return g_strcmp0(s1->alias, s2->code);

	gchar **token = g_strsplit(s1->alias, "/", -1);
	guint i = 0;
	while (token[i]) {
		if (!g_strcmp0(token[i], s2->code)) {
			g_strfreev(token);
			return 0;
		}
		i++;
	}
	g_strfreev(token);
	return 1;
}

/**
 * Catalogue must be formatted as follow:
 * This is a "; separated value" file
 * code ; ra ; dec sign ; dec ; radius ; name ; alias
 *
 */
static GSList *load_catalog(const gchar *filename, gint cat_index) {
	GFile *file;
	gchar *line;
	GSList *list = NULL;

	file = g_file_new_for_path(filename);
	GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, NULL);

	if (input_stream == NULL) {
		siril_log_message(_("File [%s] does not exist\n"), g_file_peek_path(file));
		g_object_unref(file);
		return NULL;
	}

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL,
				NULL, NULL))) {
		if (g_str_has_prefix (line, "Code")) {
			g_free(line);
			continue;
		}
		gchar **token = g_strsplit(line, ";", -1);
		guint nargs = g_strv_length(token);

		const gchar *code = NULL, *name = NULL, *alias = NULL;
		gdouble ra, dec, radius;

		/* mandatory tokens */
		code = token[0];
		ra = g_ascii_strtod(token[1], NULL) * 15.0;
		dec = g_strcmp0(token[2], "-") ? g_ascii_strtod(token[3], NULL) : g_ascii_strtod(token[3], NULL) * -1.0;
		radius = g_ascii_strtod(token[4], NULL) * 0.5;

		/* optional tokens */
		if (nargs > 5) {
			name = token[6];
			if (nargs > 6) {
				alias = token[7];
			}
		}

		CatalogObjects *object = new_catalog_object(code, ra, dec, radius, name, alias, cat_index);

		list = g_slist_prepend(list, (gpointer) object);

		g_strfreev(token);
		g_free(line);
	}
	list = g_slist_reverse(list);
	siril_debug_print("loaded %d objects from annotations catalogue %s\n", g_slist_length(list), filename);

	g_object_unref(data_input);
	g_object_unref(input_stream);
	g_object_unref(file);
	return list;
}

// for the show command, load a list of objects from a file into the temp cat
void load_csv_targets_to_temp(const gchar *filename) {
	GFile *file;
	gchar *line;
	GSList *list = NULL;

	file = g_file_new_for_path(filename);
	GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, NULL);

	if (input_stream == NULL) {
		siril_log_message(_("File [%s] does not exist\n"), g_file_peek_path(file));
		g_object_unref(file);
		return;
	}

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL,
				NULL, NULL))) {
		if (line[0] == '#') {
			g_free(line);
			continue;
		}
		gchar **token = g_strsplit(line, ",", -1);
		guint nargs = g_strv_length(token);
		if (nargs < 3) {
			siril_log_message(_("Malformed line: %s\n"), line);
			g_strfreev(token);
			g_free(line);
			continue;
		}

		/* mandatory tokens */
		const gchar *code = token[0];
		double ra, dec;
		if (strchr(token[1], ' ') || strchr(token[1], ':')) {
			ra = parse_ra_hms(token[1]);	// in hours
			dec = parse_dec_dms(token[2]);
		} else {
			ra = g_ascii_strtod(token[1], NULL);	// in degrees
			dec = g_ascii_strtod(token[2], NULL);
		}
		if (isnan(ra) || isnan(dec) || (ra == 0.0 && dec == 0.0)) {
			siril_log_message(_("Malformed line: %s\n"), line);
		} else {
			CatalogObjects *object = new_catalog_object(code, ra, dec, 0.0, NULL, NULL, USER_TEMP_CAT_INDEX);
			list = g_slist_prepend(list, (gpointer) object);

		}
		g_strfreev(token);
		g_free(line);
	}
	//list = g_slist_reverse(list);
	siril_debug_print("loaded %d objects from CSV temporary annotation %s\n", g_slist_length(list), filename);
	if (list)
		siril_catalogue_list = g_slist_concat(list, siril_catalogue_list);

	g_object_unref(data_input);
	g_object_unref(input_stream);
	g_object_unref(file);
}

static void reload_all_catalogues() {
	siril_debug_print("reloading annotation catalogues\n");
	if (siril_catalogue_list) {
		g_slist_free(com.found_object);
		com.found_object = NULL;

		g_slist_free_full(siril_catalogue_list, (GDestroyNotify)free_catalogue_object);
		siril_catalogue_list = NULL;
	}

	int cat_size = G_N_ELEMENTS(cat);
	for (int i = 0; i < cat_size; i++) {
		gchar *filename;
		if (i < USER_DSO_CAT_INDEX)
			filename = g_build_filename(siril_get_system_data_dir(), "catalogue", cat[i], NULL);
		else filename = g_build_filename(siril_get_config_dir(), PACKAGE, "catalogue", cat[i], NULL);
		if (g_file_test(filename, G_FILE_TEST_EXISTS))
			siril_catalogue_list = g_slist_concat(siril_catalogue_list, load_catalog(filename, i));
		g_free(filename);
	}
}

static GSList *get_siril_catalogue_list() {
	return siril_catalogue_list;
}

static gboolean is_catalogue_loaded() {
	return siril_catalogue_list != NULL;
}

typedef struct {
	const char *greek;		// Greek letter of stars
	const char *latin;		// Greek letter written in Latin
} GreekLetters;

static GreekLetters convert_to_greek[] = {
	{ "\u03b1", "alf" },
	{ "\u03b2", "bet" },
	{ "\u03b3", "gam" },
	{ "\u03b4", "del" },
	{ "\u03b5", "eps" },
	{ "\u03b6", "zet" },
	{ "\u03b7", "eta" },
	{ "\u03b8", "tet" },
	{ "\u03b9", "iot" },
	{ "\u03ba", "kap" },
	{ "\u03bb", "lam" },
	{ "\u03bc", "mu." },
	{ "\u03bd", "nu." },
	{ "\u03be", "ksi" },
	{ "\u03bf", "omi" },
	{ "\u03c0", "pi." },
	{ "\u03c1", "rho" },
	{ "\u03c3", "sig" },
	{ "\u03c4", "tau" },
	{ "\u03c5", "ups" },
	{ "\u03c6", "phi" },
	{ "\u03c7", "chi" },
	{ "\u03c8", "psi" },
	{ "\u03c9", "ome" },
	{ NULL, NULL }
};

static gchar* replace_str(const gchar *s, const gchar *old, const gchar *new) {
	gchar *result;
	int i, cnt = 0;
	int newlen = strlen(new);
	int oldlen = strlen(old);

	// Counting the number of times old word
	// occur in the string
	for (i = 0; s[i] != '\0'; i++) {
		if (g_strstr_len(&s[i], -1, old) == &s[i]) {
			cnt++;

			// Jumping to index after the old word.
			i += oldlen - 1;
		}
	}

	// Making new string of enough length
	result = malloc(i + cnt * (newlen - oldlen) + 1);

	i = 0;
	while (*s) {
		// compare the substring with the result
		if (g_strstr_len(s, -1, old) == s) {
			strcpy(&result[i], new);
			i += newlen;
			s += oldlen;
		} else
			result[i++] = *s++;
	}

	result[i] = '\0';
	return result;
}

gchar *retrieve_site_coord(fits *fit) {
	//  In case no localisation data is available, set the reference to geocentric, ie @500 in the request
	if (!fit->sitelat || !fit->sitelong)
		return g_strdup("@500");

	GString *formatted_site_coord = g_string_new("");

	if (fit->sitelat < 0.)
		formatted_site_coord = g_string_append(formatted_site_coord, "%2D");
	else formatted_site_coord = g_string_append(formatted_site_coord, "%2B");
	g_string_append_printf(formatted_site_coord, "%f,", fabs(fit->sitelat));

	if (fit->sitelong < 0.)
		formatted_site_coord = g_string_append(formatted_site_coord, "%2D");
	else formatted_site_coord = g_string_append(formatted_site_coord, "%2B");
	g_string_append_printf(formatted_site_coord, "%f,", fabs(fit->sitelong));

	g_string_append_printf(formatted_site_coord, "%f", fabs(fit->siteelev));

	return g_string_free(formatted_site_coord, FALSE);
}

static void write_in_user_catalogue(CatalogObjects *object, int cat_index) {
	// the temporary objects catalogue has no disk existence
	if (cat_index >= USER_TEMP_CAT_INDEX) return;

	/* First we test if root directory already exists */
	gchar *root = g_build_filename(siril_get_config_dir(), PACKAGE, "catalogue", NULL);

	if (!g_file_test(root, G_FILE_TEST_EXISTS)) {
		if (g_mkdir_with_parents(root, 0755) < 0) {
			siril_log_color_message(_("Cannot create output folder: %s\n"), "red", root);
			g_free(root);
			return;
		}
	}

	GFile *file = g_file_new_build_filename(root, cat[cat_index], NULL);
	g_free(root);

	GError *error = NULL;
	GOutputStream *output_stream;
	output_stream = (GOutputStream *)g_file_append_to(file, G_FILE_CREATE_NONE, NULL, &error);
	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
		}
		g_object_unref(file);
		return;
	}

	gchar sign = object->dec < 0 ? '-' : '+';
	gchar *output_line = g_strdup_printf("%s;%lf;%c;%lf;;;;\n", object->code, object->ra / 15.0, sign, fabs(object->dec));

	if (g_output_stream_write_all(output_stream, output_line, strlen(output_line), NULL, NULL, NULL) == FALSE) {
		siril_log_color_message(_("Error: failed to write output to user catalogue.\n"), "red");
	}

	g_free(output_line);
	g_object_unref(output_stream);
	g_object_unref(file);
}

/* get a list of objects from all catalogues (= from siril_catalogue_list) that
 * are framed in the passed plate solved image */
GSList *find_objects(fits *fit) {
	if (!has_wcs(fit)) return NULL;
	GSList *targets = NULL;

	if (!is_catalogue_loaded())
		reload_all_catalogues();
	GSList *list = get_siril_catalogue_list();

	for (GSList *l = list; l; l = l->next) {
		CatalogObjects *cur = (CatalogObjects *)l->data;
		if (!show_catalog(cur->catalogue)) continue;

		/* Search for objects in the image */
		if (is_inside(fit, cur->ra, cur->dec)) {
			if (!g_slist_find_custom(targets, cur, (GCompareFunc) object_compare)) {
				targets = g_slist_prepend(targets, cur);
			}
		}
	}

	if (targets) {
		targets = g_slist_reverse(targets);
	}
	return targets;
}

void add_object_in_catalogue(gchar *code, SirilWorldCS *wcs, gboolean is_solar_system, gboolean is_in_field) {
	int cat_index;

	if (!is_catalogue_loaded())
		reload_all_catalogues();
	/* check for the object first to avoid duplicates, not for solar system objects */
	if (!is_solar_system) {
		GSList *cur = siril_catalogue_list;
		double ra = siril_world_cs_get_alpha(wcs);
		double dec = siril_world_cs_get_delta(wcs);
		while (cur) {
			CatalogObjects *obj = cur->data;
			if (fabs(obj->ra - ra) < CATALOG_DIST_EPSILON && fabs(obj->dec - dec) < CATALOG_DIST_EPSILON) {
				siril_log_message(_("The object was already found in the %s catalog under the name %s, not adding it again\n"), cat_index_to_name(obj->catalogue), obj->code);
				return;
			}
			cur = cur->next;
		}
	}

	if (is_in_field) cat_index = USER_TEMP_CAT_INDEX;
	else {
		if (is_solar_system) cat_index = USER_SSO_CAT_INDEX;
		else cat_index = USER_DSO_CAT_INDEX;
	}

	CatalogObjects *new_object = new_catalog_object(code,
			siril_world_cs_get_alpha(wcs), siril_world_cs_get_delta(wcs), 0,
			NULL, NULL, cat_index);

	siril_catalogue_list = g_slist_prepend(siril_catalogue_list, new_object);
	write_in_user_catalogue(new_object, cat_index);
}

gchar *get_catalogue_object_code(CatalogObjects *object) {
	gboolean found = FALSE;
	int i = 0;

	/* in case of stars we want to convert to Greek letter */
	while (convert_to_greek[i].latin) {
		gchar *latin_code = g_strstr_len(object->code, -1, convert_to_greek[i].latin);
		if (latin_code) {
			found = TRUE;
			break;
		}
		i++;
	}
	if (found) {
		gchar *tmp = replace_str(object->code, convert_to_greek[i].latin, convert_to_greek[i].greek);
		gchar *code = g_strdup(tmp);
		g_free(tmp);
		g_free(object->code);
		object->code = code;
	}
	return object->code;
}

gchar *get_catalogue_object_name(CatalogObjects *object) {
	return object->name;
}

guint get_catalogue_object_cat(CatalogObjects *object) {
	return object->catalogue;
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

void free_catalogue_object(CatalogObjects *object) {
	g_free(object->code);
	g_free(object->name);
	g_free(object);
}

void refresh_found_objects() {
	if (has_wcs(&gfit)) {
		if (com.found_object) {
			g_slist_free(com.found_object);
		}
		com.found_object = find_objects(&gfit);
	}
}

static gboolean show_catalog(int catalog) {
	return com.pref.gui.catalog[catalog];
}

static void remove_temp_from_found(gpointer data, gpointer user_data) {
	CatalogObjects *obj = (CatalogObjects *) data;
	if (obj->catalogue == USER_TEMP_CAT_INDEX)
		com.found_object = g_slist_remove(com.found_object, data);
	// this is a very inefficient way to remove objects from a list,
	// but ok here becase there are not many
}

void purge_temp_user_catalogue() {
	g_slist_foreach(com.found_object, remove_temp_from_found, NULL);

	GSList *cur = siril_catalogue_list, *prev = NULL;
	while (cur) {
		CatalogObjects *obj = cur->data;
		if (obj->catalogue == USER_TEMP_CAT_INDEX) {
			// remove and free
			if (cur == siril_catalogue_list) {
				siril_catalogue_list = cur->next;
				free_catalogue_object(obj);
				g_slist_free_1(cur);
				cur = siril_catalogue_list;
			} else {
				GSList *tmp = cur;
				prev->next = cur->next;
				free_catalogue_object(obj);
				g_slist_free_1(tmp);
			}
		}
		else {
			prev = cur;
			cur = cur->next;
		}
	}
}

void on_purge_SSO_user_catalogue_clicked(GtkButton *button, gpointer user_data) {
	int confirm = siril_confirm_dialog(_("Catalog deletion"),
			_("You are about to purge your SOLAR SYSTEM user catalog. This means the "
				"file containing the manually added objects will be deleted. "
				"This operation cannot be undone."), _("Purge Catalogue"));
	if (!confirm) {
		return;
	}

	GFile *file = g_file_new_build_filename(siril_get_config_dir(),
			PACKAGE, "catalogue", cat[USER_SSO_CAT_INDEX], NULL);

	g_autoptr(GError) local_error = NULL;
	if (g_file_query_exists(file, NULL) && !g_file_delete(file, NULL, &local_error)) {
		g_warning("Failed to delete %s: %s", g_file_peek_path(file), local_error->message);
	} else {
		reload_all_catalogues();
	}

	g_object_unref(file);
}

void on_purge_DSO_user_catalogue_clicked(GtkButton *button, gpointer user_data) {
	int confirm = siril_confirm_dialog(_("Catalogue deletion"),
			_("You are about to purge your DEEP SKY user catalogue. This means the "
				"file containing the manually added objects will be deleted. "
				"This operation cannot be undone."), _("Purge Catalogue"));
	if (!confirm) {
		return;
	}

	GFile *file = g_file_new_build_filename(siril_get_config_dir(), PACKAGE, "catalogue", cat[USER_DSO_CAT_INDEX], NULL);

	g_autoptr(GError) local_error = NULL;
	if (g_file_query_exists(file, NULL) && !g_file_delete(file, NULL, &local_error)) {
		g_warning("Failed to delete %s: %s", g_file_peek_path(file), local_error->message);
	} else {
		reload_all_catalogues();
	}

	g_object_unref(file);
}

