/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_world_cs.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "algos/siril_wcs.h"
#include "io/siril_catalogues.h"

#include "annotation_catalogues.h"

#define CATALOG_DIST_EPSILON (1/3600.0)	// 1 arcsec or 1s in hrs

static GSList *siril_annot_catalogue_list = NULL; // loaded data from all annotation catalogues
static gboolean get_annotation_visibility(siril_cat_index cat_index);

static const gchar *cat[] = {
	"messier.csv",
	"ngc.csv",
	"ic.csv",
	"ldn.csv",
	"sh2.csv",
	"stars.csv",
	"constellations.csv",
	"constellationsnames.csv",
	/* below this line, user catalogues are in the user directory, not the
	 * siril install dir.  */
	"user-DSO-catalogue.csv",
	"user-SSO-catalogue.csv"
};
// make sure com.pref.gui.catalog matches this too

static const gchar *get_cat_filename_by_index(siril_cat_index cat_index) {
	if (cat_index < CAT_AN_MESSIER || cat_index > CAT_AN_USER_SSO)
		return NULL;
	return cat[cat_index - CAT_AN_INDEX_OFFSET];
}

struct _CatalogObjects {
	gchar *code;	// displayed name
	gchar *pretty_code;	// name with greek characters
	double x;	// in fits_display coords
	double y;	// in fits_display coords
	double x1;	// in fits_display coords
	double y1;	// in fits_display coords
	gdouble radius;	// in degrees but in the files it's the diameter. 0 for point-like,
			// negative for no accurate size
	gchar *alias;
	siril_cat_index catalogue; // index from the list of catalogues above
};

static CatalogObjects* new_catalog_object(const gchar *name, double x,
		double y, double x1,
		double y1, double radius, const gchar *alias,
		siril_cat_index catalogue) {
	CatalogObjects *object = g_new(CatalogObjects, 1);

	object->code = (name) ? g_strdup(name) : NULL;
	object->pretty_code = NULL;
	object->x = x;
	object->y = y;
	object->x1 = x1;
	object->y1 = y1;
	object->radius = radius;
	object->alias = (alias) ? g_strdup(alias) : NULL;
	object->catalogue = catalogue;
	return object;
}

static int compare_names(gchar *s1, gchar *s2) {
	if (!s1)
		return 1;
	if (!s2)
		return -1;
	return g_ascii_strcasecmp(s1, s2);
}

/* compare two objects, looking for duplicates based on the alias names against the name of the object to search */
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

// returns true if it was added
static gboolean add_alias_to_item(cat_item *item, gchar *name) {
	if (!compare_names(item->name, name))
		return FALSE;
	if (item->alias && item->alias[0] != '\0') {
		if (!strchr(item->alias, '/')) {
			if (!compare_names(item->alias, name))
				return FALSE;
		} else {
			gchar **token = g_strsplit(item->alias, "/", -1);
			guint i = 0;
			while (token[i]) {
				if (!compare_names(token[i], name)) {
					g_strfreev(token);
					return FALSE;
				}
				i++;
			}
			g_strfreev(token);
		}
		gchar *oldalias = item->alias;
		item->alias = g_strjoin("/", item->alias, name, NULL);
		g_free(oldalias);
		siril_debug_print("new alias: '%s'\n", item->alias);
		return TRUE;
	}
	item->alias = g_strdup(name);
	return TRUE;
}

static GSList *get_siril_annot_catalogue_list() {
	return siril_annot_catalogue_list;
}

static gboolean is_catalogue_loaded() {
	return siril_annot_catalogue_list != NULL;
}

static void free_catalogue_object(CatalogObjects *object) {
	g_free(object->code);
	g_free(object->alias);
	g_free(object);
}

gchar *get_annotation_catalog_filename(siril_cat_index cat_index, gboolean for_reading) {
	if (cat_index < CAT_AN_MESSIER || cat_index > CAT_AN_USER_SSO) {
		siril_debug_print("Wrong catalog index\n");
		return NULL;
	}
	gchar *filename = NULL;
	if (cat_index < CAT_AN_USER_DSO)
		filename = g_build_filename(siril_get_system_data_dir(), "catalogue", get_cat_filename_by_index(cat_index), NULL);
	else
		filename = g_build_filename(siril_get_config_dir(), PACKAGE, "catalogue", get_cat_filename_by_index(cat_index), NULL);
	if (for_reading && !g_file_test(filename, G_FILE_TEST_EXISTS)) {
		siril_debug_print("Catalog file %s does not exist\n", get_cat_filename_by_index(cat_index));
		g_free(filename);
		filename = NULL;
	}
	return filename;
}
/*
 * Loads a csv catalogue using generic csv parser
 */
static annotations_catalogue_t *load_catalog(siril_cat_index cat_index, const gchar *filename) {
	gboolean islocal = !filename;
	siril_catalogue *siril_cat = siril_catalog_new(cat_index);
	if (islocal)
		filename = get_annotation_catalog_filename(cat_index, TRUE);
	if (!filename || siril_catalog_load_from_file(siril_cat, filename)) {// use the generic csv parser
		siril_debug_print("Could not load the catalog %s\n", (islocal) ? get_cat_filename_by_index(cat_index) : filename);
		siril_catalog_free(siril_cat);
		return NULL;
	}
	siril_debug_print("loaded %d objects from annotations catalogue %s\n", siril_cat->nbitems, filename);
	annotations_catalogue_t *annot_cat = g_new(annotations_catalogue_t, 1);
	annot_cat->cat = siril_cat;
	annot_cat->show = get_annotation_visibility(cat_index);
	return annot_cat;
}

// loads all catalogues from files
static void load_all_catalogues() {
	siril_debug_print("reloading annotation catalogues\n");
	int cat_size = G_N_ELEMENTS(cat);
	for (int i = 0; i < cat_size; i++) {
		annotations_catalogue_t *newcat = load_catalog(i + CAT_AN_INDEX_OFFSET, NULL);
		if (newcat)
			siril_annot_catalogue_list = g_slist_prepend(siril_annot_catalogue_list, newcat);
	}
	siril_annot_catalogue_list = g_slist_reverse(siril_annot_catalogue_list);
}

static GSList *find_catalogue_by_index(siril_cat_index cat_index) {
	if (cat_index < CAT_AN_MESSIER|| cat_index > CAT_AN_USER_TEMP) {
		siril_debug_print("Wrong catalog index\n");
		return NULL;
	}
	if (!is_catalogue_loaded())
		load_all_catalogues();
	GSList *cur = siril_annot_catalogue_list;
	while (cur) {
		annotations_catalogue_t *curcat = cur->data;
		if (curcat->cat->cat_index == cat_index)
			return cur;
		cur = cur->next;
	}
	return NULL;
}

void set_annotation_visibility(siril_cat_index cat_index, gboolean visible) {
	if (cat_index < CAT_AN_MESSIER || cat_index > CAT_AN_USER_TEMP)
		return;
	com.pref.gui.catalog[cat_index - CAT_AN_INDEX_OFFSET] = visible;
}

static gboolean get_annotation_visibility(siril_cat_index cat_index) {
	return com.pref.gui.catalog[cat_index - CAT_AN_INDEX_OFFSET];
}

void refresh_annotation_visibility() {
	GSList *cur = siril_annot_catalogue_list;
	while (cur) {
		annotations_catalogue_t *curcat = cur->data;
		curcat->show = get_annotation_visibility(curcat->cat->cat_index);
		cur = cur->next;
	}
}

void refresh_annotation_to_temp() {
	GSList *cur = siril_annot_catalogue_list;
	while (cur) {
		annotations_catalogue_t *curcat = cur->data;
		curcat->show = (curcat->cat->cat_index == CAT_AN_USER_TEMP);
		cur = cur->next;
	}
}

int load_siril_cat_to_temp(siril_catalogue *siril_cat) {
	if (!siril_cat)
		return 1;
	if (!is_catalogue_loaded())
		load_all_catalogues();
	annotations_catalogue_t *annot_cat = g_new(annotations_catalogue_t, 1);
	annot_cat->cat = siril_cat;
	annot_cat->show = TRUE;
	annot_cat->cat->cat_index = CAT_AN_USER_TEMP;
	siril_annot_catalogue_list = g_slist_append(siril_annot_catalogue_list, annot_cat);
	return 0;
}

// for the show command, loads a list of objects from a file into the temp cat
// The file must have at least name, ra and dec columns
int load_csv_targets_to_temp(const gchar *filename) {
	if (!is_catalogue_loaded())
		load_all_catalogues();
	annotations_catalogue_t *annot_cat = load_catalog(CAT_AN_USER_TEMP, filename);
	if (annot_cat) {
		siril_annot_catalogue_list = g_slist_append(siril_annot_catalogue_list, annot_cat);
		siril_debug_print("loaded %d objects from CSV temporary annotation %s\n", annot_cat->cat->nbitems, filename);
		return 0;
	}
	return 1;
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

//TODO: may be much simpler to use g_string_replace
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
	result = g_malloc(i + cnt * (newlen - oldlen) + 1);

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

gchar *get_catalogue_object_code_pretty(CatalogObjects *object) {
	gboolean found = FALSE;
	int i = 0;
	if (!object->pretty_code) {
		if (!object->code)
			return NULL;
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
			object->pretty_code = replace_str(object->code,
					convert_to_greek[i].latin, convert_to_greek[i].greek);
		} else {
			object->pretty_code = g_strdup(object->code);
		}
	}
	return object->pretty_code;
}

static void write_in_user_catalogue(siril_cat_index cat_index) {;
	// only these two catalogs are writable
	if (cat_index != CAT_AN_USER_DSO && cat_index != CAT_AN_USER_SSO)
		return;
	GSList *cur = find_catalogue_by_index(cat_index);
	if (!cur || !cur->data)
		return;
	gchar *filename = get_annotation_catalog_filename(cat_index, FALSE);
	annotations_catalogue_t *annot_cat = cur->data; 
	siril_catalogue *siril_cat = annot_cat->cat;
	if (!siril_catalog_write_to_file(siril_cat, filename)) {
		siril_log_color_message(_("Could not write the updated catalogue %s\n"), "red", get_cat_filename_by_index(cat_index));
	}
}


gboolean is_inside(fits *fit, double ra, double dec) {
	return wcs2pix(fit, ra, dec, NULL, NULL) == 0;
}

/* get a list of objects from all catalogues (= from siril_annot_catalogue_list) that
 * are framed in the passed plate solved image */
GSList *find_objects_in_field(fits *fit) {
	if (!has_wcs(fit))
		return NULL;
	GSList *targets = NULL;
	if (!is_catalogue_loaded())
		load_all_catalogues();
	gdouble starradius = get_wcs_image_resolution(fit) * 6. * 60.0;
	GSList *list = get_siril_annot_catalogue_list();
	double tref = date_time_to_Julian(fit->keywords.date_obs);
	for (GSList *l = list; l; l = l->next) {
		annotations_catalogue_t *curcat = l->data;
		if (!curcat->show) // the catalog show member is set at read-out
			continue;
		siril_catalogue *siril_cat = curcat->cat;
		gboolean is_star_cat = is_star_catalogue(siril_cat->cat_index);
		if (siril_cat->projected != CAT_PROJ_WCS) {
			siril_catalog_project_with_WCS(siril_cat,fit, TRUE, TRUE); // sanity check will be done during the projection
		}
		for (int i = 0; i < siril_cat->nbitems; i++) {
			if (siril_cat->cat_index == CAT_AN_USER_SSO) { // we need to check the record is from the same night and same location
				if (tref == 0. || fabs(siril_cat->cat_items[i].dateobs - tref) > 0.75 || // 18hrs
					fabs(siril_cat->cat_items[i].sitelat - fit->keywords.sitelat) > 1.e-4 ||
					fabs(siril_cat->cat_items[i].sitelon - fit->keywords.sitelong) > 1.e-4 ||
					fabs(siril_cat->cat_items[i].siteelev - fit->keywords.siteelev) > 1.)
					continue;
			}

			if (siril_cat->cat_items[i].included) { //included means it is within the bounds of the image after projection
				double x = 0., y = 0., x1 = 0., y1 = 0.;
				// we write directly in display coordinates to avoid the flip at every redraw
				siril_to_display(siril_cat->cat_items[i].x, siril_cat->cat_items[i].y, &x, &y, fit->ry);
				x = siril_cat->cat_items[i].x;
				y = fit->ry - siril_cat->cat_items[i].y;
				if (siril_cat->cat_index == CAT_AN_CONST) {
					siril_to_display(siril_cat->cat_items[i].x1, siril_cat->cat_items[i].y1, &x1, &y1, fit->ry);
					x1 = siril_cat->cat_items[i].x1;
					y1 = fit->ry - siril_cat->cat_items[i].y1;
				}
				CatalogObjects *cur = new_catalog_object(
					siril_cat->cat_items[i].name,
					x,
					y,
					x1,
					y1,
					(is_star_cat) ? starradius : .5 * siril_cat->cat_items[i].diameter,
					siril_cat->cat_items[i].alias,
					siril_cat->cat_index
				);
				if (!g_slist_find_custom(targets, cur, (GCompareFunc)object_compare))
					targets = g_slist_prepend(targets, cur);
				else
					free_catalogue_object(cur);
			}
		}
	}
	if (targets) {
		targets = g_slist_reverse(targets);
	}
	return targets;
}

// compares (ra/dec) positions for all catalogues
// except for sso where name, dateobs and observation sites are compared
static gboolean is_same_item(cat_item *item1, cat_item *item2, siril_cat_index cat_index) {
	if (!item1 || !item2 || cat_index < CAT_AN_MESSIER || cat_index > CAT_AN_USER_TEMP) {
		siril_debug_print("problem when comparing two items\n");
		return FALSE;
	}
	switch (cat_index) {
		case CAT_AN_USER_SSO: // for SSO we need to check on more criteria than just the position
			return (!compare_names(item1->name, item2->name) || // same name
					!compare_names(item1->alias, item2->alias) || // same alias
					!compare_names(item1->name, item2->alias) || // or a mix'n'match
					!compare_names(item1->alias, item2->name))
					// same night (time diff < 18hrs), we can then use vra/vdec to compute position throughout the night
					&& fabs(item1->dateobs - item2->dateobs) < 0.75
					// same observation site
					&& fabs(item1->sitelat - item2->sitelat) < 1.e-4
					&& fabs(item1->sitelon - item2->sitelon) < 1.e-4
					&& fabs(item1->siteelev - item2->siteelev) < 1.;
		default:
			return fabs(item1->ra - item2->ra) < CATALOG_DIST_EPSILON
				&& fabs(item1->dec - item2->dec) < CATALOG_DIST_EPSILON;
	}
}

// adds an item in one of the catalogues of the static list
// and writes the updated file (for user catalogues only)
void add_item_in_catalogue(cat_item *item, siril_cat_index cat_index, gboolean check_duplicates) {
	if (!is_catalogue_loaded())
		load_all_catalogues();
	if (check_duplicates) {	// duplicate check based on coordinates within 1"
		GSList *cur = siril_annot_catalogue_list;
		while (cur) {
			annotations_catalogue_t *curcat = cur->data;
			if (curcat->cat->cat_index == CAT_AN_USER_TEMP) {
				cur = cur->next;
				continue;
			}
			for (int i = 0; i < curcat->cat->nbitems; i++) {
				if (is_same_item(item, &curcat->cat->cat_items[i], cat_index)) {
					siril_log_message(_("The object was already found in the %s catalog "
								"under the name '%s', not adding it again\n"),
								catalog_to_str(curcat->cat->cat_index), curcat->cat->cat_items[i].name);
					if (add_alias_to_item(&curcat->cat->cat_items[i], item->name)) {
						if (!(cat_index == CAT_AN_USER_DSO || cat_index == CAT_AN_USER_SSO)) {
							siril_debug_print("alias %s added to %s in the runtime catalogue\n",
							item->name, curcat->cat->cat_items[i].name);
							/* we add it to the catalogue in memory, but it will be lost on
							* purge/reload or when restarting siril. It's a bit complicated
							* to add to siril native catalogues as they should not be considered writeable.
							* at least it will help for consecutive requests
						*/
						} else { // we update the user catalogues if required
							write_in_user_catalogue(cat_index);
						}
					}
					return;
				}
			}
			cur = cur->next;
		}
	}
	siril_catalogue *siril_cat = NULL;
	GSList *cur = find_catalogue_by_index(cat_index);
	if (!cur || !cur->data) {// the catalog does not exist yet
		annotations_catalogue_t *annot_cat = calloc(1, sizeof(annotations_catalogue_t));
		siril_cat = siril_catalog_new(cat_index);
		annot_cat->cat = siril_cat;
		annot_cat->show = TRUE;
		siril_annot_catalogue_list = g_slist_append(siril_annot_catalogue_list, annot_cat);
	} else {
		siril_cat = ((annotations_catalogue_t *)cur->data)->cat;
	}
	if (siril_catalog_append_item(siril_cat, item) && (cat_index == CAT_AN_USER_DSO || cat_index == CAT_AN_USER_SSO))
		write_in_user_catalogue(cat_index);
}

// search DSO by name in all annotation catalogues, but remember the result does not
// contain magnitudes, only name and J2000 equatorial coordinates
// if *cat_index is non NULL, fills it with the catalogue type enum
// The caller must free the result
cat_item *search_in_annotations_by_name(const char *input, siril_cat_index *cat_index) {
	if (!is_catalogue_loaded())
		load_all_catalogues();
	if (!input || input[0] == '\0')
		return NULL;
	int len = strlen(input);
	gchar *target;
	// no space in our catalogues' codes
	if (len >= 3 && !g_ascii_strncasecmp(input, "M ", 2) && input[2] >= '0' && input[2] <= '9')
		target = g_strdup_printf("M%s", input+2);
	else if (len >= 4 && !g_ascii_strncasecmp(input, "IC ", 3) && input[3] >= '0' && input[3] <= '9')
		target = g_strdup_printf("IC%s", input+3);
	else if (len >= 5 && !g_ascii_strncasecmp(input, "NGC ", 4) && input[4] >= '0' && input[4] <= '9')
		target = g_strdup_printf("NGC%s", input+4);
	else {
		target = g_strdup(input);
		target[0] = g_ascii_toupper(target[0]);	// in the catalogues, they're all starting with a caps
	}
	siril_debug_print("target name after transformation: %s\n", target);

	cat_item *found = NULL, *probable = NULL;
	GSList *cur = siril_annot_catalogue_list;
	gboolean bingo = FALSE;
	while (cur) {
		annotations_catalogue_t *curcat = cur->data;
		if (curcat->cat->cat_index == CAT_AN_USER_SSO || curcat->cat->cat_index == CAT_AN_USER_TEMP) {
			cur = cur->next;
			continue;
		}
		siril_catalogue *siril_cat = curcat->cat;
		for (int i = 0; i < siril_cat->nbitems; i++) {
			cat_item *item = &siril_cat->cat_items[i];
			if (!compare_names(item->name, target)) {
				found = item;
				bingo = TRUE;
				if (cat_index)
					*cat_index = siril_cat->cat_index;
				break;
			}
			if (item->alias) {
				if (!strchr(item->alias, '/')) {
					if (!compare_names(item->alias, target)) {
						found = item;
						bingo = TRUE;
						if (cat_index)
							*cat_index = siril_cat->cat_index;
						break;
					}
				} else {
					gchar **token = g_strsplit(item->alias, "/", -1);
					guint i = 0;
					while (token[i]) {
						if (!compare_names(token[i], target)) {
							found = item;
							bingo = TRUE;
							if (cat_index)
								*cat_index = siril_cat->cat_index;
							break;
						}
						if (!probable && strstr(token[i], target)) {
							probable = item;
							if (cat_index)
								*cat_index = siril_cat->cat_index;
						}
						i++;
					}
					g_strfreev(token);
				}
			}
		}
		if (bingo)
			break;
		cur = cur->next;
	}

	g_free(target);

	if (found) {
		siril_debug_print("object found in annotation catalogues: %s\n", found->name);
		cat_item *output = calloc(1, sizeof(cat_item));
		siril_catalogue_copy_item(found, output);
		return output;
	}
	if (probable) {
		siril_debug_print("probable object found in annotation catalogues: %s\n", probable->name);
		cat_item *output = calloc(1, sizeof(cat_item));
		siril_catalogue_copy_item(probable, output);
		return output;
	}
	siril_debug_print("object %s not found in annotation catalogues\n", input);
	return NULL;
}

// search SSO user catalog by name /datobs/ site position
// The caller must free the result
cat_item *search_in_solar_annotations(sky_object_query_args *args) {
	if (!is_catalogue_loaded())
		load_all_catalogues();
	GSList *solar_an_cat = find_catalogue_by_index(CAT_AN_USER_SSO);
	if (!solar_an_cat) {
		siril_debug_print("no SSO catalogue\n");
		return NULL;
	}
	siril_catalogue *solar_cat = ((annotations_catalogue_t *)solar_an_cat->data)->cat;

	cat_item *ref_item = calloc(1, sizeof(cat_item));
	double ref = date_time_to_Julian(args->fit->keywords.date_obs);
	ref_item->name = g_strdup(args->name);
	ref_item->dateobs = ref;
	ref_item->sitelat = args->fit->keywords.sitelat;
	ref_item->sitelon = args->fit->keywords.sitelong;
	ref_item->siteelev = args->fit->keywords.siteelev;
	for (int i = 0; i < solar_cat->nbitems; i++) {
		if (is_same_item(ref_item, &solar_cat->cat_items[i], CAT_AN_USER_SSO)) {
			siril_catalogue_copy_item(&solar_cat->cat_items[i], ref_item);
			siril_log_message(_("Object %s record found in the local SSO catalog\n"), ref_item->name);
			GDateTime *timerecord = g_date_time_add_seconds(args->fit->keywords.date_obs, (ref_item->dateobs - ref) * 3600. * 24.); // doing this to avoid writing a converter julian to fits date...
			gchar *dt = date_time_to_FITS_date(timerecord);
			gchar *alpha = siril_world_cs_alpha_format_from_double(ref_item->ra, " %02dh%02dm%04.1lfs");
			gchar *delta = siril_world_cs_delta_format_from_double(ref_item->dec, "%c%02d°%02d\'%04.1lf\"");
			siril_log_message(_("at coordinates: %s, %s (on DATE-OBS:%s)\n"), alpha, delta, dt);
			g_free(alpha);
			g_free(delta);
			g_free(dt);
			g_date_time_unref(timerecord);
			return ref_item;
		}
	}
	siril_catalog_free_item(ref_item);
	free(ref_item);
	return NULL;
}


gchar *get_catalogue_object_code(const CatalogObjects *object) {
	return object->code;
}

siril_cat_index get_catalogue_object_cat(const CatalogObjects *object) {
	return object->catalogue;
}

gdouble get_catalogue_object_x(const CatalogObjects *object) {
	return object->x;
}

gdouble get_catalogue_object_y(const CatalogObjects *object) {
	return object->y;
}

gdouble get_catalogue_object_x1(const CatalogObjects *object) {
	return object->x1;
}

gdouble get_catalogue_object_y1(const CatalogObjects *object) {
	return object->y1;
}

gdouble get_catalogue_object_radius(const CatalogObjects *object) {
	return object->radius;
}

void refresh_found_objects() {
	if (has_wcs(&gfit)) {
		if (com.found_object) {
			g_slist_free_full(com.found_object, (GDestroyNotify)free_catalogue_object);
			com.found_object = NULL;
		}
		com.found_object = find_objects_in_field(&gfit);
	}
}

static void remove_user_cat_from_found(gpointer data, gpointer user_data) {
	siril_cat_index *cat_index = (siril_cat_index *)user_data;
	CatalogObjects *obj = (CatalogObjects *) data;
	if (obj->catalogue == *cat_index) {
		com.found_object = g_slist_remove(com.found_object, data);
	}
	// this is a very inefficient way to remove objects from a list,
	// but ok here because there are not many
}

void purge_user_catalogue(siril_cat_index cat_index) {
	// we remove the objects of type cat_index from found_objects
	g_slist_foreach(com.found_object, remove_user_cat_from_found, &cat_index);
	// we remove the catalog from the static list
	GSList *cur = siril_annot_catalogue_list;
	printf("nb_elts_at_start: %d\n", g_slist_length(siril_annot_catalogue_list));
	while (cur) {
		annotations_catalogue_t *curcat = cur->data;
		if (curcat->cat->cat_index == cat_index) { // we remove all catalogues of same index (there may be multiple temp cats)
			GSList *next = cur->next;
			siril_annot_catalogue_list = g_slist_remove_link(siril_annot_catalogue_list, cur);
			siril_catalog_free(curcat->cat);
			g_slist_free(cur);
			printf("nb_elts_after: %d\n", g_slist_length(siril_annot_catalogue_list));
			cur = next;
		} else {
			cur = cur->next;
		}
	}
	printf("nb_elts_on_exit: %d\n", g_slist_length(siril_annot_catalogue_list));
}

// to be called when we close an image with flag set to TRUE
void cleanup_annotation_catalogues(gboolean purge_temp) {
	g_slist_free_full(com.found_object, (GDestroyNotify)free_catalogue_object);
	com.found_object = NULL;
	if (purge_temp)
		purge_user_catalogue(CAT_AN_USER_TEMP);
	GSList *cur = siril_annot_catalogue_list;
	while (cur) {
		annotations_catalogue_t *curcat = cur->data;
		siril_catalog_reset_projection(curcat->cat);
		cur = cur->next;
	}
}

// to be called when the current image has been changed (geometry/undo/redo)
void refresh_annotations(gboolean purge_temp) {
	if (com.found_object) {
		cleanup_annotation_catalogues(purge_temp);
		refresh_found_objects();
	} else {
		cleanup_annotation_catalogues(purge_temp);
	}
}
