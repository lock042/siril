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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <unistd.h>
#endif
#endif

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "algos/search_objects.h"
#include "algos/siril_wcs.h"
#include "io/annotation_catalogues.h"
#include "algos/astrometry_solver.h"
#include "algos/comparison_stars.h"
#include "io/remote_catalogues.h"
#include "io/local_catalogues.h"
#include "registration/matching/misc.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"
#include "gui/utils.h"

// This list defines the columns that can possibly be found in any catalogue
const gchar *cat_columns[] = { 
	[CAT_FIELD_RA] = "ra",
	[CAT_FIELD_DEC] = "dec",
	[CAT_FIELD_PMRA] = "pmra",
	[CAT_FIELD_PMDEC] = "pmdec",
	[CAT_FIELD_MAG] = "mag",
	[CAT_FIELD_BMAG] = "bmag",
	[CAT_FIELD_E_MAG] = "e_mag",
	[CAT_FIELD_E_BMAG] = "e_bmag",
	[CAT_FIELD_NAME] = "name",
	[CAT_FIELD_DIAMETER] = "diameter",
	[CAT_FIELD_ALIAS] = "alias",
	[CAT_FIELD_DATEOBS] = "date-obs",
	[CAT_FIELD_SITELAT] = "sitelat",
	[CAT_FIELD_SITELON] = "sitelon",
	[CAT_FIELD_SITEELEV] = "siteelev",
	[CAT_FIELD_VRA] = "vra",
	[CAT_FIELD_VDEC] = "vdec",
	[CAT_FIELD_TYPE] = "type"
};

const gchar **get_cat_colums_names() {
	return cat_columns;
}

// This function returns the column index from a string
static int get_column_index(gchar *field) {
	for (int i = 0; i < MAX_CAT_COLUMNS; i++) {
		if (!strcasecmp(field, cat_columns[i])) // case insensitive version
			return i;
	}
	return -1;
}

static gchar *get_field_to_str(cat_item item, cat_fields field) {
	switch (field) {
		case CAT_FIELD_RA:
			return (item.ra) ? g_strdup_printf("%.6f", item.ra) : "";
		case CAT_FIELD_DEC:
			return (item.dec) ? g_strdup_printf("%.6f", item.dec) : "";
		case CAT_FIELD_PMRA:
			return (item.pmra) ? g_strdup_printf("%g", item.pmra) : "";
		case CAT_FIELD_PMDEC:
			return (item.pmdec) ? g_strdup_printf("%g", item.pmdec) : "";
		case CAT_FIELD_MAG:
			return (item.mag) ? g_strdup_printf("%g", item.mag) : "";
		case CAT_FIELD_BMAG:
			return (item.bmag) ? g_strdup_printf("%g", item.bmag) : "";
		case CAT_FIELD_E_MAG:
			return (item.e_mag) ? g_strdup_printf("%g", item.e_mag) : "";
		case CAT_FIELD_E_BMAG:
			return (item.e_bmag) ? g_strdup_printf("%g", item.e_bmag) : "";
		case CAT_FIELD_DIAMETER:
			return (item.diameter) ? g_strdup_printf("%g", item.diameter) : "";
		case CAT_FIELD_DATEOBS:
			return (item.dateobs) ? g_strdup_printf("%.12f", item.dateobs) : "";
		case CAT_FIELD_SITELAT:
			return (item.sitelat) ? g_strdup_printf("%g", item.sitelat) : "";
		case CAT_FIELD_SITELON:
			return (item.sitelon) ? g_strdup_printf("%g", item.sitelon) : "";
		case CAT_FIELD_SITEELEV:
			return (item.siteelev) ? g_strdup_printf("%g", item.siteelev) : "";
		case CAT_FIELD_VRA:
			return (item.vra) ? g_strdup_printf("%g", item.vra) : "";
		case CAT_FIELD_VDEC:
			return (item.vdec) ? g_strdup_printf("%g", item.vdec) : "";
		case CAT_FIELD_NAME:
			return (item.name) ? g_strdup(item.name) : "";
		case CAT_FIELD_ALIAS:
			return (item.alias) ? g_strdup(item.alias) : "";
		case CAT_FIELD_TYPE:
			return (item.type) ? g_strdup(item.type) : "";
		default:
			return NULL;
	}
}
// This function defines the fields that should be present in each catalog
// To be used as sanity check
uint32_t siril_catalog_columns(object_catalog cat) {
	switch (cat) {
		case CAT_TYCHO2:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG);
		case CAT_NOMAD:
		case CAT_GAIADR3:
		case CAT_SIMBAD:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_NAME);
		case CAT_PPMXL:
		case CAT_BSC:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME);
		case CAT_APASS:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_E_MAG) | (1 << CAT_FIELD_E_BMAG);
		case CAT_GCVS:
		case CAT_VSX:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME);
		case CAT_PGC:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME);
		case CAT_EXOPLANETARCHIVE:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME);
		case CAT_AAVSO_CHART:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_E_MAG) | (1 << CAT_FIELD_E_BMAG) | (1 << CAT_FIELD_NAME);
		case CAT_IMCCE:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_VRA) | (1 << CAT_FIELD_VDEC) | (1 << CAT_FIELD_TYPE);
		case CAT_AN_MESSIER:
		case CAT_AN_NGC:
		case CAT_AN_IC:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_ALIAS) | (1 << CAT_FIELD_DIAMETER) | (1 << CAT_FIELD_MAG);
		case CAT_AN_LDN:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME);
		case CAT_AN_SH2:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_ALIAS);
		case CAT_AN_STARS:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME); //| (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) //TODO:Add pm
		case CAT_AN_USER_DSO:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_ALIAS);
		case CAT_AN_USER_SSO:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_ALIAS) | (1 << CAT_FIELD_DATEOBS) | (1 << CAT_FIELD_SITELAT) | (1 << CAT_FIELD_SITELON) | (1 << CAT_FIELD_SITEELEV) | (1 << CAT_FIELD_VRA) | (1 << CAT_FIELD_VDEC) | (1 << CAT_FIELD_TYPE);
		case CAT_COMPSTARS:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_TYPE);
		case CAT_LOCAL:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG);
		case CAT_AN_USER_TEMP:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME);
		default:
			return 0;
	}
}

// This function compares two cat_item objects and return their order by mag
static int compare_items_by_mag(const void* item1, const void* item2) {
	cat_item *i1 = (cat_item*) item1;
	cat_item *i2 = (cat_item*) item2;
	if (i1->mag < i2->mag)
		return -1;
	if (i1->mag > i2->mag)
		return 1;
	return 0;
}

// This function sorts the *cat_item list of a siril_catalogue by magnitude 
void sort_cat_items_by_mag(siril_catalogue *siril_cat) {
	if (siril_cat && siril_cat->nbitems > 0 && (siril_cat->columns&(1 << CAT_FIELD_MAG)))
		qsort(siril_cat->cat_items, siril_cat->nbitems, sizeof(cat_item), compare_items_by_mag);
}

// returns the default magnitude to be used in conesearch
static float siril_catalog_get_default_limit_mag(object_catalog cat) {
	switch(cat) {
		case CAT_IMCCE:
			return 20.f;
		case CAT_AAVSO_CHART:
			return 14.5f;
		case CAT_PGC:
		case CAT_EXOPLANETARCHIVE:
			return 0.f;
		default:
			return 13.f;
	}
}

// returns the pretty name of a catalogue
const char *catalog_to_str(object_catalog cat) {
	switch (cat) {
		case CAT_TYCHO2:
			return _("Tycho-2");
		case CAT_NOMAD:
			return _("NOMAD");
		case CAT_GAIADR3:
			return _("Gaia DR3");
		case CAT_PPMXL:
			return _("PPMXL");
		case CAT_BSC:
			return _("bright stars");
		case CAT_APASS:
			return _("APASS");
		case CAT_VSX:
			return _("AAVSO Variable stars");
		case CAT_SIMBAD:
			return _("SIMBAD");
		case CAT_PGC:
			return _("PGC");
		case CAT_EXOPLANETARCHIVE:
			return _("Exoplanet archive");
		case CAT_IMCCE:
			return _("IMCCE solar system");
		case CAT_AAVSO_CHART:
			return _("AAVSO VSP Chart");
		case CAT_LOCAL:
			return _("local Tycho-2+NOMAD");
		case CAT_ASNET:
			return _("local astrometry.net");
		case CAT_AN_MESSIER:
			return "Messier";
		case CAT_AN_NGC:
			return "NGC";
		case CAT_AN_IC:
			return "IC";
		case CAT_AN_LDN:
			return "LDN";
		case CAT_AN_SH2:
			return "Sh2";
		case CAT_AN_STARS:
			return "stars";
		case CAT_AN_USER_DSO:
			return "user-DSO";
		case CAT_AN_USER_SSO:
			return "user-SSO";
		case CAT_AN_USER_TEMP:
			return "user-temp";
		default:
			return _("unknown");
	}
}

gboolean is_star_catalogue(object_catalog Catalog) {
	switch (Catalog) {
		case CAT_TYCHO2 ...	CAT_SIMBAD:
		case CAT_EXOPLANETARCHIVE:
		case CAT_AAVSO_CHART:
		case CAT_AN_STARS:
		case CAT_LOCAL:
		case CAT_AN_USER_SSO:
		case CAT_AN_USER_TEMP:
			return TRUE;
	default:
		return FALSE;
	}
}

static gboolean find_and_check_cat_columns(gchar **fields, int nbcols, object_catalog Catalog, int *indexes, uint32_t *collist) {
	if (!nbcols)
		return FALSE;

	uint32_t catspec = siril_catalog_columns(Catalog);
	uint32_t res = 0;
	for (int i = 0; i < nbcols; i++) {
		int val = get_column_index(fields[i]);
		if (val < 0) {
			siril_debug_print("Unknown column %s found in the catalog, ignoring\n", fields[i]);
			continue;
		} 
		indexes[i] = val;
		res |= (1 << val);
	}
	if (res == catspec) {
		siril_debug_print("Found same columns as in the catalog spec\n");
		*collist = res;
		return TRUE;
	}
	if ((res & catspec) == catspec) {
		siril_debug_print("Found more columns than in the catalog spec, keeping them\n");
		*collist = res;
		return TRUE;
	}
	siril_log_color_message(_("Did not find the minimal set of columns necessary for the catalog %s, aborting\n"), "red", catalog_to_str(Catalog));
	return FALSE;
}

void siril_catalog_free_item(cat_item *item) {
	if (!item)
		return;
	g_free(item->name);
	g_free(item->alias);
	g_free(item->type);
}

static void fill_cat_item(cat_item *item, const gchar *input, cat_fields index) {
	switch (index) {
		case CAT_FIELD_RA:
			item->ra = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_DEC:
			item->dec = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_PMRA:
			item->pmra = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_PMDEC:
			item->pmdec = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_MAG:
			item->mag = (float)g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_BMAG:
			item->bmag = (float)g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_E_MAG:
			item->e_mag = (float)g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_E_BMAG:
			item->e_bmag = (float)g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_NAME:
			item->name = g_shell_unquote(input, NULL); // TAP queries return quoted names
			break;
		case CAT_FIELD_DIAMETER:
			item->diameter = (float)g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_ALIAS:
			item->alias = g_strdup(input);
			break;
		case CAT_FIELD_DATEOBS:
			item->dateobs = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_SITELAT:
			item->sitelat = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_SITELON:
			item->sitelon = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_SITEELEV:
			item->siteelev = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_VRA:
			item->vra = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_VDEC:
			item->vdec = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_TYPE:
			item->type = g_strdup(input);
			break;
		case CAT_FIELD_UNDEF: // columns with unknown headers
		default:
			break;
	}
}

// copies all the data from an item to another item
void siril_catalogue_copy_item(cat_item *from, cat_item *to) {
	if (!from || !to) {
		siril_debug_print("no item to copy from or to\n");
		return;
	}
	memcpy(to, from, sizeof(cat_item));
	if (from->name)
		to->name = g_strdup(from->name);
	if (from->alias)
		to->alias = g_strdup(from->alias);
	if (from->type)
		to->type = g_strdup(from->type);
}

// frees the member cat_items of a catalogue
// This is useful to launch the same query again, updating simply the catalog center for instance
void siril_catalog_free_items(siril_catalogue *siril_cat) {
	if (!siril_cat || !siril_cat->cat_items)
		return;
	for (int i = 0; i < siril_cat->nbitems; i++)
		siril_catalog_free_item(&siril_cat->cat_items[i]);
	siril_cat->cat_items = NULL;
	siril_cat->nbitems = 0;
	siril_cat->nbincluded = -1;
	siril_cat->projected = CAT_PROJ_NONE;
}

// frees a siril_catalogue
void siril_catalog_free(siril_catalogue *siril_cat) {
	if (!siril_cat)
		return;
	siril_catalog_free_items(siril_cat);
	g_free(siril_cat->IAUcode);
	g_free(siril_cat->header);
	free(siril_cat);
}

void siril_catalog_reset_projection(siril_catalogue *siril_cat) {
	if (!siril_cat)
		return;
	for (int i = 0; i < siril_cat->nbitems; i++) {
		siril_cat->cat_items[i].x = 0.;
		siril_cat->cat_items[i].y = 0.;
		siril_cat->cat_items[i].included = TRUE;
	}
	siril_cat->nbincluded = siril_cat->nbitems;
	siril_cat->projected = CAT_PROJ_NONE;
}

// returns a siril_catalogue structure with center
siril_catalogue *siril_catalog_fill_from_fit(fits *fit, object_catalog cat, float limit_mag) {
#ifndef HAVE_WCSLIB
	return NULL;
#endif
	if (!fit) {
		return NULL;
	}
	if (!has_wcs(fit)) {
		siril_debug_print("This only works on plate solved images\n");
		return NULL;
	}
	siril_catalogue *siril_cat = calloc(1, sizeof(siril_catalogue));
	double ra, dec;
	center2wcs(fit, &ra, &dec);
	double resolution = get_wcs_image_resolution(fit) * 3600.;
	if (limit_mag == -1)
		limit_mag = siril_catalog_get_default_limit_mag(cat);
	// Preparing the catalogue query
	siril_cat->cattype = cat;
	siril_cat->columns = siril_catalog_columns(cat);
	siril_cat->center_ra = ra;
	siril_cat->center_dec = dec;
	siril_cat->radius = get_radius_deg(resolution, fit->rx, fit->ry) * 60.;
	siril_cat->limitmag = limit_mag;
	if (fit->date_obs) {
		siril_cat->dateobs = fit->date_obs;
	}
	return siril_cat;
}

/* This is the entry point to query the catalogues
 * It will call the necessary functions whether the query 
 * is for a local catalogue or an online one
 *
 * Then these raw catalogues can be used in different ways:
 * - The astrometric solver reads them and projects stars on a 2D image
 *   centered on the catalogue center and saves that to another file (the
 *   function project_catalog does that) called catalog.proj. Then the function
 *   read_projected_catalog below reads this intermediary projection catalogue
 *   and stores that in an array of psf_star objects. The obtained stars can be
 *   used for registration, but do not correspond to image coordinates.
 *
 * - The PCC reads them and projects stars on a plate-solved image using WCS
 *   and stores them in condensed form (pcc_star struct containing only
 *   x,y,b,v), done in the function project_catalog_with_WCS
 *
 * - Comparison star list creation needs equatorial coordinates and B-V
 *   magnitudes, projection is also used but only to check if a star is inside
 *   the image and far from its borders. The required object is psf_star.
 */

int siril_catalog_conesearch(siril_catalogue *siril_cat) {
	int nbstars = 0;
	if (siril_cat->cattype < CAT_AN_MESSIER) // online
#ifndef HAVE_NETWORKING
	siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
	return 0;
#endif
		nbstars = siril_catalog_get_stars_from_online_catalogues(siril_cat);
	else if (siril_cat->cattype == CAT_LOCAL)
		nbstars = siril_catalog_get_stars_from_local_catalogues(siril_cat);
	else
		siril_debug_print("trying to conesearch an invalid catalog type");
	return nbstars;
}

// This is the generic parser for all csv catalogues used by Siril
// (annotation, downloaded, nina lists...)
// Loads the csv file and fills the cat_items members of the catalogue in entry
int siril_catalog_load_from_file(siril_catalogue *siril_cat, const gchar *filename) {
	GError *error = NULL;
	GFile *catalog_file = g_file_new_for_path(filename);
	GInputStream *input_stream = (GInputStream*) g_file_read(catalog_file, NULL, &error);
	if (!input_stream) {
		if (error != NULL) {
			siril_log_message(_("Could not load the star catalog (%s)."), error->message);
			g_clear_error(&error);
		} else
			siril_log_message(_("Could not load the star catalog (%s)."), "generic error");
		return 1;
	}

	int nb_alloc = 1200, nb_items = 0;
	cat_item *cat_items = calloc(nb_alloc, sizeof(cat_item));
	if (!cat_items) {
		PRINT_ALLOC_ERR;
		g_object_unref(input_stream);
		siril_catalog_free(siril_cat);
		return 1;
	}

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	gchar *line;
	gboolean header_read = FALSE, has_error = FALSE;
	int *indexes = NULL;
	int nbcols;
	GString *header = NULL;
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		if (line[0] == '\0') {
			g_free(line);
			continue;
		}
		remove_trailing_eol(line);
		if (line[0] == COMMENT_CHAR) { // storing header lines
			if (!header)
				header = g_string_new(line);
			else
				g_string_append_printf(header, "\n%s", line);
			g_free(line);
			continue;
		}
		if (!header_read) { // reading the first line which holds the columns names
			gchar **fields = g_strsplit(line, ",", -1);
			nbcols = g_strv_length(fields);
			if (!nbcols) {
				siril_log_color_message(_("No columns found in the catalogue\n"), "red");
				has_error = TRUE;
			} else {
				indexes = malloc(nbcols * sizeof(int));
				for (int i = 0; i < nbcols; i++) {
					indexes[i] = -1;
				}
				has_error = !find_and_check_cat_columns(fields, nbcols, siril_cat->cattype, indexes, &siril_cat->columns);
			}
			g_strfreev(fields);
			if (has_error) {
				g_object_unref(data_input);
				g_object_unref(input_stream);
				siril_catalog_free(siril_cat);
				g_free(line);
				return 1;
			}
			header_read = TRUE;
			g_free(line);
			continue;
		}
		gchar **vals = g_strsplit(line, ",", -1);
		int size = g_strv_length(vals);
		if (size != nbcols) { // checking that current line has a number of columns consistent with the headers
			siril_log_color_message(_("Malformed line found %s\n"), "red", line);
			g_object_unref(data_input);
			g_object_unref(input_stream);
			siril_catalog_free(siril_cat);
			g_free(line);
			g_strfreev(vals);
			return 1;
		}
		if (nb_items >= nb_alloc) { // re-allocating if there is more to read
			nb_alloc *= 2;
			cat_item *new_array = realloc(cat_items, nb_alloc * sizeof(cat_item));
			if (!new_array) {
				PRINT_ALLOC_ERR;
				g_object_unref(data_input);
				g_object_unref(input_stream);
				siril_catalog_free(siril_cat);
				g_free(line);
				g_strfreev(vals);
				return 1;
			}
			cat_items = new_array;
		}
		memset(&cat_items[nb_items], 0, sizeof(cat_item));
		for (int i = 0; i < nbcols; i++) {
			fill_cat_item(&cat_items[nb_items], vals[i], indexes[i]);
		}
		// magnitudes above 30 are a code for 'undefined'
		if (siril_cat->phot && (cat_items[nb_items].mag == 0. || cat_items[nb_items].mag > 30. ||
					 cat_items[nb_items].bmag == 0. || cat_items[nb_items].bmag > 30.)) {
			// we reset the values and skip incrementing
			cat_items[nb_items].mag = 0.;
			cat_items[nb_items].bmag = 0.;
		} else {
			cat_items[nb_items].included = TRUE;
			nb_items++;
		}
		g_free(line);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);
	if (nb_items == 0) {
		free(cat_items);
		siril_log_color_message(_("Catalog %s was read but no items were found, nothing to show\n"), "salmon", filename);
		return 1;
	}
	cat_item *final_array = realloc(cat_items, nb_items * sizeof(cat_item));
	siril_cat->cat_items = final_array;
	siril_cat->nbitems = nb_items;
	siril_cat->nbincluded = nb_items;
	if (header)
		siril_cat->header = g_string_free(header, FALSE);
	siril_debug_print("read %d%s items from catalogue\n", nb_items, siril_cat->phot ? " photometric" : "");
	return 0;
}

// Writes the catalogue to the given filepath
// Optionally writes the header if provided
gboolean siril_catalog_write_to_file(siril_catalogue *siril_cat, const gchar *filename, gchar *header) {
	if (!siril_cat || !filename)
		return FALSE;

	/* First we test if root directory already exists */
	gchar *root = g_path_get_dirname(filename);
	if (!g_file_test(root, G_FILE_TEST_IS_DIR)) {
		if (g_mkdir_with_parents(root, 0755) < 0) {
			siril_log_color_message(_("Cannot create output folder: %s\n"), "red", root);
			g_free(root);
			return FALSE;
		}
	}

	// Then we delete it if it exists
	GError *error = NULL;
	GFile *file = g_file_new_for_path(filename);
	if (g_file_test(filename, G_FILE_TEST_EXISTS)) { // file already exists, we removed it
		if (!g_file_delete(file, NULL, &error)) {
			siril_log_color_message(_("File % cannot be deleted (%s), aborting\n"), "red", filename, (error) ? error->message : "unknown error");
			g_object_unref(file);
			return FALSE;
		}
	}

	// and we open the stream to write to
	GOutputStream *output_stream = (GOutputStream*) g_file_create(file, G_FILE_CREATE_NONE, NULL, &error);
	if (!output_stream) {
		siril_log_color_message(_("Cannot create catalogue file %s (%s)\n"), "red", filename, (error) ? error->message : "unknown error");
		g_clear_error(&error);
		g_object_unref(file);
		return FALSE;
	}
	gsize n;
	if (header)
		g_output_stream_printf(output_stream, &n, NULL, NULL, "%s\n", header);
	
	// we write the line containing the columns names based on catalog spec
	GString *columns = NULL;
	int nbcols = 0;
	int index[MAX_CAT_COLUMNS];
	for (int i = 0; i < MAX_CAT_COLUMNS; i++) {
		if (siril_cat->columns & (1 << i)) {
			if (!columns)
				columns = g_string_new(cat_columns[i]);
			else
				g_string_append_printf(columns, ",%s", cat_columns[i]);
			index[nbcols++] = i;
		}
	}
	gchar *columns_str = g_string_free(columns, FALSE);
	if (columns_str) {
		g_output_stream_printf(output_stream, &n, NULL, NULL, "%s", columns_str);
	} else {
		siril_debug_print("no columns to write");
		g_object_unref(file);
		g_object_unref(output_stream);
		return FALSE;
	}
	g_free(columns_str);
	for (int j = 0; j < siril_cat->nbitems; j++) {
		gchar **tokens = calloc(nbcols + 1, sizeof(gchar *));
		for (int i = 0; i < nbcols; i++) {
			tokens[i] = get_field_to_str(siril_cat->cat_items[j], index[i]);
		}
		gchar *newline = g_strjoinv(",", tokens);
		g_output_stream_printf(output_stream, &n, NULL, NULL, "\n%s", newline);
		g_free(newline);
	}
	g_object_unref(output_stream);
	return TRUE;
}

// appends an item at the end of a siril_catalogue
gboolean siril_catalog_append_item(siril_catalogue *siril_cat, cat_item *item) {
	if (!siril_cat || !item)
		return FALSE;
	cat_item *new_array = realloc(siril_cat->cat_items, (siril_cat->nbitems + 1) * sizeof(cat_item));
	if (!new_array) {
		PRINT_ALLOC_ERR;
		return FALSE;
	}
	siril_cat->cat_items = new_array;
	siril_catalogue_copy_item(item, &siril_cat->cat_items[siril_cat->nbitems]);
	siril_cat->nbitems++;
	// we can't have a catalogue only partially projected so we reset its projection if any
	if (siril_cat->projected > CAT_PROJ_NONE) 
		siril_catalog_reset_projection(siril_cat);
	return TRUE;
}

gboolean can_use_proper_motion(fits *fit, siril_catalogue *siril_cat) {
	if (!fit || !siril_cat)
		return FALSE;
	if (!fit->date_obs) {
		// siril_debug_print("This image does not have any DATE-OBS information, cannot account for stars proper motions\n");
		return FALSE;
	}
	if (!(siril_cat->columns & (1 << CAT_FIELD_PMRA)) || !(siril_cat->columns & (1 << CAT_FIELD_PMDEC))) {
		// siril_debug_print("This catalog does not have proper motion info, will not be computed\n");
		return FALSE;
	}
	siril_debug_print("Catalogue %s will account for proper motions\n", catalog_to_str(siril_cat->cattype));
	return TRUE;
}

gboolean can_use_velocity(fits *fit, siril_catalogue *siril_cat) {
	if (!fit || !siril_cat)
		return FALSE;
	if (!fit->date_obs) {
		// siril_debug_print("This image does not have any DATE-OBS information, cannot account for velocity\n");
		return FALSE;
	}
	if (!(siril_cat->columns & (1 << CAT_FIELD_VRA)) || !(siril_cat->columns & (1 << CAT_FIELD_VDEC))) {
		// siril_debug_print("This catalog does not have velocity info, will not be computed\n");
		return FALSE;
	}
	siril_debug_print("Catalogue %s will account for velocities\n", catalog_to_str(siril_cat->cattype));
	return TRUE;
}

// projects passed catalogue using the wcs data contained in the fit
// corrects for proper motions if the flag is TRUE and the necessary data is included
// in the fit (dateobs) and in the catalogue (pmra and pmdec fields)
// corrects for object velocity if the flag is true and if necessary data is included
// in the catalogue (vra and vdec fields)
int siril_catalog_project_with_WCS(siril_catalogue *siril_cat, fits *fit, gboolean use_proper_motion, gboolean use_velocity) {
#ifndef HAVE_WCSLIB
	return 1
#endif
	if (!(siril_cat->columns & (1 << CAT_FIELD_RA)) || !(siril_cat->columns & (1 << CAT_FIELD_DEC))) {
		siril_debug_print("catalogue %s does not have the necessary columns\n");
		return 1;
	}
	int nbincluded = 0;
	double jyears = 0.;
	double *world = NULL, *x = NULL, *y = NULL;
	int *status = NULL;
	double tobs = 0.;
	if (use_proper_motion && can_use_proper_motion(fit, siril_cat)) {
		GDateTime *dt = g_date_time_ref(fit->date_obs);
		gdouble jd = date_time_to_Julian(dt);
		g_date_time_unref(dt);
		double J2000 = 2451545.0;
		jyears = (jd - J2000) / 365.25;
	}
	if (use_velocity && can_use_velocity(fit, siril_cat)) {
		GDateTime *dt = g_date_time_ref(fit->date_obs);
		tobs = date_time_to_Julian(dt);
		g_date_time_unref(dt);
	}
	world = malloc( 2 * siril_cat->nbitems * sizeof(double));
	x = malloc(siril_cat->nbitems * sizeof(double));
	y = malloc(siril_cat->nbitems * sizeof(double));
	if (!world || !x || !y) {
		PRINT_ALLOC_ERR;
		goto clean_and_exit;
	}
	int ind = 0;
	for (int i = 0; i < siril_cat->nbitems; i++) {
		double ra = siril_cat->cat_items[i].ra;
		double dec = siril_cat->cat_items[i].dec;
		double decrad = dec * DEGTORAD;
		if (use_proper_motion) {
			ra += siril_cat->cat_items[i].pmra / cos(decrad) * jyears * 2.77777778e-7;
			dec += siril_cat->cat_items[i].pmdec * jyears * 2.77777778e-7;
		}
		if (use_velocity) {
			double deltahours = (tobs - siril_cat->cat_items[i].dateobs) * 24.;
			ra += siril_cat->cat_items[i].vra / cos(decrad) * deltahours * 2.77777778e-4;
			dec += siril_cat->cat_items[i].vdec * deltahours * 2.77777778e-4;

		}
		world[ind++] = ra;
		world[ind++] = dec;
	}
	status = wcs2pix_array(fit, siril_cat->nbitems, world, x, y);
	if (!status)
		goto clean_and_exit;
	for (int i = 0; i < siril_cat->nbitems; i++) {
		if (!status[i]) {
			siril_cat->cat_items[i].x = x[i];
			siril_cat->cat_items[i].y = y[i];
			siril_cat->cat_items[i].included = TRUE;
			nbincluded++;
		} else {
			siril_cat->cat_items[i].included = FALSE;
		}
	}
clean_and_exit:
	siril_cat->nbincluded = nbincluded;
	free(world);
	free(x);
	free(y);
	siril_cat->projected = CAT_PROJ_WCS;
	return !(nbincluded > 0);
} 

// projects passed catalogue wrt to the center ra0 and dec0 coordinates
// corrects for proper motions if the flag is TRUE and the necessary data is passed
// (dateobs) and found in the catalogue (pmra and pmdec fields)
int siril_catalog_project_at_center(siril_catalogue *siril_cat, double ra0, double dec0, gboolean use_proper_motion, GDateTime *date_obs) {
	if (!(siril_cat->columns & (1 << CAT_FIELD_RA)) || !(siril_cat->columns & (1 << CAT_FIELD_DEC)))
		return 1;
	double jyears = 0.;
	if (use_proper_motion) {
		if (!date_obs) {
			siril_log_color_message(_("no DATE-OBS information, cannot account for stars proper motions\n"), "salmon");
			use_proper_motion = FALSE;
		} else if (!(siril_cat->columns & (1 << CAT_FIELD_PMRA)) || !(siril_cat->columns & (1 << CAT_FIELD_PMDEC))) {
			siril_log_color_message(_("This catalog does not have proper motion info, will not be computed\n"), "salmon");
			use_proper_motion = FALSE;
		} else {
			GDateTime *dt = g_date_time_ref(date_obs);
			gdouble jd = date_time_to_Julian(dt);
			g_date_time_unref(dt);
			double J2000 = 2451545.0;
			jyears = (jd - J2000) / 365.25;
		}
	}
	dec0 *= DEGTORAD;
	double sindec0 = sin(dec0);
	double cosdec0 = cos(dec0);
	for (int i = 0; i < siril_cat->nbitems; i++) {
		double ra = siril_cat->cat_items[i].ra;
		double dec = siril_cat->cat_items[i].dec;
		double decrad = dec * DEGTORAD;
		if (use_proper_motion) {
			ra += siril_cat->cat_items[i].pmra / cos(decrad) * jyears * 2.77777778e-7;
			dec += siril_cat->cat_items[i].pmdec * jyears * 2.77777778e-7;
		}
		dec *= DEGTORAD;
		double delta_ra = (ra - ra0) * DEGTORAD;
		double xx = cos(dec) * sin(delta_ra);
		double yy = sindec0 * sin(dec) + cosdec0 * cos(dec) * cos(delta_ra);
		double xi = (xx / yy);
		xx = cosdec0 * sin(dec) - sindec0 * cos(dec) * cos(delta_ra);
		double eta = (xx / yy);
		siril_cat->cat_items[i].x = xi * RADtoASEC;
		siril_cat->cat_items[i].y = eta * RADtoASEC;
		siril_cat->cat_items[i].included = TRUE;
	}
	siril_cat->nbincluded = siril_cat->nbitems;
	siril_cat->projected = CAT_PROJ_PLATE;
	return 0;
}

// TODO: move to a file for callbacks and remove gtk include
static gboolean end_conesearch(gpointer p) {
	siril_catalogue *temp_cat = (siril_catalogue *) p;
	if (temp_cat) {
		purge_user_catalogue(CAT_AN_USER_TEMP);
		if (!load_siril_cat_to_temp(temp_cat)) {
			GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button"));
			refresh_found_objects();
			if (!gtk_toggle_tool_button_get_active(button)) {
				gtk_toggle_tool_button_set_active(button, TRUE);
			} else {
				refresh_found_objects();
				redraw(REDRAW_OVERLAY);
			}
		}
		// siril_catalog_free(temp_cat);
	}
	return end_generic(NULL);
}

// worker for the command conesearch
gpointer conesearch_worker(gpointer p) {
	siril_catalogue *siril_cat = (siril_catalogue *) p;
	siril_catalogue *temp_cat = NULL;
	int retval = -1;
	if (!siril_cat) {
		siril_debug_print("no query passed");
		goto exit_conesearch;
	}

	// launching the query
	if (!siril_catalog_conesearch(siril_cat)) {// returns the nb of stars
		goto exit_conesearch;
	}
	if (siril_cat->cattype != CAT_LOCAL)
		siril_log_message(_("The %s catalog has been successfully downloaded.\n"), catalog_to_str(siril_cat->cattype));

	/* project using WCS */
	gboolean use_proper_motion = can_use_proper_motion(&gfit, siril_cat);
	gboolean use_velocity = can_use_velocity(&gfit, siril_cat);
	if (siril_catalog_project_with_WCS(siril_cat, &gfit, use_proper_motion, use_velocity)) { // TODO: pass *fit instead of gfit
		goto exit_conesearch;
	}
	int nb_stars = siril_cat->nbitems;
	sort_cat_items_by_mag(siril_cat);
	int j = 0;
	// preparing the output catalog
	temp_cat = calloc(1, sizeof(siril_catalogue));
	temp_cat->cattype = CAT_AN_USER_TEMP;
	temp_cat->columns = siril_catalog_columns(CAT_AN_USER_TEMP);
	for (int i = 0; i < nb_stars; i++) {
		if (!siril_cat->cat_items[i].included || siril_cat->cat_items[i].mag > siril_cat->limitmag)
			continue;
		if (!com.script && !temp_cat->cat_items)
			temp_cat->cat_items = calloc(siril_cat->nbincluded, sizeof(cat_item));
		if (!com.script) {
			temp_cat->cat_items[j].ra = siril_cat->cat_items[i].ra;
			temp_cat->cat_items[j].dec = siril_cat->cat_items[i].dec;
			if (siril_cat->cat_items[i].name)
				temp_cat->cat_items[j].name = g_strdup(siril_cat->cat_items[i].name);
		}
		// some catalogues display the list of found objects
		if (siril_cat->cattype == CAT_IMCCE) // classes are defined at https://vo.imcce.fr/webservices/skybot/?documentation#field_1
			siril_log_message("%s (%s) - mag:%3.1f\n", siril_cat->cat_items[i].name, siril_cat->cat_items[i].type, siril_cat->cat_items[i].mag);
		if (siril_cat->cattype == CAT_AAVSO_CHART) // https://www.aavso.org/api-vsp
			siril_log_message("%s - V:%3.1f [%5.3f]- B:%3.1f [%5.3f] - RA:%g - DEC: %g\n",
			siril_cat->cat_items[i].name, 
			siril_cat->cat_items[i].mag, 
			siril_cat->cat_items[i].e_mag, 
			siril_cat->cat_items[i].bmag,
			siril_cat->cat_items[i].e_bmag, 
			siril_cat->cat_items[i].ra,
			siril_cat->cat_items[i].dec);
		if (siril_cat->cattype == CAT_EXOPLANETARCHIVE)
			siril_log_message("%s - mag:%3.1f\n", siril_cat->cat_items[i].name, siril_cat->cat_items[i].mag);
		j++;
	}
	siril_log_message("%d objects found%s in the image (mag limit %.2f)\n", j,
			siril_cat->phot ? " with valid photometry data" : "", siril_cat->limitmag);
	retval = 0;
exit_conesearch:
	siril_catalog_free(siril_cat);
	if (!com.script) {
		if (temp_cat)
			temp_cat->nbitems = j;
		siril_add_idle(end_conesearch, temp_cat); // temp_cat will be freed in the idle
	} else {
		end_generic(NULL);
	}
	return GINT_TO_POINTER(retval);
}

// Haversine formula on unit sphere
// https://en.wikipedia.org/wiki/Haversine_formula
// dec is phi, ra is lambda
// in degrees
double compute_coords_distance(double ra1, double dec1, double ra2, double dec2) {
	double dec1_r = dec1 * DEGTORAD, dec2_r = dec2 * DEGTORAD;
	double dra_2 = 0.5 * (ra2 - ra1) * DEGTORAD;
	double ddec_2 = 0.5 * (dec2_r - dec1_r);
	double sin_ddec = sin(ddec_2), sin_dra = sin(dra_2);
	double h = sin_ddec * sin_ddec + cos(dec1_r) * cos(dec2_r) * sin_dra * sin_dra;
	if (h > 1.)
		return 180.0;   // h = 1, asin(1) is pi/2
	return 2.0 * asin(sqrt(h)) * RADTODEG;
}
// TODO: using this for the moment to avoid chaging too many files
// This copies the info contained in the catalogue to a psf_star** list
// only the included items are copied over
psf_star **convert_siril_cat_to_psf_stars(siril_catalogue *siril_cat, int *nbstars) {
	*nbstars = 0;
	if (!siril_cat)
		return NULL;
	if (siril_cat->projected == CAT_PROJ_NONE) {
		siril_debug_print("Catalog has not been projected\n");
	}
	if (!(siril_cat->columns & (1 << CAT_FIELD_RA)) || !(siril_cat->columns & (1 << CAT_FIELD_DEC)) || !(siril_cat->columns & (1 << CAT_FIELD_MAG)))
		return NULL;
	psf_star **results = new_fitted_stars(siril_cat->nbincluded);

	int n = 0;
	for (int i = 0; i < siril_cat->nbitems; i++) {
		if (n > siril_cat->nbincluded) {
			siril_debug_print("problem when converting siril_cat to psf_stars, more than allocated");
		}
		if (siril_cat->cat_items[i].included) {
			results[n] = new_psf_star();
			results[n]->xpos = siril_cat->cat_items[i].x;
			results[n]->ypos = siril_cat->cat_items[i].y;
			results[n]->mag = siril_cat->cat_items[i].mag;
			results[n]->Bmag = siril_cat->cat_items[i].bmag;
			results[n]->ra = siril_cat->cat_items[i].ra;
			results[n]->dec = siril_cat->cat_items[i].dec;
			if (siril_cat->cat_items[i].name)
				results[n]->star_name = g_strdup(siril_cat->cat_items[i].name);
			n++;
		}
	}
	results[n] = NULL;
	if (n != siril_cat->nbincluded) {
		siril_debug_print("problem when converting siril_cat to psf_stars, number differs from catalogue info");
		free_fitted_stars(results);
		return NULL;
	}
	*nbstars = n;
	return results;
}

sky_object_query_args *init_sky_object_query() {
	sky_object_query_args *new_query = calloc(1, sizeof(sky_object_query_args));
	if (!new_query) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	new_query->server = -1;
	return new_query;
}

void free_sky_object_query(sky_object_query_args *args) {
	if (!args)
		return;
	g_free(args->name);
	g_free(args->prefix);
	siril_catalog_free_item(args->item);
	free(args);
}

