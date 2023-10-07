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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "algos/search_objects.h"
#include "algos/siril_wcs.h"
#include "algos/annotate.h"
#include "algos/astrometry_solver.h"
#include "algos/comparison_stars.h"
#include "io/remote_catalogues.h"
#include "io/local_catalogues.h"
#include "registration/matching/misc.h"

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
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME)| (1 << CAT_FIELD_VRA) | (1 << CAT_FIELD_VDEC) | (1 << CAT_FIELD_TYPE);
		case CAT_AN_MESSIER:
		case CAT_AN_NGC:
		case CAT_AN_IC:
		case CAT_AN_LDN:
		case CAT_AN_SH2:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_ALIAS);
		case CAT_AN_STARS:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME);
		case CAT_AN_USER_DSO:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME);
		case CAT_AN_USER_SSO:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_DATEOBS) | (1 << CAT_FIELD_SITELAT) | (1 << CAT_FIELD_SITELON) | (1 << CAT_FIELD_SITEELEV);
		case CAT_COMPSTARS:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_TYPE);
		case CAT_LOCAL:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG);
		case CAT_AN_USER_TEMP:
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

//TODO: to be rewritten
/* There are several ways of obtaining star catalogue data in Siril.
 * The raw catalogue contain in general star name, RA and dec coords, V and B
 * magnitudes. The download_catalog function makes a request to an online
 * source with a centre, a radius and a limit magnitude and stores that in a
 * raw cache file named 'cat-type-ra-dec-radius-mag.cat'.
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
		default:
			return _("unknown");
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
	siril_log_color_message(_("Did not find the minimal set of columns necessary for this catalog, aborting\n"), "red");
	return FALSE;
}

// not used for now
// static int siril_cat_init_from_filename(siril_catalogue *siril_cat, const gchar *filename) {
// 	if (!siril_cat) {
// 		PRINT_ALLOC_ERR;
// 		return 1;
// 	}
// 	if (!filename)
// 		return 1;

// 	// remove the extension
// 	gchar *basename = g_path_get_basename(filename);
// 	GString *name = g_string_new(basename);
// 	g_string_replace(name, ".csv", "", 0);
// 	gchar *namewoext = g_string_free_and_steal(name);

// 	// if a filename is passed, we will parse its name to fill the structure
// 	if (g_str_has_prefix(namewoext, "cat")) {
// 		gchar **fields = g_strsplit(namewoext, "_", -1);
// 		int n = g_strv_length(fields);
// 		if (n < 6 && n > 7) {
// 			siril_log_color_message(_("Could not parse the catalogue name %s, aborting\n"), "red", namewoext);
// 			g_strfreev(fields);
// 			return NULL;
// 		}
// 		siril_cat->cattype = (int)g_ascii_strtoll(fields[1], NULL, 10);
// 		siril_cat->center_ra = g_ascii_strtod(fields[2], NULL);
// 		siril_cat->center_dec = g_ascii_strtod(fields[3], NULL);
// 		siril_cat->radius = g_ascii_strtod(fields[4], NULL);
// 		if (n < 7) {
// 			siril_cat->limitmag = g_ascii_strtod(fields[5], NULL);
// 		} else {
// 			GDateTime *dt = FITS_date_to_date_time(fields[5]);
// 			siril_cat->dateobs = dt;
// 			siril_cat->IAUcode = g_strdup(fields[6]);
// 		}
// 		g_strfreev(fields);
// 		return 0;
// 	}
// 	return 1;
// }

static void siril_catalog_free_item(cat_item *item) {
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
			item->name = g_strdup(input);
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

void siril_catalog_free_items(siril_catalogue *siril_cat) {
	if (!siril_cat || !siril_cat->cat_items)
		return;
	for (int i = 0; i < siril_cat->nbitems; i++)
		siril_catalog_free_item(&siril_cat->cat_items[i]);
	siril_cat->cat_items = NULL;
}

void siril_catalog_free(siril_catalogue *siril_cat) {
	if (!siril_cat)
		return;
	siril_catalog_free_items(siril_cat);
	g_free(siril_cat->IAUcode);
	free(siril_cat);
}

int siril_catalog_conesearch(siril_catalogue *siril_cat) {
	int nbstars = 0;
	if (siril_cat->cattype < CAT_AN_MESSIER) // online
		nbstars = siril_catalog_get_stars_from_online_catalogues(siril_cat);
	else if (siril_cat->cattype == CAT_LOCAL)
		nbstars = siril_catalog_get_stars_from_local_catalogues(siril_cat);
	else
		siril_debug_print("trying to conesearch an invalid catalog type");
	return nbstars;
}

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
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		if (line[0] == COMMENT_CHAR) { // skipping comments
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
			nb_items++;
		}
		g_free(line);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);
	if (nb_items == 0) {
		free(cat_items);
		siril_log_color_message(_("Catalog %s was read but no items were found\n"), "red", filename);
		return 1;
	}
	cat_item *final_array = realloc(cat_items, nb_items * sizeof(cat_item));
	siril_cat->cat_items = final_array;
	siril_cat->nbitems = nb_items;
	siril_debug_print("read %d%s items from catalogue\n", nb_items, siril_cat->phot ? " photometric" : "");
	return 0;
}

int siril_catalog_project_with_WCS(siril_catalogue *siril_cat, fits *fit, gboolean use_proper_motion) {
#ifndef HAVE_WCSLIB
	return 1
#endif
	if (!(siril_cat->columns & (1 << CAT_FIELD_RA)) || !(siril_cat->columns & (1 << CAT_FIELD_DEC)) || !(siril_cat->columns & (1 << CAT_FIELD_MAG))) {
		siril_debug_print("catalogue %s does not have the necessary columns\n");
		return 1;
	}
	int nbincluded = 0;
	double jyears = 0.;
	if (use_proper_motion) {
		if (!fit->date_obs) {
			siril_log_color_message(_("This image does not have any DATE-OBS information, cannot account for stars proper motions\n"), "salmon");
			use_proper_motion = FALSE;
		} else if (!(siril_cat->columns & (1 << CAT_FIELD_PMRA)) || !(siril_cat->columns & (1 << CAT_FIELD_PMDEC))) {
			siril_log_color_message(_("This catalog does not have proper motion info, will not be computed\n"), "salmon");
			use_proper_motion = FALSE;
		} else {
			GDateTime *dt = g_date_time_ref(fit->date_obs);
			gdouble jd = date_time_to_Julian(dt);
			g_date_time_unref(dt);
			double J2000 = 2451545.0;
			jyears = (jd - J2000) / 365.25;
		}
	}
	double x, y;
	for (int i = 0; i < siril_cat->nbitems; i++) {
		double ra = siril_cat->cat_items[i].ra;
		double dec = siril_cat->cat_items[i].dec;
		if (use_proper_motion) {
			ra += siril_cat->cat_items[i].pmra / cos(dec) * jyears * 2.77777778e-7;
			dec += siril_cat->cat_items[i].pmdec * jyears * 2.77777778e-7;
		}
		if (!wcs2pix(fit, ra, dec, &x, &y)) {
			siril_cat->cat_items[i].x = x;
			siril_cat->cat_items[i].y = y;
			siril_cat->cat_items[i].included = TRUE;
			nbincluded++;
		}
	}
	siril_cat->nbincluded = nbincluded;
	return 0;
} 

int siril_catalog_project_at_center(siril_catalogue *siril_cat, double ra0, double dec0, gboolean use_proper_motion, GDateTime *date_obs) {
	if (!(siril_cat->columns & (1 << CAT_FIELD_RA)) || !(siril_cat->columns & (1 << CAT_FIELD_DEC)) || !(siril_cat->columns & (1 << CAT_FIELD_MAG)))
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
	ra0 *= DEGTORAD;
	double sindec0 = sin(dec0);
	double cosdec0 = cos(dec0);
	for (int i = 0; i < siril_cat->nbitems; i++) {
		double ra = siril_cat->cat_items[i].ra;
		double dec = siril_cat->cat_items[i].dec * DEGTORAD;
		if (use_proper_motion) {
			ra += siril_cat->cat_items[i].pmra / cos(dec) * jyears * 2.77777778e-7;
			dec += siril_cat->cat_items[i].pmdec * jyears * 2.77777778e-7;
		}
		double delta_ra = ra * DEGTORAD - ra0;
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
	return 0;
}

// TODO: using this for the moment to avoid chaging too many files
psf_star **convert_siril_cat_to_psf_stars(siril_catalogue *siril_cat, int *nbstars) {
	*nbstars = 0;
	if (!siril_cat)
		return NULL;
	if (siril_cat->projected == CAT_PROJ_NONE) {
		siril_debug_print("Catalog has not been projected\n");
		return NULL;
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

