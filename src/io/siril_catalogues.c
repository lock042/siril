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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtk/gtk.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "core/siril_date.h"
#include "core/processing.h"
#include "core/command_line_processor.h"
#include "algos/PSF.h"
#include "algos/search_objects.h"
#include "algos/siril_wcs.h"
#include "algos/photometric_cc.h"
#include "io/annotation_catalogues.h"
#include "algos/astrometry_solver.h"
#include "algos/comparison_stars.h"
#include "io/remote_catalogues.h"
#include "io/local_catalogues.h"
#include "registration/matching/misc.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"
#include "gui/utils.h"
#include "gui/siril_plot.h"


static void free_conesearch_params(conesearch_params *params);
static void free_conesearch_args(conesearch_args *args);

// This list defines the columns that can possibly be found in any catalogue
const gchar *cat_columns[] = {
	[CAT_FIELD_RA] = "ra",
	[CAT_FIELD_DEC] = "dec",
	[CAT_FIELD_RA1] = "ra1",
	[CAT_FIELD_DEC1] = "dec1",
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
	[CAT_FIELD_TYPE] = "type",
	[CAT_FIELD_TEFF] = "teff",
	[CAT_FIELD_XPSAMP] = "xpsamp",
	[CAT_FIELD_GAIASOURCEID] = "source_id"
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

static gchar *get_field_to_str(cat_item *item, cat_fields field) {
	switch (field) {
		case CAT_FIELD_RA:
			return (item->ra) ? g_strdup_printf("%.6f", item->ra) : "";
		case CAT_FIELD_DEC:
			return (item->dec) ? g_strdup_printf("%.6f", item->dec) : "";
		case CAT_FIELD_RA1:
			return (item->ra1) ? g_strdup_printf("%.6f", item->ra1) : "";
		case CAT_FIELD_DEC1:
			return (item->dec1) ? g_strdup_printf("%.6f", item->dec1) : "";
		case CAT_FIELD_PMRA:
			return (item->pmra) ? g_strdup_printf("%g", item->pmra) : "";
		case CAT_FIELD_PMDEC:
			return (item->pmdec) ? g_strdup_printf("%g", item->pmdec) : "";
		case CAT_FIELD_MAG:
			return (item->mag) ? g_strdup_printf("%g", item->mag) : "";
		case CAT_FIELD_BMAG:
			return (item->bmag) ? g_strdup_printf("%g", item->bmag) : "";
		case CAT_FIELD_E_MAG:
			return (item->e_mag) ? g_strdup_printf("%g", item->e_mag) : "";
		case CAT_FIELD_E_BMAG:
			return (item->e_bmag) ? g_strdup_printf("%g", item->e_bmag) : "";
		case CAT_FIELD_DIAMETER:
			return (item->diameter) ? g_strdup_printf("%g", item->diameter) : "";
		case CAT_FIELD_DATEOBS:
			return (item->dateobs) ? g_strdup_printf("%.12f", item->dateobs) : "";
		case CAT_FIELD_SITELAT:
			return (item->sitelat) ? g_strdup_printf("%g", item->sitelat) : "";
		case CAT_FIELD_SITELON:
			return (item->sitelon) ? g_strdup_printf("%g", item->sitelon) : "";
		case CAT_FIELD_SITEELEV:
			return (item->siteelev) ? g_strdup_printf("%g", item->siteelev) : "";
		case CAT_FIELD_VRA:
			return (item->vra) ? g_strdup_printf("%g", item->vra) : "";
		case CAT_FIELD_VDEC:
			return (item->vdec) ? g_strdup_printf("%g", item->vdec) : "";
		case CAT_FIELD_NAME:
			return (item->name) ? g_strdup(item->name) : "";
		case CAT_FIELD_ALIAS:
			return (item->alias) ? g_strdup(item->alias) : "";
		case CAT_FIELD_TYPE:
			return (item->type) ? g_strdup(item->type) : "";
		case CAT_FIELD_TEFF:
			return (item->teff) ? g_strdup_printf("%.6f", item->teff) : "";
		case CAT_FIELD_GAIASOURCEID:
			return (item->gaiasourceid) ? g_strdup_printf("%" G_GUINT64_FORMAT, item->gaiasourceid) : "";
		default:
			return NULL;
	}
}
// This function defines the fields that should be present in each catalog
// To be used as sanity check
uint32_t siril_catalog_columns(siril_cat_index cat) {
	switch (cat) {
		case CAT_TYCHO2:
		case CAT_NOMAD:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG);
		case CAT_GAIADR3:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_TEFF);
		case CAT_GAIADR3_DIRECT:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_TEFF) |
			(1 << CAT_FIELD_GAIASOURCEID);
		case CAT_PPMXL:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG);
		case CAT_BSC:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME);
		case CAT_APASS:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_E_MAG) | (1 << CAT_FIELD_E_BMAG);
		case CAT_GCVS:
		case CAT_VSX:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME);
		case CAT_SIMBAD:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_NAME);
		case CAT_VARISUM:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME);
		case CAT_PGC:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_DIAMETER);
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
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_ALIAS);
		case CAT_AN_CONST:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_RA1) | (1 << CAT_FIELD_DEC1);
		case CAT_AN_CONST_NAME:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_ALIAS);
		case CAT_AN_USER_DSO:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_ALIAS);
		case CAT_AN_USER_SSO:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_ALIAS) | (1 << CAT_FIELD_DATEOBS) | (1 << CAT_FIELD_SITELAT) | (1 << CAT_FIELD_SITELON) | (1 << CAT_FIELD_SITEELEV) | (1 << CAT_FIELD_VRA) | (1 << CAT_FIELD_VDEC) | (1 << CAT_FIELD_TYPE);
		case CAT_COMPSTARS:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_NAME) | (1 << CAT_FIELD_TYPE);
		case CAT_LOCAL:
		case CAT_LOCAL_TRIX:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_BMAG);
		case CAT_AN_USER_TEMP:
		case CAT_SHOW:
		case CAT_LOCAL_GAIA_ASTRO:
			return (1 << CAT_FIELD_GAIASOURCEID) | (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC) | (1 << CAT_FIELD_PMRA) | (1 << CAT_FIELD_PMDEC) | (1 << CAT_FIELD_MAG) | (1 << CAT_FIELD_TEFF);
		default:
			return (1 << CAT_FIELD_RA) | (1 << CAT_FIELD_DEC);
	}
}

// This function returns the epoch of the catalog
static double siril_catalog_epoch(siril_cat_index cat) {
	if ((cat == CAT_GAIADR3_DIRECT) || (cat == CAT_LOCAL_GAIA_ASTRO) || (cat == CAT_LOCAL_GAIA_XPSAMP))
		return J2016;
	return J2000;
}

static double siril_catalog_ra_multiplier(siril_cat_index cat) {
	if (cat == CAT_LOCAL_GAIA_ASTRO || cat == CAT_LOCAL_GAIA_XPSAMP)
		return 360.0 / (double) INT32_MAX;
	return 0.000001;
}

static double siril_catalog_dec_multiplier(siril_cat_index cat) {
	if (cat == CAT_LOCAL_GAIA_ASTRO || cat == CAT_LOCAL_GAIA_XPSAMP)
		return 360.0 / (double) INT32_MAX;
	return 0.00001;
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
float siril_catalog_get_default_limit_mag(siril_cat_index cat) {
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
const char *catalog_to_str(siril_cat_index cat) {
	switch (cat) {
		case CAT_TYCHO2:
			return _("Tycho-2");
		case CAT_NOMAD:
			return _("NOMAD");
		case CAT_GAIADR3:
			return _("Gaia DR3 (via Vizier)");
		case CAT_GAIADR3_DIRECT:
			return _("Gaia DR3 (direct)");
		case CAT_PPMXL:
			return _("PPMXL");
		case CAT_BSC:
			return _("bright stars");
		case CAT_APASS:
			return _("APASS");
		case CAT_VSX:
			return _("AAVSO Variable stars");
		case CAT_GCVS:
			return _("GCVS");
		case CAT_SIMBAD:
			return _("SIMBAD");
		case CAT_VARISUM:
			return _("Gaia DR3 Variability");
		case CAT_PGC:
			return _("PGC");
		case CAT_EXOPLANETARCHIVE:
			return _("Exoplanet archive");
		case CAT_IMCCE:
			return _("IMCCE solar system");
		case CAT_AAVSO_CHART:
			return _("AAVSO VSP Chart");
		case CAT_LOCAL:
			return _("Tycho-2+NOMAD");
		case CAT_LOCAL_GAIA_ASTRO:
			return _("Gaia DR3 astrometry");
		case CAT_LOCAL_GAIA_XPSAMP:
			return _("Gaia DR3 xp_sampled");
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
		case CAT_AN_CONST:
			return "IAU constellations";
		case CAT_AN_CONST_NAME:
			return "IAU constellations names";
		case CAT_AN_USER_DSO:
			return "user-DSO";
		case CAT_AN_USER_SSO:
			return "user-SSO";
		case CAT_AN_USER_TEMP:
			return "user-temp";
		case CAT_SHOW:
			return "input";
		case CAT_LOCAL_TRIX:
			return "trixel";
		default:
			return _("unknown");
	}
}

// returns TRUE if the catalog is a star catalog
// Used to set display diameters
gboolean is_star_catalogue(siril_cat_index Catalog) {
	switch (Catalog) {
		case CAT_TYCHO2 ...	CAT_VARISUM:
		case CAT_EXOPLANETARCHIVE:
		case CAT_AAVSO_CHART:
		case CAT_AN_STARS:
		case CAT_LOCAL:
		case CAT_LOCAL_TRIX:
		case CAT_AN_USER_SSO:
		case CAT_LOCAL_GAIA_ASTRO:
		case CAT_LOCAL_GAIA_XPSAMP:
			return TRUE;
	default:
		return FALSE;
	}
}

// returns TRUE if the catalog should display name information
// do not add unless the number of objects returned is limited
// used to set default for conesearch command
gboolean display_names_for_catalogue(siril_cat_index Catalog) {
	switch (Catalog) {
		case CAT_BSC:
		case CAT_GCVS:
		case CAT_PGC:
		case CAT_EXOPLANETARCHIVE:
		case CAT_AAVSO_CHART:
		case CAT_IMCCE:
			return TRUE;
	default:
		return FALSE;
	}
}

static gboolean find_and_check_cat_columns(gchar **fields, int nbcols, siril_cat_index Catalog, int *indexes, uint32_t *collist) {
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

static void fill_cat_item(cat_item *item, const gchar *input, cat_fields index) {
	switch (index) {
		case CAT_FIELD_RA:
			item->ra = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_DEC:
			item->dec = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_RA1:
			item->ra1 = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_DEC1:
			item->dec1 = g_ascii_strtod(input, NULL);
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
		case CAT_FIELD_TEFF:
			item->teff = g_ascii_strtod(input, NULL);
			break;
		case CAT_FIELD_GAIASOURCEID:
			item->gaiasourceid = g_ascii_strtoull(input, NULL, 10);
			break;
		case CAT_FIELD_UNDEF: // columns with unknown headers
		default:
			break;
	}
}

// copies all the data from a catalogue to another
void siril_catalogue_copy(siril_catalogue *from, siril_catalogue *to, gboolean metadata_only) {
	if (!from || !to) {
		siril_debug_print("no catalogue to copy from or to\n");
		return;
	}
	memcpy(to, from, sizeof(siril_catalogue));
	if (from->header)
		to->header = g_strdup(from->header);
	if (from->IAUcode)
		to->IAUcode = g_strdup(from->IAUcode);
	if (from->dateobs)
		to->dateobs = g_date_time_add(from->dateobs, 0); // makes a copy
	if (!metadata_only && from->cat_items) {
		to->cat_items = calloc(to->nbitems, sizeof(cat_item));
		for (int i = 0; i < to->nbitems; i++ )
			siril_catalogue_copy_item(from->cat_items + i, to->cat_items + i);
	} else {
		to->cat_items = NULL;
		to->nbitems = 0;
		to->nbincluded = -1;
		to->projected = CAT_PROJ_NONE;
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

// allocates a new siril_catalogue and initializes it
siril_catalogue *siril_catalog_new(siril_cat_index Catalog) {
	siril_catalogue *siril_cat = calloc(1, sizeof(siril_catalogue));
	siril_cat->cat_index = Catalog;
	siril_cat->columns = siril_catalog_columns(siril_cat->cat_index);
	siril_cat->epoch = siril_catalog_epoch(siril_cat->cat_index);
	siril_cat->ra_multiplier = siril_catalog_ra_multiplier(siril_cat->cat_index);
	siril_cat->dec_multiplier = siril_catalog_dec_multiplier(siril_cat->cat_index);
	return siril_cat;
}

// frees a siril_catalogue
void siril_catalog_free(siril_catalogue *siril_cat) {
	if (!siril_cat)
		return;
	siril_catalog_free_items(siril_cat);
	g_free(siril_cat->IAUcode);
	g_free(siril_cat->header);
	free(siril_cat);
	siril_cat = NULL;
}

// frees the member cat_items of a catalogue
// This is useful to launch the same query again, updating simply the catalog center for instance
void siril_catalog_free_items(siril_catalogue *siril_cat) {
	if (!siril_cat || !siril_cat->cat_items)
		return;
	for (int i = 0; i < siril_cat->nbitems; i++)
		siril_catalog_free_item(&siril_cat->cat_items[i]);
	free(siril_cat->cat_items);
	siril_cat->cat_items = NULL;
	siril_cat->nbitems = 0;
	siril_cat->nbincluded = -1;
	siril_cat->projected = CAT_PROJ_NONE;
}

// frees only one cat_items
void siril_catalog_free_item(cat_item *item) {
	if (!item)
		return;
	g_free(item->name);
	g_free(item->alias);
	g_free(item->type);
	free(item->xp_sampled);
}

void siril_catalog_reset_projection(siril_catalogue *siril_cat) {
	if (!siril_cat || siril_cat->projected == CAT_PROJ_NONE)
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
siril_catalogue *siril_catalog_fill_from_fit(fits *fit, siril_cat_index cat, float limit_mag) {
	if (!fit) {
		return NULL;
	}
	if (!has_wcs(fit)) {
		siril_debug_print("This only works on plate solved images\n");
		return NULL;
	}
	siril_catalogue *siril_cat = siril_catalog_new(cat);
	double ra, dec;
	center2wcs(fit, &ra, &dec);
	double resolution = get_wcs_image_resolution(fit) * 3600.;
	if (limit_mag == -1)
		limit_mag = siril_catalog_get_default_limit_mag(cat);
	// Preparing the catalogue query
	siril_cat->center_ra = ra;
	siril_cat->center_dec = dec;
	siril_cat->radius = get_radius_deg(resolution, fit->rx, fit->ry) * 60.;
	siril_cat->limitmag = limit_mag;
	if (fit->keywords.date_obs) {
		siril_cat->dateobs = fit->keywords.date_obs;
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
 * - The PCC reads them and projects stars on a plate-solved image using WCS.
 *
 * - Comparison star list creation needs equatorial coordinates and B-V
 *   magnitudes, projection is also used but only to check if a star is inside
 *   the image and far from its borders. The required object is psf_star.
 *
 *   It returns 0 on failure
 *   -1 on an empty (but successful) query
 *   nb of stars otherwise
 */

int siril_catalog_conesearch(siril_catalogue *siril_cat) {
	int nbstars = 0;
	if (siril_cat->cat_items && siril_cat->cat_index != CAT_SHOW) {
		siril_debug_print("trying to fetch a catalog while a list already exists, should not happen\n");
		return 0;
	}
	if (siril_cat->cat_index < CAT_AN_MESSIER) {
#ifndef HAVE_LIBCURL
		siril_log_color_message(_("Siril was compiled without networking support, cannot do this operation\n"), "red");
		return 0;
#else
		nbstars = siril_catalog_get_stars_from_online_catalogues(siril_cat);
		return nbstars;
#endif
	} else if (siril_cat->cat_index == CAT_LOCAL || siril_cat->cat_index == CAT_LOCAL_GAIA_ASTRO || siril_cat->cat_index == CAT_LOCAL_GAIA_XPSAMP || siril_cat->cat_index == CAT_LOCAL_TRIX) {
		nbstars = siril_catalog_get_stars_from_local_catalogues(siril_cat);
	} else if (siril_cat->cat_index == CAT_SHOW) { // for the show command
		nbstars = siril_cat->nbitems;
	} else {
		siril_debug_print("trying to conesearch an invalid catalog type\n");
	}

	return nbstars;
}


// This is the generic parser for all csv catalogues used by Siril
// (annotation, downloaded, nina lists...)
// Loads the csv file and fills the cat_items members of the catalogue in entry
int siril_catalog_load_from_file(siril_catalogue *siril_cat, const gchar *filename) {
	GError *error = NULL;
	GFile *catalog_file = g_file_new_for_path(filename);
	GInputStream *input_stream = (GInputStream*) g_file_read(catalog_file, NULL, &error);
	int retval = 1;
	if (!input_stream) {
		if (error != NULL) {
			siril_log_message(_("Could not open the catalog (%s)\n"), error->message);
			g_clear_error(&error);
		} else
			siril_log_message(_("Could not open the catalog (%s)\n"), "generic error");
		return retval;
	}

	int nb_alloc = 1200, nb_items = 0;
	cat_item *cat_items = calloc(nb_alloc, sizeof(cat_item));
	if (!cat_items) {
		PRINT_ALLOC_ERR;
		g_object_unref(input_stream);
		return retval;
	}

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	gchar *line;
	gboolean header_read = FALSE, has_error = FALSE;
	int *indexes = NULL;
	int nbcols = 0;
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
			int offset = 0;
			if (g_str_has_prefix(line, "\uFEFF")) { // see https://en.wikipedia.org/wiki/Byte_order_mark
				printf("BOM detected\n");
				offset = 3;
			}
			gchar **fields = g_strsplit(line + offset, ",", -1);
			nbcols = g_strv_length(fields);
			if (!nbcols) {
				siril_log_color_message(_("No columns found in the catalogue\n"), "red");
				has_error = TRUE;
			} else {
				indexes = malloc(nbcols * sizeof(int));
				for (int i = 0; i < nbcols; i++) {
					indexes[i] = -1;
				}
				has_error = !find_and_check_cat_columns(fields, nbcols, siril_cat->cat_index, indexes, &siril_cat->columns);
			}
			g_strfreev(fields);
			if (has_error) {
				goto siril_catalog_load_from_file_exit_on_error;
			}
			header_read = TRUE;
			g_free(line);
			continue;
		}
		gchar **vals = g_strsplit(line, ",", -1);
		int size = g_strv_length(vals);
		if (size != nbcols) { // checking that current line has a number of columns consistent with the headers
			siril_log_color_message(_("Malformed line found %s\n"), "red", line);
			g_strfreev(vals);
			goto siril_catalog_load_from_file_exit_on_error;
		}
		if (nb_items >= nb_alloc) { // re-allocating if there is more to read
			nb_alloc *= 2;
			cat_item *new_array = realloc(cat_items, nb_alloc * sizeof(cat_item));
			if (!new_array) {
				PRINT_ALLOC_ERR;
				g_strfreev(vals);
				goto siril_catalog_load_from_file_exit_on_error;
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
		line = NULL;
		g_strfreev(vals);
	}
	if (nb_items == 0) {
		siril_log_color_message(_("Catalog %s was read but no items were found in the view cone, nothing to show\n"), "salmon", filename);
		retval = -1;
		goto siril_catalog_load_from_file_exit_on_error;
	}
	cat_item *final_array = realloc(cat_items, nb_items * sizeof(cat_item));
	siril_cat->cat_items = final_array;
	siril_cat->nbitems = nb_items;
	siril_cat->nbincluded = nb_items;
	if (header)
		siril_cat->header = g_string_free(header, FALSE);
	siril_debug_print("read %d%s items from catalogue\n", nb_items, siril_cat->phot ? " photometric" : "");
	g_object_unref(data_input);
	g_object_unref(input_stream);

	return 0;
siril_catalog_load_from_file_exit_on_error:
	free(cat_items);
	siril_cat->cat_items = NULL;
	siril_cat->nbitems = 0;
	free(indexes);
	g_object_unref(data_input);
	g_object_unref(input_stream);
	g_free(line);
	if (header)
		g_string_free(header, TRUE);
	return retval;

}

// Writes the catalogue to the given output_stream
gboolean siril_catalog_write_to_output_stream(siril_catalogue *siril_cat, GOutputStream *output_stream) {
	gsize n;
	GError *error = NULL;
	if (siril_cat->header && !g_output_stream_printf(output_stream, &n, NULL, NULL, "%s\n", siril_cat->header)) {
		g_error_free(error);
		return FALSE;
	}

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
		if (!g_output_stream_printf(output_stream, &n, NULL, &error, "%s", columns_str)) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			g_free(columns_str);
			return FALSE;
		}
		g_free(columns_str);
	} else {
		siril_debug_print("no columns to write");
		g_error_free(error);
		return FALSE;
	}
	for (int j = 0; j < siril_cat->nbitems; j++) {
		gchar **tokens = calloc(nbcols + 1, sizeof(gchar *));
		for (int i = 0; i < nbcols; i++) {
			tokens[i] = get_field_to_str(&siril_cat->cat_items[j], index[i]);
		}
		gchar *newline = g_strjoinv(",", tokens);
		if (!g_output_stream_printf(output_stream, &n, NULL, &error, "\n%s", newline)) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			g_free(newline);
			return FALSE;
		}
		g_free(newline);
	}
	g_clear_error(&error);
	return TRUE;
}

// Writes the catalogue to the given filepath
gboolean siril_catalog_write_to_file(siril_catalogue *siril_cat, const gchar *filename) {
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
	g_free(root);

	// Then we delete it if it exists
	GError *error = NULL;
	GFile *file = g_file_new_for_path(filename);
	if (g_file_test(filename, G_FILE_TEST_EXISTS)) { // file already exists, we removed it
		if (!g_file_delete(file, NULL, &error)) {
			siril_log_color_message(_("File %s cannot be deleted (%s), aborting\n"), "red", filename, (error) ? error->message : "unknown error");
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
	gboolean success = siril_catalog_write_to_output_stream(siril_cat, output_stream);
	g_object_unref(file);
	g_object_unref(output_stream);
	return success;
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

static gboolean can_use_proper_motion(fits *fit, siril_catalogue *siril_cat) {
	if (!fit || !siril_cat)
		return FALSE;
	if (!fit->keywords.date_obs) {
		// siril_debug_print("This image does not have any DATE-OBS information, cannot account for stars proper motions\n");
		return FALSE;
	}
	if (!has_field(siril_cat, PMRA) || !has_field(siril_cat, PMDEC)) {
		// siril_debug_print("This catalog does not have proper motion info, will not be computed\n");
		return FALSE;
	}
	siril_debug_print("Catalogue %s will account for proper motions\n", catalog_to_str(siril_cat->cat_index));
	return TRUE;
}

static gboolean can_use_velocity(fits *fit, siril_catalogue *siril_cat) {
	if (!fit || !siril_cat)
		return FALSE;
	if (!fit->keywords.date_obs) {
		// siril_debug_print("This image does not have any DATE-OBS information, cannot account for velocity\n");
		return FALSE;
	}
	if (!has_field(siril_cat, VRA) || !has_field(siril_cat, VDEC)) {
		// siril_debug_print("This catalog does not have velocity info, will not be computed\n");
		return FALSE;
	}
	siril_debug_print("Catalogue %s will account for velocities\n", catalog_to_str(siril_cat->cat_index));
	return TRUE;
}

// projects passed catalogue using the wcs data contained in the fit
// corrects for proper motions if the flag is TRUE and the necessary data is included
// in the fit (dateobs) and in the catalogue (pmra and pmdec fields)
// corrects for object velocity if the flag is true and if necessary data is included
// in the catalogue (vra and vdec fields)
int siril_catalog_project_with_WCS(siril_catalogue *siril_cat, fits *fit, gboolean use_proper_motion, gboolean use_velocity) {
	if (!has_field(siril_cat, RA) || !has_field(siril_cat, DEC)) {
		siril_debug_print("catalogue %s does not have the necessary columns\n", catalog_to_str(siril_cat->cat_index));
		return 1;
	}
	int nbincluded = 0;
	double jyears = 0.;
	double *world = NULL, *x = NULL, *y = NULL;
	int *status = NULL;
	double tobs = 0., deltahours = 0.;
	gboolean use_tcat = FALSE;
	use_proper_motion = use_proper_motion && can_use_proper_motion(fit, siril_cat);
	use_velocity = use_velocity && can_use_velocity(fit, siril_cat);
	gboolean has_second_star = siril_cat->cat_index == CAT_AN_CONST;
	if (use_proper_motion) {
		GDateTime *dt = g_date_time_ref(fit->keywords.date_obs);
		gdouble jd = date_time_to_Julian(dt);
		g_date_time_unref(dt);
		jyears = (jd - siril_cat->epoch) / 365.25;
	}
	if (use_velocity) {
		GDateTime *dt = g_date_time_ref(fit->keywords.date_obs);
		tobs = date_time_to_Julian(dt);
		g_date_time_unref(dt);
		// for IMCCE conesearch, the dateobs is common to the whole catalogue
		// the time of the catalogue will be used instead of individual records
		if (!has_field(siril_cat, DATEOBS)) {
			if (siril_cat->dateobs) {
				deltahours = (tobs - date_time_to_Julian(siril_cat->dateobs)) * 24.;
				use_tcat = TRUE;
			} else {
				siril_log_color_message(_("Cannot correct for velocities, will use ra/dec as is\n"), "salmon");
				siril_log_color_message(_("To tag solar system objects, please refer to conesearch command\n"), "salmon");
				use_velocity = FALSE;
			}
		}
	}
	int nbinit = siril_cat->nbitems;
	if (has_second_star)
		nbinit *= 2;
	world = malloc(2 * nbinit * sizeof(double));
	x = malloc(nbinit * sizeof(double));
	y = malloc(nbinit * sizeof(double));
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
			if (!use_tcat)
				deltahours = (tobs - siril_cat->cat_items[i].dateobs) * 24.;
			ra += siril_cat->cat_items[i].vra / cos(decrad) * deltahours * 2.77777778e-4;
			dec += siril_cat->cat_items[i].vdec * deltahours * 2.77777778e-4;
		}
		world[ind++] = ra;
		world[ind++] = dec;
		if (has_second_star) {
			world[ind++] = siril_cat->cat_items[i].ra1;
			world[ind++] = siril_cat->cat_items[i].dec1;
		}
	}
	status = wcs2pix_array(fit, nbinit, world, x, y);
	if (!status)
		goto clean_and_exit;
	if (!has_second_star) {
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
	} else {
		int j = 0;
		for (int i = 0; i < siril_cat->nbitems; i++) {
			if (!status[j] || !status[j + 1]) {
				siril_cat->cat_items[i].x = x[j];
				siril_cat->cat_items[i].y = y[j];
				siril_cat->cat_items[i].x1 = x[j + 1];
				siril_cat->cat_items[i].y1 = y[j + 1];
				siril_cat->cat_items[i].included = TRUE;
				nbincluded++;
			} else {
				siril_cat->cat_items[i].included = FALSE;
			}
			j += 2;
		}
	}
clean_and_exit:
	siril_cat->nbincluded = nbincluded;
	free(status);
	free(world);
	free(x);
	free(y);
	siril_cat->projected = CAT_PROJ_WCS;
	return !(nbincluded > 0);
}

// projects passed catalogue wrt to the point ra0 and dec0 coordinates
// according to gnomonic (a.k.a TAN) projection
// corrects for proper motions if the flag is TRUE and the necessary data is passed
// (dateobs) and found in the catalogue (pmra and pmdec fields)
int siril_catalog_project_gnomonic(siril_catalogue *siril_cat, double ra0, double dec0, gboolean use_proper_motion, GDateTime *date_obs) {
	if (!has_field(siril_cat, RA) || !has_field(siril_cat, DEC))
		return 1;
	double jyears = 0.;
	if (use_proper_motion) {
		if (!date_obs) {
			siril_log_color_message(_("no DATE-OBS information, cannot account for stars proper motions\n"), "salmon");
			use_proper_motion = FALSE;
		} else if (!has_field(siril_cat, PMRA) || !has_field(siril_cat, PMDEC)) {
			siril_log_color_message(_("This catalog does not have proper motion info, will not be computed\n"), "salmon");
			use_proper_motion = FALSE;
		} else {
			GDateTime *dt = g_date_time_ref(date_obs);
			gdouble jd = date_time_to_Julian(dt);
			g_date_time_unref(dt);
			jyears = (jd - siril_cat->epoch) / 365.25;
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
	siril_cat->projected = CAT_PROJ_TAN;
	return 0;
}

// TODO: move to a file for callbacks and remove gtk include
static gboolean end_conesearch(gpointer p) {
	siril_catalogue *temp_cat = (siril_catalogue *) p;
	if (temp_cat) {
		// purge_user_catalogue(CAT_AN_USER_TEMP); // we don't clear so as to accumulate displays of various conesearches/show commands
		if (!load_siril_cat_to_temp(temp_cat)) {
			GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button"));
			// refresh_annotation_to_temp(); // may need to make that an option later on
			refresh_annotation_visibility();
			if (!gtk_toggle_tool_button_get_active(button)) {
				gtk_toggle_tool_button_set_active(button, TRUE);
			} else {
				refresh_found_objects();
				redraw(REDRAW_OVERLAY);
			}
		}
	}
	return end_generic(NULL);
	// we don't free temp_cat as it is passed as the new CAT_AN_USER_TEMP
}

// Conesearch command related functions

// Sanity checks for conesearch args
int check_conesearch_args(conesearch_args *args) {
	if (args->display_log && !has_field(args->siril_cat, NAME)) {
		siril_log_color_message(_("Won't log list of objects for catalog %s which does not have name information\n"), "salmon", catalog_to_str(args->siril_cat->cat_index));
		args->display_log = FALSE;
	}
	if (args->display_tag && !args->has_GUI) {
		siril_log_color_message(_("Won't show objects tags as there is no display, ignoring\n"), "salmon", catalog_to_str(args->siril_cat->cat_index));
		args->display_tag = FALSE;
	}
	if (args->display_tag && !has_field(args->siril_cat, NAME)) {
		siril_log_color_message(_("Won't display objects tags for catalog %s which does not have name information\n"), "salmon", catalog_to_str(args->siril_cat->cat_index));
		args->display_tag = FALSE;
	}
	return 0;
}

int execute_conesearch(conesearch_params *params) {
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		free_conesearch_params(params);
		return CMD_FOR_PLATE_SOLVED;
	}

	// Preparing the catalogue query
	siril_catalogue *siril_cat = siril_catalog_fill_from_fit(&gfit, params->cat, params->limit_mag);
	siril_cat->phot = params->photometric;
	if (params->cat == CAT_IMCCE) {
		if (params->obscode) {
			siril_cat->IAUcode = g_strdup(params->obscode);
			if (params->default_obscode_used) {
				siril_log_message(_("Using default observatory code %s\n"), params->obscode);
			}
		} else {
			siril_cat->IAUcode = g_strdup("500");
			siril_log_color_message(_("Did not specify an observatory code, using geocentric by default, positions may not be accurate\n"), "salmon");
		}
	} else if (params->obscode) {
		g_free(params->obscode);
		params->obscode = NULL;
	}
	if (params->cat == CAT_LOCAL_TRIX)
		siril_cat->trixel = params->trixel;

	siril_debug_print("centre coords: %f, %f, radius: %f arcmin\n", siril_cat->center_ra, siril_cat->center_dec, siril_cat->radius);
	conesearch_args *args = init_conesearch_args();
	args->fit = &gfit;
	args->siril_cat = siril_cat;
	args->has_GUI = !com.script;
	args->display_log = (params->display_log == BOOL_NOT_SET) ? display_names_for_catalogue(params->cat) : (gboolean) params->display_log;
	args->display_tag = (params->display_tag == BOOL_NOT_SET) ? display_names_for_catalogue(params->cat) : (gboolean) params->display_tag;
	args->outfilename = g_strdup(params->outfilename);
	args->compare = params->compare;
	free_conesearch_params(params);
	if (check_conesearch_args(args)) { // can't fail for now
		free_conesearch_args(args);
		return CMD_GENERIC_ERROR;
	}

	if (!start_in_new_thread(conesearch_worker, args)) {
		free_conesearch_args(args);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}

int execute_show_command(show_params *params) {
	if (!has_wcs(&gfit)) {
		siril_log_color_message(_("This command only works on plate solved images\n"), "red");
		return CMD_FOR_PLATE_SOLVED;
	}

	if (params->clear) { // this is safe even in script/headless,used by the command, not the UI
		purge_user_catalogue(CAT_AN_USER_TEMP);
		redraw(REDRAW_OVERLAY);
	}

	siril_catalogue *siril_cat = siril_catalog_new(CAT_SHOW);	
	conesearch_args *args = init_conesearch_args();
	args->siril_cat = siril_cat;
	args->has_GUI = TRUE;
	args->fit = &gfit;

	if (params->file) {
		int check = siril_catalog_load_from_file(siril_cat, params->file);
		if (check > 0) {
			free_conesearch_args(args);
			return CMD_ARG_ERROR;
		}
		if (check == -1) {
			free_conesearch_args(args);
			return CMD_OK;
		}
		args->display_log =
				(params->display_log == BOOL_NOT_SET) ?
						(gboolean) has_field(siril_cat, NAME) :
						(gboolean) params->display_log;
		args->display_tag =
				(params->display_tag == BOOL_NOT_SET) ?
						(gboolean) has_field(siril_cat, NAME) :
						(gboolean) params->display_tag;
		if (!start_in_new_thread(conesearch_worker, args)) {
			free_conesearch_args(args);
			return CMD_GENERIC_ERROR;
		}
	return CMD_OK;
	}

	if (params->coords) {
		cat_item *item = calloc(1, sizeof(cat_item));
		item->ra = siril_world_cs_get_alpha(params->coords);
		item->dec = siril_world_cs_get_delta(params->coords);
		siril_world_cs_unref(params->coords);
		item->name = g_strdup(params->name);
		siril_cat->columns |= (1 << CAT_FIELD_NAME);
		args->display_log = TRUE;
		args->display_tag = TRUE;

		siril_catalog_append_item(siril_cat, item);
		siril_catalog_free_item(item);
		free(item);
		if (!start_in_new_thread(conesearch_worker, args)) {
			free_conesearch_args(args);
			return CMD_GENERIC_ERROR;
		}
		return CMD_OK;
	}

	free_conesearch_args(args);
	return CMD_ARG_ERROR;
}

// worker for the command conesearch
gpointer conesearch_worker(gpointer p) {
	conesearch_args *args = (conesearch_args *)p;
	siril_catalogue *siril_cat = args->siril_cat;
	siril_catalogue *temp_cat = NULL;
	double *dx = NULL, *dy = NULL, *dxf = NULL, *dyf = NULL;
	siril_plot_data *spl_data = NULL;
	int retval = 1;
	double stardiam = 0.0;
	gboolean hide_display_tag = FALSE;

	// Check initial args
	if (!siril_cat) {
		siril_debug_print("no query passed");
		goto exit_conesearch;
	}
	if (args->compare && !args->has_GUI) {
		args->compare = FALSE;
		siril_debug_print("Compare available only with GUI\n");
	}

	// Launch the query
	int check = siril_catalog_conesearch(siril_cat);
	if (!check) {  // conesearch failed
		goto exit_conesearch;
	}
	if (siril_cat->cat_index != CAT_LOCAL &&
		siril_cat->cat_index != CAT_LOCAL_GAIA_ASTRO &&
		siril_cat->cat_index != CAT_LOCAL_TRIX)
	{
		siril_log_message(_("The %s catalog has been successfully downloaded\n"),
						  catalog_to_str(siril_cat->cat_index));
	}
	if (check == -1) {  // conesearch succeeded but the field was empty
		retval = -1;
		goto exit_conesearch;
	}

	// Project using WCS
	if (siril_catalog_project_with_WCS(siril_cat, args->fit, TRUE, TRUE)) {
		if (siril_cat->projected == CAT_PROJ_WCS)
			siril_log_color_message(_("No item found in the image\n"), "salmon");
		goto exit_conesearch;
	}
	int nb_stars = siril_cat->nbitems;
	if (siril_cat->cat_index != CAT_SHOW)
		sort_cat_items_by_mag(siril_cat);
	int j = 0, k = 0;

	// Prepare the temporary annotation catalogue if GUI or output file is requested
	if (args->has_GUI || args->outfilename) {
		temp_cat = siril_catalog_new(CAT_AN_USER_TEMP);
		if (!temp_cat) {
			PRINT_ALLOC_ERR;
			retval = 1;
			goto exit_conesearch;
		}
		temp_cat->projected = CAT_PROJ_WCS;
		if (is_star_catalogue(siril_cat->cat_index))
			stardiam = 0.2;  // in arcmin => 12"
		if (!args->display_tag && has_field(siril_cat, NAME))
			hide_display_tag = TRUE;
		temp_cat->cat_items = calloc(siril_cat->nbincluded, sizeof(cat_item));
		if (!temp_cat->cat_items) {
			PRINT_ALLOC_ERR;
			retval = 1;
			goto exit_conesearch;
		}
	}

	// Allocate memory for compare mode if needed
	if (args->compare) {
		dx = calloc(nb_stars, sizeof(double));
		if (!dx) {
			PRINT_ALLOC_ERR;
			retval = 1;
			goto exit_conesearch;
		}
		dy = calloc(nb_stars, sizeof(double));
		if (!dy) {
			PRINT_ALLOC_ERR;
			free(dx);
			dx = NULL;
			retval = 1;
			goto exit_conesearch;
		}
	}

	// Iterate over cat_items
	for (int i = 0; i < nb_stars; i++) {
		if (!siril_cat->cat_items[i].included ||
			(siril_cat->limitmag && siril_cat->cat_items[i].mag > siril_cat->limitmag))
		{
			continue;
		}
		// If GUI is active, copy the item into the temporary catalogue
		if (args->has_GUI) {
			siril_catalogue_copy_item(&siril_cat->cat_items[i],
									  &temp_cat->cat_items[j]);
			if (stardiam)
				temp_cat->cat_items[j].diameter = stardiam;
			if (hide_display_tag) {
				g_free(temp_cat->cat_items[j].name);
				temp_cat->cat_items[j].name = NULL;
			} else {
				if (siril_cat->cat_index == CAT_PGC) {
					g_free(temp_cat->cat_items[j].name);
					temp_cat->cat_items[j].name = g_strdup_printf("PGC %s",
																  siril_cat->cat_items[i].name);
				}
			}
		}
		if (args->display_log) {
			gchar *ra = siril_world_cs_alpha_format_from_double(siril_cat->cat_items[i].ra,
																"%02d %02d %04.1lf");
			gchar *dec = siril_world_cs_delta_format_from_double(siril_cat->cat_items[i].dec,
																 "%c%02d %02d %04.1lf");
			if (siril_cat->cat_index == CAT_AAVSO_CHART) {
				siril_log_message("AUID:%s - V:%3.1f [%5.3f] - B:%3.1f [%5.3f] - RA: %s, DEC: %s\n",
								  siril_cat->cat_items[i].name,
					  siril_cat->cat_items[i].mag,
					  siril_cat->cat_items[i].e_mag,
					  siril_cat->cat_items[i].bmag,
					  siril_cat->cat_items[i].e_bmag,
					  ra,
					  dec);
			} else {
				GString *msg = g_string_new("");
				g_string_append_printf(msg, "%s%s",
									   (siril_cat->cat_index == CAT_PGC) ? "PGC " : "",
									   siril_cat->cat_items[i].name);
				if (has_field(siril_cat, TYPE))
					g_string_append_printf(msg, " (%s)",
										   siril_cat->cat_items[i].type);
				g_string_append_printf(msg, ", ");
				g_string_append_printf(msg, "RA: %s, DEC: %s", ra, dec);
				if (has_field(siril_cat, MAG))
					g_string_append_printf(msg, " , mag:%3.1f",
										   siril_cat->cat_items[i].mag);
				g_string_append_printf(msg, "\n");
				gchar *printout = g_string_free(msg, FALSE);
				siril_log_message(printout);
				g_free(printout);
			}
			g_free(ra);
			ra = NULL;
			g_free(dec);
			dec = NULL;
		}
		if (args->compare) {
			double scale = 1800.0 * (fabs(args->fit->keywords.wcslib->cdelt[0]) +
			fabs(args->fit->keywords.wcslib->cdelt[1]));
			double x = siril_cat->cat_items[i].x;
			double y = siril_cat->cat_items[i].y;
			rectangle area = { 0 };
			cat_item tmp = { .x = x, .y = y };
			if (make_selection_around_a_star(&tmp, &area, args->fit)) {
				siril_debug_print("star %d is outside image or too close to border\n", i);
				continue;
			}
			psf_error error = PSF_NO_ERR;
			int layer = (args->fit->naxes[2] == 3) ? GLAYER : RLAYER;
			psf_star *star = psf_get_minimisation(args->fit, layer, &area, FALSE, FALSE, NULL,
												  FALSE, com.pref.starfinder_conf.profile, &error);
			if (star && !error) {
				dx[k] = area.x + star->x0;
				dy[k] = area.y + area.h - star->y0;
				display_to_siril(dx[k], dy[k], &dx[k], &dy[k], args->fit->ry);
				dx[k] = scale * (dx[k] - x);
				dy[k] = scale * (dy[k] - y);
				free_psf(star);
				star = NULL;
				k++;
			}
		}
		j++;
	}

	// Summary log message
	siril_log_message(_("%d objects found%s in the image (mag limit %.2f) using %s catalogue\n"),
					  j,
				   siril_cat->phot ? " with valid photometry data" : "",
				   siril_cat->limitmag,
				   catalog_to_str(siril_cat->cat_index));
	if (!j) {
		retval = -1;
		goto exit_conesearch;
	}

	// Resize temp catalogue to the actual number of items
	if (args->has_GUI || args->outfilename) {
		cat_item *final_items = realloc(temp_cat->cat_items, j * sizeof(cat_item));
		if (!final_items) {
			PRINT_ALLOC_ERR;
			retval = 1;
			goto exit_conesearch;
		}
		temp_cat->cat_items = final_items;
		temp_cat->nbitems = j;
	}

	// Memory reallocation to fit actual number of objects
	if (args->compare && k > 0) {
		double *tmp_dxf = realloc(dx, k * sizeof(double));
		if (!tmp_dxf) {
			PRINT_ALLOC_ERR;
			free(dx);
			dx = NULL;
			free(dy);
			dy = NULL;
			retval = 1;
			goto exit_conesearch;
		}
		dxf = tmp_dxf;
		double *tmp_dyf = realloc(dy, k * sizeof(double));
		if (!tmp_dyf) {
			PRINT_ALLOC_ERR;
			free(dxf);
			dxf = NULL;
			free(dy);
			dy = NULL;
			retval = 1;
			goto exit_conesearch;
		}
		dyf = tmp_dyf;
	} else {
		args->compare = FALSE;
		free(dx);
		dx = NULL;
		free(dy);
		dy = NULL;
	}

	// Write catalogue if required
	if (args->outfilename) {
		if (siril_catalog_write_to_file(temp_cat, args->outfilename)) {
			siril_log_message(_("List saved to %s\n"), args->outfilename);
		} else {
			siril_log_message(_("Failed to save list to %s\n"), args->outfilename);
		}
	}

	if (args->has_GUI && args->compare) {
		spl_data = init_siril_plot_data();
		if (!spl_data) {
			retval = 1;
		} else {
			siril_plot_set_title(spl_data, "Detected vs astrometric position");
			siril_plot_set_xlabel(spl_data, "dx [\"]");
			siril_plot_set_ylabel(spl_data, "dy [\"]");
			siril_plot_add_xydata(spl_data, NULL, k, dxf, dyf, NULL, NULL);
			siril_plot_set_savename(spl_data, "diffpos");
			siril_plot_set_nth_plot_type(spl_data, 1, KPLOT_POINTS);
		}
		free(dxf);
		dxf = NULL;
		free(dyf);
		dyf = NULL;
	}

	retval = 0;

	exit_conesearch:
	{
		gboolean go_idle = args->has_GUI;
		if ((retval || !args->has_GUI) && temp_cat) {
			siril_catalog_free(temp_cat);
			temp_cat = NULL;
		}
		free_conesearch_args(args);
		args = NULL;

		if (go_idle) {
			if (spl_data)
				siril_add_pythonsafe_idle(create_new_siril_plot_window, spl_data);
			siril_add_idle(end_conesearch, temp_cat);
		} else {
			end_generic(NULL);
		}
		if (retval == -1)  // success but empty field
			retval = 0;
		free(dx); // may still be NULL if the if (!j) conditional bails out
		free(dy); // may still be NULL if the if (!j) conditional bails out
	}

	return GINT_TO_POINTER(retval);
}

// Haversine formula on unit sphere
// https://en.wikipedia.org/wiki/Haversine_formula
// dec is phi, ra is lambda
// in degrees

double compute_coords_distance_h(double ra1, double dec1, double ra2, double dec2) {
	double dec1_r = dec1 * DEGTORAD, dec2_r = dec2 * DEGTORAD;
	double dra_2 = 0.5 * (ra2 - ra1) * DEGTORAD;
	double ddec_2 = 0.5 * (dec2_r - dec1_r);
	double sin_ddec = sin(ddec_2), sin_dra = sin(dra_2);
	double h = sin_ddec * sin_ddec + cos(dec1_r) * cos(dec2_r) * sin_dra * sin_dra;
	if (h > 1.)
		return 1.;   // h = 1, asin(1) is pi/2
	return h;
}

double compute_coords_distance(double ra1, double dec1, double ra2, double dec2) {
	double h = compute_coords_distance_h(ra1, dec1, ra2, dec2);
	if (h > 1.)
		return 180.0;   // h = 1, asin(1) is pi/2
	return 2.0 * asin(sqrt(h)) * RADTODEG;
}

int siril_catalog_inner_conesearch(siril_catalogue *siril_cat_in, siril_catalogue *siril_cat_out) {
	if (!siril_cat_in)
		return 0;
	int nb_alloc = 1200, nb_items = 0;
	cat_item *cat_items = calloc(nb_alloc, sizeof(cat_item));
	if (!cat_items) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	double ra = siril_cat_out->center_ra;
	double dec = siril_cat_out->center_dec;

	double radius_h = pow(sin(0.5 * siril_cat_out->radius / 60. * DEGTORAD), 2);
	for (int i = 0; i < siril_cat_in->nbitems; i++) {
		if (nb_items >= nb_alloc) { // re-allocating if there is more to read
			nb_alloc *= 2;
			cat_item *new_array = realloc(cat_items, nb_alloc * sizeof(cat_item));
			if (!new_array) {
				PRINT_ALLOC_ERR;
				return -1;
			}
			cat_items = new_array;
		}
		double dist_h = compute_coords_distance_h(ra, dec, siril_cat_in->cat_items[i].ra, siril_cat_in->cat_items[i].dec);
		if (dist_h <= radius_h) {
			siril_catalogue_copy_item(siril_cat_in->cat_items + i, cat_items + nb_items);
			nb_items++;
		}
	}
	cat_item *final_array = realloc(cat_items, nb_items * sizeof(cat_item));
	siril_cat_out->cat_items = final_array;
	siril_cat_out->nbitems = nb_items;
	siril_cat_out->nbincluded = nb_items;
	return nb_items;
}

// TODO: using this for the moment to avoid chaging too many files
// This copies the info contained in the catalogue to a psf_star** list
// only the included items are copied over
psf_star **convert_siril_cat_to_psf_stars(siril_catalogue *siril_cat) {
	if (!siril_cat)
		return NULL;
	if (siril_cat->projected == CAT_PROJ_NONE) {
		siril_debug_print("Catalog has not been projected\n");
	}
	if (!has_field(siril_cat, RA) || !has_field(siril_cat, DEC) || !has_field(siril_cat, MAG))
		return NULL;
	psf_star **results = new_fitted_stars(siril_cat->nbincluded + 1);

	int n = 0;
	for (int i = 0; i < siril_cat->nbitems; i++) {
		if (siril_cat->cat_items[i].included) {
			if (n >= siril_cat->nbincluded) {
				siril_debug_print("problem when converting siril_cat to psf_stars, more than allocated");
				break;
			}
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

conesearch_args *init_conesearch_args() {
	conesearch_args *args = calloc(1, sizeof(conesearch_args));
	if (!args) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	return args;
}

static void free_conesearch_args(conesearch_args *args) {
	if (!args)
		return;
	siril_catalog_free(args->siril_cat);
	g_free(args->outfilename);
	free(args);
}

conesearch_params *init_conesearch_params() {
	conesearch_params *params = g_new0(conesearch_params, 1);
	params->limit_mag = -1.0f;
	params->photometric = FALSE;
	params->display_tag = BOOL_NOT_SET;
	params->display_log = BOOL_NOT_SET;
	params->cat = CAT_AUTO;
	params->obscode = NULL;
	params->default_obscode_used = FALSE;
	params->trixel = -1;
	params->outfilename = NULL;
	if (com.pref.astrometry.default_obscode != NULL) {
		params->obscode = g_strdup(com.pref.astrometry.default_obscode);
		params->default_obscode_used = TRUE;
	}
	return params;
}

static void free_conesearch_params(conesearch_params *params) {
	if (!params)
		return;
	g_free(params->obscode);
	g_free(params->outfilename);
	g_free(params); // was alloced with gnew0
}

