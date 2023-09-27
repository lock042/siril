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
#ifndef _SIRIL_CATALOGUES_H
#define _SIRIL_CATALOGUES_H

#include <glib.h>
#include <gio/gio.h>
#include "core/siril_world_cs.h"

// number of columns that can be defined in a catalogue
#define MAX_CAT_COLUMNS 18

// all catalogues that can be used
// < 60: online
// < 40: TAP
// if value < 20: online stars catalogues (through TAP Vizier)
// if 20 <= value < 30: online other catalogues (through TAP Vizier) - galaxies etc
// if 30 <= value < 40: online exoplanets (only through TAP for now - no parser needed)
// if 40 <= value < 50: AAVSO charts (parser to csv needed)
// if 50 <= value < 60: online solar system (parser to csv needed)
// if 60 <= value < 70: local annotation catalogues
// if 70 <= value < 90: local annotation catalogues (user-defined)
// if 90 <= value <= 100: special use cases
typedef enum {
// TAP Queries from vizier
	CAT_TYCHO2,
	CAT_NOMAD,
	CAT_GAIADR3,
	CAT_PPMXL,
	CAT_BSC,
	CAT_APASS,
	CAT_GCVS,
	CAT_VSX,
	CAT_SIMBAD,
	CAT_PGC = 20,
// Other TAP Queries
	CAT_EXOPLANETARCHIVE = 30,
// Non TAP Queries (stars)
	CAT_AAVSO_CHART = 40,
// Non TAP Queries (others)
	CAT_IMCCE = 50,
// Local annotations - same as annotations_cat + 60
	CAT_AN_MESSIER = 60,
	CAT_AN_NGC = 61,
	CAT_AN_IC = 62,
	CAT_AN_LDN = 63,
	CAT_AN_SH2 = 64,
	CAT_AN_STARS = 65,
	CAT_AN_USER_DSO = 66,
	CAT_AN_USER_SSO = 67,
	CAT_AN_USER_TEMP = 68,
// Special
	CAT_COMPSTARS = 97,
	CAT_AUTO = 98,
	CAT_LOCAL = 99,		// siril local (KStars Tycho-2 and NOMAD)
	CAT_ASNET = 100,	// solve-field local (astrometry.net)
} object_catalog;

typedef enum {
	CAT_FIELD_UNDEF = -1,
	CAT_FIELD_RA,
	CAT_FIELD_DEC,
	CAT_FIELD_PMRA,
	CAT_FIELD_PMDEC,
	CAT_FIELD_MAG,
	CAT_FIELD_BMAG,
	CAT_FIELD_E_MAG,
	CAT_FIELD_E_BMAG,
	CAT_FIELD_NAME,
	CAT_FIELD_DIAMETER,
	CAT_FIELD_ALIAS,
	CAT_FIELD_DATEOBS,
	CAT_FIELD_SITELAT,
	CAT_FIELD_SITELON,
	CAT_FIELD_SITEELEV,
	CAT_FIELD_VRA,
	CAT_FIELD_VDEC,
	CAT_FIELD_TYPE
} cat_fields;

typedef struct {
	// filled from catalogue
	double ra, dec;	// celestial coordinates
	double pmra, pmdec;	// proper motions
	float mag;	// visible magnitude (V filter), for sorting and debug
	float bmag;	// B magnitude
	float e_mag, e_bmag;	// mag errors for comp_stars
	float diameter; // to be used with annotations catalogues
	double dateobs; // for solar system (in JD)
	double sitelat, sitelon, siteelev; // obs site for solar system
	double vra, vdec; // speed vector in ra/dec for solar system
	gchar *name;  // name of the object
	gchar *alias; // aliases given in annotation catalogues, tab-separated
	gchar *type; // type of the object, for solsys and compstars

	// computed
	float x, y;	// image coordinates
	gboolean included; // flag to remove items from the list without deleting them (to be used by platesolve/pcc)
} cat_item;

typedef struct {
	object_catalog cattype;
	double catalog_center_ra;
	double catalog_center_dec;
	double radius; // fov radius (unit TBD)
	double limitmag; // limiting magnitude
	double dateobs; // date-obs in JD
	double sitelat, sitelon, siteelev; // obs site
	gchar *IAUcode; // observatory code
	gboolean phot; // TRUE if can be used for photometry
	cat_item *cat_items;
	int nbitems; // the number of items stored
	int nbincluded; // the number of items included after projection
	uint32_t columns; // the list of columns which where parsed when read
} siril_catalogue;


uint32_t siril_catalog_colums(object_catalog cat);
void sort_cat_items_by_mag(siril_catalogue *siril_cat);
const char *catalog_to_str(object_catalog cat);
const gchar **get_cat_colums_names();

siril_catalogue *siril_catalog_load_from_file(const gchar *filename, gboolean phot);
int siril_catalog_project_with_WCS(siril_catalogue *siril_cat, fits *fit, gboolean use_proper_motion);
int siril_catalog_project_at_center(siril_catalogue *siril_cat, double ra0, double dec0, gboolean use_proper_motion, GDateTime *date_obs);

psf_star **convert_siril_cat_to_psf_stars(siril_catalogue *siril_cat, int *nbstars);

#endif
