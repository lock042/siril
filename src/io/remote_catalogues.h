#ifndef _REMOTE_CATALOGUES_H
#define _REMOTE_CATALOGUES_H

#include <glib.h>
#include <gio/gio.h>
#include "core/siril_world_cs.h"
#include "algos/PSF.h"


#define VIZIER_QUERY "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source="

// new queries
#define VIZIER_TAP_QUERY "http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=SELECT+"
#define EXOPLANET_TAP_QUERY "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=select+"
#define IMCCE_QUERY "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?&-mime=text&-output=basic&-filter=0&-objFilter=111&-refsys=EQJ2000&-from=Siril"

// #define AAVSOCHART_QUERY "https://app.aavso.org/vsp/api/chart/?"

// only the first 9 columns are valid TAP queries
// fields after this are used in other catalogues
#define MAX_TAP_QUERY_COLUMNS 9
#define MAX_CAT_COLUMNS 18
#define MAX_CATCODE_LEN 30


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

// The structure used to declare the columns to be queried from the tables
// for TAP queries only
typedef struct {
	gchar catcode[MAX_CATCODE_LEN];
	gchar *tap_columns[MAX_TAP_QUERY_COLUMNS];
} cat_tap_query_fields;

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
	float BV; // B magnitude - V magnitude, -99.9 if not available - could be read directly from APASS
	gboolean included; // flag to remove items from the list without deleting them (to be used by platesolve/pcc)
} cat_item;

typedef struct {
	object_catalog cattype;
	double catalog_center_ra;
	double catalog_center_dec;
	double radius; // fov radius (unit TBD)
	double limitmag; // limiting magnitude
	double dateobs; // date-obs in JD
	double sitelat, sitelon, siteelev;
	gchar *IAUcode; // observatory code
	gboolean phot; // TRUE if can be used for photometry
	cat_item *cat_items;
	int nbitems; // the number of items stored
	uint32_t columns; // the list of columns which where parsed when read
} siril_catalog;


const char *catalog_to_str(object_catalog cat);

GFile *download_catalog(object_catalog onlineCatalog, SirilWorldCS *catalog_center, double radius_arcmin, double mag, gchar *obscode, GDateTime *date_obs);
gchar *get_catalog_url(SirilWorldCS *center, double mag_limit, double dfov, int type);

gchar *fetch_url(const gchar *url);
void free_fetch_result(gchar *result);

int read_projected_catalog(GInputStream *stream, psf_star **cstars, object_catalog cat);

gpointer search_in_online_conesearch(gpointer p);
gpointer catsearch_worker(gpointer p);

int load_catalog(GFile *catalog_file, gboolean phot, psf_star **ret_stars, int *ret_nb_stars);
siril_catalog *siril_catalog_load_from_file(const gchar *filename, gboolean phot);
int siril_catalog_project_with_WCS(siril_catalog *siril_cat, fits *fit, gboolean use_proper_motion);
uint32_t siril_catalog_colums(object_catalog cat);
void sort_cat_items_by_mag(siril_catalog *siril_cat);


#endif
