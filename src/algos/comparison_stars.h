#ifndef _COMPARISON_STARS_H
#define _COMPARISON_STARS_H

#include <glib.h>
#include "io/remote_catalogues.h"
#include "algos/PSF.h"
#include "algos/photometry.h"
#include "io/siril_catalogues.h"

struct compstars_arg {
	// feature input
	fits *fit; // the fits image
	gchar *target_name;	// star name to be looked-up on online sources
	gboolean narrow_fov;	// limit to image height as diameter
	siril_cat_index cat;	// the catalogue where comparison stars will be queried
	double delta_Vmag, delta_BV, max_emag;	// comparison stars filtering criteria
	gchar *nina_file;	// optional output NINA-type file name

	// for internal use
	cat_item *target_star;	// the considered variable star
	siril_catalogue *cat_stars;	// the list of stars for the field
	siril_catalogue *comp_stars;	// the list of photometric comparison stars
	GList *var_stars_cat; // a GList holding one or more variable star catalogues to discard
	
	int nb_comp_stars;
	gchar *AAVSO_chartid;
	gchar *AAVSO_uri;
	gboolean has_GUI;
	int retval;
};

typedef struct {
	int index;
	double dist;
} compstar_dist;

struct compstars_arg* init_compstars_arg();
void free_compstars_arg(gpointer p);

gpointer compstars_worker(gpointer arg);

int parse_nina_stars_file_using_WCS(struct light_curve_args *args, const char *file_path,
		gboolean use_comp1, gboolean use_comp2, fits *first);

void write_nina_file(struct compstars_arg *args);

void fill_compstar_item(cat_item *item, double ra, double dec, float mag, gchar *name, const gchar *type);

gchar *generate_lc_subtitle(struct compstars_arg *metadata, gboolean for_plot);

#endif
