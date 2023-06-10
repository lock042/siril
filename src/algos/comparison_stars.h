#ifndef _COMPARISON_STARS_H
#define _COMPARISON_STARS_H

#include <glib.h>
#include "io/remote_catalogues.h"
#include "algos/PSF.h"
#include "algos/photometry.h"

struct compstars_arg {
	// feature input
	gchar *target_name;	// star name to be looked-up on online sources
	gboolean narrow_fov;	// limit to image height as diameter
	online_catalog cat;	// the catalogue where comparison stars will be queried
	double delta_Vmag, delta_BV;	// comparison stars filtering criteria
	gchar *nina_file;	// optional output NINA-type file name

	// for internal use
	psf_star *target_star;	// the considered variable star
	psf_star *cat_stars;	// the list of stars for the field
	int nb_cat_stars;
	psf_star **comp_stars;	// the list of photometric comparison stars
	int nb_comp_stars;
	gchar *AAVSO_chartid;
	gchar *AAVSO_uri;
	gboolean has_GUI;
	int retval;
};

gpointer compstars_worker(gpointer arg);

int parse_nina_stars_file_using_WCS(struct light_curve_args *args, const char *file_path,
		gboolean use_comp1, gboolean use_comp2, fits *first);

gchar *generate_lc_subtitle(struct compstars_arg *metadata, gboolean for_gnuplot);

#endif
