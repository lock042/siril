#ifndef SRC_ALGOS_PHOTOMETRY_H_
#define SRC_ALGOS_PHOTOMETRY_H_

#include <glib.h>
#include <gsl/gsl_matrix.h>
#include "core/siril.h"
#include "core/settings.h"
#include "algos/astrometry_solver.h"
#include "algos/PSF.h"
#include "io/siril_plot.h"

struct photometry_struct {
	double mag;	// magnitude
	double s_mag;	// magnitude uncertainty
	gboolean valid;	// TRUE if no pixel outside of the range
	double SNR;	// SNR estimation
};

struct phot_config *phot_set_adjusted_for_image(const fits *fit);

photometry *getPhotometryData(gsl_matrix* z, const psf_star *psf,
		struct phot_config *phot_set, gboolean verbose, psf_error *error);

void initialize_photometric_param();

const char *psf_error_to_string(psf_error err);
void print_psf_error_summary(gint *code_sums);

/* light curves */
/*struct light_curve_metadata {
	double delta_Vmag, delta_BV;
	int nb_comp_stars;
	gchar *AAVSO_chartid;
	gchar *AAVSO_uri;
};*/

struct compstars_arg;

struct light_curve_args {
	rectangle *areas;	// the first is the variable star's area
	int nb;			// number of areas
	sequence *seq;
	int layer;
	char *target_descr;	// the description to put in the data file and graph
	gboolean display_graph;	// if true, show it, if false, generate png
	gboolean force_rad;	// if true, fixed aperture radius, if false, dynamic aperture radius

	// metadata from the NINA file created by Siril
	struct compstars_arg *metadata;

	// spl_data for siril_plot if the light curve is displayed
	siril_plot_data *spl_data;
};

void free_light_curve_args(struct light_curve_args *args);

gpointer light_curve_worker(gpointer arg);

int new_light_curve(const char *filename, struct light_curve_args *lcargs);


struct catmag_data {
	siril_cat_index catalogue;
	gboolean limit_BV;	// NOMAD
	float refBV, dBV;
	gboolean limit_temperature; // Gaia
	float refT, dT;
	fits *fit;
};
gpointer catmag_mono_worker(gpointer arg);

#endif /* SRC_ALGOS_PHOTOMETRY_H_ */
