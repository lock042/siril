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

typedef struct {
	float x, y;	// in FITS/WCS coordinates
	float mag;	// visible magnitude (V filter), for sorting and debug
	float BV;	// B magnitude - V magnitude, -99.9 if not available
	float teff; // Gaia Teff
	uint64_t index; // Order in the Gaia results table. This is used to match
					// HDUs in the FITS, in case of excluded stars
} pcc_star;

struct phot_config *phot_set_adjusted_for_image(const fits *fit);

void fluxCut_factors (const psf_star *psf, double beta_ref, double fwhm_ref, double* in_rad, double* out_rad, double* ap_rad);
struct radii_set *radii_strat (struct phot_config *phot_set, const psf_star *psf);

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

struct radii_set {
	double ape_Radius;
	double in_Radius;
	double out_Radius;
};

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


#endif /* SRC_ALGOS_PHOTOMETRY_H_ */
