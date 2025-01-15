#ifndef SRC_ALGOS_PHOTOMETRIC_CC_H_
#define SRC_ALGOS_PHOTOMETRIC_CC_H_

#include <stdio.h>
#include <glib.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/PSF.h"
#include "algos/photometry.h"
#include "algos/astrometry_solver.h"

__inline __attribute__((always_inline)) int xisnanf(float x) { return x != x; }

typedef enum {
	CHANNEL_HIGHEST,
	CHANNEL_MIDDLE,
	CHANNEL_LOWEST
} normalization_channel;

struct photometric_cc_data {
	fits *fit;			// the image to process
	gboolean bg_auto;		// automatically select an area for bkg neutralization
	rectangle bg_area;		// the area for background if not bg_auto
	siril_catalogue *ref_stars;
	siril_cat_index catalog;		// catalog used for photometry
	limit_mag_mode mag_mode;	// automatically limit magnitude of the catalog
	double magnitude_arg;		// if not automatic, use this limit magnitude
	gboolean atmos_corr;		// whether to correct for wavelength dependent atmospheric effects
	double atmos_obs_height;	// observer height for atmosphere model
	double atmos_pressure;		// atmospheric pressure
	gboolean atmos_pressure_is_slp;	// is atmospheric pressure sea level pressure?
	pcc_star *stars;		// the list of stars with BV index in the image
	int nb_stars;			// the number of stars in the array
	float fwhm;			// representative FWHM for stars
	float t0; // lower background tolerance, in sigma units
	float t1; // upper background tolerance, in sigma units
	gchar *datalink_path;	// to hold the datalink path for SPCC
	gboolean spcc;			// set if doing SPCC
	gboolean spcc_mono_sensor; // for SPCC
	gboolean is_dslr; // for SPCC
	int selected_sensor_osc; // for SPCC
	int selected_sensor_m; // for SPCC
	int selected_filter_osc; // for SPCC
	int selected_filter_r; // for SPCC
	int selected_filter_g; // for SPCC
	int selected_filter_b; // for SPCC
	int selected_filter_lpf; // for SPCC
	int selected_white_ref; // for SPCC
	gboolean do_plot; // for SPCC
	gboolean nb_mode; // for SPCC
	double nb_center[3]; // for SPCC
	double nb_bandwidth[3]; // for SPCC
	cmsCIExyYTRIPLE primaries; // used for SPCC source profile
};

int apply_photometric_color_correction(fits *fit, const float *kw, const float *bg);
int get_stats_coefficients(fits *fit, rectangle *area, float *bg, float t0, float t1);
int photometric_cc(struct photometric_cc_data *args);
gpointer photometric_cc_standalone(gpointer p);
pcc_star *convert_siril_cat_to_pcc_stars(siril_catalogue *siril_cat, int *nbstars);
int get_favourite_spccobject(GList *list, const gchar *favourite);
int get_favourite_oscsensor(GList *list, const gchar *favourite);
int make_selection_around_a_star(cat_item *star, rectangle *area, fits *fit);

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */

