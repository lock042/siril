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

typedef struct struct_coeff {
	float value;
	int channel;
} coeff;

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

	pcc_star *stars;		// the list of stars with BV index in the image
	int nb_stars;			// the number of stars in the array
	int max_spcc_stars;			// the maximum number of stars to request from Gaia DR3
	float fwhm;			// representative FWHM for stars
	gchar *datalink_path;	// to hold the datalink path for SPCC
	gboolean spcc;			// set if doing SPCC
	gboolean spcc_mono_sensor; // for SPCC
	int selected_sensor_osc; // for SPCC
	int selected_sensor_m; // for SPCC
	int selected_filter_osc; // for SPCC
	int selected_filter_r; // for SPCC
	int selected_filter_g; // for SPCC
	int selected_filter_b; // for SPCC
	gboolean use_osc_filter; // for SPCC
	cmsCIExyYTRIPLE primaries; // used for SPCC source profile
	cmsCIExyY whitepoint; // used for SPCC source profile
};

int apply_photometric_color_correction(fits *fit, const float *kw, const coeff *bg, const float *mins, const float *maxs, int norm_channel);
int get_stats_coefficients(fits *fit, rectangle *area, coeff *bg, float *mins, float *maxs, int *norm_channel);
int photometric_cc(struct photometric_cc_data *args);
gpointer photometric_cc_standalone(gpointer p);
pcc_star *convert_siril_cat_to_pcc_stars(siril_catalogue *siril_cat, int *nbstars);

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */

