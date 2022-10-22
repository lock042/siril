#ifndef PSF_H_
#define PSF_H_

#include <gsl/gsl_matrix.h>
#include "algos/photometry.h"
#include "algos/star_finder.h"

#define _2_SQRT_2_LOG2 2.35482004503
#define INV_4_LOG2 0.360673760222241

//in siril.h: typedef struct fwhm_struct psf_star;

struct fwhm_struct {
	double B; /* average sky background value */
	double A; /* amplitude */
	double x0, y0; /* coordinates of the peak */
	double sx, sy; /* Size of the fitted function on the x and y axis in PSF coordinates */
	double fwhmx, fwhmy; /* FWHM in x and y axis */
	double fwhmx_arcsec, fwhmy_arcsec; /* FWHM in x and y axis in arc second */
	double angle; /* angle of the axis x,y with respect to the image's */
	double rmse; /* RMSE of the minimization */
	double sat; /* Level above which pixels have satured (defined by peaker) */
	int R; /* Optimized box sixe to enclose sufficient pixels in the background */
	gboolean has_saturated;

	// Moffat parameters
	double beta; /* Moffat equation beta parameter */
	starprofile profile; // Whether the profile is Gaussian or Moffat with beta {free|fixed}
	// The other parameters B, A, x0, y0, angle, rmse, sat... are the same as for Gaussian

	double xpos, ypos; /* position of the star in the image, not set by Minimization */

	/* photometry data - mag, s_mag and SNR are copied from phot if phot_is_valid,
	 * otherwise mag is approximate (computed from G and B) and the others are not set */
	double mag;	/* magnitude, approximate or accurate depending on phot_is_valid */
	double s_mag;	/* error on the magnitude, defaults to 9.9999 */
	double SNR;	/* SNR of the star, defaults to 0 */
	photometry *phot; /* photometry data */
	gboolean phot_is_valid; /* valid if computed by photometry and no saturated pixel detected */
	double BV; /* only used to pass data in photometric color calibration */

	/* uncertainties */
	double B_err;
	double A_err;
	double x_err, y_err;
	double sx_err, sy_err;
	double ang_err;
	double beta_err;

	int layer;
	char* units;
};

struct PSF_data {
	size_t n;
	double *y;
	size_t NbRows;
	size_t NbCols;
	double rmse;
	gboolean *mask;
	gboolean betafree; // not used untill we implement PSF_MOFFAT_BFIXED
	double beta; //  not used untill we implement PSF_MOFFAT_BFIXED
};

double psf_get_fwhm(fits *fit, int layer, rectangle *selection, double *roundness);

psf_star *psf_get_minimisation(fits *fit, int layer, rectangle *area,
		gboolean for_photometry, struct phot_config *phot_set, gboolean verbose,
		starprofile profile, psf_error *error);

psf_star *psf_global_minimisation(gsl_matrix* z, double bg, double sat, int convergence, 
		gboolean from_peaker, gboolean for_photometry, struct phot_config *phot_set, gboolean verbose,
		starprofile profile, psf_error *error);

void psf_display_result(psf_star *, rectangle *);
void fwhm_to_arcsec_if_needed(fits*, psf_star*);
gboolean get_fwhm_as_arcsec_if_possible(psf_star *star, double *fwhmx, double *fwhmy, char **unit);
gboolean convert_single_fwhm_to_arcsec_if_possible(double fwhm, double bin, double px_size, double flength, double *result);

psf_star *new_psf_star();
psf_star *duplicate_psf(psf_star *);
void free_psf(psf_star *psf);

#endif
