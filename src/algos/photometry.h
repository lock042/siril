#ifndef SRC_ALGOS_PHOTOMETRY_H_
#define SRC_ALGOS_PHOTOMETRY_H_

#include <glib.h>
#include <gsl/gsl_matrix.h>
#include "core/siril.h"

struct photometry_struct {
	double mag; // magnitude
	double s_mag; // magnitude uncertainty
	gboolean valid; // TRUE if no pixel outside of the range
	double SNR; // SNR estimation
};
typedef struct photometry_struct photometry;

typedef struct {
	float x, y;// in image pixels coordinates
	float mag; // visible magnitude (V filter), for sorting and debug
	float BV;
} pcc_star;

typedef enum {
	PSF_NO_ERR = 0,
	PSF_ERR_ALLOC = 3,
	PSF_ERR_UNSUPPORTED = 4,
	PSF_ERR_DIVERGED = 5,
	PSF_ERR_OUT_OF_WINDOW = 6,
	PSF_ERR_INNER_TOO_SMALL = 7,
	PSF_ERR_APERTURE_TOO_SMALL = 8,
	PSF_ERR_TOO_FEW_BG_PIX = 9,
	PSF_ERR_MEAN_FAILED = 10,
	PSF_ERR_INVALID_STD_ERROR = 11,
	PSF_ERR_INVALID_PIX_VALUE = 12,
	PSF_ERR_WINDOW_TOO_SMALL = 13,
	PSF_ERR_INVALID_IMAGE = 14,
	PSF_ERR_OUT_OF_IMAGE = 15,
	PSF_ERR_MAX_VALUE = 16	// keep last
} psf_error;

double get_camera_gain(fits *fit);

photometry *getPhotometryData(gsl_matrix* z, psf_star *psf, double gain,
		gboolean force_radius, gboolean verbose, psf_error *error);

void initialize_photometric_param();

void print_psf_error_summary(gint *code_sums);

#endif /* SRC_ALGOS_PHOTOMETRY_H_ */
