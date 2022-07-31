#ifndef SRC_ALGOS_PHOTOMETRY_H_
#define SRC_ALGOS_PHOTOMETRY_H_

#include <glib.h>
#include <gsl/gsl_matrix.h>
#include "core/siril.h"
#include "core/settings.h"

typedef struct {
	double mag; // magnitude
	double s_mag; // magnitude uncertainty
	gboolean valid; // TRUE if no pixel outside of the range
	double SNR; // SNR estimation
} photometry;

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

struct phot_config *phot_set_adjusted_for_image(fits *fit);

rectangle compute_dynamic_area_for_psf(psf_star *psf, struct phot_config *original, struct phot_config *phot_set, Homography H, Homography Href);

photometry *getPhotometryData(gsl_matrix* z, psf_star *psf,
		struct phot_config *phot_set, gboolean verbose, psf_error *error);

void initialize_photometric_param();

const char *psf_error_to_string(psf_error err);
void print_psf_error_summary(gint *code_sums);

int new_light_curve(sequence *seq, const char *filename, const char *target_descr, gboolean display_graph);

gpointer crazy_photo_worker(gpointer arg);

#endif /* SRC_ALGOS_PHOTOMETRY_H_ */
