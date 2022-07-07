#ifndef FINDER_H_
#define FINDER_H_

#include "core/siril.h"

#define SF_ERRMSG_LEN 256

typedef struct {
	fits *fit;
	sequence *from_seq;
	int index_in_seq;
} image;

struct starfinder_data {
	image im;
	int layer;
	int max_stars_fitted;
};

struct star_candidate_struct {
	int x, y;
	float mag_est;
	float bg;
	float B, sx, sy;
	int R;
};
typedef struct star_candidate_struct starc;

// criteria above 10 are mandatory
typedef enum {
	SF_OK = 0,
	SF_FWHM_TOO_LARGE = 2,
	SF_RMSE_TOO_LARGE = 3,
	SF_CENTER_OFF = 10,
	SF_NO_FWHM = 11,
	SF_NO_POS = 12,
	SF_NO_MAG = 13,
	SF_FWHM_TOO_SMALL = 14,
	SF_FWHM_NEG = 15,
	SF_ROUNDNESS_BELOW_CRIT = 16
} sf_errors;

void update_peaker_GUI();
void confirm_peaker_GUI();
psf_star **peaker(image *image, int layer, star_finder_params *sf, int *nb_stars, rectangle *area, gboolean showtime, gboolean limit_nbstars, int maxstars, int threads);
psf_star *add_star(fits *fit, int layer, int *index);
int remove_star(int index);
void sort_stars(psf_star **stars, int total);
psf_star **new_fitted_stars(size_t n);
void free_fitted_stars(psf_star **stars);
int count_stars(psf_star **stars);
void FWHM_average(psf_star **stars, int nb, float *FWHMx, float *FWHMy, char **units, float *B);
float filtered_FWHM_average(psf_star **stars, int nb);
gpointer findstar(gpointer p);

#endif
