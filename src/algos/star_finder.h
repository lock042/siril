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

typedef enum {
	SF_OK = 0,
	SF_NO_FWHM = 1,
	SF_NO_POS = 2,
	SF_NO_MAG = 3,
	SF_CENTER_OFF = 4,
	SF_FWHM_TOO_LARGE = 5,
	SF_FWHM_TOO_SMALL = 6,
	SF_FWHM_NEG = 7,
	SF_ROUNDNESS_BELOW_CRIT = 8,
	SF_RMSE_TOO_LARGE = 9,
} sf_errors;

void init_peaker_GUI();
void init_peaker_default();
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
