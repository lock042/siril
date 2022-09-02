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
	gboolean save_to_file;	// generate starfile if TRUE and in sequence
	gchar *starfile;	// save to file if not NULL
	psf_star ***stars;	// save to pointer if not NULL
	int *nb_stars;		// number of stars in stars if not NULL
	threading_type threading;
	gboolean update_GUI;	// FALSE for sequence operation
	gboolean process_all_images;	// for sequence operation
	gboolean already_in_thread;
};

struct star_candidate_struct {
	int x, y;
	float mag_est;
	float sat;
	float sx, sy;
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
void sort_stars_by_mag(psf_star **stars, int total);
psf_star **new_fitted_stars(size_t n);
void free_fitted_stars(psf_star **stars);
int count_stars(psf_star **stars);
void FWHM_stats(psf_star **stars, int nb, int bitpix, float *FWHMx, float *FWHMy, char **units, float *B, float *Acut, double Acutp) ;
psf_star **filter_stars_by_amplitude(psf_star **stars, float threshold, int *nbfilteredstars);
float filtered_FWHM_average(psf_star **stars, int nb);
int apply_findstar_to_sequence(struct starfinder_data *findstar_args);
gpointer findstar_worker(gpointer p);

#endif
