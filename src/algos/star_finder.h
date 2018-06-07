#ifndef FINDER_H_
#define FINDER_H_

void init_peaker_GUI();
void init_peaker_default();
fitted_PSF **peaker(fits *fit, int layer, star_finder_params *sf, int *nb_stars, rectangle *area, gboolean showtime);
fitted_PSF *add_star(fits *fit, int layer, int *index);
int remove_star(int index);
void sort_stars(fitted_PSF **stars, int total);
int count_stars(fitted_PSF **stars);
void FWHM_average(fitted_PSF **stars, float *FWHMx, float *FWHMy, int max);

#endif
