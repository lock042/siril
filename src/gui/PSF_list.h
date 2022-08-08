#ifndef FWHM_LIST_H_
#define FWHM_LIST_H_

#include "core/siril.h"
#include "algos/PSF.h"

void refresh_star_list(psf_star **);
void clear_stars_list(gboolean refresh_GUI);
void pick_a_star();
int save_list(gchar *filename, gboolean forcepx, int max_stars_fitted, psf_star **stars, int nbstars, star_finder_params *sf, gboolean verbose);

void popup_psf_result(psf_star *result, rectangle *area);

#endif
