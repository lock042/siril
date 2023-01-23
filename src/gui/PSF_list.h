#ifndef FWHM_LIST_H_
#define FWHM_LIST_H_

#include "core/siril.h"
#include "algos/PSF.h"

void refresh_star_list(psf_star **);
void clear_stars_list(gboolean refresh_GUI);
void pick_a_star();
void set_iter_of_clicked_psf(double x, double y);

void popup_psf_result(psf_star *result, rectangle *area, fits *fit);

#endif
