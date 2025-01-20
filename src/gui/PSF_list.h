#ifndef FWHM_LIST_H_
#define FWHM_LIST_H_

#include "core/siril.h"
#include "algos/PSF.h"

void refresh_star_list();
void clear_stars_list(gboolean refresh_GUI);
void update_star_list(psf_star **new_stars, gboolean update_PSF_list, gboolean wait_for_update);
void pick_a_star();
void set_iter_of_clicked_psf(double x, double y);

void popup_psf_result(psf_star *result, rectangle *area, fits *fit);

int get_ra_and_dec_from_star_pos(psf_star *star, gdouble *alpha, gdouble *delta);

#endif
