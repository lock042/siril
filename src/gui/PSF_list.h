#ifndef FWHM_LIST_H_
#define FWHM_LIST_H_

#include "core/siril.h"
#include "algos/PSF.h"

void refresh_star_list();
/* clear_stars_list and clear_stars_list_as_idle are in algos/PSF.h */
void clear_psf_list_display(void); /* GTK star list widget reset, called from impl_clear_star_list */
void update_star_list(psf_star **new_stars, gboolean update_PSF_list, gboolean wait_for_update);
void pick_a_star();
void set_iter_of_clicked_psf(double x, double y);

void popup_psf_result(psf_star *result, rectangle *area, fits *fit);

int get_ra_and_dec_from_star_pos(psf_star *star, gdouble *alpha, gdouble *delta);

#endif
