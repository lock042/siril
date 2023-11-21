/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_ALGOS_SIRIL_WCS_H_
#define SRC_ALGOS_SIRIL_WCS_H_

#include <wcslib.h>
#include <wcsfix.h>

/* we force naxis to 2 */
#define NAXIS 2

typedef struct wcsprm wcsprm_t;

gboolean has_wcs(fits *fit);
gboolean has_wcsdata(fits *fit);
void free_wcs(fits *fit);
wcsprm_t *wcs_deepcopy(wcsprm_t *wcssrc, int *status);
wcsprm_t *load_WCS_from_hdr(char *header, int nkeyrec);
gboolean load_WCS_from_fits(fits* fit);
// this one directly uses the WCSLIB struct
void pix2wcs2(wcsprm_t *wcslib, double x, double y, double *r, double *d);
void pix2wcs(fits *fit, double pixel_x, double pixel_y, double *world_x, double *world_y);
int wcs2pix(fits *fit, double world_x, double world_y, double *pixel_x, double *pixel_y);
int *wcs2pix_array(fits *fit, int n, double *world, double *x, double *y);
void center2wcs(fits *fit, double *r, double *d);
double get_wcs_image_resolution(fits *fit);

void wcs_cdelt2unity(wcsprm_t *prm);
void wcs_pc2mat(wcsprm_t *prm, double pc[NAXIS][NAXIS]);
void wcs_cd2mat(wcsprm_t *prm, double cd[NAXIS][NAXIS]);
void wcs_mat2pc(wcsprm_t *prm, double pc[NAXIS][NAXIS]);
void wcs_mat2cd(wcsprm_t *prm, double cd[NAXIS][NAXIS]);
void wcs_print(wcsprm_t *prm);

#endif /* SRC_ALGOS_SIRIL_WCS_H_ */
