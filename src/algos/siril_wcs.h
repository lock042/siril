/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
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
#define MAX_SIP_ORDER 6 // see also MAX_DISTO_SIZE when changing this value
#define MAX_SIP_SIZE MAX_SIP_ORDER + 1

typedef struct wcsprm wcsprm_t;

gboolean has_wcs(fits *fit);
gboolean has_wcsdata(fits *fit);
void reset_wcsdata(fits *fit);
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
void center2wcs2(struct wcsprm *wcs, int width, int height, double *r, double *d);
double get_wcs_image_resolution(fits *fit);

void wcs_pc2mat(wcsprm_t *prm, double pc[NAXIS][NAXIS]);
void wcs_cd2mat(wcsprm_t *prm, double cd[NAXIS][NAXIS]);
void wcs_mat2pc(wcsprm_t *prm, double pc[NAXIS][NAXIS]);
void wcs_mat2cd(wcsprm_t *prm, double cd[NAXIS][NAXIS]);
void wcs_mat2cdelt(wcsprm_t *prm, double cdelt[NAXIS]);
void wcs_decompose_cd(wcsprm_t *prm, double cd[NAXIS][NAXIS]);
gboolean image_is_flipped_from_wcs(struct wcsprm *wcslib);

int extract_SIP_order_and_matrices(struct disprm *dis, 
		double A[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double B[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double AP[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double BP[MAX_SIP_SIZE][MAX_SIP_SIZE]);
void update_SIP_keys(struct disprm *dis, 
		double A[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double B[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double AP[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double BP[MAX_SIP_SIZE][MAX_SIP_SIZE]);
void wcs_print(wcsprm_t *prm);

void remove_dis_from_wcs(wcsprm_t *prm);
void create_wcs(double ra0, double dec0, double scale, double framing_angle, int rx, int ry, struct wcsprm *prm);

#endif /* SRC_ALGOS_SIRIL_WCS_H_ */
