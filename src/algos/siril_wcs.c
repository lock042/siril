/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "core/siril.h"
#include "algos/astrometry_solver.h"
#include "io/image_format_fits.h"

#include "siril_wcs.h"


// Use this flag to print wcslib related verbose - not for production
// #define DEBUG_WCS

gboolean has_wcs(fits *fit) {
	return fit->keywords.wcslib != NULL;
}

// deal with cases where wcsdata is not NULL but members are set to 0
gboolean has_wcsdata(fits *fit) {
	return fit->keywords.wcsdata.pltsolvd_comment[0] != '\0';
}

void reset_wcsdata(fits *fit) {
	fit->keywords.wcsdata.pltsolvd = FALSE;
	memset(&fit->keywords.wcsdata.pltsolvd_comment, 0, sizeof(fit->keywords.wcsdata.pltsolvd_comment));
}


void free_wcs(fits *fit) {
	if (fit->keywords.wcslib) {
		if (!wcsfree(fit->keywords.wcslib))
			free(fit->keywords.wcslib);
		fit->keywords.wcslib = NULL;
	}
}

wcsprm_t *wcs_deepcopy(wcsprm_t *wcssrc, int *status) {
	if (status)
		*status = 1;
	if (!wcssrc)
		return NULL;
	wcsprm_t *wcsdst = NULL;
	int axes[2], nsub;
	nsub = 2;
	axes[0] = WCSSUB_LONGITUDE;
	axes[1] = WCSSUB_LATITUDE;
	wcsdst = calloc(1, sizeof(wcsprm_t));
	if (!wcsdst) {
		PRINT_ALLOC_ERR;
		if (status)
			*status = WCSERR_MEMORY;
		return NULL;
	}
	wcsdst->flag = -1;
	int statuscpy = wcscopy(0, wcssrc, wcsdst);
	if (statuscpy) {
		if (status)
			*status = statuscpy;
		wcsfree(wcsdst);
		return NULL;
	}
	if (status)
		*status = 0;
	return wcsdst;
}


wcsprm_t *load_WCS_from_hdr(char *header, int nkeyrec) {
	wcsprm_t *data = NULL, *wcs = NULL;
	int nreject, nwcs;
	int ctrl = 0;
#ifdef DEBUG_WCS
	ctrl = 2; // writes rejected WCS cards
#endif
	/** There was a bug with wcspih that it is not really thread-safe for wcslib version < 7.5.
	 * We now force to have 7.12 at least */
	int wcs_status = wcspih(header, nkeyrec, 0, ctrl, &nreject, &nwcs, &data);

	if (wcs_status == 0) {
		for (int i = 0; i < nwcs; i++) {
			/* Find the master celestial WCS coordinates */
			wcsprm_t *prm = data + i;
			wcsset(prm); // is it necessary?
			if (prm->lng >= 0 && prm->lat >= 0
					&& (prm->alt[0] == '\0' || prm->alt[0] == ' ')) {
				int status = -1;
				wcs = wcs_deepcopy(prm, &status);
				if (!status) {
					if (wcs->altlin & 2) { // header contains CD info
						double cd[2][2] = {{ 0. }};
						wcs_cd2mat(wcs, cd);
						wcs_decompose_cd(wcs, cd);
						wcs->altlin = 1;
						wcs->flag = 0;
						wcsset(wcs);
						siril_debug_print("contains CD\n");
					} else if (wcs->altlin & 1) { // header contains PC info
						double pc[2][2] = {{ 0. }}, cd[2][2] = {{ 0. }};
						wcs_pc2mat(wcs, pc);
						wcs_pc_to_cd(pc, wcs->cdelt, cd);
						wcs_mat2cd(wcs, cd);
						wcs->flag = 0;
						wcsset(wcs);
						// siril_debug_print("contains PC\n");
					} else { // contained some keywords but not enough to define at least a linear projection
						siril_debug_print("wcs did not contain enough info\n");
						free(wcs);
						wcs = NULL;
						break;
					}
					// siril_debug_print("at header readout\n");
					// wcs_print(wcs);
					break;
				} else {
					siril_debug_print("wcssub error %d: %s.\n", status, wcs_errmsg[status]);
					wcsfree(wcs);
					wcs = NULL;
				}
			}
		}
		wcsvfree(&nwcs, &data);
	}
	return wcs;
}


gboolean load_WCS_from_fits(fits* fit) {
	int status = 0;
	char *header;
	struct wcsprm *wcs = NULL;
	int nkeyrec;
	if (fit->keywords.wcslib) {
		free_wcs(fit);
		reset_wcsdata(fit);
	}
	ffhdr2str(fit->fptr, 1, NULL, 0, &header, &nkeyrec, &status);
	if (status) {
		report_fits_error(status);
		return FALSE;
	}

	wcs = load_WCS_from_hdr(header, nkeyrec);
	free(header);

	if (!wcs) {
		siril_debug_print("No world coordinate systems found.\n");
		wcsfree(wcs);
		return FALSE;
	}
	fit->keywords.wcslib = wcs;
	return TRUE;
}

void pix2wcs2(struct wcsprm *wcslib, double x, double y, double *r, double *d) {
	*r = 0.0;
	*d = 0.0;
	int status, stat[NWCSFIX];
	double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];

	// In WCS convention, origin of the grid is at (-0.5, -0.5) wrt siril grid
	pixcrd[0] = x + 0.5;
	pixcrd[1] = y + 0.5;

	status = wcsp2s(wcslib, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
	if (status != 0)
		return;

	*r = world[0];
	*d = world[1];
}

void pix2wcs(fits *fit, double x, double y, double *r, double *d) {
	*r = 0.0;
	*d = 0.0;
	if (fit->keywords.wcslib)
		pix2wcs2(fit->keywords.wcslib, x, y, r, d);
}

// ra in degrees
int wcs2pix(fits *fit, double ra, double dec, double *x, double *y) {
	if (x) *x = -1.0;
	if (y) *y = -1.0;
	int status, stat[NWCSFIX];
	double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];
	world[0] = ra;
	world[1] = dec;

	status = wcss2p(fit->keywords.wcslib, 1, 2, world, &phi, &theta, imgcrd, pixcrd, stat);

	if (!status) {
		double xx = pixcrd[0];
		double yy = pixcrd[1];
		if (isnan(xx) || isnan(yy))
			return WCSERR_NO_SOLUTION;
		// return values even if outside (required for celestial grid display)
		// In WCS convention, origin of the grid is at (-0.5, -0.5) wrt siril grid
		if (x) *x = xx - 0.5;
		if (y) *y = yy - 0.5;
		if (xx < 0.0 || yy < 0.0 || xx > (double)fit->rx || yy > (double)fit->ry) {
			//siril_debug_print("outside image but valid return\n");
			// wcss2p returns values between 0 and 9, picking a new one
			status = 10;
		}
	}
	return status;
}

// same as wcs2pix except it takes a world array as input
// world is an array with [ra1, dec1, ra2, dec2...ran, decn], i.e 2n elements (row major)
// it returns an allocated array of statuses (instead of a single status), which must be freed
int *wcs2pix_array(fits *fit, int n, double *world, double *x, double *y) {
	if (x) {
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
	if (y) {
		for (int i = 0; i < n; i++)
			y[i] = -1.0;
	}
	// can't pass NULL to the values we don't want to retrieve (intcrd, phi, theta)
	double *intcrd = malloc((2 * n) * sizeof(double));
	double *pixcrd = malloc((2 * n) * sizeof(double));
	double *phi = malloc(n * sizeof(double));
	double *theta = malloc(n * sizeof(double));
	int c = 0;
	int *status = calloc((unsigned)n , sizeof(int));
	int globstatus = wcss2p(fit->keywords.wcslib, n, 2, world, phi, theta, intcrd, pixcrd, status);
	if (globstatus == WCSERR_SUCCESS || globstatus == WCSERR_BAD_WORLD) {// we accept BAD_WORLD as it does not mean all of the conversions failed
		for (int i = 0; i < n; i++) {
			if (!status[i]) {
				double xx = pixcrd[c++];
				double yy = pixcrd[c++];
				// return values even if outside (required for celestial grid display)
				// In WCS convention, origin of the grid is at (-0.5, -0.5) wrt siril grid
				if (x) x[i] = xx - 0.5;
				if (y) y[i] = yy - 0.5;
				if (xx < 0.0 || yy < 0.0 || xx > (double)fit->rx || yy > (double)fit->ry) {
					//siril_debug_print("outside image but valid return\n");
					// wcss2p returns values between 0 and 9, picking a new one
					status[i] = 10;
				}
			} else {
				c += 2;
			}
		}
	} else {
		free(status);
		status = NULL;
	}
	free(intcrd);
	free(pixcrd);
	free(phi);
	free(theta);
	return status;
}


/* get image center celestial coordinates */
void center2wcs(fits *fit, double *r, double *d) {
	*r = -1.0;
	*d = -1.0;
	int status, stat[NWCSFIX];
	double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];

	// In WCS convention, origin of the grid is at (-0.5, -0.5) wrt siril grid
	pixcrd[0] = (double)(fit->rx) * 0.5 + 0.5;
	pixcrd[1] = (double)(fit->ry) * 0.5 + 0.5;

	status = wcsp2s(fit->keywords.wcslib, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
	if (status != 0)
		return;

	*r = world[0];
	*d = world[1];
}

/* get image center celestial coordinates - lower level from wcs and image width/height */
void center2wcs2(struct wcsprm *wcs, int width, int height, double *r, double *d) {
	*r = -1.0;
	*d = -1.0;
	int status, stat[NWCSFIX];
	double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];

	// In WCS convention, origin of the grid is at (-0.5, -0.5) wrt siril grid
	pixcrd[0] = (double)(width) * 0.5 + 0.5;
	pixcrd[1] = (double)(height) * 0.5 + 0.5;

	status = wcsp2s(wcs, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
	if (status != 0)
		return;

	*r = world[0];
	*d = world[1];
}

void wcs_pc2mat(wcsprm_t *prm, double pc[NAXIS][NAXIS]) {
	if (!prm || !prm->pc)
		return;
	double *pcij = prm->pc;
	for (int i = 0; i < NAXIS; i++) {
		for (int j = 0; j < NAXIS; j++) {
			pc[i][j] = *(pcij++);
		}
	}
}
void wcs_cd2mat(wcsprm_t *prm, double cd[NAXIS][NAXIS]) {
	if (!prm || !prm->cd)
		return;
	double *cdij = prm->cd;
	for (int i = 0; i < NAXIS; i++) {
		for (int j = 0; j < NAXIS; j++) {
			cd[i][j] = *(cdij++);
		}
	}
}

void wcs_mat2pc(wcsprm_t *prm, double pc[NAXIS][NAXIS]) {
	if (!prm || !prm->pc)
		return;
	double *pcij = prm->pc;
	for (int i = 0; i < NAXIS; i++) {
		for (int j = 0; j < NAXIS; j++) {
			*(pcij++) = pc[i][j];
		}
	}
}
void wcs_mat2cd(wcsprm_t *prm, double cd[NAXIS][NAXIS]) {
	if (!prm || !prm->cd)
		return;
	double *cdij = prm->cd;
	for (int i = 0; i < NAXIS; i++) {
		for (int j = 0; j < NAXIS; j++) {
			*(cdij++) = cd[i][j];
		}
	}
}

void wcs_mat2cdelt(wcsprm_t *prm, double cdelt[NAXIS]) {
	if (!prm || !prm->cdelt)
		return;
	double *cdelti = prm->cdelt;
	for (int i = 0; i < NAXIS; i++) {
		*(cdelti++) = cdelt[i];
	}
}

void wcs_decompose_cd(wcsprm_t *prm, double cd[NAXIS][NAXIS]) {
	if (!prm)
		return;
	double cdelt[NAXIS];
	double pc[NAXIS][NAXIS];
	cdelt[0] = -sqrt(cd[0][0] * cd[0][0] + cd[0][1] * cd[0][1]); // by convention cdelt1 is < 0
	cdelt[1] =  sqrt(cd[1][0] * cd[1][0] + cd[1][1] * cd[1][1]); // by convention cdelt2 is > 0
	pc[0][0] = cd[0][0] / cdelt[0];
	pc[0][1] = cd[0][1] / cdelt[0];
	pc[1][0] = cd[1][0] / cdelt[1];
	pc[1][1] = cd[1][1] / cdelt[1];
	wcs_mat2cd(prm, cd);
	wcs_mat2pc(prm, pc);
	wcs_mat2cdelt(prm, cdelt);
}

gboolean image_is_flipped_from_wcs(struct wcsprm *wcslib) {
	double cd[2][2] = {{ 0. }};
	wcs_cd2mat(wcslib, cd);
	double det = (cd[0][0] * cd[1][1] - cd[1][0] * cd[0][1]); // determinant of rotation matrix (ad - bc)
	return det > 0; // convention is that angles are positive clockwise when image is not flipped
}

/* get resolution in degree/pixel */
double get_wcs_image_resolution(fits *fit) {
	double resolution = -1.0;
	if (fit->keywords.wcslib) {
		resolution = (fabs(fit->keywords.wcslib->cdelt[0]) + fabs(fit->keywords.wcslib->cdelt[1])) * 0.5;
	}
	if (resolution <= 0.0) {
		if (fit->keywords.focal_length >= 0.0 && fit->keywords.pixel_size_x >= 0.0 && fit->keywords.pixel_size_y == fit->keywords.pixel_size_x)
			resolution = (RADCONV / fit->keywords.focal_length * fit->keywords.pixel_size_x) / 3600.0;
		// what about pix size x != y?
	}
	return resolution;
}

// return the order of the SIP polynomials and fills the coeffs matrices (if first matrix A is not NULL)
int extract_SIP_order_and_matrices(struct disprm *dis,
		double A[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double B[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double AP[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double BP[MAX_SIP_SIZE][MAX_SIP_SIZE]) {
	int order = 0;
	if (!dis)
		return 0;
	for (int n = 0; n < dis->ndp; n++) {
		if (!g_str_has_prefix(dis->dp[n].field + 4, "SIP"))
			continue;
		int mat  = dis->dp[n].j; // if 1, it's the A terms, if 2, it's the B terms
		int fwd = (g_str_has_prefix(dis->dp[n].field + 8, "FWD")) ? 1 : 2; // if 1, it's A/B, if 2, it's AP/BP
		int i, j;
		sscanf(dis->dp[n].field + 12, "%d_%d", &i, &j);
		if (fabs(dis->dp[n].value.f) > 0) {
			order = max(order, i);
			order = max(order, j);
		}
		if (!A) // Warning: for brevity, we only test for the first matrix for brevity...
			continue;
		if (mat == 1) {
			if (fwd == 1) {
				A[i][j] = dis->dp[n].value.f;
			} else {
				AP[i][j] = dis->dp[n].value.f;
			}
		} else {
			if (fwd == 1) {
				B[i][j] = dis->dp[n].value.f;
			} else {
				BP[i][j] = dis->dp[n].value.f;
			}
		}
	}
	return order;
}

void update_SIP_keys(struct disprm *dis,
		double A[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double B[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double AP[MAX_SIP_SIZE][MAX_SIP_SIZE],
		double BP[MAX_SIP_SIZE][MAX_SIP_SIZE]) {
	if (!dis)
		return;
	for (int n = 0; n < dis->ndp; n++) {
		if (!g_str_has_prefix(dis->dp[n].field + 4, "SIP"))
			continue;
		int mat  = dis->dp[n].j; // if 1, it's the A terms, if 2, it's the B terms
		int fwd = (g_str_has_prefix(dis->dp[n].field + 8, "FWD")) ? 1 : 2; // if 1, it's A/B, if 2, it's AP/BP
		int i, j;
		sscanf(dis->dp[n].field + 12, "%d_%d", &i, &j);
		if (mat == 1) {
			dis->dp[n].value.f = (fwd == 1) ? A[i][j] : AP[i][j];
		} else {
			dis->dp[n].value.f = (fwd == 1) ? B[i][j] : BP[i][j];
		}
	}
}

void wcs_print(wcsprm_t *prm) {
#ifdef DEBUG_WCS
	printf("CRPIX\n");
	int c = 0;
	for (int i = 0; i < NAXIS; i++) {
			printf("%g ", prm->crpix[c++]);
	}
	printf("\n");
	printf("CRVAL\n");
	c = 0;
	for (int i = 0; i < NAXIS; i++) {
			printf("%g ", prm->crval[c++]);
	}
	printf("\n");
	printf("PC\n");
	c = 0;
	for (int i = 0; i < NAXIS; i++) {
		for (int j = 0; j < NAXIS; j++) {
			printf("%g ", prm->pc[c++]);
		}
		printf("\n");
	}
	printf("CDELT\n");
	c = 0;
	for (int i = 0; i < NAXIS; i++) {
			printf("%g ", prm->cdelt[c++]);
	}
	printf("\n");
	printf("CD\n");
	c = 0;
	for (int i = 0; i < NAXIS; i++) {
		for (int j = 0; j < NAXIS; j++) {
			printf("%g ", prm->cd[c++]);
		}
		printf("\n");
	}
	printf("\n");
	if (prm->lin.dispre) {
		struct disprm *dis = prm->lin.dispre;
		disset(dis);
		for (int j = 0; j < dis->ndp; j++) {
			printf("%s %d", dis->dp[j].field, dis->dp[j].j);
			if (!dis->dp[j].type) //int
				printf(" %d\n", dis->dp[j].value.i);
			else //float
				printf(" %g\n", dis->dp[j].value.f);
		}
	}
#endif
}

void remove_dis_from_wcs(wcsprm_t *prm) {
	disfree(prm->lin.dispre);
	prm->lin.dispre = NULL;
	prm->flag = 0;
	wcsset(prm);
}

void create_wcs(double ra0, double dec0, double scale, double framing_angle, int rx, int ry, struct wcsprm *prm) {
	// preparing pc matrix
	double ca, sa;
	sincosd(framing_angle, &sa, &ca);
	double pc[2][2];
	pc[0][0] = ca;
	pc[1][1] = ca;
	pc[0][1] = -sa;
	pc[1][0] = sa;

	/**** Fill wcslib structure ***/
	prm->flag = -1;
	wcsinit(1, 2, prm, 0, 0, 0);
	prm->equinox = 2000.0;
	prm->crpix[0] = (double)rx * 0.5 + 0.5;
	prm->crpix[1] = (double)ry * 0.5 + 0.5;
	prm->crval[0] = ra0;
	prm->crval[1] = dec0;
	prm->cdelt[0] = -scale;
	prm->cdelt[1] =  scale;
	prm->lonpole = 180.;
	wcs_mat2pc(prm, pc);

	const char CTYPE[2][9] = { "RA---TAN", "DEC--TAN" };
	const char CUNIT[2][4] = { "deg", "deg" };
	for (int i = 0; i < NAXIS; i++) {
		strncpy(prm->cunit[i], &CUNIT[i][0], 71); // 72 char fixed buffer, keep 1 for the NULL
	}
	for (int i = 0; i < NAXIS; i++) {
		strncpy(prm->ctype[i], &CTYPE[i][0], 71); // 72 byte buffer, leave 1 byte for the NULL
	}
	wcsset(prm);
}
