/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_WCSLIB
#include <wcslib.h>
#endif

#include "core/siril.h"
#include "io/image_format_fits.h"

#include "siril_wcs.h"

static struct wcsprm *wcs = NULL;

gboolean load_WCS_from_memory(fits *fit) {
#ifdef HAVE_WCSLIB
	if (wcs) {
		free_wcs();
	}
	wcs = calloc(1, sizeof(struct wcsprm));
	wcs->flag = -1;
	wcsinit(1, fit->naxis, wcs, -1, -1, -1);
	wcs->cd[0] = fit->wcs.cd1_1;
	wcs->cd[1] = fit->wcs.cd1_2;
	wcs->cd[2] = fit->wcs.cd2_1;
	wcs->cd[3] = fit->wcs.cd2_2;

	wcs->crval[0] = fit->wcs.crval1;
	wcs->crval[1] = fit->wcs.crval2;

	wcs->crota[0] = fit->wcs.crota1;
	wcs->crota[1] = fit->wcs.crota2;

	wcs->crpix[0] = fit->wcs.crpix1;
	wcs->crpix[1] = fit->wcs.crpix2;

	wcs->cdelt[0] = fit->wcs.cdelt1;
	wcs->cdelt[1] = fit->wcs.cdelt2;

	wcs->equinox = (double) fit->wcs.equinox;

	return TRUE;
#else
	return FALSE;
#endif
}

gboolean load_WCS(fits* fit) {
#ifdef HAVE_WCSLIB
	int status = 0;
	char *header;
	int nkeyrec, nreject, nwcs;

	if (wcs) {
		free_wcs();
	}

	if (fits_hdr2str(fit->fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
		report_fits_error(status);
		return FALSE;
	}

	if ((status = wcspih(header, nkeyrec, WCSHDR_all, -3, &nreject, &nwcs, &wcs))
			!= 0) {
		free(header);
		wcsvfree(&nwcs, &wcs);
		wcs = NULL;
		siril_debug_print("wcspih error %d: %s.\n", status, wcshdr_errmsg[status]);

		return FALSE;
	}

	free(header);

	if (wcs == NULL) {
	siril_debug_print("No world coordinate systems found.\n");
		return FALSE;
	}

	// FIXME: Call above goes through EVEN if no WCS is present, so we're adding this to return for now.
	if (wcs->crpix[0] == 0) {
		wcsvfree(&nwcs, &wcs);
		wcs = NULL;
		siril_debug_print("No world coordinate systems found.\n");
		return FALSE;
	}

	int stat[2], naxis[2];
	naxis[0] = fit->rx;
	naxis[1] = fit->ry;
	wcsfix(7, naxis, wcs, &stat[0]);
	if ((status = wcsset(wcs)) != 0) {
		wcsvfree(&nwcs, &wcs);
		wcs = NULL;
		siril_debug_print("wcsset error %d: %s.\n", status, wcs_errmsg[status]);
		return FALSE;
	}

	return TRUE;
#else
	return FALSE;
#endif
}

void pix2wcs(double pixel_x, double pixel_y, double *world_x, double *world_y) {
	*world_x = -1.0;
	*world_y = -1.0;
	if (wcs == NULL) return;
#ifdef HAVE_WCSLIB
	double phi = 0, theta = 0, world[2], pixcrd[2], imgcrd[2];
	int status, stat[2];

	pixcrd[0] = pixel_x;
	pixcrd[1] = pixel_y;

	if ((status = wcsp2s(wcs, 1, 2, &pixcrd[0], &imgcrd[0], &phi, &theta,
			&world[0], &stat[0])) != 0) {
	} else {
		*world_x = world[0];
		*world_y = world[1];
	}
#endif
}

void free_wcs() {
	// Clean up.
#ifdef HAVE_WCSLIB
	wcsfree(wcs);
	free(wcs);
	wcs = NULL;
#endif
}
