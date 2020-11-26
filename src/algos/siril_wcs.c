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


#include <stdio.h>
#include <wcslib.h>

#include "core/siril.h"

#include "siril_wcs.h"

static struct wcsprm *wcs;

gboolean load_WCS(fits* fit) {
#ifdef HAVE_WCSLIB
	int status = 0;
	char *header;
	int nkeyrec, nreject, nwcs;

	if (fits_hdr2str(fit->fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
		char errmsg[512];
		fits_get_errstatus(status, errmsg);
//        lastError = errmsg;
		return FALSE;
	}

	if ((status = wcspih(header, nkeyrec, WCSHDR_all, -3, &nreject, &nwcs, &wcs))
			!= 0) {
		free(header);
		wcsvfree(&nwcs, &wcs);
		wcs = NULL;
//        lastError = QString("wcspih ERROR %1: %2.").arg(status).arg(wcshdr_errmsg[status]);
		return FALSE;
	}

	free(header);

	if (wcs == NULL) {
//        lastError = i18n("No world coordinate systems found.");
		return FALSE;
	}

	// FIXME: Call above goes through EVEN if no WCS is present, so we're adding this to return for now.
	if (wcs->crpix[0] == 0) {
		wcsvfree(&nwcs, &wcs);
		wcs = NULL;
//        lastError = i18n("No world coordinate systems found.");
		return FALSE;
	}

	int stat1[2], naxis[2];
	naxis[0] = fit->rx;
	naxis[1] = fit->ry;
	wcsfix(7, naxis, wcs, &stat1[0]);
	if ((status = wcsset(wcs)) != 0) {
		wcsvfree(&nwcs, &wcs);
		wcs = NULL;
//        lastError = QString("wcsset error %1: %2.").arg(status).arg(wcs_errmsg[status]);
		return FALSE;
	}
//
	return TRUE;
#else
	return FALSE:
#endif
}

void pix2wcs(double pixel_x, double pixel_y, double *world_x, double *world_y) {
	*world_x = -1.0;
	*world_y = -1.0;
#ifdef HAVE_WCSLIB
	double phi = 0, theta = 0, world[2], pixcrd[2], imgcrd[2];
	int status, stat[2];

	pixcrd[0] = pixel_x;
	pixcrd[1] = pixel_y;

	if ((status = wcsp2s(wcs, 1, 2, &pixcrd[0], &imgcrd[0], &phi, &theta,
			&world[0], &stat[0])) != 0) {
	} else {
		//printf("(%lf / %lf)\n\n\n\n", world[0], world[1]);
		*world_x = world[0];
		*world_y = world[1];
	}
#endif
}

void free_wcs() {
	// Clean up.
	wcsfree(wcs);
	free(wcs);
	wcs = NULL;
}
