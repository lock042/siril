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

#include "core/siril.h"
#include "core/proto.h"

#include "io/image_format_fits.h"
#include "io/fits_keywords.h"

gchar *AstroTiff_build_header(fits *fit) {
	void *memptr;
	size_t memsize = IOBUFLEN;
	int status = 0;
	fitsfile *fptr = NULL;
	memptr = malloc(memsize);
	if (!memptr) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	fits_create_memfile(&fptr, &memptr, &memsize, IOBUFLEN, realloc, &status);
	if (status) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return NULL;
	}

	if (fits_create_img(fptr, fit->bitpix, fit->naxis, fit->naxes, &status)) {
		report_fits_error(status);
		if (fptr)
			fits_close_file(fptr, &status);
		free(memptr);
		return NULL;
	}

	fits tmpfit = { 0 };
	set_all_keywords_default(&tmpfit);
	copy_fits_metadata(fit, &tmpfit);
	tmpfit.fptr = fptr;
	save_fits_header(&tmpfit);

	if (fit->header)
		free(fit->header);
	fit->header = copy_header(&tmpfit);

	fits_close_file(tmpfit.fptr, &status);
	free(memptr);

	return g_strdup(fit->header);
}
