/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#include "../core/siril.h"
#include "../core/proto.h"
#include "../io/image_format_fits.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
	fits fits1 = {0}, fits2 = {0};

	if (argc < 3) {
		fprintf(stderr, "Usage: %s image1.fit image2.fit\n", *argv);
		exit(2);
	}

	if (readfits(argv[1], &fits1, NULL, FALSE) || readfits(argv[2], &fits2, NULL, FALSE)) {
		exit(2);
	}

	if (fits1.header && fits2.header) {
		if (strcmp(fits1.header, fits2.header))
			fprintf(stdout, "headers differ\n");
	}
	if (fits1.naxis != fits2.naxis) {
		fprintf(stdout, "number of axis differ\n");
		exit(1);
	}
	if (fits1.naxes[0] != fits2.naxes[0] || fits1.naxes[1] != fits2.naxes[1] ||
			fits1.naxes[2] != fits2.naxes[2]) {
		fprintf(stdout, "image axis differ\n");
		exit(1);
	}

	if (fits1.type != fits2.type) {
		fprintf(stdout, "image type differ\n");
		exit(1);
	}

	if (fits1.type == DATA_USHORT) {
		if (memcmp(fits1.data, fits2.data, fits1.naxes[0] * fits1.naxes[1] * fits1.naxes[2] * sizeof(WORD))) {
			fprintf(stdout, "image data differ\n");
			exit(1);
		}
	} else if (fits1.type == DATA_FLOAT) {
		if (memcmp(fits1.fdata, fits2.fdata, fits1.naxes[0] * fits1.naxes[1] * fits1.naxes[2] * sizeof(float))) {
			fprintf(stdout, "image data differ\n");
			exit(1);
		}
	}


	fprintf(stdout, "images are identical\n");
	return 0;
}
