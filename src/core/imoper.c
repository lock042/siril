/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "imoper.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "io/single_image.h"	// fit_get_max

/* set a minimum level for pixel values in an image */
int threshlo(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			*buf = max(level, *buf);
			buf++;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/* set a maximum level for pixel values in an image */
int threshhi(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			*buf = min(level, *buf);
			buf++;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/* replace all zero pixels in an image by the provided level */
int nozero(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			if (*buf == 0)
				*buf = level;
			buf++;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/* add an offset to all pixel values in an image */
int off(fits *fit, int level) {
	WORD *buf[3] =
			{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	int i, layer;
	assert(fit->naxes[2] <= 3);
	if (level == 0)
		return 0;
	if (level < -USHRT_MAX)
		level = -USHRT_MAX;
	else if (level > USHRT_MAX)
		level = USHRT_MAX;
	for (i = 0; i < fit->rx * fit->ry; ++i) {
		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			WORD val = buf[layer][i];
			if ((level < 0 && val < -level))
				buf[layer][i] = 0;
			else if (level > 0 && val > USHRT_MAX - level)
				buf[layer][i] = USHRT_MAX;
			else
				buf[layer][i] = val + level;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/*****************************************************************************
 *       S I R I L      A R I T H M E T I C      O P E R A T I O N S         *
 ****************************************************************************/

/* apply an operation with a scalar argument to all pixels of an image */
int soper(fits *a, double scalar, char oper) {
	WORD *gbuf;
	int i, layer;
	int n = a->rx * a->ry;

	assert(n > 0);

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		gbuf = a->pdata[layer];
		switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] + scalar);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] - scalar);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] * scalar);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] / scalar);
			}
			break;
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* applies operation of image a with image b, for all their layers:
 * a = a oper b
 * returns 0 on success */
int imoper(fits *a, fits *b, char oper) {
	int i, layer;

	if (a->rx != b->rx || a->ry != b->ry) {
		siril_log_message(
				_("imoper: images don't have the same size (w = %u|%u, h = %u|%u)\n"),
				a->rx, b->rx, a->ry, b->ry);
		return 1;
	}
	for (layer = 0; layer < a->naxes[2]; ++layer) {
		WORD *buf = b->pdata[layer];
		WORD *gbuf = a->pdata[layer];
		int n = a->rx * a->ry;
		switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] + buf[i]);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] - buf[i]);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] * buf[i]);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				gbuf[i] = (buf[i] == 0) ? 0 : round_to_WORD(gbuf[i] / buf[i]);
			}
			break;
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* This function applies a subtraction but contrary to Sub in imoper
 * it will use double type format.
 */
int sub_background(fits* image, fits* background, int layer) {
	double *pxl_image, *pxl_bkg;
	WORD *image_buf = image->pdata[layer];
	WORD *bkg_buf = background->pdata[layer];
	size_t i, ndata;
	double median;

	if ((image->rx) != (background->rx) || ((image->ry) != (background->ry))) {
		char *msg = siril_log_message(
				_("Images don't have the same size (w = %d|%d, h = %d|%d)\n"),
				image->rx, background->rx, image->ry, background->ry);
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
		return 1;
	}
	ndata = image->rx * image->ry;

	/* First step we convert data, apply the subtraction, normalize with median,
	 * and re-convert data to USHORT
	 */
	imstats *stat = statistics(NULL, -1, image, layer, NULL, STATS_BASIC);
	median = stat->median / USHRT_MAX_DOUBLE;
	pxl_image = malloc(sizeof(double) * ndata);
	pxl_bkg = malloc(sizeof(double) * ndata);

	for (i = 0; i < ndata; i++) {
		pxl_image[i] = (double) image_buf[i] / USHRT_MAX_DOUBLE;
		pxl_bkg[i] = (double) bkg_buf[i] / USHRT_MAX_DOUBLE;
		pxl_image[i] -= pxl_bkg[i];
		pxl_image[i] += median;
		image_buf[i] = round_to_WORD(pxl_image[i] * USHRT_MAX_DOUBLE);
	}

	invalidate_stats_from_fit(image);
	// We free memory
	free_stats(stat);
	free(pxl_image);
	free(pxl_bkg);
	return 0;
}

int addmax(fits *a, fits *b) {
	WORD *gbuf[3] = { a->pdata[RLAYER], a->pdata[GLAYER], a->pdata[BLAYER] };
	WORD *buf[3] = { b->pdata[RLAYER], b->pdata[GLAYER], b->pdata[BLAYER] };
	gint i, layer;

	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		siril_log_message(
				_("addmax: images don't have the same size (w = %d|%d, h = %d|%d, layers = %d|%d)\n"),
				a->rx, b->rx, a->ry, b->ry, a->naxes[2], b->naxes[2]);
		return 1;
	}
	assert(a->naxes[2] == 1 || a->naxes[2] == 3);

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		for (i = 0; i < a->ry * a->rx; ++i) {
			if (buf[layer][i] > gbuf[layer][i])
				gbuf[layer][i] = buf[layer][i];
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* If siril_fdiv is ok, function returns 0. If overflow, siril_fdiv returns 1*/
int siril_fdiv(fits *a, fits *b, float coef) {
	int i, layer;
	int retvalue = 0;
	double temp;

	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		fprintf(stderr, "Wrong size or channel count: %u=%u? / %u=%u?\n", a->rx,
				b->rx, a->ry, b->ry);
		return -1;
	}
	for (layer = 0; layer < a->naxes[2]; ++layer) {
		WORD *buf = b->pdata[layer];
		WORD *gbuf = a->pdata[layer];
		for (i = 0; i < b->rx * b->ry; ++i) {
			if (buf[i] == 0)
				buf[i] = 1;		// avoid division by 0
			temp = ((double) coef * ((double) gbuf[i] / (double) buf[i]));
			if (temp > USHRT_MAX_DOUBLE)
				retvalue = 1;
			gbuf[i] = round_to_WORD(temp);
		}
	}
	invalidate_stats_from_fit(a);
	return retvalue;
}

/* normalized division a/b, stored in a, with max value equal to the original
 * max value of a, for each layer. */
int siril_ndiv(fits *a, fits *b) {
	double *div;
	int layer, i, nb_pixels;
	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		fprintf(stderr,
				"Wrong size or channel count: %u=%u? / %u=%u?, %ld=%ld?\n",
				a->rx, b->rx, a->ry, b->ry, a->naxes[2], b->naxes[2]);
		return 1;
	}
	nb_pixels = a->rx * a->ry;
	div = malloc(nb_pixels * sizeof(double));

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		double max = 0, norm;
		for (i = 0; i < nb_pixels; ++i) {
			if (!b->pdata[layer][i])
				div[i] = (double) a->pdata[layer][i];
			else
				div[i] = (double) a->pdata[layer][i]
						/ (double) b->pdata[layer][i];
			max = max(div[i], max);
		}
		norm = max / fit_get_max(a, layer);
		for (i = 0; i < nb_pixels; ++i) {
			a->pdata[layer][i] = round_to_WORD(div[i] / norm);
		}
	}

	invalidate_stats_from_fit(a);
	free(div);
	return 0;
}

int loglut(fits *fit) {
	// This function maps fit with a log LUT
	int i, layer;
	WORD *buf[3] = { fit->pdata[RLAYER],
			fit->pdata[GLAYER], fit->pdata[BLAYER] };

	double norm = USHRT_MAX_DOUBLE / log(USHRT_MAX_DOUBLE);

	for (i = 0; i < fit->ry * fit->rx; i++) {
		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			double px = (double)buf[layer][i];
			px = (px == 0.0) ? 1.0 : px;
			buf[layer][i] = round_to_WORD(norm * log(px));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int asinhlut(fits *fit, double beta, double offset, gboolean RGBspace) {
	int i, layer;
	WORD *buf[3] = { fit->pdata[RLAYER],
			fit->pdata[GLAYER], fit->pdata[BLAYER] };
	double norm;

	norm = get_normalized_value(fit);

	for (i = 0; i < fit->ry * fit->rx; i++) {
		double x, k;
		if (fit->naxes[2] > 1) {
			double r, g, b;

			r = (double) buf[RLAYER][i] / norm;
			g = (double) buf[GLAYER][i] / norm;
			b = (double) buf[BLAYER][i] / norm;
			/* RGB space */
			if (RGBspace)
				x = 0.2126 * r + 0.7152 * g + 0.0722 * b;
			else
				x = 0.3333 * r + 0.3333 * g + 0.3333 * b;
		} else {
			x = buf[RLAYER][i] / norm;
		}

		k = asinh(beta * x) / (x * asinh(beta));

		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			double px = (double) buf[layer][i] / norm;
			px -= offset;
			px *= k;
			buf[layer][i] = round_to_WORD(px * norm);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int ddp(fits *a, int level, float coeff, float sigma) {
	fits fit;
	memset(&fit, 0, sizeof(fits));
	copyfits(a, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
	unsharp(&fit, sigma, 0, FALSE);
	soper(&fit, (double) level, OPER_ADD);
	nozero(&fit, 1);
	siril_fdiv(a, &fit, level);
	soper(a, (double) coeff, OPER_MUL);
	clearfits(&fit);
	invalidate_stats_from_fit(a);
	return 0;
}

