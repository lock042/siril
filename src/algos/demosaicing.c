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

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "io/fits_keywords.h"
#include "algos/demosaicing.h"
#include "algos/extraction.h"

#define USE_SIRIL_DEBAYER FALSE

const char *filter_pattern[] = {
	"RGGB",
	"BGGR",
	"GBRG",
	"GRBG",

/* XTRANS */
	"GGRGGB" // ----> XTRANS_1
	"GGBGGR"
	"BRGRBG"
	"GGBGGR"
	"GGRGGB"
	"RBGBRG",

	"RBGBRG" // ----> XTRANS_2
	"GGRGGB"
	"GGBGGR"
	"BRGRBG"
	"GGBGGR"
	"GGRGGB",

	"GRGGBG"
	"BGBRGR"
	"GRGGBG"
	"GBGGRG"
	"RGRBGB"
	"GBGGRG",

	"GBGGRG"
	"RGRBGB"
	"GBGGRG"
	"GRGGBG"
	"BGBRGR"
	"GRGGBG"
};
const size_t num_filter_patterns = G_N_ELEMENTS(filter_pattern);

/** Calculate the bayer pattern color from the row and column **/
static inline int FC(const size_t row, const size_t col, const uint32_t filters) {
	return filters >> (((row << 1 & 14) + (col & 1)) << 1) & 3;
}

/* width and height are sizes of the original image */
static void super_pixel_ushort(const WORD *buf, WORD *newbuf, int width, int height,
		sensor_pattern pattern) {
	long i = 0;
	for (int row = 0; row < height - 1; row += 2) {
		for (int col = 0; col < width - 1; col += 2) {
			float tmp;
			switch (pattern) {
			default:
			case BAYER_FILTER_RGGB:
				newbuf[i + 0] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[col + (1 + row) * width];
				newbuf[i + 1] = round_to_WORD(tmp * 0.5f);
				newbuf[i + 2] = buf[1 + col + (1 + row) * width];
				break;
			case BAYER_FILTER_BGGR:
				newbuf[i + 2] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[(col + row * width) + width];
				newbuf[i + 1] = round_to_WORD(tmp * 0.5f);
				newbuf[i + 0] = buf[(1 + col + row * width) + width];
				break;
			case BAYER_FILTER_GBRG:
				newbuf[i + 2] = buf[1 + col + row * width];
				newbuf[i + 0] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				newbuf[i + 1] = round_to_WORD(tmp * 0.5f);
				break;
			case BAYER_FILTER_GRBG:
				newbuf[i + 0] = buf[1 + col + row * width];
				newbuf[i + 2] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				newbuf[i + 1] = round_to_WORD(tmp * 0.5f);
				break;
			}
			i += 3;
		}
	}
}

/* width and height are sizes of the original image */
static void super_pixel_float(const float *buf, float *newbuf, int width, int height,
		sensor_pattern pattern) {
	long i = 0;
	for (int row = 0; row < height - 1; row += 2) {
		for (int col = 0; col < width - 1; col += 2) {
			float tmp;
			switch (pattern) {
			default:
			case BAYER_FILTER_RGGB:
				newbuf[i + 0] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[col + (1 + row) * width];
				newbuf[i + 1] = tmp * 0.5f;
				newbuf[i + 2] = buf[1 + col + (1 + row) * width];
				break;
			case BAYER_FILTER_BGGR:
				newbuf[i + 2] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[(col + row * width) + width];
				newbuf[i + 1] = tmp * 0.5f;
				newbuf[i + 0] = buf[(1 + col + row * width) + width];
				break;
			case BAYER_FILTER_GBRG:
				newbuf[i + 2] = buf[1 + col + row * width];
				newbuf[i + 0] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				newbuf[i + 1] = tmp * 0.5f;
				break;
			case BAYER_FILTER_GRBG:
				newbuf[i + 0] = buf[1 + col + row * width];
				newbuf[i + 2] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				newbuf[i + 1] = tmp * 0.5f;
				break;
			}
			i += 3;
		}
	}
}

/***************************************************
 *
 * Written by Damien Douxchamps and Frederic Devernay
 * The original VNG and AHD Bayer decoding are from Dave Coffin's DCRAW.
 * https://code.google.com/p/gst-plugins-elphel/
 *
 * *************************************************/

static void ClearBorders(WORD *rgb, int sx, int sy, int w) {
	int i, j;

	/* black edges: */
	i = 3 * sx * w - 1;
	j = 3 * sx * sy - 1;
	while (i >= 0) {
		rgb[i--] = 0;
		rgb[j--] = 0;
	}

	int low = sx * (w - 1) * 3 - 1 + w * 3;
	i = low + sx * (sy - w * 2 + 1) * 3;
	while (i > low) {
		j = 6 * w;
		while (j > 0) {
			rgb[i--] = 0;
			j--;
		}
		i -= (sx - 2 * w) * 3;
	}
}

/* OpenCV's Bayer decoding */
static int bayer_Bilinear(const WORD *bayer, WORD *rgb, int sx, int sy,
		sensor_pattern tile) {
	const int bayerStep = sx;
	const int rgbStep = 3 * sx;
	int width = sx;
	int height = sy;
	int blue = tile == BAYER_FILTER_BGGR || tile == BAYER_FILTER_GBRG ? -1 : 1;
	int start_with_green = tile == BAYER_FILTER_GBRG
			|| tile == BAYER_FILTER_GRBG;

	if (tile > BAYER_FILTER_MAX || tile < BAYER_FILTER_MIN)
		return -1;

	ClearBorders(rgb, sx, sy, 1);
	rgb += rgbStep + 3 + 1;
	height -= 2;
	width -= 2;

	for (; height--; bayer += bayerStep, rgb += rgbStep) {
		int t0, t1;
		const WORD *bayerEnd = bayer + width;

		if (start_with_green) {
			t0 = (bayer[1] + bayer[bayerStep * 2 + 1] + 1) >> 1;
			t1 = (bayer[bayerStep] + bayer[bayerStep + 2] + 1) >> 1;
			rgb[-blue] = round_to_WORD(t0);
			rgb[0] = bayer[bayerStep + 1];
			rgb[blue] = round_to_WORD(t1);
			bayer++;
			rgb += 3;
		}

		if (blue > 0) {
			for (; bayer <= bayerEnd - 2; bayer += 2, rgb += 6) {
				t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
						+ bayer[bayerStep * 2 + 2] + 2) >> 2;
				t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
						+ bayer[bayerStep * 2 + 1] + 2) >> 2;
				rgb[-1] = round_to_WORD(t0);
				rgb[0] = round_to_WORD(t1);
				rgb[1] = bayer[bayerStep + 1];

				t0 = (bayer[2] + bayer[bayerStep * 2 + 2] + 1) >> 1;
				t1 = (bayer[bayerStep + 1] + bayer[bayerStep + 3] + 1) >> 1;
				rgb[2] = round_to_WORD(t0);
				rgb[3] = bayer[bayerStep + 2];
				rgb[4] = round_to_WORD(t1);
			}
		} else {
			for (; bayer <= bayerEnd - 2; bayer += 2, rgb += 6) {
				t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
						+ bayer[bayerStep * 2 + 2] + 2) >> 2;
				t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
						+ bayer[bayerStep * 2 + 1] + 2) >> 2;
				rgb[1] = round_to_WORD(t0);
				rgb[0] = round_to_WORD(t1);
				rgb[-1] = bayer[bayerStep + 1];

				t0 = (bayer[2] + bayer[bayerStep * 2 + 2] + 1) >> 1;
				t1 = (bayer[bayerStep + 1] + bayer[bayerStep + 3] + 1) >> 1;
				rgb[4] = round_to_WORD(t0);
				rgb[3] = bayer[bayerStep + 2];
				rgb[2] = round_to_WORD(t1);
			}
		}

		if (bayer < bayerEnd) {
			t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
					+ bayer[bayerStep * 2 + 2] + 2) >> 2;
			t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
					+ bayer[bayerStep * 2 + 1] + 2) >> 2;
			rgb[-blue] = round_to_WORD(t0);
			rgb[0] = round_to_WORD(t1);
			rgb[blue] = bayer[bayerStep + 1];
			bayer++;
			rgb += 3;
		}

		bayer -= width;
		rgb -= width * 3;

		blue = -blue;
		start_with_green = !start_with_green;
	}

	return 0;
}

#define ABSOLU(x) (((int)(x) ^ ((int)(x) >> 31)) - ((int)(x) >> 31))

static int bayer_VNG(const WORD *bayer, WORD *dst, int sx, int sy,
		sensor_pattern pattern) {
	const signed char bayervng_terms[] = { -2, -2, +0, -1, 0, 0x01, -2, -2,
			+0, +0, 1, 0x01, -2, -1, -1, +0, 0, 0x01, -2, -1, +0, -1, 0, 0x02, -2,
			-1, +0, +0, 0, 0x03, -2, -1, +0, +1, 1, 0x01, -2, +0, +0, -1, 0, 0x06,
			-2, +0, +0, +0, 1, 0x02, -2, +0, +0, +1, 0, 0x03, -2, +1, -1, +0, 0,
			0x04, -2, +1, +0, -1, 1, 0x04, -2, +1, +0, +0, 0, 0x06, -2, +1, +0, +1,
			0, 0x02, -2, +2, +0, +0, 1, 0x04, -2, +2, +0, +1, 0, 0x04, -1, -2, -1,
			+0, 0, 0x80, -1, -2, +0, -1, 0, 0x01, -1, -2, +1, -1, 0, 0x01, -1, -2,
			+1, +0, 1, 0x01, -1, -1, -1, +1, 0, 0x88, -1, -1, +1, -2, 0, 0x40, -1,
			-1, +1, -1, 0, 0x22, -1, -1, +1, +0, 0, 0x33, -1, -1, +1, +1, 1, 0x11,
			-1, +0, -1, +2, 0, 0x08, -1, +0, +0, -1, 0, 0x44, -1, +0, +0, +1, 0,
			0x11, -1, +0, +1, -2, 1, 0x40, -1, +0, +1, -1, 0, 0x66, -1, +0, +1, +0,
			1, 0x22, -1, +0, +1, +1, 0, 0x33, -1, +0, +1, +2, 1, 0x10, -1, +1, +1,
			-1, 1, 0x44, -1, +1, +1, +0, 0, 0x66, -1, +1, +1, +1, 0, 0x22, -1, +1,
			+1, +2, 0, 0x10, -1, +2, +0, +1, 0, 0x04, -1, +2, +1, +0, 1, 0x04, -1,
			+2, +1, +1, 0, 0x04, +0, -2, +0, +0, 1, 0x80, +0, -1, +0, +1, 1, 0x88,
			+0, -1, +1, -2, 0, 0x40, +0, -1, +1, +0, 0, 0x11, +0, -1, +2, -2, 0,
			0x40, +0, -1, +2, -1, 0, 0x20, +0, -1, +2, +0, 0, 0x30, +0, -1, +2, +1,
			1, 0x10, +0, +0, +0, +2, 1, 0x08, +0, +0, +2, -2, 1, 0x40, +0, +0, +2,
			-1, 0, 0x60, +0, +0, +2, +0, 1, 0x20, +0, +0, +2, +1, 0, 0x30, +0, +0,
			+2, +2, 1, 0x10, +0, +1, +1, +0, 0, 0x44, +0, +1, +1, +2, 0, 0x10, +0,
			+1, +2, -1, 1, 0x40, +0, +1, +2, +0, 0, 0x60, +0, +1, +2, +1, 0, 0x20,
			+0, +1, +2, +2, 0, 0x10, +1, -2, +1, +0, 0, 0x80, +1, -1, +1, +1, 0,
			0x88, +1, +0, +1, +2, 0, 0x08, +1, +0, +2, -1, 0, 0x40, +1, +0, +2, +1,
			0, 0x10 }, bayervng_chood[] = { -1, -1, -1, 0, -1, +1, 0, +1, +1, +1,
			+1, 0, +1, -1, 0, -1 };
	const int height = sy, width = sx;
	const signed char *cp;
	WORD (*brow[5])[3], *pix; /* [FD] */
	int code[8][2][320], *ip, gval[8], gmin, gmax, sum[4];
	int row, col, x, y, x1, x2, y1, y2, t, weight, grads, color, diag;
	int g, diff, thold, num, c;
	unsigned long int filters; /* [FD] */

	/* first, use bilinear bayer decoding */

	bayer_Bilinear(bayer, dst, sx, sy, pattern);

	switch (pattern) {
	case BAYER_FILTER_BGGR:
		filters = 0x16161616;
		break;
	case BAYER_FILTER_GRBG:
		filters = 0x61616161;
		break;
	case BAYER_FILTER_RGGB:
		filters = 0x94949494;
		break;
	case BAYER_FILTER_GBRG:
		filters = 0x49494949;
		break;
	default:
		return -1;
	}

	for (row = 0; row < 8; row++) { /* Precalculate for VNG */
		for (col = 0; col < 2; col++) {
			ip = code[row][col];
			for (cp = bayervng_terms, t = 0; t < 64; t++) {
				y1 = *cp++;
				x1 = *cp++;
				y2 = *cp++;
				x2 = *cp++;
				weight = *cp++;
				grads = *cp++;
				color = FC(row + y1, col + x1, filters);
				if (FC(row+y2,col+x2, filters) != (unsigned int) color)
					continue;
				diag = (FC(row,col+1, filters) == (unsigned int) color
						&& FC(row+1,col, filters) == (unsigned int) color) ? 2 : 1;
				if (abs(y1 - y2) == diag && abs(x1 - x2) == diag)
					continue;
				*ip++ = (y1 * width + x1) * 3 + color; /* [FD] */
				*ip++ = (y2 * width + x2) * 3 + color; /* [FD] */
				*ip++ = weight;
				for (g = 0; g < 8; g++)
					if (grads & 1 << g)
						*ip++ = g;
				*ip++ = -1;
			}
			*ip++ = INT_MAX;
			for (cp = bayervng_chood, g = 0; g < 8; g++) {
				y = *cp++;
				x = *cp++;
				*ip++ = (y * width + x) * 3; /* [FD] */
				color = FC(row, col, filters);
				if (FC(row+y,col+x, filters) != (unsigned int) color
						&& FC(row+y*2,col+x*2, filters) == (unsigned int) color)
					*ip++ = (y * width + x) * 6 + color; /* [FD] */
				else
					*ip++ = 0;
			}
		}
	}
	brow[4] = calloc(width * 3, sizeof **brow);
	//merror (brow[4], "vng_interpolate()");
	for (row = 0; row < 3; row++)
		brow[row] = brow[4] + row * width;
	for (row = 2; row < height - 2; row++) { /* Do VNG interpolation */
		for (col = 2; col < width - 2; col++) {
			pix = dst + (row * width + col) * 3; /* [FD] */
			ip = code[row & 7][col & 1];
			memset(gval, 0, sizeof gval);
			while ((g = ip[0]) != INT_MAX) { /* Calculate gradients */
				diff = ABSOLU(pix[g] - pix[ip[1]]) << ip[2];
				gval[ip[3]] += diff;
				ip += 5;
				if ((g = ip[-1]) == -1)
					continue;
				gval[g] += diff;
				while ((g = *ip++) != -1)
					gval[g] += diff;
			}
			ip++;
			gmin = gmax = gval[0]; /* Choose a threshold */
			for (g = 1; g < 8; g++) {
				if (gmin > gval[g])
					gmin = gval[g];
				if (gmax < gval[g])
					gmax = gval[g];
			}
			if (gmax == 0) {
				memcpy(brow[2][col], pix, 3 * sizeof *dst); /* [FD] */
				continue;
			}
			thold = gmin + (gmax >> 1);
			memset(sum, 0, sizeof sum);
			color = FC(row, col, filters);
			for (num = g = 0; g < 8; g++, ip += 2) { /* Average the neighbors */
				if (gval[g] <= thold) {
					for (c = 0; c < 3; c++) /* [FD] */
						if (c == color && ip[1])
							sum[c] += (pix[c] + pix[ip[1]]) >> 1;
						else
							sum[c] += pix[ip[0] + c];
					num++;
				}
			}
			for (c = 0; c < 3; c++) { /* [FD] Save to buffer */
				t = pix[color];
				if (c != color)
					t += (sum[c] - sum[color]) / num;
				//~ CLIP16(t,brow[2][col][c], 16); /* [FD] */
				brow[2][col][c] = round_to_WORD(t); /* [FD] */
			}
		}
		if (row > 3) /* Write buffer to image */
			memcpy(dst + 3 * ((row - 2) * width + 2), brow[0] + 2,
					(width - 4) * 3 * sizeof *dst); /* [FD] */
		for (g = 0; g < 4; g++)
			brow[(g - 1) & 3] = brow[g];
	}
	memcpy(dst + 3 * ((row - 2) * width + 2), brow[0] + 2,
			(width - 4) * 3 * sizeof *dst);
	memcpy(dst + 3 * ((row - 1) * width + 2), brow[1] + 2,
			(width - 4) * 3 * sizeof *dst);
	free(brow[4]);

	return 0;
}

/* AHD interpolation ported from dcraw to libdc1394 by Samuel Audet */
static gboolean ahd_inited = FALSE; /* WARNING: not multi-processor safe */

#define LIM(x,min,max) MAX(min,MIN(x,max))
#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))

static const float xyz_rgb[3][3] = { /* XYZ from RGB */
	{ 0.412453f, 0.357580f, 0.180423f },
	{ 0.212671f, 0.715160f, 0.072169f },
	{ 0.019334f, 0.119193f, 0.950227f } };
/* TODO: is it wise to use a D65 here? */
static const float d65_white[3] = { 0.950456f, 1.0f, 1.088754f };
/* TODO: store the precomputation of xyz_rgb * d65 instead of running a silly init */

static void cam_to_cielab(const uint16_t cam[3], float lab[3]) /* [SA] */
{
	float xyz[3];
	static float cbrt[0x10000], xyz_cam[3][4];

	if (cam == NULL) {
		int i, j;

		for (i = 0; i < 0x10000; i++) {
			float r = i / 65535.0;
			cbrt[i] = r > 0.008856 ? pow(r, 1 / 3.0) : 7.787 * r + 16 / 116.0;
		}
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++) /* [SA] */
				xyz_cam[i][j] = xyz_rgb[i][j] / d65_white[i]; /* [SA] */
	} else {
		int c;

		xyz[0] = xyz[1] = xyz[2] = 0.5;
		for (c = 0; c < 3; c++)
		{ /* [SA] */
			xyz[0] += xyz_cam[0][c] * cam[c];
			xyz[1] += xyz_cam[1][c] * cam[c];
			xyz[2] += xyz_cam[2][c] * cam[c];
		}
		xyz[0] = cbrt[round_to_WORD(xyz[0])]; /* [SA] */
		xyz[1] = cbrt[round_to_WORD(xyz[1])]; /* [SA] */
		xyz[2] = cbrt[round_to_WORD(xyz[2])]; /* [SA] */
		lab[0] = 116 * xyz[1] - 16;
		lab[1] = 500 * (xyz[0] - xyz[1]);
		lab[2] = 200 * (xyz[1] - xyz[2]);
	}
}

/*
 Adaptive Homogeneity-Directed interpolation is based on
 the work of Keigo Hirakawa, Thomas Parks, and Paul Lee.
 */
#define TS 256 /* Tile Size */

static int bayer_AHD(const WORD *bayer, WORD *dst, int sx, int sy,
		sensor_pattern pattern) {
	int i, j, top, left, row, col, tr, tc, fc, c, d, val, hm[2];
	/* the following has the same type as the image */
	uint16_t (*pix)[3], (*rix)[3]; /* [SA] */
	static const int dir[4] = { -1, 1, -TS, TS };
	unsigned ldiff[2][4], abdiff[2][4], leps, abeps;
	float flab[3];
	uint16_t (*rgb)[TS][TS][3]; /* [SA] */
	short (*lab)[TS][TS][3];
	char (*homo)[TS][TS], *buffer;
	/* start - new code for libdc1394 */
	uint32_t filters;
	const int height = sy, width = sx;
	int x, y;

	if (ahd_inited == FALSE) {
		/* WARNING: this might not be multi-processor safe */
		cam_to_cielab(NULL, NULL);
		ahd_inited = TRUE;
	}

	switch (pattern) {
	case BAYER_FILTER_BGGR:
		filters = 0x16161616;
		break;
	case BAYER_FILTER_GRBG:
		filters = 0x61616161;
		break;
	case BAYER_FILTER_RGGB:
		filters = 0x94949494;
		break;
	case BAYER_FILTER_GBRG:
		filters = 0x49494949;
		break;
	default:
		return -1;
	}

	/* fill-in destination with known exact values */
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			int channel = FC(y, x, filters);
			dst[(y * width + x) * 3 + channel] = bayer[y * width + x];
		}
	}
	/* end - new code for libdc1394 */

	/* start - code from border_interpolate(int border) */
	{
		int border = 3;
		unsigned row, col, y, x, f, c, sum[8];

		for (row = 0; row < (unsigned int) height; row++)
			for (col = 0; col < (unsigned int) width; col++) {
				if (col == (unsigned int) border && row >= (unsigned int) border
						&& row < (unsigned int) height - (unsigned int) border)
					col = width - border;
				memset(sum, 0, sizeof sum);
				for (y = row - 1; y != row + 2; y++)
					for (x = col - 1; x != col + 2; x++)
						if (y < (unsigned int) height
								&& x < (unsigned int) width) {
							f = FC(y, x, filters);
							sum[f] += dst[(y * width + x) * 3 + f]; /* [SA] */
							sum[f + 4]++;
						}
				f = FC(row, col, filters);
				for (c = 0; c < 3; c++)
					if (c != f && sum[c + 4]) /* [SA] */
						dst[(row * width + col) * 3 + c] = sum[c] / sum[c + 4]; /* [SA] */
			}
	}
	/* end - code from border_interpolate(int border) */

	buffer = (char *) malloc(26 * TS * TS); /* 1664 kB */
	/* merror (buffer, "ahd_interpolate()"); */
	rgb = (uint16_t (*)[TS][TS][3]) buffer; /* [SA] */
	lab = (short (*)[TS][TS][3]) (buffer + 12 * TS * TS);
	homo = (char (*)[TS][TS]) (buffer + 24 * TS * TS);

	for (top = 0; top < height; top += TS - 6)
		for (left = 0; left < width; left += TS - 6) {
			memset(rgb, 0, 12 * TS * TS);

			/* Interpolate green horizontally and vertically: */
			for (row = ((top < 2) ? 2 : top);
					row < top + TS && row < height - 2; row++) {
				col = left + (FC(row,left, filters) == 1);
				if (col < 2)
					col += 2;
				for (fc = FC(row, col, filters); col < left + TS && col < width - 2;
						col += 2) {
					pix = (uint16_t (*)[3]) dst + (row * width + col); /* [SA] */
					val = ((pix[-1][1] + pix[0][fc] + pix[1][1]) * 2
							- pix[-2][fc] - pix[2][fc]) >> 2;
					rgb[0][row - top][col - left][1] = ULIM(val, pix[-1][1],
							pix[1][1]);
					val = ((pix[-width][1] + pix[0][fc] + pix[width][1]) * 2
							- pix[-2 * width][fc] - pix[2 * width][fc]) >> 2;
					rgb[1][row - top][col - left][1] = ULIM(val, pix[-width][1],
							pix[width][1]);
				}
			}
			/* Interpolate red and blue, and convert to CIELab: */
			for (d = 0; d < 2; d++)
				for (row = top + 1; row < top + TS - 1 && row < height - 1;
						row++)
					for (col = left + 1; col < left + TS - 1 && col < width - 1;
							col++) {
						pix = (uint16_t (*)[3]) dst + (row * width + col); /* [SA] */
						rix = &rgb[d][row - top][col - left];
						if ((c = 2 - FC(row, col, filters)) == 1) {
							c = FC(row + 1, col, filters);
							val = pix[0][1]
									+ ((pix[-1][2 - c] + pix[1][2 - c]
											- rix[-1][1] - rix[1][1]) >> 1);
							rix[0][2 - c] = round_to_WORD(val); /* [SA] */
							val = pix[0][1]
									+ ((pix[-width][c] + pix[width][c]
											- rix[-TS][1] - rix[TS][1]) >> 1);
						} else
							val = rix[0][1]
									+ ((pix[-width - 1][c] + pix[-width + 1][c]
											+ pix[+width - 1][c]
											+ pix[+width + 1][c]
											- rix[-TS - 1][1] - rix[-TS + 1][1]
											- rix[+TS - 1][1] - rix[+TS + 1][1]
											+ 1) >> 2);
						rix[0][c] = round_to_WORD((double) val); /* [SA] */
						c = FC(row, col, filters);
						rix[0][c] = pix[0][c];
						cam_to_cielab(rix[0], flab);
						for (c = 0; c < 3; c++)
							lab[d][row - top][col - left][c] = 64 * flab[c];
					}
			/* Build homogeneity maps from the CIELab images: */
			memset(homo, 0, 2 * TS * TS);
			for (row = top + 2; row < top + TS - 2 && row < height; row++) {
				tr = row - top;
				for (col = left + 2; col < left + TS - 2 && col < width;
						col++) {
					tc = col - left;
					for (d = 0; d < 2; d++)
						for (i = 0; i < 4; i++)
							ldiff[d][i] = ABSOLU(
									lab[d][tr][tc][0]
											- lab[d][tr][tc + dir[i]][0]);
					leps = MIN(MAX(ldiff[0][0],ldiff[0][1]),
							MAX(ldiff[1][2],ldiff[1][3]));
					for (d = 0; d < 2; d++)
						for (i = 0; i < 4; i++)
							if (i >> 1 == d || ldiff[d][i] <= leps)
								abdiff[d][i] =
										SQR(
												lab[d][tr][tc][1]
														- lab[d][tr][tc + dir[i]][1]) + SQR(lab[d][tr][tc][2] -lab[d][tr][tc+dir[i]][2]);
					abeps = MIN(MAX(abdiff[0][0],abdiff[0][1]),
							MAX(abdiff[1][2],abdiff[1][3]));
					for (d = 0; d < 2; d++)
						for (i = 0; i < 4; i++)
							if (ldiff[d][i] <= leps && abdiff[d][i] <= abeps)
								homo[d][tr][tc]++;
				}
			}
			/* Combine the most homogenous pixels for the final result: */
			for (row = top + 3; row < top + TS - 3 && row < height - 3; row++) {
				tr = row - top;
				for (col = left + 3; col < left + TS - 3 && col < width - 3;
						col++) {
					tc = col - left;
					for (d = 0; d < 2; d++)
						for (hm[d] = 0, i = tr - 1; i <= tr + 1; i++)
							for (j = tc - 1; j <= tc + 1; j++)
								hm[d] += homo[d][i][j];
					if (hm[0] != hm[1])
						for (c = 0; c < 3; c++)
							dst[(row * width + col) * 3 + c] = round_to_WORD(
									rgb[hm[1] > hm[0]][tr][tc][c]); /* [SA] */
					else
						for (c = 0; c < 3; c++)
							dst[(row * width + col) * 3 + c] = round_to_WORD(
									(rgb[0][tr][tc][c] + rgb[1][tr][tc][c])
											>> 1); /* [SA] */
				}
			}
		}
	free(buffer);

	return 0;
}

#define fcol(row, col) xtrans[(row) % 6][(col) % 6]

/* Code from RAWTherapee:
 * It is a simple algorithm. Certainly not the best (probably the worst) but it works yet. */
static int fast_xtrans_interpolate(const WORD *bayer, WORD *dst, int sx, int sy,
		unsigned int xtrans[6][6]) {
	uint32_t filters = 9;
	const int height = sy, width = sx;
	int row;

	/* start - code from border_interpolate(int border) */
	{
		int border = 1;
		unsigned row, col, y, x, f, c, sum[8];

		for (row = 0; row < (unsigned int) height; row++)
			for (col = 0; col < (unsigned int) width; col++) {
				if (col == (unsigned int) border && row >= (unsigned int) border
						&& row < (unsigned int) height - (unsigned int) border)
					col = width - border;
				memset(sum, 0, sizeof sum);
				for (y = row - 1; y != row + 2; y++)
					for (x = col - 1; x != col + 2; x++)
						if (y < (unsigned int) height
								&& x < (unsigned int) width) {
							f = FC(y, x, filters);
							sum[f] += dst[(y * width + x) * 3 + f]; /* [SA] */
							sum[f + 4]++;
						}
				f = FC(row, col, filters);
				for (c = 0; c < 3; c++)
					if (c != f && sum[c + 4]) /* [SA] */
						dst[(row * width + col) * 3 + c] = sum[c] / sum[c + 4]; /* [SA] */
			}
	}
	/* end - code from border_interpolate(int border) */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (row = 1; row < (height - 1); row ++) {
		int col;
		for (col = 1; col < (width - 1); col ++) {
			float sum[3] = { 0.f };

			int v, h;
			for (v = -1; v <= 1; v++) {
				for (h = -1; h <= 1; h++) {
					sum[fcol(row + v, col + h)] += bayer[(col + h) + (row + v) * width];
				}
			}

			switch (fcol(row, col)) {
			case 0: /* Red */
				dst[(row * width + col) * 3 + 0] = bayer[col + row * width];
				dst[(row * width + col) * 3 + 1] = sum[1] * 0.2f;
				dst[(row * width + col) * 3 + 2] = sum[2] * 0.33333333f;
				break;

			case 1: /* Green */
				dst[(row * width + col) * 3 + 0] = sum[0] * 0.5f;
				dst[(row * width + col) * 3 + 1] = bayer[col + row * width];
				dst[(row * width + col) * 3 + 2] = sum[2] * 0.5f;
				break;

			case 2: /* Blue */
				dst[(row * width + col) * 3 + 0] = sum[0] * 0.33333333f;
				dst[(row * width + col) * 3 + 1] = sum[1] * 0.2f;
				dst[(row * width + col) * 3 + 2] = bayer[col + row * width];
				break;
			}

		}
	}
	return 0;
}

#undef fcol

static WORD *debayer_buffer_siril(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6], int bit_depth) {
	WORD *newbuf;
	size_t npixels;
	int retval;
	switch (interpolation) {
		default:
			npixels = (*width) * (*height);
			break;
		case BAYER_SUPER_PIXEL:
			npixels = (*width / 2 + *width % 2) * (*height / 2 + *height % 2);
			break;
	}
	newbuf = calloc(3, npixels * sizeof(WORD));
	if (!newbuf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	switch (interpolation) {
	case BAYER_BILINEAR:
		retval = bayer_Bilinear(buf, newbuf, *width, *height, pattern);
		break;
//	case BAYER_NEARESTNEIGHBOR:
//		retval = bayer_NearestNeighbor(buf, newbuf, *width, *height, pattern);
//		break;
	default:
	case BAYER_VNG:
		retval = bayer_VNG(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_AHD:
		retval = bayer_AHD(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_SUPER_PIXEL:
		super_pixel_ushort(buf, newbuf, *width, *height, pattern);
		*width = *width / 2 + *width % 2;
		*height = *height / 2 + *height % 2;
		retval = 0;
		break;
	case XTRANS:
		if (!xtrans) {
			free(newbuf);
			return NULL;
		}
		retval = fast_xtrans_interpolate(buf, newbuf, *width, *height, xtrans);
		break;
	}
	if (retval) {
		free(newbuf);
		return NULL;
	}
	return newbuf;
}

int adjust_Bayer_pattern(fits *fit, sensor_pattern *pattern) {
	int xbayeroff = 0, ybayeroff = 0;
	gboolean top_down = com.pref.debayer.top_down;

	if (com.pref.debayer.use_bayer_header) {
		if (!g_strcmp0(fit->keywords.row_order, "TOP-DOWN")) {
			top_down = TRUE;
		} else if (!g_strcmp0(fit->keywords.row_order, "BOTTOM-UP")) {
			top_down = FALSE;
		}
	}

	siril_debug_print("Debayer will be done %s\n", top_down ? "top-down" : "bottom-up");

	if (!com.pref.debayer.use_bayer_header) {
		xbayeroff = com.pref.debayer.xbayeroff;
		ybayeroff = com.pref.debayer.ybayeroff;
	} else {
		xbayeroff = fit->keywords.bayer_xoffset;
		ybayeroff = fit->keywords.bayer_yoffset;
	}

	if (xbayeroff == 1) {
		switch (*pattern) {
		case BAYER_FILTER_RGGB:
			*pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_BGGR:
			*pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_GBRG:
			*pattern = BAYER_FILTER_BGGR;
			break;
		case BAYER_FILTER_GRBG:
			*pattern = BAYER_FILTER_RGGB;
			break;
		default:
			return 1;
		}
	}

	if (ybayeroff == 1) {
		switch (*pattern) {
		case BAYER_FILTER_RGGB:
			*pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_BGGR:
			*pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_GBRG:
			*pattern = BAYER_FILTER_RGGB;
			break;
		case BAYER_FILTER_GRBG:
			*pattern = BAYER_FILTER_BGGR;
			break;
		default:
			return 1;
		}
	}

	/* read bottom-up */
	if (!top_down && !(fit->ry % 2)) {
		switch (*pattern) {
		case BAYER_FILTER_RGGB:
			*pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_BGGR:
			*pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_GBRG:
			*pattern = BAYER_FILTER_RGGB;
			break;
		case BAYER_FILTER_GRBG:
			*pattern = BAYER_FILTER_BGGR;
			break;
		default:
			return 1;
		}
	}
	return 0;
}

static void adjust_XTrans_pattern(fits *fit, unsigned int xtrans[6][6]) {
	gboolean top_down = com.pref.debayer.top_down;

	if (com.pref.debayer.use_bayer_header) {
		if (!g_strcmp0(fit->keywords.row_order, "TOP-DOWN")) {
			top_down = TRUE;
		} else if (!g_strcmp0(fit->keywords.row_order, "BOTTOM-UP")) {
			top_down = FALSE;
		}
	}
	siril_debug_print("X-Trans pattern will be used %s\n", top_down ? "top-down" : "bottom-up (inverted)");

	if (!top_down) {
		int offset = fit->ry % 6;
		if (offset)
			siril_debug_print("Image with an X-Trans sensor doesn't have a height multiple of 6\n");
		unsigned int orig[6][6];
		memcpy(orig, xtrans, 36 * sizeof(unsigned int));
		for (int i = 0; i < 6; i++) {
			int y = (5 - i + offset) % 6;
			for (int j = 0; j < 6; j++) {
				xtrans[i][j] = orig[y][j];
			}
		}
	}
}

/* convert the string-described X-Trans pattern into an int array with value corresponding to filter */
static int compile_XTrans_pattern(const char *bayer, unsigned int xtrans[6][6]) {
	int i = 0;

	if (strlen(bayer) != 36) {
		siril_log_color_message(_("FITS header does not contain a proper XTRANS pattern, demosaicing cannot be done\n"), "red");
		return 1;
	}

	for (int y = 0; y < 6; y++) {
		for (int x = 0; x < 6; x++) {
			switch (bayer[i]) {
			case 'R':
				xtrans[y][x] = 0;
				break;
			case 'G':
				xtrans[y][x] = 1;
				break;
			case 'B':
				xtrans[y][x] = 2;
				break;
			default:
				siril_log_color_message(_("Invalid character in X-Trans filter pattern: %c\n"), "red", bayer[i]);
				return 1;
			}
			i++;
		}
	}
	return 0;
}

WORD *debayer_buffer_superpixel_ushort(WORD *buf, int *width, int *height, sensor_pattern pattern) {
	int new_rx = *width / 2 + *width % 2;
	int new_ry = *height / 2 + *height % 2;
	size_t npixels = new_rx * new_ry;
	WORD *newbuf = malloc(3 * npixels * sizeof(WORD));
	if (!newbuf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	super_pixel_ushort(buf, newbuf, *width, *height, pattern);
	*width = new_rx;
	*height = new_ry;
	return newbuf;
}

/**
 * debayer a buffer of a given size into a newly allocated and returned buffer,
 * using the given bayer pattern and interpolation (only used for SER demosaicing)
 *
 * @param buf original RAW data
 * @param width width of image
 * @param height height of image
 * @param interpolation type of interpolation used for demosaicing algorithm
 * @param pattern type of pattern used for demosaicing algorithm
 * @return a new buffer of demosaiced data
 */
WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, int bit_depth) {
	if (USE_SIRIL_DEBAYER)
		return debayer_buffer_siril(buf, width, height, interpolation, pattern, NULL, bit_depth);
	return debayer_buffer_new_ushort(buf, width, height, interpolation, pattern, NULL, bit_depth);
}

float *debayer_buffer_superpixel_float(float *buf, int *width, int *height, sensor_pattern pattern) {
	int new_rx = *width / 2 + *width % 2;
	int new_ry = *height / 2 + *height % 2;
	size_t npixels = new_rx * new_ry;
	float *newbuf = malloc(3 * npixels * sizeof(float));
	if (!newbuf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	super_pixel_float(buf, newbuf, *width, *height, pattern);
	*width = new_rx;
	*height = new_ry;
	return newbuf;
}

/* From an area, get the area corresponding to the debayer data for all colors,
 * the dashed area below.
 * 0 1 2 3 4 5
 * - - - - - -
 * - - - - - -
 * - - G R - -
 * - - B G - -
 * - - - - - -
 * - - - - - -
 *
 * area is the requested area of an image (simplified as GRBG above)
 * debayer_area is the result of this function, the area with enough pixels to
 *	have a valid debayer
 * image_area is the size of the image, to avoid going out of bounds
 * debayer_offset_x and y are the offset that need to be applied to the debayer
 *	data to find the original area (between 0 and 3).
 */
void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y) {
	int right, bottom;	// temp debayer negative offsets

	/* left side */
	if (area->x & 1)
		*debayer_offset_x = 11;
	else
		*debayer_offset_x = 10;
	if (area->x - *debayer_offset_x < 0) {
		debayer_area->x = 0;
		*debayer_offset_x = area->x;
	} else {
		debayer_area->x = area->x - *debayer_offset_x;
	}

	/* right side */
	int xend = area->x + area->w - 1;
	if (xend & 1)
		right = 10;
	else
		right = 11;
	if (xend + right >= image_area->w) {
		right = image_area->w - xend - 1;
	}
	debayer_area->w = area->w + (area->x - debayer_area->x) + right;

	/* top */
	if (area->y & 1)
		*debayer_offset_y = 11;
	else
		*debayer_offset_y = 10;
	if (area->y - *debayer_offset_y < 0) {
		debayer_area->y = 0;
		*debayer_offset_y = area->y;
	} else {
		debayer_area->y = area->y - *debayer_offset_y;
	}

	/* bottom */
	int yend = area->y + area->h - 1;
	if (yend & 1)
		bottom = 10;
	else
		bottom = 11;
	if (yend + bottom >= image_area->h) {
		bottom = image_area->h - yend - 1;
	}
	debayer_area->h = area->h + (area->y - debayer_area->y) + bottom;

	assert(debayer_area->x < image_area->w);
	assert(debayer_area->y < image_area->h);
	assert(debayer_area->h > 2);
	assert(debayer_area->w > 2);
}


static int debayer_ushort(fits *fit, interpolation_method interpolation, sensor_pattern pattern) {
	size_t npixels = fit->naxes[0] * fit->naxes[1];
	int width = fit->rx;
	int height = fit->ry;
	WORD *buf = fit->data;
	gboolean top_down = com.pref.debayer.top_down;

	unsigned int xtrans[6][6];
	if (interpolation == XTRANS) {
		if (com.pref.debayer.use_bayer_header) {
			if (!g_strcmp0(fit->keywords.row_order, "TOP-DOWN")) {
				top_down = TRUE;
			} else if (!g_strcmp0(fit->keywords.row_order, "BOTTOM-UP")) {
				top_down = FALSE;
			}
		}
		if (!top_down)
			fits_flip_top_to_bottom(fit); // TODO: kind of ugly but not easy with xtrans
		compile_XTrans_pattern(fit->keywords.bayer_pattern, xtrans);
	} else {
		adjust_Bayer_pattern(fit, &pattern);
	}

	if (USE_SIRIL_DEBAYER) {
		WORD *newbuf = debayer_buffer_siril(buf, &width, &height, interpolation, pattern, xtrans, fit->bitpix);
		if (!newbuf)
			return 1;

		fit_debayer_buffer(fit, newbuf);
		// size might have changed, in case of superpixel
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->rx = width;
		fit->ry = height;
		fit->bitpix = fit->orig_bitpix;
		// color RGBRGB format to fits RRGGBB format
		for (size_t i = 0, j = 0; j < npixels; i += 3, j++) {
			WORD r = newbuf[i + RLAYER];
			WORD g = newbuf[i + GLAYER];
			WORD b = newbuf[i + BLAYER];
			fit->pdata[RLAYER][j] = (fit->bitpix == 8) ? truncate_to_BYTE(r) : r;
			fit->pdata[GLAYER][j] = (fit->bitpix == 8) ? truncate_to_BYTE(g) : g;
			fit->pdata[BLAYER][j] = (fit->bitpix == 8) ? truncate_to_BYTE(b) : b;
		}
	} else {
		// use librtprocess debayer
		WORD *newbuf = debayer_buffer_new_ushort(buf, &width, &height, interpolation, pattern, xtrans, fit->bitpix);
		if (!newbuf)
			return 1;

		fit_debayer_buffer(fit, newbuf);
		if (interpolation == XTRANS && !top_down) {
			fits_flip_top_to_bottom(fit);
		}
	}
	/* we remove Bayer header because not needed now */
	clear_Bayer_information(fit);

	/* The data is no longer mono. It's almost certainly linear, but we will
	 * not assign a color profile to it at this stage, it is likely one of
	 * many subs to be sequence processed along with others and it is more
	 * efficient to ignore color managemet until the final stacked image is
	 * available, and the user can then assign a profile as they choose.
	 */
	if (fit->icc_profile)
		cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = NULL;
	color_manage(fit, FALSE);
	return 0;
}

static int debayer_float(fits* fit, interpolation_method interpolation, sensor_pattern pattern) {
	int width = fit->rx;
	int height = fit->ry;
	float *buf = fit->fdata;
	gboolean top_down = com.pref.debayer.top_down;

	unsigned int xtrans[6][6];
	if (interpolation == XTRANS) {
		if (com.pref.debayer.use_bayer_header) {
			if (!g_strcmp0(fit->keywords.row_order, "TOP-DOWN")) {
				top_down = TRUE;
			} else if (!g_strcmp0(fit->keywords.row_order, "BOTTOM-UP")) {
				top_down = FALSE;
			}
		}
		if (!top_down)
			fits_flip_top_to_bottom(fit); // TODO: kind of ugly but not easy with xtrans
		compile_XTrans_pattern(fit->keywords.bayer_pattern, xtrans);
	} else {
		adjust_Bayer_pattern(fit, &pattern);
	}

	float *newbuf = debayer_buffer_new_float(buf, &width, &height, interpolation, pattern, xtrans);
	if (!newbuf)
		return 1;

	fit_debayer_buffer(fit, newbuf);
	if (interpolation == XTRANS && !top_down) {
		fits_flip_top_to_bottom(fit);
	}
	/* we remove Bayer header because not needed now */
	clear_Bayer_information(fit);

	/* The data is no longer mono. It's almost certainly linear, but we will
	 * not assign a color profile to it at this stage, it is likely one of
	 * many subs to be sequence processed along with others and it is more
	 * efficient to ignore color managemet until the final stacked image is
	 * available, and the user can then assign a profile as they choose.
	 */
	if (fit->icc_profile)
		cmsCloseProfile(fit->icc_profile);
	fit->icc_profile = NULL;
	color_manage(fit, FALSE);

	return 0;
}

int debayer(fits *fit, interpolation_method interpolation, sensor_pattern pattern) {
	if (fit->type == DATA_USHORT)
		return debayer_ushort(fit, interpolation, pattern);
	else if (fit->type == DATA_FLOAT)
		return debayer_float(fit, interpolation, pattern);
	else return -1;
}

// gets the index in filter_pattern
int get_cfa_pattern_index_from_string(const char *bayer) {
	for (int i = 0; i < G_N_ELEMENTS(filter_pattern); i++) {
		if (g_ascii_strcasecmp(bayer, filter_pattern[i]) == 0) {
			return i;
		}
	}
	return BAYER_FILTER_NONE;
}

// debayers the image if it's a FITS image and if debayer is activated globally
// or if the force argument is passed
int debayer_if_needed(image_type imagetype, fits *fit, gboolean force_debayer) {
	if ((imagetype != TYPEFITS && imagetype != TYPETIFF && imagetype != TYPEXISF)|| (!com.pref.debayer.open_debayer && !force_debayer))
		return 0;

	/* Siril's FITS are stored bottom-up, debayering will give wrong results.
	 * So before demosaicing we need to flip the image with fits_flip_top_to_bottom().
	 * But sometimes FITS are created by acquisition software top-down, in that case
	 * the user can indicate it ('compatibility') and we don't flip for debayer.
	 */
	if (fit->naxes[2] != 1) {
		siril_log_message(_("Cannot perform debayering on image with more than one channel\n"));
		return 0;
	}

	/* Get Bayer informations from header if available */
	sensor_pattern tmp_pattern = com.pref.debayer.bayer_pattern;
	interpolation_method tmp_algo = com.pref.debayer.bayer_inter;
	if (com.pref.debayer.use_bayer_header) {
		sensor_pattern bayer;
		bayer = get_cfa_pattern_index_from_string(fit->keywords.bayer_pattern);

		if (bayer <= BAYER_FILTER_MAX) {
			if (bayer != tmp_pattern) {
				if (bayer == BAYER_FILTER_NONE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
				}
				else {
					siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
							"blue", filter_pattern[bayer], filter_pattern[com.pref.debayer.bayer_pattern]);
					tmp_pattern = bayer;
				}
			}
		} else {
			tmp_pattern = bayer;
			tmp_algo = XTRANS;
			siril_log_message(_("XTRANS Sensor detected. Using special algorithm.\n"));
		}
	}
	if (tmp_pattern >= BAYER_FILTER_MIN && tmp_pattern <= BAYER_FILTER_MAX) {
		siril_log_message(_("Filter Pattern: %s\n"),
				filter_pattern[tmp_pattern]);
	}

	int retval = 0;
	if (debayer(fit, tmp_algo, tmp_pattern)) {
		siril_log_message(_("Cannot perform debayering\n"));
		retval = -1;
	}
	return retval;
}

static int mergecfa_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_image, MB_avail, required;
	int limit = compute_nb_images_fit_memory(args->seq, 1.0, FALSE, &MB_per_image, NULL, &MB_avail);
	required = 8 * MB_per_image;
	if (limit > 0) {
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
				thread_limit = com.max_thread;
		limit = thread_limit;
		if (for_writer)
			return 1;
	}
	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
#endif
	}
	return limit;
}

void update_bayer_pattern_information(fits *fit, sensor_pattern pattern) {
	switch (pattern) {
		case BAYER_FILTER_RGGB:;
			sprintf(fit->keywords.bayer_pattern, "RGGB");
			break;
		case BAYER_FILTER_BGGR:;
			sprintf(fit->keywords.bayer_pattern, "BGGR");
			break;
		case BAYER_FILTER_GBRG:;
			sprintf(fit->keywords.bayer_pattern, "GBRG");
			break;
		case BAYER_FILTER_GRBG:;
			sprintf(fit->keywords.bayer_pattern, "GRBG");
			break;
		case XTRANS_FILTER_1:;
			sprintf(fit->keywords.bayer_pattern, "GGRGGBGGBGGRBRGRBGGGBGGRGGRGGBRBGBRG");
			break;
		case XTRANS_FILTER_2:;
			sprintf(fit->keywords.bayer_pattern, "RBGBRGGGRGGBGGBGGRBRGRBGGGBGGRGGRGGB");
			break;
		case XTRANS_FILTER_3:;
			sprintf(fit->keywords.bayer_pattern, "GRGGBGBGBRGRGRGGBGGBGGRGRGRBGBGBGGRG");
			break;
		case XTRANS_FILTER_4:;
			sprintf(fit->keywords.bayer_pattern, "GBGGRGRGRBGBGBGGRGGRGGBGBGBRGRGRGGBG");
			break;
		default:;
			break;
	}
}

gint64 mergecfa_compute_size_hook(struct generic_seq_args *args, int nb_frames) {
	double ratio = 4.;
	double fullseqsize = seq_compute_size(args->seq, nb_frames, args->output_type);
	return (gint64)(fullseqsize * ratio);
}

int mergecfa_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads) {
	fits metadata = { 0 };
	if (seq_read_frame_metadata(args->seq, out_index, &metadata)) {
		siril_log_message(_("Could not load metadata\n"));
		return 1;
	}
	int retval = 0;
	struct merge_cfa_data *merge_cfa_args = (struct merge_cfa_data*) args->user;
	fits cfa1 = { 0 };
	fits cfa2 = { 0 };
	fits cfa3 = { 0 };
	fits *out = { 0 };
	retval = seq_read_frame(merge_cfa_args->seq1, out_index, &cfa1, args->seq->bitpix == 16, -1);
	if(retval != 0) {
		siril_log_message(_("Image %d: error opening CFA1 file\n"), args->seq->current);
		goto CLEANUP_MERGECFA;
	}
	retval = seq_read_frame(merge_cfa_args->seq2, out_index, &cfa2, args->seq->bitpix == 16, -1);
	if(retval != 0) {
		siril_log_message(_("Image %d: error opening CFA2 file\n"), args->seq->current);
		goto CLEANUP_MERGECFA;
	}
	retval = seq_read_frame(merge_cfa_args->seq3, out_index, &cfa3, args->seq->bitpix == 16, -1);
	if(retval != 0) {
		siril_log_message(_("Image %d: error opening CFA3 file\n"), args->seq->current);
		goto CLEANUP_MERGECFA;
	}
	out = merge_cfa(fit, &cfa1, &cfa2, &cfa3, merge_cfa_args->pattern);
	if (out != NULL) {
		clearfits(fit);
		copyfits(out, fit, (CP_ALLOC | CP_COPYA | CP_FORMAT), -1);
		clearfits(out);
		free(out);
	}
	copy_fits_metadata(&metadata, fit);
	update_sampling_information(fit, 0.5f);
	update_bayer_pattern_information(fit, merge_cfa_args->pattern);
	clearfits(&metadata);

CLEANUP_MERGECFA:
	clearfits(&cfa1);
	clearfits(&cfa2);
	clearfits(&cfa3);

	return retval;
}

int mergecfa_finalize_hook(struct generic_seq_args *args) {
	struct merge_cfa_data *data = (struct merge_cfa_data *) args->user;
	int retval = seq_finalize_hook(args);
	if (data->seq0 != args->seq && !check_seq_is_comseq(data->seq0)) {
		free_sequence(data->seq0, TRUE);
		data->seq0 = NULL;
	}
	if (data->seq1 != args->seq && !check_seq_is_comseq(data->seq1)) {
		free_sequence(data->seq1, TRUE);
		data->seq1 = NULL;
	}
	if (data->seq2 != args->seq && !check_seq_is_comseq(data->seq2)) {
		free_sequence(data->seq2, TRUE);
		data->seq2 = NULL;
	}
	if (data->seq3 != args->seq && !check_seq_is_comseq(data->seq3)) {
		free_sequence(data->seq3, TRUE);
		data->seq3 = NULL;
	}
	free(data);
	return retval;
}

void apply_mergecfa_to_sequence(struct merge_cfa_data *merge_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(merge_cfa_args->seq0);
	args->seq = merge_cfa_args->seq0;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = merge_cfa_args->seq0->selnum;
	args->compute_mem_limits_hook = mergecfa_compute_mem_limits;
	args->compute_size_hook = mergecfa_compute_size_hook;
	args->prepare_hook = seq_prepare_hook;
	args->image_hook = mergecfa_image_hook;
	args->finalize_hook = mergecfa_finalize_hook;
	args->description = _("Merge CFA");
	args->has_output = TRUE;
	args->new_seq_prefix = strdup(merge_cfa_args->seqEntryOut);
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = merge_cfa_args;

	if (!start_in_new_thread(generic_sequence_worker, args)) {
		if (!check_seq_is_comseq(merge_cfa_args->seq0))
			free_sequence(merge_cfa_args->seq0, TRUE);
		if (!check_seq_is_comseq(merge_cfa_args->seq1))
			free_sequence(merge_cfa_args->seq1, TRUE);
		if (!check_seq_is_comseq(merge_cfa_args->seq2))
			free_sequence(merge_cfa_args->seq2, TRUE);
		if (!check_seq_is_comseq(merge_cfa_args->seq3))
			free_sequence(merge_cfa_args->seq3, TRUE);
		free(merge_cfa_args);
		free_generic_seq_args(args, TRUE);
	}
}

//
// Re-mosaic 4 subpattern images previously separated with split_cfa
// This routine is Bayer pattern agnostic: split_cfa generates files
// CFA0, CFA1, CFA2 and CFA3 which can be RGGB, BGGR etc. This routine
// will assemble them in the same pattern as long as the files are
// provided in the correct order CFA0 to CFA3.
//
// (returns 0 on success, -1 on failure)
//

fits* merge_cfa (fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3, sensor_pattern pattern) {
	fits *out = NULL;

	// Check input files are compatible
	gboolean x_compat = (cfa0->naxes[0] == cfa1->naxes[0] && cfa1->naxes[0] == cfa2->naxes[0] && cfa2->naxes[0] == cfa3->naxes[0]);
	gboolean y_compat = (cfa0->naxes[1] == cfa1->naxes[1] && cfa1->naxes[1] == cfa2->naxes[1] && cfa2->naxes[1] == cfa3->naxes[1]);
	gboolean c_compat = (cfa0->naxes[2] == cfa1->naxes[2] && cfa1->naxes[2] == cfa2->naxes[2] && cfa2->naxes[2] == cfa3->naxes[2] && cfa3->naxes[2] == 1);
	gboolean t_compat = (cfa0->type == cfa1->type && cfa1->type == cfa2->type && cfa2->type == cfa3->type);
	if (!(x_compat && y_compat && c_compat && t_compat)) {
		siril_log_color_message(_("Input files are incompatible (all must be mono with the same size and bit depth). Aborting...\n"), "red");
		if(!x_compat)
			siril_log_message(_("X dimensions incompatible\n"));
		if(!y_compat)
			siril_log_message(_("Y dimensions incompatible\n"));
		if(!c_compat)
			siril_log_message(_("Channels not all mono\n"));
		if(!t_compat)
			siril_log_message(_("Input files not all the same bit depth\n"));
		return NULL;
	}
	int datatype = cfa0->type;

	// Create output fits twice the width and height of the cfa fits files
	if (new_fit_image(&out, cfa0->rx << 1, cfa0->ry << 1, 1, datatype)) {
		siril_log_color_message(_("Error creating output image\n"), "red");
		return NULL;
	}

//	out->header = copy_header(cfa0);
	copy_fits_metadata(cfa0, out);
	if (out->header)
		fprintf(stdout, "header ok\n");
	for (size_t outx = 0 ; outx < out->rx; outx += 2) {
		for(size_t outy = 0 ; outy < out->ry ; outy += 2) {
			size_t cfax = outx >> 1;
			size_t cfay = outy >> 1;
			size_t indexcfa = cfax + cfay * (size_t) cfa0->rx;
			// Note order because of the read orientation
			size_t indexout1 = outx + outy * out->rx;
			size_t indexout3 = (outx + 1) + outy * out->rx;
			size_t indexout0 = outx + (outy + 1) * out->rx;
			size_t indexout2 = (outx + 1) + (outy + 1) * out->rx;
			switch (datatype) {
				case DATA_FLOAT:
					out->fdata[indexout0] = cfa0->fdata[indexcfa];
					out->fdata[indexout1] = cfa1->fdata[indexcfa];
					out->fdata[indexout2] = cfa2->fdata[indexcfa];
					out->fdata[indexout3] = cfa3->fdata[indexcfa];
					break;
				case DATA_USHORT:
					out->data[indexout0] = cfa0->data[indexcfa];
					out->data[indexout1] = cfa1->data[indexcfa];
					out->data[indexout2] = cfa2->data[indexcfa];
					out->data[indexout3] = cfa3->data[indexcfa];
					break;
			}
		}
	}
	// Set Bayer pattern in FITS header
	switch (pattern) {
		case BAYER_FILTER_RGGB:
			strcpy(out->keywords.bayer_pattern, "RGGB");
			break;
		case BAYER_FILTER_BGGR:
			strcpy(out->keywords.bayer_pattern, "BGGR");
			break;
		case BAYER_FILTER_GRBG:
			strcpy(out->keywords.bayer_pattern, "GRBG");
			break;
		case BAYER_FILTER_GBRG:
			strcpy(out->keywords.bayer_pattern, "GBRG");
			break;
		default:
			break;
	}
	clearfits (cfa0);
	clearfits (cfa1);
	clearfits (cfa2);
	clearfits (cfa3);
	siril_debug_print("Merge CFA complete\n");
	return out;
}

static unsigned int compile_bayer_pattern(sensor_pattern pattern) {
	/* red is 0, green is 1, blue is 2 */
	switch (pattern) {
		case BAYER_FILTER_RGGB:
			return cpu_to_be32(0x00010102);
		case BAYER_FILTER_BGGR:
			return cpu_to_be32(0x02010100);
		case BAYER_FILTER_GBRG:
			return cpu_to_be32(0x01020001);
		case BAYER_FILTER_GRBG:
			return cpu_to_be32(0x01000201);
		default:
			return 0x03030303;
	}
}

/* from the header description of the color filter array, we create a mask that
 * contains filter values for top-down reading */
static int get_compiled_pattern(fits *fit, int layer, BYTE pattern[36], int *pattern_size) {
	g_assert(layer >= 0 && layer < 3);
	sensor_pattern idx = get_cfa_pattern_index_from_string(fit->keywords.bayer_pattern);

	if (idx >= BAYER_FILTER_MIN && idx <= BAYER_FILTER_MAX) {
		/* 2x2 Bayer matrix */
		// reorient it correctly for the image
		adjust_Bayer_pattern(fit, &idx);
		*((unsigned int *)pattern) = compile_bayer_pattern(idx);
		*pattern_size = 2;
		siril_debug_print("extracting CFA data for filter %d on Bayer pattern\n", layer);
		return 0;
	}
	else if (idx != BAYER_FILTER_NONE) {
		/* 6x6 X-Trans matrix */
		unsigned int xtrans[6][6];
		compile_XTrans_pattern(fit->keywords.bayer_pattern, xtrans);
		adjust_XTrans_pattern(fit, xtrans);
		for (int i = 0; i < 36; i++)
			pattern[i] = (BYTE)(((unsigned int *)xtrans)[i]);
		*pattern_size = 6;
		siril_debug_print("extracting CFA data for filter %d on X-Trans pattern\n", layer);
		return 0;
	}
	return 1;
}

WORD *extract_CFA_buffer_ushort(fits *fit, int layer, size_t *newsize) {
	BYTE pattern[36];	// red is 0, green is 1, blue is 2
	int pattern_size;	// 2 or 6
	if (get_compiled_pattern(fit, layer, pattern, &pattern_size))
		return NULL;

	// alloc buffer
	size_t npixels = fit->naxes[0] * fit->naxes[1];
	WORD *buf = malloc(npixels * sizeof(WORD));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	WORD *input = fit->data;

	// and we copy the pixels top-down
	size_t i = 0, j = 0;
	int x, y, pattern_y = 0;
	for (y = 0; y < fit->ry; y++) {
		int pattern_idx_y = pattern_y * pattern_size;
		for (x = 0; x < fit->rx; x++) {
			int pattern_idx = pattern_idx_y + x % pattern_size;
			if (pattern[pattern_idx] == layer)
				buf[j++] = input[i];
			i++;
		}
		pattern_y = (pattern_y + 1) % pattern_size;
	}

	void * tmp = realloc(buf, j * sizeof(WORD));
	if (!tmp) {
		PRINT_ALLOC_ERR;
		free(buf);
		return NULL;
	} else {
		buf = tmp;
	}
	*newsize = j;
	return buf;
}

WORD *extract_CFA_buffer_area_ushort(fits *fit, int layer, rectangle *bounds, size_t *newsize) {
	BYTE pattern[36];	// red is 0, green is 1, blue is 2
	int pattern_size;	// 2 or 6
	if (get_compiled_pattern(fit, layer, pattern, &pattern_size))
		return NULL;
	siril_debug_print("CFA buffer extraction with area\n");

	// alloc buffer
	size_t npixels = bounds->w * bounds->h;
	WORD *buf = malloc(npixels * sizeof(WORD));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	WORD *input = fit->data;

	// copy the pixels top-down
	size_t j = 0;
	int start_y = fit->ry - bounds->y - bounds-> h;
	int pattern_y = start_y % pattern_size;
	for (int y = start_y; y < fit->ry - bounds->y; y++) {
		int pattern_idx_y = pattern_y * pattern_size;
		size_t i = y * fit->rx + bounds->x;
		for (int x = bounds->x; x < bounds->x + bounds->w; x++) {
			int pattern_idx = pattern_idx_y + x % pattern_size;
			if (pattern[pattern_idx] == layer)
				buf[j++] = input[i];
			i++;
		}
		pattern_y = (pattern_y + 1) % pattern_size;
	}

	void *tmp = realloc(buf, j * sizeof(WORD));
	if (!tmp) {
		free(buf);
		PRINT_ALLOC_ERR;
		return NULL;
	} else {
		buf = tmp;
	}
	*newsize = j;
	return buf;
}

float *extract_CFA_buffer_float(fits *fit, int layer, size_t *newsize) {
	BYTE pattern[36];	// red is 0, green is 1, blue is 2
	int pattern_size;	// 2 or 6
	if (get_compiled_pattern(fit, layer, pattern, &pattern_size))
		return NULL;

	// alloc buffer
	size_t npixels = fit->naxes[0] * fit->naxes[1];
	float *buf = malloc(npixels * sizeof(float));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	float *input = fit->fdata;

	// copy the pixels top-down
	size_t i = 0, j = 0;
	int x, y, pattern_y = 0;
	for (y = 0; y < fit->ry; y++) {
		int pattern_idx_y = pattern_y * pattern_size;
		for (x = 0; x < fit->rx; x++) {
			int pattern_idx = pattern_idx_y + x % pattern_size;
			if (pattern[pattern_idx] == layer)
				buf[j++] = input[i];
			i++;
		}
		pattern_y = (pattern_y + 1) % pattern_size;
	}

	void *tmp = realloc(buf, j * sizeof(float));
	if (!tmp) {
		free(buf);
		PRINT_ALLOC_ERR;
		return NULL;
	} else {
		buf = tmp;
	}
	*newsize = j;
	return buf;
}

float *extract_CFA_buffer_area_float(fits *fit, int layer, rectangle *bounds, size_t *newsize) {
	BYTE pattern[36];	// red is 0, green is 1, blue is 2
	int pattern_size;	// 2 or 6
	if (get_compiled_pattern(fit, layer, pattern, &pattern_size))
		return NULL;
	siril_debug_print("CFA buffer extraction with area\n");

	// alloc buffer
	size_t npixels = bounds->w * bounds->h;
	float *buf = malloc(npixels * sizeof(float));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	float *input = fit->fdata;

	// copy the pixels top-down
	size_t j = 0;
	int start_y = fit->ry - bounds->y - bounds-> h;
	int pattern_y = start_y % pattern_size;
	for (int y = start_y; y < fit->ry - bounds->y; y++) {
		int pattern_idx_y = pattern_y * pattern_size;
		size_t i = y * fit->rx + bounds->x;
		for (int x = bounds->x; x < bounds->x + bounds->w; x++) {
			int pattern_idx = pattern_idx_y + x % pattern_size;
			if (pattern[pattern_idx] == layer)
				buf[j++] = input[i];
			i++;
		}
		pattern_y = (pattern_y + 1) % pattern_size;
	}

	void *tmp = realloc(buf, j * sizeof(float));
	if (!tmp) {
		free(buf);
		PRINT_ALLOC_ERR;
		return NULL;
	} else {
		buf = tmp;
	}
	*newsize = j;
	return buf;
}
