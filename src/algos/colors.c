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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "algos/cie_standard_observer.h" // Do not include this from anywhere else
#include "core/siril.h"
#include "core/proto.h"
#include "core/icc_profile.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "gui/progress_and_log.h"
#include "gui/histogram.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "algos/extraction.h"

/******************************************************************************
 * Note for maintainers: do not use the translation macro on the following    *
 * string. Color management relies on being able to detect "Extraction"       *
 * in FITS HISTORY header.                                                    *
 ******************************************************************************/
const gchar *extractionstring = "Extraction";

static gchar *add_filter_str[] = { "R", "G", "B"};
/*
 * A Fast HSL-to-RGB Transform
 * by Ken Fishkin
 * from "Graphics Gems", Academic Press, 1990
 * */
/*
 *  * given h,s,l on [0..1],
 *   * return r,g,b on [0..1]
 *    */
void hsl_to_rgb_float_sat(float h, float sl, float l, float * r, float * g,
		float * b) {
	float v;

	h = h >= 6.f ? h - 6.f : h;

	v = (l <= 0.5f) ? (l * (1.f + sl)) : (l + sl - l * sl);
	if (v <= 0.f) {
		*r = *g = *b = 0.f;
	} else {
		float m;
		float sv;
		int sextant;
		float fract, vsf, mid1, mid2;

		m = l + l - v;
		sv = (v - m) / v;
		sextant = h;
		fract = h - sextant;
		vsf = v * sv * fract;
		mid1 = m + vsf;
		mid2 = v - vsf;
		switch (sextant) {
			case 0:
				*r = v;
				*g = mid1;
				*b = m;
				break;
			case 1:
				*r = mid2;
				*g = v;
				*b = m;
				break;
			case 2:
				*r = m;
				*g = v;
				*b = mid1;
				break;
			case 3:
				*r = m;
				*g = mid2;
				*b = v;
				break;
			case 4:
				*r = mid1;
				*g = m;
				*b = v;
				break;
			case 5:
				*r = v;
				*g = m;
				*b = mid2;
				break;
		}
	}
}
/*
 *  * RGB-HSL transforms.
 *   * Ken Fishkin, Pixar Inc., January 1989.
 *    */

/*
 *  * given r,g,b on [0 ... 1],
 *   * return (h,s,l) on [0 ... 1]
 *    */
void rgb_to_hsl_float_sat(float r, float g, float b, float low, float *h, float *s, float *l) {
	float v;
	float m;
	float vm;
	float r2, g2, b2;

	v = max(r, g);
	v = max(v, b);
	m = min(r, g);
	m = min(m, b);

	if (m + v < low + low) {
		*l = 0.f;
		return;
	}
	*l = (m + v) / 2.f;
	*h = 0.f;
	*s = 0.f;	// init values

	if ((*s = vm = v - m) > 0.f) {
		*s /= (*l <= 0.5f) ? (v + m) : (2.f - v - m);
	} else
		return;

	if (r == v) {
		g2 = (v - g) / vm;
		b2 = (v - b) / vm;
		*h = (g == m ? 5.f + b2 : 1.f - g2);
	}else if (g == v) {
		r2 = (v - r) / vm;
		b2 = (v - b) / vm;
		*h = (b == m ? 1.f + r2 : 3.f - b2);
	} else {
		r2 = (v - r) / vm;
		g2 = (v - g) / vm;
		*h = (r == m ? 3.f + g2 : 5.f - r2);
	}

}

/*
 * A Fast HSL-to-RGB Transform
 * by Ken Fishkin
 * from "Graphics Gems", Academic Press, 1990
 * */
/*
 *  * given h,s,l on [0..1],
 *   * return r,g,b on [0..1]
 *    */
void hsl_to_rgb(double h, double sl, double l, double * r, double * g,
		double * b) {
	double v;

	assert(h >= 0.0 && h <= 1.0);
	if (h >= 1.0) h -= 1.0;		// this code doesn't work for h = 1
	v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);
	if (v <= 0) {
		*r = *g = *b = 0.0;
	} else {
		double m;
		double sv;
		int sextant;
		double fract, vsf, mid1, mid2;

		m = l + l - v;
		sv = (v - m) / v;
		h *= 6.0;
		sextant = h;
		fract = h - sextant;
		vsf = v * sv * fract;
		mid1 = m + vsf;
		mid2 = v - vsf;
		switch (sextant) {
			case 0:
				*r = v;
				*g = mid1;
				*b = m;
				break;
			case 1:
				*r = mid2;
				*g = v;
				*b = m;
				break;
			case 2:
				*r = m;
				*g = v;
				*b = mid1;
				break;
			case 3:
				*r = m;
				*g = mid2;
				*b = v;
				break;
			case 4:
				*r = mid1;
				*g = m;
				*b = v;
				break;
			case 5:
				*r = v;
				*g = m;
				*b = mid2;
				break;
		}
	}
}

/* RGB-HSL transforms.
 * Ken Fishkin, Pixar Inc., January 1989.
 *
 * given r,g,b on [0 ... 1],
 * return (h,s,l) on [0 ... 1]
 */
void rgb_to_hsl(double r, double g, double b, double *h, double *s, double *l) {
	double v;
	double m;
	double vm;
	double r2, g2, b2;

	v = max(r, g);
	v = max(v, b);
	m = min(r, g);
	m = min(m, b);
	*h = 0.0;
	*s = 0.0;	// init values

	if ((*l = (m + v) / 2.0) <= 0.0) {
		*l = 0.0;
		return;
	}
	if ((*s = vm = v - m) > 0.0) {
		*s /= (*l <= 0.5) ? (v + m) : (2.0 - v - m);
	} else
		return;

	r2 = (v - r) / vm;
	g2 = (v - g) / vm;
	b2 = (v - b) / vm;

	if (r == v)
		*h = (g == m ? 5.0 + b2 : 1.0 - g2);
	else if (g == v)
		*h = (b == m ? 1.0 + r2 : 3.0 - b2);
	else
		*h = (r == m ? 3.0 + g2 : 5.0 - r2);

	*h /= 6;
}

// Single precision versions of the above
void hsl_to_rgbf(float h, float s, float l, float * r, float * g,
		float * b) {
	float v;

	assert(h >= 0.f && h <= 1.f);
	if (h >= 1.f) h -= 1.f;		// this code doesn't work for h = 1
	v = (l <= 0.5f) ? (l * (1.f + s)) : (l + s - l * s);
	if (v <= 0.f) {
		*r = *g = *b = 0.f;
	} else {
		float m;
		float sv;
		int sextant;
		float fract, vsf, mid1, mid2;

		m = l + l - v;
		sv = (v - m) / v;
		h *= 6.f;
		sextant = h;
		fract = h - sextant;
		vsf = v * sv * fract;
		mid1 = m + vsf;
		mid2 = v - vsf;
		switch (sextant) {
			case 0:
				*r = v;
				*g = mid1;
				*b = m;
				break;
			case 1:
				*r = mid2;
				*g = v;
				*b = m;
				break;
			case 2:
				*r = m;
				*g = v;
				*b = mid1;
				break;
			case 3:
				*r = m;
				*g = mid2;
				*b = v;
				break;
			case 4:
				*r = mid1;
				*g = m;
				*b = v;
				break;
			case 5:
				*r = v;
				*g = m;
				*b = mid2;
				break;
		}
	}
}

void rgb_to_hslf(float r, float g, float b, float *h, float *s, float *l) {
	float v;
	float m;
	float vm;
	float r2, g2, b2;

	v = max(r, g);
	v = max(v, b);
	m = min(r, g);
	m = min(m, b);
	*h = 0.f;
	*s = 0.f;	// init values

	if ((*l = (m + v) / 2.f) <= 0.f) {
		*l = 0.f;
		return;
	}
	if ((*s = vm = v - m) > 0.f) {
		*s /= (*l <= 0.5f) ? (v + m) : (2.f - v - m);
	} else
		return;

	r2 = (v - r) / vm;
	g2 = (v - g) / vm;
	b2 = (v - b) / vm;

	if (r == v)
		*h = (g == m ? 5.f + b2 : 1.f - g2);
	else if (g == v)
		*h = (b == m ? 1.f + r2 : 3.f - b2);
	else
		*h = (r == m ? 3.f + g2 : 5.f - r2);

	*h /= 6.f;
}

/* all variables are between 0 and 1. h takes 0 for grey */
void rgb_to_hsv(double r, double g, double b, double *h, double *s, double *v) {
	double cmax, cmin, delta;

	cmax = max(r, g);
	cmax = max(cmax, b);
	cmin = min(r, g);
	cmin = min(cmin, b);
	delta = cmax - cmin;
	*v = cmax;
	if (delta == 0.0) {
		*s = 0.0;
		*h = 0.0;
		return;
	}
	*s = delta / cmax;

	if (cmax == r)
		*h = (((g - b) / delta)) / 6.0;
	else if (cmax == g)
		*h = (((b - r) / delta) + 2.0) / 6.0;
	else
		*h = (((r - g) / delta) + 4.0) / 6.0;

	if (*h < 0.0)
		*h += 1.0;
}

void hsv_to_rgb(double h, double s, double v, double *r, double *g, double *b) {
	double p, q, t, f;
	int i;

	if (h >= 1.0)
		h -= 1.0;
	h *= 6.0;
	i = (int)h;
	f = h - (double)i;
	p = v * (1.0 - s);
	q = v * (1.0 - (s * f));
	t = v * (1.0 - (s * (1.0 - f)));

	switch (i) {
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;
		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		case 5:
		default:
			*r = v;
			*g = p;
			*b = q;
			break;
	}
}

void rgb_to_xyz(double r, double g, double b, double *x, double *y, double *z) {
	r = (r <= 0.04045) ? r / 12.92 : pow(((r + 0.055) / 1.055), 2.4);
	g = (g <= 0.04045) ? g / 12.92 : pow(((g + 0.055) / 1.055), 2.4);
	b = (b <= 0.04045) ? b / 12.92 : pow(((b + 0.055) / 1.055), 2.4);

	r *= 100;
	g *= 100;
	b *= 100;

	*x = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b;
	*y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b;
	*z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b;
}

void rgb_to_xyzf(float r, float g, float b, float *x, float *y, float *z) {
	r = (r <= 0.04045f) ? r / 12.92f : powf(((r + 0.055f) / 1.055f), 2.4f);
	g = (g <= 0.04045f) ? g / 12.92f : powf(((g + 0.055f) / 1.055f), 2.4f);
	b = (b <= 0.04045f) ? b / 12.92f : powf(((b + 0.055f) / 1.055f), 2.4f);

	r *= 100.f;
	g *= 100.f;
	b *= 100.f;

	*x = 0.4124564f * r + 0.3575761f * g + 0.1804375f * b;
	*y = 0.2126729f * r + 0.7151522f * g + 0.0721750f * b;
	*z = 0.0193339f * r + 0.1191920f * g + 0.9503041f * b;
}

/* These functions have a scale argument as sometimes it is convenient to have
 * the Y channel output in the range [0..1] instead of [0..100]
 */

void linrgb_to_xyz(double r, double g, double b, double *x, double *y, double *z, gboolean scale) {
	if (scale) {
		r *= 100;
		g *= 100;
		b *= 100;
	}

	*x = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b;
	*y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b;
	*z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b;
}

void linrgb_to_xyzf(float r, float g, float b, float *x, float *y, float *z, gboolean scale) {
	if (scale) {
		r *= 100.f;
		g *= 100.f;
		b *= 100.f;
	}

	*x = 0.4124564f * r + 0.3575761f * g + 0.1804375f * b;
	*y = 0.2126729f * r + 0.7151522f * g + 0.0721750f * b;
	*z = 0.0193339f * r + 0.1191920f * g + 0.9503041f * b;
}

void xyz_to_linrgb(double x, double y, double z, double *r, double *g, double *b, gboolean scale) {
	if (scale) {
		x /= 100.0;
		y /= 100.0;
		z /= 100.0;
	}

	*r =  3.2404542 * x - 1.5371385 * y - 0.4985314 * z;
	*g = -0.9692660 * x + 1.8760108 * y + 0.0415560 * z;
	*b =  0.0556434 * x - 0.2040259 * y + 1.0572252 * z;
}

void xyz_to_linrgbf(float x, float y, float z, float *r, float *g, float *b, gboolean scale) {
	if (scale) {
		x /= 100.f;
		y /= 100.f;
		z /= 100.f;
	}

	*r =  3.2404542f * x - 1.5371385f * y - 0.4985314f * z;
	*g = -0.9692660f * x + 1.8760108f * y + 0.0415560f * z;
	*b =  0.0556434f * x - 0.2040259f * y + 1.0572252f * z;
}

void rgb_to_yuvf(float red, float green, float blue, float *y, float *u, float *v) {
	const float a = 1.f / sqrtf(3.f);
	const float b = 1.f / sqrtf(2.f);
	const float c = 2.f * a * sqrtf(2.f);
	*y = a * (red + green + blue);
	*u = b * (red - blue);
	*v = c * (0.25f * red - 0.5f * green + 0.25f * blue);
}

void yuv_to_rgbf(float y, float u, float v, float *red, float *green, float *blue) {
	const float a = 1.f / sqrtf(3.f);
	const float b = 1.f / sqrtf(2.f);
	const float c = a / b;
	*red = (a * y) + (b * u) + (c * 0.5f * v);
	*green = (a * y) - (c * v);
	*blue = (a * y) - (b * u) + (c * 0.5f * v);
}


void xyz_to_LAB(double x, double y, double z, double *L, double *a, double *b) {
	x /= 95.047;
	y /= 100.000;
	z /= 108.883;

	x = (x > 0.008856452) ? pow(x, 1 / 3.0) : (7.787037037 * x) + (16. / 116.);
	y = (y > 0.008856452) ? pow(y, 1 / 3.0) : (7.787037037 * y) + (16. / 116.);
	z = (z > 0.008856452) ? pow(z, 1 / 3.0) : (7.787037037 * z) + (16. / 116.);

	*L = (116.0 * y) - 16.0;
	*a = 500.0 * (x - y);
	*b = 200.0 * (y - z);
}

void xyz_to_LABf(float x, float y, float z, float *L, float *a, float *b) {
	x /= 95.047f;
	y /= 100.000f;
	z /= 108.883f;

	x = (x > 0.008856452f) ? powf(x, 1.f / 3.0f) : (7.787037037f * x) + (16.f / 116.f);
	y = (y > 0.008856452f) ? powf(y, 1.f / 3.0f) : (7.787037037f * y) + (16.f / 116.f);
	z = (z > 0.008856452f) ? powf(z, 1.f / 3.0f) : (7.787037037f * z) + (16.f / 116.f);

	*L = (116.0f * y) - 16.0f;
	*a = 500.0f * (x - y);
	*b = 200.0f * (y - z);
}

void LAB_to_xyz(double L, double a, double b, double *x, double *y, double *z) {
	*y = (L + 16.0) / 116.0;
	*x = a / 500.0 + (*y);
	*z = *y - b / 200.0;

	double x3, y3, z3;
	x3 = (*x) * (*x) * (*x);
	y3 = (*y) * (*y) * (*y);
	z3 = (*z) * (*z) * (*z);

	*x = (x3 > 0.008856452) ? x3 : (*x - 16. / 116.) / 7.787037037;
	*y = (y3 > 0.008856452) ? y3 : (*y - 16. / 116.) / 7.787037037;
	*z = (z3 > 0.008856452) ? z3 : (*z - 16. / 116.) / 7.787037037;

	*x *= 95.047;
	*y *= 100.000;
	*z *= 108.883;
}

void LAB_to_xyzf(float L, float a, float b, float *x, float *y, float *z) {
	*y = (L + 16.0f) / 116.0f;
	*x = a / 500.0f + (*y);
	*z = *y - b / 200.0f;

	float x3, y3, z3;
	x3 = (*x) * (*x) * (*x);
	y3 = (*y) * (*y) * (*y);
	z3 = (*z) * (*z) * (*z);

	*x = (x3 > 0.008856452f) ? x3 : (*x - 16.f / 116.f) / 7.787037037f;
	*y = (y3 > 0.008856452f) ? y3 : (*y - 16.f / 116.f) / 7.787037037f;
	*z = (z3 > 0.008856452f) ? z3 : (*z - 16.f / 116.f) / 7.787037037f;

	*x *= 95.047f;
	*y *= 100.000f;
	*z *= 108.883f;
}

void xyz_to_rgb(double x, double y, double z, double *r, double *g, double *b) {
	x /= 100.0;
	y /= 100.0;
	z /= 100.0;

	*r =  3.2404542 * x - 1.5371385 * y - 0.4985314 * z;
	*g = -0.9692660 * x + 1.8760108 * y + 0.0415560 * z;
	*b =  0.0556434 * x - 0.2040259 * y + 1.0572252 * z;

	*r = (*r > 0.0031308) ? 1.055 * (pow(*r, (1 / 2.4))) - 0.055 : 12.92 * (*r);
	*g = (*g > 0.0031308) ? 1.055 * (pow(*g, (1 / 2.4))) - 0.055 : 12.92 * (*g);
	*b = (*b > 0.0031308) ? 1.055 * (pow(*b, (1 / 2.4))) - 0.055 : 12.92 * (*b);

}

void xyz_to_rgbf(float x, float y, float z, float *r, float *g, float *b) {
	x /= 100.0f;
	y /= 100.0f;
	z /= 100.0f;

	*r =  3.2404542f * x - 1.5371385f * y - 0.4985314f * z;
	*g = -0.9692660f * x + 1.8760108f * y + 0.0415560f * z;
	*b =  0.0556434f * x - 0.2040259f * y + 1.0572252f * z;

	*r = (*r > 0.0031308f) ? 1.055f * (powf(*r, (1.f / 2.4f))) - 0.055f : 12.92f * (*r);
	*g = (*g > 0.0031308f) ? 1.055f * (powf(*g, (1.f / 2.4f))) - 0.055f : 12.92f * (*g);
	*b = (*b > 0.0031308f) ? 1.055f * (powf(*b, (1.f / 2.4f))) - 0.055f : 12.92f * (*b);
}

// Reference: https://en.wikipedia.org/wiki/Color_index and https://arxiv.org/abs/1201.1809 (Ballesteros, F. J., 2012)
// Uses Ballesteros' formula based on considering stars as black bodies
double BV_to_T(double BV) {
	double T;

	// make sure BV is within its bounds [-0.4, 2] otherwise the math doesnt work
	if (BV < -0.4) {
		BV = -0.4;
	} else if (BV > 2) {
		BV = 2;
	}

	// http://www.wikiwand.com/en/Color_index
	T = 4600 * ((1 / ((0.92 * BV) + 1.7)) + (1 / ((0.92 * BV) + 0.62)));

	return T;
}

// CIE XYZ Color Matching Functions
float x1931(float w) {
	int index = w - 360;
	float w2_w = w - (float) w;
	float w1_w = 1.f - w2_w;
	index = max(min(index, 470), 0);
	float retval = xyz1931_1nm[index][0] * w1_w + xyz1931_1nm[index+1][0] * w2_w;
	return retval;
}

float y1931(float w) {
	int index = w - 360;
	float w2_w = w - (float) w;
	float w1_w = 1.f - w2_w;
	index = max(min(index, 470), 0);
	float retval = xyz1931_1nm[index][1] * w1_w + xyz1931_1nm[index+1][1] * w2_w;
	return retval;
}

float z1931(float w) {
	int index = w - 360;
	float w2_w = w - (float) w;
	float w1_w = 1.f - w2_w;
	index = max(min(index, 470), 0);
	float retval = xyz1931_1nm[index][2] * w1_w + xyz1931_1nm[index+1][2] * w2_w;
	return retval;
}

float x1964(float w) {
	int index = w - 360;
	float w2_w = w - (float) w;
	float w1_w = 1.f - w2_w;
	index = max(min(index, 470), 0);
	float retval = xyz1964_1nm[index][0] * w1_w + xyz1964_1nm[index+1][0] * w2_w;
	return retval;
}

float y1964(float w) {
	int index = w - 360;
	float w2_w = w - (float) w;
	float w1_w = 1.f - w2_w;
	index = max(min(index, 470), 0);
	float retval = xyz1964_1nm[index][1] * w1_w + xyz1964_1nm[index+1][1] * w2_w;
	return retval;
}

float z1964(float w) {
	int index = w - 360;
	float w2_w = w - (float) w;
	float w1_w = 1.f - w2_w;
	index = max(min(index, 470), 0);
	float retval = xyz1964_1nm[index][2] * w1_w + xyz1964_1nm[index+1][2] * w2_w;
	return retval;
}

cmsCIExyY wl_to_xyY(double wl) {
	cmsCIEXYZ XYZ;
	cmsCIExyY xyY;
	if (com.pref.icc.cmf == CMF_1931_2DEG) {
		XYZ = (cmsCIEXYZ) { x1931(wl), y1931(wl), z1931(wl) };
	} else {
		XYZ = (cmsCIEXYZ) { x1964(wl), y1964(wl), z1964(wl) };
	}
	cmsXYZ2xyY(&xyY, &XYZ);
	return xyY;
}

// Returns the emittance of a Planckian black body spectrum at wavelength
// wl and temperature bbTemp
// from "Colour Rendering of Spectra", John Walker, Fourmilab. Public domain code, last updated March 9 2003
/*
float bb_spectrum(float wl, float bbTemp) {
	float wlm = wl * 1e-9;   // Wavelength in meters
	return (3.74183e-16f * pow(wlm, -5.f)) / (expf(1.4388e-2f / (wlm * bbTemp)) - 1.f);
}
*/

int equalize_cfa_fit_with_coeffs(fits *fit, float coeff1, float coeff2, const char *cfa_string) {
	unsigned int row, col, pat_cell;
	/* compute width of the (square) CFA pattern */
	/* added 0.1 in case the result of sqrt is something like 5.999999 and it gets casted to int as 5 */
	unsigned int pat_width = (unsigned int) (sqrt(strlen(cfa_string)) + 0.1);
	if (fit->type == DATA_USHORT) {
		WORD *data = fit->data;
		for (row = 0; row < fit->ry; row ++) {
			for (col = 0; col < fit->rx; col++) {
				pat_cell = (row % pat_width) * pat_width + col % pat_width;
				switch (cfa_string[pat_cell]) {
					case 'R':
						data[col + row * fit->rx] = round_to_WORD(data[col + row * fit->rx] / coeff1);
						break;
					case 'B':
						data[col + row * fit->rx] = round_to_WORD(data[col + row * fit->rx] / coeff2);
						break;
				}
			}
		}
	}
	else if (fit->type == DATA_FLOAT) {
		float *data = fit->fdata;
		for (row = 0; row < fit->ry; row ++) {
			for (col = 0; col < fit->rx; col++) {
				pat_cell = (row % pat_width) * pat_width + col % pat_width;
				switch (cfa_string[pat_cell]) {
					case 'R':
						data[col + row * fit->rx] /= coeff1;
						break;
					case 'B':
						data[col + row * fit->rx] /= coeff2;
						break;
				}
			}
		}
	}
	else return 1;
	return 0;
}

static gpointer extract_channels_ushort(gpointer p) {
	struct extract_channels_data *args = (struct extract_channels_data *) p;
	WORD *buf[3] = { args->fit->pdata[RLAYER], args->fit->pdata[GLAYER],
		args->fit->pdata[BLAYER] };
	size_t n = args->fit->naxes[0] * args->fit->naxes[1];
	struct timeval t_start, t_end;

	if (args->fit->naxes[2] != 3) {
		siril_log_message(
				_("Siril cannot extract layers. Make sure your image is in RGB mode.\n"));
		return GINT_TO_POINTER(1);
	}

	siril_log_color_message(_("%s channel extraction: processing...\n"), "green",
			args->str_type);
	gettimeofday(&t_start, NULL);
	gchar *histstring = NULL;
	cmsHPROFILE cielab_profile = NULL, image_profile = NULL;
	cmsColorSpaceSignature sig;
	cmsUInt32Number trans_type, lab_type;
	gboolean threaded;
	cmsHTRANSFORM transform = NULL;
	cmsUInt32Number datasize;
	cmsUInt32Number bytesperline;
	cmsUInt32Number bytesperplane;
	gchar *desc = siril_color_profile_get_description(args->fit->icc_profile);
	if(args->fit->icc_profile)
		cmsCloseProfile(args->fit->icc_profile);
	/* The extracted channels are considered raw data, and are not color
		* managed. It is up to the user to ensure that future use of them is
		* with similar data and an appropriate color profile is assigned.
		* See also the HSV and CIELAB cases below.*/
	args->fit->icc_profile = NULL;
	color_manage(args->fit, FALSE);

	switch (args->type) {
	case EXTRACT_RGB:
		histstring = g_strdup_printf(_("%s: extract RGB channel"), extractionstring);
		break;
	case EXTRACT_HSL:
		histstring = g_strdup_printf(_("%s: extract HSL channel"), extractionstring);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < n; i++) {
			double h, s, l;
			double r = (double) buf[RLAYER][i] / USHRT_MAX_DOUBLE;
			double g = (double) buf[GLAYER][i] / USHRT_MAX_DOUBLE;
			double b = (double) buf[BLAYER][i] / USHRT_MAX_DOUBLE;
			// RGB to HSL is a coordinate transform and does not require lcms
			rgb_to_hsl(r, g, b, &h, &s, &l);
			buf[RLAYER][i] = round_to_WORD(h * 360.0);
			buf[GLAYER][i] = round_to_WORD(s * USHRT_MAX_DOUBLE);
			buf[BLAYER][i] = round_to_WORD(l * USHRT_MAX_DOUBLE);
		}
		break;
	case EXTRACT_HSV:
		histstring = g_strdup_printf(_("%s: extract HSV channel"), extractionstring);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < n; i++) {
			double h, s, v;
			double r = (double) buf[RLAYER][i] / USHRT_MAX_DOUBLE;
			double g = (double) buf[GLAYER][i] / USHRT_MAX_DOUBLE;
			double b = (double) buf[BLAYER][i] / USHRT_MAX_DOUBLE;
			// RGB to HSV is a coordinate transform and does not require lcms
			rgb_to_hsv(r, g, b, &h, &s, &v);
			buf[RLAYER][i] = round_to_WORD(h * 360.0);
			buf[GLAYER][i] = round_to_WORD(s * USHRT_MAX_DOUBLE);
			buf[BLAYER][i] = round_to_WORD(v * USHRT_MAX_DOUBLE);
		}
		break;
	case EXTRACT_CIELAB:
		histstring = g_strdup_printf(_("%s: extract LAB channel"), extractionstring);
		cielab_profile = cmsCreateLab4Profile(NULL);
		if (args->fit->icc_profile) {
			image_profile = copyICCProfile(args->fit->icc_profile);
		} else {
			siril_log_message(_("Image is not color managed. Assuming sRGB.\n"));
			image_profile = srgb_trc();
		}
		sig = cmsGetColorSpace(image_profile);
		trans_type = get_planar_formatter_type(sig, args->fit->type, FALSE);
		lab_type = TYPE_Lab_16_PLANAR;
		threaded = !get_thread_run();
		// We use sRGB as the fallback for non-color managed images
		transform = cmsCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), image_profile, trans_type, cielab_profile, lab_type, INTENT_PERCEPTUAL, com.icc.rendering_flags);
		cmsCloseProfile(cielab_profile);
		cmsCloseProfile(image_profile);
		datasize = sizeof(WORD);
		bytesperline = args->fit->rx * datasize;
		bytesperplane = args->fit->rx * args->fit->ry * datasize;
		cmsDoTransformLineStride(transform, args->fit->data, args->fit->data, args->fit->rx, args->fit->ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
		cmsDeleteTransform(transform);
	}
	gchar *fitfilter = g_strdup(args->fit->filter);
	if (desc) {
		args->fit->history = g_slist_append(args->fit->history, g_strdup_printf(_("Channel extraction from 3-channel image with ICC profile:")));
		args->fit->history = g_slist_append(args->fit->history, g_strdup_printf("%s", desc));
	}
	for (int i = 0; i < 3; i++) {
		if (args->channel[i]) {
			update_filter_information(args->fit, add_filter_str[i], TRUE);
			if (i > 0) {
				GSList *current = args->fit->history;
				while (current->next != NULL && current->next->next != NULL) {
					current = current->next;
				}
				// Check if there is only one element in the list.
				if (current->next == NULL) {
					g_slist_free_full(args->fit->history, g_free);
					args->fit->history = NULL;
				} else {
					// Remove the last element.
					GSList *last = current->next;
					current->next = NULL;
					g_free(last->data);
					g_slist_free_1(last);
				}
			}
			args->fit->history = g_slist_append(args->fit->history, g_strdup_printf("%s %d", histstring, i));
			save1fits16(args->channel[i], args->fit, i);
			update_filter_information(args->fit, fitfilter, FALSE); //reinstate original filter name
		}
	}
	g_free(fitfilter);
	g_free(histstring);
	g_free(desc);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	return GINT_TO_POINTER(0);
}

static gpointer extract_channels_float(gpointer p) {
	struct extract_channels_data *args = (struct extract_channels_data *) p;
	float *buf[3] = { args->fit->fpdata[RLAYER], args->fit->fpdata[GLAYER],
		args->fit->fpdata[BLAYER] };
	struct timeval t_start, t_end;
	size_t n = args->fit->naxes[0] * args->fit->naxes[1];

	if (args->fit->naxes[2] != 3) {
		siril_log_message(
				_("Siril cannot extract channels. Make sure your image is in RGB mode.\n"));
		return GINT_TO_POINTER(1);
	}

	siril_log_color_message(_("%s channel extraction: processing...\n"), "green",
			args->str_type);
	gettimeofday(&t_start, NULL);
	gchar *histstring = NULL;
	cmsHPROFILE cielab_profile = NULL, image_profile = NULL;
	cmsColorSpaceSignature sig;
	cmsUInt32Number trans_type, lab_type;
	gboolean threaded;
	cmsHTRANSFORM transform = NULL;
	cmsUInt32Number datasize;
	cmsUInt32Number bytesperline;
	cmsUInt32Number bytesperplane;
	gchar *desc = NULL;
	if(args->fit->icc_profile) {
		desc = siril_color_profile_get_description(args->fit->icc_profile);
		cmsCloseProfile(args->fit->icc_profile);
	}
	/* The extracted channels are considered raw data, and are not color
		* managed. It is up to the user to ensure that future use of them is
		* with similar data and an appropriate color profile is assigned.
		* See also the HSV and CIELAB cases below.*/
	args->fit->icc_profile = NULL;
	color_manage(args->fit, FALSE);

	switch (args->type) {
	case EXTRACT_RGB:
		histstring = g_strdup_printf(_("%s: extract RGB channel"), extractionstring);
		break;
	case EXTRACT_HSL:
		histstring = g_strdup_printf(_("%s: extract HSL channel"), extractionstring);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < n; i++) {
			double h, s, l;
			double r = (double) buf[RLAYER][i];
			double g = (double) buf[GLAYER][i];
			double b = (double) buf[BLAYER][i];
			rgb_to_hsl(r, g, b, &h, &s, &l);
			buf[RLAYER][i] = (float) h;
			buf[GLAYER][i] = (float) s;
			buf[BLAYER][i] = (float) l;
		}
		break;
	case EXTRACT_HSV:
		histstring = g_strdup_printf(_("%s: extract HSV channel"), extractionstring);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < n; i++) {
			double h, s, v;
			double r = (double) buf[RLAYER][i];
			double g = (double) buf[GLAYER][i];
			double b = (double) buf[BLAYER][i];
			rgb_to_hsv(r, g, b, &h, &s, &v);
			buf[RLAYER][i] = (float) h;
			buf[GLAYER][i] = (float) s;
			buf[BLAYER][i] = (float) v;
		}
		break;
	case EXTRACT_CIELAB:
		histstring = g_strdup_printf(_("%s: extract LAB channel"), extractionstring);
		cielab_profile = cmsCreateLab4Profile(NULL);
		if (args->fit->icc_profile) {
			image_profile = copyICCProfile(args->fit->icc_profile);
		} else {
			siril_log_message(_("Image is not color managed. Assuming sRGB.\n"));
			image_profile = srgb_trc();
		}
		sig = cmsGetColorSpace(image_profile);
		trans_type = get_planar_formatter_type(sig, args->fit->type, FALSE);
		lab_type = TYPE_Lab_FLT_PLANAR;
		threaded = !get_thread_run();
		transform = cmsCreateTransformTHR((threaded ? com.icc.context_threaded : com.icc.context_single), image_profile, trans_type, cielab_profile, lab_type, com.pref.icc.processing_intent, com.icc.rendering_flags);
		cmsCloseProfile(cielab_profile);
		cmsCloseProfile(image_profile);
		datasize = sizeof(float);
		bytesperline = args->fit->rx * datasize;
		bytesperplane = args->fit->rx * args->fit->ry * datasize;
		/* Note this output is in CIE La*b* ranges (ie L [0..100] etc, not Siril's
		 * usual [0..1] range.
		 * TODO: convert to Siril ranges?
		 */
		cmsDoTransformLineStride(transform, args->fit->fdata, args->fit->fdata, args->fit->rx, args->fit->ry, bytesperline, bytesperline, bytesperplane, bytesperplane);
		cmsDeleteTransform(transform);
	}
	gchar *fitfilter = g_strdup(args->fit->filter);
	if (desc) {
		args->fit->history = g_slist_append(args->fit->history, g_strdup_printf(_("Channel extraction from 3-channel image with ICC profile:")));
		args->fit->history = g_slist_append(args->fit->history, g_strdup_printf("%s", desc));
	}
	for (int i = 0; i < 3; i++) {
		if (args->channel[i]) {
			update_filter_information(args->fit, add_filter_str[i], TRUE);
			if (i > 0) {
				GSList *current = args->fit->history;
				while (current->next != NULL && current->next->next != NULL) {
					current = current->next;
				}
				// Check if there is only one element in the list.
				if (current->next == NULL) {
					g_slist_free_full(args->fit->history, g_free);
					args->fit->history = NULL;
				} else {
					// Remove the last element.
					GSList *last = current->next;
					current->next = NULL;
					g_free(last->data);
					g_slist_free_1(last);
				}
			}
			args->fit->history = g_slist_append(args->fit->history, g_strdup_printf("%s %d", histstring, i));
			save1fits32(args->channel[i], args->fit, i);
			update_filter_information(args->fit, fitfilter, FALSE); //reinstate original filter name
		}
	}
	g_free(desc);
	g_free(histstring);
	g_free(fitfilter);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	return GINT_TO_POINTER(0);
}

// channels extraction computed in-place
gpointer extract_channels(gpointer p) {
	struct extract_channels_data *args = (struct extract_channels_data *)p;
	gpointer retval = GINT_TO_POINTER(1);

	if (args->fit->type == DATA_USHORT)
		retval = extract_channels_ushort(p);
	else if (args->fit->type == DATA_FLOAT)
		retval = extract_channels_float(p);

	clearfits(args->fit);
	free(args->channel[0]);
	free(args->channel[1]);
	free(args->channel[2]);
	free(args);
	siril_add_idle(end_generic, NULL);
	return retval;
}

/****************** Color calibration ************************/

/* This function equalize the background by giving equal value for all layers */
void background_neutralize(fits* fit, rectangle black_selection) {
	int chan;
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	imstats* stats[3];
	double ref = 0;

	assert(fit->naxes[2] == 3);

	for (chan = 0; chan < 3; chan++) {
		stats[chan] = statistics(NULL, -1, fit, chan, &black_selection, STATS_BASIC, MULTI_THREADED);
		if (!stats[chan]) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return;
		}
		ref += stats[chan]->median;
	}
	ref /= 3.0;

	if (fit->type == DATA_USHORT) {
		for (chan = 0; chan < 3; chan++) {
			double offset = stats[chan]->mean - ref;
			WORD *buf = fit->pdata[chan];
			for (i = 0; i < n; i++) {
				buf[i] = round_to_WORD((double)buf[i] - offset);
			}
			free_stats(stats[chan]);
		}
	}
	else if (fit->type == DATA_FLOAT) {
		for (chan = 0; chan < 3; chan++) {
			float offset = stats[chan]->mean - ref;
			float *buf = fit->fpdata[chan];
			for (i = 0; i < n; i++) {
				buf[i] = buf[i] - offset;
			}
			free_stats(stats[chan]);
		}
	}

	invalidate_stats_from_fit(fit);
	invalidate_gfit_histogram();
}

void get_coeff_for_wb(fits *fit, rectangle white, rectangle black,
		double kw[], double bg[], double norm, double low, double high) {
	int chan, i, j, n;
	double tmp[3] = { 0.0, 0.0, 0.0 };

	assert(fit->naxes[2] == 3);

	if (fit->type == DATA_USHORT) {
		WORD lo = round_to_WORD(low * (norm));
		WORD hi = round_to_WORD(high * (norm));

		for (chan = 0; chan < 3; chan++) {
			n = 0;
			WORD *from = fit->pdata[chan] + (fit->ry - white.y - white.h) * fit->rx
				+ white.x;
			int stridefrom = fit->rx - white.w;

			for (i = 0; i < white.h; i++) {
				for (j = 0; j < white.w; j++) {
					if (*from > lo && *from < hi ) {
						kw[chan] += (double)*from / norm;
						n++;
					}
					from++;
				}
				from += stridefrom;
			}
			if (n > 0)
				kw[chan] /= (double)n;
		}
	}
	else if (fit->type == DATA_FLOAT) {
		for (chan = 0; chan < 3; chan++) {
			n = 0;
			float *from = fit->fpdata[chan] + (fit->ry - white.y - white.h) * fit->rx
				+ white.x;
			int stridefrom = fit->rx - white.w;

			for (i = 0; i < white.h; i++) {
				for (j = 0; j < white.w; j++) {
					double f = (double)*from;
					if (f > low && f < high) {
						kw[chan] += f;
						n++;
					}
					from++;
				}
				from += stridefrom;
			}
			if (n > 0)
				kw[chan] /= (double)n;
		}
	}
	else return;

	siril_log_message(_("Background reference:\n"));
	for (chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, &black, STATS_BASIC, MULTI_THREADED);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return;
		}
		bg[chan] = stat->median / stat->normValue;
		siril_log_message("B%d: %.5e\n", chan, bg[chan]);
		free_stats(stat);

	}

	siril_log_message(_("White reference:\n"));
	for (chan = 0; chan < 3; chan++) {
		siril_log_message("W%d: %.5e\n", chan, kw[chan]);
		kw[chan] = fabs(kw[chan] - bg[chan]);
	}

	int rc = (kw[0] > kw[1]) ? ((kw[0] > kw[2]) ? 0 : 2) :
		((kw[1] > kw[2]) ? 1 : 2);
	for (chan = 0; chan < 3; chan++) {
		if (chan == rc)
			tmp[chan] = 1.0;
		else
			tmp[chan] = kw[rc] / kw[chan];
	}

	siril_log_message(_("Color calibration factors:\n"));
	for (chan = 0; chan < 3; chan++) {
		kw[chan] = tmp[chan];
		siril_log_message("K%d: %5.3lf\n", chan, kw[chan]);
	}
}

int calibrate(fits *fit, int layer, double kw, double bg, double norm) {
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->type == DATA_USHORT) {
		double bgNorm = bg * norm;
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < n; ++i) {
			buf[i] = round_to_WORD((buf[i] - bgNorm) * kw + bgNorm);
		}
	}
	else if (fit->type == DATA_FLOAT) {
		float *buf = fit->fpdata[layer];
		for (i = 0; i < n; ++i) {
			buf[i] = (float)(((double)buf[i] - bg) * kw + bg);
		}
	}
	else return 1;
	return 0;
}

int pos_to_neg(fits *fit) {
	size_t i, n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
	if (fit->type == DATA_USHORT) {
		WORD norm = (WORD)get_normalized_value(fit);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
		for (i = 0; i < n; i++) {
			fit->data[i] = norm - fit->data[i];
		}
	}
	else if (fit->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
		for (i = 0; i < n; i++) {
			fit->fdata[i] = 1.0f - fit->fdata[i];
		}
	}
	else return 1;

	return 0;
}
