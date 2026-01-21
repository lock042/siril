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

#ifndef SIRIL_CURVE_TRANSFORM_H
#define SIRIL_CURVE_TRANSFORM_H

#include "core/siril.h"
#include <glib.h>

#define MAX_POINTS 100
#define LUT_SIZE 65536

enum curve_algorithm {
	CUBIC_SPLINE,
	LINEAR,
	AKIMA_SPLINE
};

enum curve_channel {
	CHAN_RGB_K,
	CHAN_R, CHAN_G, CHAN_B,
	CHAN_L, CHAN_C, CHAN_S,
	CHAN_COUNT
};

// Configuration for a single channel (Points + Range Masking)
typedef struct {
	GList *points;              // List of point* (x,y)
	gboolean active;            // Is this channel effectively modified?
	
	// Luminance Range Masking
	gboolean range_enabled;
	float lum_min;              // 0.0 - 1.0
	float lum_max;              // 0.0 - 1.0
	float feather;              // 0.0 - 1.0 (Sigma equivalent)
} curve_channel_config;

struct curve_params {
	destructor destroy_fn;
	
	// Array of configurations for ALL channels
	curve_channel_config channels[CHAN_COUNT];
	
	enum curve_algorithm algorithm;
	enum curve_channel target_channel; 
	
	fits *fit;
	gboolean verbose;
	gboolean for_preview;

	long *clipped_count;
};

// Math structs (Internal use for fitting)
typedef struct {
	double a[MAX_POINTS], b[MAX_POINTS], c[MAX_POINTS], d[MAX_POINTS];
	double x_values[MAX_POINTS], y_values[MAX_POINTS];
	int n;
} cubic_spline_data;

typedef struct {
	double b[MAX_POINTS], c[MAX_POINTS], d[MAX_POINTS];
	double x_values[MAX_POINTS], y_values[MAX_POINTS];
	int n;
} akima_spline_data;

// Runtime Cache (LUT + Mask buffer)
typedef struct {
	gboolean active;
	float lut[LUT_SIZE]; 
	gboolean use_mask;
	float *mask_buffer;  
} channel_runtime_cache;

struct curve_params *new_curve_params();
void free_curve_params(void *args);

void apply_curve(fits *from, fits *to, struct curve_params *params, gboolean multithreaded);
int curve_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

gchar *curves_log_hook(gpointer p, log_hook_detail detail);
gboolean curve_preview_idle(gpointer p);
gboolean curve_apply_idle(gpointer p);

// EXPOSED MATH FUNCTIONS (Crucial for GUI)
void linear_fit(GList *points, double *slopes);
float linear_interpolate(float x, GList *points, double *slopes);

void cubic_spline_fit(GList *points, cubic_spline_data *cspline_data);
float cubic_spline_interpolate(float x, cubic_spline_data *cspline_data);

void akima_spline_fit(GList *points, akima_spline_data *akima_data);
float akima_spline_interpolate(float x, akima_spline_data *akima_data);

#endif