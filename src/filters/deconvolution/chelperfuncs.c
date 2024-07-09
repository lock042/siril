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

#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "chelperfuncs.h"
#include "gui/progress_and_log.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "core/siril_log.h"

#ifndef RT_INCLUDE
#undef max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#undef min
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

// Helpful wrapper functions to call Siril functions without causing clashes with some of the img_expr functionality

void updateprogress(const char *text, double percent) {
	set_progress_bar_data(text, percent);
}

int is_thread_stopped() {
	return (get_thread_run() ? 0 : 1);
}

void sirillog(const char* text) {
	siril_log_message(_(text));
}

int updatenoise(float *array, int nx, int ny, int nchans, double *noise) {
	return sos_update_noise_float(array, (long) nx, (long) ny, (long) nchans, noise);
}
