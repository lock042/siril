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

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/optimize_utils.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "algos/statistics.h"
#include "algos/sorting.h"

#include "median.h"
#include "algos/median_fast.h"

gchar* median_log_hook(gpointer p, log_hook_detail detail) {
	struct median_filter_data *params = (struct median_filter_data *) p;
	if (detail == SUMMARY)
		return g_strdup_printf(_("Median filter: ksize %d, iters %d, amount %.3f"),
				params->ksize, params->iterations, params->amount);
	return g_strdup_printf(_("Median filter: kernel size %d, iters %d, amount %.3f"),
			params->ksize, params->iterations, params->amount);
}

/*****************************************************************************
 *                      M E D I A N     I M A G E     F I L T E R S          *
 ****************************************************************************/

double get_median_ushort(const WORD *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	WORD *values;
	double median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(WORD));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				if (include_self || x != xx || y != yy) {
					values[n++] = buf[x + y * w];
				}
			}
		}
	}
	median = quickmedian(values, n);
	free(values);
	return median;
}

double get_median_float(const float *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	float *values;
	double median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(float));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				if (include_self || x != xx || y != yy) {
					values[n++] = buf[x + y * w];
				}
			}
		}
	}
	median = quickmedian_float(values, n);
	free(values);
	return median;
}

static float get_median_float_fast(const float *buf, const int xx, const int yy, const int w,
		const int h, int radius) {

	int ksize = radius * 2 + 1;
	float values[ksize * ksize];

	int ystart = (yy - radius) < 0 ? 0 : yy - radius;
	int yend = (yy + radius) >= h ? h - 1 : yy + radius;
	int xstart = (xx - radius) < 0 ? 0 : xx - radius;
	int xend = (xx + radius) >= w ? w - 1 : xx + radius;
	int n = 0;
	for (int y = ystart; y <= yend; ++y) {
		for (int x = xstart; x <= xend; ++x) {
			values[n++] = buf[x + y * w];
		}
	}
	return quickmedian_float(values, n);
}

static float get_median_ushort_fast(const WORD *buf, const int xx, const int yy, const int w,
		const int h, int radius) {

	int ksize = radius * 2 + 1;
	WORD values[ksize * ksize];

	int ystart = (yy - radius) < 0 ? 0 : yy - radius;
	int yend = (yy + radius) >= h ? h - 1 : yy + radius;
	int xstart = (xx - radius) < 0 ? 0 : xx - radius;
	int xend = (xx + radius) >= w ? w - 1 : xx + radius;
	int n = 0;
	for (int y = ystart; y <= yend; ++y) {
		for (int x = xstart; x <= xend; ++x) {
			values[n++] = buf[x + y * w];
		}
	}
	return quickmedian(values, n);
}

double get_median_gsl(gsl_matrix *mat, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	double *values, median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(double));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				if (include_self || x != xx || y != yy) {
					values[n++] = gsl_matrix_get(mat, y, x);
				}
			}
		}
	}
	median = quickmedian_double(values, n);
	free(values);
	return median;
}


/*****************************************************************************
 *                      M E D I A N     F I L T E R                          *
 ****************************************************************************/

static gpointer median_filter_ushort(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *)p;
	int progress = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	double total;
	int radius = (args->ksize - 1) / 2;

	g_assert(args->ksize % 2 == 1 && args->ksize > 1);
	g_assert(nx > 0 && ny > 0);
	total = ny * args->fit->naxes[2] * args->iterations;

	size_t alloc_size = args->fit->naxes[0] * args->fit->naxes[1] * sizeof(WORD);
	WORD *temp = calloc(1, alloc_size);
	if (!temp) {
		PRINT_ALLOC_ERR;
		return GINT_TO_POINTER(-1);
	}
	float amountf = args->amount;
	for (int layer = 0; layer < args->fit->naxes[2]; layer++) {
		for (int iter = 0; iter < args->iterations; ++iter) {
			WORD *dst = (iter % 2) ? args->fit->pdata[layer] : temp;
			WORD *src = (iter % 2) ? temp : args->fit->pdata[layer];
			for (int y = 0; y < ny; y++) {
				if (y < radius || y >= ny - radius) {
					for (int x = 0; x < nx; x++) {
						if (x < radius || x >= nx - radius) {
							int pix_idx = y * nx + x;
							float median = get_median_ushort_fast(src, x, y, nx, ny, radius);
							dst[pix_idx] = roundf_to_WORD(intpf(amountf, median, (float) src[pix_idx]));
						}
					}
				}
			}
			if (args->ksize == 3) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread)
#endif
				for (int y = 1; y < ny - 1; y++) {
#ifdef _OPENMP
#pragma omp simd
#endif
					for (int x = 1; x < nx - 1; x++) {
						float median = median9f((float)src[(y - 1) * nx + x - 1],
								(float)src[(y - 1) * nx + x],
								(float)src[(y - 1) * nx + x + 1],
								(float)src[y * nx + x - 1],
								(float)src[y * nx + x],
								(float)src[y * nx + x + 1],
								(float)src[(y + 1) * nx + x - 1],
								(float)src[(y + 1) * nx + x],
								(float)src[(y + 1) * nx + x + 1]);
						dst[y * nx + x] = roundf_to_WORD(intpf(amountf, median, (float) src[y * nx + x]));
					}
#ifdef _OPENMP
#pragma omp atomic
#endif
					++progress;
					if (!(progress % 32))
						gui_iface.set_progress((double)progress / total, NULL);
				}
			} else if (args->ksize == 5) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[25];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
					for (int y = 2; y < ny - 2; y++) {
						for (int x = 2; x < nx - 2; x++) {
							for (int i = -2; i <= 2; ++i)
								for (int j = -2; j <= 2; ++j)
									medbuf[(i + 2) * 5 + (j + 2)] = (float) src[(y + i) * nx + x + j];
							float median = median5x5(medbuf);
							dst[y * nx + x] = roundf_to_WORD(intpf(amountf, median, (float) src[y * nx + x]));
						}
#ifdef _OPENMP
#pragma omp atomic
#endif
						++progress;
						if (!(progress % 32))
							gui_iface.set_progress((double)progress / total, NULL);
					}
				}
			} else if (args->ksize == 7) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[49];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
					for (int y = 3; y < ny - 3; y++) {
						for (int x = 3; x < nx - 3; x++) {
							for (int i = -3; i <= 3; ++i)
								for (int j = -3; j <= 3; ++j)
									medbuf[(i + 3) * 7 + (j + 3)] = (float) src[(y + i) * nx + x + j];
							float median = median7x7(medbuf);
							dst[y * nx + x] = roundf_to_WORD(intpf(amountf, median, (float) src[y * nx + x]));
						}
#ifdef _OPENMP
#pragma omp atomic
#endif
						++progress;
						if (!(progress % 32))
							gui_iface.set_progress((double)progress / total, NULL);
					}
				}
			} else if (args->ksize == 9) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[81];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
					for (int y = 4; y < ny - 4; y++) {
						for (int x = 4; x < nx - 4; x++) {
							for (int i = -4; i <= 4; ++i)
								for (int j = -4; j <= 4; ++j)
									medbuf[(i + 4) * 9 + (j + 4)] = (float) src[(y + i) * nx + x + j];
							float median = median9x9(medbuf);
							dst[y * nx + x] = roundf_to_WORD(intpf(amountf, median, (float) src[y * nx + x]));
						}
#ifdef _OPENMP
#pragma omp atomic
#endif
						++progress;
						if (!(progress % 32))
							gui_iface.set_progress((double)progress / total, NULL);
					}
				}
			} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
				for (int y = 0; y < ny; y++) {
					for (int x = 0; x < nx; x++) {
						float median = get_median_ushort_fast(src, x, y, nx, ny, radius);
						dst[y * nx + x] = roundf_to_WORD(intpf(amountf, median, (float) src[y * nx + x]));
					}
#ifdef _OPENMP
#pragma omp atomic
#endif
					++progress;
					if (!(progress % 32))
						gui_iface.set_progress((double)progress / total, NULL);
				}
			}
		}
		if (args->iterations % 2) {
			WORD *dst = args->fit->pdata[layer];
			WORD *src = temp;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
			for (int y = 0; y < ny; y++) {
				for (int x = 0; x < nx; x++) {
					dst[y * nx + x] = src[y * nx + x];
				}
			}
		}
	}
	free(temp);
	invalidate_stats_from_fit(args->fit);

	return GINT_TO_POINTER(0);
}

static gpointer median_filter_float(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *)p;
	int progress = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	double total;
	int radius = (args->ksize - 1) / 2;

	g_assert(args->ksize % 2 == 1 && args->ksize > 1);
	g_assert(nx > 0 && ny > 0);
	total = ny * args->fit->naxes[2] * args->iterations;

	size_t alloc_size = args->fit->naxes[0] * args->fit->naxes[1] * sizeof(float);
	float *temp = calloc(1, alloc_size);
	if (!temp) {
		PRINT_ALLOC_ERR;
		return GINT_TO_POINTER(-1);
	}
	float amountf = args->amount;
	for (int layer = 0; layer < args->fit->naxes[2]; layer++) {
		for (int iter = 0; iter < args->iterations; ++iter) {
			float *dst = (iter % 2) ? args->fit->fpdata[layer] : temp;
			float *src = (iter % 2) ? temp : args->fit->fpdata[layer];
			for (int y = 0; y < ny; y++) {
				if (y < radius || y >= ny - radius) {
					for (int x = 0; x < nx; x++) {
						if (x < radius || x >= nx - radius) {
							int pix_idx = y * nx + x;
							float median = get_median_float_fast(src, x, y, nx, ny, radius);
							dst[pix_idx] = intpf(amountf, median, src[pix_idx]);
						}
					}
				}
			}
			if (args->ksize == 3) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(com.max_thread)
#endif
				for (int y = 1; y < ny - 1; y++) {
#ifdef _OPENMP
#pragma omp simd
#endif
					for (int x = 1; x < nx - 1; x++) {
						float median = median9f(src[(y - 1) * nx + x - 1],
								src[(y - 1) * nx + x],
								src[(y - 1) * nx + x + 1],
								src[y * nx + x - 1],
								src[y * nx + x],
								src[y * nx + x + 1],
								src[(y + 1) * nx + x - 1],
								src[(y + 1) * nx + x],
								src[(y + 1) * nx + x + 1]);
						dst[y * nx + x] = intpf(amountf, median, src[y * nx + x]);
					}
#ifdef _OPENMP
#pragma omp atomic
#endif
					++progress;
					if (!(progress % 32))
						gui_iface.set_progress((double)progress / total, NULL);
				}
			} else if (args->ksize == 5) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[25];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
					for (int y = 2; y < ny - 2; y++) {
						for (int x = 2; x < nx - 2; x++) {
							for (int i = -2; i <= 2; ++i)
								for (int j = -2; j <= 2; ++j)
									medbuf[(i + 2) * 5 + (j + 2)] = src[(y + i) * nx + x + j];
							float median = median5x5(medbuf);
							dst[y * nx + x] = intpf(amountf, median, src[y * nx + x]);
						}
#ifdef _OPENMP
#pragma omp atomic
#endif
						++progress;
						if (!(progress % 32))
							gui_iface.set_progress((double)progress / total, NULL);
					}
				}
			} else if (args->ksize == 7) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[49];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
					for (int y = 3; y < ny - 3; y++) {
						for (int x = 3; x < nx - 3; x++) {
							for (int i = -3; i <= 3; ++i)
								for (int j = -3; j <= 3; ++j)
									medbuf[(i + 3) * 7 + (j + 3)] = src[(y + i) * nx + x + j];
							float median = median7x7(medbuf);
							dst[y * nx + x] = intpf(amountf, median, src[y * nx + x]);
						}
#ifdef _OPENMP
#pragma omp atomic
#endif
						++progress;
						if (!(progress % 32))
							gui_iface.set_progress((double)progress / total, NULL);
					}
				}
			} else if (args->ksize == 9) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[81];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
					for (int y = 4; y < ny - 4; y++) {
						for (int x = 4; x < nx - 4; x++) {
							for (int i = -4; i <= 4; ++i)
								for (int j = -4; j <= 4; ++j)
									medbuf[(i + 4) * 9 + (j + 4)] = src[(y + i) * nx + x + j];
							float median = median9x9(medbuf);
							dst[y * nx + x] = intpf(amountf, median, src[y * nx + x]);
						}
#ifdef _OPENMP
#pragma omp atomic
#endif
						++progress;
						if (!(progress % 32))
							gui_iface.set_progress((double)progress / total, NULL);
					}
				}
			} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
				for (int y = 0; y < ny; y++) {
					for (int x = 0; x < nx; x++) {
						float median = get_median_float_fast(src, x, y, nx, ny, radius);
						dst[y * nx + x] = intpf(amountf, median, src[y * nx + x]);
					}
#ifdef _OPENMP
#pragma omp atomic
#endif
					++progress;
					if (!(progress % 32))
						gui_iface.set_progress((double)progress / total, NULL);
				}
			}
		}
		if (args->iterations % 2) {
			float *dst = args->fit->fpdata[layer];
			float *src = temp;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
			for (int y = 0; y < ny; y++) {
				for (int x = 0; x < nx; x++) {
					dst[y * nx + x] = src[y * nx + x];
				}
			}
		}
	}
	free(temp);
	invalidate_stats_from_fit(args->fit);

	return GINT_TO_POINTER(0);
}

int median_image_hook(struct generic_img_args *args, fits *fit, int nb_threads) {
	struct median_filter_data *params = (struct median_filter_data*) args->user;
	if (!params)
		return 1;
	params->fit = fit;
	gpointer result;
	if (fit->type == DATA_USHORT)
		result = median_filter_ushort(params);
	else if (fit->type == DATA_FLOAT)
		result = median_filter_float(params);
	else
		return 1;
	return GPOINTER_TO_INT(result);
}

/* GUI callbacks moved to src/gui/median.c */
