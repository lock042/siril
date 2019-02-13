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

/* This file provides ways to cache the operations used in multi-point
 * planetary algorithms: conversion to monochrome, gaussian filtering and
 * laplacian operator.
 * Data is computed from input images, can be saved in the SER cache or
 * obtained instead from the cache.
 */

#include "planetary_caching.h"
#include "io/sequence.h"
#include "core/proto.h"
#include "opencv/opencv.h"
#include "io/ser.h"

static void disable_cache_reading(struct planetary_cache *args);
static void disable_cache_writing(struct planetary_cache *args);

void init_caching(const char *seqname, struct planetary_cache *args, int kernel_size) {
	args->kernel_size = kernel_size;

	if (args->use_cached) {
		char *gaussian_seqname = malloc(strlen(seqname) + 20);
		sprintf(gaussian_seqname, "%s-gaussian%d.seq", seqname, args->kernel_size);

		if (!(args->seq_gaussian = readseqfile(gaussian_seqname)) ||
				seq_check_basic_data(args->seq_gaussian, FALSE) == -1) {
			fprintf(stderr, "failed not load sequence %s, not using cached data\n", gaussian_seqname);
			disable_cache_reading(args);
			free(gaussian_seqname);
			goto cache_writing;
		}
		siril_log_message("using cache (%s) to get gaussian filtered images\n", gaussian_seqname);
		free(gaussian_seqname);

		/* same as above for the laplacian */
		char *laplacian_seqname = malloc(strlen(seqname) + 20);
		sprintf(laplacian_seqname, "%s-laplacian%d.seq", seqname, args->kernel_size);

		if (!(args->seq_laplacian = readseqfile(laplacian_seqname)) ||
				seq_check_basic_data(args->seq_laplacian, FALSE) == -1) {
			fprintf(stderr, "failed not load sequence %s, not using cached data\n", laplacian_seqname);
			disable_cache_reading(args);
			free(laplacian_seqname);
			goto cache_writing;
		}
		siril_log_message("using cache (%s) to get laplacian filtered images\n", laplacian_seqname);
		free(laplacian_seqname);
	}

cache_writing:
	if (args->cache_data) {
		char *gaussian_seqname = malloc(strlen(seqname) + 20);
		sprintf(gaussian_seqname, "%s-gaussian%d.ser", seqname, args->kernel_size);

		// basically do the same as ser_prepare_hook
		args->ser_gaussian = malloc(sizeof(struct ser_struct));
		if (ser_create_file(gaussian_seqname, args->ser_gaussian, TRUE, NULL)) {
			free(gaussian_seqname);
			disable_cache_writing(args);
			return;
		}
		siril_log_message("saving gaussian filtered images to cache (%s)\n", gaussian_seqname);
		free(gaussian_seqname);

		/* same as above for the laplacian */
		char *laplacian_seqname = malloc(strlen(seqname) + 20);
		sprintf(laplacian_seqname, "%s-laplacian%d.ser", seqname, args->kernel_size);

		args->ser_laplacian = malloc(sizeof(struct ser_struct));
		if (ser_create_file(laplacian_seqname, args->ser_laplacian, TRUE, NULL)) {
			free(laplacian_seqname);
			disable_cache_writing(args);
			return;
		}
		siril_log_message("saving laplacian filtered images to cache (%s)\n", laplacian_seqname);
		free(laplacian_seqname);
	}
}

static void finalize_cache_reading(struct planetary_cache *args) {
	if (args->seq_gaussian) {
		free_sequence(args->seq_gaussian, TRUE);
	}
	args->seq_gaussian = NULL;

	if (args->seq_laplacian) {
		free_sequence(args->seq_laplacian, TRUE);
	}
	args->seq_laplacian = NULL;
}

static void finalize_cache_writing(struct planetary_cache *args) {
	gboolean had_some_seq = FALSE;
	if (args->ser_gaussian) {
		ser_write_and_close(args->ser_gaussian);
		free(args->ser_gaussian);
		had_some_seq = TRUE;
	}
	args->ser_gaussian = NULL;

	if (args->ser_laplacian) {
		ser_write_and_close(args->ser_laplacian);
		free(args->ser_laplacian);
		had_some_seq = TRUE;
	}
	args->ser_laplacian = NULL;

	if (had_some_seq)
		check_seq(FALSE);
}

void finalize_caching(struct planetary_cache *args) {
	finalize_cache_reading(args);
	finalize_cache_writing(args);
}

static void disable_cache_reading(struct planetary_cache *args) {
	args->use_cached = FALSE;
	finalize_cache_reading(args);
}

static void disable_cache_writing(struct planetary_cache *args) {
	args->cache_data = FALSE;
	finalize_cache_writing(args);
}

/* get gaussian filtered image data for the monochrome image given as index in
 * the sequence if cached or as fits */
WORD * get_gaussian_data_for_image(int index, fits *fit, struct planetary_cache *args) {
	if (args->use_cached && args->seq_gaussian) {
		// getting the data from cache
		fits cachedfit = { 0 };
		if (seq_read_frame(args->seq_gaussian, index, &cachedfit)) {
			disable_cache_reading(args);
			goto computing_it;
		}
		WORD *buf = cachedfit.data;
		cachedfit.data = NULL;
		clearfits(&cachedfit);
		return buf;
	}
	else {
		// compute gaussian for fit
computing_it:
		if (!fit)
			return NULL;
		// convert to monochrome first
		WORD *monochrome_data;
		if (fit->naxes[2] == 3) {
			monochrome_data = malloc(fit->rx * fit->ry * sizeof(WORD));
			cvToMonochrome(fit->pdata, fit->rx, fit->ry, monochrome_data);
		} else {
			monochrome_data = fit->data;
		}

		WORD *buf = malloc(fit->rx * fit->ry * sizeof(WORD));
		cvGaussian(monochrome_data, fit->rx, fit->ry, args->kernel_size, buf);

		// store the data in cache, with a fits envelope
		if (args->cache_data && args->ser_gaussian) {
			fits cachefit = { 0 }, *cachefitptr = &cachefit;
			if (new_fit_image(&cachefitptr, fit->rx, fit->ry, 1, buf)) {
				disable_cache_writing(args);
			}
			else ser_write_frame_from_fit(args->ser_gaussian, cachefitptr, index);
			cachefit.data = NULL;
			clearfits(cachefitptr);
		}
		return buf;
	}
}

WORD * get_laplacian_data_for_image(int index, WORD *gaussian_data,
		int width, int height, struct planetary_cache *args) {
	if (args->use_cached && args->seq_laplacian) {
		// getting the data from cache
		fits cachedfit = { 0 };
		if (seq_read_frame(args->seq_laplacian, index, &cachedfit)) {
			disable_cache_reading(args);
			goto computing_it;
		}
		WORD *buf = cachedfit.data;
		cachedfit.data = NULL;
		clearfits(&cachedfit);
		return buf;
	}
	else {
		// compute laplacian for gaussian_data
computing_it:
		{
			WORD *buf = malloc(width * height * sizeof(WORD));
			cvLaplacian(gaussian_data, width, height, args->kernel_size, buf);

			// store the data in cache, with a fits envelope
			if (args->cache_data && args->ser_laplacian) {
				fits cachefit = { 0 }, *cachefitptr = &cachefit;
				if (new_fit_image(&cachefitptr, width, height, 1, buf)) {
					disable_cache_writing(args);
				}
				else ser_write_frame_from_fit(args->ser_laplacian, cachefitptr, index);
				cachefit.data = NULL;
				clearfits(cachefitptr);
			}
			return buf;
		}
	}
}


