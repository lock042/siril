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
/*
 The following code is based on algorithms written by Richard White at STScI and made
 available for use in CFITSIO in July 1999 and updated in January 2008.
 The code has been updated by Cyril Richard (2016) in order to work with ushort data used by Siril
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "core/proto.h"
#include "core/siril.h"
#include "core/processing.h"
#include "sorting.h"

/* more than this many standard deviations from the mean is an outlier */
#define SIGMA_CLIP     5.
#define NITER          3	/* number of sigma-clipping iterations */

static int FnMeanSigma_ushort(WORD *array, long npix, 
		long *ngoodpix, double *mean, double *sigma,
		int *status);
static int FnMeanSigma_float(float *array, long npix, 
		long *ngoodpix, double *mean, double *sigma,
		int *status);

static int FnMeanSigma_int(int *array, long npix,
		long *ngoodpix, double *mean, double *sigma, int *status);

static int FnNoise1_ushort(WORD *array, long nx, long ny, 
		double *noise, threading_type threads, int *status);

static int FnNoise1_float(float *array, long nx, long ny, 
			double *noise, threading_type threads, int *status);


static int FnNoise5_ushort(WORD *array, long nx, long ny, 
		long *ngood, WORD *minval, WORD *maxval, double *n2,
		double *n3, double *n5, int *status);

static int FnNoise5_float(float *array, long nx, long ny, 
		long *ngood, float *minval, float *maxval,
		double *n2, double *n3, double *n5, int *status);

static int FnCompare_double(const void *, const void *);

/*--------------------------------------------------------------------------*/
int siril_fits_img_stats_ushort(WORD *array, /*  2 dimensional array of image pixels */
long nx, /* number of pixels in each row of the image */
long ny, /* number of rows in the image */
/* (if this is a 3D image, then ny should be the */
/* product of the no. of rows times the no. of planes) */

/* returned parameters (if the pointer is not null)  */
long *ngoodpix, /* number of non-null pixels in the image */
WORD *minvalue, /* returned minimum non-null value in the array */
WORD *maxvalue, /* returned maximum non-null value in the array */
double *mean, /* returned mean value of all non-null pixels */
double *sigma, /* returned R.M.S. value of all non-null pixels */
double *noise1, /* 1st order estimate of noise in image background level */
double *noise2, /* 2nd order estimate of noise in image background level */
double *noise3, /* 3rd order estimate of noise in image background level */
double *noise5, /* 5th order estimate of noise in image background level */
threading_type threads,
int *status) /* error status */

/*
 Compute statistics of the input ushort integer image.
 */
{
	long ngood = 0;
	WORD minval = 0, maxval = 0;
	double xmean = 0., xsigma = 0., xnoise = 0., xnoise2 = 0., xnoise3 = 0.,
			xnoise5 = 0.;

	/* need to calculate mean and/or sigma and/or limits? */
	if (mean || sigma) {
		FnMeanSigma_ushort(array, nx * ny, &ngood, &xmean,
				&xsigma, status);

		if (ngoodpix)
			*ngoodpix = ngood;
		if (mean)
			*mean = xmean;
		if (sigma)
			*sigma = xsigma;
	}

	if (noise1) {
		FnNoise1_ushort(array, nx, ny, &xnoise, threads, status);

		*noise1 = xnoise;
	}

	/* this following code is not used for now */
	if (minvalue || maxvalue || noise2 || noise3 || noise5) {
		FnNoise5_ushort(array, nx, ny, &ngood, &minval,
				&maxval, &xnoise2, &xnoise3, &xnoise5, status);

		if (ngoodpix)
			*ngoodpix = ngood;
		if (minvalue)
			*minvalue = minval;
		if (maxvalue)
			*maxvalue = maxval;
		if (noise2)
			*noise2 = xnoise2;
		if (noise3)
			*noise3 = xnoise3;
		if (noise5)
			*noise5 = xnoise5;
	}
	return (*status);
}

int siril_fits_img_stats_float(float *array, /*  2 dimensional array of image pixels */
long nx, /* number of pixels in each row of the image */
long ny, /* number of rows in the image */
/* (if this is a 3D image, then ny should be the */
/* product of the no. of rows times the no. of planes) */

/* returned parameters (if the pointer is not null)  */
long *ngoodpix, /* number of non-null pixels in the image */
float *minvalue, /* returned minimum non-null value in the array */
float *maxvalue, /* returned maximum non-null value in the array */
double *mean, /* returned mean value of all non-null pixels */
double *sigma, /* returned R.M.S. value of all non-null pixels */
double *noise1, /* 1st order estimate of noise in image background level */
double *noise2, /* 2nd order estimate of noise in image background level */
double *noise3, /* 3rd order estimate of noise in image background level */
double *noise5, /* 5th order estimate of noise in image background level */
threading_type threads,
int *status) /* error status */

/*
 Compute statistics of the input float image.
 */
{
	long ngood = 0L;
	float minval = 0.f, maxval = 0.f;
	double xmean = 0., xsigma = 0., xnoise = 0., xnoise2 = 0., xnoise3 = 0.,
			xnoise5 = 0.;

	/* need to calculate mean and/or sigma and/or limits? */
	if (mean || sigma) {
		FnMeanSigma_float(array, nx * ny, &ngood, &xmean,
				&xsigma, status);

		if (ngoodpix)
			*ngoodpix = ngood;
		if (mean)
			*mean = xmean;
		if (sigma)
			*sigma = xsigma;
	}

	if (noise1) {
		FnNoise1_float(array, nx, ny, &xnoise, threads, status);

		*noise1 = xnoise;
	}

	if (minvalue || maxvalue || noise3) {
		FnNoise5_float(array, nx, ny, &ngood, &minval,
				&maxval, &xnoise2, &xnoise3, &xnoise5, status);

		if (ngoodpix)
			*ngoodpix = ngood;
		if (minvalue)
			*minvalue = minval;
		if (maxvalue)
			*maxvalue = maxval;
		if (noise2)
			*noise2 = xnoise2;
		if (noise3)
			*noise3 = xnoise3;
		if (noise5)
			*noise5 = xnoise5;
	}
	return (*status);
}

/*--------------------------------------------------------------------------*/
static int FnMeanSigma_ushort(WORD *array, /*  2 dimensional array of image pixels */
long npix, /* number of pixels in the image */

/* returned parameters */

long *ngoodpix, /* number of non-null pixels in the image */
double *mean, /* returned mean value of all non-null pixels */
double *sigma, /* returned R.M.S. value of all non-null pixels */
int *status) /* error status */

/*
 Compute mean and RMS sigma of the non-null pixels in the input array.
 */
{
	long ii, ngood = 0;
	WORD *value;
	double sum = 0., sum2 = 0., xtemp;

	value = array;


	for (ii = 0; ii < npix; ii++, value++) {
		if (*value != 0) {
			ngood++;
			xtemp = (double) *value;
			sum += xtemp;
			sum2 += (xtemp * xtemp);
		}
	}

	if (ngood > 1) {
		if (ngoodpix)
			*ngoodpix = ngood;
		xtemp = sum / ngood;
		if (mean)
			*mean = xtemp;
		if (sigma)
			*sigma = sqrt((sum2 / ngood) - (xtemp * xtemp));
	} else if (ngood == 1) {
		if (ngoodpix)
			*ngoodpix = 1;
		if (mean)
			*mean = sum;
		if (sigma)
			*sigma = 0.0;
	} else {
		if (ngoodpix)
			*ngoodpix = 0;
		if (mean)
			*mean = 0.;
		if (sigma)
			*sigma = 0.;
	}
	return (*status);
}

/*--------------------------------------------------------------------------*/
static int FnMeanSigma_float(float *array, /*  2 dimensional array of image pixels */
long npix, /* number of pixels in the image */

/* returned parameters */

long *ngoodpix, /* number of non-null pixels in the image */
double *mean, /* returned mean value of all non-null pixels */
double *sigma, /* returned R.M.S. value of all non-null pixels */
int *status) /* error status */

/*
 Compute mean and RMS sigma of the non-null pixels in the input array.
 */
{
	long ii, ngood = 0;
	float *value;
	double sum = 0., sum2 = 0., xtemp;

	value = array;

	for (ii = 0; ii < npix; ii++, value++) {
		if (*value != 0.f && !isnan(*value)) {
			ngood++;
			xtemp = (double) *value;
			sum += xtemp;
			sum2 += (xtemp * xtemp);
		}
	}

	if (ngood > 1) {
		if (ngoodpix)
			*ngoodpix = ngood;
		xtemp = sum / ngood;
		if (mean)
			*mean = xtemp;
		if (sigma)
			*sigma = sqrt((sum2 / ngood) - (xtemp * xtemp));
	} else if (ngood == 1) {
		if (ngoodpix)
			*ngoodpix = 1;
		if (mean)
			*mean = sum;
		if (sigma)
			*sigma = 0.0;
	} else {
		if (ngoodpix)
			*ngoodpix = 0;
		if (mean)
			*mean = 0.;
		if (sigma)
			*sigma = 0.;
	}
	return (*status);
}


/*--------------------------------------------------------------------------*/
static int FnMeanSigma_int(int *array, /*  2 dimensional array of image pixels */
long npix, /* number of pixels in the image */

/* returned parameters */

long *ngoodpix, /* number of non-null pixels in the image */
double *mean, /* returned mean value of all non-null pixels */
double *sigma, /* returned R.M.S. value of all non-null pixels */
int *status) /* error status */

/*
 Compute mean and RMS sigma of the non-null pixels in the input array.
 */
{
	long ii, ngood = 0;
	int *value;
	double sum = 0., sum2 = 0., xtemp;

	value = array;


	for (ii = 0; ii < npix; ii++, value++) {
		if (*value != 0) {
			ngood++;
			xtemp = (double) *value;
			sum += xtemp;
			sum2 += (xtemp * xtemp);
		}
	}

	if (ngood > 1) {
		if (ngoodpix)
			*ngoodpix = ngood;
		xtemp = sum / ngood;
		if (mean)
			*mean = xtemp;
		if (sigma)
			*sigma = sqrt((sum2 / ngood) - (xtemp * xtemp));
	} else if (ngood == 1) {
		if (ngoodpix)
			*ngoodpix = 1;
		if (mean)
			*mean = sum;
		if (sigma)
			*sigma = 0.0;
	} else {
		if (ngoodpix)
			*ngoodpix = 0;
		if (mean)
			*mean = 0.;
		if (sigma)
			*sigma = 0.;
	}
	return (*status);
}
/*--------------------------------------------------------------------------*/

static int FnNoise5_ushort(WORD *array, /*  2 dimensional array of image pixels */
long nx, /* number of pixels in each row of the image */
long ny, /* number of rows in the image */
/* returned parameters */
long *ngood, /* number of good, non-null pixels? */
WORD *minval, /* minimum non-null value */
WORD *maxval, /* maximum non-null value */
double *noise2, /* returned 2nd order MAD of all non-null pixels */
double *noise3, /* returned 3rd order MAD of all non-null pixels */
double *noise5, /* returned 5th order MAD of all non-null pixels */
int *status) /* error status */

/*
 Estimate the median and background noise in the input image using 2nd, 3rd and 5th
 order Median Absolute Differences.

 The noise in the background of the image is calculated using the MAD algorithms
 developed for deriving the signal to noise ratio in spectra
 (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)

 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))

 The returned estimates are the median of the values that are computed for each
 row of the image.
 */
{
	long ii, jj, nrows = 0, nrows2 = 0, nvals, nvals2, ngoodpix = 0;
	int *differences2, *differences3, *differences5;
	WORD *rowpix, v1, v2, v3, v4, v5, v6, v7, v8, v9;
	WORD xminval = USHRT_MAX, xmaxval = 0;
	int do_range = 0;
	double *diffs2, *diffs3, *diffs5;
	double xnoise2 = 0, xnoise3 = 0, xnoise5 = 0;

	if (nx < 9) {
		/* treat entire array as an image with a single row */
		nx = nx * ny;
		ny = 1;
	}

	/* rows must have at least 9 pixels */
	if (nx < 9) {

		for (ii = 0; ii < nx; ii++) {
			if (array[ii] == 0)
				continue;
			else {
				if (array[ii] < xminval)
					xminval = array[ii];
				if (array[ii] > xmaxval)
					xmaxval = array[ii];
				ngoodpix++;
			}
		}
		if (minval)
			*minval = xminval;
		if (maxval)
			*maxval = xmaxval;
		if (ngood)
			*ngood = ngoodpix;
		if (noise2)
			*noise2 = 0.;
		if (noise3)
			*noise3 = 0.;
		if (noise5)
			*noise5 = 0.;
		return (*status);
	}

	/* do we need to compute the min and max value? */
	if (minval || maxval)
		do_range = 1;

	/* allocate arrays used to compute the median and noise estimates */
	differences2 = calloc(nx, sizeof(int));
	if (!differences2) {
		*status = MEMORY_ALLOCATION;
		return (*status);
	}
	differences3 = calloc(nx, sizeof(int));
	if (!differences3) {
		free(differences2);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}
	differences5 = calloc(nx, sizeof(int));
	if (!differences5) {
		free(differences2);
		free(differences3);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	diffs2 = calloc(ny, sizeof(double));
	if (!diffs2) {
		free(differences2);
		free(differences3);
		free(differences5);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	diffs3 = calloc(ny, sizeof(double));
	if (!diffs3) {
		free(differences2);
		free(differences3);
		free(differences5);
		free(diffs2);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	diffs5 = calloc(ny, sizeof(double));
	if (!diffs5) {
		free(differences2);
		free(differences3);
		free(differences5);
		free(diffs2);
		free(diffs3);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	/* loop over each row of the image */
	for (jj = 0; jj < ny; jj++) {

		rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;

		while (ii < nx && rowpix[ii] == 0)
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v1 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v1 < xminval)
				xminval = v1;
			if (v1 > xmaxval)
				xmaxval = v1;
		}

		/***** find the 2nd valid pixel in row (which we will skip over) */
		ii++;

		while (ii < nx && rowpix[ii] == 0)
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v2 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v2 < xminval)
				xminval = v2;
			if (v2 > xmaxval)
				xmaxval = v2;
		}

		/***** find the 3rd valid pixel in row */
		ii++;
		while (ii < nx && rowpix[ii] == 0)
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v3 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v3 < xminval)
				xminval = v3;
			if (v3 > xmaxval)
				xmaxval = v3;
		}

		/* find the 4nd valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && rowpix[ii] == 0)
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v4 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v4 < xminval)
				xminval = v4;
			if (v4 > xmaxval)
				xmaxval = v4;
		}

		/* find the 5th valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && rowpix[ii] == 0)
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v5 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v5 < xminval)
				xminval = v5;
			if (v5 > xmaxval)
				xmaxval = v5;
		}

		/* find the 6th valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && rowpix[ii] == 0)
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v6 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v6 < xminval)
				xminval = v6;
			if (v6 > xmaxval)
				xmaxval = v6;
		}

		/* find the 7th valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && rowpix[ii] == 0)
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v7 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v7 < xminval)
				xminval = v7;
			if (v7 > xmaxval)
				xmaxval = v7;
		}

		/* find the 8th valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && rowpix[ii] == 0)
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v8 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v8 < xminval)
				xminval = v8;
			if (v8 > xmaxval)
				xmaxval = v8;
		}
		/* now populate the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		nvals2 = 0;
		for (ii++; ii < nx; ii++) {

			/* find the next valid pixel in row */
			while (ii < nx && rowpix[ii] == 0)
				ii++;

			if (ii == nx)
				break; /* hit end of row */
			v9 = rowpix[ii]; /* store the good pixel value */

			if (do_range) {
				if (v9 < xminval)
					xminval = v9;
				if (v9 > xmaxval)
					xmaxval = v9;
			}

			/* construct array of absolute differences */

			if (!(v5 == v6 && v6 == v7)) {
				differences2[nvals2] = abs((int) v5 - (int) v7);
				nvals2++;
			}

			if (!(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7)) {
				differences3[nvals] = abs((2 * (int) v5) - (int) v3 - (int) v7);
				differences5[nvals] = abs(
						(6 * (int) v5) - (4 * (int) v3) - (4 * (int) v7)
								+ (int) v1 + (int) v9);
				nvals++;
			} else {
				/* ignore constant background regions */
				ngoodpix++;
			}

			/* shift over 1 pixel */
			v1 = v2;
			v2 = v3;
			v3 = v4;
			v4 = v5;
			v5 = v6;
			v6 = v7;
			v7 = v8;
			v8 = v9;
		} /* end of loop over pixels in the row */

		/* compute the median diffs */
		/* Note that there are 8 more pixel values than there are diffs values. */
		ngoodpix += nvals;

		if (nvals == 0) {
			continue; /* cannot compute medians on this row */
		} else if (nvals == 1) {
			if (nvals2 == 1) {
				diffs2[nrows2] = differences2[0];
				nrows2++;
			}

			diffs3[nrows] = differences3[0];
			diffs5[nrows] = differences5[0];
		} else {
			// quickmedian returns the median MUCH faster than using quicksort
			if (nvals2 > 1) {
				diffs2[nrows2] = quickmedian_int(differences2, nvals);
				nrows2++;
			}

			diffs3[nrows] = quickmedian_int(differences3, nvals);
			diffs5[nrows] = quickmedian_int(differences5, nvals);
		}

		nrows++;
	} /* end of loop over rows */

	/* compute median of the values for each row */
	if (nrows == 0) {
		xnoise3 = 0;
		xnoise5 = 0;
	} else if (nrows == 1) {
		xnoise3 = diffs3[0];
		xnoise5 = diffs5[0];
	} else {
		qsort(diffs3, nrows, sizeof(double), FnCompare_double);
		qsort(diffs5, nrows, sizeof(double), FnCompare_double);
		xnoise3 = (diffs3[(nrows - 1) / 2] + diffs3[nrows / 2]) / 2.;
		xnoise5 = (diffs5[(nrows - 1) / 2] + diffs5[nrows / 2]) / 2.;
	}

	if (nrows2 == 0) {
		xnoise2 = 0;
	} else if (nrows2 == 1) {
		xnoise2 = diffs2[0];
	} else {
		qsort(diffs2, nrows2, sizeof(double), FnCompare_double);
		xnoise2 = (diffs2[(nrows2 - 1) / 2] + diffs2[nrows2 / 2]) / 2.;
	}

	if (ngood)
		*ngood = ngoodpix;
	if (minval)
		*minval = xminval;
	if (maxval)
		*maxval = xmaxval;
	if (noise2)
		*noise2 = 1.0483579 * xnoise2;
	if (noise3)
		*noise3 = 0.6052697 * xnoise3;
	if (noise5)
		*noise5 = 0.1772048 * xnoise5;

	free(diffs5);
	free(diffs3);
	free(diffs2);
	free(differences5);
	free(differences3);
	free(differences2);

	return (*status);
}

static int FnNoise5_float(float *array, /*  2 dimensional array of image pixels */
long nx, /* number of pixels in each row of the image */
long ny, /* number of rows in the image */
/* returned parameters */
long *ngood, /* number of good, non-null pixels? */
float *minval, /* minimum non-null value */
float *maxval, /* maximum non-null value */
double *noise2, /* returned 2nd order MAD of all non-null pixels */
double *noise3, /* returned 3rd order MAD of all non-null pixels */
double *noise5, /* returned 5th order MAD of all non-null pixels */
int *status) /* error status */

/*
 Estimate the median and background noise in the input image using 2nd, 3rd and 5th
 order Median Absolute Differences.

 The noise in the background of the image is calculated using the MAD algorithms
 developed for deriving the signal to noise ratio in spectra
 (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)

 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))

 The returned estimates are the median of the values that are computed for each
 row of the image.
 */
{
	long ii, jj, nrows = 0, nrows2 = 0, nvals, nvals2, ngoodpix = 0;
	float *differences2, *differences3, *differences5;
	float *rowpix, v1, v2, v3, v4, v5, v6, v7, v8, v9;
	float xminval = FLT_MAX, xmaxval = -FLT_MAX;
	int do_range = 0;
	double *diffs2, *diffs3, *diffs5;
	double xnoise2 = 0, xnoise3 = 0, xnoise5 = 0;

	if (nx < 9) {
		/* treat entire array as an image with a single row */
		nx = nx * ny;
		ny = 1;
	}

	/* rows must have at least 9 pixels */
	if (nx < 9) {

		for (ii = 0; ii < nx; ii++) {
			if (array[ii] == 0.f || isnan(array[ii]))
				continue;
			else {
				if (array[ii] < xminval)
					xminval = array[ii];
				if (array[ii] > xmaxval)
					xmaxval = array[ii];
				ngoodpix++;
			}
		}
		if (minval)
			*minval = xminval;
		if (maxval)
			*maxval = xmaxval;
		if (ngood)
			*ngood = ngoodpix;
		if (noise2)
			*noise2 = 0.;
		if (noise3)
			*noise3 = 0.;
		if (noise5)
			*noise5 = 0.;
		return (*status);
	}

	/* do we need to compute the min and max value? */
	if (minval || maxval)
		do_range = 1;

	/* allocate arrays used to compute the median and noise estimates */
	differences2 = calloc(nx, sizeof(float));
	if (!differences2) {
		*status = MEMORY_ALLOCATION;
		return (*status);
	}
	differences3 = calloc(nx, sizeof(float));
	if (!differences3) {
		free(differences2);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}
	differences5 = calloc(nx, sizeof(float));
	if (!differences5) {
		free(differences2);
		free(differences3);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	diffs2 = calloc(ny, sizeof(double));
	if (!diffs2) {
		free(differences2);
		free(differences3);
		free(differences5);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	diffs3 = calloc(ny, sizeof(double));
	if (!diffs3) {
		free(differences2);
		free(differences3);
		free(differences5);
		free(diffs2);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	diffs5 = calloc(ny, sizeof(double));
	if (!diffs5) {
		free(differences2);
		free(differences3);
		free(differences5);
		free(diffs2);
		free(diffs3);
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	/* loop over each row of the image */
	for (jj = 0; jj < ny; jj++) {

		rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;

		while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v1 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v1 < xminval)
				xminval = v1;
			if (v1 > xmaxval)
				xmaxval = v1;
		}

		/***** find the 2nd valid pixel in row (which we will skip over) */
		ii++;
		while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v2 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v2 < xminval)
				xminval = v2;
			if (v2 > xmaxval)
				xmaxval = v2;
		}

		/***** find the 3rd valid pixel in row */
		ii++;
		while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v3 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v3 < xminval)
				xminval = v3;
			if (v3 > xmaxval)
				xmaxval = v3;
		}

		/* find the 4nd valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v4 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v4 < xminval)
				xminval = v4;
			if (v4 > xmaxval)
				xmaxval = v4;
		}

		/* find the 5th valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v5 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v5 < xminval)
				xminval = v5;
			if (v5 > xmaxval)
				xmaxval = v5;
		}

		/* find the 6th valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v6 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v6 < xminval)
				xminval = v6;
			if (v6 > xmaxval)
				xmaxval = v6;
		}

		/* find the 7th valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v7 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v7 < xminval)
				xminval = v7;
			if (v7 > xmaxval)
				xmaxval = v7;
		}

		/* find the 8th valid pixel in row (to be skipped) */
		ii++;
		while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
			ii++;

		if (ii == nx)
			continue; /* hit end of row */
		v8 = rowpix[ii]; /* store the good pixel value */
		ngoodpix++;

		if (do_range) {
			if (v8 < xminval)
				xminval = v8;
			if (v8 > xmaxval)
				xmaxval = v8;
		}
		/* now populate the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		nvals2 = 0;
		for (ii++; ii < nx; ii++) {

			/* find the next valid pixel in row */
			while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
				ii++;

			if (ii == nx)
				break; /* hit end of row */
			v9 = rowpix[ii]; /* store the good pixel value */

			if (do_range) {
				if (v9 < xminval)
					xminval = v9;
				if (v9 > xmaxval)
					xmaxval = v9;
			}

			/* construct array of absolute differences */

			if (!(v5 == v6 && v6 == v7)) {
				differences2[nvals2] = (float) fabs(v5 - v7);
				nvals2++;
			}

			if (!(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7)) {
				differences3[nvals] = (float) fabs((2 * v5) - v3 - v7);
				differences5[nvals] = (float) fabs(
						(6 * v5) - (4 * v3) - (4 * v7) + v1 + v9);
				nvals++;
			} else {
				/* ignore constant background regions */
				ngoodpix++;
			}

			/* shift over 1 pixel */
			v1 = v2;
			v2 = v3;
			v3 = v4;
			v4 = v5;
			v5 = v6;
			v6 = v7;
			v7 = v8;
			v8 = v9;
		} /* end of loop over pixels in the row */

		/* compute the median diffs */
		/* Note that there are 8 more pixel values than there are diffs values. */
		ngoodpix += nvals;

		if (nvals == 0) {
			continue; /* cannot compute medians on this row */
		} else if (nvals == 1) {
			if (nvals2 == 1) {
				diffs2[nrows2] = differences2[0];
				nrows2++;
			}

			diffs3[nrows] = differences3[0];
			diffs5[nrows] = differences5[0];
		} else {
			/* quick_select returns the median MUCH faster than using qsort */
			if (nvals2 > 1) {
				diffs2[nrows2] = quickmedian_float(differences2, nvals);
				nrows2++;
			}

			diffs3[nrows] = quickmedian_float(differences3, nvals);
			diffs5[nrows] = quickmedian_float(differences5, nvals);
		}

		nrows++;
	} /* end of loop over rows */

	/* compute median of the values for each row */
	if (nrows == 0) {
		xnoise3 = 0;
		xnoise5 = 0;
	} else if (nrows == 1) {
		xnoise3 = diffs3[0];
		xnoise5 = diffs5[0];
	} else {
		qsort(diffs3, nrows, sizeof(double), FnCompare_double);
		qsort(diffs5, nrows, sizeof(double), FnCompare_double);
		xnoise3 = (diffs3[(nrows - 1) / 2] + diffs3[nrows / 2]) / 2.;
		xnoise5 = (diffs5[(nrows - 1) / 2] + diffs5[nrows / 2]) / 2.;
	}

	if (nrows2 == 0) {
		xnoise2 = 0;
	} else if (nrows2 == 1) {
		xnoise2 = diffs2[0];
	} else {
		qsort(diffs2, nrows2, sizeof(double), FnCompare_double);
		xnoise2 = (diffs2[(nrows2 - 1) / 2] + diffs2[nrows2 / 2]) / 2.;
	}

	if (ngood)
		*ngood = ngoodpix;
	if (minval)
		*minval = xminval;
	if (maxval)
		*maxval = xmaxval;
	if (noise2)
		*noise2 = 1.0483579 * xnoise2;
	if (noise3)
		*noise3 = 0.6052697 * xnoise3;
	if (noise5)
		*noise5 = 0.1772048 * xnoise5;

	free(diffs5);
	free(diffs3);
	free(diffs2);
	free(differences5);
	free(differences3);
	free(differences2);

	return (*status);
}

/*--------------------------------------------------------------------------*/
static int FnNoise1_ushort(WORD *array, /*  2 dimensional array of image pixels */
long nx, /* number of pixels in each row of the image */
long ny, /* number of rows in the image */
/* returned parameters */
double *noise, /* returned R.M.S. value of all non-null pixels */
threading_type threads,
int *status) /* error status */
/*
 Estimate the background noise in the input image using sigma of 1st order differences.

 noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])

 The returned estimate is the median of the values that are computed for each
 row of the image.
 */
{
	long jj, nrows = 0;
	double *diffs, xnoise;
	/* rows must have at least 3 pixels to estimate noise */
	if (nx < 3) {
		*noise = 0;
		return (*status);
	}

	/* allocate arrays used to compute the median and noise estimates */
	diffs = calloc(ny, sizeof(double));
	if (!diffs) {
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	/* loop over each row of the image */
#ifdef _OPENMP
	threads = limit_threading(&threads, 200000, nx*ny);
	//if (threads > 1)
	//	siril_debug_print("threading FnNoise1_ushort with %d threads\n", threads);
#pragma omp parallel num_threads(threads) if (threads>1)
#endif
	{
#ifdef _OPENMP
		if (threads > 1 && omp_get_num_threads() != threads)
			siril_debug_print("actual number of threads: %d of %d requested (level %d)\n", omp_get_num_threads(), threads, omp_get_level());
#endif
		WORD *rowpix, v1;
		double mean, stdev;
		int *differences;
		differences = calloc(nx, sizeof(int)); // no check here at the moment, allocation is small, should be no problem

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (jj = 0; jj < ny; jj++) {
				long ii, kk, nvals;
				rowpix = array + (jj * nx); /* point to first pixel in the row */
				int iter;

				/***** find the first valid pixel in row */
				ii = 0;

				while (ii < nx && rowpix[ii] == 0)
					ii++;

				if (ii == nx)
					continue; /* hit end of row */
				v1 = rowpix[ii]; /* store the good pixel value */

				/* now continue populating the differences arrays */
				/* for the remaining pixels in the row */
				nvals = 0;
				for (ii++; ii < nx; ii++) {

					/* find the next valid pixel in row */
					while (ii < nx && rowpix[ii] == 0)
						ii++;

					if (ii == nx)
						break; /* hit end of row */

					/* construct array of 1st order differences */
					differences[nvals] = v1 - rowpix[ii];

					nvals++;
					/* shift over 1 pixel */
					v1 = rowpix[ii];
				} /* end of loop over pixels in the row */

				if (nvals < 2)
					continue;
				else {

					FnMeanSigma_int(differences, nvals, NULL, &mean, &stdev,
							status);

					if (stdev > 0.) {
						for (iter = 0; iter < NITER; iter++) {
							kk = 0;
							for (ii = 0; ii < nvals; ii++) {
								if (fabs(differences[ii] - mean)
										< SIGMA_CLIP * stdev) {
									if (kk < ii)
										differences[kk] = differences[ii];
									kk++;
								}
							}
							if (kk == nvals)
								break;

							nvals = kk;
							FnMeanSigma_int(differences, nvals, NULL, &mean,
									&stdev, status);
						}
					}

#ifdef _OPENMP
#pragma omp critical
#endif
					{
						diffs[nrows] = stdev;
						nrows++;
					}
				}
			} /* end of loop over rows */
		free(differences);

	}
	/* compute median of the values for each row */
	if (nrows == 0) {
		xnoise = 0;
	} else if (nrows == 1) {
		xnoise = diffs[0];
	} else {
		xnoise = quickmedian_double(diffs, nrows);
	}

	*noise = .70710678 * xnoise;

	free(diffs);

	return (*status);
}

static int FnNoise1_float(float *array, /*  2 dimensional array of image pixels */
long nx, /* number of pixels in each row of the image */
long ny, /* number of rows in the image */
/* returned parameters */
double *noise, /* returned R.M.S. value of all non-null pixels */
threading_type threads,
int *status)        /* error status */
/*
Estimate the background noise in the input image using sigma of 1st order differences.

  noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])

The returned estimate is the median of the values that are computed for each
row of the image.
*/
{
	long jj, nrows = 0;
	double *diffs, xnoise;

	/* rows must have at least 3 pixels to estimate noise */
	if (nx < 3) {
		*noise = 0;
		return (*status);
	}

//	/* allocate arrays used to compute the median and noise estimates */
//	differences = calloc(nx, sizeof(float));
//	if (!differences) {
//		*status = MEMORY_ALLOCATION;
//		return (*status);
//	}

	diffs = calloc(ny, sizeof(double));
	if (!diffs) {
		*status = MEMORY_ALLOCATION;
		return (*status);
	}

	/* loop over each row of the image */
#ifdef _OPENMP
	threads = limit_threading(&threads, 200000, nx*ny);
	//if (threads > 1)
		//siril_debug_print("threading FnNoise1_float with %d threads\n", threads);
#pragma omp parallel num_threads(threads) if (threads>1)
#endif
	{
#ifdef _OPENMP
		if (threads > 1 && omp_get_num_threads() != threads)
			siril_debug_print("actual number of threads: %d of %d requested (level %d)\n", omp_get_num_threads(), threads, omp_get_level());
#endif
		float *rowpix, v1;
		double mean, stdev;
		float *differences;
		differences = calloc(nx, sizeof(float)); // no check here at the moment, allocation is small, should be no problem
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
		for (jj = 0; jj < ny; jj++) {
			long ii, kk, nvals;
			rowpix = array + (jj * nx); /* point to first pixel in the row */
			int iter;

			/***** find the first valid pixel in row */
			ii = 0;

			while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
				ii++;

			if (ii == nx)
				continue; /* hit end of row */
			v1 = rowpix[ii]; /* store the good pixel value */

			/* now continue populating the differences arrays */
			/* for the remaining pixels in the row */
			nvals = 0;
			for (ii++; ii < nx; ii++) {

				/* find the next valid pixel in row */
				while (ii < nx && (rowpix[ii] == 0.f || isnan(rowpix[ii])))
					ii++;

				if (ii == nx)
					break; /* hit end of row */

				/* construct array of 1st order differences */
				differences[nvals] = v1 - rowpix[ii];

				nvals++;
				/* shift over 1 pixel */
				v1 = rowpix[ii];
			} /* end of loop over pixels in the row */

			if (nvals < 2)
				continue;
			else {

				FnMeanSigma_float(differences, nvals, NULL, &mean, &stdev,
						status);

				if (stdev > 0.) {
					for (iter = 0; iter < NITER; iter++) {
						kk = 0;
						for (ii = 0; ii < nvals; ii++) {
							if (fabs(differences[ii] - (float)mean)
									< SIGMA_CLIP * stdev) {
								if (kk < ii)
									differences[kk] = differences[ii];
								kk++;
							}
						}
						if (kk == nvals)
							break;

						nvals = kk;
						FnMeanSigma_float(differences, nvals, NULL, &mean,
								&stdev, status);
					}
				}

#ifdef _OPENMP
#pragma omp critical
#endif
				{
					diffs[nrows] = stdev;
					nrows++;
				}
			}
		} /* end of loop over rows */
		free(differences);
	}
	/* compute median of the values for each row */
	if (nrows == 0) {
		xnoise = 0;
	} else if (nrows == 1) {
		xnoise = diffs[0];
	} else {
		qsort(diffs, nrows, sizeof(double), FnCompare_double);
		xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;
	}

	*noise = .70710678 * xnoise;

	free(diffs);

	return (*status);
}
/*--------------------------------------------------------------------------*/

static int FnCompare_double(const void *v1, const void *v2) {
	const double *i1 = v1;
	const double *i2 = v2;

	if (*i1 < *i2)
		return (-1);
	else if (*i1 > *i2)
		return (1);
	else
		return (0);
}
