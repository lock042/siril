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

/* This files contains all sorting algorithms of siril.
 * See src/test/sorting.c for testing and metrics.
 */

#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"

#include "rt/rt_algo.h"
#include "sorting.h"

/**
 * In-place insertion sort of array of double a of size n
 * @param a array to sort
 * @param n size of the array
 */
static void insertionSort_d(double a[], size_t n) {
	for (long i = 1; i < n; i++) {
		const double val = a[i];
		long j = i - 1;

		/* Move elements of a[0..i-1], that are greater than val, to one position ahead of their current position */
		while (j >= 0 && a[j] > val) {
			a[j + 1] = a[j];
			--j;
		}
		a[j + 1] = val;
	}
}

/**
 * In-place quick sort of array of double a of size n
 * @param a array to sort
 * @param n size of the array
 */
void quicksort_d (double *a, size_t n) {
	if (n <= 32) {
		return insertionSort_d(a, n);
	}

	double pivot = a[n / 2];
	double *left = a;
	double *right = a + n - 1;
	register double t;

	while (left <= right) {
		if (*left < pivot) {
			left++;
			continue;
		}
		if (*right > pivot) {
			right--;
			continue;
		}
		t = *left;
		*left++ = *right;
		*right-- = t;
	}
	quicksort_d(a, right - a + 1);
	quicksort_d(left, a + n - left);
}

/**
 * In-place insertion sort of array of float a of size n
 * @param a array to sort
 * @param n size of the array
 */
 static void insertionSort_f(float a[], size_t n) {
	for (long i = 1; i < n; i++) {
		const float val = a[i];
		long j = i - 1;

		/* Move elements of a[0..i-1], that are greater than val, to one position ahead of their current position */
		while (j >= 0 && a[j] > val) {
			a[j + 1] = a[j];
			--j;
		}
		a[j + 1] = val;
	}
}

/**
 * In-place quick sort of array of float a of size n
 * @param a array to sort
 * @param n size of the array
 */
void quicksort_f (float *a, size_t n) {
	if (n <= 32) {
		return insertionSort_f(a, n);
	}

	float pivot = a[n / 2];
	float *left = a;
	float *right = a + n - 1;
	register float t;

	while (left <= right) {
		if (*left < pivot) {
			left++;
			continue;
		}
		if (*right > pivot) {
			right--;
			continue;
		}
		t = *left;
		*left++ = *right;
		*right-- = t;
	}
	quicksort_f(a, right - a + 1);
	quicksort_f(left, a + n - left);
}

/**
 * In-place insertion sort of array of WORD a of size n
 * @param a array to sort
 * @param n size of the array
 */
static void insertionSort_s(WORD a[], size_t n) {
	for (long i = 1; i < n; i++) {
		const WORD val = a[i];
		long j = i - 1;

		/* Move elements of a[0..i-1], that are greater than val, to one position ahead of their current position */
		while (j >= 0 && a[j] > val) {
			a[j + 1] = a[j];
			--j;
		}
		a[j + 1] = val;
	}
}

/**
 * In-place quick sort of array of WORD a of size n
 * @param a array to sort
 * @param n size of the array
 */
void quicksort_s(WORD *a, size_t n) {
	if (n <= 32) {
		return insertionSort_s(a, n);
	}
	WORD pivot = a[n / 2];
	WORD *left = a;
	WORD *right = a + n - 1;
	register WORD t;

	while (left <= right) {
		if (*left < pivot) {
			left++;
			continue;
		}
		if (*right > pivot) {
			right--;
			continue;
		}
		t = *left;
		*left++ = *right;
		*right-- = t;
	}
	quicksort_s(a, right - a + 1);
	quicksort_s(left, a + n - left);
}

/* quickmedian returns the median from array of length n
 * Derived from the original quickselect algorithm from Hoare
 * warning: data are sorted in place
 * non recurssive version modified to return median value
 * @param a array of WORD to search
 * @param n size of the array
 * @return median as double for even size average the middle two elements
*/
double quickmedian(WORD *a, size_t n) {
	// Use faster and robust sorting network for small size array
	if (n < 9)
		return sortnet_median(a, n);

	WORD pivot, tmp;
	size_t k = n / 2;	// size to sort
	size_t left = 0; 	// left index
	size_t right = n - 1; 	// right index

	while (left < right) { //we stop when our indicies have crossed
		size_t pindex = (left + right) / 2; // pivot selection, this can be whatever
		pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP(pivot,right)

		for (size_t i = pindex = left; i < right; i++) {
			if (a[i] < pivot) { // SWAP
				tmp = a[pindex];
				a[pindex] = a[i];
				a[i] = tmp;
				pindex++;
			}
		}
		a[right] = a[pindex];
		a[pindex] = pivot; // SWAP

		if (pindex < k)
			left = pindex + 1;
		else
			// pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ?
			((double) a[k - 1] + (double) a[k]) / 2.0 : (double) a[k];
}

/* quickmedian_float returns the median from array of length n
 * Derived from the original quickselect algorithm from Hoare
 * warning: data are sorted in place
 * non recurssive version modified to return median value
 * @param a array of float to search
 * @param n size of the array
 * @return median as double for even size average the middle two elements
 */
double quickmedian_float(float *a, size_t n) {
	// Use faster and robust sorting network for small size array
	if (n < 9)
		return sortnet_median_float(a, n);
	size_t k = n / 2;	// size to sort
	size_t left = 0; 	// left index
	size_t right = n - 1; 	// right index
	float tmp;

	while (left < right) { //we stop when our indicies have crossed
		size_t pindex = (left + right) / 2; // pivot selection, this can be whatever
		float pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP(pivot,right)

		for (size_t i = pindex = left; i < right; i++) {
			if (a[i] < pivot) { // SWAP
				tmp = a[pindex];
				a[pindex] = a[i];
				a[i] = tmp;
				pindex++;
			}
		}
		a[right] = a[pindex];
		a[pindex] = pivot; // SWAP(right,j)

		if (pindex < k)
			left = pindex + 1;
		else
			// pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ? ((double) a[k - 1] + a[k]) / 2.0 : (double) a[k];
}

/* quickmedian_double returns the median from array of length n
 * Derived from the original quickselect algorithm from Hoare
 * warning: data are sorted in place
 * non recurssive version modified to return median value
 * @param a array of double to search
 * @param n size of the array
 * @return median as double for even size average the middle two elements
 */
double quickmedian_double(double *a, size_t n) {
	// Use faster and robust sorting network for small size array
	if (n < 9)
		return sortnet_median_double(a, n);
	size_t k = n / 2;	// size to sort
	size_t left = 0; 	// left index
	size_t right = n - 1; 	// right index
	double tmp;

	while (left < right) { //we stop when our indicies have crossed
		size_t pindex = (left + right) / 2; // pivot selection, this can be whatever
		double pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP(pivot,right)

		for (size_t i = pindex = left; i < right; i++) {
			if (a[i] < pivot) { // SWAP
				tmp = a[pindex];
				a[pindex] = a[i];
				a[i] = tmp;
				pindex++;
			}
		}
		a[right] = a[pindex];
		a[pindex] = pivot; // SWAP(right,j)

		if (pindex < k)
			left = pindex + 1;
		else
			// pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2.0 : a[k];
}

/*
 * quickmedian_int returns the median from array of int of of length n
 * Derived from original quickselect algorithm from Hoare
 * @param a array to search
 * @param n size of the array
 * @return median as double
 */
double quickmedian_int(int *a, size_t n) {
	size_t k = n / 2;	// size to sort
	size_t left = 0; 	// left index
	size_t right = n - 1; 	// right index
	int tmp;

	while (left < right) { //we stop when our indicies have crossed
		size_t pindex = (left + right) / 2; // pivot selection, this can be whatever
		int pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP

		for (size_t i = pindex = left; i < right; i++) {
			if (a[i] < pivot) { // SWAP
				tmp = a[pindex];
				a[pindex] = a[i];
				a[i] = tmp;
				pindex++;
			}
		}
		a[right] = a[pindex];
		a[pindex] = pivot; // SWAP

		if (pindex < k)
			left = pindex + 1;
		else
			// pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ?
			((double) a[k - 1] + (double) a[k]) / 2.0 : (double) a[k];
}

/*
 * Optimal sorting network array of size [2,9] to retrieve median value
 * (C) Emmanuel Brandt 2019-02
 * @param a array of double to sort (warning in place sorting)
 * @param n size of the array to sort [2,9]
 * warning in-place sorting
 */
#define sw(i,j) if(a[i] > a[j]) { register WORD t=a[i]; a[i]=a[j]; a[j]=t; }
double sortnet_median (WORD *a, size_t n) {
	size_t k = n / 2;

	switch (n) {
		case 1: return a[0]; break;

		case 2: sw(0,1); break;

		case 3: sw(0,1); sw(1,2); sw(0,1); break;

		case 4: sw(0,1); sw(2,3); sw(0,2); sw(1,3); sw(1,2); break;

		case 5: sw(0,1); sw(2,3); sw(1,3); sw(2,4); sw(0,2);
			sw(1,4); sw(1,2); sw(3,4); sw(2,3); break;

		case 6: sw(0,1); sw(2,3); sw(4,5); sw(0,2); sw(3,5);
			sw(1,4); sw(0,1); sw(2,3); sw(4,5); sw(1,2); sw(3,4);
			sw(2,3); break;

		case 7: sw(1,2); sw(3,4); sw(5,6); sw(0,2); sw(4,6); sw(3,5);
			sw(2,6); sw(1,5); sw(0,4); sw(2,5); sw(0,3); sw(2,4);
			sw(1,3); sw(0,1); sw(2,3); sw(4,5); break;

		case 8: sw(0,1); sw(2,3); sw(4,5); sw(6,7);
			sw(0,2); sw(1,3); sw(4,6); sw(5,7);
			sw(1,2); sw(5,6);
			sw(0,4); sw(1,5); sw(2,6); sw(3,7);
			sw(2,4); sw(3,5);
			sw(1,2); sw(3,4); sw(5,6); break;

		case 9: sw(1,8); sw(2,7); sw(3,6); sw(4,5);
			sw(1,4); sw(5,8);
			sw(0,2); sw(6,7);
			sw(2,6); sw(7,8);
			sw(0,3); sw(4,5);
			sw(0,1); sw(3,5); sw(6,7);
			sw(2,4);
			sw(1,3); sw(5,7);
			sw(4,6);
			sw(1,2); sw(3,4); sw(5,6); sw(7,8);
			sw(2,3); sw(4,5); break;

		default: return 0.0; break; // no sort
	}
	return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2.0 : a[k];
}
#undef sw

/*
 */
#define sw(i,j) if(a[i] > a[j]) { register double t=a[i]; a[i]=a[j]; a[j]=t; }
double sortnet_median_double(double *a, size_t n) {
	size_t k = n / 2;

	switch (n) {
		case 1: return a[0]; break;

		case 2: sw(0,1); break;

		case 3: sw(0,1); sw(1,2); sw(0,1); break;

		case 4: sw(0,1); sw(2,3); sw(0,2); sw(1,3); sw(1,2); break;

		case 5: sw(0,1); sw(2,3); sw(1,3); sw(2,4); sw(0,2);
			sw(1,4); sw(1,2); sw(3,4); sw(2,3); break;

		case 6: sw(0,1); sw(2,3); sw(4,5); sw(0,2); sw(3,5);
			sw(1,4); sw(0,1); sw(2,3); sw(4,5); sw(1,2); sw(3,4);
			sw(2,3); break;

		case 7: sw(1,2); sw(3,4); sw(5,6); sw(0,2); sw(4,6); sw(3,5);
			sw(2,6); sw(1,5); sw(0,4); sw(2,5); sw(0,3); sw(2,4);
			sw(1,3); sw(0,1); sw(2,3); sw(4,5); break;

		case 8: sw(0,1); sw(2,3); sw(4,5); sw(6,7);
			sw(0,2); sw(1,3); sw(4,6); sw(5,7);
			sw(1,2); sw(5,6);
			sw(0,4); sw(1,5); sw(2,6); sw(3,7);
			sw(2,4); sw(3,5);
			sw(1,2); sw(3,4); sw(5,6); break;

		case 9: sw(1,8); sw(2,7); sw(3,6); sw(4,5);
			sw(1,4); sw(5,8);
			sw(0,2); sw(6,7);
			sw(2,6); sw(7,8);
			sw(0,3); sw(4,5);
			sw(0,1); sw(3,5); sw(6,7);
			sw(2,4);
			sw(1,3); sw(5,7);
			sw(4,6);
			sw(1,2); sw(3,4); sw(5,6); sw(7,8);
			sw(2,3); sw(4,5); break;

		default: return 0.0; break; // no sort
	}
	return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2.0 : a[k];
}
#undef sw

/*
 */
#define sw(i,j) if(a[i] > a[j]) { register float t=a[i]; a[i]=a[j]; a[j]=t; }
double sortnet_median_float(float *a, size_t n) {
	size_t k = n / 2;

	switch (n) {
		case 1: return a[0]; break;

		case 2: sw(0,1); break;

		case 3: sw(0,1); sw(1,2); sw(0,1); break;

		case 4: sw(0,1); sw(2,3); sw(0,2); sw(1,3); sw(1,2); break;

		case 5: sw(0,1); sw(2,3); sw(1,3); sw(2,4); sw(0,2);
			sw(1,4); sw(1,2); sw(3,4); sw(2,3); break;

		case 6: sw(0,1); sw(2,3); sw(4,5); sw(0,2); sw(3,5);
			sw(1,4); sw(0,1); sw(2,3); sw(4,5); sw(1,2); sw(3,4);
			sw(2,3); break;

		case 7: sw(1,2); sw(3,4); sw(5,6); sw(0,2); sw(4,6); sw(3,5);
			sw(2,6); sw(1,5); sw(0,4); sw(2,5); sw(0,3); sw(2,4);
			sw(1,3); sw(0,1); sw(2,3); sw(4,5); break;

		case 8: sw(0,1); sw(2,3); sw(4,5); sw(6,7);
			sw(0,2); sw(1,3); sw(4,6); sw(5,7);
			sw(1,2); sw(5,6);
			sw(0,4); sw(1,5); sw(2,6); sw(3,7);
			sw(2,4); sw(3,5);
			sw(1,2); sw(3,4); sw(5,6); break;

		case 9: sw(1,8); sw(2,7); sw(3,6); sw(4,5);
			sw(1,4); sw(5,8);
			sw(0,2); sw(6,7);
			sw(2,6); sw(7,8);
			sw(0,3); sw(4,5);
			sw(0,1); sw(3,5); sw(6,7);
			sw(2,4);
			sw(1,3); sw(5,7);
			sw(4,6);
			sw(1,2); sw(3,4); sw(5,6); sw(7,8);
			sw(2,3); sw(4,5); break;

		default: return 0.0; break; // no sort
	}
	return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2.0 : a[k];
}
#undef sw

/*
 * Optimal sorting network [2,9]
 * (C) Emmanuel Brandt 2019-02
 * @param a array of double to sort (warning in place sorting)
 * @param n size of the array to sort [2,9]
 * warning in-place sorting
 */
#define sw(i,j) if(a[i] > a[j]) { register WORD t=a[i]; a[i]=a[j]; a[j]=t; }
void sortnet (WORD *a, size_t n) {

	switch (n) {
		case 2: sw(0,1); break;

		case 3: sw(0,1); sw(1,2); sw(0,1); break;

		case 4: sw(0,1); sw(2,3); sw(0,2); sw(1,3); sw(1,2); break;

		case 5: sw(0,1); sw(2,3); sw(1,3); sw(2,4); sw(0,2);
			sw(1,4); sw(1,2); sw(3,4); sw(2,3); break;

		case 6: sw(0,1); sw(2,3); sw(4,5); sw(0,2); sw(3,5);
			sw(1,4); sw(0,1); sw(2,3); sw(4,5); sw(1,2); sw(3,4);
			sw(2,3); break;

		case 7: sw(1,2); sw(3,4); sw(5,6); sw(0,2); sw(4,6); sw(3,5);
			sw(2,6); sw(1,5); sw(0,4); sw(2,5); sw(0,3); sw(2,4);
			sw(1,3); sw(0,1); sw(2,3); sw(4,5); break;

		case 8: sw(0,1); sw(2,3); sw(4,5); sw(6,7);
			sw(0,2); sw(1,3); sw(4,6); sw(5,7);
			sw(1,2); sw(5,6);
			sw(0,4); sw(1,5); sw(2,6); sw(3,7);
			sw(2,4); sw(3,5);
			sw(1,2); sw(3,4); sw(5,6); break;

		case 9: sw(1,8); sw(2,7); sw(3,6); sw(4,5);
			sw(1,4); sw(5,8);
			sw(0,2); sw(6,7);
			sw(2,6); sw(7,8);
			sw(0,3); sw(4,5);
			sw(0,1); sw(3,5); sw(6,7);
			sw(2,4);
			sw(1,3); sw(5,7);
			sw(4,6);
			sw(1,2); sw(3,4); sw(5,6); sw(7,8);
			sw(2,3); sw(4,5); break;

		default: break; // no sort
	}
}
#undef sw

/*
 * Histogram median for very large array of unsigned short
 * (C) Emmanuel Brandt 2019-02
 * @param a array of unsigned short to search
 * @param n size of the array
 * @return median as a double (for n odd)
 * Use temp storage h to build the histogram. Complexity O(2*N) in time and
 * O(N) in memory if single threaded or O(N+N*nb_threads) if multi-threaded
 */
double histogram_median(WORD *a, size_t n, threading_type threading) {
	// For arrays n < 10 histogram is use fast and simple sortnet_median
	if (n < 10)
		return sortnet_median(a, n);

	threading = limit_threading(&threading, 2000000, n);
#ifndef _OPENMP
	threading = SINGLE_THREADED;
#endif
	const size_t s = sizeof(unsigned int);
	unsigned int *h = (unsigned int*) calloc(USHRT_MAX + 1, s);
	if (!h) {
		PRINT_ALLOC_ERR;
		return -1.0;
	}

	if (threading == SINGLE_THREADED) {
		for (size_t i = 0; i < n; i++) {
			h[a[i]]++;
		}
	}
#ifdef _OPENMP
	else {
		/* parallel version */
#pragma omp parallel num_threads(threading) if (threading > 1)
		{
			unsigned int *hthr = (unsigned int*) calloc(USHRT_MAX + 1, s);
			if (!hthr) {
				PRINT_ALLOC_ERR;
			}
			else {
#pragma omp for schedule(static) nowait
				for (size_t i = 0; i < n; i++) {
					hthr[a[i]]++;
				}
				/* this critical block doesn't only apply to this parallel group,
				 * but to the callers as well, so it prevents this function from
				 * running in parallel with one thread for each call */
#pragma omp critical
				{
					// add per-thread histogram to main histogram
#pragma omp simd
					for (int ii = 0; ii <= USHRT_MAX; ++ii) {
						h[ii] += hthr[ii];
					}
				}
				free(hthr);
			}
		}
	}
#endif

	unsigned int i = 0, j = 0, k = n / 2;
	unsigned int sum = 0;
	if (n % 2 == 0) {
		for (; sum <= k - 1; j++)
			sum += h[j];
		i = j;
	}

	for (; sum <= k; i++)
		sum += h[i];

	free(h);
	return (n % 2 == 0) ? (double) (i + j - 2) / 2.0 : (double) (i - 1);
}

double histogram_median_float(float *a, size_t n, threading_type threading) {
	float median;
	int threads = check_threading(&threading);
	findMinMaxPercentile(a, n, 0.5f, &median, 0.5f, &median, threads);
	return median;
}

/**
 * Compares a and b like strcmp()
 * @param a a gconstpointer
 * @param b a gconstpointer
 * @return an integer less than, equal to, or greater than zero, if a is than b .
 */
gint strcompare(gconstpointer *a, gconstpointer *b) {
	gchar *collate_key1, *collate_key2;
	gint result;

	const gchar *s1 = (const gchar *) a;
	const gchar *s2 = (const gchar *) b;

	collate_key1 = g_utf8_collate_key_for_filename(s1, strlen(s1));
	collate_key2 = g_utf8_collate_key_for_filename(s2, strlen(s2));

	result = g_strcmp0(collate_key1, collate_key2);
	g_free(collate_key1);
	g_free(collate_key2);

	return result;
}

/*
 * Optimal sorting network array of size [2,9] to retrieve median value
 * (C) Emmanuel Brandt 2019-02 (modified)
 * @param a array of float to sort (warning in place sorting)
 * *param b array of float weights (warning in place sorting)
 * @param n size of the array to sort [2,9]
 * warning in-place sorting
 */

#define sw(i, j)                 \
    if (a[i] > a[j]) {           \
        register float t = a[i]; \
        a[i] = a[j];             \
        a[j] = t;                \
        t = b[i];                \
        b[i] = b[j];             \
        b[j] = t;                \
    }

double sortnet_weighted_median_float(float *a, float *b, size_t n) {
    size_t k = n / 2;

    switch (n) {
        case 1:
            return a[0];
            break;

        case 2:
            sw(0, 1);
            break;

        case 3:
            sw(0, 1); sw(1, 2); sw(0, 1);
            break;

        case 4:
            sw(0, 1); sw(2, 3); sw(0, 2); sw(1, 3);
			sw(1, 2);
            break;

        case 5:
            sw(0, 1); sw(2, 3); sw(1, 3); sw(2, 4);
            sw(0, 2); sw(1, 4); sw(1, 2); sw(3, 4);
            sw(2, 3);
            break;

        case 6:
            sw(0, 1); sw(2, 3); sw(4, 5); sw(0, 2);
			sw(3, 5); sw(1, 4); sw(0, 1); sw(2, 3);
            sw(4, 5); sw(1, 2); sw(3, 4); sw(2, 3);
            break;

        case 7:
            sw(1, 2); sw(3, 4); sw(5, 6); sw(0, 2);
            sw(4, 6); sw(3, 5); sw(2, 6); sw(1, 5);
            sw(0, 4); sw(2, 5); sw(0, 3); sw(2, 4);
            sw(1, 3); sw(0, 1); sw(2, 3); sw(4, 5);
            break;

        case 8:
            sw(0, 1); sw(2, 3); sw(4, 5); sw(6, 7);
            sw(0, 2); sw(1, 3); sw(4, 6); sw(5, 7);
            sw(1, 2); sw(5, 6); sw(0, 4); sw(1, 5);
            sw(2, 6); sw(3, 7); sw(2, 4); sw(3, 5);
            sw(1, 2); sw(3, 4); sw(5, 6);
            break;

        case 9:
            sw(1, 8); sw(2, 7); sw(3, 6); sw(4, 5);
            sw(1, 4); sw(5, 8); sw(0, 2); sw(6, 7);
            sw(2, 6); sw(7, 8); sw(0, 3); sw(4, 5);
            sw(0, 1); sw(3, 5); sw(6, 7); sw(2, 4);
            sw(1, 3); sw(5, 7); sw(4, 6); sw(1, 2);
            sw(3, 4); sw(5, 6); sw(7, 8); sw(2, 3);
            sw(4, 5);
            break;

        default:
            return 0.0;
            break;  // no sort
    }
    return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2.0 : a[k];
}

#undef sw

double weighted_median_float(float *a, float *b, size_t n) {
    // Use faster and robust sorting network for small size array
    if (n < 9)
        return sortnet_weighted_median_float(a, b, n); // Implement sortnet_weighted_median_float function

    float total_weight = 0.0;
    for (size_t i = 0; i < n; ++i)
        total_weight += b[i];

    float half_total_weight = total_weight / 2.0;
    size_t k = 0;
    float cumulative_weight = 0.0;

    size_t left = 0;
    size_t right = n - 1;
    float tmp;

    while (left < right) {
        size_t pindex = (left + right) / 2;
        float pivot = a[pindex];
        float pivot_weight = b[pindex];
        a[pindex] = a[right];
        a[right] = pivot;
        b[pindex] = b[right];
        b[right] = pivot_weight;

        for (size_t i = pindex = left; i < right; ++i) {
            if (a[i] < pivot) {
                tmp = a[pindex];
                a[pindex] = a[i];
                a[i] = tmp;
                tmp = b[pindex];
                b[pindex] = b[i];
                b[i] = tmp;
                pindex++;
            }
        }
        a[right] = a[pindex];
        a[pindex] = pivot;
        b[right] = b[pindex];
        b[pindex] = pivot_weight;

        if (pindex < k)
            left = pindex + 1;
        else
            right = pindex;
    }

    // Find the index k such that sum of weights from 0 to k-1 is less than half of total weight
    while (cumulative_weight < half_total_weight && k < n) {
        cumulative_weight += b[k];
        ++k;
    }

    float median;
    if (cumulative_weight == half_total_weight || (cumulative_weight < half_total_weight && k > 0))
        median = a[k - 1]; // Weighted median lies within the same element
    else
        median = a[k]; // Weighted median lies between two elements

    return median;
}

#define sw_WORD(i, j)             \
    if (a[i] > a[j]) {            \
        WORD t = a[i];            \
        a[i] = a[j];              \
        a[j] = t;                 \
        t = b[i];                 \
        b[i] = b[j];              \
        b[j] = t;                 \
    }

double sortnet_weighted_median_WORD(WORD *a, WORD *b, size_t n) {
    size_t k = n / 2;

    switch (n) {
        case 1:
            return a[0];
            break;

        case 2:
            sw_WORD(0, 1);
            break;

        case 3:
            sw_WORD(0, 1); sw_WORD(1, 2); sw_WORD(0, 1);
            break;

        case 4:
            sw_WORD(0, 1); sw_WORD(2, 3); sw_WORD(0, 2); sw_WORD(1, 3);
			sw_WORD(1, 2);
            break;

        case 5:
            sw_WORD(0, 1); sw_WORD(2, 3); sw_WORD(1, 3); sw_WORD(2, 4);
            sw_WORD(0, 2); sw_WORD(1, 4); sw_WORD(1, 2); sw_WORD(3, 4);
            sw_WORD(2, 3);
            break;

        case 6:
            sw_WORD(0, 1); sw_WORD(2, 3); sw_WORD(4, 5); sw_WORD(0, 2);
			sw_WORD(3, 5); sw_WORD(1, 4); sw_WORD(0, 1); sw_WORD(2, 3);
            sw_WORD(4, 5); sw_WORD(1, 2); sw_WORD(3, 4); sw_WORD(2, 3);
            break;

        case 7:
            sw_WORD(1, 2); sw_WORD(3, 4); sw_WORD(5, 6); sw_WORD(0, 2);
            sw_WORD(4, 6); sw_WORD(3, 5); sw_WORD(2, 6); sw_WORD(1, 5);
            sw_WORD(0, 4); sw_WORD(2, 5); sw_WORD(0, 3); sw_WORD(2, 4);
            sw_WORD(1, 3); sw_WORD(0, 1); sw_WORD(2, 3); sw_WORD(4, 5);
            break;

        case 8:
            sw_WORD(0, 1); sw_WORD(2, 3); sw_WORD(4, 5); sw_WORD(6, 7);
            sw_WORD(0, 2); sw_WORD(1, 3); sw_WORD(4, 6); sw_WORD(5, 7);
            sw_WORD(1, 2); sw_WORD(5, 6); sw_WORD(0, 4); sw_WORD(1, 5);
            sw_WORD(2, 6); sw_WORD(3, 7); sw_WORD(2, 4); sw_WORD(3, 5);
            sw_WORD(1, 2); sw_WORD(3, 4); sw_WORD(5, 6);
            break;

        case 9:
            sw_WORD(1, 8); sw_WORD(2, 7); sw_WORD(3, 6); sw_WORD(4, 5);
            sw_WORD(1, 4); sw_WORD(5, 8); sw_WORD(0, 2); sw_WORD(6, 7);
            sw_WORD(2, 6); sw_WORD(7, 8); sw_WORD(0, 3); sw_WORD(4, 5);
            sw_WORD(0, 1); sw_WORD(3, 5); sw_WORD(6, 7); sw_WORD(2, 4);
            sw_WORD(1, 3); sw_WORD(5, 7); sw_WORD(4, 6); sw_WORD(1, 2);
            sw_WORD(3, 4); sw_WORD(5, 6); sw_WORD(7, 8); sw_WORD(2, 3);
            sw_WORD(4, 5);
            break;

        default:
            return 0.0;
            break;  // no sort
    }
    return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2.0 : a[k];
}

#undef sw_WORD

double weighted_median_WORD(WORD *a, WORD *b, size_t n) {
    if (n < 9)
        return sortnet_weighted_median_WORD(a, b, n);

    uint32_t total_weight = 0;
    for (size_t i = 0; i < n; ++i)
        total_weight += b[i];

    uint32_t half_total_weight = total_weight / 2;
    size_t k = 0;
    uint32_t cumulative_weight = 0;

    size_t left = 0;
    size_t right = n - 1;
    WORD tmp;

    while (left < right) {
        size_t pindex = (left + right) / 2;
        WORD pivot = a[pindex];
        WORD pivot_weight = b[pindex];
        a[pindex] = a[right];
        a[right] = pivot;
        b[pindex] = b[right];
        b[right] = pivot_weight;

        for (size_t i = pindex = left; i < right; ++i) {
            if (a[i] < pivot) {
                tmp = a[pindex];
                a[pindex] = a[i];
                a[i] = tmp;
                tmp = b[pindex];
                b[pindex] = b[i];
                b[i] = tmp;
                pindex++;
            }
        }
        a[right] = a[pindex];
        a[pindex] = pivot;
        b[right] = b[pindex];
        b[pindex] = pivot_weight;

        if (pindex < k)
            left = pindex + 1;
        else
            right = pindex;
    }

    while (cumulative_weight < half_total_weight && k < n) {
        cumulative_weight += b[k];
        ++k;
    }

    double median;
    if (cumulative_weight == half_total_weight || (cumulative_weight < half_total_weight && k > 0))
        median = a[k - 1]; // Weighted median lies within the same element
    else
        median = a[k]; // Weighted median lies between two elements

    return median;
}
