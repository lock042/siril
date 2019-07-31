#ifndef SORTING_H_
#define SORTING_H_

#include "core/siril.h"	// for types

/* the quicksorts */
void quicksort_d (double *a, int n);
void quicksort_s (WORD *a, int n);

/* fast median for array of 3 elements */
WORD minmax_med3 (WORD *a);

/* median for insertion sort */
WORD insertion_sort_median (WORD *a, int n);

/* Quick median based on quick select */
double quickmedian (WORD *a, int n);
double quickmedian_double(double *a, int n);
double quickmedian_int (int *a, int n);

/* Search for maximum value */
WORD maximum (WORD *a, int n);
double maximum_d (double *a, int n);

/* Quick select based on Three Way Partition */
double twpqs (WORD *a, int n);
double twpqs_d (double *a, int n);

/* Histogram median for very large array of unsigned short */
double histogram_median (WORD *a, int n);
double histogram_median_double (double *a, int n);

/* Sorting netnork */
double sortnet_median (WORD *a, int n);
void sortnet (WORD *a, int n);

gint strcompare(gconstpointer *a, gconstpointer *b);

#endif
