/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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
 *
 * This files contains all sorting algorithms for siril.
 * See src/test/sorting.c for testing and metrics.
 */
 
#include "sorting.h"
#include <string.h>
#include <math.h>

/*
 * In-place quick sort of array of double of size n
 * @param a array to sort
 * @param n array size
 */
void quicksort_d (double *a, int n) 
{
	if (n < 2) return;
	
	double pivot = a[n / 2];
	double *left = a;
	double *right = a + n - 1;
	register double t; 	// temporary value used for swap

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

/*
 * In-place quick sort of array of WORD of size n
 * @param a array of WORD
 * @param n array size
 */
void quicksort_s (WORD *a, int n) 
{
	if (n < 2) return;
	
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

/* quickmedian returns the median of an array of WORD of size n
 * Non recursive version of the original quickselect algorithm from
 * Tony Hoare modified to return the median value
 * Warning: values are shuffled in place
 * Written by Emmanuel Brandt 2019-02
 * @param a array of WORD
 * @param n array size
 * @return median as double for even size average the middle two elements
*/
double quickmedian (WORD *a, int n) 
{	// Use fast and robust sorting network for small size array
	if (n < 10) return sortnet_median (a, n);
	
	int i;
	int k = n / 2;		// size to sort
	int pindex;			// pivot index
	int left = 0; 		// left index
	int right = n - 1; 	// right index
	WORD pivot, tmp;

	while (left < right) { //we stop when our indicies have crossed
		pindex = (left + right) / 2; // pivot selection, this can be whatever
		pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP(pivot,right)

		for (i = pindex = left; i < right; i++) {
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
		else // pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ? ((double) a[k - 1] + (double) a[k]) / 2.0 : (double) a[k];
}

/* 
 * quickmedian_double returns the median of an array of double of size n
 * Non recursive version of the original quickselect algorithm from
 * Tony Hoare modified to return the median value
 * Warning: values are shuffled in place
 * Written by Emmanuel Brandt 2019-02
 * @param a array of double 
 * @param n array size 
 * @return median as double for even size average the middle two elements
 */
double quickmedian_double (double *a, int n) 
{
	int i;
	int k = n / 2;		// size to sort
	int pindex;			// pivot index
	int left = 0; 		// left index
	int right = n - 1; 	// right index
	double pivot, tmp;

	while (left < right) { //we stop when our indicies have crossed
		pindex = (left + right) / 2; // pivot selection, this can be whatever
		pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP(pivot,right)

		for (i = pindex = left; i < right; i++) {
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
		else // pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2 : a[k];
}

/*
 * quickmedian_int returns the median from an array of int of size n
 * Derived from original quickselect algorithm from Tony Hoare
 * non recurssive version modified to return the median value
 * Written by Emmanuel Brandt 2019-02
 * @param a array of int
 * @param n array size
 * @return median value as double
 */
double quickmedian_int (int *a, int n) {
	int i;
	int k = n / 2;		// size to sort
	int pindex;			// pivot index
	int left = 0; 		// left index
	int right = n - 1; 	// right index
	int pivot, tmp;

	while (left < right) { //we stop when our indicies have crossed
		pindex = (left + right) / 2; // pivot selection, this can be whatever
		pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP

		for (i = pindex = left; i < right; i++) {
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
		else // pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ?
			((double) a[k - 1] + (double) a[k]) / 2.0 : (double) a[k];
}

/*
 * sortnet_median in-place sort an array of size [2,9] using fast and robust 
 * sorting network and returns the median as double
 * Written by Emmanuel Brandt 2019-02
 * @param a array of WORD
 * @param n array size in range [2,9]
 * @return median value as double
 */
#define sw(i,j) if(a[i] > a[j]) { register WORD t=a[i]; a[i]=a[j]; a[j]=t; }
double sortnet_median (WORD *a, int n) 
{
	int k = n / 2;

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
	return (n % 2 == 0) ?
			((double) a[k - 1] + (double) a[k]) / 2.0 : (double) a[k];
}
#undef sw

/*
 * sortnet in-place sort of an array of WORD of size [2,9] using 
 * fast and robust sorting network.
 * Written by Emmanuel Brandt 2019-02
 * @param a array of WORD to sort
 * @param n array size in range [2,9]
 */
#define sw(i,j) if(a[i] > a[j]) { register WORD t=a[i]; a[i]=a[j]; a[j]=t; }
void sortnet (WORD *a, int n) 
{
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
 * histogram_median returns the median of large array of unsigned short
 * Use temp storage for histogram
 * Written by Emmanuel Brandt 2019-02
 * @param a array of unsigned short to search
 * @param n array size
 * @return median as a double
 */
double histogram_median (WORD *a, int n) 
{
	// For arrays size n < 10 use sortnet_median
	if (n < 10) return sortnet_median(a, n);

	unsigned int i, j, k = n / 2;
	size_t s = sizeof(unsigned int);
	unsigned int *h = (unsigned int *) calloc(USHRT_MAX + 1, s);

	for (i = 0; i < n; i++)
		h[a[i]]++;

	i = j = 0;
	if (n % 2 == 0) {
		for (; h[j] <= k - 1; j++)
			h[j + 1] += h[j];
		i = j;
	}

	for (; h[i] <= k; i++)
		h[i + 1] += h[i];

	free(h);
	return (n % 2 == 0) ? (double) (i + j) / 2.0 : (double) i;
}

/*
 * histogram_median_double return median of an array of double of size n
 * Developp for large array of values in range [0,1] 
 * Use temp storage for the histogram
 * It is faster than quickmedian_double for array > 20000
 * Written by Emmanuel Brandt 2019-02
 * @param a array of double
 * @param n array size
 * @return median as a double
 */
double histogram_median_double (double *a, int n) 
{
	unsigned int i, j, k = n / 2, nb_bins = 100000;
	size_t s = sizeof(unsigned int);
	unsigned int *h = (unsigned int *) calloc(nb_bins + 1, s);

	for (i = 0; i < n; i++)
		h[(unsigned int) (a[i] * nb_bins)]++;

	i = j = 0;
	if (n % 2 == 0) {
		for (; h[j] <= k - 1; j++)
			h[j + 1] += h[j];
		i = j;
	}

	for (; h[i] <= k; i++)
		h[i + 1] += h[i];

	free(h);
	return (n % 2 == 0) ? (double) (i + j) / ( 2.0 * (double) nb_bins) : (double) i / (double) nb_bins;
}

/*
 * maximum returns the maximum value of an array of WORD of size n
 * Written by Emmanuel Brandt 2019-03
 * @param a array of WORD
 * @param n array size
 * @return the maximum value as WORD
 */
WORD maximum (WORD *a, int n) 
{
	WORD m = a[0];
	for (int i=1; i<n; i++) 
		if (a[i]>m) m=a[i];
	
	return m;
}


// Same as maximum for double values
double maximum_d (double *a, int n) 
{
	double m = a[0];
	for (int i=1; i<n; i++) 
		if (a[i]>m) m=a[i];
	
	return m;
}
/*
 * Three Way Partition Quick Select returns the median of an array of WORD of size n
 * Combines the three way partition algorithm or "Dutch Flag" from 
 * Edsger Dijkstra and the Quick Select from Tony Hoare.
 * Splits the array in 3 partitions: < = > of pivot and iterate on partition 
 * containing the median. Outperforms QuickSelect algorithm for arrays with
 * duplicates or pre sorted
 * Warning: in place value shuffling
 * Written by Emmanuel Brandt 2019-04
 * @param a array of WORD
 * @param n array size
 * @return median as double
 */
double twpqs (WORD *a, int n) 
{
	// For arrays size n < 10 use sortnet_median
	if (n < 10) return sortnet_median(a, n);

	int left, l, j, right, r, k;
	WORD p, tmp;
	left=l=j=0;
	right=r=n-1;
	k=n/2; // pivot index
	p=a[k]; // pivot value
	
	while (left<=right) {
		while (j<=r) { // 3 way partition < = >
			if (a[j]<p) {
				tmp=a[j]; a[j]=a[l]; a[l]=tmp; // swap a[j] a[l]
				l+=1;
				j+=1;
			} else if (a[j]>p) {
				tmp=a[j]; a[j]=a[r]; a[r]=tmp; // swap a[j] a[r]
				r-=1;
			} else 
				j+=1;
		}

		if (k>r) {
			left=l=j=r+1;
			r=right;
		} else if (k<l) {
			r=right=l-1;
			l=j=left;	
		} else {
			if (n%2==0 && k<=l)
				return ((double) p + maximum (a,k)) / 2.0;
			else
				return (double) p;
		}
		p=a[(l+r)/2];
	}
}
// Same as twpqs for double values
double twpqs_d (double *a, int n) 
{
	int left, l, j, right, r, k=n/2;
	double p=a[k], tmp;
	
	left=l=j=0;
	right=r=n-1;
	
	while (left<=right) {
		while (j<=r) { // 3 way partition < = >
			if (a[j]<p) {
				tmp=a[j]; a[j]=a[l]; a[l]=tmp; // swap a[j] a[l]
				l+=1;
				j+=1;
			} else if (a[j]>p) {
				tmp=a[j]; a[j]=a[r]; a[r]=tmp; // swap a[j] a[r]
				r-=1;
			} else 
				j+=1;
		}

		if (k>r) {
			left=l=j=r+1;
			r=right;
		} else if (k<l) {
			r=right=l-1;
			l=j=left;	
		} else {
			if (n%2==0 && k<=l)
				return (p + maximum_d (a,k)) / 2.0;
			else
				return p;
		}
		
		p=a[(l+r)/2];
	}
}

/**
 * Compares a and b like strcmp()
 * @param a a gconstpointer
 * @param b a gconstpointer
 * @return an integer less than, equal to, or greater than zero, if a is than b .
 */
gint strcompare (gconstpointer *a, gconstpointer *b) 
{
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