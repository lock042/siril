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

#include <criterion/criterion.h>
#include "../core/siril.h"
#include "../algos/sorting.h"

#include <stdlib.h>
#include <time.h>

/* This program tests the new implementation of the quickselect from Emmanuel
 * and the old from statistics.
 * It compares the results with the quicksort*/

#define NBTRIES 200	// for result checking, unit test of implementations

cominfo com;	// the core data struct
guiinfo gui;	// the gui data struct
fits gfit;	// currently loaded image

double median_from_sorted_array(WORD *arr, int size)
{
	if (size % 2)
		return arr[(size-1)/2];
	int sum = (int)arr[(size-1)/2] + (int)arr[size/2];
	return (double)sum/2.0;
}

int compare_median_algos(int datasize, int threads)
{
	WORD *data, *data_backup;
	double result_qsel1, result_qsel2, result_qsort;
	int i, retval = 0;

	data = malloc(datasize * sizeof(WORD));
	data_backup = malloc(datasize * sizeof(WORD));
	for (i=0; i<datasize; i++) {
		int val = rand() % USHRT_MAX;
		data[i] = (WORD)val;
		data_backup[i] = (WORD)val;
	}

	quicksort_s(data, datasize);
	result_qsort = median_from_sorted_array(data, datasize);
	memcpy(data_backup, data, datasize * sizeof(WORD));

	result_qsel1 = quickmedian(data, datasize);
	memcpy(data_backup, data, datasize * sizeof(WORD));

	result_qsel2 = histogram_median(data, datasize, threads);
	memcpy(data_backup, data, datasize * sizeof(WORD));

	if (result_qsel1 != result_qsort || result_qsel2 != result_qsort) {
		cr_log_error("got %g (quickmedian), %g (histogram_median) and %g (qsort)\n",
					 result_qsel1, result_qsel2, result_qsort);
		retval = 1;
	}

	free(data);
	free(data_backup);
	return retval;
}

void common_setup()
{
	srand(time(NULL));
	com.max_thread = g_get_num_processors();
}

TestSuite(Sorting, .init=common_setup);

Test(Sorting, Median)
{
	int size = 1;
	for (int i = 0; i < NBTRIES; i++, size++) {
		cr_assert(compare_median_algos(size, 1) == 0, "Failed at size=%u", size);
	}

	for (int i = 0; i < NBTRIES; i++, size++) {
		cr_assert(compare_median_algos(size, 2) == 0, "Failed at size=%u", size);
	}
}
