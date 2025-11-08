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

#include "../core/siril.h"
#include "../algos/sorting.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define THREADING_TYPE MULTI_THREADED

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

double _siril_qsort(WORD *data, size_t datasize)
{
	quicksort_s(data, datasize);
	return median_from_sorted_array(data, datasize);
}

double _histogram_sort(WORD *data, size_t datasize)
{
	return histogram_median(data, datasize, THREADING_TYPE);
}

clock_t perf_test(double (*function)(WORD *data, size_t datasize),
				  int datasize,
				  int nb_draws,
				  int nb_times_each)
{
	WORD *data = malloc(datasize * sizeof(WORD));
	WORD *data_backup = malloc(datasize * sizeof(WORD));
	int i, draws, times;

	clock_t t_start = clock();
	for (draws = 0; draws < nb_draws; draws++) {
		for (i=0; i<datasize; i++) {
			int val = rand() % USHRT_MAX;
			data[i] = (WORD)val;
			data_backup[i] = (WORD)val;
		}

		for (times = 0; times < nb_times_each; times++) {
			function(data, datasize);

			memcpy(data, data_backup, datasize * sizeof(WORD));
			data[times % datasize] = times % USHRT_MAX;
		}
	}
	clock_t t_end = clock();

	free(data);
	free(data_backup);

	return t_end - t_start;
}

/******
 *
 * As clock_t has an unspecified type, we follow the advice on this link to output results:
 * https://stackoverflow.com/questions/1083142/what-s-the-correct-way-to-use-printf-to-print-a-clock-t
 *
 *
 *  */

void MeasureSmall()
{
	int datasize = 8;
	int nb_draws = 100;
	int nb_times_each = 200000;

	fprintf(stdout, "== small dataset (%d elements, %d different draws run %d times)\n",
			datasize, nb_draws, nb_times_each);

	clock_t t_siril = perf_test(_siril_qsort, datasize, nb_draws, nb_times_each);
	clock_t t_quick = perf_test(quickmedian, datasize, nb_draws, nb_times_each);
	clock_t t_hist = perf_test(_histogram_sort, datasize, nb_draws, nb_times_each);

	fprintf(stdout, "siril quicksort time:\t%.0Lf\n", (long double) t_siril);
	fprintf(stdout, "quickmedian time:\t%.0Lf\n", (long double) t_quick);
	fprintf(stdout, "histogram_median time:\t%.0Lf\n", (long double) t_hist);
}

void MeasureBig()
{
	int datasize = 30000000;

	fprintf(stdout, "== large dataset (%d elements, same for each)\n", datasize);

	clock_t t_siril = perf_test(_siril_qsort, datasize, 1, 1);
	clock_t t_quick = perf_test(quickmedian, datasize, 1, 1);
	clock_t t_hist = perf_test(_histogram_sort, datasize, 1, 1);

	fprintf(stdout, "siril quicksort time:\t%.0Lf\n", (long double) t_siril);
	fprintf(stdout, "quickmedian time:\t%.0Lf\n", (long double) t_quick);
	fprintf(stdout, "histogram_median time:\t%.0Lf\n", (long double) t_hist);
}

int main()
{
	srand(time(NULL));
	com.max_thread = g_get_num_processors();

	fputc('\n', stdout);
	MeasureBig();
	fputc('\n', stdout);
	MeasureSmall();
	fputc('\n', stdout);
	return 0;
}
