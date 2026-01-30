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

/* The goal of this file is to test the pixel rejection functions. They cannot
 * be called directly from the main function apply_rejection_float() because
 * there are too many things to declare and types to include to do so, so we
 * copy as closely as possible what is done in this function here.
 */

#include <criterion/criterion.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_cdf.h>

#include "core/siril.h"
#include "stacking/rejection_float.c"

cominfo com;	// the core data struct
guiinfo gui;	// the gui data struct
fits gfit;	// currently loaded image

float set1[] = { 145, 125, 190, 135, 220, 130, 210, 3, 165, 165, 150, 350, 170, 180, 195, 440, 215, 135, 410, 40, 140, 175 };

float set2[] = { 7.7110e-2f, 4.7330e-1f, 5.7340e-1f, 3.3310e-1f, 5.3160e-1f, 3.6550e-1f, 3.1900e-1f, 3.4650e-1f, 2.2340e-1f, 5.3680e-1f, 4.8200e-1f, 4.8150e-1f, 2.5420e-1f, 7.3770e-1f, 6.6930e-1f, 3.8980e-1f, 5.8780e-1f, 6.6680e-1f, 6.9580e-1f, 3.6260e-1f, 7.1870e-1f, 2.6420e-1f, 5.2890e-1f, 6.1350e-1f, 2.4980e-1f, 2.7930e-1f, 7.9300e-1f, 6.6690e-1f, 5.9180e-1f, 6.5240e-1f, 8.4440e-2f, 8.1500e-1f, 3.5880e-1f, 3.7450e-1f, 5.6660e-1f, 2.5050e-1f, 5.6520e-1f, 4.6880e-1f, 9.7020e-2f, 4.9380e-1 };


static float *duplicate(const float *f, int size) {
	float *d = malloc(size * sizeof(float));
	memcpy(d, f, size * sizeof(float));
	return d;
}

static float compute_mean(const float *f, int size) {
	double sum = 0.0;
	for (int i = 0; i < size; i++)
		sum += f[i];
	return sum / size;
}

/*static void print_outliers(struct ESD_outliers *rej, int N) {
	fprintf(stderr, "outliers are: ");
	for (int i = 0; i < N; i++) {
		if (rej[i].out) fprintf(stderr, "%f ", rej[i].x);
	}
	fprintf(stderr, "\n");
}*/

static float calculate_critical_value(int size, float alpha) {
	float t_dist = gsl_cdf_tdist_Pinv(1 - alpha / (2 * size), size - 2);
	float numerator = (size - 1) * t_dist;
	float denominator = sqrtf(size) * sqrtf(size - 2 + (t_dist * t_dist));
	return numerator / denominator;
}

static float ESD_test(float *stack, int size, float alpha, int max_outliers) {
	struct ESD_outliers *out = malloc(max_outliers * sizeof(struct ESD_outliers));

	quicksort_f(stack, size);
	double median = gsl_stats_float_median_from_sorted_data(stack, 1, size);
	float *w_stack = duplicate(stack, size);
	int *rejected = calloc(size, sizeof(int));
	int w_size = size;
	int cold = 0;

	for (int iter = 0; iter < max_outliers; iter++) {
		float Gstat, Gcritical;
		int max_index = 0;

		Gcritical = calculate_critical_value(w_size, alpha);
		grubbs_stat(w_stack, w_size, &Gstat, &max_index);
		out[iter].out = check_G_values(Gstat, Gcritical);
		out[iter].x = w_stack[max_index];
		out[iter].i = (max_index == 0) ? cold++ : max_index;
		remove_element(w_stack, max_index, w_size);
		w_size--;
	}
	int count[2] = { 0, 0 };
	confirm_outliers(out, max_outliers, median, rejected, count);
	//print_outliers(out, max_outliers);

	cr_expect_eq(count[0], 2);
	cr_expect_eq(count[1], 3);

	cr_expect_float_eq(out[0].x, 440.0f, 1e-6);
	cr_expect_float_eq(out[1].x, 410.0f, 1e-6);
	cr_expect_float_eq(out[2].x, 350.0f, 1e-6);
	cr_expect_float_eq(out[3].x, 3.000f, 1e-6);
	cr_expect_float_eq(out[4].x, 40.00f, 1e-6);

	int kept = 0;
	double sum = 0.0;
	for (int frame = 0; frame < size; frame++) {
		if (!rejected[frame]) {
			sum += stack[frame];
			kept++;
		}
		//else printf("rejected %f\n", stack[frame]);
	}
	return sum / kept;
}

void test_GESDT_float() {
	float mean = ESD_test(set1, G_N_ELEMENTS(set1), 0.05, 7);
	cr_expect_float_eq(mean, 167.352936, 1e-6);
}

Test(rejection, GESDT) { test_GESDT_float(); }

/* NO_REJEC, PERCENTILE, SIGMA, MAD, SIGMEDIAN, WINSORIZED, LINEARFIT, GESDT */

/* test functions below take a stack and its size, populates the rej[2] and returns the mean */

static float percentile_test(float *stack, int size, float sig[2], int rej[2]) {
	float median = (float)quickmedian_float(stack, size);
	int kept = 0;
	double sum = 0.0;
	for (int frame = 0; frame < size; frame++) {
		if (!percentile_clipping(stack[frame], sig, median, rej)) {
			sum += stack[frame];
			kept++;
		}
		//else printf("rejected %f\n", stack[frame]);
	}
	return sum / kept;
}

static void test_percentile_float() {
	int rej[] = { 0, 0 };
	float sig[] = { 0.3, 0.4 };
	float mean = percentile_test(duplicate(set1, G_N_ELEMENTS(set1)), G_N_ELEMENTS(set1), sig, rej);
	//printf("percentile: %f\trej[0] = %d, rej[1] = %d\n", mean, rej[0], rej[1]);
	cr_expect_eq(rej[0], 2);
	cr_expect_eq(rej[1], 3);
	cr_expect_float_eq(mean, 167.352936f, 1e-6);

	rej[0] = 0; rej[1] = 0;
	sig[0] = 1.0; sig[1] = 1.0;
	mean = percentile_test(duplicate(set1, G_N_ELEMENTS(set1)), G_N_ELEMENTS(set1), sig, rej);
	//printf("percentile: %f\trej[0] = %d, rej[1] = %d\n", mean, rej[0], rej[1]);
	cr_expect_eq(rej[0], 0);
	cr_expect_eq(rej[1], 3);
	cr_expect_float_eq(mean, 152.0f, 1e-6);
}

Test(rejection, percentile) { test_percentile_float(); }


static float linearfit_test(float *stack, int size, float sig[2], int rej[2]) {
	float *xf = malloc(sizeof(float) * size);
	float m_x = (size - 1) * 0.5f;
	float m_dx2 = 0.f;
	for (int j = 0; j < size; ++j) {
		const float dx = j - m_x;
		xf[j] = 1.f / (j + 1);
		m_dx2 += (dx * dx - m_dx2) * xf[j];
	}
	m_dx2 = 1.f / m_dx2;
	int *rejected = calloc(size, sizeof(int));

	int changed, r = 0, N = size;
	do {
		quicksort_f(stack, N);

		float a, b;
		siril_fit_linear(xf, stack, m_x, m_dx2, N, &b, &a);
		float sigma = 0.f;
		for (int frame = 0; frame < N; frame++)
			sigma += fabsf(stack[frame] - (a * frame + b));
		sigma /= (float)N;
		//printf("sigma = %f\n", sigma);

		for (int frame = 0; frame < N; frame++) {
			if (N - r <= 4) {
				// no more rejections
				rejected[frame] = 0;
			} else {
				rejected[frame] = line_clipping(stack[frame], sig, sigma, frame, a, b, rej);
				if (rejected[frame] != 0) {
					r++;
					//printf("rejected %f\n", stack[frame]);
				}
			}
		}
		int output = 0;
		for (int pixel = 0; pixel < N; pixel++) {
			if (!rejected[pixel]) {
				// copy only if there was a rejection
				if (pixel != output)
					stack[output] = stack[pixel];
				output++;
			}
		}
		changed = N != output;
		N = output;
	} while (changed && N > 3);
	return compute_mean(stack, N);
}

static void test_linearfit_float() {
	int rej[] = { 0, 0 };
	float sig[] = { 2.5f, 2.5f };
	float mean = linearfit_test(duplicate(set2, G_N_ELEMENTS(set2)), G_N_ELEMENTS(set2), sig, rej);
	//printf("linearfit: %f\trej[0] = %d, rej[1] = %d\n", mean, rej[0], rej[1]);
	cr_expect_eq(rej[0], 3);
	cr_expect_eq(rej[1], 2);
	cr_expect_float_eq(mean, 0.476394f, 1e-6);

	rej[0] = 0; rej[1] = 0;
	sig[0] = 1.0f; sig[1] = 1.0f;
	mean = linearfit_test(duplicate(set2, G_N_ELEMENTS(set2)), G_N_ELEMENTS(set2), sig, rej);
	//printf("linearfit: %f\trej[0] = %d, rej[1] = %d\n", mean, rej[0], rej[1]);
	cr_expect_eq(rej[0], 7);
	cr_expect_eq(rej[1], 12);
	cr_expect_float_eq(mean, 0.4966f, 1e-5);
}

Test(rejection, linearfit) { test_linearfit_float(); }

