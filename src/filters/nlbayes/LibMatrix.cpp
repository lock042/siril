/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file LibMatrix.cpp
 * @brief Tools for matrix manipulation, based on ccmath functions
 *        by Daniel A. Atkinson.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

 #include "LibMatrix.h"
#include <iostream>
 #include <math.h>
 #include <stdlib.h>

 using namespace std;

/**
 * @brief Invert (in place) a symmetric real matrix, V -> Inv(V).
 *			The input matrix V is symmetric (V[i,j] = V[j,i]).
 *
 * @param io_mat : array containing a symmetric input matrix. This is converted to the inverse
			matrix;
 * @param p_N : dimension of the system (dim(v)=n*n)
 *
 * @return	EXIT_SUCCESS -> normal exit
 *			EXIT_FAILURE -> input matrix not positive definite
 **/
int inverseMatrix(
	vector<float> &io_mat
,	const unsigned p_N
){
	//! Initializations
	float z;
	unsigned p, q, r, s, t, j, k;

	for (j = 0, p = 0; j < p_N ; j++, p += p_N + 1) {
		for (q = j * p_N; q < p ; q++) {
			io_mat[p] -= io_mat[q] * io_mat[q];
		}

		if (io_mat[p] <= 0.f) {
			return EXIT_FAILURE;
		}

		io_mat[p] = sqrtf(io_mat[p]);

		for (k = j + 1, q = p + p_N; k < p_N ; k++, q += p_N) {
			for (r = j * p_N, s = k * p_N, z = 0.f; r < p; r++, s++) {
				z += io_mat[r] * io_mat[s];
			}

			io_mat[q] -= z;
			io_mat[q] /= io_mat[p];
		}
	}

	transposeMatrix(io_mat, p_N);

	for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
		io_mat[p] = 1.f / io_mat[p];

		for (q = j, t = 0; q < p; t += p_N + 1, q += p_N) {
			for (s = q, r = t, z = 0.f; s < p; s += p_N, r++) {
				z -= io_mat[s] * io_mat[r];
			}

			io_mat[q] = z * io_mat[p];
		}
	}

	for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
		for (q = j, t = p - j; q <= p; q += p_N, t++) {
			for (k = j, r = p, s = q, z = 0.f; k < p_N; k++, r++, s++) {
				z += io_mat[r] * io_mat[s];
			}

			io_mat[t] = io_mat[q] = z;
		}
	}

	return EXIT_SUCCESS;
}

/**
 * @brief Transpose a real square matrix in place io_mat -> io_mat~.
 *
 * @param io_mat : pointer to an array of p_N by p_N input matrix io_mat. This is overloaded by
			the transpose of io_mat;
 * @param p_N : dimension (dim(io_mat) = p_N * p_N).
 *
 * @return none.
 **/
void transposeMatrix(
	vector<float> &io_mat
,	const unsigned p_N
){
	for (unsigned i = 0; i < p_N - 1; i++) {
		unsigned p = i * (p_N + 1) + 1;
		unsigned q = i * (p_N + 1) + p_N;

		for (unsigned j = 0; j < p_N - 1 - i; j++, p++, q += p_N) {
			const float s = io_mat[p];
			io_mat[p] = io_mat[q];
			io_mat[q] = s;
		}
	}
}

/**
 * @brief Compute the covariance matrix.
 *
 * @param i_patches: set of patches of size (nb x N);
 * @param o_covMat: will contain patches' * patches;
 * @param p_N : size of a patch;
 * @param p_nb: number of similar patches in the set of patches.
 *
 * @return none.
 **/
void covarianceMatrix(
	vector<float> const& i_patches
,	vector<float> &o_covMat
,	const unsigned p_nb
,	const unsigned p_N
){
	const float coefNorm = 1.f / (float) (p_nb);

	for (unsigned i = 0; i < p_N; i++) {
		for (unsigned j = 0; j < i + 1; j++) {
			float val = 0.f;

			for (unsigned k = 0; k < p_nb; k++) {
				val += i_patches[i * p_nb + k] * i_patches[j * p_nb + k];
			}
			o_covMat[i * p_N + j] = val * coefNorm;
			o_covMat[i + j * p_N] = val * coefNorm;
		}
	}
}

/**
 * @brief Multiply two matrix A * B. (all matrices stored in row order).
 *
 * @param o_mat = array containing n by l product matrix at exit;
 * @param i_A = input array containing n by m matrix;
 * @param i_B = input array containing m by l matrix;
 * @param p_n, p_m, p_l = dimension parameters of arrays.
 *
 * @return  none.
 **/
void productMatrix(
	vector<float> &o_mat
,	vector<float> const& i_A
,	vector<float> const& i_B
,	const unsigned p_n
,	const unsigned p_m
,	const unsigned p_l
){
	vector<float> q0(p_m, 0.f);

	for (unsigned i = 0; i < p_l; i++) {
		for (unsigned k = 0; k < p_m; k++) {
			q0[k] = i_B[i + k * p_l];
		}

		for (unsigned j = 0; j < p_n; j++) {
			float z = 0.f;

			for (unsigned k = 0; k < p_m; k++) {
				z += i_A[j * p_m + k] * q0[k];
			}

			o_mat[i + j * p_l] = z;
		}
	}
}
