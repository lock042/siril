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
#ifndef LIB_MATRIX_H_INCLUDED
#define LIB_MATRIX_H_INCLUDED

/*
 * Patch-matrix linear algebra for NL-Bayes. Templated on the element type and
 * container-agnostic (raw pointer + dimension), so the same routines operate on
 * std::vector buffers and on img_t-owned buffers with no copy. (std::span would
 * be the natural type but the project is C++17.) Inner reductions accumulate in
 * double and, where the loop is canonical, use `omp simd`; thread-level
 * parallelism is intentionally NOT added because NL-Bayes already saturates
 * threads at the sub-image level, so these run nested.
 *
 * Explicit float and double instantiations are provided in LibMatrix.cpp.
 */

/**
 * @brief Invert (in place) a symmetric positive-definite matrix, V -> Inv(V).
 * @param io_mat : p_N x p_N symmetric matrix, overwritten with its inverse.
 * @param p_N : dimension.
 * @return EXIT_SUCCESS, or EXIT_FAILURE if not positive definite.
 */
template <typename T>
int inverseMatrix(T* io_mat, const unsigned p_N);

/**
 * @brief Transpose a p_N x p_N square matrix in place.
 */
template <typename T>
void transposeMatrix(T* io_mat, const unsigned p_N);

/**
 * @brief Compute the covariance matrix o_covMat = i_patches^T * i_patches / nb.
 * @param i_patches : p_N x p_nb matrix (row = patch pixel, col = patch).
 * @param o_covMat : p_N x p_N output.
 * @param p_nb : number of patches.
 * @param p_N : patch size.
 */
template <typename T>
void covarianceMatrix(const T* i_patches, T* o_covMat, const unsigned p_nb,
                      const unsigned p_N);

/**
 * @brief Multiply two row-major matrices: o_mat (n x l) = i_A (n x m) * i_B (m x l).
 */
template <typename T>
void productMatrix(T* o_mat, const T* i_A, const T* i_B, const unsigned p_n,
                   const unsigned p_m, const unsigned p_l);

#endif // LIB_MATRIX_H_INCLUDED
