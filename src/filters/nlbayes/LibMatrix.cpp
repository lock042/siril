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
#include <cmath>
#include <cstdlib>
#include <vector>

template <typename T>
int inverseMatrix(T* io_mat, const unsigned p_N) {
    // Cholesky-based symmetric positive-definite inverse (ccmath psinv).
    // The pointer-walk inner loops are kept as-is (canonicalising them for
    // omp simd is high risk for low gain: p_N is small and the factorisation
    // is sequential). Accumulations are widened to double for precision.
    unsigned p, q, r, s, t, j, k;

    for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
        double diag = io_mat[p];
        for (q = j * p_N; q < p; q++)
            diag -= static_cast<double>(io_mat[q]) * io_mat[q];

        if (diag <= 0.0)
            return EXIT_FAILURE;

        io_mat[p] = static_cast<T>(std::sqrt(diag));

        for (k = j + 1, q = p + p_N; k < p_N; k++, q += p_N) {
            double z = 0.0;
            for (r = j * p_N, s = k * p_N; r < p; r++, s++)
                z += static_cast<double>(io_mat[r]) * io_mat[s];

            io_mat[q] = static_cast<T>((io_mat[q] - z) / io_mat[p]);
        }
    }

    transposeMatrix(io_mat, p_N);

    for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
        io_mat[p] = static_cast<T>(1.0 / io_mat[p]);

        for (q = j, t = 0; q < p; t += p_N + 1, q += p_N) {
            double z = 0.0;
            for (s = q, r = t; s < p; s += p_N, r++)
                z -= static_cast<double>(io_mat[s]) * io_mat[r];

            io_mat[q] = static_cast<T>(z * io_mat[p]);
        }
    }

    for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
        for (q = j, t = p - j; q <= p; q += p_N, t++) {
            double z = 0.0;
            for (k = j, r = p, s = q; k < p_N; k++, r++, s++)
                z += static_cast<double>(io_mat[r]) * io_mat[s];

            io_mat[t] = io_mat[q] = static_cast<T>(z);
        }
    }

    return EXIT_SUCCESS;
}

template <typename T>
void transposeMatrix(T* io_mat, const unsigned p_N) {
    for (unsigned i = 0; i < p_N - 1; i++) {
        unsigned p = i * (p_N + 1) + 1;
        unsigned q = i * (p_N + 1) + p_N;

        for (unsigned j = 0; j < p_N - 1 - i; j++, p++, q += p_N) {
            const T s = io_mat[p];
            io_mat[p] = io_mat[q];
            io_mat[q] = s;
        }
    }
}

template <typename T>
void covarianceMatrix(const T* i_patches, T* o_covMat, const unsigned p_nb,
                      const unsigned p_N) {
    const double coefNorm = 1.0 / static_cast<double>(p_nb);

    for (unsigned i = 0; i < p_N; i++) {
        const T* pi = i_patches + i * p_nb;
        for (unsigned j = 0; j < i + 1; j++) {
            const T* pj = i_patches + j * p_nb;
            double val = 0.0;
#ifdef _OPENMP
#pragma omp simd reduction(+:val)
#endif
            for (unsigned k = 0; k < p_nb; k++)
                val += static_cast<double>(pi[k]) * pj[k];

            o_covMat[i * p_N + j] = static_cast<T>(val * coefNorm);
            o_covMat[i + j * p_N] = static_cast<T>(val * coefNorm);
        }
    }
}

template <typename T>
void productMatrix(T* o_mat, const T* i_A, const T* i_B, const unsigned p_n,
                   const unsigned p_m, const unsigned p_l) {
    std::vector<T> q0(p_m, T(0));

    for (unsigned i = 0; i < p_l; i++) {
        for (unsigned k = 0; k < p_m; k++)
            q0[k] = i_B[i + k * p_l];

        for (unsigned j = 0; j < p_n; j++) {
            const T* a = i_A + j * p_m;
            double z = 0.0;
#ifdef _OPENMP
#pragma omp simd reduction(+:z)
#endif
            for (unsigned k = 0; k < p_m; k++)
                z += static_cast<double>(a[k]) * q0[k];

            o_mat[i + j * p_l] = static_cast<T>(z);
        }
    }
}

// Explicit instantiations.
template int inverseMatrix<float>(float*, const unsigned);
template int inverseMatrix<double>(double*, const unsigned);
template void transposeMatrix<float>(float*, const unsigned);
template void transposeMatrix<double>(double*, const unsigned);
template void covarianceMatrix<float>(const float*, float*, const unsigned, const unsigned);
template void covarianceMatrix<double>(const double*, double*, const unsigned, const unsigned);
template void productMatrix<float>(float*, const float*, const float*, const unsigned, const unsigned, const unsigned);
template void productMatrix<double>(double*, const double*, const double*, const unsigned, const unsigned, const unsigned);
