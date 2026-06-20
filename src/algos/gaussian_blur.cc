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

/*
 * Separable Young-van Vliet recursive Gaussian blur (mono, float).
 *
 * The algorithm is derived from RawTherapee's src/rt/gauss.cc (Ingo Weyrich,
 * GPLv3+). The vertical and small-sigma passes there are already plain
 * #pragma omp simd loops and are reproduced here as-is. The original
 * horizontal pass (gaussHorizontalSse) was the codebase's last hand-written
 * SSE-intrinsic routine: it gained SIMD by stuffing four image rows into the
 * four lanes of an __m128. That trick is reproduced here as a portable
 * #pragma omp simd over an 8-row block, which the compiler widens to the
 * target's native width (AVX2 at the x86-64-v3 baseline, NEON on aarch64,
 * etc.) instead of being pinned to 128-bit SSE. No x86 intrinsics remain,
 * so this no longer needs src/rt's helpersse2.h / opthelper.h.
 */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstddef>
#include <vector>
#include <array>
#include <memory>
#include <new>
#include <algorithm>

#include "algos/gaussian_blur.h"

namespace {

template <typename T, std::size_t Alignment>
class AlignedAllocator {
public:
    using value_type = T;

    AlignedAllocator() noexcept = default;

    template <typename U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    T* allocate(std::size_t n) {
        if (n == 0)
            return nullptr;
        if (n > static_cast<std::size_t>(-1) / sizeof(T))
            throw std::bad_alloc();
        void* ptr = nullptr;
#ifdef _WIN32
        ptr = _aligned_malloc(n * sizeof(T), Alignment);
        if (!ptr)
            throw std::bad_alloc();
#else
        if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0)
            throw std::bad_alloc();
#endif
        return static_cast<T*>(ptr);
    }

    void deallocate(T* p, std::size_t) noexcept {
#ifdef _WIN32
        _aligned_free(p);
#else
        free(p);
#endif
    }

    template <typename U>
    struct rebind { using other = AlignedAllocator<U, Alignment>; };

    bool operator==(const AlignedAllocator&) const noexcept { return true; }
    bool operator!=(const AlignedAllocator& o) const noexcept { return !(*this == o); }
};

// Young-van Vliet / Triggs-Sdika 3rd-order recursive Gaussian coefficients.
void calculateYvVFactors(const double sigma, double &b1, double &b2, double &b3, double &B, double M[3][3])
{
    double q;
    if (sigma < 2.5)
        q = 3.97156 - 4.14554 * sqrt(1.0 - 0.26891 * sigma);
    else
        q = 0.98711 * sigma - 0.96330;

    double b0 = 1.57825 + 2.44413 * q + 1.4281 * q * q + 0.422205 * q * q * q;
    b1 = 2.44413 * q + 2.85619 * q * q + 1.26661 * q * q * q;
    b2 = -1.4281 * q * q - 1.26661 * q * q * q;
    b3 = 0.422205 * q * q * q;
    B = 1.0 - (b1 + b2 + b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    M[0][0] = -b3 * b1 + 1.0 - b3 * b3 - b2;
    M[0][1] = (b3 + b1) * (b2 + b3 * b1);
    M[0][2] = b3 * (b1 + b3 * b2);
    M[1][0] = b1 + b3 * b2;
    M[1][1] = -(b2 - 1.0) * (b2 + b3 * b1);
    M[1][2] = -(b3 * b1 + b3 * b3 + b2 - 1.0) * b3;
    M[2][0] = b3 * b1 + b2 + b1 * b1 - b2 * b2;
    M[2][1] = b1 * b2 + b3 * b2 * b2 - b1 * b3 * b3 - b3 * b3 * b3 - b3 * b2 + b3;
    M[2][2] = b3 * (b1 + b3 * b2);
}

// classical 3x3 filtering for very small sigma when src != dst
void gauss3x3(float** src, float** dst, const int W, const int H,
              const float c0, const float c1, const float c2, const float b0, const float b1)
{
#ifdef _OPENMP
    #pragma omp single nowait
#endif
    {
        dst[0][0] = src[0][0];
        for (int j = 1; j < W - 1; j++)
            dst[0][j] = b1 * (src[0][j - 1] + src[0][j + 1]) + b0 * src[0][j];
        dst[0][W - 1] = src[0][W - 1];
    }

#ifdef _OPENMP
    #pragma omp for nowait
#endif
    for (int i = 1; i < H - 1; i++) {
        dst[i][0] = b1 * (src[i - 1][0] + src[i + 1][0]) + b0 * src[i][0];
        for (int j = 1; j < W - 1; j++) {
            dst[i][j] = c2 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1])
                      + c1 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j])
                      + c0 * src[i][j];
        }
        dst[i][W - 1] = b1 * (src[i - 1][W - 1] + src[i + 1][W - 1]) + b0 * src[i][W - 1];
    }

#ifdef _OPENMP
    #pragma omp single
#endif
    {
        dst[H - 1][0] = src[H - 1][0];
        for (int j = 1; j < W - 1; j++)
            dst[H - 1][j] = b1 * (src[H - 1][j - 1] + src[H - 1][j + 1]) + b0 * src[H - 1][j];
        dst[H - 1][W - 1] = src[H - 1][W - 1];
    }
}

void gaussHorizontal3(float** src, float** dst, int W, int H, const float c0, const float c1)
{
    std::vector<float, AlignedAllocator<float, 16>> tempv(W);
    float *temp = tempv.data();
#ifdef _OPENMP
    #pragma omp for
#endif
    for (int i = 0; i < H; i++) {
        for (int j = 1; j < W - 1; j++)
            temp[j] = c1 * (src[i][j - 1] + src[i][j + 1]) + c0 * src[i][j];
        dst[i][0] = src[i][0];
        memcpy(dst[i] + 1, temp + 1, (W - 2) * sizeof(float));
        dst[i][W - 1] = src[i][W - 1];
    }
}

void gaussVertical3(float** src, float** dst, int W, int H, const float c0, const float c1)
{
    // process BLKSIZE columns at a time for cache efficiency; keep prev/curr/next
    // rows in local arrays so in-place (src==dst) is safe without aliasing
    constexpr int BLKSIZE = 8;
    std::vector<float, AlignedAllocator<float, 16>> prev(BLKSIZE), curr(BLKSIZE), next_row(BLKSIZE);

#ifdef _OPENMP
    #pragma omp for nowait
#endif
    for (int i = 0; i < W - BLKSIZE + 1; i += BLKSIZE) {
#pragma omp simd simdlen(8)
        for (int k = 0; k < BLKSIZE; k++) prev[k] = src[0][i + k];
#pragma omp simd simdlen(8)
        for (int k = 0; k < BLKSIZE; k++) dst[0][i + k] = prev[k];

        if (H > 1) {
#pragma omp simd simdlen(8)
            for (int k = 0; k < BLKSIZE; k++) curr[k] = src[1][i + k];

            for (int j = 1; j < H - 1; j++) {
#pragma omp simd simdlen(8)
                for (int k = 0; k < BLKSIZE; k++) next_row[k] = src[j + 1][i + k];
#pragma omp simd simdlen(8)
                for (int k = 0; k < BLKSIZE; k++)
                    dst[j][i + k] = c1 * (prev[k] + next_row[k]) + c0 * curr[k];
                std::swap(prev, curr);
                std::swap(curr, next_row);
            }
#pragma omp simd simdlen(8)
            for (int k = 0; k < BLKSIZE; k++) dst[H - 1][i + k] = src[H - 1][i + k];
        }
    }

    // remaining columns — buffer entire column before writing to handle src==dst
    std::vector<float, AlignedAllocator<float, 16>> colbuf(H);
#ifdef _OPENMP
    #pragma omp single
#endif
    for (int i = W - (W % BLKSIZE); i < W; i++) {
        for (int j = 0; j < H; j++) colbuf[j] = src[j][i];
        dst[0][i] = colbuf[0];
        for (int j = 1; j < H - 1; j++)
            dst[j][i] = c1 * (colbuf[j - 1] + colbuf[j + 1]) + c0 * colbuf[j];
        dst[H - 1][i] = colbuf[H - 1];
    }
}

// Horizontal recursive pass, vectorised across a block of rows (one row per
// SIMD lane). Portable reimplementation of the former gaussHorizontalSse:
// identical float arithmetic and Triggs-Sdika boundary, but expressed with
// #pragma omp simd so the width follows the target ISA.
void gaussHorizontalOmp(float** src, float** dst, const int W, const int H, const double sigma)
{
    double b1, b2, b3, B, M[3][3];
    calculateYvVFactors(sigma, b1, b2, b3, B, M);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            M[i][j] *= (1.0 + b2 + (b1 - b3) * b3);
            M[i][j] /= (1.0 + b1 - b2 + b3) * (1.0 - b1 - b2 - b3);
        }

    constexpr int NR = 8; // rows per SIMD block
    const float Bf = (float)B, b1f = (float)b1, b2f = (float)b2, b3f = (float)b3;
    float Mf[3][3];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            Mf[a][b] = (float)M[a][b];

    // tmp[column][row-lane]
    std::vector<std::array<float, NR>, AlignedAllocator<std::array<float, NR>, 64>> tmpv(W);
    auto *tmp = tmpv.data();

#ifdef _OPENMP
    #pragma omp for nowait
#endif
    for (int r = 0; r < H - NR + 1; r += NR) {
        // causal (forward) pass, one independent recursion per row-lane
#pragma omp simd simdlen(NR)
        for (int k = 0; k < NR; k++)
            tmp[0][k] = src[r + k][0] * (Bf + b1f + b2f + b3f);
#pragma omp simd simdlen(NR)
        for (int k = 0; k < NR; k++)
            tmp[1][k] = Bf * src[r + k][1] + b1f * tmp[0][k] + (b2f + b3f) * src[r + k][0];
#pragma omp simd simdlen(NR)
        for (int k = 0; k < NR; k++)
            tmp[2][k] = Bf * src[r + k][2] + b1f * tmp[1][k] + b2f * tmp[0][k] + b3f * src[r + k][0];

        for (int j = 3; j < W; j++) {
#pragma omp simd simdlen(NR)
            for (int k = 0; k < NR; k++)
                tmp[j][k] = Bf * src[r + k][j] + b1f * tmp[j - 1][k] + b2f * tmp[j - 2][k] + b3f * tmp[j - 3][k];
        }

        // Triggs-Sdika anti-causal boundary
        float tW[NR], tWp1[NR];
#pragma omp simd simdlen(NR)
        for (int k = 0; k < NR; k++) {
            const float sW = src[r + k][W - 1];
            const float u1 = tmp[W - 1][k] - sW;
            const float u2 = tmp[W - 2][k] - sW;
            const float u3 = tmp[W - 3][k] - sW;
            tW[k]   = sW + Mf[1][0] * u1 + Mf[1][1] * u2 + Mf[1][2] * u3;
            tWp1[k] = sW + Mf[2][0] * u1 + Mf[2][1] * u2 + Mf[2][2] * u3;
            tmp[W - 1][k] = sW + Mf[0][0] * u1 + Mf[0][1] * u2 + Mf[0][2] * u3;
        }
#pragma omp simd simdlen(NR)
        for (int k = 0; k < NR; k++) {
            tmp[W - 2][k] = Bf * tmp[W - 2][k] + b1f * tmp[W - 1][k] + b2f * tW[k] + b3f * tWp1[k];
            tmp[W - 3][k] = Bf * tmp[W - 3][k] + b1f * tmp[W - 2][k] + b2f * tmp[W - 1][k] + b3f * tW[k];
        }

        // anti-causal (backward) pass
        for (int j = W - 4; j >= 0; j--) {
#pragma omp simd simdlen(NR)
            for (int k = 0; k < NR; k++)
                tmp[j][k] = Bf * tmp[j][k] + b1f * tmp[j + 1][k] + b2f * tmp[j + 2][k] + b3f * tmp[j + 3][k];
        }

        // write back, one contiguous destination row per lane (safe in-place:
        // each row was fully read into tmp before any store)
        for (int k = 0; k < NR; k++) {
            float *drow = dst[r + k];
            for (int j = 0; j < W; j++)
                drow[j] = tmp[j][k];
        }
    }

    // remaining < NR rows, scalar in double precision (matches the original
    // gaussHorizontalSse border path)
#ifdef _OPENMP
    #pragma omp single
#endif
    {
        std::vector<double> tv(W);
        double *t = tv.data();
        for (int i = H - (H % NR); i < H; i++) {
            t[0] = src[i][0] * (B + b1 + b2 + b3);
            t[1] = B * src[i][1] + b1 * t[0] + (b2 + b3) * src[i][0];
            t[2] = B * src[i][2] + b1 * t[1] + b2 * t[0] + b3 * src[i][0];
            for (int j = 3; j < W; j++)
                t[j] = B * src[i][j] + b1 * t[j - 1] + b2 * t[j - 2] + b3 * t[j - 3];

            const double sW = src[i][W - 1];
            const double wm1   = sW + M[0][0] * (t[W - 1] - sW) + M[0][1] * (t[W - 2] - sW) + M[0][2] * (t[W - 3] - sW);
            const double tWs   = sW + M[1][0] * (t[W - 1] - sW) + M[1][1] * (t[W - 2] - sW) + M[1][2] * (t[W - 3] - sW);
            const double tWp1s = sW + M[2][0] * (t[W - 1] - sW) + M[2][1] * (t[W - 2] - sW) + M[2][2] * (t[W - 3] - sW);
            t[W - 1] = wm1;
            t[W - 2] = B * t[W - 2] + b1 * t[W - 1] + b2 * tWs + b3 * tWp1s;
            t[W - 3] = B * t[W - 3] + b1 * t[W - 2] + b2 * t[W - 1] + b3 * tWs;
            for (int j = W - 4; j >= 0; j--)
                t[j] = B * t[j] + b1 * t[j + 1] + b2 * t[j + 2] + b3 * t[j + 3];

            for (int j = 0; j < W; j++)
                dst[i][j] = (float)t[j];
        }
    }
}

// scalar horizontal pass (double precision) for large sigma; matches RT
void gaussHorizontal(float** src, float** dst, const int W, const int H, const double sigma)
{
    double b1, b2, b3, B, M[3][3];
    calculateYvVFactors(sigma, b1, b2, b3, B, M);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            M[i][j] /= (1.0 + b1 - b2 + b3) * (1.0 + b2 + (b1 - b3) * b3);

    std::vector<double, AlignedAllocator<double, 16>> tempv(W);
    double *temp2 = tempv.data();

#ifdef _OPENMP
    #pragma omp for
#endif
    for (int i = 0; i < H; i++) {
        temp2[0] = B * src[i][0] + b1 * src[i][0] + b2 * src[i][0] + b3 * src[i][0];
        temp2[1] = B * src[i][1] + b1 * temp2[0]  + b2 * src[i][0] + b3 * src[i][0];
        temp2[2] = B * src[i][2] + b1 * temp2[1]  + b2 * temp2[0]  + b3 * src[i][0];

        for (int j = 3; j < W; j++)
            temp2[j] = B * src[i][j] + b1 * temp2[j - 1] + b2 * temp2[j - 2] + b3 * temp2[j - 3];

        double temp2Wm1 = src[i][W - 1] + M[0][0] * (temp2[W - 1] - src[i][W - 1]) + M[0][1] * (temp2[W - 2] - src[i][W - 1]) + M[0][2] * (temp2[W - 3] - src[i][W - 1]);
        double temp2W   = src[i][W - 1] + M[1][0] * (temp2[W - 1] - src[i][W - 1]) + M[1][1] * (temp2[W - 2] - src[i][W - 1]) + M[1][2] * (temp2[W - 3] - src[i][W - 1]);
        double temp2Wp1 = src[i][W - 1] + M[2][0] * (temp2[W - 1] - src[i][W - 1]) + M[2][1] * (temp2[W - 2] - src[i][W - 1]) + M[2][2] * (temp2[W - 3] - src[i][W - 1]);

        temp2[W - 1] = temp2Wm1;
        temp2[W - 2] = B * temp2[W - 2] + b1 * temp2[W - 1] + b2 * temp2W + b3 * temp2Wp1;
        temp2[W - 3] = B * temp2[W - 3] + b1 * temp2[W - 2] + b2 * temp2[W - 1] + b3 * temp2W;

        for (int j = W - 4; j >= 0; j--)
            temp2[j] = B * temp2[j] + b1 * temp2[j + 1] + b2 * temp2[j + 2] + b3 * temp2[j + 3];

        for (int j = 0; j < W; j++)
            dst[i][j] = (float)temp2[j];
    }
}

// vertical recursive pass (double precision), 8 columns per SIMD lane group
void gaussVertical(float** src, float** dst, const int W, const int H, const double sigma)
{
    double b1, b2, b3, B, M[3][3];
    calculateYvVFactors(sigma, b1, b2, b3, B, M);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            M[i][j] /= (1.0 + b1 - b2 + b3) * (1.0 + b2 + (b1 - b3) * b3);

    constexpr int numcols = 8;
    auto temp2_storage = std::make_unique<std::vector<double, AlignedAllocator<double, 16>>>(H * numcols);
    std::unique_ptr<double*[]> temp2(new double*[H]);
    for (int i = 0; i < H; ++i)
        temp2[i] = &(*temp2_storage)[i * numcols];

    std::vector<double> temp2Hm1v(numcols), temp2Hv(numcols), temp2Hp1v(numcols);
    double *temp2Hm1 = temp2Hm1v.data(), *temp2H = temp2Hv.data(), *temp2Hp1 = temp2Hp1v.data();

#ifdef _OPENMP
    #pragma omp for nowait
#endif
    for (unsigned int i = 0; i < static_cast<unsigned>(std::max(0, W - numcols + 1)); i += numcols) {
#pragma omp simd simdlen(8)
        for (int k = 0; k < numcols; k++) {
            temp2[0][k] = B * src[0][i + k] + b1 * src[0][i + k] + b2 * src[0][i + k] + b3 * src[0][i + k];
            temp2[1][k] = B * src[1][i + k] + b1 * temp2[0][k] + b2 * src[0][i + k] + b3 * src[0][i + k];
            temp2[2][k] = B * src[2][i + k] + b1 * temp2[1][k] + b2 * temp2[0][k] + b3 * src[0][i + k];
        }

        for (int j = 3; j < H; j++) {
#pragma omp simd simdlen(8)
            for (int k = 0; k < numcols; k++)
                temp2[j][k] = B * src[j][i + k] + b1 * temp2[j - 1][k] + b2 * temp2[j - 2][k] + b3 * temp2[j - 3][k];
        }

#pragma omp simd simdlen(8)
        for (int k = 0; k < numcols; k++) {
            temp2Hm1[k] = src[H - 1][i + k] + M[0][0] * (temp2[H - 1][k] - src[H - 1][i + k]) + M[0][1] * (temp2[H - 2][k] - src[H - 1][i + k]) + M[0][2] * (temp2[H - 3][k] - src[H - 1][i + k]);
            temp2H[k]   = src[H - 1][i + k] + M[1][0] * (temp2[H - 1][k] - src[H - 1][i + k]) + M[1][1] * (temp2[H - 2][k] - src[H - 1][i + k]) + M[1][2] * (temp2[H - 3][k] - src[H - 1][i + k]);
            temp2Hp1[k] = src[H - 1][i + k] + M[2][0] * (temp2[H - 1][k] - src[H - 1][i + k]) + M[2][1] * (temp2[H - 2][k] - src[H - 1][i + k]) + M[2][2] * (temp2[H - 3][k] - src[H - 1][i + k]);
        }

#pragma omp simd simdlen(8)
        for (int k = 0; k < numcols; k++) {
            dst[H - 1][i + k] = temp2[H - 1][k] = temp2Hm1[k];
            dst[H - 2][i + k] = temp2[H - 2][k] = B * temp2[H - 2][k] + b1 * temp2[H - 1][k] + b2 * temp2H[k] + b3 * temp2Hp1[k];
            dst[H - 3][i + k] = temp2[H - 3][k] = B * temp2[H - 3][k] + b1 * temp2[H - 2][k] + b2 * temp2[H - 1][k] + b3 * temp2H[k];
        }

        for (int j = H - 4; j >= 0; j--) {
#pragma omp simd simdlen(8)
            for (int k = 0; k < numcols; k++)
                dst[j][i + k] = temp2[j][k] = B * temp2[j][k] + b1 * temp2[j + 1][k] + b2 * temp2[j + 2][k] + b3 * temp2[j + 3][k];
        }
    }

#ifdef _OPENMP
    #pragma omp single
#endif
    for (int i = W - (W % numcols); i < W; i++) {
        temp2[0][0] = B * src[0][i] + b1 * src[0][i] + b2 * src[0][i] + b3 * src[0][i];
        temp2[1][0] = B * src[1][i] + b1 * temp2[0][0]  + b2 * src[0][i] + b3 * src[0][i];
        temp2[2][0] = B * src[2][i] + b1 * temp2[1][0]  + b2 * temp2[0][0]  + b3 * src[0][i];

        for (int j = 3; j < H; j++)
            temp2[j][0] = B * src[j][i] + b1 * temp2[j - 1][0] + b2 * temp2[j - 2][0] + b3 * temp2[j - 3][0];

        double temp2Hm1s = src[H - 1][i] + M[0][0] * (temp2[H - 1][0] - src[H - 1][i]) + M[0][1] * (temp2[H - 2][0] - src[H - 1][i]) + M[0][2] * (temp2[H - 3][0] - src[H - 1][i]);
        double temp2Hs   = src[H - 1][i] + M[1][0] * (temp2[H - 1][0] - src[H - 1][i]) + M[1][1] * (temp2[H - 2][0] - src[H - 1][i]) + M[1][2] * (temp2[H - 3][0] - src[H - 1][i]);
        double temp2Hp1s = src[H - 1][i] + M[2][0] * (temp2[H - 1][0] - src[H - 1][i]) + M[2][1] * (temp2[H - 2][0] - src[H - 1][i]) + M[2][2] * (temp2[H - 3][0] - src[H - 1][i]);

        dst[H - 1][i] = temp2[H - 1][0] = temp2Hm1s;
        dst[H - 2][i] = temp2[H - 2][0] = B * temp2[H - 2][0] + b1 * temp2[H - 1][0] + b2 * temp2Hs + b3 * temp2Hp1s;
        dst[H - 3][i] = temp2[H - 3][0] = B * temp2[H - 3][0] + b1 * temp2[H - 2][0] + b2 * temp2[H - 1][0] + b3 * temp2Hs;

        for (int j = H - 4; j >= 0; j--)
            dst[j][i] = temp2[j][0] = B * temp2[j][0] + b1 * temp2[j + 1][0] + b2 * temp2[j + 2][0] + b3 * temp2[j + 3][0];
    }
}

void gaussianBlurImpl(float** src, float** dst, const int W, const int H, const double sigma)
{
    static constexpr double GAUSS_SKIP = 0.25;
    static constexpr double GAUSS_3X3_LIMIT = 0.6;
    static constexpr double GAUSS_DOUBLE = 25.0;

    if (sigma < GAUSS_SKIP) {
#ifdef _OPENMP
        #pragma omp single
#endif
        if (src != dst) {
            for (int i = 0; i < H; ++i)
                memcpy(dst[i], src[i], W * sizeof(float));
        }
    } else if (sigma < GAUSS_3X3_LIMIT) {
        if (src != dst) {
            double c0 = 1.0;
            double c1 = exp(-0.5 * (1.0 / (sigma * sigma)));
            double c2 = exp(-(1.0 / (sigma * sigma)));
            double sum = c0 + 4.0 * (c1 + c2);
            c0 /= sum; c1 /= sum; c2 /= sum;
            double b1 = exp(-1.0 / (2.0 * sigma * sigma));
            double bsum = 2.0 * b1 + 1.0;
            b1 /= bsum;
            double b0 = 1.0 / bsum;
            gauss3x3(src, dst, W, H, c0, c1, c2, b0, b1);
        } else {
            double c1 = exp(-1.0 / (2.0 * sigma * sigma));
            double csum = 2.0 * c1 + 1.0;
            c1 /= csum;
            double c0 = 1.0 / csum;
            gaussHorizontal3(src, dst, W, H, c0, c1);
            gaussVertical3(dst, dst, W, H, c0, c1);
        }
    } else {
        if (sigma < GAUSS_DOUBLE) {
            gaussHorizontalOmp(src, dst, W, H, sigma);
            gaussVertical(dst, dst, W, H, sigma);
        } else {
            gaussHorizontal(src, dst, W, H, sigma);
            gaussVertical(dst, dst, W, H, sigma);
        }
    }
}

} // namespace

extern "C" void gaussian_blur_mono(float** src, float** dst, const int W, const int H, const double sigma, int threads)
{
#ifdef _OPENMP
    #pragma omp parallel num_threads(threads)
#endif
    gaussianBlurImpl(src, dst, W, H, sigma);
}
