/*
MIT License

Copyright (c) 2018 Jérémy Anger, Gabriele Facciolo, Mauricio Delbracio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

#include "image.hpp"

/// project the gradients by shearing + accumulation
template <typename T>
void projectImage(img_t<T>& projections, const img_t<T>& u_x, const img_t<T>& u_y,
                const std::vector<angle_t>& angleSet)
{
    assert(u_x.w == u_y.w);
    assert(u_x.h == u_y.h);

    int w = u_x.w;
    int h = u_x.h;
    int maxSize = w + h;

    projections.ensure_size(maxSize, angleSet.size());

    // Transpose both images outside the parallel region
    img_t<T> u_yt = u_y;
    u_yt.transpose();
    img_t<T> u_xt = u_x;
    u_xt.transpose();

#ifdef _OPENMP
    #pragma omp parallel
    {
#endif
        std::vector<T> accumulationLine(maxSize);
        std::vector<int> countLine(maxSize);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic)
#endif
        for (unsigned a = 0; a < angleSet.size(); a++) {
            // Reset the accumulators
            std::fill(accumulationLine.begin(), accumulationLine.end(), 0.);
            std::fill(countLine.begin(), countLine.end(), 0);

            double factor;
            bool horizontalShear;
            double cos = std::cos(angleSet[a].angle);
            double sin = std::sin(angleSet[a].angle);
            if (angleSet[a].angle >= -M_PI/4 && angleSet[a].angle <= M_PI/4) {
                factor = std::tan(angleSet[a].angle);
                horizontalShear = true;
            } else {
                factor = 1. / std::tan(angleSet[a].angle);
                horizontalShear = false;
            }

            if (horizontalShear) {
                int start = (maxSize - w - factor*h) / 2;
#ifdef _OPENMP
                #pragma omp simd
#endif
                for (int y = 0; y < h; y++) {
                    int offset = start + round(factor * y);
                    for (int x = 0; x < w; x++) {
                        accumulationLine[x + offset] += u_x(x, y) * cos + u_y(x, y) * sin;
                        countLine[x + offset]++;
                    }
                }
            } else {
                int start = (maxSize - h - factor*w) / 2;
#ifdef _OPENMP
                #pragma omp simd
#endif
                for (int x = 0; x < w; x++) {
                    int offset = start + round(factor * x);
                    for (int y = 0; y < h; y++) {
                        accumulationLine[y + offset] += u_xt(y, x) * cos + u_yt(y, x) * sin;
                        countLine[y + offset]++;
                    }
                }
            }

            // Replace values that didn't get any samples by NAN
#ifdef _OPENMP
            #pragma omp simd
#endif
            for (int i = 0; i < maxSize; i++) {
                if (!countLine[i])
                    accumulationLine[i] = NAN;
            }

            // Copy to the resulting image
            std::copy(accumulationLine.begin(), accumulationLine.end(), &projections(0, a));
        }
#ifdef _OPENMP
    }
#endif
}

/// project the intensity by shearing + accumulation
template <typename T>
void projectImage(img_t<T>& projections, const img_t<T>& u,
                const std::vector<angle_t>& angleSet)
{
    int w = u.w;
    int h = u.h;
    int maxSize = w + h;

    projections.ensure_size(maxSize, angleSet.size());

    // Transpose the image outside the parallel region
    img_t<T> ut = u;
    ut.transpose();

#ifdef _OPENMP
    #pragma omp parallel
    {
#endif
        std::vector<T> accumulationLine(maxSize);
        std::vector<int> countLine(maxSize);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic)
#endif
        for (unsigned a = 0; a < angleSet.size(); a++) {
            // Reset the accumulators
            std::fill(accumulationLine.begin(), accumulationLine.end(), 0.);
            std::fill(countLine.begin(), countLine.end(), 0);

            double factor;
            bool horizontalShear;
            if (angleSet[a].angle >= -M_PI/4 && angleSet[a].angle <= M_PI/4) {
                factor = std::tan(angleSet[a].angle);
                horizontalShear = true;
            } else {
                factor = 1. / std::tan(angleSet[a].angle);
                horizontalShear = false;
            }

            if (horizontalShear) {
                int start = (maxSize - w - factor*h) / 2;
#ifdef _OPENMP
                #pragma omp simd
#endif
                for (int y = 0; y < h; y++) {
                    int offset = start + round(factor * y);
                    for (int x = 0; x < w; x++) {
                        accumulationLine[x + offset] += u(x, y);
                        countLine[x + offset]++;
                    }
                }
            } else {
                int start = (maxSize - h - factor*w) / 2;
#ifdef _OPENMP
                #pragma omp simd
#endif
                for (int x = 0; x < w; x++) {
                    int offset = start + round(factor * x);
                    for (int y = 0; y < h; y++) {
                        accumulationLine[y + offset] += ut(y, x);
                        countLine[y + offset]++;
                    }
                }
            }

            // Replace values that didn't get any samples by NAN
#ifdef _OPENMP
            #pragma omp simd
#endif
for (int i = 0; i < maxSize; i++) {
                if (!countLine[i])
                    accumulationLine[i] = NAN;
            }

            // Copy to the resulting image
            std::copy(accumulationLine.begin(), accumulationLine.end(), &projections(0, a));
        }
#ifdef _OPENMP
    }
#endif
}
