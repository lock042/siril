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

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric> // For std::gcd in C++17
#ifdef _OPENMP
#include <omp.h>
#endif

#include "image.hpp"
#include "angleSet.hpp"

// compute the angle set that allows reaching each pixel in a square
// of size kernelSize*kernelSize and starting at position 0,0
static void computeProjectionHalfAngleSet(std::vector<angle_t>& angles, int kernelSize)
{
    std::vector<angle_t> localAngles;

    // Reserve approximate space to reduce reallocations
    localAngles.reserve((kernelSize + 1) * (kernelSize + 1) / 2);

#ifdef _OPENMP
    int available_threads = com.max_thread - omp_get_num_threads();
#pragma omp parallel num_threads(available_threads) if (available_threads > 1)
    {
#endif
        std::vector<angle_t> privateAngles;
#ifdef _OPENMP
#pragma omp for collapse(2)
#endif
        for (int x = 0; x <= kernelSize; ++x) {
            for (int y = 0; y <= kernelSize; ++y) {
                // If gcd(x, y) is not one, skip to avoid duplicate angles
                if (std::gcd(x, y) != 1)
                    continue;

                angle_t angle;
                angle.angle = std::atan2(static_cast<double>(y), static_cast<double>(x));
                angle.x = x;
                angle.y = y;
                privateAngles.push_back(angle);
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        localAngles.insert(localAngles.end(), privateAngles.begin(), privateAngles.end());
#ifdef _OPENMP
    }
#endif

    // Sort by angle in descending order
    std::sort(localAngles.begin(), localAngles.end(), [](const angle_t& a, const angle_t& b) {
        return a.angle > b.angle;
    });

    angles.swap(localAngles);
}

// same as computeProjectionHalfAngleSet but with additional mirrored angles
void computeProjectionAngleSet(std::vector<angle_t>& angles, int kernelSize)
{
    // Get angles from pi/2 to 0
    computeProjectionHalfAngleSet(angles, kernelSize);

    // Calculate the new size of the angles vector
    int s = angles.size();
    angles.reserve(s * 2 - 2); // Reserve space to avoid reallocations

    // Add mirrored angles
    std::vector<angle_t> mirroredAngles(s - 2);

#ifdef _OPENMP
    int available_threads = com.max_thread - omp_get_num_threads();
#pragma omp parallel for schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
    for (int i = 1; i < s - 1; ++i) { // Start from 1 to skip theta = pi/2 and go to s-1 to skip theta = 0
        angle_t ref = angles[s - 1 - i];
        mirroredAngles[i - 1].angle = -ref.angle;
        mirroredAngles[i - 1].x = ref.x;
        mirroredAngles[i - 1].y = -ref.y;
    }

    angles.insert(angles.end(), mirroredAngles.begin(), mirroredAngles.end());
}
