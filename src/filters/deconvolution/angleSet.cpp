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

#define dontneedcppfftwmultithreaded
#define dontneedcppfftwflags
#include "image.hpp"
#include "angleSet.hpp"
#undef dontneedcppfftwmultithreaded
#undef dontneedcppfftwflags


// greatest common factor
static int gcd(int a, int b)
{
    return b ? gcd(b, a % b) : a;
}
/*
int touch_the_vars() {
    return cppfftwmultithreaded * cppfftwflags;
}
*/
/// compute the angle set that allows to reach each pixels in a square
/// of size kernelSize*kernelSize and starting at position 0,0
static void computeProjectionHalfAngleSet(std::vector<angle_t>& angles, int kernelSize)
{
    for (int x = 0; x < kernelSize+1; x++) {
        for (int y = 0; y < kernelSize+1; y++) {
            // if gcd(x,y) is not one, this means that we could find (x',y')
            // such that (x,y) = (i*x',i*y'), thus we would have the same angle twice
            if (gcd(x, y) != 1)
                continue;

            angle_t angle;
            angle.angle = std::atan2(y, x);
            angle.x = x;
            angle.y = y;
            angles.push_back(angle);
        }
    }

    struct angle_sort_by_angle {
        bool operator() (const angle_t& a, const angle_t& b) {
            return a.angle > b.angle;
        }
    };
    std::sort(angles.begin(), angles.end(), angle_sort_by_angle());
}

/// same as computeProjectionHalfAngleSet but with additional mirrored angles
void computeProjectionAngleSet(std::vector<angle_t>& angles, int kernelSize)
{
    // get angles from pi/2 to 0
    computeProjectionHalfAngleSet(angles, kernelSize);

    // copy the angles in reverse order (so that the final set goes from pi/2 to -pi/2)
    int s = angles.size();

    // don't copy the first one and the last one (theta=pi/2 and theta=0)
    angles.resize(s * 2 - 2);

    for (int i = 0; i < s-2; i++) {
        angle_t ref = angles[s-1 - i - 1];
        angles[s + i] = ref;
        angles[s + i].angle = - ref.angle;
        angles[s + i].y = - ref.y;
    }
}

