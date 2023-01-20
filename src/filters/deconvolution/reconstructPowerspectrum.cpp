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
#include <complex>
#include <cassert>

#include "image.hpp"
#include "angleSet.hpp"

/// reconstruct the power spectrum from a set of autocorrelations of projections
/// each projection is used to reconstruct one or more coefficients
template <typename T>
void reconstructPowerspectrum(img_t<T>& powerSpectrum, const img_t<T> acProjections,
                              const std::vector<angle_t>& angleSet, int psSize)
{
    powerSpectrum.ensure_size(psSize*2+1, psSize*2+1);
    powerSpectrum.set_value(0.);

    img_t<std::complex<T>> ftAutocorrelation(acProjections.w, 1);
    img_t<T> powerSpectrumSlice(acProjections.w, 1);

    for (unsigned j = 0; j < angleSet.size(); j++) {
        // compute the discrete Fourier transform of the autocorrelation
        // (= power spectrum by the Wiener-Khinchin theorem)
        for (int x = 0; x < ftAutocorrelation.w; x++)
            ftAutocorrelation[x] = acProjections(x, j);
        ftAutocorrelation.fft(ftAutocorrelation);

        for (int x = 0; x < acProjections.w; x++) {
            powerSpectrumSlice[x] = std::abs(ftAutocorrelation[x]);
        }

        T normalize = powerSpectrumSlice[0];
        for (int x = 0; x < acProjections.w; x++) {
            powerSpectrumSlice[x] /= normalize;
        }

        powerSpectrumSlice.fftshift();

        // extract and place back the coefficient that intersect the grid
        for (int i = 1; i < psSize + 1; i++) {
            int xOffset = i * angleSet[j].x;
            int yOffset = i * angleSet[j].y;

            if (std::abs(xOffset) > psSize || std::abs(yOffset) > psSize)
                break;

            // place the sample in the 2D power spectrum
            int sliceOffset = std::max(std::abs(xOffset), std::abs(yOffset));
            powerSpectrum(psSize + xOffset, psSize + yOffset) = powerSpectrumSlice[psSize + sliceOffset];
            powerSpectrum(psSize - xOffset, psSize - yOffset) = powerSpectrumSlice[psSize + sliceOffset];
        }
    }

    // the DC value of the kernel is 1
    powerSpectrum(psSize, psSize) = 1.;
}

