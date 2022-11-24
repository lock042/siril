#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fftw3.h>
#include "fftw_allocator.hpp"
#include <functional>
#include "image.hpp"
#include <limits>
#include <numeric>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

struct options {
    std::string input;
    int kernelSize;
    std::string out_kernel;
    std::string out_deconv;

    int Ninner;
    int Ntries;
    int Nouter;
    float compensationFactor;
    int medianFilter;

    float finalDeconvolutionWeight;
    float intermediateDeconvolutionWeight;
    int seed;
};

struct angle_t {
    int x, y;
    double angle;
    double dist;
    bool valid;
};

void computeProjectionAngleSet(std::vector<angle_t>& angles, int kernelSize);

template <typename T>
void reconstructPowerspectrum(img_t<T>& powerSpectrum, const img_t<T> acProjections,
                              const std::vector<angle_t>& angleSet, int psSize);

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

template <typename T>
void projectImage(img_t<T>& projections, const img_t<T>& imgX, const img_t<T>& imgY,
                const std::vector<angle_t>& angleSet);

template <typename T>
void projectImage(img_t<T>& projections, const img_t<T>& img,
                const std::vector<angle_t>& angleSet);

/// project the gradients by shearing + accumulation
template <typename T>
void projectImage(img_t<T>& projections, const img_t<T>& u_x, const img_t<T>& u_y,
                const std::vector<angle_t>& angleSet)
{
    assert(u_x.w == u_y.w);
    assert(u_x.h == u_y.h);

    // transpose both images to speedup pixel access for vertical shear
    // it can lead to a measured speedup of 10x
    img_t<T> u_yt = u_y;
    u_yt.transpose();
    img_t<T> u_xt = u_x;
    u_xt.transpose();

    int w = u_x.w;
    int h = u_x.h;
    int maxSize = w + h;

    projections.ensure_size(maxSize, angleSet.size());

    // parallelize per orientation
#pragma omp parallel
    {
        std::vector<T> accumulationLine(maxSize);
        std::vector<int> countLine(maxSize);
#pragma omp for
        for (unsigned a = 0; a < angleSet.size(); a++) {
            // reset the accumulators
            for (int i = 0; i < maxSize; i++) {
                accumulationLine[i] = 0.;
                countLine[i] = 0;
            }

            // compute the shearing factor and orientation
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
                for (int y = 0; y < h; y++) {
                    int offset = start + round(factor * y);
                    for (int x = 0; x < w; x++) {
                        accumulationLine[x + offset] += u_x(x, y) * cos + u_y(x, y) * sin;
                        countLine[x + offset]++;
                    }
                }
            } else {
                int start = (maxSize - h - factor*w) / 2;
                for (int x = 0; x < w; x++) {
                    int offset = start + round(factor * x);
                    for (int y = 0; y < h; y++) {
                        // the next line is the original one,
                        //accumulationLine[y + offset] += u_x(x, y) * cos + u_y(x, y) * sin;
                        // and below is the one sped up by the transposition
                        accumulationLine[y + offset] += u_xt(y, x) * cos + u_yt(y, x) * sin;
                        countLine[y + offset]++;
                    }
                }
            }

            // replace values that didn't get any samples by NAN
            // we do so in order to extract the valid values for the autocorrelation
            for (int i = 0; i < maxSize; i++) {
                if (!countLine[i])
                    accumulationLine[i] = NAN;
            }

            // copy to the resulting image
            std::copy(accumulationLine.begin(), accumulationLine.end(), &projections(0, a));
        }
    }
}

/// project the intensity by shearing + accumulation
template <typename T>
void projectImage(img_t<T>& projections, const img_t<T>& u,
                const std::vector<angle_t>& angleSet)
{
    // transpose both images to speedup pixel access for vertical shear
    // it can lead to a measured speedup of 10x
    img_t<T> ut = u;
    ut.transpose();

    int w = u.w;
    int h = u.h;
    int maxSize = w + h;

    projections.ensure_size(maxSize, angleSet.size());

    // parallelize per orientation
#pragma omp parallel
    {
        std::vector<T> accumulationLine(maxSize);
        std::vector<int> countLine(maxSize);
#pragma omp for
        for (unsigned a = 0; a < angleSet.size(); a++) {
            // reset the accumulators
            for (int i = 0; i < maxSize; i++) {
                accumulationLine[i] = 0.;
                countLine[i] = 0;
            }

            // compute the shearing factor and orientation
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
                for (int y = 0; y < h; y++) {
                    int offset = start + round(factor * y);
                    for (int x = 0; x < w; x++) {
                        accumulationLine[x + offset] += u(x, y);
                        countLine[x + offset]++;
                    }
                }
            } else {
                int start = (maxSize - h - factor*w) / 2;
                for (int x = 0; x < w; x++) {
                    int offset = start + round(factor * x);
                    for (int y = 0; y < h; y++) {
                        // the next line is the original equation
                        //accumulationLine[y + offset] += u_x(x, y) * cos + u_y(x, y) * sin;
                        // and below is the one sped up by the transposition
                        accumulationLine[y + offset] += ut(y, x);
                        countLine[y + offset]++;
                    }
                }
            }

            // replace values that didn't get any samples by NAN
            // we do so in order to extract the valid values for the autocorrelation
            for (int i = 0; i < maxSize; i++) {
                if (!countLine[i])
                    accumulationLine[i] = NAN;
            }

            // copy to the resulting image
            std::copy(accumulationLine.begin(), accumulationLine.end(), &projections(0, a));
        }
    }
}

template <typename T>
void phaseRetrieval(img_t<T>& kernel, const img_t<T>& blurredPatch,
                    const img_t<T>& powerSpectrum, int kernelSize,
                    const options& opts);

/// Algorithm 6
template <typename T>
static void singlePhaseRetrieval(img_t<T>& kernel, const img_t<T>& magnitude,
                                 int kernelSize, int nbIterations)
{
    using complex = std::complex<T>;
    static const complex I(0, 1);

    kernel.ensure_size(kernelSize, kernelSize);
    kernel.set_value(0);

    const T alpha = 0.95;
    const T beta0 = 0.75;

    img_t<complex> ftkernel(magnitude.w, magnitude.h);
    img_t<T> g(magnitude.w, magnitude.h);
    img_t<complex> gft(g.w, g.h, g.d);
    img_t<T> g2(g);
    img_t<T> R(g);
    img_t<char> omega(g.w, g.h); // can't use bool because of std::vector

    for (int i = 0; i < ftkernel.size; i++) {
        T phase = ((T)rand()/RAND_MAX) * M_PI * 2 - M_PI;
        ftkernel[i] = magnitude[i] * std::exp(I * phase);
    }
    ftkernel.ifft(ftkernel);
    for (int i = 0; i < g.size; i++) {
        g[i] = std::real(ftkernel[i]);
    }

    for (int m = 0; m < nbIterations; m++) {
        T beta = beta0 + (T(1.) - beta0) * (T(1.) - std::exp(- std::pow(m / T(7.), T(3.))));

        for (int i = 0; i < g.size; i++) {
            gft[i] = g[i];
        }
        gft.fft(gft);

        for (int i = 0; i < gft.size; i++) {
            gft[i] = (alpha * magnitude[i] + (T(1.) - alpha) * std::abs(gft[i]))
                     * std::exp(I * std::arg(gft[i]));
        }

        gft.ifft(gft);
        for (int i = 0; i < g.size; i++) {
            g2[i] = std::real(gft[i]);
        }

        for (int i = 0; i < R.size; i++) {
            R[i] = T(2.) * g2[i] - g[i];
        }
        for (int i = 0; i < omega.size; i++) {
            omega[i] = R[i] < T(0.);
        }

        for (int y = 0; y < magnitude.h; y++)
        for (int x = kernelSize; x < magnitude.w; x++) {
            omega(x, y) = true;
        }
        for (int y = kernelSize; y < magnitude.h; y++)
        for (int x = 0; x < magnitude.w; x++) {
            omega(x, y) = true;
        }

        for (int i = 0; i < g.size; i++) {
            g[i] = omega[i] ? beta * g[i] + (T(1.) - T(2.)*beta) * g2[i] : g2[i];
        }
    }

    for (int y = 0; y < kernelSize; y++)
    for (int x = 0; x < kernelSize; x++)
        kernel(x, y) = g2(x, y) >= T(0.) ? g2(x, y) : T(0.);
    kernel.normalize();

    // apply the thresholding of 1/255
    for (int i = 0; i < kernel.size; i++) {
        kernel[i] = kernel[i] < T(1./255.) ? T(0.) : kernel[i];
    }
    kernel.normalize();
}

/// center the kernel at the center of the image
template <typename T>
static void centerKernel(img_t<T>& kernel)
{
    // compute its barycenter
    float dx = 0.f;
    float dy = 0.f;
    for (int y = 0; y < kernel.h; y++) {
        for (int x = 0; x < kernel.w; x++) {
            dx += kernel(x, y) * x;
            dy += kernel(x, y) * y;
        }
    }
    dx = std::round(dx);
    dy = std::round(dy);

    // center the kernel
    img_t<T> copy(kernel);
    for (int y = 0; y < kernel.h; y++) {
        for (int x = 0; x < kernel.w; x++) {
            int nx = (x + (int)dx + (kernel.w/2+1)) % kernel.w;
            int ny = (y + (int)dy + (kernel.h/2+1)) % kernel.h;
            kernel(x, y) = copy(nx, ny);
        }
    }
}

/// evaluate a kernel on a given blurry subimage
template <typename T>
static T evaluateKernel(const img_t<T>& kernel, const img_t<T>& blurredPatch, T deconvLambda)
{
    assert(blurredPatch.d == 1);

    // pad and deconvolve the patch
    img_t<T> paddedBlurredPatch;
    pad_and_taper(paddedBlurredPatch, blurredPatch, kernel);
    img_t<T> deconvPadded;
    deconvBregman(deconvPadded, paddedBlurredPatch, kernel, 10, deconvLambda);
    img_t<T> deconv;
    unpad(deconv, deconvPadded, kernel);

    // compute the l1 and l2 norm of the gradient of the deconvolved patch
    T normL1 = 0.;
    T normL2p2 = 0.;
    for (int y = 1; y < deconv.h; y++) {
        for (int x = 1; x < deconv.w; x++) {
            T dx = deconv(x, y) - deconv(x - 1, y);
            T dy = deconv(x, y) - deconv(x, y - 1);
            T norm = std::sqrt(dx*dx + dy*dy);
            normL1 += norm;
            normL2p2 += norm*norm;
        }
    }

    // returns the score of the kernel
    return normL1 / std::sqrt(normL2p2);
}

/// Algorithm 5
template <typename T>
void phaseRetrieval(img_t<T>& outkernel, const img_t<T>& blurredPatch,
                    const img_t<T>& powerSpectrum, int kernelSize,
                    const options& opts)
{
    img_t<T> magnitude(powerSpectrum.w, powerSpectrum.h);
    for (int i = 0; i < powerSpectrum.size; i++)
        magnitude[i] = std::sqrt(powerSpectrum[i]);
    magnitude.ifftshift(); // unshift the magnitude

    T globalCurrentScore = std::numeric_limits<T>::max();
#pragma omp parallel
    {
        img_t<T> kernel;
        img_t<T> kernel_mirror;
        T currentScore = std::numeric_limits<T>::max();
        img_t<T> bestKernel;
#pragma omp for nowait
        for (int k = 0; k < opts.Ntries; k++) {
            // retrieve one possible kernel
            singlePhaseRetrieval(kernel, magnitude, kernelSize, opts.Ninner);
            centerKernel(kernel);

            // mirror the kernel (because the phase retrieval can't distinguish between the kernel and its mirror)
            kernel_mirror.ensure_size(kernel.w, kernel.h);
            for (int y = 0; y < kernel.h; y++) {
                for (int x = 0; x < kernel.w; x++) {
                    kernel_mirror(x, y) = kernel(kernel.w-1 - x, kernel.h-1 - y);
                }
            }

            // evaluate the two kernels
            img_t<T>* kernels[2] = {&kernel, &kernel_mirror};
            T scores[2];
            for (int i = 0; i < 2; i++) {
                scores[i] = evaluateKernel(*(kernels[i]), blurredPatch, opts.intermediateDeconvolutionWeight);
            }

            // keep the best one
            if (scores[1] < scores[0]) {
                scores[0] = scores[1];
                kernels[0] = kernels[1];
            }

            // if the best of two is better than the current best, keep it
            if (scores[0] < currentScore) {
                currentScore = scores[0];
                bestKernel = *(kernels[0]);
            }
        }

        // aggregate results by keeping the best kernel
#pragma omp critical
        {
            if (currentScore < globalCurrentScore) {
                globalCurrentScore = currentScore;
                outkernel = bestKernel;
            }
        }
    }
}

template <typename T>
using linear_map_t = void(*)(T *y, T *x, int n, void *e);

template <typename T>
static T scalar_product(T *x, T *y, int n)
{
	T r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * y[i];
	return r;
}

#define FOR(i,n) for(int i = 0; i < n; i++)

template <typename T>
static void fancy_conjugate_gradient(T *x,
		linear_map_t<T> A, const T *b, int n, void *e,
		T *x0, int max_iter, T min_residual)
{
	T *r  = (T*) malloc(n * sizeof(T));
	T *p  = (T*) malloc(n * sizeof(T));
	T *Ap = (T*) malloc(n * sizeof(T));

	A(Ap, x0, n, e);

	FOR(i,n) x[i] = x0[i];
	FOR(i,n) r[i] = b[i] - Ap[i];
	FOR(i,n) p[i] = r[i];

	for (int iter = 0; iter < max_iter; iter++) {
		A(Ap, p, n, e);
		T   App    = scalar_product(Ap, p, n);
		T   rr_old = scalar_product(r, r, n);
		T   alpha  = rr_old / App;
		FOR(i,n) x[i]   = x[i] + alpha * p[i];
		FOR(i,n) r[i]   = r[i] - alpha * Ap[i];
		T   rr_new = scalar_product(r, r, n);
		/*fprintf(stderr, "iter=%d, rr_new=%g\n", iter, rr_new);*/
		if (sqrt(rr_new) < min_residual)
			break;
		T   beta   = rr_new / rr_old;
		FOR(i,n) p[i]   = r[i] + beta * p[i];
	}

	free(r);
	free(p);
	free(Ap);
}

template <typename T>
void conjugate_gradient(T *x, linear_map_t<T> A, const T *b, int n, void *e)
{
	for (int i = 0; i < n; i++)
		x[i] = 0;
	int max_iter = n;
	T min_residual = 1e-10;

	fancy_conjugate_gradient(x, A, b, n, e, x, max_iter, min_residual);
}

template <typename T>
void computeProjectionsAutocorrelation(img_t<T>& acProjections, const img_t<T>& imgBlur,
                                       const std::vector<angle_t>& angleSet,
                                       int psSize, T compensationFactor);

/// computed the autocorrelation of a signal up to a given window size
/// values at a radius windowRadius of the border of data are not used
template <typename T>
void computeAutocorrelation(std::vector<T>& out, const std::vector<T>& data,
                            int windowRadius)
{
    out.resize(windowRadius*2+1);
    std::fill(out.begin(), out.end(), 0.);

    int len = data.size();
    for (int i = 0; i <= windowRadius; i++) {
        T suma = 0.;
        T sumb = 0.;
        for (int j = 0; j < len - (windowRadius*2 + 1); j++) {
            suma += data[windowRadius + j] * data[windowRadius - i + j];
            sumb += data[windowRadius + j] * data[windowRadius + i + j];
        }

        T res = (suma + sumb) / 2.;
        out[windowRadius-i] = res;
        out[windowRadius+i] = res;
    }

    for (int i = 0; i < windowRadius*2+1; i++) {
        out[i] /= len - (windowRadius*2+1);
    }
}

/// whiten an image by convolving it with a 9 points 1D differentiation filter
/// the filter is applied to rows and columns (returns two images)
template <typename T>
static void whitenImage(img_t<T>& imgBlurX, img_t<T>& imgBlurY,
                        const img_t<T>& imgBlur)
{
    const T filter[] = { 3/840., -32/840., 168/840., -672/840.,
                         0, 672/840., -168/840., 32/840., -3/840. };
    int filterSize = sizeof(filter) / sizeof(*filter);
    int w = imgBlur.w;
    int h = imgBlur.h;

    imgBlurX.ensure_size(w, h);
    imgBlurX.set_value(0);
    imgBlurY.ensure_size(w, h);
    imgBlurY.set_value(0);

    // apply the filter vertically and horizontally
    for (int y = filterSize/2; y < h - filterSize/2; y++) {
        for (int x = filterSize/2; x < w - filterSize/2; x++) {
            for (int i = 0; i < filterSize; i++) {
                imgBlurX(x, y) += filter[filterSize-1 - i] * imgBlur(x + i - filterSize/2, y);
                imgBlurY(x, y) += filter[filterSize-1 - i] * imgBlur(x, y + i - filterSize/2);
            }
        }
    }
}

template <typename T>
void conv(T *y, T *x, int n, void *e)
{
    std::vector<T>& filter = *(std::vector<T>*)e;
    assert((int)filter.size() == n);
    for (int i = 0; i < n; i++) {
        y[i] = 0;
        for (int j = -n/2; j < n/2; j++) {
            int jj = j + n/2;
            // repeating boundaries
            int idx = std::max(0, std::min(n-1, i - j));
            y[i] += filter[jj] * x[idx];
        }
    }
}

template <typename T>
static void deconvolveAutocorrelation(std::vector<T>& acCompensatedRow,
                                      const std::vector<T>& acRow,
                                      const std::vector<T>& compensationFilter)
{
    auto deconvRow = acRow;
    // solve `acRow = convolve(compensationFilter, deconvRow)` to find deconvRow
    conjugate_gradient<T>(&deconvRow[0], conv, &acRow[0], acRow.size(), (void*)&compensationFilter);

    // detect negatives values in the center of the row
    bool hasNegatives = false;
    for (unsigned x = deconvRow.size() / 2 - 2; x <= deconvRow.size() / 2 + 2; x++)
        hasNegatives |= deconvRow[x] < 0.;

    // if there are some negative values, revert the changes
    if (hasNegatives) {
        acCompensatedRow = acRow;
    } else {
        acCompensatedRow = deconvRow;
    }
}

/// computes the autocorrelation of the projections of the whitened image
template <typename T>
void computeProjectionsAutocorrelation(img_t<T>& acProjections, const img_t<T>& imgBlur,
                                       const std::vector<angle_t>& angleSet,
                                       int psSize, T compensationFactor)
{
    img_t<T> projections;

    // whiten the image (horizontal and vertical filtering with the filter 'd')
    img_t<T> imgBlurX;
    img_t<T> imgBlurY;
    whitenImage(imgBlurX, imgBlurY, imgBlur);

    // compute the projections of the derivative
    projectImage(projections, imgBlurX, imgBlurY, angleSet);

    acProjections.ensure_size(psSize*2 + 1, projections.h);

    // build the compensation filter
    // k(x) = 1 / x^compensationFactor
    std::vector<T> compensationFilter(acProjections.w);
    if (compensationFactor > 0.) {
        int center = compensationFilter.size() / 2;
        T sum = 0.;
        for (int i = 0; i < (int) compensationFilter.size(); i++) {
            compensationFilter[i] = 1. / std::pow(std::abs(i - center) + 1, compensationFactor);
            sum += compensationFilter[i];
        }
        for (unsigned int i = 0; i < compensationFilter.size(); i++) {
            compensationFilter[i] /= sum;
        }
    }

#pragma omp parallel
    {
        std::vector<T> proj1d(projections.h);
#pragma omp for
        for (int j = 0; j < projections.h; j++) {
            // extract meaningful values of the projections
            proj1d.resize(0);
            for (int i = 0; i < projections.w; i++) {
                if (!std::isnan(projections(i, j))) {
                    proj1d.push_back(projections(i, j));
                }
            }

            // compute the mean
            T mean = 0.;
            for (T v : proj1d)
                mean += v;
            mean /= proj1d.size();

            // center
            for (T& v : proj1d)
                v -= mean;

            // compute the norm
            T norm = 0.;
            for (T v : proj1d)
                norm += v * v;
            norm = std::sqrt(norm);

            // normalize
            for (T& v : proj1d)
                v /= norm;

            // compute autocorrelation of the projection
            std::vector<T> autocorrelation;
            computeAutocorrelation(autocorrelation, proj1d, psSize);

            // apply the compensation filter
            if (compensationFactor > 0.) {
                // deconvolve the autocorrelation with the compensation filter
                deconvolveAutocorrelation(autocorrelation, autocorrelation, compensationFilter);
            }

            // extract the central part of the autocorrelation
            auto start = &autocorrelation[0] + autocorrelation.size() / 2 - acProjections.w / 2;
            std::copy(start, start + acProjections.w, &acProjections(0, j));
        }
    }
}

template <typename T>
void estimateKernel(img_t<T>& kernel, const img_t<T>& img,
                    int kernelSize, const options& opts);

/// search a patch with high variance in the greyscale blurred image
template <typename T>
static void searchBlurredPatch(img_t<T>& window, const img_t<T>& blurredImage,
                               int windowSize, int searchSamples)
{
    assert(blurredImage.d == 1);
    T best = 0;
    int best_x = 0;
    int best_y = 0;
    for (int i = 0; i < searchSamples; i++) {
        // draw a random position
        int x = (blurredImage.w - windowSize) * (T)rand() / RAND_MAX;
        int y = (blurredImage.h - windowSize) * (T)rand() / RAND_MAX;

        // compute the mean
        T mean = 0.;
        for (int yy = 0; yy < windowSize; yy++) {
            for (int xx = 0; xx < windowSize; xx++) {
                mean += blurredImage(x + xx, y + yy);
            }
        }
        mean /= windowSize * windowSize;

        // compute the variance
        T var = 0.;
        for (int yy = 0; yy < windowSize; yy++) {
            for (int xx = 0; xx < windowSize; xx++) {
                var += std::pow(blurredImage(x + xx, y + yy) - mean, 2);
            }
        }
        var /= windowSize * windowSize;

        // if the variance is greater, we keep the position
        if (var > best) {
            best = var;
            best_x = x;
            best_y = y;
        }
    }

    // copy the chosen window to the output
    window.ensure_size(windowSize, windowSize);
    for (int yy = 0; yy < windowSize; yy++) {
        for (int xx = 0; xx < windowSize; xx++) {
            window(xx, yy) = blurredImage(best_x + xx, best_y + yy);
        }
    }
}

/// apply a circular median filter on each column
template <typename T>
static void circularMedianFilter(img_t<T>& img, int size)
{
    img_t<T> copy(img);
    std::vector<T> pixels(size);
    for (int y = 0; y < img.h; y++) {
        for (int x = 0; x < img.w; x++) {
            // extract 'size' pixels of the column
            // centered around (x,y)
            for (int s = 0; s < size; s++) {
                int yy = (y + s - size/2 + img.h) % img.h;
                pixels[s] = copy(x, yy);
            }

            // sort the value
            std::sort(pixels.begin(), pixels.end());

            // save the median
            if (size % 2) {
                img(x, y) = pixels[size/2];
            } else {
                img(x, y) = (pixels[size/2-1] + pixels[size/2]) / 2.;
            }
        }
    }
}

/// adjust the autocorrelations using the given support (substract mu_theta and clip to 0, normalize, and median filter)
/// solve eq. (39), see Algorithm 4
template <typename T>
static void adjustAutocorrelations(img_t<T>& acProjectionsCorrected,
                                   const img_t<T>& acProjections,
                                   const std::vector<int>& support,
                                   bool medianFilter)
{
    int h = acProjections.h;
    int w = acProjections.w;

    acProjectionsCorrected.ensure_size(w, h);

    // for each orientation, cut and normalize the autocorrelation
    // according to the estimated support
    for (int y = 0; y < h; y++) {
        int cut = w / 2 - support[y];
        T cutVal = acProjections(cut, y);
        for (int x = 0; x < w; x++) {
            if (x <= cut || x >= w-1 - cut) {
                // set to 0 outside the support
                acProjectionsCorrected(x, y) = 0.;
            } else {
                // subtract mu_theta and clip to 0
                acProjectionsCorrected(x, y) = std::max(T(0.), acProjections(x, y) - cutVal);
            }
        }

        // normalize (compensate c_theta)
        T sum = 0.;
        for (int x = 0; x < w; x++) {
            sum += acProjectionsCorrected(x, y);
        }
        for (int x = 0; x < w; x++) {
            acProjectionsCorrected(x, y) /= sum;
        }
    }

    // apply a median filter across angles to smooth the coefficients
    if (medianFilter) {
        circularMedianFilter(acProjectionsCorrected, std::round(2.*std::sqrt(acProjections.h)));
    }
}

/// computed the autocorrelation of a signal up to a given window size
/// uses valid boundary condition
/// (this function differs from the computeAutocorrelation in the file computeProjectionsAutocorrelation.cpp
///  by the fact that this one is adapted to small support of data)
template <typename T>
void computeAutocorrelationSmallSupport(std::vector<T>& out, const std::vector<T>& data,
                                        int windowRadius)
{
    out.resize(windowRadius*2+1);
    std::fill(out.begin(), out.end(), 0.);

    int len = data.size();
    for (int i = -windowRadius; i <= windowRadius; i++) {
        T sum = 0.;

        for (int j = 0; j < len; j++) {
            if (j - i >= 0 && j - i < len) {  // valid boundary condition
                sum += data[j] * data[j - i];
            }
        }

        out[windowRadius + i] = sum;
    }
}

/// Reevaluate the support of the projection of the kernel
template <typename T>
static void reestimateKernelSupport(std::vector<int>& support, const img_t<T>& kernel,
                                    const std::vector<angle_t>& angleSet,
                                    int acRadius, T threshold=5e-2)
{
    std::vector<T> angles(angleSet.size());
    for (unsigned j = 0; j < angleSet.size(); j++) {
        angles[j] = angleSet[j].angle;
    }

    // compute the shear projections of the kernel
    img_t<T> shearProjections;
    projectImage(shearProjections, kernel, angleSet);
    shearProjections.transpose();

    support.resize(angleSet.size());

    img_t<T> ac(acRadius*2+1, angles.size());

    std::vector<T> autocorrelation;
    std::vector<T> proj(shearProjections.h);
    // for each orientation, compute the autocorrelation of the estimated kernel,
    // and estimate its support
    for (unsigned j = 0; j < angles.size(); j++) {
        // extract the projection
        for (int i = 0; i < shearProjections.h; i++) {
            proj[i] = shearProjections(j, i);
            if (std::isnan(proj[i]))
                proj[i] = 0.;
        }

        // compute the autocorrelation of the projection
        computeAutocorrelationSmallSupport(autocorrelation, proj, acRadius);

        // find the max of the autocorrelation (for the adaptative threshold)
        T max = *std::max_element(autocorrelation.begin(), autocorrelation.end());

        // search for the first value exceeding the threshold on the right side of the autocorrelation
        unsigned cursor = autocorrelation.size() - 1;
        while (cursor > autocorrelation.size()/2 && autocorrelation[cursor] < max*threshold)
            cursor--;
        // go back to the value that was below the threshold
        cursor++;
        // and substract half of the autocorrelation size to be within [0, acRadius]
        cursor -= autocorrelation.size() / 2;

        // set the support
        support[j] = cursor;
    }
}


/// estimate the kernel support from the autocorrelation of whitened projections
/// Algorithm 3
template <typename T>
static void initialSupportEstimation(std::vector<int>& support,
                                     const img_t<T> acProjections,
                                     T maxSlope=20./700.)
{
    int h = acProjections.h;
    int w = acProjections.w;

    support.resize(h);

    // local minima of the autocorrelations
    for (int j = 0; j < h; j++) {
        const T* start = &acProjections(w / 2, j);
        support[j] = std::min_element(start, start + w / 2) - start;
    }

    // apply the Lipschitz continuity constraint on the support
    std::vector<T> currentMinimums(h);
    std::fill(currentMinimums.begin(), currentMinimums.end(), w / 2);

    for (int j = 0; j < h; j++) {
        // if the support is lower than the current maximum, we apply the constraint
        if (support[j] < currentMinimums[j]) {
            currentMinimums[j] = support[j];

            // propagate the constraint to other values (to both sides at the same time)
            for (int yy = 1; yy <= h / 2; yy++) {
                // minimum value that is allowed at a distance 'yy' from the current position
                T minimumValue = support[j] + yy * maxSlope;

                int prevRow = (j - yy + h) % h;
                currentMinimums[prevRow] = std::min(currentMinimums[prevRow], minimumValue);

                int nextRow = (j + yy) % h;
                currentMinimums[nextRow] = std::min(currentMinimums[nextRow], minimumValue);
            }
        }
    }

    for (int j = 0; j < h; j++) {
        support[j] = std::ceil(currentMinimums[j]);
    }
}

template <typename T>
void goldstein_fattal(img_t<T>& kernel, const img_t<T>& img,
                    int kernelSize, const options& opts);
