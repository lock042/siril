#include "goldfat.hpp"

/// estimate the kernel from a blurred image and a kernel size
/// Algorithm 1 of the paper
template <typename T>
void goldstein_fattal(img_t<T>& kernel, const img_t<T>& img,
                    int kernelSize, const options& opts)
{
    kernel.ensure_size(kernelSize, kernelSize);

    // convert the image to greyscale
    img_t<T> grey(img.w, img.h);
    grey.greyfromcolor(img);

    // search a blurred patch which will be used for kernel evaluation
    img_t<T> blurredPatch;
    searchBlurredPatch(blurredPatch, grey, 150, 100);

    // compute the angle set
    std::vector<angle_t> angleSet;
    computeProjectionAngleSet(angleSet, kernelSize*2);

    // compute the autocorrelation of the projection of the whitened image
    img_t<T> acProjections;
    computeProjectionsAutocorrelation(acProjections, grey, angleSet, kernelSize*2, opts.compensationFactor);
    int acRadius = acProjections.w / 2;

    // initial support estimation
    std::vector<int> support;
    initialSupportEstimation(support, acProjections);

    // iterative estimation
    img_t<T> powerSpectrum;
    img_t<T> acProjectionsCorrected;
    for (int i = 0; i < opts.Nouter; i++) {
        // adjust the autocorrelation using the estimated support
        adjustAutocorrelations(acProjectionsCorrected, acProjections, support, opts.medianFilter);

        // compute the power spectrum from the autocorrelation
        reconstructPowerspectrum(powerSpectrum, acProjectionsCorrected, angleSet, acRadius);

        // retrieve a kernel in spatial domain using the power spectrum
        phaseRetrieval(kernel, blurredPatch, powerSpectrum, kernelSize, opts);

        // reestimate the kernel support
        reestimateKernelSupport(support, kernel, angleSet, acRadius);
    }
}
