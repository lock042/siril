#pragma once

#include <string>

struct options {
    std::string input;
    int kernelSize;
    std::string out_kernel;
    std::string out_deconv;

    int Ninner;
    int Ntries;
    int Nouter;
    flt compensationFactor;
    int medianFilter;

    flt finalDeconvolutionWeight;
    flt intermediateDeconvolutionWeight;
    int seed;
};

