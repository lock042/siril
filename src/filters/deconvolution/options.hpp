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
    float compensationFactor;
    int medianFilter;

    float finalDeconvolutionWeight;
    float intermediateDeconvolutionWeight;
    int seed;
};

