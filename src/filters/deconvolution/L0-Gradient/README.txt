
This code is based on the code from "Blind Image Deblurring using the l0 Gradient Prior".

Jérémy Anger, Gabriele Facciolo, Mauricio Delbracio

This program is part of the IPOL publication:
    https://www.ipol.im/pub/art/2019/243/

Version 20190110

Compilation:
    run "make" to produce executables named "estimate-kernel".
    requires a C++11 compatible compiler and the following libraries: libpng, libtiff, libjpeg, libfftw3

Usage:
    ./estimate-kernel KERNEL_SIZE BLURRY_IMAGE KERNEL_OUTPUT [options]

        KERNEL_SIZE: should be an odd integer large enough to contains the actual estimated kernel
        BLURRY_IMAGE: should be a tiff, png or jpeg file.
        KERNEL_OUTPUT: output file for the estimated kernel, should be a .tif in order to keep floating point values
        [options]: allows to change the parameters of the method, here is the list:
            --lambda=[lambda]                 L0 regularization weight
            --lambda-ratio=[lambda-ratio]     decay of lambda
            --lambda-min=[lambda-min]         L0 regularization weight minimum value
            --gamma=[gamma]                   kernel regularization weight
            --iterations=[iterations]         number of iterations per scale
            --no-multiscale                   disable the multiscale scheme
            --scale-factor=[scale-factor]     downsampling factor
            --kernel-threshold-max=[kernel-threshold-max] threshold the kernel at max(kernel)*kernel-threshold-max
            --remove-isolated=[remove-isolated] remove isolated connected component of the kernel
            --output-sharp=[output-sharp]     output the sharp image to file
            --debug=[debug]                   output all kernels, sharp and blurry images
            --verbose                         output more information

    For more info, use "--help"

Example:
    ./estimate-kernel 15 hollywood.jpg kernel.tif

Credits:
    iio and imscript: from https://github.com/mnhrdt/imscript
    args.hxx: from https://github.com/Taywee/args

