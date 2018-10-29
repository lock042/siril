#ifndef SIRIL_OPENCV_H_
#define SIRIL_OPENCV_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef __cplusplus

// This is a C++ function reused in another C++ file, ecc.cpp.
#include <opencv2/core/core.hpp>
void convert_MatH_to_H(cv::Mat from, Homography *to);

extern "C" {
#endif

#include <stdint.h>
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "gui/progress_and_log.h"


WORD *cvResizeGaussian_buf(WORD *buf, int srcx, int srcy, int destx, int desty);
int cvResizeGaussian(fits *, int, int, int);
int cvResizeGaussian_data8(uint8_t *dataIn, int rx, int ry, uint8_t *dataOut,
		int toX, int toY, int chan, int interpolation);
int cvRotateImage(fits *, double, int, int);
int cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *H);
int cvApplyScaleToH(Homography *H1, double scale);
void cvTransformBuf(WORD *image, int size, Homography *Hom);
int cvTransformImage(fits *, point, Homography, int);
int cvUnsharpFilter(fits*, double, double);
int cvComputeFinestScale(fits *image);
int cvLucyRichardson(fits *image, double sigma, int iterations);
int cvLaplacian(fits *image);
#ifdef __cplusplus
}
#endif

#endif
