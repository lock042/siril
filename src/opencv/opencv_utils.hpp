#pragma once

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <opencv2/core/core.hpp>
#include "io/image_format_fits.h"

/* functions from opencv.cpp that cannot be declared in a C header, but useful for a siril C++ library */

int image_to_Mat(fits *image, cv::Mat *in, cv::Mat *out, void **bgr, int target_rx, int target_ry);
int Mat_to_image(fits *image, cv::Mat *in, cv::Mat *out, void *bgr, int target_rx, int target_ry);
void init_guide(fits *image, unsigned int target_rx, unsigned int target_ry, cv::Mat *guide);

