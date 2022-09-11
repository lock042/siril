/*
 * main.cpp
 *
 *  Created on: 24/mar/2015
 *      Author: nicola pierazzo <nicola.pierazzo@cmla.ens-cachan.fr>
 */

#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <utility>
#include "Image.hpp"
#include "Utils.hpp"
#include "DemoUtils.hpp"
#include "DA3D.hpp"
extern "C" {
#include "core/proto.h"
}

using std::cerr;
using std::endl;

using utils::pick_option;
using utils::read_image;
using utils::save_image;
using utils::isMonochrome;
using utils::makeMonochrome;

using da3d::Image;
using da3d::DA3D;

float* call_da3d(float *in, float *gd, int width, int height, int nchans, float sigma) {

#ifndef _OPENMP
  siril_log_message(_("OpenMP not available. The DA3D algorithm will run in a single thread.\n"));
#endif

  Image input(in, height, width, nchans);
  Image guide(gd, height, width, nchans);
  // DA3D doesn't work if a color image has monochromatic noise
  if (input.channels()>1 && isMonochrome(input)) {
     siril_log_color_message(_("Warning: input color image has monochromatic noise! Converting to monochrome."), "red");
     input = makeMonochrome(input);
     guide = makeMonochrome(guide);
  }
  Image output = DA3D(input, guide, sigma);
  float *out = output.data();
  return out;
}
