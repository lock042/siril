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

int call_da3d(float *in, float *gd, float *out, unsigned width, unsigned height, unsigned nchans, float sigma) {

#ifndef _OPENMP
  siril_log_message(_("OpenMP not available. The DA3D algorithm will run in a single thread.\n"));
#endif

  Image input(in, (int) height, (int) width, (int) nchans);
  Image guide(gd, (int) height, (int) width, (int) nchans);
//  float sigma = atof(argv[3]);
  // DA3D doesn't work if a color image has monochromatic noise
  if (input.channels()>1 && isMonochrome(input)) {
     siril_log_color_message(_("Warning: input color image has monochromatic noise! Converting to monochrome."), "red");
     input = makeMonochrome(input);
     guide = makeMonochrome(guide);
  }
  Image output = DA3D(input, guide, sigma);
  out = output.data();
  return EXIT_SUCCESS;
}
