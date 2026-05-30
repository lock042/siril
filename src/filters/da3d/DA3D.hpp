/*
 * da3d.hpp
 *
 *  Created on: 24/mar/2015
 *      Author: nicola pierazzo <nicola.pierazzo@cmla.ens-cachan.fr>
 */

#ifndef DA3D_DA3D_HPP_
#define DA3D_DA3D_HPP_

#include "algos/img_t/image.hpp"

namespace da3d {

img_t<float> DA3D(int &retval, const img_t<float> &noisy, const img_t<float> &guide, float sigma,
           int nthreads = 0, int r = 31, float sigma_s = 14.f,
           float gamma_r = .7f, float threshold = 2.f);

}  // namespace da3d

#endif  // DA3D_DA3D_HPP_
