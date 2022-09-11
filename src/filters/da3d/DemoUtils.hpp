/*
 * DemoUtils.cpp
 *
 *  Created on: 10/apr/2017
 *      Author: gabriele facciolo <facciolo@cmla.ens-cachan.fr>
 */

#ifndef UTILS_DEMOUTILS_HPP
#define UTILS_DEMOUTILS_HPP

#include "Image.hpp"
#include <string>

namespace utils {

const char *pick_option(int *c, char **v, const char *o, const char *d);
da3d::Image read_image(const std::string& filename);
void save_image(const da3d::Image& image, const std::string& filename);

}  // namespace utils

#endif //UTILS_DEMOUTILS_HPP
