/*
 * Utils.hpp
 *
 *  Created on: 23/mar/2015
 *      Author: nicola pierazzo <nicola.pierazzo@cmla.ens-cachan.fr>
 */

#ifndef DA3D_UTILS_HPP_
#define DA3D_UTILS_HPP_

#include <cstring>
#include <string>
#include <vector>
#include <utility>
#include "Image.hpp"
#include "WeightMap.hpp"

namespace utils {

// number of bits used to represent n
inline int NumberOfBits(int n) {
  int ans = 0;
  while (n) {
    ans += 1;
    n >>= 1;
  }
  return ans;
}

// http://graphics.stanford.edu/~seander/bithacks.html
inline int NextPowerOf2(int n) {
  --n;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  ++n;
  return n;
}

inline float fastexp(float x) {
  int result = static_cast<int>(12102203 * x) + 1065353216;
  result *= result > 0;
  std::memcpy(&x, &result, sizeof(result));
  return x;
}

inline int SymmetricCoordinate(int pos, int size) {
  if (pos < 0) pos = -pos - 1;
  if (pos >= 2 * size) pos %= 2 * size;
  if (pos >= size) pos = 2 * size - 1 - pos;
  return pos;
}

bool isMonochrome (const da3d::Image &u);

da3d::Image makeMonochrome (const da3d::Image &u);

std::pair<int, int> ComputeTiling(int rows, int columns, int tiles);
std::vector<da3d::Image> SplitTiles(const da3d::Image &src, int pad_before,
                                    int pad_after, std::pair<int, int> tiling);
da3d::Image MergeTiles(const std::vector<std::pair<da3d::Image, da3d::Image>> &src,
                       std::pair<int, int> shape, int pad_before, int pad_after,
                       std::pair<int, int> tiling);
}  // namespace utils

#endif  // DA3D_UTILS_HPP_
