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

}  // namespace utils

#endif  // DA3D_UTILS_HPP_
