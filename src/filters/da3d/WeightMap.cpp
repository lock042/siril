/*
 * WeightMap.cpp
 *
 *  Created on: 12/feb/2015
 *      Author: nicola pierazzo <nicola.pierazzo@cmla.ens-cachan.fr>
 */

#include <cassert>
#include <limits>
#include <cstdlib>
#include <utility>
#include <tuple>
#include <algorithm>
#include "WeightMap.hpp"
#include "Image.hpp"
#include "Utils.hpp"

using std::max;
using std::min;
using std::pair;
using std::make_pair;
using std::tie;

namespace da3d {

WeightMap::WeightMap(int rows, int columns) {
  Init(rows, columns);
}

void WeightMap::Init(int rows, int columns) {
  assert (rows > 0 && columns > 0);
  int height_rows = utils::NumberOfBits(rows - 1) + 1;
  int height_columns = utils::NumberOfBits(columns - 1) + 1;
  num_levels_ = max(height_rows, height_columns);
  width_ = columns;
  height_ = rows;
  rows_.resize(num_levels_);
  columns_.resize(num_levels_);
  data_.resize(num_levels_);

  // initialization of layers
  int rows_rounded = utils::NextPowerOf2(rows);
  int cols_rounded = utils::NextPowerOf2(columns);
  for (int l = 0; l < num_levels_; ++l) {
    rows_[l] = rows_rounded;
    columns_[l] = cols_rounded;
    data_[l].resize(rows_rounded * cols_rounded);
    // zeros in the good area, MAXFLT elsewhere
    for (int row = 0; row < rows; ++row) {
      for (int col = 0; col < columns; ++col)
        val(col, row, l) = 0.f;
      for (int col = columns; col < cols_rounded; ++col) {
        val(col, row, l) = std::numeric_limits<float>::infinity();
      }
    }
    for (int row = rows; row < rows_rounded; ++row) {
      for (int col = 0; col < cols_rounded; ++col) {
        val(col, row, l) = std::numeric_limits<float>::infinity();
      }
    }
    rows = (rows + 1) >> 1;  // it stays at least 1
    columns = (columns + 1) >> 1;
    rows_rounded = ((rows_rounded + 3) >> 2) << 1;  // it stays at least 2
    cols_rounded = ((cols_rounded + 3) >> 2) << 1;
  }
  assert(rows == 1);
  assert(columns == 1);
}

float WeightMap::Minimum() const {
  return val(0, 0, num_levels_ - 1);
}

pair<int, int> WeightMap::FindMinimum() const {
  int row = 0;
  int col = 0;
  for (int l = num_levels_ - 2; l >= 0; --l) {
    row <<= 1;
    col <<= 1;
    // This finds the position of the minimum in the 2x2 square
    tie(row, col) = min({make_pair(row, col),     make_pair(row + 1, col),
                         make_pair(row, col + 1), make_pair(row + 1, col + 1)},
                        [this, &l](const pair<int, int> &a,
                                   const pair<int, int> &b) {
                          return val(a.second, a.first, l)
                              < val(b.second, b.first, l);
                        });
  }
  return {row, col};
}

void WeightMap::IncreaseWeights(const Image &weights, int row0, int col0) {
  assert(weights.channels() == 1);
  int firstrow = max(0, row0);
  int lastrow = min(height(), row0 + weights.rows()) - 1;
  int firstcol = max(0, col0);
  int lastcol = min(width(), col0 + weights.columns()) - 1;

  // Updates the level zero
  for (int row = firstrow; row <= lastrow; ++row) {
    for (int col = firstcol; col <= lastcol; ++col) {
        val(col, row) += weights.val(col - col0, row - row0);
    }
  }
  // Updates the other levels
  for (int l = 1; l < num_levels_; ++l) {
    for (int row = firstrow >> l; row <= lastrow >> l; ++row) {
      for (int col = firstcol >> l; col <= lastcol >> l; ++col) {
        int dc = col << 1, dr = row << 1;
        val(col, row, l) = min({val(dc,     dr,     l - 1),
                                val(dc + 1, dr,     l - 1),
                                val(dc,     dr + 1, l - 1),
                                val(dc + 1, dr + 1, l - 1)});
      }
    }
  }
}

}  // namespace da3d
