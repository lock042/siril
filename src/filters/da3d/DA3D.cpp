/*
 * DA3D.cpp
 *
 *  Created on: 24/mar/2015
 *      Author: nicola pierazzo <nicola.pierazzo@cmla.ens-cachan.fr>
 */

#include <cmath>
#include <tuple>
#include <utility>
#include <algorithm>
#include <numeric>
#include "Image.hpp"
#include "DA3D.hpp"
#include "WeightMap.hpp"
#include "Utils.hpp"
#include "DftPatch.hpp"

#include "core/processing.h"
#include "gui/progress_and_log.h"


#ifdef _OPENMP
#include <omp.h>
#endif  // OPENMP

using std::max;
using std::min;
using std::vector;
using std::pair;
using std::tie;
using std::move;
using std::sqrt;
using std::log2;
using std::floor;
using std::modf;
using std::abs;
using std::accumulate;
using std::norm;
using utils::fastexp;
using utils::NextPowerOf2;
using utils::ComputeTiling;
using utils::SplitTiles;
using utils::MergeTiles;

namespace da3d {

namespace {

Image ColorTransform(Image&& src) {
  Image img = move(src);
  if (img.channels() == 3) {
    for (int row = 0; row < img.rows(); ++row) {
      for (int col = 0; col < img.columns(); ++col) {
        float r, g, b;
        r = img.val(col, row, 0);
        g = img.val(col, row, 1);
        b = img.val(col, row, 2);
        img.val(col, row, 0) = (r + g + b) / sqrt(3.f);
        img.val(col, row, 1) = (r - b) / sqrt(2.f);
        img.val(col, row, 2) = (r - 2 * g + b) / sqrt(6.f);
      }
    }
  }
  return img;
}

Image ColorTransformInverse(Image&& src) {
  Image img = move(src);
  if (img.channels() == 3) {
    for (int row = 0; row < img.rows(); ++row) {
      for (int col = 0; col < img.columns(); ++col) {
        float y, u, v;
        y = img.val(col, row, 0);
        u = img.val(col, row, 1);
        v = img.val(col, row, 2);
        img.val(col, row, 0) = (sqrt(2.f) * y + sqrt(3.f) * u + v) / sqrt(6.f);
        img.val(col, row, 1) = (y - sqrt(2.f) * v) / sqrt(3.f);
        img.val(col, row, 2) = (sqrt(2.f) * y - sqrt(3.f) * u + v) / sqrt(6.f);
      }
    }
  }
  return img;
}

void ExtractPatch(const Image &src, int pr, int pc, Image *dst) {
  // src is padded, so (pr, pc) becomes the upper left pixel
  int i = 0, j = (pr * src.columns() + pc) * src.channels();
  for (int row = 0; row < dst->rows(); ++row) {
    for (int el = 0; el < dst->columns() * dst->channels(); ++el) {
      dst->val(i) = src.val(j);
      ++i;
      ++j;
    }
    j += (src.columns() - dst->columns()) * src.channels();
  }
}

void BilateralWeight(const Image &g, Image *k, int r, float gamma_r_sigma2,
                     float sigma_s2) {
  for (int row = 0; row < g.rows(); ++row) {
    for (int col = 0; col < g.columns(); ++col) {
      float x = 0.f;
      for (int chan = 0; chan < g.channels(); ++chan) {
        float y = g.val(col, row, chan) - g.val(r, r, chan);
        x += y * y;
      }
      x /= gamma_r_sigma2;
      x += ((row - r) * (row - r) + (col - r) * (col - r)) / (2 * sigma_s2);
      k->val(col, row) = utils::fastexp(-x);
    }
  }
}

/* void ComputeRegressionPlaneIRLS(const Image &y,
                                const Image &g,
                                const Image &k,
                                int r,
                                vector<pair<float, float>> *reg_plane,
                                int iterations = 2) {
  constexpr float delta = 0.0001f;
  for (pair<float, float> &v : *reg_plane) v = {0.f, 0.f};
  for (int chan = 0; chan < y.channels(); ++chan) {
    float x1 = 0.f, x2 = 0.f;
    for (int t = 0; t < iterations; ++t) {
      float a = 0.f, b = 0.f, c = 0.f, d = 0.f, e = 0.f;
      float central = g.val(r, r, chan);
      for (int row = 0; row < y.rows(); ++row) {
        for (int col = 0; col < y.columns(); ++col) {
          int cc = col - r, rr = row - r;
          float w = k.val(col, row)
              / max(delta, (y.val(col, row, chan) - central - x1 * rr - x2 * cc));
          a += rr * rr * w;
          b += rr * cc * w;
          c += cc * cc * w;
          d += rr * (y.val(col, row, chan) - central) * w;
          e += cc * (y.val(col, row, chan) - central) * w;
        }
      }
      float det = a * c - b * b;
      if (abs(det) < delta) break;

      // Solves the system
      // |a   b| |x1|   |d|
      // |     | |  | = | |
      // |b   c| |x2|   |e|
      x1 = (c * d - b * e) / det;
      x2 = (a * e - b * d) / det;
    }
    (*reg_plane)[chan] = {x1, x2};
  }
} */

void ComputeRegressionPlane(const Image &y, const Image &g, const Image &k,
                            int r, vector<pair<float, float>> *reg_plane) {
  constexpr float epsilon = 0.0001f;
  float a = 0.f, b = 0.f, c = 0.f;
  for (int row = 0; row < y.rows(); ++row) {
    for (int col = 0; col < y.columns(); ++col) {
      a += (row - r) * (row - r) * k.val(col, row);
      b += (row - r) * (col - r) * k.val(col, row);
      c += (col - r) * (col - r) * k.val(col, row);
    }
  }
  float det = a * c - b * b;
  if (abs(det) < epsilon) {
    for (int chan = 0; chan < y.channels(); ++chan) {
      (*reg_plane)[chan] = {0.f, 0.f};
    }
  } else {
    for (int chan = 0; chan < y.channels(); ++chan) {
      float d = 0.f, e = 0.f;
      float central = g.val(r, r, chan);
      for (int row = 0; row < y.rows(); ++row) {
        for (int col = 0; col < y.columns(); ++col) {
          d += (row - r) * (y.val(col, row, chan) - central) * k.val(col, row);
          e += (col - r) * (y.val(col, row, chan) - central) * k.val(col, row);
        }
      }
      // Solves the system
      // |a   b| |x1|   |d|
      // |     | |  | = | |
      // |b   c| |x2|   |e|
      (*reg_plane)[chan] = {(c * d - b * e) / det, (a * e - b * d) / det};
    }
  }
}

void SubtractPlane(int r, vector<pair<float, float>> reg_plane, Image *y) {
  for (int row = 0; row < y->rows(); ++row) {
    for (int col = 0; col < y->columns(); ++col) {
      for (int chan = 0; chan < y->channels(); ++chan) {
        y->val(col, row, chan) -= reg_plane[chan].first * (row - r) +
                                  reg_plane[chan].second * (col - r);
      }
    }
  }
}

void ModifyPatch(const Image &patch, const Image &k, DftPatch *modified,
                 float *average = nullptr) {
  // compute the total weight of the mask
  float weight = accumulate(k.begin(), k.end(), 0.f);

  for (int chan = 0; chan < patch.channels(); ++chan) {
    float avg = 0.f;
    for (int row = 0; row < patch.rows(); ++row) {
      for (int col = 0; col < patch.columns(); ++col) {
        avg += k.val(col, row) * patch.val(col, row, chan);
      }
    }
    avg /= weight;
    for (int row = 0; row < patch.rows(); ++row) {
      for (int col = 0; col < patch.columns(); ++col) {
        modified->space(col, row, chan) =
            k.val(col, row) * patch.val(col, row, chan) +
            (1.f - k.val(col, row)) * avg;
      }
    }
    if (average) average[chan] = avg;
  }
}

pair<Image, Image> DA3D_block(int &retval, const Image &noisy, const Image &guide,
                              float sigma, int r, float sigma_s,
                              float gamma_r, float threshold, unsigned nThreads, unsigned thread) {
  // useful values
  const int s = utils::NextPowerOf2(2 * r + 1);
  const float sigma2 = sigma * sigma;
  const float gamma_r_sigma2 = gamma_r * sigma2;
  const float sigma_s2 = sigma_s * sigma_s;

  // Work estimate for progress bar. Empirical based on a sample of images, may need tuning
  const unsigned predicted_loops = noisy.pixels() * 5 / 1000;
  unsigned loop = 0;
  unsigned lastupdate = 0;

  // regression parameters
  const float gamma_rr_sigma2 = gamma_r_sigma2 * 10.f;
  const float sigma_sr2 = sigma_s2 * 2.f;

  // declaration of internal variables
  Image y(s, s, guide.channels());
  Image g(s, s, guide.channels());
  Image k_reg(s, s);
  Image k(s, s);
  DftPatch y_m(s, s, guide.channels());
  DftPatch g_m(s, s, guide.channels());
  int pr, pc;  // coordinates of the central pixel
  vector<pair<float, float>> reg_plane(guide.channels());  // parameters of the regression plane
  float yt[guide.channels()];  // weighted average of the patch
  WeightMap agg_weights(guide.rows() - s + 1, guide.columns() - s + 1);  // line 1

  Image output(guide.rows(), guide.columns(), guide.channels());
  Image weights(guide.rows(), guide.columns());

  // Instrumentation
  // main loop
  while (agg_weights.Minimum() < threshold) {  // line 4
    if (!get_thread_run()) {
      retval++;
      break;
    }
    if (thread == 0) {
      loop++;
      if (loop - lastupdate > (predicted_loops >> 3)) {
        lastupdate = loop;
        if ((double) loop / (nThreads * predicted_loops) < 1.0)
          set_progress_bar_data("DA3D denoising...", (double) loop / (nThreads * predicted_loops));
        else
          set_progress_bar_data("DA3D denoising...", PROGRESS_PULSATE);
      }

    }


    tie(pr, pc) = agg_weights.FindMinimum();  // line 5
    ExtractPatch(noisy, pr, pc, &y);  // line 6
    ExtractPatch(guide, pr, pc, &g);  // line 7
    BilateralWeight(g, &k_reg, r, gamma_rr_sigma2, sigma_sr2);  // line 8
    ComputeRegressionPlane(y, g, k_reg, r, &reg_plane);  // line 9
    SubtractPlane(r, reg_plane, &y);  // line 10
    SubtractPlane(r, reg_plane, &g);  // line 11
    BilateralWeight(g, &k, r, gamma_r_sigma2, sigma_s2);  // line 12
    if (accumulate(k.begin(), k.end(), 0.f) < 10.f) {     // line 13
      for (float& v : k) v *= v;  // Square the weights   // line 14
      for (int row = 0; row < s; ++row) {                 // line 15-16
        for (int col = 0; col < s; ++col) {
          for (int chan = 0; chan < output.channels(); ++chan) {
            output.val(col + pc, row + pr, chan) +=
                (g.val(col, row, chan) + reg_plane[chan].first * (row - r) +
                reg_plane[chan].second * (col - r)) * k.val(col, row);
          }
          weights.val(col + pc, row + pr) += k.val(col, row);
        }
      }
    } else {
      ModifyPatch(y, k, &y_m, yt);  // line 18
      ModifyPatch(g, k, &g_m);  // line 19
      y_m.ToFreq();  // line 20
      g_m.ToFreq();  // line 21
      float sigma_f2 = 0.f;
      for (int row = 0; row < s; ++row) {
        for (int col = 0; col < s; ++col) {
          sigma_f2 += k.val(col, row) * k.val(col, row);
        }
      }
      sigma_f2 *= sigma2;  // line 22
      for (int row = 0; row < y_m.frows(); ++row) {
        for (int col = 0; col < y_m.fcolumns(); ++col) {
          for (int chan = 0; chan < y_m.channels(); ++chan) {
            if (row || col) {
              float x = norm(g_m.freq(col, row, chan)) / sigma_f2;
              float K;
              K = utils::fastexp(-.8f / x);  // line 23
              y_m.freq(col, row, chan) *= K;
            }
          }
        }
      }
      y_m.ToSpace();  // line 24

      // lines 25,26,30
      // col and row are the "internal" indexes (with respect to the patch).
      for (int row = 0; row < s; ++row) {
        for (int col = 0; col < s; ++col) {
          for (int chan = 0; chan < output.channels(); ++chan) {
            float pij = (row - r) * reg_plane[chan].first +
                        (col - r) * reg_plane[chan].second;
            float kij = k.val(col, row);
            output.val(col + pc, row + pr, chan) +=
                (y_m.space(col, row, chan) - (1.f - kij) * yt[chan] +
                pij * kij) * kij;
          }
          k.val(col, row) *= k.val(col, row);  // line 27
          weights.val(col + pc, row + pr) += k.val(col, row);
        }
      }
    }
    agg_weights.IncreaseWeights(k, pr - r, pc - r);  // line 29
  }
  return {move(output), move(weights)};
}

}  // namespace

Image DA3D(int &retval, const Image &noisy, const Image &guide, float sigma,
           int nthreads, int r, float sigma_s, float gamma_r,
           float threshold) {
  // padding and color transformation
  const int s = utils::NextPowerOf2(2 * r + 1);

#ifdef _OPENMP
//  if (!nthreads) nthreads = omp_get_max_threads();  // number of threads
  if (!nthreads) nthreads = com.max_thread;  // obey Siril max number of threads
#else
  nthreads = 1;
#endif  // _OPENMP

  pair<int, int> tiling = ComputeTiling(guide.rows(), guide.columns(),
                                        nthreads);
  vector<Image> noisy_tiles = SplitTiles(ColorTransform(noisy.copy()), r,
                                         s - r - 1, tiling);
  vector<Image> guide_tiles = SplitTiles(ColorTransform(guide.copy()), r,
                                         s - r - 1, tiling);
  vector<pair<Image, Image>> result_tiles(nthreads);

  unsigned thread;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif  // _OPENMP
  for (int i = 0; i < nthreads; ++i) {
#ifdef _OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif

    result_tiles[i] = DA3D_block(retval, noisy_tiles[i], guide_tiles[i], sigma,
                                 r, sigma_s, gamma_r, threshold, nthreads, thread);
  }
  return ColorTransformInverse(MergeTiles(result_tiles, guide.shape(), r,
                                          s - r - 1, tiling));
}

}  // namespace da3d
