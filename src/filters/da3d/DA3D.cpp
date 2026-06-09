/*
 * DA3D.cpp
 *
 *  Created on: 24/mar/2015
 *      Author: nicola pierazzo <nicola.pierazzo@cmla.ens-cachan.fr>
 *
 *  Migrated onto the shared img_t image type (branch harmonize_img_t):
 *  da3d::Image -> img_t<float> (planar), DftPatch -> imgops::dft_patch,
 *  and tiling/colour transforms routed through imgops. img_t is planar, which
 *  matches Siril's fits layout, so colour images are now processed correctly
 *  (the old da3d::Image read planar input as interleaved - a latent bug).
 *
 *  Coordinate mapping from the old da3d::Image API:
 *    Image(rows, cols, ch)  -> img_t<float>(cols, rows, ch)   (w=cols, h=rows)
 *    img.rows()    -> img.h        img.columns() -> img.w
 *    img.channels()-> img.d        img.pixels()  -> img.w*img.h
 *    img.val(col,row,chan) -> img(col,row,chan)
 */

#include <cmath>
#include <tuple>
#include <utility>
#include <algorithm>
#include <numeric>
#include "DA3D.hpp"
#include "WeightMap.hpp"
#include "Utils.hpp"

#include "algos/img_t/image.hpp"
#include "algos/img_t/image_ops.hpp"
#include "algos/img_t/image_dft.hpp"

#include "core/processing.h"
#include "core/gui_iface.h"

#ifdef _OPENMP
#include <omp.h>
#endif  // OPENMP

using std::max;
using std::min;
using std::vector;
using std::pair;
using std::tie;
using std::move;
using std::abs;
using std::accumulate;
using std::norm;

namespace da3d {

namespace {

void ExtractPatch(const img_t<float> &src, int pr, int pc, img_t<float> *dst) {
  // src is padded, so (pr, pc) becomes the upper left pixel (row pr, col pc).
  for (int chan = 0; chan < dst->d; ++chan)
    for (int row = 0; row < dst->h; ++row)
      for (int col = 0; col < dst->w; ++col)
        (*dst)(col, row, chan) = src(pc + col, pr + row, chan);
}

void BilateralWeight(const img_t<float> &g, img_t<float> *k, int r, float gamma_r_sigma2,
                     float sigma_s2) {
  const float inv_gamma_r_sigma2 = 1.f / gamma_r_sigma2;
  const float inv_2sigma_s2      = 0.5f / sigma_s2;
  const int channels = g.d;
  // Cache center-pixel values to avoid repeated accesses in the inner loop
  float center[4];  // channels is at most 4
  for (int chan = 0; chan < channels; ++chan)
    center[chan] = g(r, r, chan);
  for (int row = 0; row < g.h; ++row) {
    const float dr = static_cast<float>(row - r);
    const float row_dist2 = dr * dr * inv_2sigma_s2;
    for (int col = 0; col < g.w; ++col) {
      float x = 0.f;
      for (int chan = 0; chan < channels; ++chan) {
        float y = g(col, row, chan) - center[chan];
        x += y * y;
      }
      x *= inv_gamma_r_sigma2;
      const float dc = static_cast<float>(col - r);
      x += row_dist2 + dc * dc * inv_2sigma_s2;
      (*k)(col, row) = utils::fastexp(-x);
    }
  }
}

void ComputeRegressionPlane(const img_t<float> &y, const img_t<float> &g, const img_t<float> &k,
                            int r, vector<pair<float, float>> *reg_plane) {
  constexpr float epsilon = 0.0001f;
  float a = 0.f, b = 0.f, c = 0.f;
  for (int row = 0; row < y.h; ++row) {
    for (int col = 0; col < y.w; ++col) {
      a += (row - r) * (row - r) * k(col, row);
      b += (row - r) * (col - r) * k(col, row);
      c += (col - r) * (col - r) * k(col, row);
    }
  }
  float det = a * c - b * b;
  if (abs(det) < epsilon) {
    for (int chan = 0; chan < y.d; ++chan) {
      (*reg_plane)[chan] = {0.f, 0.f};
    }
  } else {
    for (int chan = 0; chan < y.d; ++chan) {
      float d = 0.f, e = 0.f;
      float central = g(r, r, chan);
      for (int row = 0; row < y.h; ++row) {
        for (int col = 0; col < y.w; ++col) {
          d += (row - r) * (y(col, row, chan) - central) * k(col, row);
          e += (col - r) * (y(col, row, chan) - central) * k(col, row);
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

void SubtractPlane(int r, const vector<pair<float, float>> &reg_plane, img_t<float> *y) {
  for (int row = 0; row < y->h; ++row) {
    for (int col = 0; col < y->w; ++col) {
      for (int chan = 0; chan < y->d; ++chan) {
        (*y)(col, row, chan) -= reg_plane[chan].first * (row - r) +
                                reg_plane[chan].second * (col - r);
      }
    }
  }
}

void ModifyPatch(const img_t<float> &patch, const img_t<float> &k, imgops::dft_patch *modified,
                 float *average = nullptr) {
  // compute the total weight of the mask
  float weight = accumulate(k.data.begin(), k.data.end(), 0.f);
  if (!weight)
    return;

  for (int chan = 0; chan < patch.d; ++chan) {
    float avg = 0.f;
    for (int row = 0; row < patch.h; ++row) {
      for (int col = 0; col < patch.w; ++col) {
        avg += k(col, row) * patch(col, row, chan);
      }
    }
    avg /= weight;
    for (int row = 0; row < patch.h; ++row) {
      for (int col = 0; col < patch.w; ++col) {
        modified->space(col, row, chan) =
            k(col, row) * patch(col, row, chan) +
            (1.f - k(col, row)) * avg;
      }
    }
    if (average) average[chan] = avg;
  }
}

pair<img_t<float>, img_t<float>> DA3D_block(int &retval, const img_t<float> &noisy,
                              const img_t<float> &guide,
                              float sigma, int r, float sigma_s,
                              float gamma_r, float threshold, unsigned nThreads, unsigned thread) {
  // useful values
  const int s = utils::NextPowerOf2(2 * r + 1);
  const float sigma2 = sigma * sigma;
  const float gamma_r_sigma2 = gamma_r * sigma2;
  const float sigma_s2 = sigma_s * sigma_s;

  // Work estimate for progress bar. Empirical based on a sample of images, may need tuning
  const unsigned predicted_loops = noisy.w * noisy.h * 5 / 1000;
  unsigned loop = 0;
  unsigned lastupdate = 0;

  // regression parameters
  const float gamma_rr_sigma2 = gamma_r_sigma2 * 10.f;
  const float sigma_sr2 = sigma_s2 * 2.f;

  // declaration of internal variables (img_t: w=cols, h=rows)
  img_t<float> y(s, s, guide.d);
  img_t<float> g(s, s, guide.d);
  img_t<float> k_reg(s, s);
  img_t<float> k(s, s);
  imgops::dft_patch y_m(s, s, guide.d);
  imgops::dft_patch g_m(s, s, guide.d);
  int pr, pc;  // coordinates of the central pixel
  vector<pair<float, float>> reg_plane(guide.d);  // parameters of the regression plane
  vector<float> ytv(guide.d);
  float* yt = ytv.data();  // weighted average of the patch
  WeightMap agg_weights(guide.h - s + 1, guide.w - s + 1);  // line 1

  img_t<float> output(guide.w, guide.h, guide.d);
  img_t<float> weights(guide.w, guide.h);
  output.set_value(0.f);
  weights.set_value(0.f);

  // main loop
  while (agg_weights.Minimum() < threshold) {  // line 4
    if (!processing_should_continue()) {
      retval++;
      break;
    }
    if (thread == 0) {
      loop++;
      if (loop - lastupdate > (predicted_loops >> 3)) {
        lastupdate = loop;
        if ((double) loop / (nThreads * predicted_loops) < 1.0)
          gui_iface.set_progress((double) loop / (nThreads * predicted_loops), _("DA3D denoising..."));
        else
          gui_iface.set_progress(PROGRESS_PULSATE, "DA3D denoising...");
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
    if (accumulate(k.data.begin(), k.data.end(), 0.f) < 10.f) {  // line 13
      for (float& v : k.data) v *= v;  // Square the weights   // line 14
      for (int row = 0; row < s; ++row) {                 // line 15-16
        for (int col = 0; col < s; ++col) {
          for (int chan = 0; chan < output.d; ++chan) {
            output(col + pc, row + pr, chan) +=
                (g(col, row, chan) + reg_plane[chan].first * (row - r) +
                reg_plane[chan].second * (col - r)) * k(col, row);
          }
          weights(col + pc, row + pr) += k(col, row);
        }
      }
    } else {
      ModifyPatch(y, k, &y_m, yt);  // line 18
      ModifyPatch(g, k, &g_m);  // line 19
      y_m.to_freq();  // line 20
      g_m.to_freq();  // line 21
      float sigma_f2 = 0.f;
      for (int row = 0; row < s; ++row) {
        for (int col = 0; col < s; ++col) {
          sigma_f2 += k(col, row) * k(col, row);
        }
      }
      sigma_f2 *= sigma2;  // line 22
      for (int row = 0; row < y_m.h(); ++row) {
        for (int col = 0; col < y_m.fw(); ++col) {
          for (int chan = 0; chan < y_m.d(); ++chan) {
            if (row || col) {
              float x = norm(g_m.freq(col, row, chan)) / sigma_f2;
              float K;
              K = utils::fastexp(-.8f / x);  // line 23
              y_m.freq(col, row, chan) *= K;
            }
          }
        }
      }
      y_m.to_space();  // line 24

      // lines 25,26,30
      // col and row are the "internal" indexes (with respect to the patch).
      for (int row = 0; row < s; ++row) {
        for (int col = 0; col < s; ++col) {
          for (int chan = 0; chan < output.d; ++chan) {
            float pij = (row - r) * reg_plane[chan].first +
                        (col - r) * reg_plane[chan].second;
            float kij = k(col, row);
            output(col + pc, row + pr, chan) +=
                (y_m.space(col, row, chan) - (1.f - kij) * yt[chan] +
                pij * kij) * kij;
          }
          k(col, row) *= k(col, row);  // line 27
          weights(col + pc, row + pr) += k(col, row);
        }
      }
    }
    agg_weights.IncreaseWeights(k, pr - r, pc - r);  // line 29
  }
  return {std::move(output), std::move(weights)};
}

}  // namespace

img_t<float> DA3D(int &retval, const img_t<float> &noisy, const img_t<float> &guide, float sigma,
           int nthreads, int r, float sigma_s, float gamma_r,
           float threshold) {
  // padding and color transformation
  const int s = utils::NextPowerOf2(2 * r + 1);

#ifdef _OPENMP
  if (!nthreads) nthreads = com.max_thread;  // obey Siril max number of threads
#else
  nthreads = 1;
#endif  // _OPENMP

  imgops::tiling_t tiling = imgops::compute_tiling(guide.h, guide.w, nthreads);

  img_t<float> noisy_t = noisy;
  imgops::color_transform(noisy_t, true);
  img_t<float> guide_t = guide;
  imgops::color_transform(guide_t, true);
  vector<img_t<float>> noisy_tiles = imgops::split_tiles(noisy_t, r, s - r - 1, tiling);
  vector<img_t<float>> guide_tiles = imgops::split_tiles(guide_t, r, s - r - 1, tiling);
  vector<pair<img_t<float>, img_t<float>>> result_tiles(nthreads);

  unsigned thread;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) private(thread)
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

  img_t<float> merged = imgops::merge_tiles(result_tiles, guide.w, guide.h, r, s - r - 1, tiling);
  imgops::color_transform(merged, false);
  return merged;
}

}  // namespace da3d
