/*
 * DftPatch.hpp
 *
 *  Created on: 11/feb/2015
 *      Author: nicola pierazzo <nicola.pierazzo@cmla.ens-cachan.fr>
 */

#ifndef DA3D_DFTPATCH_HPP_
#define DA3D_DFTPATCH_HPP_

#include <fftw3.h>
#include <cassert>
#include <complex>

namespace da3d {

class DftPatch {
 public:
  DftPatch(int rows, int columns, int channels = 1);
  ~DftPatch();
  void ToFreq();
  void ToSpace();
  int rows() const { return rows_; }
  int columns() const { return columns_; }
  int frows() const { return rows_; }
  int fcolumns() const { return fcolumns_; }
  int channels() const { return channels_; }
  float& space(int col, int row, int chan = 0);
  std::complex<float>& freq(int col, int row, int chan = 0);

 private:
  float *space_;
  std::complex<float> *freq_;
  fftwf_plan plan_forward_;
  fftwf_plan plan_backward_;
  int rows_, columns_, fcolumns_, channels_;
};

inline float& DftPatch::space(int col, int row, int chan) {
  assert(0 <= col && col < columns_);
  assert(0 <= row && row < rows_);
  assert(0 <= chan && chan < channels_);
  return space_[row * columns_ * channels_ + col * channels_ + chan];
}

inline std::complex<float>& DftPatch::freq(int col, int row, int chan) {
  assert(0 <= col && col < fcolumns_);
  assert(0 <= row && row < rows_);
  assert(0 <= chan && chan < channels_);
  return freq_[row * fcolumns_ * channels_ + col * channels_ + chan];
}

inline DftPatch::DftPatch(int rows, int columns, int channels)
    : rows_(rows), columns_(columns), fcolumns_(columns / 2 + 1), channels_(channels) {
  int N = rows * columns * channels;
  int N_half = rows * fcolumns_ * channels;
  space_ = reinterpret_cast<float *>(fftwf_malloc(sizeof(float) * N));
  freq_ = reinterpret_cast<std::complex<float> *>(fftwf_malloc(
      sizeof(fftwf_complex) * N_half));
  int n[] = {rows, columns};
#pragma omp critical
  {
    plan_forward_ = fftwf_plan_many_dft_r2c(2, n, channels, space_, NULL,
                                            channels, 1,
                                            reinterpret_cast<fftwf_complex *>(freq_),
                                            NULL, channels, 1, FFTW_MEASURE);
    plan_backward_ = fftwf_plan_many_dft_c2r(2, n, channels,
                                             reinterpret_cast<fftwf_complex *>(freq_),
                                             NULL, channels, 1, space_, NULL,
                                             channels, 1, FFTW_MEASURE);
  }
}

inline DftPatch::~DftPatch() {
  fftwf_free(space_);
  fftwf_free(freq_);
  fftwf_destroy_plan(plan_forward_);
  fftwf_destroy_plan(plan_backward_);
}

inline void DftPatch::ToFreq() {
  fftwf_execute(plan_forward_);
}

inline void DftPatch::ToSpace() {
  fftwf_execute(plan_backward_);
  for (int i = 0; i < rows_ * columns_ * channels_; ++i) {
    space_[i] /= rows_ * columns_;
  }
}

} /* namespace da3d */

#endif  // DA3D_DFTPATCH_HPP_
