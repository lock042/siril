/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Per-patch real<->complex DFT helper for patch-based denoisers (da3d),
 * backed by img_t storage (branch harmonize_img_t).
 *
 * Why this exists rather than reusing img_t's fft.hpp: fft::r2c there runs a
 * FULL complex FFT (building a complex copy of the input) and allocates a new
 * buffer on every call. A patch denoiser transforms thousands of small patches
 * in a hot loop, so it needs (a) the real-optimized half-spectrum transforms
 * and (b) persistent buffers + reused plans. This mirrors da3d::DftPatch but
 * uses planar img_t storage, which is exactly FFTW's 2D row-major layout
 * (within a channel, index = y*w + x), and exposes the buffers as img_t so the
 * rest of the pipeline stays in the shared type.
 */

#pragma once

#include <complex>

#include "image.hpp"

namespace imgops {

class dft_patch {
public:
    //! Allocate a patch of w x h with `d` channels and build the r2c/c2r plans.
    //! Plan construction (FFTW_MEASURE) scribbles on the buffers, so create the
    //! patch first, then load data.
    dft_patch(int w, int h, int d = 1)
        : space_(w, h, d), freq_(w / 2 + 1, h, d),
          w_(w), h_(h), d_(d), fw_(w / 2 + 1) {
        float* in = space_.data.data();
        fftwf_complex* out = reinterpret_cast<fftwf_complex*>(freq_.data.data());
        int n[2] = {h_, w_};
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
        {
            // planar batch over channels: istride/ostride 1, dist = channel size
            plan_forward_ = fftwf_plan_many_dft_r2c(2, n, d_, in, nullptr, 1, w_ * h_,
                                                    out, nullptr, 1, fw_ * h_,
                                                    FFTW_MEASURE);
            plan_backward_ = fftwf_plan_many_dft_c2r(2, n, d_, out, nullptr, 1, fw_ * h_,
                                                     in, nullptr, 1, w_ * h_,
                                                     FFTW_MEASURE);
        }
    }

    ~dft_patch() {
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
        {
            if (plan_forward_) fftwf_destroy_plan(plan_forward_);
            if (plan_backward_) fftwf_destroy_plan(plan_backward_);
        }
    }

    dft_patch(const dft_patch&) = delete;
    dft_patch& operator=(const dft_patch&) = delete;

    //! Real-space buffer (w x h x d). Write input here, read output after to_space().
    img_t<float>& space() { return space_; }
    const img_t<float>& space() const { return space_; }
    float& space(int x, int y, int c = 0) { return space_(x, y, c); }

    //! Half-spectrum frequency buffer ((w/2+1) x h x d).
    img_t<std::complex<float>>& freq() { return freq_; }
    const img_t<std::complex<float>>& freq() const { return freq_; }
    std::complex<float>& freq(int x, int y, int c = 0) { return freq_(x, y, c); }

    int w() const { return w_; }
    int h() const { return h_; }
    int d() const { return d_; }
    int fw() const { return fw_; }  // half-spectrum width = w/2 + 1

    //! Forward real-to-complex transform (space -> freq).
    void to_freq() { fftwf_execute(plan_forward_); }

    //! Inverse complex-to-real transform (freq -> space), normalized so that
    //! to_freq() followed by to_space() is the identity.
    void to_space() {
        fftwf_execute(plan_backward_);
        const float inv_n = 1.f / (static_cast<float>(w_) * h_);
        for (float& v : space_.data)
            v *= inv_n;
    }

private:
    img_t<float> space_;
    img_t<std::complex<float>> freq_;
    fftwf_plan plan_forward_ = nullptr;
    fftwf_plan plan_backward_ = nullptr;
    int w_, h_, d_, fw_;
};

}  // namespace imgops
