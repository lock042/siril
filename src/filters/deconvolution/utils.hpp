/*
MIT License

Copyright (c) 2018 Jérémy Anger, Gabriele Facciolo, Mauricio Delbracio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#pragma once

#include "core/siril.h"
#include "core/arithm.h"
#include "image.hpp"
#include "image_expr.hpp"
#include "labeling.hpp"
#include "vec2.hpp"

namespace utils {

    inline void blur(img_t<float>& out, const img_t<float>& in, float sigma) {
        assert(in.d == 1);
        out.resize(in.w, in.h);
        gaussblur(&out[0], (float*) &in[0], in.w, in.h, sigma);
    }

    inline void upsample(img_t<float>& out, const img_t<float>& _in, float factor,
                         int targetw, int targeth) {
        img_t<float> in = _in; // copy input
        out.resize(targetw, targeth, in.d);
        magnify(&out[0], &in[0], out.w, out.h, out.d, in.w, in.h, factor);
    }

    template <typename T>
    void transpose(img_t<T>& out, const img_t<T>& in)
    {
        if (&in == &out) {
            auto copy = in;
            return transpose(out, copy);
        }

        out.resize(in.h, in.w, in.d);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(3) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int d = 0; d < in.d; d++) {
            for (int y = 0; y < in.h; y++) {
                for (int x = 0; x < in.w; x++) {
                    out(y, x, d) = in(x, y, d);
                }
            }
        }
    }

    template <typename T>
    img_t<T> add_padding(const img_t<T>& _f, int hw, int hh)
    {
        img_t<T> f(_f.w + hw*2, _f.h + hh*2, _f.d);
        f.set_value(T(0));
        slice(f, _sl(hw, -hw-1), _sl(hh, -hh-1)).map(_f);
        // replicate borders
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
{
        #pragma omp for collapse(3) schedule(static)
#endif
        for (int y = 0; y < hh; y++) {
            for (int x = 0; x < f.w; x++) {
                for (int l = 0; l < f.d; l++) {
                    f(x, y, l) = f(x, 2*hh - y, l);
                    f(x, f.h-1-y, l) = f(x, f.h-1-2*hh+y, l);
                }
            }
        }
#ifdef _OPENMP
        #pragma omp for collapse(3) schedule(static)
#endif
        for (int y = 0; y < f.h; y++) {
            for (int x = 0; x < hw; x++) {
                for (int l = 0; l < f.d; l++) {
                    f(x, y, l) = f(2*hw - x, y, l);
                    f(f.w-1-x, y, l) = f(f.w-1-2*hw+x, y, l);
                }
            }
        }
#ifdef _OPENMP
}
#endif
        return f;
    }

    template <typename T>
    img_t<T> add_padding(const img_t<T>& f, const img_t<T>& K)
    {
        return add_padding(f, K.w/2, K.h/2);
    }

    template <typename T>
    img_t<T> remove_padding(const img_t<T>& f, int hw, int hh)
    {
        return to_img(slice(f, _sl(hw, -hw-1), _sl(hh, -hh-1)));
    }

    template <typename T>
    img_t<T> remove_padding(const img_t<T>& f, const img_t<T>& K)
    {
        return remove_padding(f, K.w/2, K.h/2);
    }

    template <typename T>
    void center_kernel(img_t<T>& kernel) {
        T dx = 0.f;
        T dy = 0.f;
        T sum = kernel.sum();
        if (!sum)
            return;
    #ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
            T local_dx = 0.f;
            T local_dy = 0.f;

            #pragma omp for schedule(static) nowait
            for (int y = 0; y < kernel.h; y++) {
                for (int x = 0; x < kernel.w; x++) {
                    local_dx += kernel(x, y) * x;
                    local_dy += kernel(x, y) * y;
                }
            }

            #pragma omp atomic
            dx += local_dx;
            #pragma omp atomic
            dy += local_dy;

            #pragma omp barrier

            #pragma omp single
            {
                dx = std::round(dx / sum);
                dy = std::round(dy / sum);
            }

            img_t<T> copy(kernel);

            #pragma omp for collapse(2) schedule(static)
            for (int y = 0; y < kernel.h; y++) {
                for (int x = 0; x < kernel.w; x++) {
                    int nx = x + static_cast<int>(dx) - kernel.w / 2;
                    int ny = y + static_cast<int>(dy) - kernel.h / 2;
                    if (nx >= 0 && nx < kernel.w && ny >= 0 && ny < kernel.h) {
                        kernel(x, y) = copy(nx, ny);
                    }
                }
            }
        }
    #else
        for (int y = 0; y < kernel.h; y++) {
            for (int x = 0; x < kernel.w; x++) {
                dx += kernel(x, y) * x;
                dy += kernel(x, y) * y;
            }
        }

        dx = std::round(dx / sum);
        dy = std::round(dy / sum);

        img_t<T> copy(kernel);
        kernel.set_value(0);

        for (int y = 0; y < kernel.h; y++) {
            for (int x = 0; x < kernel.w; x++) {
                int nx = x + static_cast<int>(dx) - kernel.w / 2;
                int ny = y + static_cast<int>(dy) - kernel.h / 2;
                if (nx >= 0 && nx < kernel.w && ny >= 0 && ny < kernel.h) {
                    kernel(x, y) = copy(nx, ny);
                }
            }
        }
    #endif
    }

    template <typename T>
    void circular_gradients(vec2<img_t<T>>& out, const img_t<T>& in) {
        out[0].resize(in);
        out[1].resize(in);

        int w = in.w;
        int h = in.h;
        int d = in.d;
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(3) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int l = 0; l < d; l++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    out[0](x, y, l) = in((x+1)%w, y, l) - in(x, y, l);
                    out[1](x, y, l) = in(x, (y+1)%h, l) - in(x, y, l);
                }
            }
        }
    }

    template <typename T>
    void circular_divergence(img_t<T>& out, const vec2<img_t<T>>& in) {
        out.resize(in[0]);

        int w = out.w;
        int h = out.h;
        int d = out.d;
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(3) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int l = 0; l < d; l++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    out(x, y, l) = in[0](x, y, l) - in[0]((x-1+w)%w, y, l)
                                 + in[1](x, y, l) - in[1](x, (y-1+h)%h, l);
                }
            }
        }
    }

    template <typename T>
    void downsa2(img_t<T>& out, const img_t<T>& in) {
        if (out.size == 0)
            out.resize(in.w/2, in.h/2, in.d);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(3) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int d = 0; d < out.d; d++) {
            for (int j = 0; j < out.h; j++) {
                for (int i = 0; i < out.w; i++) {
                    T m = getpixel_1(in, 2*i, 2*j, d)
                        + getpixel_1(in, 2*i+1, 2*j, d)
                        + getpixel_1(in, 2*i, 2*j+1, d)
                        + getpixel_1(in, 2*i+1, 2*j+1, d);
                    out(i, j, d) = m / T(4);
                }
            }
        }
    }

    template <typename T>
    void upsa2(img_t<T>& out, const img_t<T>& in) {
        if (out.size == 0)
            out.resize(in.w*2, in.h*2, in.d);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(3) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int d = 0; d < out.d; d++) {
            for (int j = 0; j < out.h; j++) {
                for (int i = 0; i < out.w; i++) {
                    float x = (i - 0.5) / 2;
                    float y = (j - 0.5) / 2;
                    out(i, j, d) = bilinear_interpolation(in, x, y, d);
                }
            }
        }
    }

    template <typename T>
    void remove_isolated_cc(img_t<T>& k) {
        T sum = k.sum();
        if (!sum) return; // avoid div/0
        img_t<int> lab;
        std::vector<T> sums;
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
            #pragma omp for schedule(static)
#endif
            // First loop: normalize k
            for (int i = 0; i < k.size; i++)
                k[i] /= sum;

            // Synchronize threads before sequential operations
#ifdef _OPENMP
            #pragma omp single
            {
#endif
                labeling::labels(lab, k);
                std::map<int, T> sum_map = labeling::sum(lab, k);
                sums.reserve(sum_map.size());
                for (const auto& pair : sum_map) {
                    sums.push_back(pair.second);
                }
#ifdef _OPENMP
            }

            #pragma omp for schedule(static)
#endif
            // Second loop: remove low-intensity components
            for (int i = 0; i < k.size; i++) {
                if (sums[lab[i]] < T(0.1))
                    k[i] = T(0);
            }
#ifdef _OPENMP
        }
#endif
    }
}
