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

#include <memory>

#ifndef IMG_NO_FFTW
#include <fftw3.h>
#endif

#include "image.hpp"
#include "image_expr.hpp"

#include <map>
#include <unordered_map>

struct dim_t {
    int h, w;
    int d;

    bool operator==(const dim_t &other) const {
        return h == other.h && w == other.w && d == other.d;
    }
};

namespace std {
    template <>
    struct hash<dim_t>
    {
        std::size_t operator()(const dim_t& k) const
        {
            return ((std::hash<int>()(k.h)
                  ^ (std::hash<int>()(k.w) << 1)) >> 1)
                  ^ (std::hash<int>()(k.d) << 1);
        }
    };
}

template <typename T>
struct plan_t {
};

template <>
struct plan_t<float> {
    typedef fftwf_plan plan_type;
    typedef fftwf_complex value_type;
    plan_type plan_forward = nullptr;
    plan_type plan_backward = nullptr;

    plan_t(plan_t<float>&& p) {
        std::swap(plan_forward, p.plan_forward);
        std::swap(plan_backward, p.plan_backward);
    }

    plan_t(dim_t dim, int flags) {
        img_t<std::complex<float>> img(dim.w, dim.h, dim.d);
        auto out = reinterpret_cast<value_type*>(&img[0]);
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
        {
            plan_forward = fftwf_plan_many_dft(2, &dim.h, dim.d, out, &dim.h, dim.d, 1, out,
                                               &dim.h, dim.d, 1, FFTW_FORWARD, flags);
            plan_backward = fftwf_plan_many_dft(2, &dim.h, dim.d, out, &dim.h, dim.d, 1, out,
                                                &dim.h, dim.d, 1, FFTW_BACKWARD, flags);
        }
    }

    ~plan_t() {
        return; // don't free the plan...
    }

    void execute_forward(std::complex<float>* _out) const {
        auto out = reinterpret_cast<value_type*>(_out);
        fftwf_execute_dft(plan_forward, out, out);
    }

    void execute_backward(std::complex<float>* _out) const {
        auto out = reinterpret_cast<value_type*>(_out);
        fftwf_execute_dft(plan_backward, out, out);
    }

};

template <typename T>
inline plan_t<T>* make_plan(dim_t dim, int flags)
{
    static std::unordered_map<dim_t, plan_t<T>> cache;
    auto [it, inserted] = cache.emplace(dim, plan_t<T>(dim, flags));
    return &it->second;
}

namespace fft {

    template <typename E>
    auto c2c(const E& o, bool fast=true) {
        using T = typename E::value_type::value_type;
        dim_t dim = {.h=o.h, .w=o.w, .d=o.d};
        const auto plan = make_plan<T>(dim, fast ? FFTW_ESTIMATE : com.pref.fftw_conf.strategy);
        auto tmp = to_img(o);
        plan->execute_forward(&tmp[0]);
        return tmp;
    }

    template <typename E>
    auto c2r(const E& o, bool fast=true) {
        return to_img(std::real(c2c(to_img(o), fast)));
    }

    template <typename E>
    auto r2c(const E& o, bool fast=true) {
        return c2c(to_img<std::complex<typename E::value_type>>(o), fast);
    }

    template <typename E>
    auto r2r(const E& o, bool fast=true) {
        return c2r(to_img<std::complex<typename E::value_type>>(o), fast);
    }

    template <typename T>
    auto shift(const img_t<T>& in) {
        img_t<T> out;
        out.resize(in);

        int halfw = (in.w + 1) / 2.;
        int halfh = (in.h + 1) / 2.;
        int ohalfw = in.w - halfw;
        int ohalfh = in.h - halfh;
        for (int l = 0; l < in.d; l++) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.fftw_max_thread) if(in.size > 40000 && com.fftw_max_thread > 1)
{
#pragma omp for simd schedule(static) collapse(2)
#endif
            for (int y = 0; y < halfh; y++) {
                for (int x = 0; x < ohalfw; x++) {
                    out(x, y + ohalfh, l) = in(x + halfw, y, l);
                }
            }
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
            for (int y = 0; y < halfh; y++) {
                for (int x = 0; x < halfw; x++) {
                    out(x + ohalfw, y + ohalfh, l) = in(x, y, l);
                }
            }
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
            for (int y = 0; y < ohalfh; y++) {
                for (int x = 0; x < ohalfw; x++) {
                    out(x, y, l) = in(x + halfw, y + halfh, l);
                }
            }
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
            for (int y = 0; y < ohalfh; y++) {
                for (int x = 0; x < halfw; x++) {
                    out(x + ohalfw, y, l) = in(x, y + halfh, l);
                }
            }
#ifdef _OPENMP
}
#endif
        }

        return out;
    }

    inline const int* get_optimal_table(int& _lut_size) {
        // based on the matlab code of Sunghyun Cho
        const int lut_size = 16383;
        static int is_optimal[lut_size] = {0};
        if (!is_optimal[1]) {
            for (int e2 = 1; e2 < lut_size; e2 *= 2)
            for (int e3 = e2; e3 < lut_size; e3 *= 3)
            for (int e5 = e3; e5 < lut_size; e5 *= 5)
            for (int e7 = e5; e7 < lut_size; e7 *= 7) {
                is_optimal[e7] = true;
                if (e7 * 11 < lut_size)
                    is_optimal[e7*11] = true;
                if (e7 * 13 < lut_size)
                    is_optimal[e7*13] = true;
            }
        }
        _lut_size = lut_size;
        return is_optimal;
    }

    inline int get_optimal_size_up(int size) {
        int lut_size;
        const int* is_optimal = get_optimal_table(lut_size);
        for (int i = size; i < lut_size; i++) {
            if (is_optimal[i]) {
                return i;
            }
        }
        return size;
    }

    inline int get_optimal_size_down(int size) {
        int lut_size;
        const int* is_optimal = get_optimal_table(lut_size);
        for (int i = size; i > 0; i--) {
            if (is_optimal[i]) {
                return i;
            }
        }
        return size;
    }

    template <typename T>
    void psf2otf(img_t<std::complex<T>>& out, const img_t<T>& k, int w, int h, int d=1)
    {
        out.resize(w, h, d);
        out.padcirc(k);
        out = c2c(out);
    }

}

namespace ifft {

    template <typename E>
    img_t<typename E::value_type> c2c(const E& o, bool fast=true) {
        using T = typename E::value_type::value_type;
        dim_t dim = {.h=o.h, .w=o.w, .d=o.d};
        const auto plan = make_plan<T>(dim, fast ? FFTW_ESTIMATE : com.pref.fftw_conf.strategy);
        auto tmp = to_img(o);
        plan->execute_backward(&tmp[0]);

        T norm = tmp.w * tmp.h;
        tmp.map(tmp / norm);
        return tmp;
    }

    template <typename E>
    auto c2r(const E& o, bool fast=true) {
        return to_img(img::real(c2c(to_img(o), fast)));
    }

    template <typename E>
    auto r2c(const E& o, bool fast=true) {
        return c2c(to_img<std::complex<typename E::value_type>>(o), fast);
    }

    template <typename E>
    auto r2r(const E& o, bool fast=true) {
        return c2r(to_img<std::complex<typename E::value_type>>(o), fast);
    }

    template <typename T>
    auto shift(const img_t<T>& in) {
        img_t<T> out;
        out.resize(in);

        int halfw = (in.w + 1) / 2.;
        int halfh = (in.h + 1) / 2.;
        int ohalfw = in.w - halfw;
        int ohalfh = in.h - halfh;
        for (int l = 0; l < in.d; l++) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) if(in.size > 40000 && com.max_thread > 1)
{
#pragma omp for simd schedule(static) collapse(2)
#endif
            for (int y = 0; y < ohalfh; y++) {
                for (int x = 0; x < halfw; x++) {
                    out(x, y + halfh, l) = in(x + ohalfw, y, l);
                }
            }
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
            for (int y = 0; y < ohalfh; y++) {
                for (int x = 0; x < ohalfw; x++) {
                    out(x + halfw, y + halfh, l) = in(x, y, l);
                }
            }
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
            for (int y = 0; y < halfh; y++) {
                for (int x = 0; x < halfw; x++) {
                    out(x, y, l) = in(x + ohalfw, y + ohalfh, l);
                }
            }
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
            for (int y = 0; y < halfh; y++) {
                for (int x = 0; x < ohalfw; x++) {
                    out(x + halfw, y, l) = in(x, y + ohalfh, l);
                }
            }
#ifdef _OPENMP
}
#endif
        }
        return out;
    }
}

