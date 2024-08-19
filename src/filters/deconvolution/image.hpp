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

/*

Modifications and extensions for use in Siril by Adrian Knagg-Baugh (c) team
free-astro 2022-2023.

*/
#pragma once

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

// All other files get the definition of HAVE_FFTW3F_MULTITHREAD from
// siril.h. This file doesn't include siril.h so we define it again here.
#if defined (HAVE_FFTW3F_OMP) || defined (HAVE_FFTW3F_THREADS)
#define HAVE_FFTW3F_MULTITHREAD
#endif

#include <unordered_map>

#include <cassert>
#include <complex>
#include <limits>
#include <functional>
#include <vector>
#include <numeric>
#include <memory>
#include "core/siril.h"
#include <fftw3.h>
#include "fftw_allocator.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

extern int sequence_is_running;

template <typename T>
class img_t {
public:
    // Data
    typedef T value_type;
    int size, w, h, d;
    std::vector<T, fftw_alloc<T>> data;
    fftwf_plan forwardplanf = nullptr;
    fftwf_plan backwardplanf = nullptr;

    // Methods
    static void use_threading(int n) {
        fftwf_set_timelimit(com.pref.fftw_conf.timelimit);
#ifdef HAVE_FFTW3F_MULTITHREAD
        if (n > 1) {
            fftwf_plan_with_nthreads(n);
        }
#endif
        fprintf(stdout, "fftwf initialized with %d threads, planning time limit %.1f seconds\n", n, com.pref.fftw_conf.timelimit);
    }

    auto begin() {
        return data.begin();
    }

    auto end() {
        return data.end();
    }

    img_t() : size(0), w(0), h(0), d(0) {
    }

    img_t(int w, int h, int d=1)
        : size(w*h*d), w(w), h(h), d(d), data(w*d*h) {
    }

    img_t(int w, int h, int d, T* data)
        : size(w*h*d), w(w), h(h), d(d) {
        this->data.assign(data, data+w*h*d);
    }

    img_t(const img_t<T>& o)
        : size(o.size), w(o.w), h(o.h), d(o.d), data(o.data) {
    }

    img_t& operator=(img_t<T>&& o) noexcept {
        if (this != &o) {  // Check for self-assignment
            // Destroy any existing plans
            if (forwardplanf)
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
                fftwf_destroy_plan(forwardplanf);
            if (backwardplanf)
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
                fftwf_destroy_plan(backwardplanf);

            // Transfer ownership of resources
            w = o.w;
            h = o.h;
            d = o.d;
            size = o.size;
            data = std::move(o.data);
            forwardplanf = o.forwardplanf;
            backwardplanf = o.backwardplanf;

            // Leave the source object in a valid state
            o.w = 0;
            o.h = 0;
            o.d = 0;
            o.size = 0;
            o.forwardplanf = nullptr;
            o.backwardplanf = nullptr;
        }
        return *this;
    }

    img_t(img_t<T>&& o) noexcept : size(o.size), w(o.w), h(o.h), d(o.d), data(std::move(o.data)), forwardplanf(std::move(o.forwardplanf)), backwardplanf(std::move(o.backwardplanf)) {
        o.w = 0;
        o.h = 0;
        o.d = 0;
        o.size = 0;
        o.forwardplanf = nullptr;
        o.backwardplanf = nullptr;
    }

    img_t operator=(const img_t<T>& o) {
        w = o.w;
        h = o.h;
        d = o.d;
        size = w*d*h;
        forwardplanf = o.forwardplanf;
        backwardplanf = o.backwardplanf;
        data = o.data;
        return *this;
    }

    ~img_t() {
        if (forwardplanf)
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
            fftwf_destroy_plan(forwardplanf);
        if (backwardplanf)
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
            fftwf_destroy_plan(backwardplanf);
    }

    // indexing
    inline T& operator[](int i) {
        return data[i];
    }
    inline const T& operator[](int i) const {
        return data[i];
    }
    inline T& operator()(int x, int y, int dd=0) {
        return data[dd+d*(x+y*w)];
    }
    inline const T& operator()(int x, int y, int dd=0) const {
        return data[dd+d*(x+y*w)];
    }

    void set_value(const T& v) {
        std::fill(data.begin(), data.end(), v);
    }

    template <typename T2>
    void set_value(const T2& v) {
        std::fill(data.begin(), data.end(), v);
    }

    template <typename T2>
    void fill_array(T2** array) {
        size_t size = data.size();
        T2* temp = (T2*) realloc((*array), size * sizeof(T));
        if (!temp) {
            PRINT_ALLOC_ERR;
            return;
        } else {
            (*array) = temp;
            memcpy((*array), data.data(), data.size() * sizeof(T2));
        }
    }

    T sum() const {
        return fold<T>(std::plus<T>());
    }

    T max() const {
        return fold<T>([](const T& a, const T& b) { return a > b ? a : b; });
    }

    T min() const {
        return fold<T>([](const T& a, const T& b) { return a < b ? a : b; });
    }

    template <typename T2>
    T2 fold(const std::function<T2(const T&, const T2&)>& f) const {
        return std::accumulate(data.begin(), data.end(), T2(), f);
    }

    template <typename E>
    bool similar(const E& o) const {
        // should check if the cast makes sense in case T != E::value_type
        bool similar = w == o.w && h == o.h && d == o.d;
        if (!similar)
            fprintf(stderr, "%dx%dx%d (type %s) != %dx%dx%d (type %s)\n",
                    w, h, d, typeid(*this).name(), o.w, o.h, o.d, typeid(o).name());
        return similar;
    }

    void resize(int w, int h, int d=1) {
        if (this->w != w || this->h != h || this->d != d) {
            if (forwardplanf)
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
                fftwf_destroy_plan(forwardplanf);
            if (backwardplanf)
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
                fftwf_destroy_plan(backwardplanf);

            forwardplanf = nullptr;
            backwardplanf = nullptr;
            this->w = w;
            this->h = h;
            this->d = d;
            size = w * h * d;
            data.resize(w * h * d);
        }
    }

    template <typename T2>
    void resize(const img_t<T2>& o) {
        resize(o.w, o.h, o.d);
    }

    // map (no arg, xyd arg)
    template <class E>
    void map(const E& o) {
        assert(o.similar(*this));
        int n = w * h * d;
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
#pragma omp parallel for simd schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < n; i++)
            data[i] = o[i];
    }

    // naive convolution
    template <typename T2, typename T3>
    void conv2(const img_t<T2>& x, const img_t<T3>& k) {
        assert(k.w / 2 && k.w == k.h); // Ensure kernel is square and odd
        assert(x.similar(*this));

        img_t<T> out(x.w, x.h, x.d);
        const int ix = (k.w-1)/2;

        // Main convolution (excluding borders)
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
            #pragma omp for collapse(3) schedule(guided)
#endif
            for (int c = 0; c < d; c++) {
                for (int i = ix; i < x.w-ix; i++) {
                    for (int j = ix; j < x.h-ix; j++) {
                        T val = T(0);
#ifdef _OPENMP
                        #pragma omp simd reduction(+:val)
#endif
                        for (int m = -ix; m < ix+1; m++) {
                            for (int n = -ix; n < ix+1; n++) {
                                val += x(i+m, j+n, c) * k(m+ix, n+ix, 0);
                            }
                        }
                        out(i, j, c) = val;
                    }
                }
            }

            // Left border
#ifdef _OPENMP
            #pragma omp for collapse(3) schedule(guided)
#endif
            for (int c = 0; c < d; c++) {
                for (int i = 0; i < ix; i++) {
                    for (int j = 0; j < x.h; j++) {
                        T val = T(0);
#ifdef _OPENMP
                        #pragma omp simd reduction(+:val)
#endif
                        for (int m = -ix; m < ix+1; m++) {
                            for (int n = -ix; n < ix+1; n++) {
                                int im = i + m, jn = j + n;
                                if (im >= 0 && im < x.w && jn >= 0 && jn < x.h) {
                                    val += x(im, jn, c) * k(m+ix, n+ix, 0);
                                }
                            }
                        }
                        out(i, j, c) = val;
                    }
                }
            }

            // Right border
#ifdef _OPENMP
            #pragma omp for collapse(3) schedule(guided)
#endif
            for (int c = 0; c < d; c++) {
                for (int i = x.w - ix; i < x.w; i++) {
                    for (int j = 0; j < x.h; j++) {
                        T val = T(0);
#ifdef _OPENMP
                        #pragma omp simd reduction(+:val)
#endif
                        for (int m = -ix; m < ix+1; m++) {
                            for (int n = -ix; n < ix+1; n++) {
                                int im = i + m, jn = j + n;
                                if (im >= 0 && im < x.w && jn >= 0 && jn < x.h) {
                                    val += x(im, jn, c) * k(m+ix, n+ix, 0);
                                }
                            }
                        }
                        out(i, j, c) = val;
                    }
                }
            }

            // Top border (excluding corners already handled by left/right)
#ifdef _OPENMP
            #pragma omp for collapse(3) schedule(guided)
#endif
            for (int c = 0; c < d; c++) {
                for (int i = ix; i < x.w - ix; i++) {
                    for (int j = 0; j < ix; j++) {
                        T val = T(0);
#ifdef _OPENMP
                        #pragma omp simd reduction(+:val)
#endif
                        for (int m = -ix; m < ix+1; m++) {
                            for (int n = -ix; n < ix+1; n++) {
                                int im = i + m, jn = j + n;
                                if (im >= 0 && im < x.w && jn >= 0 && jn < x.h) {
                                    val += x(im, jn, c) * k(m+ix, n+ix, 0);
                                }
                            }
                        }
                        out(i, j, c) = val;
                    }
                }
            }

            // Bottom border (excluding corners already handled by left/right)
#ifdef _OPENMP
            #pragma omp for collapse(3) schedule(guided)
#endif
            for (int c = 0; c < d; c++) {
                for (int i = ix; i < x.w - ix; i++) {
                    for (int j = x.h - ix; j < x.h; j++) {
                        T val = T(0);
#ifdef _OPENMP
                        #pragma omp simd reduction(+:val)
#endif
                        for (int m = -ix; m < ix+1; m++) {
                            for (int n = -ix; n < ix+1; n++) {
                                int im = i + m, jn = j + n;
                                if (im >= 0 && im < x.w && jn >= 0 && jn < x.h) {
                                    val += x(im, jn, c) * k(m+ix, n+ix, 0);
                                }
                            }
                        }
                        out(i, j, c) = val;
                    }
                }
            }
#ifdef _OPENMP
        }
#endif

        (*this).map(out);
    }

    // Updated mapf() function, hopefully should avoid the MacOS problems
    template <typename F>
    void mapf(F&& f) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for simd schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < size; i++) {
            (*this)[i] = f((*this)[i]);
        }
    }

    // Overload for when the function returns a different type
    template <typename F, typename R = std::invoke_result_t<F, T>>
    std::enable_if_t<!std::is_same_v<R, T>, img_t<R>> mapf(F&& f) const {
        img_t<R> result(w, h, d);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < size; i++) {
            result[i] = f((*this)[i]);
        }
        return result;
    }

/*
// Keeping these functions for later work, not needed yet
    template <typename T2>
    void pasteinto(const img_t<T2>&o, int xoff, int yoff, int width, int height) { // Could this be better optimised?
        assert (width <= o.w && height <= o.h);
        assert(width + xoff <= this->w);
        assert(height + yoff <= this->h);
        for (int i = 0 ; i < width ; i++) {
            for (int j = 0 ; j < height; j++) {
                for (int c = 0 ; c < o.d ; c++) {
                (*this)(i + xoff, j + yoff, c) = o(i, j, c);
                }
            }
        }
    }

    template <typename T2>
    void masktolist(const img_t<T>& in,
                    const img_t<T2>& mask) {
        assert(in.similar(mask));
        int index = 0;
        for (int i = 0 ; i < mask.w; i++) {
            for (int j = 0 ; j < mask.h ; j++) {
                if (mask(i,j) == T2(1)) {
                    for (int c = 0 ; c < in.d ; c++) {
                        (*this)(index, 0, c) = in(i, j, c);
                        index++;
                    }
                }
            }
        }
    }

    template <typename T2>
    void expandlisttoimage(const img_t<T>& list,
                           const img_t<T2>& mask) {
        assert ((*this).similar(mask));
        (*this).set_value(T(0));
        int index = 0;
        for (int i = 0 ; i < mask.w ; i++) {
            for (int j = 0 ; j < mask.h ; j++) {
                if (mask(i,j) == T2(1)) {
                    for (int c = 0 ; c < list.d ; c++) {
                        (*this)(i, j, c) = list(index, 0, c);
                        index++;
                    }
                }
            }
        }
    }

    template <typename T2>
    void map_channel(const img_t<T2>&src, int srcchan, int destchan) {
        for (int i = 0 ; i < src.w ; i++) {
            for (int j = 0 ; j < src.h ; j++) {
                (*this)(i, j, destchan) = src(i, j, srcchan);
            }
        }
    }
*/

    template <typename T2, class F>
    void map(const img_t<T2>& o, const F& f) {
        assert(o.similar(*this));
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < size; i++) {
            (*this)[i] = f(o[i]);
        }
    }

    template <typename T2>
    void flip(const img_t<T2>&o) {
        assert(o.similar(*this));

#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) collapse(3) num_threads(available_threads) if (available_threads > 1)
#endif
       for (int y = 0; y < o.h; y++) {
            for (int x = 0; x < o.w; x++) {
                for (int dd = 0; dd < d; dd++) {
                    (*this)(x, y, dd) = o(w-x-1, h-y-1, dd);
                }
            }
        }
    }

    void copy(const img_t<T>& o) {
        assert(o.similar(*this));
        std::copy(o.data.begin(), o.data.end(), data.begin());
    }

    template <typename T2>
    void copy(const img_t<T2>& o) {
        assert(o.similar(*this));
        std::copy(o.data.begin(), o.data.end(), data.begin());
    }

    bool inside(int x, int y, int dd=0) const {
        return x >= 0 && x < w && y >= 0 && y < h && dd >= 0 && dd < d;
    }

    // forward differences
    template <typename T2>
    void gradients(const img_t<T2>& u_) {
        const img_t<typename T::value_type>& u = *static_cast<const img_t<typename T::value_type>*>(&u_);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) collapse(3) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int l = 0; l < d; l++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    (*this)(x, y, l)[0] = (x < w - 1) ? (u(x+1, y, l) - u(x, y, l)) : 0;
                    (*this)(x, y, l)[1] = (y < h - 1) ? (u(x, y+1, l) - u(x, y, l)) : 0;
                }
            }
        }
    }

    template <typename T2>
    void circular_gradients(const img_t<T2>& u_) {
        const img_t<typename T::value_type>& u = *static_cast<const img_t<typename T::value_type>*>(&u_);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
            #pragma omp for collapse(3) schedule(guided)
#endif
            for (int l = 0; l < d; l++) {
                for (int y = 0; y < h; y++) {
                    for (int x = 0; x < w; x++) {
                        (*this)(x, y, l)[0] = u((x+1)%w, y, l) - u(x, y, l);
                        (*this)(x, y, l)[1] = u(x, (y+1)%h, l) - u(x, y, l);
                    }
                }
            }
#ifdef _OPENMP
        }
#endif
    }

    void gradientx(const img_t<T>& u) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
            #pragma omp for collapse(2) schedule(guided)
#endif
            for (int l = 0; l < d; l++) {
                for (int y = 0; y < h; y++) {
#ifdef _OPENMP
                    #pragma omp simd
#endif
                    for (int x = 0; x < w-1; x++) {
                        (*this)(x, y, l) = u(x+1, y, l) - u(x, y, l);
                    }
                    (*this)(w-1, y, l) = 0;
                }
            }
#ifdef _OPENMP
        }
#endif
    }

    void gradienty(const img_t<T>& u) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) collapse(3) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int l = 0; l < d; l++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    (*this)(x, y, l) = (y < h - 1) ? (u(x, y+1, l) - u(x, y, l)) : 0;
                }
            }
        }
    }

    void gradientxx(const img_t<T>& u) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
            #pragma omp for collapse(2) schedule(guided)
#endif
            for (int l = 0; l < d; l++) {
                for (int y = 0; y < h; y++) {
#ifdef _OPENMP
                    #pragma omp simd
#endif
                    for (int x = 1; x < w-1; x++) {
                        (*this)(x, y, l) = u(x+1, y, l) + u(x-1, y, l) - 2 * u(x, y, l);
                    }
                    (*this)(0, y, l) = 0;
                    (*this)(w-1, y, l) = 0;
                }
            }
#ifdef _OPENMP
        }
#endif
    }

    void gradientyy(const img_t<T>& u) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
            #pragma omp for collapse(2) schedule(guided)
#endif
            for (int l = 0; l < d; l++) {
                for (int x = 0; x < w; x++) {
#ifdef _OPENMP
                    #pragma omp simd
#endif
                    for (int y = 1; y < h-1; y++) {
                        (*this)(x, y, l) = u(x, y+1, l) + u(x, y-1, l) - 2 * u(x, y, l);
                    }
                    (*this)(x, 0, l) = 0;
                    (*this)(x, h-1, l) = 0;
                }
            }
#ifdef _OPENMP
        }
#endif
    }

    void gradientxy(const img_t<T>& u) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) collapse(3) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int l = 0; l < d; l++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    (*this)(x, y, l) = (x < w - 1 && y < h - 1) ?
                        (u(x+1, y+1, l) - u(x+1, y, l) - u(x, y+1, l) + u(x, y, l)) : 0;
                }
            }
        }
    }

    template <typename T2>
    void divergence(const img_t<T2>& g) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
#endif
            for (int l = 0; l < d; l++) {
#ifdef _OPENMP
                #pragma omp for schedule(guided)
#endif
                for (int y = 1; y < h-1; y++) {
#ifdef _OPENMP
                    #pragma omp simd
#endif
                    for (int x = 1; x < w-1; x++) {
                        (*this)(x, y, l) = g(x, y, l)[0] - g(x-1, y, l)[0] + g(x, y, l)[1] - g(x, y-1, l)[1];
                    }
                }

#ifdef _OPENMP
                #pragma omp single
                {
#endif
                    (*this)(0, 0, l) = g(0, 0, l)[0] + g(0, 0, l)[1];
                    (*this)(w-1, 0, l) = -g(w-2, 0, l)[0] + g(w-1, 0, l)[1];
                    (*this)(0, h-1, l) = g(0, h-1, l)[0] - g(0, h-2, l)[1];
                    (*this)(w-1, h-1, l) = -g(w-2, h-1, l)[0] - g(w-1, h-2, l)[1];
#ifdef _OPENMP
                }

                #pragma omp sections
                {
                    #pragma omp section
                    {
                        #pragma omp simd
#endif
                        for (int y = 1; y < h-1; y++)
                            (*this)(0, y, l) = g(0, y, l)[0] + g(0, y, l)[1] - g(0, y-1, l)[1];
#ifdef _OPENMP
                    }
                    #pragma omp section
                    {
                        #pragma omp simd
#endif
                        for (int x = 1; x < w-1; x++)
                            (*this)(x, 0, l) = g(x, 0, l)[0] - g(x-1, 0, l)[0] + g(x, 0, l)[1];
#ifdef _OPENMP
                    }
                    #pragma omp section
                    {
                        #pragma omp simd
#endif
                        for (int y = 1; y < h-1; y++)
                            (*this)(w-1, y, l) = -g(w-2, y, l)[0] + g(w-1, y, l)[1] - g(w-1, y-1, l)[1];
#ifdef _OPENMP
                    }
                    #pragma omp section
                    {
                        #pragma omp simd
#endif
                        for (int x = 1; x < w-1; x++)
                            (*this)(x, h-1, l) = g(x, h-1, l)[0] - g(x-1, h-1, l)[0] - g(x, h-2, l)[1];
#ifdef _OPENMP
                    }
                }
#endif
            }
#ifdef _OPENMP
        }
#endif
    }

    template <typename T2>
    void circular_divergence(const img_t<T2>& g) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) collapse(3) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int l = 0; l < d; l++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    (*this)(x, y, l) = g(x, y, l)[0] - g((x-1+w)%w, y, l)[0] + g(x, y, l)[1] - g(x, (y-1+h)%h, l)[1];
                }
            }
        }
    }

    void divergence(const img_t<T>& gx, const img_t<T>& gy) {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
#endif
            for (int l = 0; l < d; l++) {
#ifdef _OPENMP
                #pragma omp for schedule(guided)
#endif
                for (int y = 1; y < h-1; y++) {
#ifdef _OPENMP
                    #pragma omp simd
#endif
                    for (int x = 1; x < w-1; x++) {
                        (*this)(x, y, l) = gx(x, y, l) - gx(x-1, y, l) + gy(x, y, l) - gy(x, y-1, l);
                    }
                }

#ifdef _OPENMP
                #pragma omp single
                {
#endif
                    (*this)(0, 0, l) = gx(0, 0, l) + gy(0, 0, l);
                    (*this)(w-1, 0, l) = -gx(w-2, 0, l) + gy(w-1, 0, l);
                    (*this)(0, h-1, l) = gx(0, h-1, l) - gy(0, h-2, l);
                    (*this)(w-1, h-1, l) = -gx(w-2, h-1, l) - gy(w-1, h-2, l);
#ifdef _OPENMP
                }

                #pragma omp sections
                {
                    #pragma omp section
                    {
                        #pragma omp simd
#endif
                        for (int y = 1; y < h-1; y++)
                            (*this)(0, y, l) = gx(0, y, l) + gy(0, y, l) - gy(0, y-1, l);
#ifdef _OPENMP
                    }
                    #pragma omp section
                    {
                        #pragma omp simd
#endif
                        for (int x = 1; x < w-1; x++)
                            (*this)(x, 0, l) = gx(x, 0, l) - gx(x-1, 0, l) + gy(x, 0, l);
#ifdef _OPENMP
                    }
                    #pragma omp section
                    {
                        #pragma omp simd
#endif
                        for (int y = 1; y < h-1; y++)
                            (*this)(w-1, y, l) = -gx(w-2, y, l) + gy(w-1, y, l) - gy(w-1, y-1, l);
#ifdef _OPENMP
                    }
                    #pragma omp section
                    {
                        #pragma omp simd
#endif
                        for (int x = 1; x < w-1; x++)
                            (*this)(x, h-1, l) = gx(x, h-1, l) - gx(x-1, h-1, l) - gy(x, h-2, l);
#ifdef _OPENMP
                    }
                }
#endif
            }
#ifdef _OPENMP
        }
#endif
    }
    void fft(const img_t<std::complex<float> >& o) {
        static_assert(std::is_same<T, std::complex<float>>::value, "T must be complex float");
        assert(w == o.w);
        assert(h == o.h);
        assert(d == o.d);
        auto out = reinterpret_cast<fftwf_complex*>(&data[0]);
        if (!forwardplanf) {
            img_t<T> tmp(w, h, d);
            tmp.copy(*this);
            int n[] = {h, w};
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
            if (!(forwardplanf = fftwf_plan_many_dft(2, n, d, out, n, d, 1, out, n, d, 1, FFTW_FORWARD, FFTW_WISDOM_ONLY)))
                forwardplanf = fftwf_plan_many_dft(2, n, d, out, n, d, 1, out, n, d, 1, FFTW_FORWARD, FFTW_MEASURE);

            copy(tmp);
        }
        if (this != &o) {
            this->copy(o);
        }
        fftwf_execute(forwardplanf);
    }

    void ifft(const img_t<std::complex<float> >& o) {
        static_assert(std::is_same<T, std::complex<float>>::value, "T must be complex float");
        assert(w == o.w);
        assert(h == o.h);
        assert(d == o.d);
        auto out = reinterpret_cast<fftwf_complex*>(&data[0]);
        if (!backwardplanf) {
            img_t<T> tmp(w, h, d);
            tmp.copy(*this);
            int n[] = {h, w};
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
            if (!(backwardplanf = fftwf_plan_many_dft(2, n, d, out, n, d, 1, out, n, d, 1, FFTW_BACKWARD, FFTW_WISDOM_ONLY)))
                backwardplanf = fftwf_plan_many_dft(2, n, d, out, n, d, 1, out, n, d, 1, FFTW_BACKWARD, FFTW_MEASURE);
            copy(tmp);
        }
        float norm = w * h;
        if (this == &o) {
            for (auto& x : data) {
                x /= norm;
            }
        } else {
            this->map(o, [norm](T x){ return x / norm; });
        }
        fftwf_execute(backwardplanf);
    }

    void fftshift() {
        img_t<T> out;
        out.resize(*this);

        int halfw = (this->w + 1) / 2;
        int halfh = (this->h + 1) / 2;
        int ohalfw = this->w - halfw;
        int ohalfh = this->h - halfh;

#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
#endif
            for (int l = 0; l < this->d; l++) {
#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < halfh; y++) {
                    for (int x = 0; x < ohalfw; x++) {
                        out(x, y + ohalfh, l) = (*this)(x + halfw, y, l);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < halfh; y++) {
                    for (int x = 0; x < halfw; x++) {
                        out(x + ohalfw, y + ohalfh, l) = (*this)(x, y, l);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < ohalfh; y++) {
                    for (int x = 0; x < ohalfw; x++) {
                        out(x, y, l) = (*this)(x + halfw, y + halfh, l);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < ohalfh; y++) {
                    for (int x = 0; x < halfw; x++) {
                        out(x + ohalfw, y, l) = (*this)(x, y + halfh, l);
                    }
                }
            }
#ifdef _OPENMP
        }
#endif

        *this = std::move(out);
    }

    void ifftshift() {
        img_t<T> out;
        out.resize(*this);

        int halfw = (this->w + 1) / 2;
        int halfh = (this->h + 1) / 2;
        int ohalfw = this->w - halfw;
        int ohalfh = this->h - halfh;

#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
#endif
            for (int l = 0; l < this->d; l++) {
#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < ohalfh; y++) {
                    for (int x = 0; x < halfw; x++) {
                        out(x, y + halfh, l) = (*this)(x + ohalfw, y, l);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < ohalfh; y++) {
                    for (int x = 0; x < ohalfw; x++) {
                        out(x + halfw, y + halfh, l) = (*this)(x, y, l);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < halfh; y++) {
                    for (int x = 0; x < halfw; x++) {
                        out(x, y, l) = (*this)(x + ohalfw, y + ohalfh, l);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < halfh; y++) {
                    for (int x = 0; x < ohalfw; x++) {
                        out(x + halfw, y, l) = (*this)(x, y + ohalfh, l);
                    }
                }
            }
#ifdef _OPENMP
        }
#endif

        *this = std::move(out);
    }

    template <typename T2>
    void padcirc(const img_t<T2>& o) {
        set_value(0);
        int ww = o.w / 2;
        int hh = o.h / 2;

#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel num_threads(available_threads) if (available_threads > 1)
        {
#endif
            for (int dd = 0; dd < d; dd++) {
                int od = 0;
                if (d == o.d)
                    od = dd;
                else if (o.d == 1)
                    od = 0;
                else {
                    printf("Error: depths are %d and %d\n", d, o.d);
                    assert(false);
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < hh; y++) {
                    for (int x = 0; x < ww; x++) {
                        (*this)(w - ww + x, h - hh + y, dd) = o(x, y, od);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = 0; y < hh; y++) {
                    for (int x = ww; x < o.w; x++) {
                        (*this)(- ww + x, h - hh + y, dd) = o(x, y, od);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = hh; y < o.h; y++) {
                    for (int x = 0; x < ww; x++) {
                        (*this)(w - ww + x, - hh + y, dd) = o(x, y, od);
                    }
                }

#ifdef _OPENMP
                #pragma omp for collapse(2) schedule(static)
#endif
                for (int y = hh; y < o.h; y++) {
                    for (int x = ww; x < o.w; x++) {
                        (*this)(- ww + x, - hh + y, dd) = o(x, y, od);
                    }
                }
            }
#ifdef _OPENMP
        }
#endif
    }
    void ensure_size(int w, int h, int d=1) {
        assert(w > 0);
        assert(h > 0);
        assert(d > 0);
        if (this->w != w || this->h != h || this->d != d) {
            this->w = w;
            this->h = h;
            this->d = d;
            size = w * h * d;
            data.resize(size);

            if (forwardplanf) {
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
                fftwf_destroy_plan(forwardplanf);
                forwardplanf = nullptr;
            }
            if (backwardplanf) {
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
                fftwf_destroy_plan(backwardplanf);
                backwardplanf = nullptr;
            }
        }
    }

    T mean() const {
        T foo = (*this).sum();
        foo /= this->size;
        return foo;
    }

    void normalize() {
        T sum = this->sum();
        if (sum != 0.) {
            for (T& v : data) {
                v /= sum;
            }
        }
    }

    void sanitize() {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < size; i++) {
            if (((*this)[i] != (*this)[i]) || (*this)[i] == T(0))
                (*this)[i] = T(1.e-9);
        }
    }

    void greyfromcolor(const img_t<T>& color) {
        assert(d == 1);
        assert(w == color.w);
        assert(h == color.h);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < size; i++) {
            T val(0);
            for (int dd = 0; dd < color.d; dd++) {
                val += color[i * color.d + dd];
            }
            (*this)[i] = val / color.d;
        }
    }

    void desaturate() {
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(2) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                T val = T(0);
                for (int c = 0; c < d; c++) {
                    val += (*this)(i, j, c);
                }
                (*this)(i, j, 0) = val / d;
            }
        }
        d = 1;
        data.resize(w * h);
    }

    void transpose() {
        img_t<T> o(*this);
        this->w = o.h;
        this->h = o.w;
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(3) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int y = 0; y < o.h; y++) {
            for (int x = 0; x < o.w; x++) {
                for (int dd = 0; dd < d; dd++) {
                    (*this)(y, x, dd) = o(x, y, dd);
                }
            }
        }
    }

    void transposeToMatlab() {
        img_t<T> o(*this);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(3) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int dd = 0; dd < d; dd++) {
                    (*this)[y + h * (x + w * dd)] = o(x, y, dd);
                }
            }
        }
    }

    void transposeFromMatlab() {
        img_t<T> o(*this);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for collapse(3) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int dd = 0; dd < d; dd++) {
                    (*this)(x, y, dd) = o[y + h * (x + w * dd)];
                }
            }
        }
    }
};

namespace img {
    template <typename T, typename E>
    T sum(const E& img) {
        T a(0);
#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for reduction(+:a) schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < img.size; i++) {
            a += img[i];
        }
        return a;
    }

    template <typename T, typename E>
    T sumL1(const E& img) {
        return sum<T>(std::abs(img));
    }

    // Frobenius norm
    template <typename T, typename E>
    T sumL2(const E& img) {
        return std::sqrt(sum<T>(img*img));
    }

    // convert an image to YCbCr colorspace (from RGB)
    template <typename T>
    void rgb2ycbcr(img_t<T>& out, const img_t<T>& in)
    {
        assert(in.d == 3);
        out.ensure_size(in.w, in.h, in.d);

#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < out.w*out.h; i++) {
            T r = in[i*3+0];
            T g = in[i*3+1];
            T b = in[i*3+2];
            out[i*3+0] = 0.299*r + 0.587*g + 0.114*b;
            out[i*3+1] = (b - out[i*3+0]) * 0.564 + 0.5;
            out[i*3+2] = (r - out[i*3+0]) * 0.713 + 0.5;
        }
    }

    /// convert an image to RGB colorspace (from YCbCr)
    template <typename T>
    void ycbcr2rgb(img_t<T>& out, const img_t<T>& in)
    {
        assert(in.d == 3);
        out.ensure_size(in.w, in.h, in.d);

#ifdef _OPENMP
        int available_threads = com.max_thread - omp_get_num_threads();
        #pragma omp parallel for schedule(static) num_threads(available_threads) if (available_threads > 1)
#endif
        for (int i = 0; i < out.w*out.h; i++) {
            T y = in[i*3+0];
            T cb = in[i*3+1];
            T cr = in[i*3+2];
            out[i*3+0] = y + 1.403 * (cr - 0.5);
            out[i*3+1] = y - 0.714 * (cr - 0.5) - 0.344 * (cb - 0.5);
            out[i*3+2] = y + 1.773 * (cb - 0.5);
        }
    }
};

