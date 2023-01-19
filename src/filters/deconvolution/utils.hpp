/*

From https://ipolcore.ipol.im/demo/clientApp/demo.html?id=77777000086
Anger, Jérémy, Mauricio Delbracio, and Gabriele Facciolo. "Efficient Blind Deblurring under High Noise Levels.", Published at International Symposium on Image and Signal Processing and Analysis (ISPA 2019).

Licenced under the GNU AFFERO GENERAL PUBLIC LICENSE, Version 3, 19 November 2007

*/
#pragma once

#include "image.hpp"
#include "image_expr.hpp"
#include "labeling.hpp"
#include "vec2.hpp"
#include "chelperfuncs.h"

namespace utils {

    inline void blur(img_t<float>& out, const img_t<float>& in, float sigma) {
        assert(in.d == 1);
        out.resize(in.w, in.h);
        gaussblur(&out[0], (float*) &in[0], in.w, in.h, sigma);
    }

    // upsample an image
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
        for (int d = 0; d < in.d; d++) {
            for (int y = 0; y < in.h; y++) {
                for (int x = 0; x < in.w; x++) {
                    out(y, x, d) = in(x, y, d);
                }
            }
        }
    }

    // add zero padding of the size of the kernel
    template <typename T>
    img_t<T> zero_pad(const img_t<T>& _f, int hw, int hh)
    {
        img_t<T> f(_f.w + hw*2, _f.h + hh*2, _f.d);
        f.set_value(T(0));
        slice(f, _(hw, -hw-1), _(hh, -hh-1)).map(_f);
        return f;
    }

    // add symmetric padding of the size of the kernel
    template <typename T>
    img_t<T> add_padding(const img_t<T>& _f, int hw, int hh)
    {
        img_t<T> f(_f.w + hw*2, _f.h + hh*2, _f.d);
        f.set_value(T(0));
        slice(f, _(hw, -hw-1), _(hh, -hh-1)).map(_f);
        // replicate borders
        for (int y = 0; y < hh; y++) {
            for (int x = 0; x < f.w; x++) {
                for (int l = 0; l < f.d; l++) {
                    f(x, y, l) = f(x, 2*hh - y, l);
                    f(x, f.h-1-y, l) = f(x, f.h-1-2*hh+y, l);
                }
            }
        }
        for (int y = 0; y < f.h; y++) {
            for (int x = 0; x < hw; x++) {
                for (int l = 0; l < f.d; l++) {
                    f(x, y, l) = f(2*hw - x, y, l);
                    f(f.w-1-x, y, l) = f(f.w-1-2*hw+x, y, l);
                }
            }
        }
        return f;
    }

    template <typename T>
    img_t<T> add_padding(const img_t<T>& f, const img_t<T>& K)
    {
//        printf("Kernel width %d and height %d\n", K.w, K.h);
        return add_padding(f, K.w/2, K.h/2);
    }

    template <typename T>
    img_t<T> remove_padding(const img_t<T>& f, int hw, int hh)
    {
        return to_img(slice(f, _(hw, -hw-1), _(hh, -hh-1)));
    }

    template <typename T>
    img_t<T> remove_padding(const img_t<T>& f, const img_t<T>& K)
    {
        return remove_padding(f, K.w/2, K.h/2);
    }

    template <typename T>
    img_t<T> crop_to_evens(const img_t<T>& f)
    {
        int xcrop = 0, ycrop = 0;
        if (f.w % 2)
            xcrop = 1;
        if (f.h % 2)
            ycrop = 1;
        return to_img(slice(f, _(0, -xcrop-1), _(0, -ycrop-1)));
    }

    template <typename T>
    void center_kernel(img_t<T>& kernel) {
        T dx = 0.f;
        T dy = 0.f;
        T sum = kernel.sum();
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
                int nx = x + (int)dx - kernel.w/2;
                int ny = y + (int)dy - kernel.h/2;
                if (nx >= 0 && nx < kernel.w && ny >= 0 && ny < kernel.h) {
                    kernel(x, y) = copy(nx, ny);
                }
            }
        }
    }

    // compute connected component of the support of k
    // and set to zero pixels belonging to low intensity connected components
    template <typename T>
    void remove_isolated_cc(img_t<T>& k) {
        T sum = k.sum();
        for (int i = 0; i < k.size; i++)
            k[i] /= sum;
        img_t<int> lab;
        labeling::labels(lab, k);
        auto sums = labeling::sum(lab, k);
        for (int i = 0; i < k.size; i++) {
            if (sums[lab[i]] < T(0.1))
                k[i] = T(0);
        }
    }

    // compute the circular gradients by forward difference
    template <typename T>
    void circular_gradients(vec2<img_t<T>>& out, const img_t<T>& in) {
        out[0].resize(in);
        out[1].resize(in);

        int w = in.w;
        int h = in.h;
        int d = in.d;
        for (int l = 0; l < d; l++) {
            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                    out[0](x, y, l) = in((x+1)%w, y, l) - in(x, y, l);

            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                    out[1](x, y, l) = in(x, (y+1)%h, l) - in(x, y, l);
        }
    }

    // compute the circular divergence by backward difference
    template <typename T>
    void circular_divergence(img_t<T>& out, const vec2<img_t<T>>& in) {
        out.resize(in[0]);

        int w = out.w;
        int h = out.h;
        int d = out.d;
        for (int l = 0; l < d; l++) {
            for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
                out(x, y, l) = in[0](x, y, l) - in[0]((x-1+w)%w, y, l)
                             + in[1](x, y, l) - in[1](x, (y-1+h)%h, l);
        }
    }

    template <typename T>
    T getpixel_1(const img_t<T>& x, int i, int j, int d=0)
    {
        i = std::max(std::min(i, x.w - 1), 0);
        j = std::max(std::min(j, x.h - 1), 0);
        return x(i, j, d);
    }

    template <typename T>
    void downsa2(img_t<T>& out, const img_t<T>& in)
    {
        if (out.size == 0)
            out.resize(in.w/2, in.h/2, in.d);
        for (int d = 0; d < out.d; d++)
        for (int j = 0; j < out.h; j++)
        for (int i = 0; i < out.w; i++)
        {
            T m = getpixel_1(in, 2*i, 2*j, d)
                + getpixel_1(in, 2*i+1, 2*j, d)
                + getpixel_1(in, 2*i, 2*j+1, d)
                + getpixel_1(in, 2*i+1, 2*j+1, d);
            out(i, j, d) = m / T(4);
        }
    }

    template <typename T>
    T evaluate_bilinear_cell(T a[4], float x, float y)
    {
        T r = 0;
        r += a[0] * (1-x) * (1-y);
        r += a[1] * ( x ) * (1-y);
        r += a[2] * (1-x) * ( y );
        r += a[3] * ( x ) * ( y );
        return r;
    }

    template <typename T>
    T bilinear_interpolation(const img_t<T>& x, float p, float q, int d)
    {
        int ip = floor(p);
        int iq = floor(q);
        T a[4] = {
            getpixel_1(x, ip  , iq  , d),
            getpixel_1(x, ip+1, iq  , d),
            getpixel_1(x, ip  , iq+1, d),
            getpixel_1(x, ip+1, iq+1, d)
        };
        T r = evaluate_bilinear_cell(a, p-ip, q-iq);
        return r;
    }

    template <typename T>
    void upsa2(img_t<T>& out, const img_t<T>& in)
    {
        if (out.size == 0)
            out.resize(in.w*2, in.h*2, in.d);
        for (int d = 0; d < out.d; d++)
        for (int j = 0; j < out.h; j++)
        for (int i = 0; i < out.w; i++)
        {
            float x = (i - 0.5) / 2;
            float y = (j - 0.5) / 2;
            out(i, j, d) = bilinear_interpolation(in, x, y, d);
        }
    }

}

