#include <array>

#include "image.hpp"
#include "deblur.hpp"

/// pad an image using constant boundaries
template <typename T>
static void padimage_replicate(img_t<T>& out, const img_t<T>& in, int padding)
{
    out.ensure_size(in.w + padding*2, in.h + padding*2, in.d);

    for (int y = 0; y < in.h; y++) {
        for (int x = 0; x < in.w; x++) {
            for (int l = 0; l < in.d; l++) {
                out(x+padding, y+padding, l) = in(x, y, l);
            }
        }
    }

    // pad top and bottom
    for (int x = 0; x < out.w; x++) {
        int xx = std::min(std::max(0, x - padding), in.w-1);
        for (int l = 0; l < in.d; l++) {
            T val_top = in(xx, 0, l);
            T val_bottom = in(xx, in.h-1, l);
            for (int y = 0; y < padding; y++) {
                out(x, y, l) = val_top;
                out(x, out.h-1 - y, l) = val_bottom;
            }
        }
    }

    // pad left and right
    for (int y = 0; y < out.h; y++) {
        int yy = std::min(std::max(0, y - padding), in.h-1);
        for (int l = 0; l < in.d; l++) {
            T val_left = in(0, yy, l);
            T val_right = in(in.w-1, yy, l);
            for (int x = 0; x < padding; x++) {
                out(x, y, l) = val_left;
                out(out.w-1 - x, y, l) = val_right;
            }
        }
    }
}

/// remove the padding of an image
template <typename T>
static void unpadimage(img_t<T>& out, const img_t<T>& in, int padding)
{
    out.ensure_size(in.w - 2*padding, in.h - 2*padding, in.d);

    for (int y = 0; y < out.h; y++) {
        for (int x = 0; x < out.w; x++) {
            for (int l = 0; l < out.d; l++) {
                out(x, y, l) = in(x+padding, y+padding, l);
            }
        }
    }
}

template <typename T>
void pad_and_taper(img_t<T>& u, const img_t<T>& f, const img_t<T>& K)
{
    int padding = std::max(K.w, K.h);
    img_t<T> padded;
    padimage_replicate(padded, f, padding);

    img_t<T> tapered;
    edgetaper(u, padded, K, 4);
}

template <typename T>
void unpad(img_t<T>& u, const img_t<T>& f, const img_t<T>& K)
{
    int padding = std::max(K.w, K.h);
    unpadimage(u, f, padding);
}

// convert an image to YCbCr colorspace (from RGB)
template <typename T>
static void rgb2ycbcr(img_t<T>& out, const img_t<T>& in)
{
    assert(in.d == 3);

    out.ensure_size(in.w, in.h, in.d);

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
static void ycbcr2rgb(img_t<T>& out, const img_t<T>& in)
{
    assert(in.d == 3);

    out.ensure_size(in.w, in.h, in.d);

    for (int i = 0; i < out.w*out.h; i++) {
        T y = in[i*3+0];
        T cb = in[i*3+1];
        T cr = in[i*3+2];
        out[i*3+0] = y + 1.403 * (cr - 0.5);
        out[i*3+1] = y - 0.714 * (cr - 0.5) - 0.344 * (cb - 0.5);
        out[i*3+2] = y + 1.773 * (cb - 0.5);
    }
}


/// deconvolve an image using Split bregman
/// deconvolve only the luminance
/// boundaries have to be handled elsewhere
template <typename T>
void deconvBregman(img_t<T>& u, const img_t<T>& f, const img_t<T>& K,
                  int numIter, T lambda, T beta)
{
    if (f.d == 3) {
        // convert to YCbCr
        img_t<T> ycbcr;
        rgb2ycbcr(ycbcr, f);
        img_t<T> y(ycbcr.w, ycbcr.h);
        for (int i = 0; i < y.w*y.h; i++)
            y[i] = ycbcr[i*3];

        // deconvolve Y
        img_t<T> ydeconv;
        deconvBregman(ydeconv, y, K, numIter, lambda, beta);

        // convert to RGB
        for (int i = 0; i < y.w*y.h; i++)
            ycbcr[i*3] = ydeconv[i];
        ycbcr2rgb(u, ycbcr);
        return;
    }

    // reorder to planar
    img_t<T> f_planar(f.w, f.h, f.d);
    img_t<T> deconv_planar(f.w, f.h, f.d);
    if (f.d != 1) {
        for (int y = 0; y < f.h; y++) {
            for (int x = 0; x < f.w; x++) {
                for (int l = 0; l < f.d; l++) {
                    f_planar[x + f.w*(y + f.h*l)] = f(x, y, l);
                    deconv_planar[x + f.w*(y + f.h*l)] = f(x, y, l);
                }
            }
        }
    } else {
        f_planar.copy(f);
        deconv_planar.copy(f);
    }

    // deconvolve
    deblur::rof::split_continuation(deconv_planar, f_planar, K, 2.f / lambda, beta, 2.f * std::sqrt(2.f), std::pow(2.f, 8.f), 1);

    // reorder to interleaved
    u.ensure_size(deconv_planar.w, deconv_planar.h, deconv_planar.d);
    if (u.d != 1) {
        for (int y = 0; y < u.h; y++) {
            for (int x = 0; x < u.w; x++) {
                for (int l = 0; l < u.d; l++) {
                    u(x, y, l) = deconv_planar[x + u.w*(y + u.h*l)];
                }
            }
        }
    } else {
        u.copy(deconv_planar);
    }
}

