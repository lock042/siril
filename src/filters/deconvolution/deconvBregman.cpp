#include <array>

#include "image.hpp"
#include "deconvolve.hpp"

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
        img::rgb2ycbcr(ycbcr, f);
        img_t<T> y(ycbcr.w, ycbcr.h);
        for (int i = 0; i < y.w*y.h; i++)
            y[i] = ycbcr[i*3];

        // deconvolve Y
        img_t<T> ydeconv;
        deconvBregman(ydeconv, y, K, numIter, lambda, beta);

        // convert to RGB
        for (int i = 0; i < y.w*y.h; i++)
            ycbcr[i*3] = ydeconv[i];
        img::ycbcr2rgb(u, ycbcr);
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
    deconvolve::sb_deconvolve(deconv_planar, f_planar, K, 2.f / lambda, beta, 2.f * std::sqrt(2.f), std::pow(2.f, 8.f), 1);

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

