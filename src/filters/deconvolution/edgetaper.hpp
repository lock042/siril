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

// Moved into a separate file out of deconvBregman.cpp

#pragma once

#include "image.hpp"
#include <glib.h>

template <typename T>
void edgetaper(img_t<T>& out, const img_t<T>& in, const img_t<T>& kernel, int iterations=1)
{
    img_t<T> weights(in.w, in.h);
    // kind of tukey window
    for (int y = 0; y < in.h; y++) {
        T wy = 1.;
        if (y < kernel.h) {
            wy = std::pow(std::sin(y * G_PI / (kernel.h*2 - 1)), 2.);
        } else if (y > in.h - kernel.h) {
            wy = std::pow(std::sin((in.h-1 - y) * G_PI / (kernel.h*2 - 1)), 2.);
        }
        for (int x = 0; x < in.w; x++) {
            T wx = 1.;
            if (x < kernel.w) {
                wx = std::pow(std::sin(x * G_PI / (kernel.w*2 - 1)), 2.);
            } else if (x > in.w - kernel.w) {
                wx = std::pow(std::sin((in.w-1 - x) * G_PI / (kernel.w*2 - 1)), 2.);
            }
            weights(x, y) = wx * wy;
        }
    }

    // kernel's fft
    img_t<T> blurred(in.w, in.h, in.d);
    img_t<std::complex<T>> kernel_ft(in.w, in.h, in.d);
    kernel_ft.padcirc(kernel);
    kernel_ft.fft(kernel_ft);

    img_t<std::complex<T>> blurred_ft(in.w, in.h, in.d);

    out = in;
    for (int i = 0; i < iterations; i++) {
        blurred_ft.copy(out);

        blurred_ft.fft(blurred_ft);
        for (int y = 0; y < out.h; y++)
            for (int x = 0; x < out.w; x++)
                for (int l = 0; l < out.d; l++)
                    blurred_ft(x, y, l) *= kernel_ft(x, y);
        blurred_ft.ifft(blurred_ft);

        for (int i = 0; i < blurred.size; i++)
            blurred[i] = std::real(blurred_ft[i]);

        // blend the images
        for (int y = 0; y < out.h; y++) {
            for (int x = 0; x < out.w; x++) {
                T w = weights(x, y);
                for (int l = 0; l < out.d; l++) {
                    out(x, y, l) = w * out(x, y, l) + (1. - w) * blurred(x, y, l);
                }
            }
        }
    }
}


