#pragma once

#include "image.hpp"
#include "image_expr.hpp"
#include "vec2.hpp"

namespace optimization {

    template <typename T>
    struct Functional {
        typedef T value_type;
        auto gradient(const img_t<T>& x);
        auto prox(const img_t<T>& x, T tau);
        auto prox_adjoint(const img_t<T>& x, T tau);
    };

    template <typename T, typename T2>
    struct Operator {
        typedef T in_type;
        typedef T2 out_type;
        auto direct(const img_t<T>& x);
        auto adjoint(const img_t<T2>& x);
    };

    namespace operators {
        template <typename T>
        struct identity : Operator<T, T> {
            auto direct(const img_t<T>& x) {
                return to_expr(x);
            }
            auto adjoint(const img_t<T>& x) {
                return to_expr(x);
            }
        };

        template <typename T>
        struct gradient : Operator<T, vec2<T>> {
            img_t<vec2<T>> grad;
            img_t<T> div;
            gradient(const img_t<T> f) : grad(f.w, f.h, f.d), div(f.w, f.h, f.d) {}

            auto direct(const img_t<T>& x) {
                grad.gradients(x);
                return to_expr(grad);
            }

            auto adjoint(const img_t<vec2<T>>& x) {
                div.divergence(x);
                return -div;
            }
        };

        template <typename T>
        struct circular_gradient : Operator<T, vec2<T>> {
            img_t<vec2<T>> grad;
            img_t<T> div;
            circular_gradient(const img_t<T> f) : grad(f.w, f.h, f.d), div(f.w, f.h, f.d) {}

            auto direct(const img_t<T>& x) {
                grad.circular_gradients(x);
                return to_expr(grad);
            }

            auto adjoint(const img_t<vec2<T>>& x) {
                div.circular_divergence(x);
                return -div;
            }
        };

        template <typename T>
        struct circular_convolution : Operator<T, T> {
            img_t<std::complex<T>> ftK;
            img_t<std::complex<T>> ftx;

            circular_convolution(const img_t<T>& f, const img_t<T>& K) :
                    ftx(f.w, f.h, f.d),
                    ftK(f.w, f.h, f.d) {
                assert(K.w % 2);
                assert(K.h % 2);

                ftK.padcirc(K);
                ftK.map(ftK * T(K.d) / K.sum());
                ftK.fft(ftK);
            }

            auto direct(const img_t<T>& x) {
                ftx.map(x);
                ftx.fft(ftx);
                ftx.map(ftx * ftK);
                ftx.ifft(ftx);
                return std::real(ftx);
            }

            auto adjoint(const img_t<T>& x) {
                ftx.map(x);
                ftx.fft(ftx);
                ftx.map(ftx * std::conj(ftK));
                ftx.ifft(ftx);
                return std::real(ftx);
            }
        };
    }
}
