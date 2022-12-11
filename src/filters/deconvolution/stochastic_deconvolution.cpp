#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <set>
#include <map>

#include <omp.h>

#include "image.hpp"
#include "image_expr.hpp"
#include "utils.hpp"
#include "drand48.h"

int NORMALIZE_INPUT = 0;

float reset_probability = 0.005f;
float mutation_sigma = 4.f;
//float REG_WEIGHT = 0.0008f; // Uses lambda instead
int ITERATIONS = 300;
float ed_init = 0.025f;
float clip = 0.95f;
float quantize = 0.f; //1.f/65536.f;
float gamma_value = 0.f;
float gamma_tv = 0.f;
int fast = 1;
int multiscale = 1;

static std::default_random_engine generator;

void mutate(int& x, int& y, int& l, int w, int h, int d)
{
    if (drand48() > reset_probability) {
        std::normal_distribution<float> distribution(0., mutation_sigma);
        int xx = x;
        int yy = y;
        do {
            xx = std::round(x + distribution(generator));
            yy = std::round(y + distribution(generator));
        } while ((xx == x && yy == y)
                 || xx < 0 || xx >= w
                 || yy < 0 || yy >= h);
        x = xx;
        y = yy;
        // same channel
    } else {
        x = drand48() * w;
        y = drand48() * h;
        l = rand() % d;
    }
}

template <typename T, typename P>
T evaluate(img_t<T>& u, P& problem, int x, int y, int l, T& ed)
{
    T init = problem.energy(u, x, y, l, 0);
    T plus = problem.energy(u, x, y, l, ed);
    T minus = problem.energy(u, x, y, l, -ed);

    if (minus < plus) {
        ed = - ed;
        return init - minus;
    }
    return init - plus;
}

template <typename T, typename P>
void stochastic_deconvolution(img_t<T>& u, P& problem, int w, int h, int d)
{
    bool positivity = gamma_value > 0;
    T ED(ed_init);
    const int mutations = w * h * d;
    for (int i = 0; i < ITERATIONS; i++) {
        int rate = 0;

        int sections = 1;
#ifdef _OPENMP
        sections = omp_get_max_threads();
#pragma omp parallel for
#endif
        for (int p = 0; p < sections; p++) {
            int hh = h / sections;
            int offy = p * hh;
            int x = drand48() * w;
            int y = drand48() * hh;
            int l = rand() % d;
            T fx(0);

            for (int j = 0; j < mutations/sections; j++) {
                int xx = x;
                int yy = y;
                int ll = l;
                mutate(xx, yy, ll, w, hh, d);
                T ed = ED;
                T fy = evaluate(u, problem, xx, yy+offy, ll, ed);

                if (fy > T(0)) {
#ifdef _OPENMP
#pragma omp atomic
#endif
                    rate++;
                    if (u(xx, yy+offy, ll)+ed >= 0. || !positivity) {
                        u(xx, yy+offy, ll) += ed;
                        problem.on_accept(u, xx, yy+offy, ll, ed);
                    }
                }

                if ((fx <= T(0) && fy >= fx)
                    || drand48() < std::min(T(1), fy / fx)) {
                    x = xx;
                    y = yy;
                    l = ll;
                    fx = fy;
                }
            }
        }

        double accept_rate = (double) rate / mutations;
        if (accept_rate < 0.075) {
            ED *= 0.75;
            if (ED <= 1e-3)
                break;
        }
    }
}

template <typename T>
class data {
public:
    img_t<T> f;

    virtual void on_accept(const img_t<T>& u, int x, int y, int l, T ed) = 0;
    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) = 0;
    virtual ~data() { };
};

template <typename T>
class regularizer {
    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) = 0;
public:
    virtual ~regularizer() { };
};

template <typename T, typename D, typename R>
struct problem_t {
    D& data;
    R& regularizer;
    const img_t<T>& f;
    img_t<T> hits;

    problem_t(D& data, R& regularizer, const img_t<T>& f, const img_t<T>& u)
        : data(data), regularizer(regularizer), f(f), hits(u.w, u.h, u.d) {
        data.f = f;
    }

    T energy(img_t<T>& u, int x, int y, int l, T ed=0) {
        return data.evaluate(u, x, y, l, ed) + regularizer.evaluate(u, x, y, l, ed);
    }

    void on_accept(const img_t<T>& u, int x, int y, int l, T ed) {
        hits(x, y, l) += 1.f;
        data.on_accept(u, x, y, l, ed);
    }

};

template <typename T>
struct point_t {
    int x, y;
    T val;
};

template <typename T>
std::vector<point_t<T>> build_kernel(const img_t<T>& K) {
    std::vector<point_t<T>> kernel;
    for (int y = 0; y < K.h; y++) {
        for (int x = 0; x < K.w; x++) {
            if (K(x, y) != 0.)
                kernel.push_back({.x=x - K.w/2, .y=y - K.h/2, .val=K(x, y)});
        }
    }
    return kernel;
}

template <typename T>
class conv_data : public data<T> {
public:

    conv_data<T>(const img_t<T>& K, const img_t<T>& u) :
        K(K), kernel(build_kernel(K)), bounds(u.w, u.h), b(u.w, u.h, u.d)
    {
        bounds.set_value(0);
        slice(bounds, _(K.w/2, -K.w/2-1), _(K.h/2, -K.h/2-1)).map(1);

        // b allows to compute the data term without convolving the image
        reset_b(u);
    }

    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        (void) u; // use 'b' instead

        T sum(0);
        for (auto& k : kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (b.inside(xx, yy) && bounds(xx, yy)) {
                T v = b(xx, yy, l) + ed*k.val;
                sum += std::pow(v - data<T>::f(xx, yy, l), T(2));
            }
        }
        return sum;
    }

    virtual void on_accept(const img_t<T>& u, int x, int y, int l, T ed) override {
        for (auto& k : kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (b.inside(xx, yy))
                b(xx, yy, l) += ed * k.val;
        }
    }

private:
    void reset_b(const img_t<T>& u) {
        b.set_value(0);
        for (int y = 0; y < b.h; y++) {
            for (int x = 0; x < b.w; x++) {
                for (auto& k : kernel) {
                    int xx = x + k.x;
                    int yy = y + k.y;
                    if (b.inside(xx, yy))
                        for (int l = 0; l < u.d; l++)
                            b(xx, yy, l) += u(x, y, l) * k.val;
                }
            }
        }
    }

protected:
    img_t<T> K;
    std::vector<point_t<T>> kernel;
    img_t<uint8_t> bounds;
    img_t<T> b;
};

template <typename T>
class clip_conv_data : public conv_data<T> {
public:
    clip_conv_data(const img_t<T>& K, const img_t<T>& u, T clip) :
        conv_data<T>(K, u), clip(clip)
    {
    }

    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        T sum(0);
        for (auto& k : conv_data<T>::kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (conv_data<T>::b.inside(xx, yy) && conv_data<T>::bounds(xx, yy)) {
                T v = std::min(conv_data<T>::b(xx, yy, l) + ed*k.val, clip);
                sum += std::pow(v - conv_data<T>::f(xx, yy, l), T(2));
            }
        }
        return sum;
    }

private:
    T clip;
};

template <typename T>
class quantize_conv_data : public conv_data<T> {
public:
    quantize_conv_data(const img_t<T>& K, const img_t<T>& u, T quantize) :
        conv_data<T>(K, u), quantize(quantize)
    {
    }

    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        (void) u;
        T sum(0);
        for (auto& k : conv_data<T>::kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (conv_data<T>::b.inside(xx, yy) && conv_data<T>::bounds(xx, yy)) {
                T v = conv_data<T>::b(xx, yy, l) + ed*k.val;
                T e = std::max(std::abs(v - conv_data<T>::f(xx, yy, l)) - T(0.5 / quantize), T(0));
                sum += std::pow(e, T(2));
            }
        }
        return sum;
    }

private:
    T quantize;
};

template <typename T>
class quantize_clip_conv_data : public conv_data<T> {
public:
    quantize_clip_conv_data(const img_t<T>& K, const img_t<T>& u, T quantize, T clip) :
        conv_data<T>(K, u), quantize(quantize), clip(clip)
    {
    }

    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        T sum(0);
        for (auto& k : conv_data<T>::kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (conv_data<T>::b.inside(xx, yy) && conv_data<T>::bounds(xx, yy)) {
                T v = conv_data<T>::b(xx, yy, l) + ed*k.val;
                T v2 = std::min(v, clip);
                T e = std::max(std::abs(v2 - conv_data<T>::f(xx, yy, l)) - T(0.5 / quantize), T(0));
                sum += std::pow(e, T(2));
            }
        }
        return sum;
    }

private:
    T quantize;
    T clip;
};

template <typename T>
class gamma_conv_data : public conv_data<T> {
public:
    gamma_conv_data(const img_t<T>& K, const img_t<T>& u, T gamma) :
        conv_data<T>(K, u), gamma(gamma)
    {
    }

    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        (void) u;
        T sum(0);
        for (auto& k : conv_data<T>::kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (conv_data<T>::b.inside(xx, yy) && conv_data<T>::bounds(xx, yy)) {
                T v = conv_data<T>::b(xx, yy, l) + ed*k.val;
                if (v < T(0))
                    return sum + 9999999.f;
                T e = std::pow(v, gamma) - conv_data<T>::f(xx, yy, l);
                sum += std::pow(e, T(2));
            }
        }
        return sum;
    }

private:
    T gamma;
};

template <typename T>
class gamma_quantize_conv_data : public conv_data<T> {
public:
    gamma_quantize_conv_data(const img_t<T>& K, const img_t<T>& u, T quantize, T gamma) :
        conv_data<T>(K, u), quantize(quantize), gamma(gamma)
    {
    }

    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        (void) u;
        T sum(0);
        for (auto& k : conv_data<T>::kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (conv_data<T>::b.inside(xx, yy) && conv_data<T>::bounds(xx, yy)) {
                T v = conv_data<T>::b(xx, yy, l) + ed*k.val;
                if (v < T(0))
                    return sum + 9999999.f;
                T v2 = std::pow(v, gamma);
                T e = std::max(std::abs(v2 - conv_data<T>::f(xx, yy, l)) - T(0.5 / quantize), T(0));
                sum += std::pow(e, T(2));
            }
        }
        return sum;
    }

private:
    T quantize;
    T gamma;
};

template <typename T>
class gamma_clip_conv_data : public conv_data<T> {
public:
    gamma_clip_conv_data(const img_t<T>& K, const img_t<T>& u, T clip, T gamma) :
        conv_data<T>(K, u), clip(clip), gamma(gamma)
    {
    }

    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        (void) u;
        T sum(0);
        for (auto& k : conv_data<T>::kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (conv_data<T>::b.inside(xx, yy) && conv_data<T>::bounds(xx, yy)) {
                T v = conv_data<T>::b(xx, yy, l) + ed*k.val;
                if (v < T(0))
                    return sum + 9999999.f;
                T v2 = std::min(v, clip);
                T v3 = std::pow(v2, gamma);
                T e = v3 - conv_data<T>::f(xx, yy, l);
                sum += std::pow(e, T(2));
            }
        }
        return sum;
    }

private:
    T clip;
    T gamma;
};

template <typename T>
class gamma_quantize_clip_conv_data : public conv_data<T> {
public:
    gamma_quantize_clip_conv_data(const img_t<T>& K, const img_t<T>& u, T quantize, T clip, T gamma) :
        conv_data<T>(K, u), quantize(quantize), clip(clip), gamma(gamma)
    {
    }

    virtual T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        (void) u;
        T sum(0);
        for (auto& k : conv_data<T>::kernel) {
            int xx = x + k.x;
            int yy = y + k.y;
            if (conv_data<T>::b.inside(xx, yy) && conv_data<T>::bounds(xx, yy)) {
                T v = conv_data<T>::b(xx, yy, l) + ed*k.val;
                if (v < T(0))
                    return sum + 9999999.f;
                T v2 = std::min(v, clip);
                T v3 = std::pow(v2, gamma);
                T e = std::max(std::abs(v3 - conv_data<T>::f(xx, yy, l)) - T(0.5 / quantize), T(0));
                sum += std::pow(e, T(2));
            }
        }
        return sum;
    }

private:
    T quantize;
    T clip;
    T gamma;
};

template <typename T>
static inline constexpr T colorgrad(const img_t<T>& u, int x, int y) {
    T sum(0);
    for (int l = 0; l < u.d; l++) {
        T gx = x >= 1 ? u(x, y, l) - u(x - 1, y, l) : 0.;
        T gy = y >= 1 ? u(x, y, l) - u(x, y - 1, l) : 0.;
        sum += gx*gx+gy*gy;
    }
    return std::sqrt(sum);
}

template <typename T>
static inline constexpr T colorgrad_gamma(const img_t<T>& u, int x, int y, T gamma) {
    T sum(0);
    for (int l = 0; l < u.d; l++) {
        T gx = x >= 1 ? std::pow(u(x, y, l), gamma) - std::pow(u(x - 1, y, l), gamma) : 0.;
        T gy = y >= 1 ? std::pow(u(x, y, l), gamma) - std::pow(u(x, y - 1, l), gamma) : 0.;
        sum += gx*gx+gy*gy;
    }
    return std::sqrt(sum);
}

template <typename T>
class colortv_regularizer : regularizer<T> {
public:
    colortv_regularizer(T lambda, T gamma) :
        lambda(lambda), gamma(gamma)
    {
    }

    T evaluate(img_t<T>& u, int x, int y, int l, T ed=0) override {
        u(x, y, l) += ed;
        T sum(0);
        if (gamma > 0) {
            sum = (colorgrad_gamma(u, x, y, gamma)
                        + (x < u.w-1 ? colorgrad_gamma(u, x+1, y, gamma) : 0.)
                        + (y < u.h-1 ? colorgrad_gamma(u, x, y+1, gamma) : 0.));
        } else {
            sum = (colorgrad(u, x, y)
                        + (x < u.w-1 ? colorgrad(u, x+1, y) : 0.)
                        + (y < u.h-1 ? colorgrad(u, x, y+1) : 0.));
        }
        sum += - std::min(T(0), u(x, y, l)) * 10; // penalize negative intensities
        u(x, y, l) -= ed;
        return lambda * sum;
    }

private:
    T lambda;
    T gamma;
};

template <typename T>
void tv_stochastic_deconvolution(img_t<T>& u, img_t<T> f, const img_t<T>& K, T lambda)
{
    if (K.size <= 1) {
        u = f;
        return;
    }

    // initial guess of 'u' is the blurred image, if nothing else is provided
    if (u.size != f.size) {
        u = f;
    }

    // add padding to the blurred image and to the initial guess
    int hw = K.w/2;
    int hh = K.h/2;
    f = utils::add_padding(f, hw, hh);
    u = utils::add_padding(u, hw, hh);

    colortv_regularizer<T> regularizer(lambda, gamma_tv);

    if (gamma_value > 0.) {
        if (clip > 0. && quantize > 0.) {
            gamma_quantize_clip_conv_data<T> data(K, u, quantize, clip, gamma_value);
            auto p = problem_t<T, decltype(data), decltype(regularizer)>(data, regularizer, f, u);
            stochastic_deconvolution(u, p, f.w, f.h, f.d);
        } else if (quantize > 0.) {
            gamma_quantize_conv_data<T> data(K, u, quantize, gamma_value);
            auto p = problem_t<T, decltype(data), decltype(regularizer)>(data, regularizer, f, u);
            stochastic_deconvolution(u, p, f.w, f.h, f.d);
        } else if (clip > 0.) {
            gamma_clip_conv_data<T> data(K, u, clip, gamma_value);
            auto p = problem_t<T, decltype(data), decltype(regularizer)>(data, regularizer, f, u);
            stochastic_deconvolution(u, p, f.w, f.h, f.d);
        } else {
            gamma_conv_data<T> data(K, u, gamma_value);
            auto p = problem_t<T, decltype(data), decltype(regularizer)>(data, regularizer, f, u);
            stochastic_deconvolution(u, p, f.w, f.h, f.d);
        }
    } else if (clip > 0. && quantize > 0.) {
        quantize_clip_conv_data<T> data_clip(K, u, quantize, clip);
        auto p = problem_t<T, decltype(data_clip), decltype(regularizer)>(data_clip, regularizer, f, u);
        stochastic_deconvolution(u, p, f.w, f.h, f.d);
    } else if (clip > 0.) {
        clip_conv_data<T> data_clip(K, u, clip);
        auto p = problem_t<T, decltype(data_clip), decltype(regularizer)>(data_clip, regularizer, f, u);
        stochastic_deconvolution(u, p, f.w, f.h, f.d);
    } else if (quantize > 0.) {
        quantize_conv_data<T> data_quantize(K, u, quantize);
        auto p = problem_t<T, decltype(data_quantize), decltype(regularizer)>(data_quantize, regularizer, f, u);
        stochastic_deconvolution(u, p, f.w, f.h, f.d);
    } else {
        conv_data<T> data(K, u);
        auto p = problem_t<T, decltype(data), decltype(regularizer)>(data, regularizer, f, u);
        stochastic_deconvolution(u, p, f.w, f.h, f.d);
    }

    u = utils::remove_padding(u, hw, hh);
}

#include "utils.hpp"
template <typename T>
void tv_stochastic_deconvolution_rec(img_t<T>& u, const img_t<T>& f, const img_t<T>& K, T lambda)
{
    int iter = ITERATIONS;
    ITERATIONS *= 2; // use more iterations at coarser levels because it's cheap

    if ((f.w > 5 || f.h > 5) && (K.w > 3 || K.h > 3)) {
        img_t<T> fs;
        img_t<T> Ks;
        // downsample the kernel and the blurry image
        Ks.resize(K.w / 2 + !(K.w/2 % 2), K.h / 2 + !(K.w/2 % 2), K.d); // make sure it stays odd
        utils::downsa2(Ks, K);
        Ks.map(Ks / img::sum<T>(Ks));
        utils::downsa2(fs, f);

        // dive into recursion
        img_t<T> us;
        tv_stochastic_deconvolution_rec(us, fs, Ks, lambda);

        // upsample the result to the same size as f
        u.resize(f);
        utils::upsa2(u, us);
    } else {
        // at low resolution, the deconvolved image is the blurry image (k=delta or f=delta)
        u = f;
    }

    ITERATIONS = iter;

    // deconvolve
    tv_stochastic_deconvolution(u, f, K, lambda);
}

extern "C" int stochastic(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kernel, int kernelsize, float lambda) {

    if (!fast)
        omp_set_num_threads(1);

    img_t<float> f(rx, ry, nchans, fdata);
    img_t<float> K(kernelsize, kernelsize, 1, kernel);
    K.map(K / K.sum());

    float max = f.max();
    if (NORMALIZE_INPUT)
        f.map(f / max);

    if (ed_init <= 0)
        ed_init = (f.max() - f.min()) * 0.03;

    img_t<float> u;
    if (multiscale) {
        tv_stochastic_deconvolution_rec(u, f, K, lambda);
    } else {
        tv_stochastic_deconvolution(u, f, K, lambda);
    }

    if (NORMALIZE_INPUT)
        u.map(u * max);

    // copy u.data.data back to image.fdata
    for (unsigned i = 0; i < rx * ry * nchans; i++)
        fdata[i] = u.data[i];

    return 0;
}
