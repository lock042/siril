#include "image.hpp"
#include "image_expr.hpp"
#include "deconvolution.h"
#include "chelperfuncs.h"
#include "utils.hpp"

// SPA works on mono images: to process a RGB image it must be called 3 times

template <typename T>
void SPA (const img_t<std::complex<T>>& z_init,
          const img_t<T>& h,
          const img_t<T>& mask,
          const img_t<std::complex<T>>& P_x,
          const T sigma_w,
          const int noiter,
          const T reltol,
          img_t<std::complex<T>>& z_hat) {

    assert(z_init.d == 1 && "Only one channel at a time");

    printf("Calculating spectral pre-adaptation...\n");
    int nx = h.w;
    int ny = h.h;
    int NDx = z_init.w;
    int NDy = z_init.h;

    assert(h.w%2 && h.h%2); // Kernel must be odd

    int Lx = std::ceil((nx-1)/2);
    int Ly = std::ceil((ny-1)/2);
    // Algorithm 1 steps 1.1, 1.2
    // Calculate mean of observed part and subtract it
    img_t<std::complex<T>> y_ext(NDx, NDy);
    y_ext.map(z_init * mask);
    printf("Sum of image z_init: %f\n", std::real(z_init.sum()));
    std::complex<T> mu = y_ext.mean() / mask.mean();
    y_ext.map(y_ext - mu);
    printf("Sum of image y_ext2: %f\n", std::real(y_ext.sum()));

    // Algorithm 1 step 1.3
    // Map the blurring kernel into the Fourier domain
    img_t<std::complex<T>> h_ext(NDx, NDy);
    h_ext.set_value(std::complex<T>(0));
    h_ext.pasteinto(h, NDx/2+1-Lx,NDy/2+1-Ly, nx, ny);
    printf("Sum h: %f, sum h_ext: %f\n", h.sum(), std::real(h_ext.sum()));
    h_ext.fftshift();
    img_t<std::complex<T>> H(NDx, NDy);
    H.fft(h_ext);

    // Algorithm 1 step 1.4
    // Compute reference PSD P_z for z using H, P_x and sigma_w
    img_t<std::complex<T>> P_z(NDx, NDy);
    P_z.map(std::abs(H));
    P_z.map(P_z * P_z * P_x + sigma_w * sigma_w);
//checked ok with a numeric test against matlab

    // Algorithm 1 step 1.5
    // Fourier transform of the observed pixels
    img_t<std::complex<T>> Y_f(NDx, NDy);
    Y_f.fft(y_ext);

    // Algorithm 1 step 1.6
    // Then q = inverse FFT of FFT(y_ext) / P_z
    Y_f.map(Y_f / P_z);
    Y_f.ifft(Y_f);
    img_t<std::complex<T>> q(NDx, NDy);
    q.map(std::real(std::complex<T>(0)-Y_f));

    // Algorithm 1 step 1.7
    // Select interpolated pixels as b_y = -q(NOT mask)
    img_t<T> nmask(NDx, NDy);
    nmask.map(T(1) - mask);
    int masksum = static_cast<int>(mask.sum());
    int nmasksum = static_cast<int>(nmask.sum());
    img_t<std::complex<T>> b_y(nmasksum, 1);
    b_y.masktolist(q, nmask);

    // Algorithm 1 Step 2
    // Find the soluton using the Conjugate Gradient method
    // Step 2.1: set z_i and z from the initial guess
    img_t<std::complex<T>> z_i(nmasksum, 1);
    {
        img_t<std::complex<T>> z(z_init); // Temp image to hold z_init - mu
        z.map(z_init - mu);
        z_i.masktolist(z, nmask); // z_i = list of masked pixels in z_init
    }

    // Conjugate Gradient method
    T nb = std::real(std::sqrt(std::inner_product(b_y.begin(), b_y.end(), b_y.begin(), std::complex<T>(0))));

    img_t<std::complex<T>> x(z_i); // x = x0 (= z_i)
    // matrixAMultiplication: r = b - fh(x) = b - matAmult(x)
    img_t<std::complex<T>> r(nmasksum, 1);
    {
        img_t<std::complex<T>> X(NDx, NDy);
        X.expandlisttoimage(x, nmask);
        X.fft(X);
        X.map(X / P_z);
        X.ifft(X);
        r.masktolist(X, nmask);
        r.map(std::real(r));
    }
    r.map(b_y - r);
    img_t<std::complex<T>> p(r);
    std::complex<T> rzold = std::inner_product(r.begin(), r.end(), r.begin(), std::complex<T>(0));

    // Main CG loop
    //
    for (int i = 0 ; i < noiter ; i++) {
        printf("CG iter %d\n", i);

        // matrixAMultiplication: Ap = matAmult(p)
        img_t<std::complex<T>> Ap(nmasksum, 1);
        {
            img_t<std::complex<T>> X(NDx, NDy);
            X.expandlisttoimage(p, nmask);
            X.fft(X);
            X.map(X / P_z);
            X.ifft(X);
            Ap.masktolist(X, nmask);
            Ap.map(std::real(Ap));
        }
        std::complex<T> alpha = rzold / std::inner_product(p.begin(), p.end(), Ap.begin(), std::complex<T>(0));
        x.map(x+alpha*p);
        r.map(r-alpha*Ap);
        T relres = std::abs(std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), std::complex<T>(0))) / nb);
        if (relres < reltol)
            break;
        std::complex<T> rznew = std::inner_product(r.begin(), r.end(), r.begin(), std::complex<T>(0));
        p.map(r + p * (rznew / rzold));
        rzold = rznew;
    } // End conjugate gradient method

    // Populate z_hat with y_ext(mask) + x(nmask)
    //
    {
        img_t<std::complex<T>> X(NDx, NDy);
        X.expandlisttoimage(x, nmask);
        z_hat.map(y_ext + X);
    }
    z_hat.map(z_hat + mu); // Add the mean that was previously subtracted
}

template <typename T>
void extrapol(const img_t<T>& im, const img_t<T>& mask, img_t<std::complex<T>>& imx) {
    int NDx = im.w;
    int NDy = im.h;
    int Nx = mask.w;
    int Ny = mask.h;
    printf("Interpolating smooth boundaries...\n");
    img_t<std::complex<T>> img(NDx, NDy);
    img.map(im * mask);
    img_t<T> fil(NDx, NDy);
    for (int i = 0 ; i < NDx ; i++) {
        T xcoord = i - NDx / 2;
        for (int j = 0 ; j < NDy ; j++) {
            T ycoord = j - NDy / 2;
            fil(i,j) = T(1) / std::pow(((xcoord * xcoord) + (ycoord * ycoord)), T(3.5));
        }
    }
    fil.fftshift();
    fil(0,0) = T(4);
    // The value of fil(0,0) - what should this be - seems to change the strength of the noise -
    // 1000 in the matlab version but that was based on a pixel range of [0,255].
    img_t<std::complex<T>> F(NDx, NDy);
    F.map(fil);
    F.fft(F);
    F.map(std::real(F));
    img_t<std::complex<T>> num(NDx, NDy);
    img_t<std::complex<T>> denom(NDx, NDy);
    num.fft(img);
    num.map(F * num);
    num.ifft(num);
    {
        img_t<std::complex<T>> complexmask(NDx, NDy);
        complexmask.map(mask);
        denom.fft(complexmask);
    }
    denom.map(F * denom);
    denom.ifft(denom);
//    denom.sanitize();
    imx.map(num / denom);
    imx.map((imx * (std::complex<T>(1) - mask)) + (img * mask));
    printf("Smooth boundary estimation complete.\n");
}

template <typename T>
void get_x_psd(img_t<std::complex<T>>& PX_f, T im_range, int Nx, int Ny) {
    printf("Computing Power Spectrum Density...\n");
    im_range = 1.f;
    T sigma_x = (T(30) / T(256)) * im_range; // sigma = 30 is from the original code but does it make sense for astro images?
    T sig2x = sigma_x * sigma_x;
    T rho = T(0.65);
    img_t<T> xc(Nx, Ny);
    img_t<T> yc(Nx, Ny);
    for (int i = 0; i < Nx ; i++) {
        for (int j = 0 ; j < Ny ; j++) {
            xc(i,j) = (i - T(Nx)/2) / Nx; // now equivalent to u in the matlab code
            yc(i,j) = (j - T(Ny)/2) / Ny; // now equivalent to v in the matlab code
        }
    }
    T logrhosquared = std::log(rho) * std::log(rho);
    T fourpisq = T(4) * T(M_PI) * T(M_PI);
    PX_f.set_value(sig2x * T(4) * logrhosquared);
    {
        img_t<T> denom1(xc);
        {
            img_t<T> denom2(yc);
            denom1.map(denom1 * denom1); // u^2
            denom1.map(denom1 * fourpisq); // 4*pi^2*u^2
            denom1.map(denom1 + logrhosquared); // log(rho)^2 + 4*pi^2*u^2
            denom2.map(denom2 * denom2); // v^2
            denom2.map(denom2 * fourpisq); // 4*pi^2*v^2
            denom2.map(denom2 + logrhosquared); // log(rho)^2 + 4*pi^*v^2
            denom1.map(denom1 * denom2); // (log(rho)^2 + 4*pi^2*u^2) + (log(rho)^2 + 4*pi^*v^2)
//          No need to sanitize as this must always be positive because of the log (positive constant)^2 term.
        }
        PX_f.map(PX_f / denom1);
    }
    PX_f.fftshift();
    printf("PSD computation complete.\n");
}

template <typename T>
void prepare_and_run_spa(const img_t<T>& input, const img_t<T>& kernel, img_t<T>& output, T sigma_w, int noiter, int deconv_algo, T lambda, int finaliters, T stopcriterion, int fftw_max_thread) {
    float Nmul = 8.f;
    int Nx = input.w;
    int Ny = input.h;
    int chans = input.d;
    int nx = kernel.w;
    int ny = kernel.h;
    int Lx = ceil((nx-1)/T(2));
    int Ly = ceil((ny-1)/T(2));
    int Lex = ((nx * ceil(1+floor((Nx+ceil(nx/2))/nx))) - Nx) / 2;
    int Ley = ((ny * ceil(1+floor((Ny+ceil(ny/2))/ny))) - Ny) / 2;
    int NDx = Nx + 2 * Lex;
    int NDy = Ny +  2 * Ley;
    T reltol = T(1.e-12);

    // PSD calculation
    //
    T im_power = std::ceil(std::log(input.max())/std::log(T(2)));
    T im_range = std::pow(T(2), im_power);
    img_t<std::complex<T>> PX_f(NDx, NDy);
    get_x_psd(PX_f, im_range, NDx, NDy);

    // Mask setup
    //
    img_t<T> mask;
    {
        img_t<T> tmpmask(Nx, Ny);
        tmpmask.set_value(T(1));
        mask = utils::zero_pad(tmpmask, Lex, Ley);
    }
    printf("Mask populated.\n");

    //Create expanded support y and paste the image monoinput into it
    //
    img_t<T> y = utils::zero_pad(input, Lex, Ley);
    printf("Support constructed.\n");

    // Remove the mean of the masked pixels
    //
    img_t<T> z00(NDx, NDy);
    T mu;
    {
        img_t<T> y00_obsv(NDx, NDy);
        y00_obsv.map(y * mask);
        mu = y00_obsv.mean() / mask.mean();
        img_t<T> nmask(NDx, NDy);
        nmask.map(T(1) - mask);
        z00.map((y00_obsv * mask) + (mu * nmask));
    }
    // z00 has the nmasked (border) area replaced with the mean value of the rest of the image

    // Interpolate initial estimate of smoothed border
    //
    img_t<std::complex<T>> z0(NDx, NDy);
    extrapol(z00, mask, z0);

    // Call SPA
    //
    img_t<T> z_out(NDx, NDy, chans);
    {
        img_t<std::complex<T>> z_SPA(NDx, NDy);
        SPA(z0, kernel, mask, PX_f, sigma_w, noiter, reltol, z_SPA);
        z_out.map(std::real(z_SPA));
    }

    // Change back to float* so as to use the existing C deconv wrappers
    T* fdata = (T*) malloc(z_out.size * sizeof(T));
    for (int i = 0 ; i < z_out.w ; i++)
        for (int j = 0 ; j < z_out.h ; j++)
            for (int c = 0 ; c < z_out.d ; c++)
                fdata[c*z_out.w*z_out.h + j * z_out.w + i] = z_out(i, j, c);

    T* kdata = (T*) malloc(nx * nx * sizeof(T));
    for (int i = 0 ; i < kernel.w ; i++)
        for (int j = 0 ; j < kernel.h ; j++)
            for (int c = 0 ; c < kernel.d ; c++)
                kdata[c*kernel.w*kernel.h + j * kernel.w + i] = kernel(i, j, c);

    // Pass z_SPA to be deconvolved
    switch (deconv_algo) {
        case 0:
            split_bregman(fdata, z_out.w, z_out.h, z_out.d, kdata, nx, lambda, finaliters, fftw_max_thread);
            break;
        case 2:
        case 1:
            richardson_lucy(fdata, z_out.w, z_out.h, z_out.d, kdata, nx, lambda, finaliters, stopcriterion, fftw_max_thread);
            break;
    }

    free(kdata);
    img_t<T> final_result(NDx, NDy, chans, fdata);
    // Trim back to the correct size
    output = utils::remove_padding(final_result, Lex, Ley);
    free(fdata);
}

extern "C" int spectral_pre_adaption(float *fdata, unsigned rx, unsigned ry, unsigned nchans, float *kdata, int kernelsize, float lambda, int max_threads, int deconv_algo, int finaliters, float stopcriterion) {
    double sigma_w;
    if (rx %2 || ry %2) {
        printf("Error, image dimensions must be even! Aborting!\n");
        return -1;
    }
    updatenoise(fdata, rx, ry, nchans, &sigma_w);
    printf("noise est %f\n", sigma_w);
//    sigma_w *= 10; // experimental
    img_t<float> u(rx, ry, nchans, fdata);
    img_t<float> k(kernelsize, kernelsize, 1, kdata);
    img_t<float> out(rx,ry);
    int noiter = 50;
    float Nmul = 8.f;
    deconv_algo = 1;
    if (nchans == 1)
    {
        prepare_and_run_spa(u, k, out, (float) sigma_w, noiter, deconv_algo, lambda, finaliters, stopcriterion, max_threads);
        for (unsigned i = 0; i < rx * ry * nchans; i++)
            fdata[i] = out.data[i];
    } else {
        img_t<float> ycbcr;
        img::rgb2ycbcr(ycbcr, u);
        img_t<float> f(rx, ry, 1);
        for (int i = 0 ; i < f.w * f.h ; i++)
            f[i] = ycbcr[i*3];
        prepare_and_run_spa(f, k, out, (float) sigma_w, noiter, deconv_algo, lambda, finaliters, stopcriterion, max_threads);
        for (int i = 0 ; i < out.size ; i++)
            ycbcr[i*3] = out[i];
        img_t<float> final_out;
        img::ycbcr2rgb(final_out,ycbcr);
        // copy u.data.data back to image.fdata
        for (unsigned i = 0; i < rx * ry * nchans; i++)
            fdata[i] = final_out.data[i];
    }

    return 0;
}
