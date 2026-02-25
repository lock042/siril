/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stddef.h>
#include <math.h>
#include "core/siril.h"
#include "algos/anscombe.h"

void generalized_anscombe_array(float *x, const float mu, const float sigma, const float gain, const size_t ndata) {
    /*
    Compute the generalized anscombe variance stabilizing transform,
    which assumes that the data provided to it is a mixture of poisson
    and gaussian noise. Applies to an array of floats.
    The input signal  z  is assumed to follow the Poisson-Gaussian noise model
        x = gain * p + n
    where gain is the camera gain and mu and sigma are the read noise
    mean and standard deviation.
    We assume that x contains only positive values.  Values that are
    less than or equal to 0 are ignored by the transform.
    Note, this transform will show some bias for counts less than
    about 20.
    */
    float addterm = powf(gain, 2.f) * 0.375f + sigma * sigma - gain * mu;
    float factor = (2.f / gain);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
    for (size_t i = 0; i < ndata; i++) {
        float y = gain * x[i] + addterm;
        // Clamp to zero before taking the square root.
        x[i] = factor * sqrtf(max(y, 0.f));
    }
}

void inverse_generalized_anscombe_array(float *x, const float mu, const float sigma, const float gain, const size_t ndata) {
    /*
    Applies the closed-form approximation of the exact unbiased
    inverse of Generalized Anscombe variance-stabilizing
    transformation to an array of floats.
    The input signal x is transform back into a Poisson random variable
    based on the assumption that the original signal from which it was
    derived follows the Poisson-Gaussian noise model:
        x = gain * p + n
    where gain is the camera gain and mu and sigma are the read noise
    mean and standard deviation.
    Roference: M. Makitalo and A. Foi, "Optimal inversion of the
    generalized Anscombe transformation for Poisson-Gaussian noise",
    IEEE Trans. Image Process., doi:10.1109/TIP.2012.2202675
    */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < ndata ; i++) {
        float test = max(x[i], 1.0);
        float exact_inverse = ( 0.25f * powf(test, 2.f) +
                                0.25f * sqrtf(1.5f) * powf(test, -1.f) -
                                1.375f * powf(test, -2.f) +
                                0.625f * sqrtf(1.5f) * powf(test, -3.f) -
                                0.125f - powf(sigma, 2.f) );
        // Clamp to zero
        exact_inverse = max(0.f, exact_inverse);
        exact_inverse *= gain;
        exact_inverse += mu;
        if (exact_inverse != exact_inverse) // Catch NaNs
            exact_inverse = 0.f;
        x[i] = exact_inverse;
    }
}

float generalized_anscombe(const float x, const float mu, const float sigma, const float gain) {
    /*
    Compute the generalized anscombe variance stabilizing transform,
    which assumes that the data provided to it is a mixture of poisson
    and gaussian noise. Applies to a single float and returns a float.
    The input signal  z  is assumed to follow the Poisson-Gaussian noise model
        x = gain * p + n
    where gain is the camera gain and mu and sigma are the read noise
    mean and standard deviation.
    We assume that x contains only positive values.  Values that are
    less than or equal to 0 are ignored by the transform.
    Note, this transform will show some bias for counts less than
    about 20.
    */
    float y = gain * x + powf(gain, 2.f) * 0.375f + sigma*sigma - gain*mu;
    // Clamp to zero before taking the square root.
    return (2.f / gain) * sqrtf(max(y, 0.f));
}

float inverse_generalized_anscombe(const float x, const float mu, const float sigma, const float gain) {
    /*
    Applies the closed-form approximation of the exact unbiased
    inverse of Generalized Anscombe variance-stabilizing
    transformation. Applies to a single float and returns a float.
    The input signal x is transform back into a Poisson random variable
    based on the assumption that the original signal from which it was
    derived follows the Poisson-Gaussian noise model:
        x = gain * p + n
    where gain is the camera gain and mu and sigma are the read noise
    mean and standard deviation.
    Roference: M. Makitalo and A. Foi, "Optimal inversion of the
    generalized Anscombe transformation for Poisson-Gaussian noise",
    IEEE Trans. Image Process., doi:10.1109/TIP.2012.2202675
    */
    float test = max(x, 1.0);
    float exact_inverse = ( 0.25f * powf(test, 2.f) +
                            0.25f * sqrtf(1.5f)*powf(test, -1.f) -
                            1.375f * powf(test, -2.f) +
                            0.625f * sqrtf(1.5f) * powf(test, -3.f) -
                            0.125f - powf(sigma, 2.f) );
    // Clamp to zero
    exact_inverse = max(0.f, exact_inverse);
    exact_inverse *= gain;
    exact_inverse += mu;
    if (exact_inverse != exact_inverse) // Catch NaNs
        exact_inverse = 0.f;
    return exact_inverse;
}

float anscombe(float x) {
    /*
    Compute the anscombe variance stabilizing transform.
    the input   x   is noisy Poisson-distributed data
    the output  fx  has variance approximately equal to 1.
    Reference: Anscombe, F. J. (1948), "The transformation of Poisson,
    binomial and negative-binomial data", Biometrika 35 (3-4): 246-254
    */

    return 2.f * sqrtf(x + 3.f / 8.f);
}

float exact_unbiased_inverse_anscombe(float z) {
    /*
    Compute the inverse transform using an approximation of the exact unbiased inverse.
    Reference: Makitalo, M., & Foi, A. (2011). A closed-form approximation of the exact
    unbiased inverse of the Anscombe variance-stabilizing transformation.
    Image Processing.
    */

    return (0.25f * powf(z, 2.f) +
            0.25f * sqrtf(1.5f) * powf(z, -1.f) -
            1.375f * powf(z, -2.f) +
            0.625f * sqrtf(1.5f) * powf(z, -3.f) - 0.125f);
}
