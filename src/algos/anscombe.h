#ifndef ANSCOMBE_H_
#define ANSCOMBE_H_

void generalized_anscombe_array(float *x, const float mu, const float sigma, const float gain, const size_t ndata);

void inverse_generalized_anscombe_array(float *x, const float mu, const float sigma, const float gain, const size_t ndata);

float generalized_anscombe(const float x, const float mu, const float sigma, const float gain);

float inverse_generalized_anscombe(const float x, const float mu, const float sigma, const float gain);

float anscombe(float x);

float exact_unbiased_inverse_anscombe(float x);

#endif
