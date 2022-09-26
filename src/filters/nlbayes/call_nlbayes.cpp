#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>
#include <iterator>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <utility>
#include <functional>
#include "opencv/opencv.h"

#include "filters/da3d/Image.hpp"
#include "filters/da3d/Utils.hpp"
#include "filters/da3d/DA3D.hpp"
#include "Utilities.h"
#include "NlBayes.h"
#include "LibImages.h"

extern "C" {
#include "core/proto.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "io/image_format_fits.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "algos/anscombe.h"
#include "core/processing.h"
#include "filters/cosmetic_correction.h"
}

#define ACTIVATE_NULLCHECK_FLOAT 1

using namespace std;
using std::cerr;
using std::endl;
using std::transform;
using std::bind;
using std::multiplies;
using std::divides;

using utils::isMonochrome;
using utils::makeMonochrome;
using NlBayes::runNlBayes;
using da3d::Image;
using da3d::DA3D;

extern "C" int do_nlbayes(fits *fit, const float modulation, unsigned sos, int da3d, const float rho, const gboolean do_anscombe, const gboolean do_cosme) {
    // Parameters
    const unsigned width = fit->naxes[0];
    const unsigned height = fit->naxes[1];
    const unsigned nchans = fit->naxes[2];
    const unsigned npixels = width * height;
    float fSigma = 0.f;
    float lastfSigma = 0.f;
    double intermediate_bgnoise;

    float normalize = 1.f;
    if (fit->type == DATA_USHORT)
    	normalize = USHRT_MAX_SINGLE;
    float invnorm = 1 / normalize;
    if (sos < 1)
      sos = 1;

	// The useArea bools set the paste trick in both halves of the Nl-Bayes function
    const bool useArea1 = true;
    const unsigned useArea2 = true;
    const bool verbose = FALSE;

    float *bgr_f = (float*) calloc(npixels * nchans, sizeof(float));
    float *bgr_fout;

    // Carry out cosmetic correction at the start, if selected
    if (do_cosme == TRUE)
      denoise_hook_cosmetic(fit);

    if (fit->type == DATA_FLOAT) {
        for (size_t i = 0; i < npixels * nchans; i++)
          bgr_f[i] = fit->fdata[i];
    } else {
        for (size_t i = 0; i < npixels * nchans; i++)
          bgr_f[i] = (float) fit->data[i] * invnorm;
    }

    // Measure image noise using the custom wrapper to FnNoise1_float in quantize.c
    // Initial noise measurement
    sos_update_noise_float(bgr_f, width, height, nchans, &intermediate_bgnoise);
    lastfSigma = (float) intermediate_bgnoise;

    vector<float> bgr_v_orig { bgr_f, bgr_f + width * height * nchans };
    vector<float> bgr_vout = bgr_v_orig;
    vector<float> bgr_v(width * height * nchans, 0.f);
    vector<float> basic;
    ImageSize imSize;
    imSize.width = width;
    imSize.height = height;
    imSize.nChannels = nchans;
    imSize.wh = width * height;
    imSize.whc = width * height * nchans;

    set_progress_bar_data("NL-Bayes denoising...", 0.0);

    if(!get_thread_run()) {
      return EXIT_FAILURE;
    }
    if (do_anscombe) {
      sos = 1; // SOS doesn't make sense with VST as it would be adding non-Gaussian noise back into an AWGN denoising algorithm.
      da3d = false; // DA3D doesn't produce good results in combination with the Anscombe VST, so we disable it.
    }

    // SOS iteration loop
    for (unsigned iter=0; iter < sos; iter++) {

      // Strengthen result of previous iteration bgr_v by mixing back a fraction of the original noisy image
      for (unsigned i = 0; i < bgr_v_orig.size(); i++)
        bgr_v[i] = (rho * bgr_vout[i] + (1.f - rho) * bgr_v_orig[i]);
      // Update image noise measurement
      sos_update_noise_float(bgr_v.data(), width, height, nchans, &intermediate_bgnoise);
      fSigma = (float) intermediate_bgnoise;
      if (fSigma > lastfSigma * 1.01f) {
        // Note we only check this on the first iteration: if the first iteration converges then subsequent iterations appear to converge reliably, however the noise level on successive iterations does not necessarily decrease monotonically.
        siril_log_color_message(_("Error: SOS is not converging. Try a smaller value of rho.\n"),"red");
        bgr_vout = bgr_v_orig;
        break;
      }
      if (iter == 0)
        siril_log_message(_("NL-Bayes auto parametrisation: measured background noise level is %.3e\n"),fSigma);

      if (do_anscombe && iter == 0) {
        siril_log_message(_("Applying Anscombe VST\n"));
        if (nchans == 1) {
          // Mono images
          transform(bgr_v.begin(), bgr_v.end(), bgr_v.begin(),
               bind(multiplies<float>(), std::placeholders::_1, 65536.f));
          sos_update_noise_float(bgr_v.data(), width, height, nchans, &intermediate_bgnoise);
          fSigma = (float) intermediate_bgnoise;
          generalized_anscombe_array(bgr_v.data(), 0.f, fSigma, 1.f, imSize.whc);
          sos_update_noise_float(bgr_v.data(), width, height, nchans, &intermediate_bgnoise);
          fSigma = (float) intermediate_bgnoise;
        } else {
          // RGB images
          transform(bgr_v.begin(), bgr_v.end(), bgr_v.begin(),
               bind(multiplies<float>(), std::placeholders::_1, 65536.f));
          generalized_anscombe_array(bgr_v.data(), 0.f, 0.f, 1.f, imSize.whc);
          float norm = 65536.f * 2.f * sqrtf(11.f / 8.f);
          transform(bgr_v.begin(), bgr_v.end(), bgr_v.begin(),
               bind(divides<float>(), std::placeholders::_1, norm));
          sos_update_noise_float(bgr_v.data(), width, height, nchans, &intermediate_bgnoise);
          fSigma = (float) intermediate_bgnoise;
        }
      }
      // Operate the NL-Bayes algorithm
      if (runNlBayes(bgr_v, basic, bgr_vout, imSize, useArea1, useArea2, fSigma, verbose) != EXIT_SUCCESS)
        return EXIT_FAILURE;

      if (do_anscombe && iter == 0) {
        siril_log_message(_("Applying exact unbiased inverse Anscombe VST\n"));
        if (nchans == 1) {
          // Inverse transform for mono images
          inverse_generalized_anscombe_array(bgr_vout.data(), 0.f, fSigma, 1.f, imSize.whc);
          transform(bgr_vout.begin(), bgr_vout.end(), bgr_vout.begin(),
               bind(divides<float>(), std::placeholders::_1, 65536.f));

        } else {
          float norm = 65536.f * 2.f * sqrtf(11.f / 8.f);
          transform(bgr_vout.begin(), bgr_vout.end(), bgr_vout.begin(),
               bind(multiplies<float>(), std::placeholders::_1, norm));
          inverse_generalized_anscombe_array(bgr_vout.data(), 0.f, 0.f, 1.f, imSize.whc);
          transform(bgr_vout.begin(), bgr_vout.end(), bgr_vout.begin(),
               bind(divides<float>(), std::placeholders::_1, 65536.f));
        }
      }

      // Subtraction step to produce input for next iteration. Not performed on the final iteration.
      if (iter+1 < sos) {
        for (unsigned i = 0; i < bgr_v.size(); i++) {
          bgr_vout[i] *= 1.f / (1.f - rho);
          bgr_vout[i] -= bgr_v[i] * (rho / (1.f - rho));
        }
      if (sos > 1)
        siril_log_message(_("SOS iteration %d of %d complete\n"), iter+1, sos);
      }
    }

    bgr_fout = bgr_vout.data();
    float *bgr_da3dout;

    // Carry out final-stage DA3D denoising if required
    if (da3d != 0) {
      siril_log_message(_("DA3D final-stage denoising...\n"));
#ifndef _OPENMP
      siril_log_message(_("OpenMP not available. The DA3D algorithm will run in a single thread.\n"));
#endif
      Image input(bgr_f, height, width, nchans);
      Image guide(bgr_fout, height, width, nchans);
      // DA3D doesn't work if a color image has monochromatic noise
      if (input.channels()>1 && isMonochrome(input)) {
        siril_log_color_message(_("Warning: input color image has monochromatic noise! Converting to monochrome."), "red");
        input = makeMonochrome(input);
        guide = makeMonochrome(guide);
      }
      int retval = 0;
      Image output = DA3D(retval, input, guide, lastfSigma);
      if (retval != 0)
        return EXIT_FAILURE;
      bgr_da3dout = output.data();
      memcpy(bgr_fout, bgr_da3dout, height * width * nchans * sizeof(float));
    }

    // Final noise measurement
    sos_update_noise_float(bgr_fout, width, height, nchans, &intermediate_bgnoise);
    fSigma = (float) intermediate_bgnoise;
    siril_log_message(_("NL-Bayes output: measured background noise level is %.3e\n"),fSigma);

    // Convert output from bgrbgr back to planar rgb and put back into fit
    if (fit->type == DATA_FLOAT) {
        for (unsigned i = 0; i < npixels * nchans; i++)
          fit->fdata[i] = (1.f-modulation) * fit->fdata[i] + modulation * bgr_fout[i];
    } else {
       for (unsigned i = 0; i < npixels * nchans; i++)
          fit->data[i] = round_to_WORD(USHRT_MAX * ((1.f-modulation) * (fit->data[i] * invnorm) + modulation * bgr_fout[i]));
    }

    return EXIT_SUCCESS;
}

