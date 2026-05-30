/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

#include "filters/da3d/DA3D.hpp"
#include "algos/img_t/image_ops.hpp"
#include "algos/img_t/img_fits.hpp"
#include "Utilities.h"
#include "NlBayes.h"
#include "LibImages.h"

extern "C" {
#include "core/proto.h"
#include "core/siril.h"
#include "core/gui_iface.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "io/image_format_fits.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "algos/anscombe.h"
#include "core/processing.h"
#include "filters/cosmetic_correction.h"
}

using namespace std;
using std::cerr;
using std::endl;
using std::transform;
using std::bind;
using std::multiplies;
using std::divides;

using NlBayes::runNlBayes;
using da3d::DA3D;

extern "C" int do_nlbayes(fits *fit, const float modulation, unsigned sos, int da3d, const float rho, const gboolean do_anscombe) {
    // Parameters
    const unsigned width = fit->naxes[0];
    const unsigned height = fit->naxes[1];
    const unsigned nchans = fit->naxes[2];
    const unsigned whc = width * height * nchans;
    float fSigma = 0.f;
    float lastfSigma = 0.f;
    double intermediate_bgnoise;

    if (sos < 1)
      sos = 1;

    // The useArea bools set the paste trick in both halves of the Nl-Bayes function
    const bool useArea1 = true;
    const unsigned useArea2 = true;
    const bool verbose = false;

    // Build the source image once from the fits (normalise USHORT, clip to
    // [0,1]); the whole SOS loop then runs entirely in img_t<float>, so there
    // is no raw-buffer / vector / img_t juggling.
    img_t<float> orig = imgops::from_fits(fit);
    sos_update_noise_float(orig.data.data(), width, height, nchans, &intermediate_bgnoise);
    lastfSigma = (float) intermediate_bgnoise;

    img_t<float> vout = orig;                 // running denoised result
    img_t<float> v(width, height, nchans);    // SOS-strengthened input

    gui_iface.set_progress(0.0, _("NL-Bayes denoising..."));

    if (!processing_should_continue()) {
        siril_log_debug("do_nlbayes: processing_should_continue() returned FALSE\n");
        return EXIT_FAILURE;
    }
    if (do_anscombe) {
      sos = 1; // SOS doesn't make sense with VST as it would be adding non-Gaussian noise back into an AWGN denoising algorithm.
      da3d = false; // DA3D doesn't produce good results in combination with the Anscombe VST, so we disable it.
    }

    // SOS iteration loop
    for (unsigned iter = 0; iter < sos; iter++) {

      // Strengthen the previous result by mixing back a fraction of the original noisy image
      for (unsigned i = 0; i < whc; i++)
        v.data[i] = rho * vout.data[i] + (1.f - rho) * orig.data[i];
      // Update image noise measurement
      sos_update_noise_float(v.data.data(), width, height, nchans, &intermediate_bgnoise);
      fSigma = (float) intermediate_bgnoise;
      if (fSigma > lastfSigma * 1.01f) {
        // Note we only check this on the first iteration: if the first iteration converges then subsequent iterations appear to converge reliably, however the noise level on successive iterations does not necessarily decrease monotonically.
        siril_log_error(_("Error: SOS is not converging. Try a smaller value of rho.\n"));
        vout = orig;
        break;
      }
      if (iter == 0)
        siril_log_message(_("NL-Bayes auto parametrisation: measured background noise level is %.3e\n"), fSigma);

      if (do_anscombe && iter == 0) {
        siril_log_message(_("Applying Anscombe VST\n"));
        transform(v.data.begin(), v.data.end(), v.data.begin(),
              bind(multiplies<float>(), std::placeholders::_1, 65536.f));
        generalized_anscombe_array(v.data.data(), 0.f, 0.f, 1.f, whc);
        float norm = 65536.f * 2.f * sqrtf(11.f / 8.f);
        transform(v.data.begin(), v.data.end(), v.data.begin(),
              bind(divides<float>(), std::placeholders::_1, norm));
        sos_update_noise_float(v.data.data(), width, height, nchans, &intermediate_bgnoise);
        fSigma = (float) intermediate_bgnoise;
      }

      // Operate the NL-Bayes algorithm directly on the img_t buffers.
      img_t<float> imBasic, imFinal;
      if (runNlBayes(v, imBasic, imFinal, useArea1, useArea2, fSigma, verbose) != EXIT_SUCCESS) {
        siril_log_debug("do_nlbayes: runNlBayes returned an error code\n");
        return EXIT_FAILURE;
      }
      vout = std::move(imFinal);

      if (do_anscombe && iter == 0) {
        siril_log_message(_("Applying exact unbiased inverse Anscombe VST\n"));
        float norm = 65536.f * 2.f * sqrtf(11.f / 8.f);
        transform(vout.data.begin(), vout.data.end(), vout.data.begin(),
              bind(multiplies<float>(), std::placeholders::_1, norm));
        inverse_generalized_anscombe_array(vout.data.data(), 0.f, 0.f, 1.f, whc);
        transform(vout.data.begin(), vout.data.end(), vout.data.begin(),
              bind(divides<float>(), std::placeholders::_1, 65536.f));
      }

      // Subtraction step to produce input for next iteration. Not performed on the final iteration.
      if (iter + 1 < sos) {
        for (unsigned i = 0; i < whc; i++) {
          vout.data[i] *= 1.f / (1.f - rho);
          vout.data[i] -= v.data[i] * (rho / (1.f - rho));
        }
        if (sos > 1)
          siril_log_message(_("SOS iteration %d of %d complete\n"), iter + 1, sos);
      }
    }

    // Clip in case NL-Bayes pushed a saturated pixel slightly over 1.0.
    for (unsigned i = 0; i < whc; i++)
      vout.data[i] = max(min(1.f, vout.data[i]), 0.f);

    // Carry out final-stage DA3D denoising if required
    if (da3d != 0) {
      siril_log_message(_("DA3D final-stage denoising...\n"));
#ifndef _OPENMP
      siril_log_message(_("OpenMP not available. The DA3D algorithm will run in a single thread.\n"));
#endif
      img_t<float> input = orig;   // original (noisy) image
      img_t<float> guide = vout;   // NL-Bayes result, used as the DA3D guide
      // DA3D doesn't work if a color image has monochromatic noise
      if (input.d > 1 && imgops::is_monochrome(input)) {
        siril_log_error(_("Warning: input color image has monochromatic noise! Converting to monochrome."));
        input = imgops::make_monochrome(input);
        guide = imgops::make_monochrome(guide);
      }
      int retval = 0;
      img_t<float> output = DA3D(retval, input, guide, lastfSigma);
      if (retval != 0) {
        siril_log_debug("do_nlbayes: DA3D returned an error code\n");
        return EXIT_FAILURE;
      }
      vout = std::move(output);
    }

    // Final noise measurement
    sos_update_noise_float(vout.data.data(), width, height, nchans, &intermediate_bgnoise);
    fSigma = (float) intermediate_bgnoise;
    siril_log_message(_("NL-Bayes output: measured background noise level is %.3e\n"), fSigma);

    // Blend the denoised result back into the fits.
    imgops::to_fits(vout, fit, modulation);

    return EXIT_SUCCESS;
}
