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
#include "opencv/opencv.h"

#include "filters/bm3d/call_bm3d.h" // We use some of the same fits handling functions
#include "filters/da3d/Image.hpp"
#include "filters/da3d/Utils.hpp"
#include "filters/da3d/DA3D.hpp"
#include "Utilities.h"
#include "NlBayes.h"
#include "LibImages.h"

extern "C" {
#include "core/proto.h"
#include "core/OS_utils.h"
#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"
#include "algos/statistics.h"
#include "core/processing.h"
}

using namespace std;
using std::cerr;
using std::endl;

using utils::isMonochrome;
using utils::makeMonochrome;
using NlBayes::runNlBayes;
using da3d::Image;
using da3d::DA3D;

extern "C" int do_nlbayes(fits *fit, float modulation, int da3d) {
    // Parameters
    const unsigned width = fit->naxes[0];
    const unsigned height = fit->naxes[1];
    const unsigned nchans = fit->naxes[2];
    const unsigned npixels = width * height;
    float norm = 1.f;
    if (fit->type == DATA_USHORT)
      norm = USHRT_MAX_SINGLE;
    float invnorm = 1 / norm;

    // Measure background noise
    float fSigma = 0.f;
    for (size_t chan = 0 ; chan < nchans ; chan++) {
      imstats* stat = statistics(NULL, -1, fit, chan, NULL, STATS_SIGMEAN, MULTI_THREADED);
      fSigma += stat->bgnoise / norm;
      free_stats(stat);
    }
    fSigma /= nchans;
    siril_log_message(_("NL-Bayes auto parametrisation: measured background noise level is %f\n"),fSigma);

	// The useArea bools set the paste trick in both halves of the Nl-Bayes process
    const bool useArea1 = true;
    const unsigned useArea2 = true;
    const bool verbose = FALSE;

    float *bgr_f, *bgr_fout;

    if (fit->type == DATA_FLOAT) {
      if (nchans == 3) {
        bgr_f = fits_to_bgrbgr_float(fit);
      } else {
        bgr_f = (float*) calloc(npixels, sizeof(float));
        for (unsigned i=0; i<npixels; i++) {
          bgr_f[i] = fit->fdata[i];
        }
      }
    } else {
      if (nchans == 3) {
        bgr_f = fits_to_bgrbgr_wordtofloat(fit);
      } else {
        bgr_f = (float*) calloc(npixels, sizeof(float));
        for (unsigned i=0; i<npixels; i++) {
          bgr_f[i] = (float) fit->data[i] * invnorm;
        }
      }
    }

    vector<float> bgr_v { bgr_f, bgr_f + width * height * nchans };
    vector<float> basic, bgr_vout;
    ImageSize imSize;
    imSize.width = width;
    imSize.height = height;
    imSize.nChannels = nchans;
    imSize.wh = width * height;
    imSize.whc = width * height * nchans;

    if(!get_thread_run()) {
      return EXIT_FAILURE;
    }

    if (runNlBayes(bgr_v, basic, bgr_vout, imSize, useArea1, useArea2, fSigma, verbose) != EXIT_SUCCESS)
        return EXIT_FAILURE;

      set_progress_bar_data("NL-Bayes denoising...", 0.0);

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
      Image output = DA3D(input, guide, fSigma);
      bgr_da3dout = output.data();
      memcpy(bgr_fout, bgr_da3dout, height * width * nchans * sizeof(float));
    }

    // Convert output from bgrbgr back to planar rgb and put back into fit
    if (fit->type == DATA_FLOAT) {
      if (nchans == 3) {
        bgrbgr_float_to_fits(fit, bgr_fout, modulation);
      } else {
        for (unsigned i=0; i<npixels; i++) {
          fit->fdata[i] = (1.f-modulation) * fit->fdata[i] + modulation * bgr_fout[i];
        }
      }
    } else {
      if (nchans == 3) {
        bgrbgr_float_to_word_fits(fit, bgr_fout, modulation);
      } else {
         for (unsigned i=0; i<npixels; i++) {
          fit->data[i] = round_to_WORD(USHRT_MAX * ((1.f-modulation) * (fit->data[i] * invnorm) + modulation * bgr_fout[i]));
        }
      }
    }
    return EXIT_SUCCESS;
}

