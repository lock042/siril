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

void bgrbgr_float_to_fits(fits *image, float *bgrbgr, float modulation) {
  unsigned npixels = image->naxes[0] * image->naxes[1];

  for (unsigned i = 0; i < npixels; i++) {
    image->fpdata[BLAYER][i] = (1.f-modulation) * image->fpdata[BLAYER][i] + modulation * bgrbgr[i*3];
	image->fpdata[GLAYER][i] = (1.f-modulation) * image->fpdata[GLAYER][i] + modulation * bgrbgr[i*3+1];
	image->fpdata[RLAYER][i] = (1.f-modulation) * image->fpdata[RLAYER][i] + modulation * bgrbgr[i*3+2];
  }
}

void bgrbgr_float_to_word_fits(fits *image, float *bgrbgr, float modulation) {
	size_t ndata = image->rx * image->ry * 3;
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		image->pdata[BLAYER][j] = round_to_WORD(((modulation * bgrbgr[i + 0]) + (1.f - modulation) * image->pdata[BLAYER][j] / USHRT_MAX_SINGLE) * USHRT_MAX);
		image->pdata[GLAYER][j] = round_to_WORD(((modulation * bgrbgr[i + 1]) + (1.f - modulation) * image->pdata[GLAYER][j] / USHRT_MAX_SINGLE) * USHRT_MAX);
		image->pdata[RLAYER][j] = round_to_WORD(((modulation * bgrbgr[i + 2]) + (1.f - modulation) * image->pdata[RLAYER][j] / USHRT_MAX_SINGLE) * USHRT_MAX);
	}
}

float *fits_to_bgrbgr_wordtofloat(fits *image) {
    size_t ndata = image->rx * image->ry * 3;
    float invnorm = 1 / USHRT_MAX_SINGLE;
	float *bgrbgr = (float *)malloc(ndata * sizeof(float));
	if (!bgrbgr) { PRINT_ALLOC_ERR; return NULL; }
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		bgrbgr[i + 0] = (float) image->pdata[BLAYER][j] * invnorm;
		bgrbgr[i + 1] = (float) image->pdata[GLAYER][j] * invnorm;
		bgrbgr[i + 2] = (float) image->pdata[RLAYER][j] * invnorm;
	}
	return bgrbgr;
}

extern "C" int do_nlbayes(fits *fit, float modulation, unsigned sos, int da3d) {
    // Parameters
    const unsigned width = fit->naxes[0];
    const unsigned height = fit->naxes[1];
    const unsigned nchans = fit->naxes[2];
    const unsigned npixels = width * height;
    float norm = 1.f;
    if (fit->type == DATA_USHORT)
      norm = USHRT_MAX_SINGLE;
    float invnorm = 1 / norm;
    if (sos < 1)
      sos = 1;

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

    float rho = 0.3; // Proportion of original noise folded back in is (1-rho)
    // SOS iteration loop
    for (unsigned i=0; i < sos; i++) {

      // Strengthen result of previous iteration bgr_v by mixing back a fraction of the original noisy image
      for (unsigned i = 0; i < bgr_v_orig.size(); i++)
        bgr_v[i] = (rho * bgr_vout[i] + (1.f - rho) * bgr_v_orig[i]);

      // Operate the NL-Bayes algorithm
      if (runNlBayes(bgr_v, basic, bgr_vout, imSize, useArea1, useArea2, fSigma, verbose) != EXIT_SUCCESS)
        return EXIT_FAILURE;

      // Subtraction step to produce input for next iteration. Not performed on the final iteration.
      if (i+1 < sos) {
        for (unsigned i = 0; i < bgr_v.size(); i++) {
          bgr_vout[i] *= 1.f / (1.f - rho);
          bgr_vout[i] -= bgr_v[i] * (rho / (1.f - rho));
          // Measure background noise
          fSigma = 0.f;
          for (size_t chan = 0 ; chan < nchans ; chan++) {
            imstats* stat = statistics(NULL, -1, fit, chan, NULL, STATS_SIGMEAN, MULTI_THREADED);
            fSigma += stat->bgnoise / norm;
            free_stats(stat);
          }
          fSigma /= nchans;
        }
      if (sos > 1)
        siril_log_message(_("SOS iteration %d of %d...\n"), i+1, sos);
      siril_log_message(_("NL-Bayes auto parametrisation: measured background noise level is %f\n"),fSigma);
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
      Image output = DA3D(retval, input, guide, fSigma);
      if (retval != 0)
        return EXIT_FAILURE;
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

