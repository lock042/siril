#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>

#include "bm3d.h"
#include "bm3d_utilities.h"
#include "opencv/opencv.h"

extern "C" {
#include "core/proto.h"
#include "core/OS_utils.h"
#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"
#include "algos/statistics.h"
}

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3
#define DCT       4
#define BIOR      5
#define HADAMARD  6
#define NONE      7

using namespace std;

void bgrbgr_float_to_fits(fits *image, float *bgrbgr, float modulation) {
	size_t ndata = image->rx * image->ry * 3;
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		image->fpdata[BLAYER][j] = (1.f-modulation) * image->fpdata[BLAYER][j] + modulation * bgrbgr[i + 0];
		image->fpdata[GLAYER][j] = (1.f-modulation) * image->fpdata[GLAYER][j] + modulation * bgrbgr[i + 1];
		image->fpdata[RLAYER][j] = (1.f-modulation) * image->fpdata[RLAYER][j] + modulation * bgrbgr[i + 2];
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

extern "C" int do_bm3d(fits *fit, float modulation) {
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
    siril_log_message(_("Auto parametrisation: measured background noise level is %f\n"),fSigma);

	// The algorithm parameters set below are the optimal parameters discussed in
	// (Lebrun, 2012): http://www.ipol.im/pub/art/2012/l-bm3d/
	const bool useSD_1 = FALSE;
    const bool useSD_2 = FALSE;
    const unsigned tau_2D_hard = BIOR;
    const unsigned tau_2D_wien = DCT;
    const unsigned color_space = OPP; //OPP; // YUV, OPP or YCBCR
    const int patch_size = 0;
    const int nb_threads = 0;
    const bool verbose = FALSE;

    float *bgr_f;

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

    // Work out chunks required to avoid OOM
    float memGB = (float) (get_available_memory() / 1000000000);
    float imgmemMpix = npixels / 1000000;
    unsigned numchunks = (unsigned) imgmemMpix / (memGB / 6.f);
    if (numchunks < 1)
      numchunks = 1;

    // Cut into manageable chunks to avoid OOM
    vector<vector<float> > chunk_noisy(numchunks);
    vector<vector<float> > chunk_basic(numchunks);
    vector<vector<float> > chunk_denoised(numchunks);
    vector<unsigned> w_table(numchunks);
    vector<unsigned> h_table(numchunks);
    sub_divide(bgr_v, chunk_noisy, w_table, h_table, width, height, nchans, 32, TRUE); // 32 is nWien * 2;

    // Run bm3d on each chunk in turn.
    for (unsigned i = 0; i < numchunks ; i++) {

      if(!get_thread_run()) {
        return EXIT_FAILURE;
      }

      if (run_bm3d(fSigma, chunk_noisy[i], chunk_basic[i], chunk_denoised[i],
          w_table[i], h_table[i], nchans, useSD_1, useSD_2, tau_2D_hard, tau_2D_wien, color_space, patch_size, nb_threads, verbose) != EXIT_SUCCESS)
        return EXIT_FAILURE;

      set_progress_bar_data("BM3D denoising...", ((double) (i+1) / (double) numchunks));
    }

    // Reassemble chunks
    sub_divide(bgr_v, chunk_denoised, w_table, h_table, width, height, nchans, 32, FALSE);
    bgr_f = bgr_v.data();

      // Convert output from bgrbgr back to planar rgb and put back into fit
    if (fit->type == DATA_FLOAT) {
      if (nchans == 3) {
        bgrbgr_float_to_fits(fit, bgr_f, modulation);
      } else {
        for (unsigned i=0; i<npixels; i++) {
          fit->fdata[i] = (1.f-modulation) * fit->fdata[i] + modulation * bgr_f[i];
        }
      }
    } else {
      if (nchans == 3) {
        bgrbgr_float_to_word_fits(fit, bgr_f, modulation);
      } else {
         for (unsigned i=0; i<npixels; i++) {
          fit->data[i] = round_to_WORD(USHRT_MAX * ((1.f-modulation) * (fit->data[i] * invnorm) + modulation * bgr_f[i]));
        }
      }
    }
	return EXIT_SUCCESS;
}
