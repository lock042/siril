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

/*
float* cutchunk (float *img, unsigned x0, unsigned y0, unsigned xlen, unsigned ylen, unsigned width, unsigned nchans) {
	// img stored as bgrbgrbgr
	float *chunk = (float*) calloc(xlen * ylen * nchans, sizeof(float));
	for (unsigned i = y0, ii = 0; i < (y0 + ylen); i++, ii++) {
		for (unsigned j = x0, jj = 0; j < (x0 + xlen); j++, jj++) {
			for (unsigned chan=0; chan < nchans ; chan++) {
				chunk[chan + (jj + ii * xlen) * nchans] = img[chan + (j + i * width) * nchans];
			}
		}
	}
	return chunk;
}

void pastechunk (float *img, float *chunk, unsigned x0, unsigned y0, unsigned xlen, unsigned ylen, unsigned pad, unsigned width, unsigned nchans) {
	// pad is the half-window size
	for (unsigned i = y0, ii = pad; i < (y0 + ylen); i++, ii++) {
		for (unsigned j = x0, jj = pad; j < (x0 + xlen); j++, jj++) {
			for (unsigned chan=0; chan < nchans; chan++) {
				img[chan + (j + i * width) * nchans] = chunk[chan + (jj + ii * xlen) * nchans];
			}
		}
	}
}

unsigned round_up(unsigned num, unsigned factor) {
	return (num - 1 - (num - 1) % factor + factor);
}
*/

void bgrbgr_float_to_fits(fits *image, float *bgrbgr) {
	size_t ndata = image->rx * image->ry * 3;
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		image->fpdata[BLAYER][j] = bgrbgr[i + 0];
		image->fpdata[GLAYER][j] = bgrbgr[i + 1];
		image->fpdata[RLAYER][j] = bgrbgr[i + 2];
	}
}

void bgrbgr_float_to_word_fits(fits *image, float *bgrbgr) {
	size_t ndata = image->rx * image->ry * 3;
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		image->pdata[BLAYER][j] = round_to_WORD(bgrbgr[i + 0] * USHRT_MAX);
		image->pdata[GLAYER][j] = round_to_WORD(bgrbgr[i + 1] * USHRT_MAX);
		image->pdata[RLAYER][j] = round_to_WORD(bgrbgr[i + 2] * USHRT_MAX);
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

extern "C" int do_bm3d(fits *fit) {
    // Parameters
  const unsigned width = fit->naxes[0];
    const unsigned height = fit->naxes[1];
    const unsigned nchans = fit->naxes[2];
    const unsigned npixels = width * height;
    float norm = 1.f;
    if (fit->type == DATA_USHORT)
      norm = USHRT_MAX_SINGLE;
    float invnorm = 1 / norm;

    float ratio = 0.5f; // This will become a user parameter
    if (nchans == 1)
      ratio *= 0.75f;


// Measure background noise
    float fSigma = 0.f;
    for (size_t chan = 0 ; chan < nchans ; chan++) {
      imstats* stat = statistics(NULL, -1, fit, chan, NULL, STATS_SIGMEAN, MULTI_THREADED);
      fSigma += stat->bgnoise / norm;
      free_stats(stat);
    }
    fSigma /= nchans;
    siril_log_message(_("Auto parametrisation: measured background noise level is %f\n"),fSigma);
    fSigma *= ratio; // Replace 1.0f with a user-specified parameter. Reasonable default seems to be
                     // 1.f for colour images and 0.75f for mono.

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
    unsigned numchunks = (unsigned) imgmemMpix / (memGB / 5.0); // Be smarter about this considering available memory and size of image
    if (numchunks < 1)
      numchunks = 1;
    fprintf(stdout, "memory: %f GB, imgsize: %f Mpix, numchunks: %u\n", memGB, imgmemMpix, numchunks);
 	// Cut into manageable chunks to avoid OOM
    vector<vector<float> > chunk_noisy(numchunks);
    vector<vector<float> > chunk_basic(numchunks);
    vector<vector<float> > chunk_denoised(numchunks);
    vector<unsigned> w_table(numchunks);
    vector<unsigned> h_table(numchunks);
    sub_divide(bgr_v, chunk_noisy, w_table, h_table, width, height, nchans, 32, TRUE); // 32 is nWien * 2;

    // Run bm3d on each chunk in turn.
    for (unsigned i = 0; i < numchunks ; i++) {
      if (run_bm3d(fSigma, chunk_noisy[i], chunk_basic[i], chunk_denoised[i],
         w_table[i], h_table[i], nchans, useSD_1, useSD_2, tau_2D_hard, tau_2D_wien, color_space, patch_size, nb_threads, verbose) != EXIT_SUCCESS)
      return EXIT_FAILURE;

      set_progress_bar_data("BM3D denoising...", (i / numchunks));
    }

	// Reassemble chunks
    sub_divide(bgr_v, chunk_denoised, w_table, h_table, width, height, nchans, 32, FALSE);
    bgr_f = bgr_v.data();

	// Convert output from bgrbgr back to planar rgb and put back into fit
    if (fit->type == DATA_FLOAT) {
      if (nchans == 3) {
        bgrbgr_float_to_fits(fit, bgr_f);
      } else {
        for (unsigned i=0; i<npixels; i++) {
          fit->fdata[i] = bgr_f[i];
        }
      }
    } else {
      if (nchans == 3) {
        bgrbgr_float_to_word_fits(fit, bgr_f);
      } else {
         for (unsigned i=0; i<npixels; i++) {
          fit->data[i] = round_to_WORD(bgr_f[i] * USHRT_MAX);
        }
      }
    }
	return EXIT_SUCCESS;
}


/*
int fits_apply_bm3d(fits *fit) {
    // Parameters
	float fSigma = 2; // Replace this with getting the nosie level from statistics
    const bool useSD_1 = FALSE;
    const bool useSD_2 = FALSE;
    const unsigned tau_2D_hard = DCT;
    const unsigned tau_2D_wien = BIOR;
    const unsigned color_space = OPP;
    const int patch_size = 0;
    const int nb_threads = 0;
    const bool verbose = FALSE;

	//! Declarations
	vector<float> strip_basic, strip_denoised;
    const unsigned width = fit->naxes[0], height = fit->naxes[1], chnls = fit->naxes[2];
    const unsigned npixels = width * height;
    const unsigned num_strips = min(npixels / 1e6, 1);
    const unsigned ndata = npixels * chnls;
    const unsigned strip_height = height % num_strips;
    // Ensure strip height is at least big enough for the minimum dimension for the denoiser
    const unsigned strip_npixels = strip_height * width;
    const unsigned strip_ndata = strip_height * width * chnls;
    const unsigned last_strip_height = height - (num_strips - 1) * strip_height;
    // Ensure last_strip_height is at least big enough for the minimum dimension for the denoiser
    const unsigned last_strip_npixels = last_strip_height * width;
    const unsigned last_strip_ndata = last_strip_npixels * chnls;

    //! Load image
    float *bgr_f = fits_to_bgrbgr_float(fit);

   //! Denoising
   // Split image into strips
    float *buf_strip[num_strips];
    float *denoised_strip;
    for (unsigned strip = 0; strip < num_strips - 1 ; strip ++) { // All except the last strip
      buf_strip[strip] = bgr_f + (strip_ndata * strip);
      vector<float> strip_noisy {buf_strip[strip] , buf_strip[strip] + strip_ndata};
      if (run_bm3d(fSigma, strip_noisy, strip_basic, strip_denoised,
            width, strip_height, chnls, useSD_1, useSD_2, tau_2D_hard, tau_2D_wien, color_space, patch_size, nb_threads, verbose)
            != EXIT_SUCCESS)
        return EXIT_FAILURE;
      // Put the denoised strips back into the fits
      denoised_strip = strip_denoised.data();
      for (unsigned i=0; i<strip_npixels;i++) {
        fit->fpdata[RLAYER][(strip * strip_npixels) + i] = denoised_strip[3*i];
        fit->fpdata[GLAYER][(strip * strip_npixels) + i] = denoised_strip[3*i+1];
        fit->fpdata[BLAYER][(strip * strip_npixels) + i] = denoised_strip[3*i+2];
      }
      // Add a progress bar update here
    }
    // Do the last strip
    unsigned laststrip = num_strips - 1;
    buf_strip[num_strips - 1] = bgr_f + (strip_ndata * (laststrip));
    vector<float> strip_noisy {buf_strip[laststrip] , buf_strip[laststrip] + last_strip_ndata};
     if (run_bm3d(fSigma, strip_noisy, strip_basic, strip_denoised,
          width, strip_height, chnls, useSD_1, useSD_2, tau_2D_hard, tau_2D_wien, color_space, patch_size, nb_threads, verbose)
          != EXIT_SUCCESS)
      return EXIT_FAILURE;
    denoised_strip = strip_denoised.data();
    for (unsigned i=0; i<last_strip_npixels;i++) {
      fit->fpdata[RLAYER][(laststrip * strip_npixels) + i] = denoised_strip[3*i];
      fit->fpdata[GLAYER][(laststrip * strip_npixels) + i] = denoised_strip[3*i+1];
      fit->fpdata[BLAYER][(laststrip * strip_npixels) + i] = denoised_strip[3*i+2];
    }
    free(bgr_f);
/*   //! save noisy, denoised and differences images
   cout << endl << "Save images...";

   if (argc > 4)
   if (save_image(argv[4], img_basic, width, height, chnls) != EXIT_SUCCESS)
      return EXIT_FAILURE;

	if (save_image(argv[3], img_denoised, width, height, chnls) != EXIT_SUCCESS)
		return EXIT_FAILURE;

    cout << "done." << endl;

	return EXIT_SUCCESS;
}
*/
