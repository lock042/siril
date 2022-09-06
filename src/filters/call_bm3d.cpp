#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>

#include "bm3d.h"
#include "bm3d_utilities.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3
#define DCT       4
#define BIOR      5
#define HADAMARD  6
#define NONE      7

using namespace std;

// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
const char *pick_option(int *c, char **v, const char *o, const char *d) {
  int id = d ? 1 : 0;
  for (int i = 0; i < *c - id; i++) {
    if (v[i][0] == '-' && 0 == strcmp(v[i] + 1, o)) {
      char *r = v[i + id] + 1 - id;
      for (int j = i; j < *c - id; j++)
        v[j] = v[j + id + 1];
      *c -= id + 1;
      return r;
    }
  }
  return d;
}

/**
 * @file   main.cpp
 * @brief  Main executable file. Do not use lib_fftw to
 *         process DCT.
 *
 * @author MARC LEBRUN  <marc.lebrun@cmla.ens-cachan.fr>
 */
/*int main(int argc, char **argv)
{
  //! Variables initialization
  const char *_tau_2D_hard = pick_option(&argc, argv, "tau_2d_hard", "bior");
  const char *_tau_2D_wien = pick_option(&argc, argv, "tau_2d_wien", "dct");
  const char *_color_space = pick_option(&argc, argv, "color_space", "opp");
  const char *_patch_size = pick_option(&argc, argv, "patch_size", "0"); // >0: overrides default
  const char *_nb_threads = pick_option(&argc, argv, "nb_threads", "0");
  const bool useSD_1 = pick_option(&argc, argv, "useSD_hard", NULL) != NULL;
  const bool useSD_2 = pick_option(&argc, argv, "useSD_wien", NULL) != NULL;
  const bool verbose = pick_option(&argc, argv, "verbose", NULL) != NULL;

	//! Check parameters
	const unsigned tau_2D_hard  = (strcmp(_tau_2D_hard, "dct" ) == 0 ? DCT :
                                 (strcmp(_tau_2D_hard, "bior") == 0 ? BIOR : NONE));
    if (tau_2D_hard == NONE) {
        cout << "tau_2d_hard is not known." << endl;
        argc = 0; //abort
    }
	const unsigned tau_2D_wien  = (strcmp(_tau_2D_wien, "dct" ) == 0 ? DCT :
                                 (strcmp(_tau_2D_wien, "bior") == 0 ? BIOR : NONE));
    if (tau_2D_wien == NONE) {
        cout << "tau_2d_wien is not known." << endl;
        argc = 0; //abort
    };
	const unsigned color_space  = (strcmp(_color_space, "rgb"  ) == 0 ? RGB   :
                                 (strcmp(_color_space, "yuv"  ) == 0 ? YUV   :
                                 (strcmp(_color_space, "ycbcr") == 0 ? YCBCR :
                                 (strcmp(_color_space, "opp"  ) == 0 ? OPP   : NONE))));
    if (color_space == NONE) {
        cout << "color_space is not known." << endl;
        argc = 0; //abort
    };

  const int patch_size = atoi(_patch_size);
    if (patch_size < 0)
    {
      cout << "The patch_size parameter must not be negative." << endl;
      return EXIT_FAILURE;
    } else {
      const unsigned patch_size = (unsigned) patch_size;
    }
  const int nb_threads = atoi(_nb_threads);
    if (nb_threads < 0)
    {
      cout << "The nb_threads parameter must not be negative." << endl;
      return EXIT_FAILURE;
    } else {
      const unsigned nb_threads = (unsigned) nb_threads;
    }

  //! Check if there is the right call for the algorithm
  if (argc < 4) {
    cerr << "usage: " << argv[0] << " input sigma output [basic]\n\
             [-tau_2d_hard {dct,bior} (default: bior)]\n\
             [-useSD_hard]\n\
             [-tau_2d_wien {dct,bior} (default: dct)]\n\
             [-useSD_wien]\n\
             [-color_space {rgb,yuv,opp,ycbcr} (default: opp)]\n\
             [-patch_size {0,8,...} (default: 0, auto size, 8 or 12 depending on sigma)]\n\
             [-nb_threads (default: 0, auto number)]\n\
             [-verbose]" << endl;
    return EXIT_FAILURE;
  }
*/
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

int do_bm3d(fits *fit) {
    // Parameters
    const unsigned width = fit->naxes[0];
    const unsigned height = fit->naxes[1];
    const unsigned nchans = fit->naxes[2];
	float fSigma = 2; // Replace this with getting the noise level from statistics

	// The algorithm parameters set below are the optimal parameters discussed in
	// (Lebrun, 2012): http://www.ipol.im/pub/art/2012/l-bm3d/
	const bool useSD_1 = FALSE;
    const bool useSD_2 = FALSE;
    const unsigned tau_2D_hard = BIOR;
    const unsigned tau_2D_wien = DCT;
    const unsigned color_space = OPP;
    const int patch_size = 0;
    const int nb_threads = 0;
    const bool verbose = FALSE;

    float *bgr_f = fits_to_bgrbgr_float(fit);
    vector<float> bgr_v { bgr_f, bgr_f + width * height * nchans * sizeof(float) };

    // Cut into manageable chunks to avoid OOM
	unsigned numchunks = 16; // Be smarter about this considering available memory and size of image
    vector<vector<float> > chunk(numchunks);
    vector<vector<float> > chunk_basic(numchunks);
    vector<vector<float> > chunk_denoised(numchunks);
    vector<unsigned> w_table(numchunks);
    vector<unsigned> h_table(numchunks);
    sub_divide(bgr_v, chunk, w_table, h_table, width, height, nchans, 32, FALSE); // 32 is nWien * 2;

	for (unsigned i = 0; i < numchunks ; i++) {
      if (run_bm3d(fSigma, chunk[i], chunk_basic[i], chunk_denoised[i],
         w_table[i], h_table[i], nchans, useSD_1, useSD_2, tau_2D_hard, tau_2D_wien, color_space, patch_size, nb_threads, verbose) != EXIT_SUCCESS)
      return EXIT_FAILURE;
	}

	// Reassemble chunks
    sub_divide(bgr_v, chunk_denoised, w_table, h_table, width, height, nchans, 32, TRUE);
    bgr_f = bgr_v.data();

	// Convert output from bgrbgr back to planar rgb and put back into fit
    bgrbgr_float_to_fits(fit, bgr_f);

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
