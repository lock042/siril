/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibImages.h"
#include "Utilities.h"
//#include "mt19937ar.h"
#include "algos/img_t/image_ops.hpp"  // imgops::pad_symmetric / unpad

#include <iostream>
#include <sstream>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

 using namespace std;

/**
 * @brief Load image, check the number of channels.
 *
 * @param p_name : name of the image to read;
 * @param o_im : vector which will contain the image : R, G and B concatenated;
 * @param o_imSize : will contain the size of the image;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_SUCCESS if the image has been loaded, EXIT_FAILURE otherwise.
 **/
/*int loadImage(
    char* p_name
,   std::vector<float> &o_im
,   ImageSize &o_imSize
,   const bool p_verbose
){
    //! read input image
    if (p_verbose) {
        cout << endl << "Read input image...";
    }
	float *imTmp = NULL;
	size_t w, h, c;
	imTmp = read_png_f32(p_name, &w, &h, &c);
	if (!imTmp) {
		cout << "error :: " << p_name << " not found or not a correct png image" << endl;
		return EXIT_FAILURE;
	}
	if (p_verbose) {
        cout << "done." << endl;
	}

	//! test if image is really a color image and exclude the alpha channel
	if (c > 2) {
	    unsigned k = 0;
	    while (k < w * h && imTmp[k] == imTmp[w * h + k] && imTmp[k] == imTmp[2 * w * h + k]) {
            k++;
	    }
        c = (k == w * h ? 1 : 3);
	}

	//! Some image informations
	if (p_verbose) {
        cout << "image size :" << endl;
        cout << " - width          = " << w << endl;
        cout << " - height         = " << h << endl;
        cout << " - nb of channels = " << c << endl;
	}

	//! Initializations
	o_imSize.width      = w;
	o_imSize.height     = h;
	o_imSize.nChannels  = c;
	o_imSize.wh         = w * h;
	o_imSize.whc        = w * h * c;
	o_im.resize(w * h * c);
	for (unsigned k = 0; k < w * h * c; k++)
        o_im[k] = imTmp[k];

    return EXIT_SUCCESS;
}
*/
/**
 * @brief write image.
 *
 * @param p_name : path+name+extension of the image;
 * @param i_im : vector which contains the image;
 * @param p_imSize : size of the image;
 * @param p_min, p_max : range of data (usually [0, 255]).
 *
 * @return EXIT_SUCCESS if the image has been saved, EXIT_FAILURE otherwise
 **/
/*int saveImage(
    char* p_name
,   std::vector<float> const& i_im
,   const ImageSize &p_imSize
,   const float p_min
,   const float p_max
){
    //! Allocate Memory
    float* imTmp = new float[p_imSize.whc];

    //! Check for boundary problems
    for (unsigned k = 0; k < p_imSize.whc; k++) {
        imTmp[k] = clip(i_im[k], p_min, p_max);
    }

    if (write_png_f32(p_name, imTmp, p_imSize.width, p_imSize.height, p_imSize.nChannels) != 0) {
        cout << "... failed to save png image " << p_name << endl;
        return EXIT_FAILURE;
    }

    //! Free Memory
    delete[] imTmp;

    return EXIT_SUCCESS;
}
*/
/**
 * @brief add noise to img.
 *
 * @param i_im : original noise-free image;
 * @param o_imNoisy = im + noise;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return none.
 **/
/*void addNoise(
    std::vector<float> const& i_im
,   std::vector<float> &o_imNoisy
,   const float p_sigma
,   const bool p_verbose
){
    if (p_verbose) {
        cout << "Add noise [sigma = " << p_sigma << "] ...";
    }

	//! Initialization
    o_imNoisy = i_im;
    mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) getpid());

    //! Add noise
    for (unsigned k = 0; k < i_im.size(); k++) {
        const double a = mt_genrand_res53();
        const double b = mt_genrand_res53();

        o_imNoisy[k] += p_sigma * (float) (sqrtl(-2.0l * log(a)) * cos(2.0l * M_PI * b));
    }

    if (p_verbose) {
        cout << "done." << endl;
    }
}
*/
/**
 * @brief Compute PSNR and RMSE between i_im1 and i_im2
 *
 * @param i_im1 : pointer to an allocated array of pixels;
 * @param i_im2 : pointer to an allocated array of pixels;
 * @param o_psnr  : will contain the PSNR;
 * @param o_rmse  : will contain the RMSE;
 * @param p_imageName: name of the image;
 * @param p_verbose: if true, print values of PSNR and RMSE.
 *
 * @return EXIT_FAILURE if both images haven't the same size.
 **/
int computePsnr(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   float &o_psnr
,   float &o_rmse
,   const char* p_imageName
,   const bool p_verbose
){
    if (i_im1.size() != i_im2.size()) {
        cout << "Can't compute PSNR & RMSE: images have different sizes: " << endl;
        cout << "i_im1 : " << i_im1.size() << endl;
        cout << "i_im2 : " << i_im2.size() << endl;
        return EXIT_FAILURE;
    }

    float sum = 0.f;
    for (unsigned k = 0; k < i_im1.size(); k++)
        sum += (i_im1[k] - i_im2[k]) * (i_im1[k] - i_im2[k]);

    o_rmse = sqrtf(sum / (float) i_im1.size());
    o_psnr = 20.f * log10f(255.f / o_rmse);

    if (p_verbose) {
        cout << p_imageName << endl;
        cout << "PSNR = " << o_psnr << endl;
        cout << "RMSE = " << o_rmse << endl;
    }

    return EXIT_SUCCESS;
}

/**
 * @brief Compute a difference image between i_im1 and i_im2.
 *
 * @param i_im1: reference image;
 * @param i_im2: image to compare;
 * @param o_imDiff: will contain the difference;
 * @param p_sigma : standard deviation of the noise;
 * @param p_min, p_max : range of data (usually [0, 255]);
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if i_im1 and i_im2 don't have the same size.
 **/
int computeDiff(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   std::vector<float> &o_imDiff
,   const float p_sigma
,   const float p_min
,   const float p_max
,   const bool p_verbose
){
    if (i_im1.size() != i_im2.size()) {
        cout << "Can't compute difference, i_im1 and i_im2 don't have the same size" << endl;
        cout << "i_im1 : " << i_im1.size() << endl;
        cout << "i_im2 : " << i_im2.size() << endl;
        return EXIT_FAILURE;
    }

    if (p_verbose) {
        cout << "Compute difference..." << endl;
    }

    const unsigned size = i_im1.size();
    if (o_imDiff.size() != size) {
        o_imDiff.resize(size);
    }

    for (unsigned k = 0; k < size; k++) {
        float value =  (i_im1[k] - i_im2[k] + p_sigma) * p_max / (2.f * p_sigma);
        o_imDiff[k] = clip(value, p_min, p_max);
    }

    if (p_verbose) {
        cout << "done." << endl;
    }

    return EXIT_SUCCESS;
}

/**
 * @brief Add boundary by symetry.
 *
 * @param i_im : image to symetrize;
 * @param o_imSym : will contain i_img with symetrized boundaries;
 * @param p_imSize : size of i_im;
 * @param p_imSizeSym : size of o_imSym.
 *
 * @return none.
 **/
int addBoundary(
	img_t<float> const& i_im
,	img_t<float> &o_imSym
,	const unsigned p_widthSym
,	const unsigned p_heightSym
){
	//! Parameters declarations
	const unsigned width	= i_im.w;
	const unsigned height	= i_im.h;
	const unsigned chnls	= i_im.d;
	const unsigned h		= p_heightSym;
	const unsigned w		= p_widthSym;

	if (w < width || h < height) {
		cout << "o_imSym must be greater than i_im!!!" << endl;
		return EXIT_FAILURE;
	}

	if ((unsigned) o_imSym.size != chnls * h * w) {
		o_imSym.resize(w, h, chnls);
	}

	//! Declaration
	for (unsigned c = 0; c < chnls; c++) {
		const unsigned dc1 = c * width * height;
		const unsigned dc2 = c * w * h;

		//! Center of the image
		for (unsigned i = 0; i < height; i++) {
			for (unsigned j = 0; j < width; j++) {
				o_imSym[dc2 + i * w + j] = i_im[dc1 + i * width + j];
			}
		}

		//! Right
		for (unsigned i = 0; i < height; i++) {
			for (unsigned j = width; j < w; j++) {
				o_imSym[dc2 + i * w + j] = o_imSym[dc2 + i * w + 2 * width - j - 1];
			}
		}

		//! Bottom
		for (unsigned i = height; i < h; i++) {
			for (unsigned j = 0; j < w; j++) {
				o_imSym[dc2 + i * w + j] = o_imSym[dc2 + (2 * height - i - 1) * w + j];
			}
		}
	}

	return EXIT_SUCCESS;
}

/**
 * @brief Remove boundaries added with addBoundary
 *
 * @param o_im : will contain the inner image;
 * @param i_imSym : contains i_im with symetrized boundaries;
 * @param p_imSize: size of o_im;
 * @param p_imSizeSym : size of i_imSym.
 *
 * @return none.
 **/
int removeBoundary(
	img_t<float> &o_im
,	img_t<float> const& i_imSym
,	const unsigned p_width
,	const unsigned p_height
){
	//! Parameters declaration
	const unsigned width	= p_width;
	const unsigned height	= p_height;
	const unsigned chnls	= i_imSym.d;
	const unsigned h		= i_imSym.h;
	const unsigned w		= i_imSym.w;

	if (w < width || h < height) {
		cout << "i_imSym must be greater than o_im!!!" << endl;
		return EXIT_FAILURE;
	}

	if ((unsigned) o_im.size != chnls * height * width) {
		o_im.resize(width, height, chnls);
	}

	for (unsigned c = 0, k = 0; c < chnls; c++) {
		const unsigned dc = c * w * h;

		for (unsigned i = 0; i < height; i++) {
			for (unsigned j = 0; j < width; j++, k++) {
				o_im[k] = i_imSym[dc + i * w + j];
			}
		}
	}

	return EXIT_SUCCESS;
}

/**
 * @brief Add boundaries by symetry
 *
 * @param i_im1: if p_isForward, contains the original image, otherwise contains
 *      the symetrized image;
 * @param o_im2: if p_isForward, will contain i_im1 symetrized, otherwise will
 *      contain the inner image of i_im1;
 * @param p_imSize : if p_isForward, size of i_im1, otherwise size of o_im2;
 * @param p_borderSize : size of the boundary;
 * @param p_isForward: if true, build io_imSym, otherwise build io_im.
 *
 * @return none.
 **/
void symetrizeImage(
	img_t<float> const& i_im1
,	img_t<float> &o_im2
,	const unsigned p_borderSize
,	const bool p_isForward
){
	// Delegates to the shared half-sample symmetric padding. nlbayes stores
	// images planar ([c*wh + row*w + col]), which is exactly img_t's planar
	// layout, so the img_t views map faithfully.
	const int border = static_cast<int>(p_borderSize);

	if (p_isForward) {
		// input is i_im1; output grows by border on every side
		o_im2 = imgops::pad_symmetric(i_im1, border);
	} else {
		// input is the padded image; output is the inner image
		o_im2 = imgops::unpad(i_im1, border);
	}
}

/**
 * @brief Transform the color space of an image, from RGB to YUV, or vice-versa.
 *
 * @param io_im: image on which the transform will be applied;
 * @param p_imSize: size of io_im;
 * @param p_isForward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
 *
 * @return none.
 **/
void transformColorSpace(
	img_t<float> &io_im
,	const bool p_isForward
){
	//! If the image has only one channel, do nothing
	if (io_im.d == 1) {
		return;
	}

	// Delegates to the shared orthonormal opponent transform. nlbayes' planar
	// [c*wh + ...] layout is img_t's planar layout;
	// imgops::color_transform replicates the original arithmetic bit-for-bit.
	imgops::color_transform(io_im, p_isForward);
}

/**
 * @brief Subdivide an image into small sub-images
 *
 * @param i_im : image to subdivide;
 * @param o_imSub : will contain all sub-images;
 * @param p_imSize : size of i_im;
 * @param p_imSizeSub : size of sub-images;
 * @param p_N : boundary around sub-images;
 * @param p_nb : number of sub-images wanted. Need to be a power of 2.
 *
 * @return EXIT_FAILURE in case of problems.
 **/
int subDivide(
	img_t<float> const& i_im
,	std::vector<img_t<float> > &o_imSub
,	const unsigned p_N
,	const unsigned p_nb
){
	const unsigned chnls = i_im.d;

	//! Determine width and height composition
	unsigned nW, nH;
	determineFactor(p_nb, nW, nH);
	const unsigned wTmp = ceil(float(i_im.w)  / float(nW));
	const unsigned hTmp = ceil(float(i_im.h) / float(nH));

	//! Add boundaries left and bottom
	img_t<float> imTmp;
	const unsigned widthTmp  = nW * wTmp;
	const unsigned heightTmp = nH * hTmp;
	if (addBoundary(i_im, imTmp, widthTmp, heightTmp) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
	}

    //! Symetrize boundaries of image
	img_t<float> imSymTmp;
	const unsigned widthSym  = widthTmp  + 2 * p_N;
	const unsigned heightSym = heightTmp + 2 * p_N;
	const unsigned whSym     = widthSym * heightSym;
    symetrizeImage(imTmp, imSymTmp, p_N, true);

    //! Obtain sub-images
    const unsigned widthSub  = wTmp + 2 * p_N;
    const unsigned heightSub = hTmp + 2 * p_N;
	o_imSub.resize(p_nb);
	for (unsigned p = 0, n = 0; p < nH; p++) {
        for (unsigned q = 0; q < nW; q++, n++) {
            o_imSub[n].resize(widthSub, heightSub, chnls);
            for (unsigned c = 0, k = 0; c < chnls; c++) {
                const unsigned dc = c * whSym + p * hTmp * widthSym + wTmp * q;
                for (unsigned i = 0; i < heightSub; i++) {
                    for (unsigned j = 0; j < widthSub; j++, k++) {
                        o_imSub[n][k] = imSymTmp[dc + i * widthSym + j];
                    }
                }
            }
        }
	}

	return EXIT_SUCCESS;
}

/**
 * @brief Reconstruct an image from its small sub-images
 *
 * @param o_im : image to reconstruct;
 * @param i_imSub : will contain all sub-images;
 * @param p_imSize : size of o_im;
 * @param p_imSizeSub : size of sub-images;
 * @param p_N : boundary around sub-images.
 *
 * @return EXIT_FAILURE in case of problems.
 **/
int subBuild(
	img_t<float> &o_im
,	std::vector<img_t<float> > const& i_imSub
,	const unsigned p_N
){
	const unsigned chnls    = o_im.d;
	const unsigned width    = o_im.w;
	const unsigned height   = o_im.h;
	const unsigned widthSub  = i_imSub[0].w;
	const unsigned heightSub = i_imSub[0].h;
	const unsigned whSub     = widthSub * heightSub;

    //! Determine width and height composition
	unsigned nW, nH;
	determineFactor(i_imSub.size(), nW, nH);
	const unsigned hTmp = heightSub - 2 * p_N;
	const unsigned wTmp = widthSub  - 2 * p_N;

	//! Obtain inner image (containing boundaries)
	const unsigned widthTmp  = wTmp * nW;
	const unsigned heightTmp = hTmp * nH;
	const unsigned whTmp     = widthTmp * heightTmp;
	img_t<float> imTmp(widthTmp, heightTmp, chnls);

	for (unsigned p = 0, n = 0; p < nH; p++) {
        for (unsigned q = 0; q < nW; q++, n++) {
            for (unsigned c = 0; c < chnls; c++) {
                const unsigned dc   = c * whTmp + p * hTmp * widthTmp + q * wTmp;
                const unsigned dcS  = c * whSub + p_N * widthSub + p_N;
                for (unsigned i = 0; i < hTmp; i++) {
                    for (unsigned j = 0; j < wTmp; j++) {
                        imTmp[dc + i * widthTmp + j] =
                            i_imSub[n][dcS + i * widthSub + j];
                    }
                }
            }
        }
	}

	//! Remove Boundaries
	return removeBoundary(o_im, imTmp, width, height);
}

