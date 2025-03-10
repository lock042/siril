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
 * @file NlBayes.cpp
 * @brief NL-Bayes denoising functions
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <math.h>

#include "NlBayes.h"
#include "LibMatrix.h"
#include "LibImages.h"
#include "Utilities.h"
#include "algos/anscombe.h"
#include "core/processing.h"
#include "gui/progress_and_log.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace NlBayes {

using namespace std;

/**
 * @brief Initialize Parameters of the NL-Bayes algorithm.
 *
 * @param o_paramStep1 : will contain the nlbParams for the first step of the algorithm;
 * @param o_paramStep2 : will contain the nlbParams for the second step of the algorithm;
 * @param p_sigma : standard deviation of the noise;
 * @param p_imSize: size of the image;
 * @param p_useArea1 : if true, use the homogeneous area trick for the first step;
 * @param p_useArea2 : if true, use the homogeneous area trick for the second step;
 * @param p_verbose : if true, print some informations.
 *
 * @return none.
 **/
void initializeNlbParameters(
	nlbParams &o_paramStep1
,	nlbParams &o_paramStep2
,   const float p_sigma
,	const ImageSize &p_imSize
,	const bool p_useArea1
,   const bool p_useArea2
,   const bool p_verbose
){
	//! Standard deviation of the noise
	o_paramStep1.sigma = p_sigma;
	o_paramStep2.sigma = p_sigma;

	//! Size of patches
	if (p_imSize.nChannels == 1) {
		o_paramStep1.sizePatch = (p_sigma < 30.f ? 5 : 7);
		o_paramStep2.sizePatch = 5;
	}
	else {
		o_paramStep1.sizePatch = (p_sigma < 20.f ? 3 :
                                 (p_sigma < 50.f ? 5 : 7));
        o_paramStep2.sizePatch = (p_sigma < 50.f ? 3 :
                                 (p_sigma < 70.f ? 5 : 7));
	}

	//! Number of similar patches
	if (p_imSize.nChannels == 1) {
		o_paramStep1.nSimilarPatches =	(p_sigma < 10.f ? 35 :
										(p_sigma < 30.f ? 45 :
										(p_sigma < 80.f ? 90 : 100)));
		o_paramStep2.nSimilarPatches =	(p_sigma < 20.f ? 15 :
										(p_sigma < 40.f ? 25 :
										(p_sigma < 80.f ? 30 : 45)));
	}
	else {
		o_paramStep1.nSimilarPatches = o_paramStep1.sizePatch * o_paramStep1.sizePatch * 3;
		o_paramStep2.nSimilarPatches = o_paramStep2.sizePatch * o_paramStep2.sizePatch * 3;
	}

	//! Offset: step between two similar patches
	o_paramStep1.offSet = o_paramStep1.sizePatch / 2;
	o_paramStep2.offSet = o_paramStep2.sizePatch / 2;

	//! Use the homogeneous area detection trick
	o_paramStep1.useHomogeneousArea = p_useArea1;
	o_paramStep2.useHomogeneousArea = p_useArea2;

	//! Size of the search window around the reference patch (must be odd)
	o_paramStep1.sizeSearchWindow = o_paramStep1.nSimilarPatches / 2;
	if (o_paramStep1.sizeSearchWindow % 2 == 0) {
		o_paramStep1.sizeSearchWindow++;
	}
	o_paramStep2.sizeSearchWindow = o_paramStep2.nSimilarPatches / 2;
	if (o_paramStep2.sizeSearchWindow % 2 == 0) {
		o_paramStep2.sizeSearchWindow++;
	}

	//! Size of boundaries used during the sub division
	o_paramStep1.boundary = int(1.5f * float(o_paramStep1.sizeSearchWindow));
	o_paramStep2.boundary = int(1.5f * float(o_paramStep2.sizeSearchWindow));

	//! Parameter used to determine if an area is homogeneous
	o_paramStep1.gamma = 1.05f;
	o_paramStep2.gamma = 1.05f;

	//! Parameter used to estimate the covariance matrix
	if (p_imSize.nChannels == 1) {
		o_paramStep1.beta = (p_sigma < 15.f ? 1.1f :
                            (p_sigma < 70.f ? 1.f : 0.9f));
		o_paramStep2.beta = (p_sigma < 15.f ? 1.1f :
                            (p_sigma < 35.f ? 1.f : 0.9f));
	}
	else {
		o_paramStep1.beta = 1.f;
		o_paramStep2.beta = (p_sigma < 50.f ? 1.2f : 1.f);
	}

	//! Parameter used to determine similar patches
	o_paramStep2.tau = 16.f * o_paramStep2.sizePatch * o_paramStep2.sizePatch * p_imSize.nChannels;

	//! Print information?
	o_paramStep1.verbose = p_verbose;
	o_paramStep2.verbose = p_verbose;

	//! Is first step?
	o_paramStep1.isFirstStep = true;
	o_paramStep2.isFirstStep = false;

	//! Boost the paste trick
	o_paramStep1.doPasteBoost = true;
	o_paramStep2.doPasteBoost = true;
}

/**
 * @brief Sanitizes bad data (inf, nan and data <0.0 or >1.0
 *
 * @param input: contains the image data to be sanitized;
 * @param ImageSize: contains dimensional information.
 *
 * @return sanitized image data.
 **/
vector<float> sanitize(vector<float> input, const ImageSize size) {
#define EPSILON 1.e-29
	size_t error = 0;
	for (size_t i = 0 ; i < input.size() ; i++) {
		size_t x = i % size.width;
		if (input[i] != input[i]) {
			error++;
		}
		if (isinf(input[i])) {
			error++;
		}
		if (error) {
			if (i == 0)
				input[i] = EPSILON;
			if (i > 0 && x < size.width - 1)
				input[i] = input[i-1];
			else if (x == 0)
				input[i] = input[i-size.width];
		}
		if (input[i] < 0.0) {
			input[i] = fmaxf(EPSILON, input[i]);
		}
	}
	siril_debug_print("\n");
	return input;
}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * @param i_imNoisy: contains the noisy image;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_imSize: size of the image;
 * @param p_useArea1 : if true, use the homogeneous area trick for the first step;
 * @param p_useArea2 : if true, use the homogeneous area trick for the second step;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if something wrong happens during the whole process.
 **/

int runNlBayes(
	std::vector<float> const& i_imNoisy
,   std::vector<float> &o_imBasic
,	std::vector<float> &o_imFinal
,	const ImageSize &p_imSize
,	const bool p_useArea1
,	const bool p_useArea2
,	const float p_sigma
,   const bool p_verbose
){
	//! Only 1, 3 or 4-channels images can be processed.
	const unsigned chnls = p_imSize.nChannels;
	if (! (chnls == 1 || chnls == 3 || chnls == 4)) {
		cout << "Wrong number of channels. Must be 1 or 3!!" << endl;
		return EXIT_FAILURE;
	}

	//! Number of available cores
	unsigned nbThreads = 1;
#ifdef _OPENMP
//    nbThreads = omp_get_max_threads();
    nbThreads = com.max_thread;
    if (p_verbose) {
        cout << "Open MP is used" << endl;
    }
#endif

	//! Initialization
	o_imBasic.resize(i_imNoisy.size());
	o_imFinal.resize(i_imNoisy.size());

	//! Parameters Initialization
	nlbParams paramStep1, paramStep2;
	initializeNlbParameters(paramStep1, paramStep2, p_sigma, p_imSize, p_useArea1, p_useArea2,
                         p_verbose);

	//! Step 1
	if (paramStep1.verbose) {
		cout << "1st Step...";
	}

	//! RGB to YUV
	vector<float> imNoisy = sanitize(i_imNoisy, p_imSize);
	transformColorSpace(imNoisy, p_imSize, true);

	//! Divide the noisy image into sub-images in order to easier parallelize the process
	const unsigned nbParts = 2 * nbThreads;
	vector<vector<float> > imNoisySub(nbParts), imBasicSub(nbParts), imFinalSub(nbParts);
	ImageSize imSizeSub;
	if (subDivide(imNoisy, imNoisySub, p_imSize, imSizeSub, paramStep1.boundary, nbParts)
		!= EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}
	int retval = 0;
	int thread = 0;
	int pass = 0;
	bool passset = false;
	unsigned threadcomplete = 0;
	//! Process all sub-images
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, nbParts/nbThreads) \
            shared(imNoisySub, imBasicSub, imFinalSub, imSizeSub, retval, passset, threadcomplete) \
			private(thread) firstprivate (paramStep1)
#endif
	for (int n = 0; n < (int) nbParts; n++) {
#ifdef _OPENMP
#pragma omp atomic
		threadcomplete++;
		thread = omp_get_thread_num();
		if (!passset && threadcomplete > nbThreads) {
			pass = 1;
			passset = true;
		}
#endif
		retval += processNlBayes(imNoisySub[n], imBasicSub[n], imFinalSub[n], imSizeSub, paramStep1, thread, pass);
	}

	if (retval != 0)
		return EXIT_FAILURE;

	//! Get the basic estimate
	if (subBuild(o_imBasic, imBasicSub, p_imSize, imSizeSub, paramStep1.boundary)
		!= EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}
	fprintf(stdout,"ready for YUV2RGB\n");
	//! YUV to RGB
	transformColorSpace(o_imBasic, p_imSize, false);
	o_imBasic = sanitize(o_imBasic, p_imSize);

	if (paramStep1.verbose) {
		cout << "done." << endl;
	}

	//! 2nd Step
	if (paramStep2.verbose) {
		cout << "2nd Step...";
	}

	//! Divide the noisy and basic images into sub-images in order to easier parallelize the process
	if (subDivide(i_imNoisy, imNoisySub, p_imSize, imSizeSub, paramStep2.boundary, nbParts)
        != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}
	if (subDivide(o_imBasic, imBasicSub, p_imSize, imSizeSub, paramStep2.boundary, nbParts)
		!= EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}
	thread = 0;
	pass = 0;
	passset = false;
	threadcomplete = 0;
	//! Process all sub-images
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, nbParts/nbThreads) \
            shared(imNoisySub, imBasicSub, imFinalSub, threadcomplete, passset) \
			private(thread) firstprivate (paramStep2)
#endif
	for (int n = 0; n < (int) nbParts; n++) {
#ifdef _OPENMP
#pragma omp atomic
		threadcomplete++;
		thread = omp_get_thread_num();
		if (!passset && threadcomplete > nbThreads) {
			pass = 1;
			passset = true;
		}
#endif
		retval += processNlBayes(imNoisySub[n], imBasicSub[n], imFinalSub[n], imSizeSub, paramStep2, thread, pass);
	}

	if (retval != 0)
		return EXIT_FAILURE;

	//! Get the final result
	if (subBuild(o_imFinal, imFinalSub, p_imSize, imSizeSub, paramStep2.boundary)
		!= EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	if (paramStep2.verbose) {
		cout << "done." << endl << endl;
	}

	return EXIT_SUCCESS;
}

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy image;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_imSize: size of i_imNoisy;
 * @param p_params: see nlbParams.
 *
 * @return none.
 **/
int processNlBayes(
	std::vector<float> const& i_imNoisy
,	std::vector<float> &io_imBasic
,	std::vector<float> &o_imFinal
,	const ImageSize &p_imSize
,	nlbParams &p_params
,	int thread
,	int pass
){
	//! Parameters initialization
	const unsigned sW		= p_params.sizeSearchWindow;
	const unsigned sP		= p_params.sizePatch;
	const unsigned sP2		= sP * sP;
	const unsigned sPC		= sP2 * p_imSize.nChannels;
	const unsigned nSP		= p_params.nSimilarPatches;
	unsigned nInverseFailed	= 0;
	const float threshold	= p_params.sigma * p_params.sigma * p_params.gamma *
                                (p_params.isFirstStep ? p_imSize.nChannels : 1.f);

	//! Allocate Sizes
	if (p_params.isFirstStep) {
		io_imBasic.resize(p_imSize.whc);
	}
	o_imFinal.resize(p_imSize.whc);

	//! Used matrices during Bayes' estimate
	vector<vector<float> > group3d(p_imSize.nChannels, vector<float> (nSP * sP2));
	vector<float> group3dNoisy(sW * sW * sPC), group3dBasic(sW * sW * sPC);
	vector<unsigned> index(p_params.isFirstStep ? nSP : sW * sW);
	matParams mat;
	mat.group3dTranspose.resize(p_params.isFirstStep ? nSP * sP2 : sW * sW * sPC);
	mat.tmpMat          .resize(p_params.isFirstStep ? sP2 * sP2 : sPC * sPC);
	mat.baricenter      .resize(p_params.isFirstStep ? sP2 : sPC);
	mat.covMat          .resize(p_params.isFirstStep ? sP2 * sP2 : sPC * sPC);
	mat.covMatTmp       .resize(p_params.isFirstStep ? sP2 * sP2 : sPC * sPC);

	//! ponderation: weight sum per pixel
	vector<float> weight(i_imNoisy.size(), 0.f);

	//! Mask: non-already processed patches
	vector<bool> mask(p_imSize.wh, false);

	//! Only pixels of the center of the image must be processed (not the boundaries)
	for (unsigned i = sW; i < p_imSize.height - sW; i++) {
		for (unsigned j = sW; j < p_imSize.width - sW; j++) {
			mask[i * p_imSize.width + j] = true;
		}
	}
	unsigned lastupdate = 0;
	float offset = 0;
	for (unsigned ij = 0; ij < p_imSize.wh; ij += p_params.offSet) {
		if (!get_thread_run())
			return 1;
		//! Only non-seen patches are processed
		if (mask[ij]) {
			//! Search for similar patches around the reference one
			unsigned nSimP = p_params.nSimilarPatches;
			if (p_params.isFirstStep) {
				estimateSimilarPatchesStep1(i_imNoisy, group3d, index, ij, p_imSize, p_params);
			}
			else {
				nSimP = estimateSimilarPatchesStep2(i_imNoisy, io_imBasic, group3dNoisy,
					group3dBasic, index, ij, p_imSize, p_params);
				nSimP = (nSimP < 1 ? 1 : nSimP);
			}

			//! Initialization
			bool doBayesEstimate = true;

			//! If we use the homogeneous area trick
			if (p_params.useHomogeneousArea) {
				if (p_params.isFirstStep) {
					doBayesEstimate = !computeHomogeneousAreaStep1(group3d, sP, nSP,
						threshold, p_imSize);
				}
				else {
					doBayesEstimate = !computeHomogeneousAreaStep2(group3dNoisy, group3dBasic,
						sP, nSimP, threshold, p_imSize);
				}
			}

			//! Else, use Bayes' estimate
			if (doBayesEstimate) {
				if (p_params.isFirstStep) {
					computeBayesEstimateStep1(group3d, mat, nInverseFailed, p_params);
				}
				else {
					computeBayesEstimateStep2(group3dNoisy, group3dBasic, mat, nInverseFailed,
						p_imSize, p_params, nSimP);
				}
			}

			//! Aggregation
			if (p_params.isFirstStep) {
				computeAggregationStep1(io_imBasic, weight, mask, group3d, index, p_imSize,
					p_params);
			}
			else {
				computeAggregationStep2(o_imFinal, weight, mask, group3dBasic, index, p_imSize,
					p_params, nSimP);
			}
		}
		if (thread == 0) {
			if (ij - lastupdate > p_imSize.wh >> 4) {
				lastupdate = ij;
				float second;
				second = (p_params.isFirstStep ? 0.f : 2.f);
				offset = (pass + second) / 4.f;
				set_progress_bar_data("NL-Bayes denoising...", offset + ((double)ij / (4.f * p_imSize.wh)));
			}
		}
	}

	//! Weighted aggregation
	computeWeightedAggregation(i_imNoisy, io_imBasic, o_imFinal, weight, p_params, p_imSize);

	if (nInverseFailed > 0 && p_params.verbose) {
		cout << "nInverseFailed = " << nInverseFailed << endl;
	}
	return 0;
}

/**
 * @brief Estimate the best similar patches to a reference one.
 *
 * @param i_im: contains the noisy image on which distances are processed;
 * @param o_group3d: will contain values of similar patches;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_imSize: size of the image;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
void estimateSimilarPatchesStep1(
	std::vector<float> const& i_im
,	std::vector<std::vector<float> > &o_group3d
,	std::vector<unsigned> &o_index
,	const unsigned p_ij
,	const ImageSize &p_imSize
,	const nlbParams &p_params
){
	//! Initialization
	const unsigned sW		= p_params.sizeSearchWindow;
	const unsigned sP		= p_params.sizePatch;
	const unsigned width	= p_imSize.width;
	const unsigned chnls	= p_imSize.nChannels;
	const unsigned wh		= width * p_imSize.height;
	const unsigned ind		= p_ij - (sW - 1) * (width + 1) / 2;
	const unsigned nSimP	= p_params.nSimilarPatches;
	vector<pair<float, unsigned> > distance(sW * sW);

	//! Compute distance between patches
	for (unsigned i = 0; i < sW; i++) {
		for (unsigned j = 0; j < sW; j++) {
			const unsigned k = i * width + j + ind;
			float diff = 0.f;
			for (unsigned p = 0; p < sP; p++) {
				for (unsigned q = 0; q < sP; q++) {
					const float tmpValue = i_im[p_ij + p * width + q] - i_im[k + p * width + q];
					diff += tmpValue * tmpValue;
				}
			}

			//! Save all distances
			distance[i * sW + j] = make_pair(diff, k);
		}
	}

	//! Keep only the N2 best similar patches
	partial_sort(distance.begin(), distance.begin() + nSimP, distance.end(), comparaisonFirst);

	//! Register position of patches
	for (unsigned n = 0; n < nSimP; n++) {
		o_index[n] = distance[n].second;
	}

	//! Register similar patches into the 3D group
	for (unsigned c = 0; c < chnls; c++) {
		for (unsigned p = 0, k = 0; p < sP; p++) {
			for (unsigned q = 0; q < sP; q++) {
				for (unsigned n = 0; n < nSimP; n++, k++) {
					o_group3d[c][k] = i_im[o_index[n] + p * width + q + c * wh];
				}
			}
		}
	}
}

/**
 * @brief Keep from all near patches the similar ones to the reference patch for the second step.
 *
 * @param i_imNoisy: contains the original noisy image;
 * @param i_imBasic: contains the basic estimation;
 * @param o_group3dNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_group3dBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_imSize: size of images;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	std::vector<float> const& i_imNoisy
,	std::vector<float> const& i_imBasic
,	std::vector<float> &o_group3dNoisy
,	std::vector<float> &o_group3dBasic
,	std::vector<unsigned> &o_index
,	const unsigned p_ij
,	const ImageSize &p_imSize
,	const nlbParams &p_params
){
	//! Initialization
	const unsigned width	= p_imSize.width;
	const unsigned chnls	= p_imSize.nChannels;
	const unsigned wh		= width * p_imSize.height;
	const unsigned sP		= p_params.sizePatch;
	const unsigned sW		= p_params.sizeSearchWindow;
	const unsigned ind		= p_ij - (sW - 1) * (width + 1) / 2;
	vector<pair<float, unsigned> > distance(sW * sW);

	//! Compute distance between patches
	for (unsigned i = 0; i < sW; i++) {
		for (unsigned j = 0; j < sW; j++) {
			const unsigned k = i * width + j + ind;
			float diff = 0.0f;

			for (unsigned c = 0; c < chnls; c++) {
				const unsigned dc = c * wh;
				for (unsigned p = 0; p < sP; p++) {
					for (unsigned q = 0; q < sP; q++) {
						const float tmpValue = i_imBasic[dc + p_ij + p * width + q]
											- i_imBasic[dc + k + p * width + q];
						diff += tmpValue * tmpValue;
					}
				}
			}

			//! Save all distances
			distance[i * sW + j] = make_pair(diff, k);
		}
	}

	//! Keep only the nSimilarPatches best similar patches
	partial_sort(distance.begin(), distance.begin() + p_params.nSimilarPatches, distance.end(),
		comparaisonFirst);

	//! Save index of similar patches
	const float threshold = (p_params.tau > distance[p_params.nSimilarPatches - 1].first ?
							p_params.tau : distance[p_params.nSimilarPatches - 1].first);
	unsigned nSimP = 0;

	//! Register position of similar patches
	for (unsigned n = 0; n < distance.size(); n++) {
		if (distance[n].first < threshold) {
			o_index[nSimP++] = distance[n].second;
		}
	}

	//! Save similar patches into 3D groups
	for (unsigned c = 0, k = 0; c < chnls; c++) {
		for (unsigned p = 0; p < sP; p++) {
			for (unsigned q = 0; q < sP; q++) {
				for (unsigned n = 0; n < nSimP; n++, k++) {
					o_group3dNoisy[k] = i_imNoisy[c * wh + o_index[n] + p * width + q];
					o_group3dBasic[k] = i_imBasic[c * wh + o_index[n] + p * width + q];
				}
			}
		}
	}

	return nSimP;
}

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group3d: contains for each channels values of similar patches. If an homogeneous area
 *			is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_doLinearRegression: if true, apply a linear regression to average value of pixels;
 * @param p_imSize: size of the image.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep1(
	std::vector<std::vector<float> > &io_group3d
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const ImageSize &p_imSize
){
	//! Initialization
	const unsigned N = p_sP * p_sP * p_nSimP;

	//! Compute the standard deviation of the set of patches
	float stdDev = 0.f;
	for (unsigned c = 0; c < p_imSize.nChannels; c++) {
		stdDev += computeStdDeviation(io_group3d[c], p_sP * p_sP, p_nSimP, 1);
	}

	//! If we are in an homogeneous area
	if (stdDev < p_threshold) {
		for (unsigned c = 0; c < p_imSize.nChannels; c++) {
            float mean = 0.f;

            for (unsigned k = 0; k < N; k++) {
                mean += io_group3d[c][k];
            }

            mean /= (float) N;

            for (unsigned k = 0; k < N; k++) {
                io_group3d[c][k] = mean;
            }
        }
		return 1;
	}
	else {
		return 0;
	}
}

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group3dNoisy: contains values of similar patches for the noisy image;
 * @param io_group3dBasic: contains values of similar patches for the basic image. If an homogeneous
 *		area is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_imSize: size of the image.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep2(
	std::vector<float> const& i_group3dNoisy
,	std::vector<float> &io_group3dBasic
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const ImageSize &p_imSize
){
	//! Parameters
	const unsigned sP2	= p_sP * p_sP;
	const unsigned sPC = sP2 * p_imSize.nChannels;

	//! Compute the standard deviation of the set of patches
	const float stdDev = computeStdDeviation(i_group3dNoisy, sP2, p_nSimP, p_imSize.nChannels);

	//! If we are in an homogeneous area
	if (stdDev < p_threshold) {
		for (unsigned c = 0; c < p_imSize.nChannels; c++) {
            float mean = 0.f;

            for (unsigned n = 0; n < p_nSimP; n++) {
                for (unsigned k = 0; k < sP2; k++) {
                    mean += io_group3dBasic[n * sPC + c * sP2 + k];
                }
            }

            mean /= float(sP2 * p_nSimP);

            for (unsigned n = 0; n < p_nSimP; n++) {
                for (unsigned k = 0; k < sP2; k++) {
                    io_group3dBasic[n * sPC + c * sP2 + k] = mean;
                }
            }
		}
		return 1;
	}
	else {
		return 0;
	}
}

/**
 * @brief Compute the Bayes estimation.
 *
 * @param io_group3d: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *		- group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
 void computeBayesEstimateStep1(
	std::vector<std::vector<float> > &io_group3d
,	matParams &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams &p_params
){
	//! Parameters
	const unsigned chnls = io_group3d.size();
	const unsigned nSimP = p_params.nSimilarPatches;
	const unsigned sP2   = p_params.sizePatch * p_params.sizePatch;
	const float valDiag  = p_params.beta * p_params.sigma * p_params.sigma;

	//! Bayes estimate
	for (unsigned c = 0; c < chnls; c++) {

	    //! Center data around the baricenter
		centerData(io_group3d[c], i_mat.baricenter, nSimP, sP2);

		//! Compute the covariance matrix of the set of similar patches
		covarianceMatrix(io_group3d[c], i_mat.covMat, nSimP, sP2);

		//! Bayes' Filtering
		if (inverseMatrix(i_mat.covMat, sP2) == EXIT_SUCCESS) {
            productMatrix(i_mat.group3dTranspose, i_mat.covMat, io_group3d[c], sP2, sP2, nSimP);
            for (unsigned k = 0; k < sP2 * nSimP; k++) {
                io_group3d[c][k] -= valDiag * i_mat.group3dTranspose[k];
            }
		}
		else {
			io_nInverseFailed++;
		}

		//! Add baricenter
		for (unsigned j = 0, k = 0; j < sP2; j++) {
			for (unsigned i = 0; i < nSimP; i++, k++) {
			    io_group3d[c][k] += i_mat.baricenter[j];
			}
		}
	}
}

/**
 * @brief Compute the Bayes estimation.
 *
 * @param i_group3dNoisy: contains all similar patches in the noisy image;
 * @param io_group3dBasic: contains all similar patches in the basic image. Will contain estimates
 *			for all similar patches;
 * @param i_mat: contains :
 *		- group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_imSize: size of the image;
 * @param p_params: see processStep2 for more explanations;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeBayesEstimateStep2(
	std::vector<float> &i_group3dNoisy
,	std::vector<float> &io_group3dBasic
,	matParams &i_mat
,	unsigned &io_nInverseFailed
,	const ImageSize &p_imSize
,	nlbParams p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float diagVal = p_params.beta * p_params.sigma * p_params.sigma;
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch * p_imSize.nChannels;

	//! Center 3D groups around their baricenter
	centerData(io_group3dBasic, i_mat.baricenter, p_nSimP, sPC);
	centerData(i_group3dNoisy, i_mat.baricenter, p_nSimP, sPC);

	//! Compute the covariance matrix of the set of similar patches
	covarianceMatrix(io_group3dBasic, i_mat.covMat, p_nSimP, sPC);

	//! Bayes' Filtering
    for (unsigned k = 0; k < sPC; k++) {
        i_mat.covMat[k * sPC + k] += diagVal;
    }

	//! Compute the estimate
	if (inverseMatrix(i_mat.covMat, sPC) == EXIT_SUCCESS) {
        productMatrix(io_group3dBasic, i_mat.covMat, i_group3dNoisy, sPC, sPC, p_nSimP);
        for (unsigned k = 0; k < sPC * p_nSimP; k++) {
            io_group3dBasic[k] = i_group3dNoisy[k] - diagVal * io_group3dBasic[k];
        }
	}
	else {
		io_nInverseFailed++;
	}

	//! Add baricenter
	for (unsigned j = 0, k = 0; j < sPC; j++) {
		for (unsigned i = 0; i < p_nSimP; i++, k++) {
			io_group3dBasic[k] += i_mat.baricenter[j];
		}
	}
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
void computeAggregationStep1(
	std::vector<float> &io_im
,	std::vector<float> &io_weight
,	std::vector<bool> &io_mask
,	std::vector<std::vector<float> > const& i_group3d
,	std::vector<unsigned> const& i_index
,	const ImageSize &p_imSize
,	const nlbParams &p_params
){
	//! Parameters initializations
	const unsigned chnls	= p_imSize.nChannels;
	const unsigned width	= p_imSize.width;
	const unsigned height	= p_imSize.height;
	const unsigned sP		= p_params.sizePatch;
	const unsigned nSimP	= p_params.nSimilarPatches;

	//! Aggregate estimates
	for (unsigned n = 0; n < nSimP; n++) {
		const unsigned ind = i_index[n];
		for (unsigned c = 0; c < chnls; c++) {
			const unsigned ij = ind + c * width * height;
			for (unsigned p = 0; p < sP; p++) {
				for (unsigned q = 0; q < sP; q++) {
					io_im[ij + p * width + q] += i_group3d[c][(p * sP + q) * nSimP + n];
					io_weight[ij + p * width + q]++;
				}
			}
		}

		//! Use Paste Trick
		io_mask[ind] = false;

		if (p_params.doPasteBoost) {
			io_mask[ind - width ] = false;
			io_mask[ind + width ] = false;
			io_mask[ind - 1		] = false;
			io_mask[ind + 1		] = false;
		}
	}
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeAggregationStep2(
	std::vector<float> &io_im
,	std::vector<float> &io_weight
,	std::vector<bool> &io_mask
,	std::vector<float> const& i_group3d
,	std::vector<unsigned> const& i_index
,	const ImageSize &p_imSize
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls	= p_imSize.nChannels;
	const unsigned width	= p_imSize.width;
	const unsigned wh		= width * p_imSize.height;
	const unsigned sP		= p_params.sizePatch;

	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++) {
		const unsigned ind = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++) {
			const unsigned ij = ind + c * wh;
			for (unsigned p = 0; p < sP; p++) {
				for (unsigned q = 0; q < sP; q++, k++) {
					io_im[ij + p * width + q] += i_group3d[k * p_nSimP + n];
					io_weight[ij + p * width + q]++;
				}
			}
		}

		//! Apply Paste Trick
		io_mask[ind] = false;

		if (p_params.doPasteBoost) {
			io_mask[ind - width ] = false;
			io_mask[ind + width ] = false;
			io_mask[ind - 1     ] = false;
			io_mask[ind + 1     ] = false;
		}
	}
}

/**
 * @brief Compute the final weighted aggregation.
 *
 * i_imReference: image of reference, when the weight if null;
 * io_imResult: will contain the final image;
 * i_weight: associated weight for each estimate of pixels.
 *
 * @return : none.
 **/
void computeWeightedAggregation(
	std::vector<float> const& i_imNoisy
,	std::vector<float> &io_imBasic
,	std::vector<float> &io_imFinal
,	std::vector<float> const& i_weight
,	const nlbParams &p_params
,	const ImageSize &p_imSize
){
	for (unsigned c = 0, k = 0; c < p_imSize.nChannels; c++) {

		for (unsigned ij = 0; ij < p_imSize.wh; ij++, k++) {

			//! To avoid weighting problem (particularly near boundaries of the image)
			if (i_weight[k] > 0.f) {
				if (p_params.isFirstStep) {
					io_imBasic[k] /= i_weight[k];
				}
				else {
					io_imFinal[k] /= i_weight[k];
				}
			}
			else {
				if (p_params.isFirstStep) {
					io_imBasic[k] = i_imNoisy[k];
				}
				else {
					io_imFinal[k] = io_imBasic[k];
				}
			}
		}
	}
}
}
