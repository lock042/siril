/*M///////////////////////////////////////////////////////////////////////////////////////
 //
 //  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
 //
 //  By downloading, copying, installing or using the software you agree to this license.
 //  If you do not agree to this license, do not download, install,
 //  copy or use the software.
 //
 //
 //                        Intel License Agreement
 //                For Open Source Computer Vision Library
 //
 // Copyright (C) 2000, Intel Corporation, all rights reserved.
 // Third party copyrights are property of their respective owners.
 //
 // Redistribution and use in source and binary forms, with or without modification,
 // are permitted provided that the following conditions are met:
 //
 //   * Redistribution's of source code must retain the above copyright notice,
 //     this list of conditions and the following disclaimer.
 //
 //   * Redistribution's in binary form must reproduce the above copyright notice,
 //     this list of conditions and the following disclaimer in the documentation
 //     and/or other materials provided with the distribution.
 //
 //   * The name of Intel Corporation may not be used to endorse or promote products
 //     derived from this software without specific prior written permission.
 //
 // This software is provided by the copyright holders and contributors "as is" and
 // any express or implied warranties, including, but not limited to, the implied
 // warranties of merchantability and fitness for a particular purpose are disclaimed.
 // In no event shall the Intel Corporation or contributors be liable for any direct,
 // indirect, incidental, special, exemplary, or consequential damages
 // (including, but not limited to, procurement of substitute goods or services;
 // loss of use, data, or profits; or business interruption) however caused
 // and on any theory of liability, whether in contract, strict liability,
 // or tort (including negligence or otherwise) arising in any way out of
 // the use of this software, even if advised of the possibility of such damage.
 //
 //M*/

/* the start of this file comes from OpenCV 3, the end is siril's code */

/* An example of the ECC algorithm usage can be found here:
 * https://docs.opencv.org/3.3.0/d0/d7f/image_alignment_8cpp-example.html
 */

#include <cassert>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "core/siril.h"
#include "core/proto.h"
#include "opencv/ecc/ecc.h"
#include "opencv/opencv.h"

#define ECC_DEBUG

/****************************************************************************************\
*                                       Image Alignment (ECC algorithm)                  *
 \****************************************************************************************/

using namespace cv;

static void image_jacobian_homo_ECC(const Mat& src1, const Mat& src2,
		const Mat& src3, const Mat& src4, const Mat& src5, Mat& dst) {

	assert(src1.size() == src2.size());
	assert(src1.size() == src3.size());
	assert(src1.size() == src4.size());

	assert(src1.rows == dst.rows);
	assert(dst.cols == (src1.cols * 8));
	assert(dst.type() == CV_32FC1);

	assert(src5.isContinuous());

	const float* hptr = src5.ptr<float>(0);

	const float h0_ = hptr[0];
	const float h1_ = hptr[3];
	const float h2_ = hptr[6];
	const float h3_ = hptr[1];
	const float h4_ = hptr[4];
	const float h5_ = hptr[7];
	const float h6_ = hptr[2];
	const float h7_ = hptr[5];

	const int w = src1.cols;

	//create denominator for all points as a block
	Mat den_ = src3 * h2_ + src4 * h5_ + 1.0; //check the time of this! otherwise use addWeighted

	//create projected points
	Mat hatX_ = -src3 * h0_ - src4 * h3_ - h6_;
	divide(hatX_, den_, hatX_);
	Mat hatY_ = -src3 * h1_ - src4 * h4_ - h7_;
	divide(hatY_, den_, hatY_);

	//instead of dividing each block with den,
	//just pre-devide the block of gradients (it's more efficient)

	Mat src1Divided_;
	Mat src2Divided_;

	divide(src1, den_, src1Divided_);
	divide(src2, den_, src2Divided_);

	//compute Jacobian blocks (8 blocks)

	dst.colRange(0, w) = src1Divided_.mul(src3);    //1

	dst.colRange(w, 2 * w) = src2Divided_.mul(src3);    //2

	Mat temp_ = (hatX_.mul(src1Divided_) + hatY_.mul(src2Divided_));
	dst.colRange(2 * w, 3 * w) = temp_.mul(src3);    //3

	hatX_.release();
	hatY_.release();

	dst.colRange(3 * w, 4 * w) = src1Divided_.mul(src4);    //4

	dst.colRange(4 * w, 5 * w) = src2Divided_.mul(src4);    //5

	dst.colRange(5 * w, 6 * w) = temp_.mul(src4);    //6

	src1Divided_.copyTo(dst.colRange(6 * w, 7 * w));    //7

	src2Divided_.copyTo(dst.colRange(7 * w, 8 * w));    //8
}

static void image_jacobian_euclidean_ECC(const Mat& src1, const Mat& src2,
		const Mat& src3, const Mat& src4, const Mat& src5, Mat& dst) {

	assert(src1.size() == src2.size());
	assert(src1.size() == src3.size());
	assert(src1.size() == src4.size());

	assert(src1.rows == dst.rows);
	assert(dst.cols == (src1.cols * 3));
	assert(dst.type() == CV_32FC1);

	assert(src5.isContinuous());

	const float* hptr = src5.ptr<float>(0);

	const float h0 = hptr[0];    //cos(theta)
	const float h1 = hptr[3];    //sin(theta)

	const int w = src1.cols;

	//create -sin(theta)*X -cos(theta)*Y for all points as a block -> hatX
	Mat hatX = -(src3 * h1) - (src4 * h0);

	//create cos(theta)*X -sin(theta)*Y for all points as a block -> hatY
	Mat hatY = (src3 * h0) - (src4 * h1);

	//compute Jacobian blocks (3 blocks)
	dst.colRange(0, w) = (src1.mul(hatX)) + (src2.mul(hatY));    //1

	src1.copyTo(dst.colRange(w, 2 * w));    //2
	src2.copyTo(dst.colRange(2 * w, 3 * w));    //3
}

static void image_jacobian_affine_ECC(const Mat& src1, const Mat& src2,
		const Mat& src3, const Mat& src4, Mat& dst) {

	assert(src1.size() == src2.size());
	assert(src1.size() == src3.size());
	assert(src1.size() == src4.size());

	assert(src1.rows == dst.rows);
	assert(dst.cols == (6 * src1.cols));

	assert(dst.type() == CV_32FC1);

	const int w = src1.cols;

	//compute Jacobian blocks (6 blocks)

	dst.colRange(0, w) = src1.mul(src3);    //1
	dst.colRange(w, 2 * w) = src2.mul(src3);    //2
	dst.colRange(2 * w, 3 * w) = src1.mul(src4);    //3
	dst.colRange(3 * w, 4 * w) = src2.mul(src4);    //4
	src1.copyTo(dst.colRange(4 * w, 5 * w));    //5
	src2.copyTo(dst.colRange(5 * w, 6 * w));    //6
}

static void image_jacobian_translation_ECC(const Mat& src1, const Mat& src2,
		Mat& dst) {

	assert(src1.size() == src2.size());

	assert(src1.rows == dst.rows);
	assert(dst.cols == (src1.cols * 2));
	assert(dst.type() == CV_32FC1);

	const int w = src1.cols;

	//compute Jacobian blocks (2 blocks)
	src1.copyTo(dst.colRange(0, w));
	src2.copyTo(dst.colRange(w, 2 * w));
}

static void project_onto_jacobian_ECC(const Mat& src1, const Mat& src2,
		Mat& dst) {
	/* this functions is used for two types of projections. If src1.cols ==src.cols
	 it does a blockwise multiplication (like in the outer product of vectors)
	 of the blocks in matrices src1 and src2 and dst
	 has size (number_of_blcks x number_of_blocks), otherwise dst is a vector of size
	 (number_of_blocks x 1) since src2 is "multiplied"(dot) with each block of src1.

	 The number_of_blocks is equal to the number of parameters we are lloking for
	 (i.e. rtanslation:2, euclidean: 3, affine: 6, homography: 8)

	 */
	assert(src1.rows == src2.rows);
	assert((src1.cols % src2.cols) == 0);
	int w;

	float* dstPtr = dst.ptr<float>(0);

	if (src1.cols != src2.cols) {    //dst.cols==1
		w = src2.cols;
		for (int i = 0; i < dst.rows; i++) {
			dstPtr[i] = (float) src2.dot(src1.colRange(i * w, (i + 1) * w));
		}
	}

	else {
		assert(dst.cols == dst.rows); //dst is square (and symmetric)
		w = src2.cols / dst.cols;
		Mat mat;
		for (int i = 0; i < dst.rows; i++) {

			mat = Mat(src1.colRange(i * w, (i + 1) * w));
			dstPtr[i * (dst.rows + 1)] = (float) pow(norm(mat), 2); //diagonal elements

			for (int j = i + 1; j < dst.cols; j++) { //j starts from i+1
				dstPtr[i * dst.cols + j] = (float) mat.dot(
						src2.colRange(j * w, (j + 1) * w));
				dstPtr[j * dst.cols + i] = dstPtr[i * dst.cols + j]; //due to symmetry
			}
		}
	}
}

static void update_warping_matrix_ECC(Mat& map_matrix, const Mat& update,
		const int motionType) {
	assert(map_matrix.type() == CV_32FC1);
	assert(update.type() == CV_32FC1);

	assert(
			motionType == WARP_MODE_TRANSLATION
					|| motionType == WARP_MODE_EUCLIDEAN
					|| motionType == WARP_MODE_AFFINE
					|| motionType == WARP_MODE_HOMOGRAPHY);

	if (motionType == WARP_MODE_HOMOGRAPHY)
		assert(map_matrix.rows == 3 && update.rows == 8);
	else if (motionType == WARP_MODE_AFFINE)
		assert(map_matrix.rows == 2 && update.rows == 6);
	else if (motionType == WARP_MODE_EUCLIDEAN)
		assert(map_matrix.rows == 2 && update.rows == 3);
	else
		assert(map_matrix.rows == 2 && update.rows == 2);

	assert(update.cols == 1);

	assert(map_matrix.isContinuous());
	assert(update.isContinuous());

	float* mapPtr = map_matrix.ptr<float>(0);
	const float* updatePtr = update.ptr<float>(0);

	if (motionType == WARP_MODE_TRANSLATION) {
		mapPtr[2] += updatePtr[0];
		mapPtr[5] += updatePtr[1];
	}
	if (motionType == WARP_MODE_AFFINE) {
		mapPtr[0] += updatePtr[0];
		mapPtr[3] += updatePtr[1];
		mapPtr[1] += updatePtr[2];
		mapPtr[4] += updatePtr[3];
		mapPtr[2] += updatePtr[4];
		mapPtr[5] += updatePtr[5];
	}
	if (motionType == WARP_MODE_HOMOGRAPHY) {
		mapPtr[0] += updatePtr[0];
		mapPtr[3] += updatePtr[1];
		mapPtr[6] += updatePtr[2];
		mapPtr[1] += updatePtr[3];
		mapPtr[4] += updatePtr[4];
		mapPtr[7] += updatePtr[5];
		mapPtr[2] += updatePtr[6];
		mapPtr[5] += updatePtr[7];
	}
	if (motionType == WARP_MODE_EUCLIDEAN) {
		double new_theta = updatePtr[0];
		if (mapPtr[3] > 0)
			new_theta += acos(mapPtr[0]);

		if (mapPtr[3] < 0)
			new_theta -= acos(mapPtr[0]);

		mapPtr[2] += updatePtr[1];
		mapPtr[5] += updatePtr[2];
		mapPtr[0] = mapPtr[4] = (float) cos(new_theta);
		mapPtr[3] = (float) sin(new_theta);
		mapPtr[1] = -mapPtr[3];
	}
}

double findTransform_ECC(InputArray templateImage, InputArray inputImage,
		InputOutputArray warpMatrix, int motionType, TermCriteria criteria,
		InputArray inputMask) {

	Mat src = templateImage.getMat(); //template image
	Mat dst = inputImage.getMat(); //input image (to be warped)
	Mat map = warpMatrix.getMat(); //warp (transformation)

	assert(!src.empty());
	assert(!dst.empty());

	if (!(src.type() == dst.type()))
		std::cerr << "Both input images must have the same data type"
				<< std::endl;

	//accept only 1-channel images
	if (src.type() != CV_8UC1 && src.type() != CV_32FC1)
		std::cerr << "Images must have 8uC1 or 32fC1 type" << std::endl;

	if (map.type() != CV_32FC1)
		std::cerr << "warpMatrix must be single-channel floating-point matrix"
				<< std::endl;

	assert(map.cols == 3);
	assert(map.rows == 2 || map.rows == 3);

	assert(motionType == WARP_MODE_AFFINE || motionType == WARP_MODE_HOMOGRAPHY
					|| motionType == WARP_MODE_EUCLIDEAN
					|| motionType == WARP_MODE_TRANSLATION);

	if (motionType == WARP_MODE_HOMOGRAPHY) {
		assert(map.rows == 3);
	}

	assert((criteria.type & TermCriteria::COUNT)
					|| (criteria.type & TermCriteria::EPS));
	const int numberOfIterations =
			(criteria.type & TermCriteria::COUNT) ? criteria.maxCount : 200;
	const double termination_eps =
			(criteria.type & TermCriteria::EPS) ? criteria.epsilon : -1;

	int paramTemp = 6; //default: affine
	switch (motionType) {
	case WARP_MODE_TRANSLATION:
		paramTemp = 2;
		break;
	case WARP_MODE_EUCLIDEAN:
		paramTemp = 3;
		break;
	case WARP_MODE_HOMOGRAPHY:
		paramTemp = 8;
		break;
	}

	const int numberOfParameters = paramTemp;

	const int ws = src.cols;
	const int hs = src.rows;
	const int wd = dst.cols;
	const int hd = dst.rows;

	Mat Xcoord = Mat(1, ws, CV_32F);
	Mat Ycoord = Mat(hs, 1, CV_32F);
	Mat Xgrid = Mat(hs, ws, CV_32F);
	Mat Ygrid = Mat(hs, ws, CV_32F);

	float* XcoPtr = Xcoord.ptr<float>(0);
	float* YcoPtr = Ycoord.ptr<float>(0);
	int j;
	for (j = 0; j < ws; j++)
		XcoPtr[j] = (float) j;
	for (j = 0; j < hs; j++)
		YcoPtr[j] = (float) j;

	repeat(Xcoord, hs, 1, Xgrid);
	repeat(Ycoord, 1, ws, Ygrid);

	Xcoord.release();
	Ycoord.release();

	Mat templateZM = Mat(hs, ws, CV_32F); // to store the (smoothed)zero-mean version of template
	Mat templateFloat = Mat(hs, ws, CV_32F); // to store the (smoothed) template
	Mat imageFloat = Mat(hd, wd, CV_32F); // to store the (smoothed) input image
	Mat imageWarped = Mat(hs, ws, CV_32F); // to store the warped zero-mean input image
	Mat imageMask = Mat(hs, ws, CV_8U); //to store the final mask

	Mat inputMaskMat = inputMask.getMat();
	//to use it for mask warping
	Mat preMask;
	if (inputMask.empty())
		preMask = Mat::ones(hd, wd, CV_8U);
	else
		threshold(inputMask, preMask, 0, 1, THRESH_BINARY);

	//gaussian filtering is optional
	src.convertTo(templateFloat, templateFloat.type());
	GaussianBlur(templateFloat, templateFloat, Size(5, 5), 0, 0);

	Mat preMaskFloat;
	preMask.convertTo(preMaskFloat, CV_32F);
	GaussianBlur(preMaskFloat, preMaskFloat, Size(5, 5), 0, 0);
	// Change threshold.
	preMaskFloat *= (0.5 / 0.95);
	// Rounding conversion.
	preMaskFloat.convertTo(preMask, preMask.type());
	preMask.convertTo(preMaskFloat, preMaskFloat.type());

	dst.convertTo(imageFloat, imageFloat.type());
	GaussianBlur(imageFloat, imageFloat, Size(5, 5), 0, 0);

	// needed matrices for gradients and warped gradients
	Mat gradientX = Mat::zeros(hd, wd, CV_32FC1);
	Mat gradientY = Mat::zeros(hd, wd, CV_32FC1);
	Mat gradientXWarped = Mat(hs, ws, CV_32FC1);
	Mat gradientYWarped = Mat(hs, ws, CV_32FC1);

	// calculate first order image derivatives
	Matx13f dx(-0.5f, 0.0f, 0.5f);

	filter2D(imageFloat, gradientX, -1, dx);
	filter2D(imageFloat, gradientY, -1, dx.t());

	gradientX = gradientX.mul(preMaskFloat);
	gradientY = gradientY.mul(preMaskFloat);

	// matrices needed for solving linear equation system for maximizing ECC
	Mat jacobian = Mat(hs, ws * numberOfParameters, CV_32F);
	Mat hessian = Mat(numberOfParameters, numberOfParameters, CV_32F);
	Mat hessianInv = Mat(numberOfParameters, numberOfParameters, CV_32F);
	Mat imageProjection = Mat(numberOfParameters, 1, CV_32F);
	Mat templateProjection = Mat(numberOfParameters, 1, CV_32F);
	Mat imageProjectionHessian = Mat(numberOfParameters, 1, CV_32F);
	Mat errorProjection = Mat(numberOfParameters, 1, CV_32F);

	Mat deltaP = Mat(numberOfParameters, 1, CV_32F); //transformation parameter correction
	Mat error = Mat(hs, ws, CV_32F); //error as 2D matrix

	const int imageFlags = INTER_LINEAR + WARP_INVERSE_MAP;
	const int maskFlags = INTER_NEAREST + WARP_INVERSE_MAP;

	// iteratively update map_matrix
	double rho = -1;
	double last_rho = -termination_eps;
	for (int i = 1;
			(i <= numberOfIterations)
					&& (fabs(rho - last_rho) >= termination_eps); i++) {

		// warp-back portion of the inputImage and gradients to the coordinate space of the templateImage
		if (motionType != WARP_MODE_HOMOGRAPHY) {
			warpAffine(imageFloat, imageWarped, map, imageWarped.size(),
					imageFlags);
			warpAffine(gradientX, gradientXWarped, map, gradientXWarped.size(),
					imageFlags);
			warpAffine(gradientY, gradientYWarped, map, gradientYWarped.size(),
					imageFlags);
			warpAffine(preMask, imageMask, map, imageMask.size(), maskFlags);
		} else {
			warpPerspective(imageFloat, imageWarped, map, imageWarped.size(),
					imageFlags);
			warpPerspective(gradientX, gradientXWarped, map,
					gradientXWarped.size(), imageFlags);
			warpPerspective(gradientY, gradientYWarped, map,
					gradientYWarped.size(), imageFlags);
			warpPerspective(preMask, imageMask, map, imageMask.size(),
					maskFlags);
		}

		Scalar imgMean, imgStd, tmpMean, tmpStd;
		meanStdDev(imageWarped, imgMean, imgStd, imageMask);
		meanStdDev(templateFloat, tmpMean, tmpStd, imageMask);

		subtract(imageWarped, imgMean, imageWarped, imageMask); //zero-mean input
		templateZM = Mat::zeros(templateZM.rows, templateZM.cols,
				templateZM.type());
		subtract(templateFloat, tmpMean, templateZM, imageMask); //zero-mean template

		const double tmpNorm = std::sqrt(
				countNonZero(imageMask) * (tmpStd.val[0]) * (tmpStd.val[0]));
		const double imgNorm = std::sqrt(
				countNonZero(imageMask) * (imgStd.val[0]) * (imgStd.val[0]));

		// calculate jacobian of image wrt parameters
		switch (motionType) {
		case WARP_MODE_AFFINE:
			image_jacobian_affine_ECC(gradientXWarped, gradientYWarped, Xgrid,
					Ygrid, jacobian);
			break;
		case WARP_MODE_HOMOGRAPHY:
			image_jacobian_homo_ECC(gradientXWarped, gradientYWarped, Xgrid,
					Ygrid, map, jacobian);
			break;
		case WARP_MODE_TRANSLATION:
			image_jacobian_translation_ECC(gradientXWarped, gradientYWarped,
					jacobian);
			break;
		case WARP_MODE_EUCLIDEAN:
			image_jacobian_euclidean_ECC(gradientXWarped, gradientYWarped,
					Xgrid, Ygrid, map, jacobian);
			break;
		}

		// calculate Hessian and its inverse
		project_onto_jacobian_ECC(jacobian, jacobian, hessian);

		hessianInv = hessian.inv();

		const double correlation = templateZM.dot(imageWarped);

		// calculate enhanced correlation coefficiont (ECC)->rho
		last_rho = rho;
		rho = correlation / (imgNorm * tmpNorm);
		if (cvIsNaN(rho)) {
			std::cerr << "NaN encountered." << std::endl;
		}

		// project images into jacobian
		project_onto_jacobian_ECC(jacobian, imageWarped, imageProjection);
		project_onto_jacobian_ECC(jacobian, templateZM, templateProjection);

		// calculate the parameter lambda to account for illumination variation
		imageProjectionHessian = hessianInv * imageProjection;
		const double lambda_n = (imgNorm * imgNorm)
				- imageProjection.dot(imageProjectionHessian);
		const double lambda_d = correlation
				- templateProjection.dot(imageProjectionHessian);
		if (lambda_d <= 0.0) {
			rho = -1;
			std::cerr
					<< "The algorithm stopped before its convergence. The correlation is going to be minimized. Images may be uncorrelated or non-overlapped"
					<< std::endl;
			;

		}
		const double lambda = (lambda_n / lambda_d);

		// estimate the update step delta_p
		error = lambda * templateZM - imageWarped;
		project_onto_jacobian_ECC(jacobian, error, errorProjection);
		deltaP = hessianInv * errorProjection;

		// update warping matrix
		update_warping_matrix_ECC(map, deltaP, motionType);

	}

	// return final correlation coefficient
	return rho;
}


/*******************************************************
 *        S I R I L             C O D E                *
 *******************************************************/

int findTransformBuf(WORD *reference, int ref_rows, int ref_cols,
		WORD *image, int im_rows, int im_cols, reg_ecc *reg_param)
{
	Mat ref(ref_rows, ref_cols, CV_16UC1, reference);
	Mat im(im_rows, im_cols, CV_16UC1, image);
	Mat warp_matrix = Mat::eye(3, 3, CV_64F);
	WARP_MODE warp_mode = WARP_MODE_TRANSLATION;	// for tests
	int number_of_iterations = 50;
	double termination_eps = 0.002;
	int retvalue = 0;

	setIdentity(warp_matrix);
	ref.convertTo(ref, CV_8UC1);
	im.convertTo(im, CV_8UC1);

	if (warp_mode != WARP_MODE_HOMOGRAPHY)
		warp_matrix.rows = 2;

	// Define termination criteria
	TermCriteria criteria (TermCriteria::COUNT+TermCriteria::EPS, number_of_iterations, termination_eps);

	double rho = findTransform_ECC(ref, im, warp_matrix, warp_mode, criteria,
			noArray());

	if (rho > 0) {

#ifdef ECC_DEBUG
		std::cout << "rho=" << rho << std::endl;
		std::cout << "result=" << std::endl << warp_matrix << std::endl;
#endif

		switch (warp_mode) {
		case WARP_MODE_TRANSLATION:
			reg_param->dx = warp_matrix.at<float>(0, 2);
			reg_param->dy = warp_matrix.at<float>(1, 2);
			break;
		default:
			std::cout << "Not handled yet" << std::endl;
			retvalue = 1;
		}
	}
	else
		retvalue = rho;

	warp_matrix.release();
	ref.release();
	im.release();
	return retvalue;
}

int findTransform(fits *reference, fits *image, int layer, reg_ecc *reg_param)
{
	return findTransformBuf(reference->pdata[layer], reference->ry, reference->rx,
			image->pdata[layer], image->ry, image->rx, reg_param);
}

/* finds the translation with ECC between a reference image and a tested image,
 * both monochrome, the same size and square of side size. It first downscales
 * the images to estimate the translation then refines it with the full images.
 */
int ecc_find_translation_buf(WORD *reference, WORD *image, int size, double im_max, // input args
		reg_ecc *reg_param) // output args
{
	// 1. create the input images
	Mat ref(size, size, CV_16UC1, reference);
	Mat im(size, size, CV_16UC1, image);
	ref.convertTo(ref, CV_32FC1, 1.0/im_max);
	im.convertTo(im, CV_32FC1, 1.0/im_max);

	// 2. create the downscaled images
	int down_size = size / 2;
	Mat down_ref(down_size, down_size, CV_32FC1);
	Mat down_im(down_size, down_size, CV_32FC1);
	resize(im, down_im, down_im.size(), 0, 0, OPENCV_LINEAR);
	resize(ref, down_ref, down_ref.size(), 0, 0, OPENCV_LINEAR);

	// 3. find basic transform for downscaled
	Mat warp_matrix = Mat::eye(2, 3, CV_32F); // 64 not supported
	WARP_MODE motion_type = WARP_MODE_TRANSLATION;
	TermCriteria criteria(TermCriteria::COUNT+TermCriteria::EPS, 200, 0.008);
	double ecc;

	ecc = findTransform_ECC(down_ref, down_im, warp_matrix, motion_type, criteria, noArray());
#ifdef ECC_DEBUG
	std::cout << "ecc down = " << ecc << std::endl;
	std::cout << "result down = " << std::endl << warp_matrix << std::endl;
#endif
	if (ecc < 0.3) {
		// failure
		return -1;
	}
	// since we use it for an image twice as big, translation is twice as big too
	warp_matrix.at<float>(0, 2) *= 2.0f;
	warp_matrix.at<float>(1, 2) *= 2.0f;

	// 4. find the best transfrom between full-size images
	criteria.maxCount = 1200;
	criteria.epsilon = 0.002;

	ecc = findTransform_ECC(ref, im, warp_matrix, motion_type, criteria, noArray());
#ifdef ECC_DEBUG
	std::cout << "ecc = " << ecc << std::endl;
	std::cout << "result = " << std::endl << warp_matrix << std::endl;
#endif
	if (ecc < 0.8) {
		// failure
		return -1;
	}
	// 5. copy result
	reg_param->dx = warp_matrix.at<float>(0, 2);
	reg_param->dy = warp_matrix.at<float>(1, 2);
	return 0;
}

/* finds the transform with ECC between a reference image and a tested image,
 * both monochrome, the same size and square of side size. It first downscales
 * the images to estimate a translation and computes a 2D transformation
 * (homography) using the full images.
 */
int ecc_find_transform_buf(WORD *reference, WORD *image, int size, double im_max, // input args
		reg_ecc *reg_param, Homography *transform) // output args
{
	// 1. create the input images
	Mat ref(size, size, CV_16UC1, reference);
	Mat im(size, size, CV_16UC1, image);
	ref.convertTo(ref, CV_32FC1, 1.0/im_max);
	im.convertTo(im, CV_32FC1, 1.0/im_max);

	// 2. create the downscaled images
	int down_size = size / 2;
	Mat down_ref(down_size, down_size, CV_32FC1);
	Mat down_im(down_size, down_size, CV_32FC1);
	resize(im, down_im, down_im.size(), 0, 0, OPENCV_LINEAR);
	resize(ref, down_ref, down_ref.size(), 0, 0, OPENCV_LINEAR);

	// 3. find basic transform for downscaled
	Mat warp_matrix = Mat::eye(2, 3, CV_64F);
	WARP_MODE motion_type = WARP_MODE_TRANSLATION;
	TermCriteria criteria(TermCriteria::COUNT+TermCriteria::EPS, 300, 0.006);
	double ecc;

	ecc = findTransform_ECC(down_ref, down_im, warp_matrix, motion_type, criteria, noArray());
#ifdef ECC_DEBUG
	std::cout << "ecc down = " << ecc << std::endl;
	std::cout << "result down = " << std::endl << warp_matrix << std::endl;
#endif
	if (ecc < 0.3) {
		// failure
		return -1;
	}
	// since we use it for an image twice as big, translation is twice as big too
	reg_param->dx = warp_matrix.at<float>(0, 2) * 2.0f;
	reg_param->dy = warp_matrix.at<float>(1, 2) * 2.0f;

	// 4. find the best transfrom between full-size images
	criteria.maxCount = 2000;
	criteria.epsilon = 0.001;
	motion_type = WARP_MODE_HOMOGRAPHY;
	Mat warp_matrix_full = Mat::eye(3, 3, CV_64F);
	/*for (int i = 0; i < 2; i++)
		for (int j = 0; j < 3; j++)
			warp_matrix_full.at<float>(i, j) = warp_matrix.at<float>(i, j);*/
	warp_matrix_full.at<float>(0, 2) = warp_matrix.at<float>(0, 2) * 2.0f;
	warp_matrix_full.at<float>(1, 2) = warp_matrix.at<float>(1, 2) * 2.0f;

	ecc = findTransform_ECC(ref, im, warp_matrix_full, motion_type, criteria, noArray());
#ifdef ECC_DEBUG
	std::cout << "ecc = " << ecc << std::endl;
	std::cout << "result = " << std::endl << warp_matrix_full << std::endl;
#endif
	if (ecc < 0.9) {
		// failure
		return -1;
	}

	// 5. copy result
	convert_MatH_to_H(warp_matrix_full, transform);
	return 0;
}
