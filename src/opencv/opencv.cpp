/*
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2015 team free-astro (see more in AUTHORS file)
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
 *
 * Useful links about OpenCV:
 * http://docs.opencv.org/modules/core/doc/intro.html
 * http://docs.opencv.org/modules/imgproc/doc/geometric_transformations.html#resize
 */
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "opencv2/core/version.hpp"
#if CV_MAJOR_VERSION == 2
#include "opencv/findHomography/calib3d.hpp"
#else
#if CV_MAJOR_VERSION == 4
#define CV_RANSAC FM_RANSAC
#endif
#include <opencv2/calib3d.hpp>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "opencv.h"
#include "opencv/ecc/ecc.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "algos/statistics.h"
#ifdef __cplusplus
}
#endif

#define defaultRANSACReprojThreshold 3

using namespace cv;

/* TODO:
 * fix memory leak
 */

static WORD *fits_to_bgrbgr_ushort(fits *image) {
	size_t ndata = image->rx * image->ry * 3;
	WORD *bgrbgr = (WORD *)malloc(ndata * sizeof(WORD));
	if (!bgrbgr) { PRINT_ALLOC_ERR; return NULL; }
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		bgrbgr[i + 0] = image->pdata[BLAYER][j];
		bgrbgr[i + 1] = image->pdata[GLAYER][j];
		bgrbgr[i + 2] = image->pdata[RLAYER][j];
	}
	return bgrbgr;
}

static float *fits_to_bgrbgr_float(fits *image) {
	size_t ndata = image->rx * image->ry * 3;
	float *bgrbgr = (float *)malloc(ndata * sizeof(float));
	if (!bgrbgr) { PRINT_ALLOC_ERR; return NULL; }
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		bgrbgr[i + 0] = image->fpdata[BLAYER][j];
		bgrbgr[i + 1] = image->fpdata[GLAYER][j];
		bgrbgr[i + 2] = image->fpdata[RLAYER][j];
	}
	return bgrbgr;
}

static BYTE *fits8_to_bgrbgr(fits *image) {
	size_t ndata = image->rx * image->ry * 3;
	BYTE *bgrbgr = (BYTE *)malloc(ndata * sizeof(BYTE));
	if (!bgrbgr) { PRINT_ALLOC_ERR; return NULL; }
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		bgrbgr[i + 0] = (BYTE)image->pdata[BLAYER][j];
		bgrbgr[i + 1] = (BYTE)image->pdata[GLAYER][j];
		bgrbgr[i + 2] = (BYTE)image->pdata[RLAYER][j];
	}
	return bgrbgr;
}

/* this prepares input and output images, but lets the input in a non-usable state, beware! */
static int image_to_Mat(fits *image, Mat *in, Mat *out, void **bgr, int target_rx, int target_ry) {
	if (image->naxes[2] != 1 && image->naxes[2] != 3) {
		siril_log_message(_("Images with %ld channels are not supported\n"), image->naxes[2]);
		return -1;
	}
	if (image->type == DATA_USHORT) {
		if (image->naxes[2] == 1) {
			*in = Mat(image->ry, image->rx, CV_16UC1, image->data);
			WORD *newbuf = (WORD *)calloc(target_rx * target_ry, sizeof(WORD));
			if (!newbuf) {
				PRINT_ALLOC_ERR;
				return -1;
			}
			*out = Mat(target_ry, target_rx, CV_16UC1, newbuf);
		}
		else if (image->naxes[2] == 3) {
			WORD *bgr_u = fits_to_bgrbgr_ushort(image);
			if (!bgr_u) return -1;
			free(image->data);
			image->data = NULL;
			memset(image->pdata, 0, sizeof image->pdata);
			*in = Mat(image->ry, image->rx, CV_16UC3, bgr_u);
			*out = Mat(target_ry, target_rx, CV_16UC3, Scalar(0));
			*bgr = bgr_u;
		}
	}
	else if (image->type == DATA_FLOAT) {
		if (image->naxes[2] == 1) {
			*in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
			float *newbuf = (float *)calloc(target_rx * target_ry, sizeof(float));
			if (!newbuf) {
				PRINT_ALLOC_ERR;
				return -1;
			}
			*out = Mat(target_ry, target_rx, CV_32FC1, newbuf);
		}
		else if (image->naxes[2] == 3) {
			float *bgr_f = fits_to_bgrbgr_float(image);
			if (!bgr_f) return -1;
			free(image->fdata);
			image->fdata = NULL;
			memset(image->fpdata, 0, sizeof image->fpdata);
			*in = Mat(image->ry, image->rx, CV_32FC3, bgr_f);
			*out = Mat(target_ry, target_rx, CV_32FC3, Scalar(0.0f));
			*bgr = bgr_f;
		}
	}
	else return -1;
	return 0;
}

static int Mat_to_image(fits *image, Mat *in, Mat *out, void *bgr, int target_rx, int target_ry) {
	in->release();
	if (bgr) free(bgr);
	if (image->naxes[2] == 1)
		free(image->data);

	size_t ndata = target_rx * target_ry;
	if (image->type == DATA_USHORT) {
		if (image->naxes[2] == 3) {
			size_t data_size = ndata * sizeof(WORD);
			// normally image->data is NULL here, as done in image_to_Mat
			WORD *newdata = (WORD*) realloc(image->data, data_size * image->naxes[2]);
			if (!newdata) {
				PRINT_ALLOC_ERR;
				out->release();
				return 1;
			}
			image->data = newdata;
			image->pdata[RLAYER] = image->data;
			image->pdata[GLAYER] = image->data + ndata;
			image->pdata[BLAYER] = image->data + ndata * 2;

			Mat channel[3];
			channel[0] = Mat(target_ry, target_rx, CV_16UC1, image->pdata[BLAYER]);
			channel[1] = Mat(target_ry, target_rx, CV_16UC1, image->pdata[GLAYER]);
			channel[2] = Mat(target_ry, target_rx, CV_16UC1, image->pdata[RLAYER]);

			split(*out, channel);

			channel[0].release();
			channel[1].release();
			channel[2].release();
		} else {
			image->data = (WORD *)out->data;
			image->pdata[RLAYER] = image->data;
			image->pdata[GLAYER] = image->data;
			image->pdata[BLAYER] = image->data;
		}
		out->release();
	} else {
		if (image->naxes[2] == 3) {
			size_t data_size = ndata * sizeof(float);
			float *newdata = (float *) realloc(image->fdata, data_size * image->naxes[2]);
			if (!newdata) {
				PRINT_ALLOC_ERR;
				out->release();
				return 1;
			}
			image->fdata = newdata;
			image->fpdata[RLAYER] = image->fdata;
			image->fpdata[GLAYER] = image->fdata + ndata;
			image->fpdata[BLAYER] = image->fdata + ndata * 2;

			Mat channel[3];
			channel[0] = Mat(target_ry, target_rx, CV_32FC1, image->fpdata[BLAYER]);
			channel[1] = Mat(target_ry, target_rx, CV_32FC1, image->fpdata[GLAYER]);
			channel[2] = Mat(target_ry, target_rx, CV_32FC1, image->fpdata[RLAYER]);

			split(*out, channel);

			channel[0].release();
			channel[1].release();
			channel[2].release();
		} else {
			image->fdata = (float *)out->data;
			image->fpdata[RLAYER] = image->fdata;
			image->fpdata[GLAYER] = image->fdata;
			image->fpdata[BLAYER] = image->fdata;
		}
		out->release();
	}
	image->rx = target_rx;
	image->ry = target_ry;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;
	invalidate_stats_from_fit(image);
	return 0;
}

/* resizes image to the sizes toX * toY, and stores it back in image */
int cvResizeGaussian(fits *image, int toX, int toY, int interpolation) {
	Mat in, out;
	void *bgr = NULL;

	if (image_to_Mat(image, &in, &out, &bgr, toX, toY))
		return 1;

	// OpenCV function
	resize(in, out, out.size(), 0, 0, interpolation);

	return Mat_to_image(image, &in, &out, bgr, toX, toY);
}

int cvRotateImage(fits *image, point center, double angle, int interpolation, int cropped) {
	Mat in, out;
	void *bgr = NULL;
	int target_rx = image->rx, target_ry = image->ry;
	Rect frame;
	Point2f pt(center.x, center.y);

	gboolean is_fast = fmod(angle, 90.0) == 0.0 && fmod(angle, 180) != 0.0;
	if (interpolation == -1)
		assert(is_fast);

	if (is_fast && (interpolation == -1 || !cropped)) {
		target_rx = image->ry;
		target_ry = image->rx;
	}
	else if (!cropped) {
		frame = RotatedRect(pt, Size(image->rx, image->ry), angle).boundingRect();
		target_rx = frame.width;
		target_ry = frame.height;
		siril_debug_print("after rotation, new image size will be %d x %d\n", target_rx, target_ry);
	}

	if (image_to_Mat(image, &in, &out, &bgr, target_rx, target_ry))
		return 1;

	if (is_fast && (interpolation == -1 || !cropped)) {	// fast rotation
		transpose(in, out);
		/* flip third argument: how to flip the array; 0 means flipping around the
		 * x-axis and positive value (for example, 1) means flipping around y-axis.
		 * Negative value (for example, -1) means flipping around both axes. 
		 */
		if (angle == 90.0)
			flip(out, out, 0);
		else // 270, -90
			flip(out, out, 1);
	} else {
		Mat r = getRotationMatrix2D(pt, angle, 1.0);
		if (cropped == 1) {
			warpAffine(in, out, r, in.size(), interpolation);
		} else {
			// adjust transformation matrix
			r.at<double>(0, 2) += frame.width / 2.0 - pt.x;
			r.at<double>(1, 2) += frame.height / 2.0 - pt.y;

			warpAffine(in, out, r, frame.size(), interpolation);
		}
	}

	return Mat_to_image(image, &in, &out, bgr, target_rx, target_ry);
}

int cvAffineTransformation(fits *image, pointf *refpoints, pointf *curpoints, int nb_points, gboolean upscale2x, int interpolation) {
	// see https://docs.opencv.org/3.4/d4/d61/tutorial_warp_affine.html
	std::vector<Point2f> ref;
	std::vector<Point2f> cur;

	/* build vectors with lists of 3 stars. */
	for (int i = 0; i < nb_points; i++) {
		ref.push_back(Point2f(refpoints[i].x, image->ry - refpoints[i].y - 1));
		cur.push_back(Point2f(curpoints[i].x, image->ry - curpoints[i].y - 1));
	}

	Mat m = estimateAffinePartial2D(cur, ref);
	//std::cout << m << std::endl;

	/* test that m is not a zero matrix */
	if (countNonZero(m) < 1) {
		siril_log_color_message(_("Singular Matrix. Cannot compute Affine Transformation.\n"), "red");
		return -1;
	}

	Mat in, out;
	void *bgr = NULL;
	int target_rx = image->rx, target_ry = image->ry;
	if (upscale2x) {
		target_rx *= 2;
		target_ry *= 2;
		m *= 2;
	}

	if (image_to_Mat(image, &in, &out, &bgr, target_rx, target_ry))
		return 1;

	warpAffine(in, out, m, out.size(), interpolation, BORDER_TRANSPARENT);

	return Mat_to_image(image, &in, &out, bgr, target_rx, target_ry);
}

static void convert_H_to_MatH(Homography *from, Mat &to) {
	to.at<double>(0, 0) = from->h00;
	to.at<double>(0, 1) = from->h01;
	to.at<double>(0, 2) = from->h02;
	to.at<double>(1, 0) = from->h10;
	to.at<double>(1, 1) = from->h11;
	to.at<double>(1, 2) = from->h12;
	to.at<double>(2, 0) = from->h20;
	to.at<double>(2, 1) = from->h21;
	to.at<double>(2, 2) = from->h22;
}

static void convert_MatH_to_H(Mat from, Homography *to) {
	to->h00 = from.at<double>(0, 0);
	to->h01 = from.at<double>(0, 1);
	to->h02 = from.at<double>(0, 2);
	to->h10 = from.at<double>(1, 0);
	to->h11 = from.at<double>(1, 1);
	to->h12 = from.at<double>(1, 2);
	to->h20 = from.at<double>(2, 0);
	to->h21 = from.at<double>(2, 1);
	to->h22 = from.at<double>(2, 2);
}

unsigned char *cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *Hom) {

	std::vector<Point2f> ref;
	std::vector<Point2f> img;
	Mat mask;
	unsigned char *ret = NULL;
	int i;

	/* build vectors with lists of stars. */
	for (i = 0; i < n; i++) {
		ref.push_back(Point2f(star_array_ref[i].x, star_array_ref[i].y));
		img.push_back(Point2f(star_array_img[i].x, star_array_img[i].y));
	}

	Mat H = findHomography(img, ref, CV_RANSAC, defaultRANSACReprojThreshold, mask);
	if (countNonZero(H) < 1) {
		return NULL;
	}
	Hom->Inliers = countNonZero(mask);
	if (n > 0) {
		ret = (unsigned char *) malloc(n * sizeof(unsigned char));
		for (i = 0; i < n; i++) {
			ret[i] = mask.at<uchar>(i);
		}
	}

	convert_MatH_to_H(H, Hom);

	mask.release();
	return ret;
}

// transform an image using the homography.
int cvTransformImage(fits *image, long width, long height, Homography Hom, gboolean upscale2x, int interpolation) {
	Mat in, out;
	void *bgr = NULL;
	int target_rx = image->rx, target_ry = image->ry;
	if (upscale2x) {
		target_rx *= 2;
		target_ry *= 2;
	}

	if (image_to_Mat(image, &in, &out, &bgr, target_rx, target_ry))
		return 1;

	Mat H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(&Hom, H);

	/* modify matrix for reverse Y axis */
	Mat F = Mat::eye(3, 3, CV_64FC1);
	F.at<double>(1,1) = -1.0;
	F.at<double>(1,2) = image->ry - 1.0;
	H = F * H * F.inv();

	if (upscale2x) {
		Mat S = Mat::eye(3, 3, CV_64FC1);
		S.at<double>(0,0) = 2.0;
		S.at<double>(1,1) = 2.0;
		H = S * H;
	}

	// OpenCV function
	warpPerspective(in, out, H, Size(width, height), interpolation, BORDER_TRANSPARENT);

	return Mat_to_image(image, &in, &out, bgr, target_rx, target_ry);
}

int cvUnsharpFilter(fits* image, double sigma, double amount) {
	Mat in, out;
	void *bgr = NULL;
	int target_rx = image->rx, target_ry = image->ry;

	if (image_to_Mat(image, &in, &out, &bgr, target_rx, target_ry))
		return 1;

	/* 3rd argument: Gaussian kernel size. When width and height are zeros
	 * they are computed from sigma.
	 */
	GaussianBlur(in, out, Size(), sigma);
	if (fabs(amount) > 0.0) {
		out = in * (1 + amount) + out * (-amount);
	}

	return Mat_to_image(image, &in, &out, bgr, target_rx, target_ry);
}

/* Work on grey images. If image is in RGB it must be first converted
 * in CieLAB. Then, only the first channel is applied
 */
static int cvClahe_ushort(fits *image, double clip_limit, int size) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	// preparing data
	Mat in, out;

	Ptr<CLAHE> clahe = createCLAHE();
	clahe->setClipLimit(clip_limit);
	clahe->setTilesGridSize(Size(size, size));

	if (image->naxes[2] == 3) {
		Mat lab_image;
		std::vector<Mat> lab_planes(3);
		BYTE *bgrbgr8;
		WORD *bgrbgr;

		switch (image->bitpix) {
			case BYTE_IMG:
				bgrbgr8 = fits8_to_bgrbgr(image);
				in = Mat(image->ry, image->rx, CV_8UC3, bgrbgr8);
				out = Mat();
				// convert the RGB color image to Lab
				cvtColor(in, lab_image, COLOR_BGR2Lab);

				// Extract the L channel
				split(lab_image, lab_planes); // now we have the L image in lab_planes[0]

				// apply the CLAHE algorithm to the L channel
				clahe->apply(lab_planes[0], lab_planes[0]);

				// Merge the color planes back into an Lab image
				merge(lab_planes, lab_image);

				// convert back to RGB
				cvtColor(lab_image, out, COLOR_Lab2BGR);
				out.convertTo(out, CV_16UC3, 1.0);

				free(bgrbgr8);
				break;
			default:
			case USHORT_IMG:
				bgrbgr = fits_to_bgrbgr_ushort(image);
				in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr);
				in.convertTo(in, CV_32F, 1.0 / USHRT_MAX_DOUBLE);
				out = Mat();

				// convert the RGB color image to Lab
				cvtColor(in, lab_image, COLOR_BGR2Lab);

				// Extract the L channel
				split(lab_image, lab_planes); // now we have the L image in lab_planes[0]

				// apply the CLAHE algorithm to the L channel (does not work with 32F images)
				lab_planes[0].convertTo(lab_planes[0], CV_16U,	USHRT_MAX_DOUBLE / 100.0);
				clahe->apply(lab_planes[0], lab_planes[0]);
				lab_planes[0].convertTo(lab_planes[0], CV_32F, 100.0 / USHRT_MAX_DOUBLE);

				// Merge the color planes back into an Lab image
				merge(lab_planes, lab_image);

				// convert back to RGB
				cvtColor(lab_image, out, COLOR_Lab2BGR);
				out.convertTo(out, CV_16UC3, USHRT_MAX_DOUBLE);

				free(bgrbgr);
		}

	} else {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat();
		switch (image->bitpix) {
			case BYTE_IMG:
				in.convertTo(in, CV_8U, 1.0);
				clahe->apply(in, out);
				out.convertTo(out, CV_16UC3, 1.0);
				// dynamic range is important with CLAHE, use 16 bits output
				break;
			default:
			case USHORT_IMG:
				clahe->apply(in, out);
		}
	}

	std::vector<Mat> channel(3);
	split(out, channel);

	size_t nbpixels = image->naxes[0] * image->naxes[1];
	size_t ndata = nbpixels * sizeof(WORD);
	if (image->naxes[2] == 3) {
		memcpy(image->data, channel[2].data, ndata);
		memcpy(image->data + nbpixels, channel[1].data, ndata);
		memcpy(image->data + nbpixels * 2, channel[0].data, ndata);
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + nbpixels;
		image->pdata[BLAYER] = image->data + nbpixels* 2;
	} else {
		memcpy(image->data, channel[0].data, ndata);
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	}

	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	/* free data */
	in.release();
	out.release();
	channel[0].release();
	channel[1].release();
	channel[2].release();
	invalidate_stats_from_fit(image);

	return 0;
}

static int cvClahe_float(fits *image, double clip_limit, int size) {
	assert(image->fdata);
	assert(image->rx);
	assert(image->ry);

	// preparing data
	Mat in, out;

	Ptr<CLAHE> clahe = createCLAHE();
	clahe->setClipLimit(clip_limit);
	clahe->setTilesGridSize(Size(size, size));

	if (image->naxes[2] == 3) {
		Mat lab_image;
		std::vector<Mat> lab_planes(3);
		float *bgrbgr;

		bgrbgr = fits_to_bgrbgr_float(image);
		in = Mat(image->ry, image->rx, CV_32FC3, bgrbgr);
		out = Mat();

		// convert the RGB color image to Lab
		cvtColor(in, lab_image, COLOR_BGR2Lab);

		// Extract the L channel
		split(lab_image, lab_planes); // now we have the L image in lab_planes[0]

		// apply the CLAHE algorithm to the L channel (does not work with 32F images)
		lab_planes[0].convertTo(lab_planes[0], CV_16U, USHRT_MAX_DOUBLE / 100.0);
		clahe->apply(lab_planes[0], lab_planes[0]);
		lab_planes[0].convertTo(lab_planes[0], CV_32F, 100.0 / USHRT_MAX_DOUBLE);

		// Merge the color planes back into an Lab image
		merge(lab_planes, lab_image);

		// convert back to RGB
		cvtColor(lab_image, out, COLOR_Lab2BGR);
		out.convertTo(out, CV_32FC3, 1.0);

		free(bgrbgr);

	} else {
		in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
		out = Mat();

		in.convertTo(in, CV_16U, USHRT_MAX_DOUBLE);
		clahe->apply(in, out);
		out.convertTo(out, CV_32F, 1.0 / USHRT_MAX_DOUBLE);
	}

	std::vector<Mat> channel(3);
	split(out, channel);

	size_t nbpixels = image->naxes[0] * image->naxes[1];
	size_t ndata = nbpixels * sizeof(float);
	if (image->naxes[2] == 3) {
		memcpy(image->fdata, channel[2].data, ndata);
		memcpy(image->fdata + nbpixels, channel[1].data, ndata);
		memcpy(image->fdata + nbpixels * 2, channel[0].data, ndata);
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata + nbpixels;
		image->fpdata[BLAYER] = image->fdata + nbpixels* 2;
	} else {
		memcpy(image->fdata, channel[0].data, ndata);
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata;
		image->fpdata[BLAYER] = image->fdata;
	}

	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	/* free data */
	in.release();
	out.release();
	channel[0].release();
	channel[1].release();
	channel[2].release();
	invalidate_stats_from_fit(image);

	return 0;
}

int cvClahe(fits *image, double clip_limit, int size) {
	if (image->type == DATA_USHORT)
		return cvClahe_ushort(image, clip_limit, size);
	if (image->type == DATA_FLOAT)
		return cvClahe_float(image, clip_limit, size);
	return -1;
}

uint cvHoughLines(fits *image, int idx, int layer, uchar threshvalue, uint minlen) {
	Mat src, gray, thresh;
	std::vector<Vec2f> lines; // will hold the results of the detection

	if (layer < 0) layer = (image->naxes[2] == 3) ? 1 : 0; // if no layer defined, select green layer if color, first layer if mono
	if (image->type == DATA_USHORT) {
		src = Mat(image->ry, image->rx, CV_16UC1, image->pdata[layer]);
	} else {
		src = Mat(image->ry, image->rx, CV_32FC1, image->fpdata[layer]);
	}

	src.convertTo(gray, CV_8UC1); //converting to UCHAR for thresholding
	threshold(gray, thresh, threshvalue, 255, THRESH_BINARY);
    HoughLines(thresh, lines, 1, CV_PI/180, minlen, 0, 0 );

#ifdef SIRIL_OUTPUT_DEBUG
    // Output the lines if any
    for( size_t i = 0; i < lines.size(); i++ )
    {
        float rho = lines[i][0], theta = lines[i][1];
        Point pt1, pt2;
        double a = cos(theta), b = sin(theta);
        double x0 = a*rho, y0 = b*rho;
        pt1.x = cvRound(x0 + 1000*(-b));
        pt1.y = cvRound(y0 + 1000*(a));
        pt2.x = cvRound(x0 - 1000*(-b));
        pt2.y = cvRound(y0 - 1000*(a));
		siril_debug_print(_("Line detected in frame# %d. (%d,%d)->(%d,%d)\n"), idx, pt1.x, pt1.y, pt2.x, pt2.y);
    }
#endif
	src.release();
	gray.release();
	thresh.release();
	return lines.size();
}


