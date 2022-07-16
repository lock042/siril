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
#define CV_RANSAC FM_RANSAC
#include <opencv2/calib3d.hpp>

#include "core/siril.h"
#include "core/proto.h"
#include "registration/registration.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "opencv.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "algos/statistics.h"
#ifdef __cplusplus
}
#endif

#define defaultRANSACReprojThreshold 3

using namespace cv;

static void convert_H_to_MatH(Homography *from, Mat &to);
static void convert_MatH_to_H(Mat from, Homography *to);

/* TODO:
 * fix memory leak
 *
 * Notice that functions are not static anymore as they are used
 * in kombat.cpp; and "static" function are only valid within current
 * translation unit in C++.
 */
WORD *fits_to_bgrbgr_ushort(fits *image) {
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

float *fits_to_bgrbgr_float(fits *image) {
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

/* this prepares input and output images, but lets the input in a non-usable state, beware!
 * the memory consumption of the combination of this and Mat_to_image is O(n) */
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
	if (image->naxes[2] == 1) {
		free(image->data);
		free(image->fdata);
	}

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

void cvResizeArray(double *in, double *out, int inX, int inY, int outX, int outY) {
	Mat in_mat(inX, inY, CV_64F, in);
	Mat out_mat(outX, outY, CV_64F, out);
	resize(in_mat, out_mat, out_mat.size(), 0, 0);
}

int cvRotateImage(fits *image, point center, double angle, int interpolation, int cropped) {
	Mat in, out;
	void *bgr = NULL;
	int target_rx = image->rx, target_ry = image->ry;
	Rect frame;
	Point2f pt(center.x, center.y);

	/*  the angle will be mapped to the range of [-360, 360] */
	angle = (fmod((angle / 90), 4)) * 90;

	gboolean is_fast = fmod(angle, 90.0) == 0.0;
	if (interpolation == -1)
		assert(is_fast);

	if (is_fast && (interpolation == -1 || !cropped)) {
		if (fmod(angle, 180.0) != 0.0) {
			target_rx = image->ry;
			target_ry = image->rx;
		}
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
		/* flip third argument: how to flip the array; 0 means flipping around the
		 * x-axis and positive value (for example, 1) means flipping around y-axis.
		 * Negative value (for example, -1) means flipping around both axes.
		 */
		if (angle == 90 || angle == -270) {
			transpose(in, out);
			flip(out, out, 0);
		} else if (angle == 180 || angle == -180) {
			flip(in, out, -1);
		}
		else { // 270, -90
			transpose(in, out);
			flip(out, out, 1);
		}
	} else {
		Mat r = getRotationMatrix2D(pt, angle, 1.0);
		if (cropped == 1) {
			warpAffine(in, out, r, out.size(), interpolation);
		} else {
			// adjust transformation matrix
			r.at<double>(0, 2) += frame.width / 2.0 - pt.x;
			r.at<double>(1, 2) += frame.height / 2.0 - pt.y;

			warpAffine(in, out, r, frame.size(), interpolation);
		}
	}

	return Mat_to_image(image, &in, &out, bgr, target_rx, target_ry);
}

void cvRotateImageRefPoint(fits *image, point center, double angle, int cropped, point refpointin, point *refpointout) {
	Point2f pt(center.x, center.y);
	Point3d refptin(refpointin.x, refpointin.y, 1);
	Mat refptout;
	Rect frame;

	if (!cropped) {
		frame = RotatedRect(pt, Size(center.x *2, center.y *2), angle).boundingRect2f();
	}

	Mat r = getRotationMatrix2D(pt, angle, 1.0);
	if (!cropped) {
		// adjust transformation matrix
		r.at<double>(0, 2) += frame.width / 2.0 - pt.x;
		r.at<double>(1, 2) += frame.height / 2.0 - pt.y;
	}
	refptout = r * Mat(refptin);
	refpointout->x = refptout.at<double>(0, 0);
	refpointout->y = refptout.at<double>(0, 1);
}


int cvAffineTransformation(fits *image, pointf *refpoints, pointf *curpoints, int nb_points, gboolean upscale2x, int interpolation, Homography *Hom) {
	// see https://docs.opencv.org/3.4/d4/d61/tutorial_warp_affine.html
	std::vector<Point2f> ref;
	std::vector<Point2f> cur;

	/* build vectors with lists of 3 stars. */
	for (int i = 0; i < nb_points; i++) {
		ref.push_back(Point2f(refpoints[i].x, image->ry - refpoints[i].y - 1));
		cur.push_back(Point2f(curpoints[i].x, image->ry - curpoints[i].y - 1));
	}

	Mat m = estimateAffinePartial2D(cur, ref);
	// std::cout << m << std::endl;

	/* test that m is not a zero matrix */
	if (countNonZero(m) < 1) {
		siril_log_color_message(_("Singular Matrix. Cannot compute Affine Transformation.\n"), "red");
		return -1;
	}
	Mat H = Mat::eye(3, 3, CV_64FC1);
	m.copyTo(H(cv::Rect_<int>(0,0,3,2))); //slicing is (x, y, w, h)
	convert_MatH_to_H(H, Hom);

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

void cvGetEye(Homography *Hom) {
	Mat M = Mat::eye(3, 3, CV_64FC1);
	convert_MatH_to_H(M, Hom);
}

void cvTransfPoint(double *x, double *y, Homography Href, Homography Himg) {
	Mat_<double> ref(3,1);
	Mat_<double> dst;
	Mat H0 = Mat(3, 3, CV_64FC1);
	Mat H1 = Mat(3, 3, CV_64FC1);

	ref(0,0) = *x;
	ref(1,0) = *y;
	ref(2,0) = 1.;
	convert_H_to_MatH(&Href, H0);
	convert_H_to_MatH(&Himg, H1);
	if (cv::determinant(H1) == 0) return;
	dst = H1.inv() * H0 * ref;
	*x = dst(0,0);
	*y = dst(1,0);
}

void cvTransfH(Homography Href, Homography Himg, Homography *Hres) {
	Mat_<double> ref(3,1);
	Mat_<double> dst;
	Mat H0 = Mat(3, 3, CV_64FC1);
	Mat H1 = Mat(3, 3, CV_64FC1);
	Mat H2 = Mat(3, 3, CV_64FC1);

	convert_H_to_MatH(&Href, H0);
	convert_H_to_MatH(&Himg, H1);
	H2 = H1.inv() * H0;
	convert_MatH_to_H(H2, Hres);
}

unsigned char *cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *Hom, transformation_type type) {

	std::vector<Point2f> ref;
	std::vector<Point2f> img;
	// needed for shift transform which uses estimateTranslation3D
#ifdef HAVE_CV44
	std::vector<Point3f> ref3;
	std::vector<Point3f> img3;
#endif
	Mat H, a, mask, s;
	unsigned char *ret = NULL;

	/* build vectors with lists of stars. */
	switch (type) {
	case SIMILARITY_TRANSFORMATION:
	case HOMOGRAPHY_TRANSFORMATION:
	case AFFINE_TRANSFORMATION:
		for (int i = 0; i < n; i++) {
			ref.push_back(Point2f(star_array_ref[i].x, star_array_ref[i].y));
			img.push_back(Point2f(star_array_img[i].x, star_array_img[i].y));
		}
	break;
#ifdef HAVE_CV44
	case SHIFT_TRANSFORMATION:
		for (int i = 0; i < n; i++) {
			ref3.push_back(Point3f(star_array_ref[i].x, star_array_ref[i].y, 0.));
			img3.push_back(Point3f(star_array_img[i].x, star_array_img[i].y, 0.));
		}
	break;
#endif
	default:
		return NULL;
	}

	//fitting the model
	switch (type) {
#ifdef HAVE_CV44
	case SHIFT_TRANSFORMATION:
		estimateTranslation3D(img3, ref3, s, mask, CV_RANSAC, defaultRANSACReprojThreshold);
		if (!s.cols) return NULL; // exit if could not find a match at all=> s is null
		H = Mat::eye(3, 3, CV_64FC1);
		H.at<double>(0,2) = s.at<double>(0);
		H.at<double>(1,2) = s.at<double>(1);
	break;
#endif
	case SIMILARITY_TRANSFORMATION:
		a = estimateAffinePartial2D(img, ref, mask, CV_RANSAC, defaultRANSACReprojThreshold);
		if (countNonZero(a) < 1) return NULL; //must count before filling H, otherwise zero elements cannot be caught
		H = Mat::eye(3, 3, CV_64FC1);
		a.copyTo(H(cv::Rect_<int>(0,0,3,2))); //slicing is (x, y, w, h)
	break;
	case HOMOGRAPHY_TRANSFORMATION:
		H = findHomography(img, ref, CV_RANSAC, defaultRANSACReprojThreshold, mask);
		if (countNonZero(H) < 1) return NULL;
		break;
	case AFFINE_TRANSFORMATION:
		a = estimateAffine2D(img, ref, mask, CV_RANSAC, defaultRANSACReprojThreshold);
		if (countNonZero(a) < 1) return NULL; //must count before filling H, otherwise zero elements cannot be caught
		H = Mat::eye(3, 3, CV_64FC1);
		a.copyTo(H(cv::Rect_<int>(0,0,3,2))); //slicing is (x, y, w, h)
		break;
	default:
		return NULL;
	}

	Hom->Inliers = countNonZero(mask);
	if (n > 0) {
		ret = (unsigned char *) malloc(n * sizeof(unsigned char));
		for (int i = 0; i < n; i++) {
			ret[i] = mask.at<uchar>(i);
		}
	} else {
		return NULL;
	}

	convert_MatH_to_H(H, Hom);

	mask.release();
	return ret;
}

// transform an image using the homography.
int cvTransformImage(fits *image, unsigned int width, unsigned int height, Homography Hom, gboolean upscale2x, int interpolation) {
    Mat in, out;
    void *bgr = NULL;
    int target_rx = width, target_ry = height;
    int source_ry = image->ry;

    if (image_to_Mat(image, &in, &out, &bgr, target_rx, target_ry))
        return 1;

    Mat H = Mat(3, 3, CV_64FC1);
    convert_H_to_MatH(&Hom, H);

    if (upscale2x) {
        Mat S = Mat::eye(3, 3, CV_64FC1);
        S.at<double>(0,0) = 2.0;
        S.at<double>(1,1) = 2.0;
        H = S * H;
    }

    /* modify matrix for reverse Y axis */
    Mat F1 = Mat::eye(3, 3, CV_64FC1);
    F1.at<double>(1,1) = -1.0;
    F1.at<double>(1,2) = source_ry - 1.0;

    Mat F2 = Mat::eye(3, 3, CV_64FC1);
    F2.at<double>(1,1) = -1.0;
    F2.at<double>(1,2) = target_ry - 1.0;

    H = F2.inv() * H * F1;

    // OpenCV function
    warpPerspective(in, out, H, Size(target_rx, target_ry), interpolation, BORDER_TRANSPARENT);

    return Mat_to_image(image, &in, &out, bgr, target_rx, target_ry);
}

int cvUnsharpFilter(fits* image, double sigma, double amount) {
	Mat in, out;
	void *bgr = NULL;
	int target_rx = image->rx, target_ry = image->ry;

	if (image_to_Mat(image, &in, &out, &bgr, target_rx, target_ry))
		return 1;

	//setUseOptimized(false);
	//std::cout << "---- OpenCV setUseOptimize(false) ----" << std::endl;
	//std::cout << getBuildInformation();

	/* 3rd argument: Gaussian kernel size. When width and height are zeros
	 * they are computed from sigma.
	 */
	siril_debug_print("using opencv GaussianBlur (CPU)\n");
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

// https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
double cvCalculRigidTransform(s_star *star_array_in,
		struct s_star *star_array_out, int n, Homography *Hom) {

	Mat out = Mat(2, n, CV_64FC1);
	Mat in = Mat(2, n, CV_64FC1);
	Mat A = Mat(2, 2, CV_64FC1);
	Mat B = Mat(2, 2, CV_64FC1);
	Mat U = Mat(2, 2, CV_64FC1);
	Mat V = Mat(2, 2, CV_64FC1);
	Mat S = Mat(2, 2, CV_64FC1);
	Mat T = Mat(2, 2, CV_64FC1);
	Mat R = Mat(2, 2, CV_64FC1);
	Mat H = Mat(3, 3, CV_64FC1);
	Mat res = Mat(3, 1, CV_64FC1);
	Mat outC = Mat(2, 1, CV_64FC1);
	Mat inC = Mat(2, 1, CV_64FC1);
	Mat shift = Mat(2, 1, CV_64FC1);

	double outCx = 0., outCy = 0., inCx = 0., inCy = 0.;
	for (int i = 0; i < n; i++) {
		outCx += star_array_out[i].x;
		outCy += star_array_out[i].y;
		inCx += star_array_in[i].x;
		inCy += star_array_in[i].y;
		out.at<double>(0,i) = star_array_out[i].x;
		out.at<double>(1,i) = star_array_out[i].y;
		in.at<double>(0,i) = star_array_in[i].x;
		in.at<double>(1,i) = star_array_in[i].y;
	}
	outCx /= n;
	outCy /= n;
	inCx /= n;
	inCy /= n;
	outC.at<double>(0,0) = outCx;
	outC.at<double>(1,0) = outCy;
	inC.at<double>(0,0) = inCx;
	inC.at<double>(1,0) = inCy;	

	for (int i = 0; i < n; i++) {
		out.at<double>(0,i) -= outCx;
		out.at<double>(1,i) -= outCy;
		in.at<double>(0,i) -= inCx;
		in.at<double>(1,i) -= inCy;
	}
	// std::cout << "out\n" << out << std::endl;
	// std::cout << "in\n" << in << std::endl;

	A = out * Mat::eye(n, n, CV_64FC1) * in.t();
	B = A.t() * A;

	// std::cout << "A\n" << A << std::endl;
	// std::cout << "B\n" << B << std::endl;

	double b = -(B.at<double>(0,0) + B.at<double>(1,1));
	double c = cv::determinant(B);

	double l1 = (-b + sqrt(b * b - 4 *c)) * 0.5;
	double l2 = (-b - sqrt(b * b - 4 *c)) * 0.5;

	double v1 = B.at<double>(0,1) / (l1 - B.at<double>(0,0));
	double v2 = B.at<double>(0,1) / (l2 - B.at<double>(0,0));
	// std::cout << "l1\n" << l1 << std::endl;
	// std::cout << "l2\n" << l2 << std::endl;

	double n1 = sqrt( 1 + v1 * v1);
	double n2 = sqrt( 1 + v2 * v2);

	V.at<double>(0,0) = v1 / n1;
	V.at<double>(1,0) = 1. / n1;
	V.at<double>(0,1) = v2 / n2;
	V.at<double>(1,1) = 1. / n2;
	// std::cout << "V\n" << V << std::endl;

	if (n == 3) {
		S = Mat::eye(2, 2, CV_64FC1);
		S.at<double>(0,0) = 1. / sqrt(l1);
		S.at<double>(1,1) = 1. / sqrt(l2);
		// std::cout << "S\n" << S << std::endl;
		U = (S * (A * V).t()).t();
	} else {
		U = A * V / sqrt(l1);
		U.at<double>(0,1) = U.at<double>(1,0);
		U.at<double>(1,1) = -U.at<double>(0,0);
	}
	// std::cout << "U\n" << U << std::endl;

	T = Mat::eye(2, 2, CV_64FC1);
	T.at<double>(1,1) = cv:: determinant(V * U.t());
	// std::cout << "T\n" << T << std::endl;

	R = V * T * U.t();
	// std::cout << "R\n" << R << std::endl;
	shift = inC - R * outC;
	// std::cout << "s\n" << shift << std::endl;


	H = Mat::eye(3, 3, CV_64FC1);
	R.copyTo(H(cv::Rect_<int>(0,0,2,2)));
	H.at<double>(0,2) = shift.at<double>(0,0);
	H.at<double>(1,2) = shift.at<double>(1,0);
	// std::cout << "H\n" << H << std::endl;

	// computing the maximum error
	double err = 0., norm2;
	for (int i = 0; i < n; i++) {
		res = H * (Mat_<double>(3,1) << star_array_out[i].x, star_array_out[i].y, 1) - (Mat_<double>(3,1) << star_array_in[i].x, star_array_in[i].y, 1);
		norm2 = cv::norm(res);
		if (norm2 > err) err = norm2;
	}

	Hom->Inliers = n;
	convert_MatH_to_H(H, Hom);

	return err;

}

void multH(Homography H1, Homography H2, Homography *Hout) {
	Mat _H1 = Mat(3, 3, CV_64FC1);
	Mat _H2 = Mat(3, 3, CV_64FC1);
	Mat _H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(&H1, _H1);
	convert_H_to_MatH(&H2, _H2);
	_H = _H1 * _H2;
	convert_MatH_to_H(_H, Hout);
}

void cvGetMatrixReframe(int x, int y, int w, int h, double angle, Homography *Hom) {
	double dx = (double)x + (double)w * 0.5;
	double dy = (double)y + (double)h * 0.5;
	Point2f pt(dx, dy);
	std::cout << pt << std::endl;
	//get rotation matrix from orginal to rotated about the new center
	Mat r = getRotationMatrix2D(pt, angle, 1.0);
	Mat H = Mat::eye(3, 3, CV_64FC1);
	r.copyTo(H(cv::Rect_<int>(0,0,3,2))); //slicing is (x, y, w, h)
	std::cout << H << std::endl;
	// invert to get from rotated to original
	H = H.inv();
	std::cout << H << std::endl;
	// shift backwards to set the top-left point at selection x,y
	H.at<double>(0, 2) -= (double)x;
	H.at<double>(1, 2) -= (double)y;
	std::cout << H << std::endl;
	convert_MatH_to_H(H, Hom);
}