/*
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
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
#include <opencv2/core/matx.hpp>
#define CV_RANSAC FM_RANSAC
#include <opencv2/calib3d.hpp>

// for mosaics
#include <opencv2/stitching/warpers.hpp>

#include "core/siril.h"
#include "core/proto.h"
#include "core/masks.h"
#include "core/siril_log.h"
#include "core/settings.h"
#include "gui/callbacks.h"
#include "registration/registration.h"
#include "registration/matching/atpmatch.h"
#include "opencv.h"
#include "guidedfilter.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "io/image_format_fits.h"
#include "algos/statistics.h"
#include "registration/distorsion.h"
#ifdef __cplusplus
}
#endif

#define defaultRANSACReprojThreshold 3
#define CLAMPING_FACTOR 0.98

using namespace cv;

static void convert_H_to_MatH(const Homography *from, Mat &to);
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
		image->data = NULL; // This must be set to NULL as the function can return on error before it is reassigned to point at the modified data
		free(image->fdata);
		image->fdata = NULL; // This must be set to NULL as the function can return on error before it is reassigned to point at the modified data
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
			image->fdata = (float *) out->data;
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

// int guide image for clamping
static void init_guide(fits *image, unsigned int target_rx, unsigned int target_ry, Mat *guide) {
	if (image->type == DATA_USHORT) {
		*guide = Mat(target_ry, target_rx, (image->naxes[2] == 1) ? CV_16UC1 : CV_16UC3, Scalar(0));
	} else {
		*guide = Mat(target_ry, target_rx, (image->naxes[2] == 1) ? CV_32FC1 : CV_32FC3, Scalar(0.0f));
	}
}

/* resizes image to the sizes toX * toY, and stores it back in image */
int cvResizeGaussian(fits *image, int toX, int toY, int interpolation, gboolean clamp) {
	Mat in, out;
	void *bgr = NULL;
	if (image_to_Mat(image, &in, &out, &bgr, toX, toY))
		return 1;

	// OpenCV function
	resize(in, out, out.size(), 0, 0, interpolation);

	if ((interpolation == OPENCV_LANCZOS4 || interpolation == OPENCV_CUBIC) && clamp) {
		Mat guide, tmp1;
		init_guide(image, toX, toY, &guide);
		// Create guide image
		resize(in, guide, out.size(), 0, 0, OPENCV_AREA);
		tmp1 = (out < CLAMPING_FACTOR * guide);
		Mat element = getStructuringElement(MORPH_ELLIPSE, Size(3, 3), Point(1,1));
		dilate(tmp1, tmp1, element);

		copyTo(guide, out, tmp1); // Guide copied to the clamped pixels
		guide.release();
		tmp1.release();
	}
	return Mat_to_image(image, &in, &out, bgr, toX, toY);
}

void cvResizeArray(double *in, double *out, int inX, int inY, int outX, int outY) {
	Mat in_mat(inX, inY, CV_64F, in);
	Mat out_mat(outX, outY, CV_64F, out);
	resize(in_mat, out_mat, out_mat.size(), 0, 0);
}

// This function is now used only for mod90 rotations w/o interp
int cvRotateImage(fits *image, int angle) {
	Mat in, out;
	void *bgr = NULL;
	int target_rx = image->rx, target_ry = image->ry;

	gboolean is_fast = (angle % 90 == 0);
	if (!is_fast) return 1;

	if (angle % 180 != 0.0) {
		target_rx = image->ry;
		target_ry = image->rx;
	}

	if (image_to_Mat(image, &in, &out, &bgr, target_rx, target_ry))
		return 1;

	// fast rotation
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
	return Mat_to_image(image, &in, &out, bgr, target_rx, target_ry);
}

void cvTransformImageRefPoint(Homography Hom, point refpointin, point *refpointout) {
	Mat refptout;
	Point3d refptin(refpointin.x, refpointin.y, 1);
	Mat H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(&Hom, H);
	refptout = H * Mat(refptin);
	refpointout->x = refptout.at<double>(0) / refptout.at<double>(2);
	refpointout->y = refptout.at<double>(1) / refptout.at<double>(2);
}

static void convert_H_to_MatH(const Homography *from, Mat &to) {
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
	convert_MatH_to_H(std::move(M), Hom);
}

void cvTransfPoint(double *x, double *y, Homography Href, Homography Himg, double scale) {
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
	if (scale != 1.) {
		Mat S = Mat::eye(3, 3, CV_64FC1);
		S.at<double>(0,0) = scale;
		S.at<double>(1,1) = scale;
		dst = S * dst;
	}
	*x = dst(0,0) / dst(2,0);
	*y = dst(1,0) / dst(2,0);
}

void cvTransfH(Homography *Href, Homography *Himg, Homography *Hres) {
	Mat_<double> ref(3,1);
	Mat_<double> dst;
	Mat H0 = Mat(3, 3, CV_64FC1);
	Mat H1 = Mat(3, 3, CV_64FC1);
	Mat H2 = Mat(3, 3, CV_64FC1);

	convert_H_to_MatH(Href, H0);
	convert_H_to_MatH(Himg, H1);
	H2 = H1.inv() * H0;
	convert_MatH_to_H(std::move(H2), Hres);
}

unsigned char *cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *Hom, transformation_type type, float offset) {

	std::vector<Point2f> ref;
	std::vector<Point2f> img;
	// needed for shift transform which uses estimateTranslation3D
#ifdef HAVE_CV44
	std::vector<Point3f> ref3;
	std::vector<Point3f> img3;
#endif
	Mat H, a, mask, s;
	unsigned char *ret = NULL;

	/* 	build vectors with lists of stars.
		When defined with offset, the -0.5 term comes from the difference
		in convention between how we compute PSF
		and opencv (zero coordinate at edge of pixel vs at center) */
	switch (type) {
	case SIMILARITY_TRANSFORMATION:
	case HOMOGRAPHY_TRANSFORMATION:
	case AFFINE_TRANSFORMATION:
		for (int i = 0; i < n; i++) {
			ref.push_back(Point2f(star_array_ref[i].x + offset, star_array_ref[i].y + offset));
			img.push_back(Point2f(star_array_img[i].x + offset, star_array_img[i].y + offset));
		}
	break;
#ifdef HAVE_CV44
	case SHIFT_TRANSFORMATION:
		for (int i = 0; i < n; i++) {
			ref3.push_back(Point3f(star_array_ref[i].x + offset, star_array_ref[i].y + offset, 0.));
			img3.push_back(Point3f(star_array_img[i].x + offset, star_array_img[i].y + offset, 0.));
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

	convert_MatH_to_H(std::move(H), Hom);

	mask.release();
	return ret;
}

static void cvPrepareH(Mat H, double scale, int source_rx, int source_ry, int target_rx, int target_ry) {
	if (scale != 1.) {
		Mat S = Mat::eye(3, 3, CV_64FC1);
		S.at<double>(0,0) = scale;
		S.at<double>(1,1) = scale;
		// the terms below are the corrections to make the scaling about the image center (and not top left)
		// it's the result of -T2^-1*So*T1
		// where T1 and T2 are the shifts to the center in original and scaled images 
		// using opencv convention, so (rx/2-0.5;ry/2-0.5) and (target_rx/2-0.5;target_ry/2-0.5)
		// and So the scale about top left origin which is simply [s 0 0][0 s 0][0 0 1]
		// rx/ry are the sizes of the area to be projected, that is dst_rx/scale and dst_ry/scale
		// We then simplify as 0.5 * (1 - scale)
		S.at<double>(0,2) -= 0.5 * (1 - scale);
		S.at<double>(1,2) -= 0.5 * (1 - scale);
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
}

void cvPrepareDrizzleH(Homography *Hom, double scale, int source_rx, int source_ry, int target_rx, int target_ry) {
	Mat H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(Hom, H);
	cvPrepareH(H, scale, source_rx, source_ry, target_rx, target_ry);
	convert_MatH_to_H(std::move(H), Hom);
}

// transform an image using the homography.
int cvTransformImage(fits *image, unsigned int width, unsigned int height, Homography Hom, float scale, int interpolation, gboolean clamp, disto_data *disto) {
	Mat in, out;
	void *bgr = NULL;
	int target_rx = width, target_ry = height;

	if (image_to_Mat(image, &in, &out, &bgr, target_rx, target_ry))
		return 1;

	Mat H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(&Hom, H);
	cvPrepareH(H, scale, image->rx, image->ry, target_rx, target_ry);

	// no distortion case
	if (!disto || (disto->dtype != DISTO_MAP_D2S && disto->dtype != DISTO_D2S)) {
		// OpenCV function
		warpPerspective(in, out, H, Size(target_rx, target_ry), interpolation, BORDER_TRANSPARENT);
		if ((interpolation == OPENCV_LANCZOS4 || interpolation == OPENCV_CUBIC) && clamp) {
			Mat guide, tmp1;
			init_guide(image, target_rx, target_ry, &guide);
			// Create guide image
			warpPerspective(in, guide, H, Size(target_rx, target_ry), OPENCV_AREA, BORDER_TRANSPARENT);
			tmp1 = (out < guide * CLAMPING_FACTOR);
			Mat element = getStructuringElement( MORPH_ELLIPSE,
						Size(3, 3), Point(-1,-1));
			dilate(tmp1, tmp1, element);

			copyTo(guide, out, tmp1); // Guide copied to the clamped pixels
			guide.release();
			tmp1.release();
		}
		return Mat_to_image(image, &in, &out, bgr, target_rx, target_ry);
	}

	// distortion case - we need to compute maps to pass to remapping to avoid double interpolation
	float *xmap, *ymap;
	Mat Hinv = H.inv();  // dst->src
	xmap = (float *)malloc(target_rx * target_ry *sizeof(float));
	ymap = (float *)malloc(target_rx * target_ry *sizeof(float));
	if (!xmap || !ymap) {
		free(xmap);
		free(ymap);
		PRINT_ALLOC_ERR;
		return 2;
	}
	prepare_H_with_disto_4remap(Hinv.ptr<double>(0), image->rx, image->ry, target_rx, target_ry, disto, xmap, ymap);
	Mat uxmap = Mat(target_ry, target_rx, CV_32FC1, xmap);
	Mat uymap = Mat(target_ry, target_rx, CV_32FC1, ymap);
	remap(in, out, uxmap, uymap, interpolation, BORDER_TRANSPARENT);
	if ((interpolation == OPENCV_LANCZOS4 || interpolation == OPENCV_CUBIC) && clamp) {
		Mat guide, tmp1;
		init_guide(image, target_rx, target_ry, &guide);
		// Create guide image
		remap(in, guide, uxmap, uymap, OPENCV_AREA, BORDER_TRANSPARENT);
		tmp1 = (out < guide * CLAMPING_FACTOR);
		Mat element = getStructuringElement( MORPH_ELLIPSE,
					Size(3, 3), Point(-1,-1));
		dilate(tmp1, tmp1, element);

		copyTo(guide, out, tmp1); // Guide copied to the clamped pixels
		guide.release();
		tmp1.release();
	}
	free(xmap);
	free(ymap);
	return Mat_to_image(image, &in, &out, bgr, width, height);
}

void cvDownscaleBlendMask(int rx, int ry, int out_rx, int out_ry, uint8_t *maskin, float *maskout) {
	Mat _maskin = Mat(ry, rx, CV_8U, maskin);
	Mat _maskindown = Mat(out_ry + 2, out_rx + 2, CV_8U, Scalar(0));
	Mat _maskoutdown32 = Mat(out_ry + 2, out_rx + 2, CV_32F, Scalar(0.));
	Mat _maskout = Mat(out_ry, out_rx, CV_32F, maskout);
	int morph_size = 3;
	// first we dilate and erode to make sure we fill holes left by drizzling if any
	Mat element = getStructuringElement(MORPH_RECT,
			Size(2 * morph_size + 1, 2 * morph_size + 1),
			Point(morph_size, morph_size));
  	dilate(_maskin, _maskin, element);
	erode(_maskin, _maskin, element);
	// we leave a border of one pixel black to make sure the distance transform
	// will later count distance from this border (if the initial mask was entirely white for instance)
	// we get a handle to the inner part of the maskout (without the border)
	Mat _maskindownroi = _maskindown(Rect(1, 1, out_rx, out_ry));
	// we resize
	resize(_maskin, _maskindownroi, _maskindownroi.size(), 0, 0, INTER_LINEAR);
	// we compute the distances, it returns a 32b array
	distanceTransform(_maskindown, _maskoutdown32, DIST_L2, 3, CV_32F);
	Mat _maskoutdownroi = _maskoutdown32(Rect(1, 1, out_rx, out_ry));
	_maskoutdownroi.copyTo(_maskout);
}

void cvUpscaleBlendMask(int rx, int ry, int out_rx, int out_ry, float *maskin, float *maskout) {
	Mat _maskin = Mat(ry, rx, CV_32F, maskin);
	Mat _maskout = Mat(out_ry, out_rx, CV_32F, maskout);
	resize(_maskin, _maskout, _maskout.size(), 0, 0, INTER_LINEAR);
	flip(_maskout, _maskout, 0);  // don't know why we need to flip but does not work without
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

int cvBilateralFilter(fits* image, double d, double sigma_col, double sigma_space) {
	Mat in, out;
	void *bgr = NULL;
	if (image_to_Mat(image, &in, &out, &bgr, image->rx, image->ry))
		return 1;
	siril_debug_print("using opencv Bilateral Filter (CPU)\n");
	bilateralFilter(in, out, d, sigma_col, sigma_space, BORDER_DEFAULT);
	return Mat_to_image(image, &in, &out, bgr, image->rx, image->ry);
}

int cvGuidedFilter(fits* image, fits *guide, double r, double eps) {
	Mat in, out, guide_mat;
	void *bgr = NULL;
	int rx = guide->rx, ry = guide->ry;

	if (image_to_Mat(image, &in, &out, &bgr, image->rx, image->ry))
		return 1;
	if (image == guide) {
		guide_mat = in.clone();
	} else {
		if (guide->type == DATA_USHORT) {
			if (guide->naxes[2] == 1) {
				guide_mat = Mat(ry, rx, CV_16UC1, guide->data);
			} else if (guide->naxes[2] == 3) {
				WORD *bgr_u = fits_to_bgrbgr_ushort(guide);
				if (!bgr_u) return -1;
				guide_mat = Mat(ry, rx, CV_16UC3, bgr_u);
			}
		}
		else if (guide->type == DATA_FLOAT) {
			if (guide->naxes[2] == 1) {
				guide_mat = Mat(ry, rx, CV_32FC1, guide->fdata);
			}
			else if (guide->naxes[2] == 3) {
				float *bgr_f = fits_to_bgrbgr_float(guide);
				if (!bgr_f) return -1;
				guide_mat = Mat(ry, rx, CV_32FC3, bgr_f);
			}
		}
	}
	if (guide_mat.channels() != in.channels())
		cvtColor(guide_mat, guide_mat, COLOR_GRAY2BGR);
	siril_debug_print("using Guided Filter (CPU)\n");
	Mat result = guidedFilter(guide_mat, in, r, eps, -1);
	result.copyTo(out);
	guide_mat.release();
	return Mat_to_image(image, &in, &out, bgr, image->rx, image->ry);
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
				out.convertTo(out, CV_16UC1, 1.0);
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
		image->pdata[BLAYER] = image->data + nbpixels * 2;
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
		image->fpdata[BLAYER] = image->fdata + nbpixels * 2;
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
	/* the -0.5 term comes from the difference in convention between how we compute PSF
	/ and opencv (zero coordinate at edge of pixel vs at center) */
	for (int i = 0; i < n; i++) {
		outCx += star_array_out[i].x - 0.5;
		outCy += star_array_out[i].y - 0.5;
		inCx += star_array_in[i].x - 0.5;
		inCy += star_array_in[i].y - 0.5;
		out.at<double>(0,i) = star_array_out[i].x - 0.5;
		out.at<double>(1,i) = star_array_out[i].y - 0.5;
		in.at<double>(0,i) = star_array_in[i].x - 0.5;
		in.at<double>(1,i) = star_array_in[i].y - 0.5;
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
		res = H * (Mat_<double>(3,1) << star_array_out[i].x - 0.5, star_array_out[i].y - 0.5, 1) - (Mat_<double>(3,1) << star_array_in[i].x - 0.5, star_array_in[i].y - 0.5, 1);
		norm2 = cv::norm(res);
		if (norm2 > err) err = norm2;
	}

	Hom->Inliers = n;
	convert_MatH_to_H(std::move(H), Hom);

	return err;

}

void cvMultH(Homography H1, Homography H2, Homography *Hout) {
	Mat _H1 = Mat(3, 3, CV_64FC1);
	Mat _H2 = Mat(3, 3, CV_64FC1);
	Mat _H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(&H1, _H1);
	convert_H_to_MatH(&H2, _H2);
	_H = _H1 * _H2;
	convert_MatH_to_H(std::move(_H), Hout);
}

void cvGetMatrixReframe(double x, double y, int w, int h, double angle, Homography *Hom) {
	double dx = x + (double)w * 0.5;
	double dy = y + (double)h * 0.5;
	Point2f pt(0., 0.);

	// shift to get to the center of the initial image
	Mat S =  Mat::eye(3, 3, CV_64FC1);
	S.at<double>(0, 2) = dx;
	S.at<double>(1, 2) = dy;

	// shift backwards to set the top-left point of the final image
	Mat S2 =  Mat::eye(3, 3, CV_64FC1);
	S2.at<double>(0, 2) = -(double)w * 0.5;
	S2.at<double>(1, 2) = -(double)h * 0.5;

	// get rot matrix about origin {0, 0}
	Mat r = getRotationMatrix2D(pt, angle, 1.0);
	Mat H = Mat::eye(3, 3, CV_64FC1);
	r.copyTo(H(cv::Rect_<int>(0,0,3,2))); //slicing is (x, y, w, h)
	// std::cout << H << std::endl;

	H = S * H * S2;
	std::cout << H << std::endl;

	// transform is final to orginal, we need to inverse
	// to have H from original to final
	H = H.inv();
	std::cout << H << std::endl;

	convert_MatH_to_H(std::move(H), Hom);
}

void cvGetBoundingRectSize(fits *image, point center, double angle, int *w, int *h) {
	Rect frame;
	Point2f pt(center.x, center.y);
	frame = RotatedRect(pt, Size(image->rx, image->ry), angle).boundingRect2f();
	siril_debug_print("after rotation, new image size will be %d x %d\n", frame.width, frame.height);
	*w = (int)frame.width;
	*h = (int)frame.height;
}

void cvInvertH(Homography *Hom) {
	Mat H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(Hom, H);
	H = H.inv();
	convert_MatH_to_H(std::move(H), Hom);
	// std::cout << "H" << std::endl;
	// std::cout << H << std::endl;
}

void cvApplyFlips(Homography *Hom, int source_ry, int target_ry) {
	Mat H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(Hom, H);
	/* modify matrix for reverse Y axis */
	Mat F1 = Mat::eye(3, 3, CV_64FC1);
	F1.at<double>(1,1) = -1.0;
	F1.at<double>(1,2) = source_ry;

	Mat F2 = Mat::eye(3, 3, CV_64FC1);
	F2.at<double>(1,1) = -1.0;
	F2.at<double>(1,2) = target_ry;

	H = F2.inv() * H * F1;
	convert_MatH_to_H(std::move(H), Hom);
}


// Used to convert a H matrix written in display convention to opencv convention
void cvdisplay2ocv(Homography *Hom) {
	Mat H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(Hom, H);

	/* modify matrix to go to opencv convention */
	Mat S1 = Mat::eye(3, 3, CV_64FC1);
	S1.at<double>(0,2) = 0.5;
	S1.at<double>(1,2) = 0.5;

	H = S1.inv() * H * S1;
	convert_MatH_to_H(std::move(H), Hom);
}

/*
Takes an angle in degrees and a rotationtype (X, Y or Z)
Fills the 3D rotation matrix
Returns FALSE on error
*/
static gboolean cvRotMat1(double angle, rotation_type rottype, Mat &R) {
	angle *= G_PI / 180.;
	double ca = cos(angle);
	double sa = sin(angle);
	gboolean retval = TRUE;

	switch (rottype) {
		case ROTX:
			R.at<double>(1,1) = ca;
			R.at<double>(2,2) = ca;
			R.at<double>(1,2) = -sa;
			R.at<double>(2,1) = sa;
			break;
		case ROTY:
			R.at<double>(0,0) = ca;
			R.at<double>(2,2) = ca;
			R.at<double>(2,0) = -sa;
			R.at<double>(0,2) = sa;
			break;
		case ROTZ:
			R.at<double>(0,0) = ca;
			R.at<double>(1,1) = ca;
			R.at<double>(0,1) = -sa;
			R.at<double>(1,0) = sa;
			break;
		default:
			retval = FALSE;
			break;
	}
	return retval;
}


/*
Takes an angle array[3] in degrees and a rotationtype vector[3] (X, Y or Z)
Fills the 3D rotation matrix
The angles are applied in the input order
The flag W2C specifies if the matrix out is World-to-Camera or Camera-to-World
Returns FALSE on error
*/
gboolean cvRotMat3(double angles[3], rotation_type rottype[3], gboolean W2C, Homography *R) {
	Mat _R = Mat::eye(3, 3, CV_64FC1);
	gboolean retval = TRUE;
	for (int i = 0; i < 3; i++) {
		Mat M = Mat::eye(3, 3, CV_64FC1);
		if (!cvRotMat1(angles[i], rottype[i], M)) {
			retval = FALSE;
			break;
		}
		_R = M * _R; // left-multiply as the full matrix is R3*R2*R1
	}
	if (W2C)
		_R = _R.t(); // if we need W2C, we need to invert (i.e. transpose for R mats)
	convert_MatH_to_H(std::move(_R), R);
	// std::cout << R << std::endl;
	return retval;
}

// Computes the relative rotation matrix between Rref and R
void cvRelRot(Homography *Ref, Homography *R, Homography *Rout) {
	Mat _Ref = Mat(3, 3, CV_64FC1);
	Mat _R = Mat(3, 3, CV_64FC1);
	Mat _H = Mat(3, 3, CV_64FC1);
	Mat _Rout = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(Ref, _Ref);
	convert_H_to_MatH(R, _R);
	_Rout = _Ref.t() * _R;
	convert_MatH_to_H(std::move(_Rout), Rout);
}

// Computes Homography from cameras R and K
void cvcalcH_fromKKR(Homography *Kref, Homography *K, Homography *R, Homography *H) {
	Mat _Kref = Mat(3, 3, CV_64FC1);
	Mat _K = Mat(3, 3, CV_64FC1);
	Mat _R = Mat(3, 3, CV_64FC1);
	Mat _H = Mat(3, 3, CV_64FC1);
	convert_H_to_MatH(Kref, _Kref);
	convert_H_to_MatH(K, _K);
	convert_H_to_MatH(R, _R);

	//Compute H and returning
	_H = _Kref * _R * _K.inv();
	convert_MatH_to_H(std::move(_H), H);
}

int cvCalcH_from_corners(double *x_img, double *y_img, double *x_ref, double *y_ref, Homography *Hom) {
	std::vector<Point2f> ref;
	std::vector<Point2f> img;

	for (int i = 0; i < 4; i++) {
		ref.push_back(Point2f(x_ref[i], y_ref[i]));
		img.push_back(Point2f(x_img[i], y_img[i]));
	}

	Mat H = findHomography(img, ref, 0);
	if (H.empty())
		return -1;
	convert_MatH_to_H(std::move(H), Hom);
	return 0;
}

/******************************************
 * Mask functions that make use of OpenCV *
 * ***************************************/

// Apply Gaussian blur to the mask
int mask_apply_gaussian_blur(fits *fit, float radius) {
	if (!fit || !fit->mask || !fit->mask->data) {
		siril_debug_print("mask_apply_gaussian_blur: invalid mask\n");
		return 1;
	}
	if (radius <= 0.f) {
		siril_debug_print("mask_apply_gaussian_blur: radius must be positive\n");
		return 1;
	}

	size_t rx = fit->rx, ry = fit->ry;
	int cv_type;

	// Determine OpenCV type based on bitpix
	switch (fit->mask->bitpix) {
		case 8:
			cv_type = CV_8UC1;
			break;
		case 16:
			cv_type = CV_16UC1;
			break;
		case 32:
			cv_type = CV_32FC1;
			break;
		default:
			siril_debug_print("mask_apply_gaussian_blur: unsupported bitpix %d\n", fit->mask->bitpix);
			return 1;
	}

	// Create OpenCV Mat using existing data (no copy)
	cv::Mat src(ry, rx, cv_type, fit->mask->data);
	cv::Mat dst;

	// Calculate kernel size from radius (should be odd)
	// Using approximation: kernel_size = 2 * ceil(3 * sigma) + 1, where sigma = radius
	int kernel_size = 2 * (int)ceilf(3.f * radius) + 1;

	// Apply Gaussian blur
	cv::GaussianBlur(src, dst, cv::Size(kernel_size, kernel_size), radius, radius);

	// Allocate new data and copy result
	size_t data_size = rx * ry * (fit->mask->bitpix / 8);
	void *blur_data = malloc(data_size);
	if (!blur_data) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	memcpy(blur_data, dst.data, data_size);

	// Replace original data with blurred data
	free(fit->mask->data);
	fit->mask->data = blur_data;

	return 0;
}

void set_poly_in_mask(UserPolygon *poly, fits *fit, gboolean state) {
	if (!poly || !fit || !fit->mask || !fit->mask->data) return;

	int width = fit->rx;
	int height = fit->ry;
	mask_t *m = fit->mask;

	// Prepare OpenCV points array
	std::vector<cv::Point> pts;
	pts.reserve(poly->n_points);
	for (int i = 0; i < poly->n_points; i++) {
		pts.push_back(cv::Point(
			round_to_int(poly->points[i].x),
			round_to_int(fit->ry - 1 - poly->points[i].y)
		));
	}

	// Define the fill value based on bitpix and state
	cv::Scalar color;
	if (state) {
		if (m->bitpix == 8) color = cv::Scalar(255);
		else if (m->bitpix == 16) color = cv::Scalar(65535);
		else color = cv::Scalar(1.0);
	} else {
		color = cv::Scalar(0);
	}

	// Create a Mat using the existing raw data (no copy)
	int type = (m->bitpix == 8) ? CV_8UC1 : (m->bitpix == 16 ? CV_16UC1 : CV_32FC1);
	cv::Mat mat(height, width, type, m->data);

	// Perform the polygon fill
	std::vector<std::vector<cv::Point>> contours = { pts };
	cv::fillPoly(mat, contours, color, cv::LINE_8);
}

int mask_feather(fits *fit, float feather_dist, feather_mode mode) {
	if (!fit || !fit->mask || !fit->mask->data) {
		siril_debug_print("mask_feather: invalid mask\n");
		return 1;
	}

	if (feather_dist <= 0.f) {
		siril_debug_print("mask_feather: feather_dist must be positive\n");
		return 1;
	}

	size_t rx = fit->rx, ry = fit->ry;
	size_t npixels = rx * ry;
	uint8_t bitpix = fit->mask->bitpix;

	// Convert mask to binary uint8 for OpenCV processing
	cv::Mat binary(ry, rx, CV_8UC1);
	uint8_t *binary_data = binary.data;

	// Convert mask to binary uint8
	switch (bitpix) {
		case 8: {
			uint8_t *m = (uint8_t*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				binary_data[i] = (m[i] > 127) ? 255 : 0;
			}
			break;
		}
		case 16: {
			uint16_t *m = (uint16_t*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				binary_data[i] = (m[i] > 32767) ? 255 : 0;
			}
			break;
		}
		case 32: {
			float *m = (float*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				binary_data[i] = (m[i] > 0.5f) ? 255 : 0;
			}
			break;
		}
	}

	// Determine which distance transforms we need
	gboolean need_inside = (mode == FEATHER_INNER || mode == FEATHER_EDGE);
	gboolean need_outside = (mode == FEATHER_OUTER || mode == FEATHER_EDGE);

	// Adjust feather distance for EDGE mode
	float effective_dist = (mode == FEATHER_EDGE) ? feather_dist / 2.0f : feather_dist;

	// Distance maps
	cv::Mat dist_inside, dist_outside;

	// Compute inside distance transform if needed
	if (need_inside) {
		cv::distanceTransform(binary, dist_inside, cv::DIST_L2, cv::DIST_MASK_PRECISE);
	}

	// Compute outside distance transform if needed
	if (need_outside) {
		// Invert binary mask for outside distance
		cv::Mat binary_inv;
		cv::bitwise_not(binary, binary_inv);
		cv::distanceTransform(binary_inv, dist_outside, cv::DIST_L2, cv::DIST_MASK_PRECISE);
	}

	// Apply feathering based on mode
	switch (bitpix) {
		case 8: {
			uint8_t *m = (uint8_t*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float value;
				if (mode == FEATHER_INNER) {
					if (dist_inside.at<float>(i) >= feather_dist) {
						value = 1.0f;
					} else {
						value = dist_inside.at<float>(i) / feather_dist;
					}
				} else if (mode == FEATHER_OUTER) {
					if (binary_data[i] > 127) {
						value = 1.0f;
					} else if (dist_outside.at<float>(i) >= feather_dist) {
						value = 0.0f;
					} else {
						value = 1.0f - (dist_outside.at<float>(i) / feather_dist);
					}
				} else { // FEATHER_EDGE
					if (binary_data[i] > 127) {
						// Inside: feather inward
						if (dist_inside.at<float>(i) >= effective_dist) {
							value = 1.0f;
						} else {
							value = 0.5f + (dist_inside.at<float>(i) / feather_dist);
						}
					} else {
						// Outside: feather outward
						if (dist_outside.at<float>(i) >= effective_dist) {
							value = 0.0f;
						} else {
							value = 0.5f - (dist_outside.at<float>(i) / feather_dist);
						}
					}
				}
				m[i] = (uint8_t)roundf(value * 255.0f);
			}
			break;
		}
		case 16: {
			uint16_t *m = (uint16_t*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float value;
				if (mode == FEATHER_INNER) {
					if (dist_inside.at<float>(i) >= feather_dist) {
						value = 1.0f;
					} else {
						value = dist_inside.at<float>(i) / feather_dist;
					}
				} else if (mode == FEATHER_OUTER) {
					if (binary_data[i] > 127) {
						value = 1.0f;
					} else if (dist_outside.at<float>(i) >= feather_dist) {
						value = 0.0f;
					} else {
						value = 1.0f - (dist_outside.at<float>(i) / feather_dist);
					}
				} else { // FEATHER_EDGE
					if (binary_data[i] > 127) {
						// Inside: feather inward
						if (dist_inside.at<float>(i) >= effective_dist) {
							value = 1.0f;
						} else {
							value = 0.5f + (dist_inside.at<float>(i) / feather_dist);
						}
					} else {
						// Outside: feather outward
						if (dist_outside.at<float>(i) >= effective_dist) {
							value = 0.0f;
						} else {
							value = 0.5f - (dist_outside.at<float>(i) / feather_dist);
						}
					}
				}
				m[i] = (uint16_t)roundf(value * 65535.0f);
			}
			break;
		}
		case 32: {
			float *m = (float*)fit->mask->data;
			for (size_t i = 0; i < npixels; i++) {
				float value;
				if (mode == FEATHER_INNER) {
					if (dist_inside.at<float>(i) >= feather_dist) {
						value = 1.0f;
					} else {
						value = dist_inside.at<float>(i) / feather_dist;
					}
				} else if (mode == FEATHER_OUTER) {
					if (binary_data[i] > 127) {
						value = 1.0f;
					} else if (dist_outside.at<float>(i) >= feather_dist) {
						value = 0.0f;
					} else {
						value = 1.0f - (dist_outside.at<float>(i) / feather_dist);
					}
				} else { // FEATHER_EDGE
					if (binary_data[i] > 127) {
						// Inside: feather inward
						if (dist_inside.at<float>(i) >= effective_dist) {
							value = 1.0f;
						} else {
							value = 0.5f + (dist_inside.at<float>(i) / feather_dist);
						}
					} else {
						// Outside: feather outward
						if (dist_outside.at<float>(i) >= effective_dist) {
							value = 0.0f;
						} else {
							value = 0.5f - (dist_outside.at<float>(i) / feather_dist);
						}
					}
				}
				m[i] = value;
			}
			break;
		}
	}

	const char *mode_str = (mode == FEATHER_INNER) ? "inward" :
						(mode == FEATHER_OUTER) ? "outward" : "on edge";
	siril_log_message(_("Mask feathered %s with distance %.1f pixels\n"), mode_str, feather_dist);
	return 0;
}

/**
 * mask_cleanup_morphological:
 * @fit: The fits image containing the mask to clean up
 * @close_size: Size of closing kernel to remove small holes (0 to skip)
 * @open_size: Size of opening kernel to remove small noise (0 to skip)
 * @denoise_threshold: Threshold for area-based denoising (0 to skip)
 *
 * Cleans up a noisy mask using morphological operations and connected component analysis.
 * This is particularly useful for color masks where hue is poorly defined in dark areas.
 *
 * Returns: 0 on success, -1 on error
 */
int mask_cleanup_morphological(fits *fit, int close_size, int open_size, int denoise_threshold) {
	if (!fit || !fit->mask) {
		siril_log_message(_("No mask to clean up\n"));
		return -1;
	}

	int width = fit->rx;
	int height = fit->ry;
	int ndata = width * height;

	// Convert mask to OpenCV Mat
	cv::Mat mask_mat(height, width, CV_8UC1);
	uint8_t *data8 = NULL;
	uint16_t *data16 = NULL;
	float *data32 = NULL;

	switch (fit->mask->bitpix) {
		case 8:  // 8-bit
			data8 = (uint8_t *)fit->mask->data;
			memcpy(mask_mat.data, data8, ndata);
			break;
		case 16:  // 16-bit
			data16 = (uint16_t *)fit->mask->data;
			for (int i = 0; i < ndata; i++) {
				mask_mat.data[i] = data16[i] >> 8;
			}
			break;
		case 32:  // 32-bit float
			data32 = (float *)fit->mask->data;
			for (int i = 0; i < ndata; i++) {
				mask_mat.data[i] = (uint8_t)(data32[i] * 255.0f);
			}
			break;
		default:
			siril_log_message(_("Unknown mask bitpix: %d\n"), fit->mask->bitpix);
			return -1;
	}

	cv::Mat processed = mask_mat.clone();

	// Step 1: Morphological closing - fills small holes
	if (close_size > 0) {
		cv::Mat element_close = cv::getStructuringElement(
			cv::MORPH_ELLIPSE,
			cv::Size(2 * close_size + 1, 2 * close_size + 1),
			cv::Point(close_size, close_size)
		);
		cv::morphologyEx(processed, processed, cv::MORPH_CLOSE, element_close);
	}

	// Step 2: Morphological opening - removes small noise
	if (open_size > 0) {
		cv::Mat element_open = cv::getStructuringElement(
			cv::MORPH_ELLIPSE,
			cv::Size(2 * open_size + 1, 2 * open_size + 1),
			cv::Point(open_size, open_size)
		);
		cv::morphologyEx(processed, processed, cv::MORPH_OPEN, element_open);
	}

	// Step 3: Connected component analysis to remove small isolated regions
	if (denoise_threshold > 0) {
		cv::Mat labels, stats, centroids;
		int num_labels = cv::connectedComponentsWithStats(
			processed, labels, stats, centroids, 8, CV_32S
		);

		// Create output mask - start with all zeros
		cv::Mat filtered = cv::Mat::zeros(height, width, CV_8UC1);

		// Keep only components larger than threshold
		// Label 0 is background, so start from 1
		for (int label = 1; label < num_labels; label++) {
			int area = stats.at<int>(label, cv::CC_STAT_AREA);

			if (area >= denoise_threshold) {
				// Keep this component
				cv::Mat component_mask = (labels == label);
				filtered.setTo(255, component_mask);
			}
		}

		processed = filtered;
	}

	// Convert back to original mask format
	switch (fit->mask->bitpix) {
		case 8:  // 8-bit
			data8 = (uint8_t *)fit->mask->data;
			memcpy(data8, processed.data, ndata);
			break;
		case 16:  // 16-bit
			data16 = (uint16_t *)fit->mask->data;
			for (int i = 0; i < ndata; i++) {
				data16[i] = ((uint16_t)processed.data[i]) << 8;
			}
			break;
		case 32:  // 32-bit float
			data32 = (float *)fit->mask->data;
			for (int i = 0; i < ndata; i++) {
				data32[i] = processed.data[i] / 255.0f;
			}
			break;
	}

	return 0;
}

// Create a mask based on image gradient magnitude.
// The gradient magnitude is normalized so the maximum gradient
// maps to the maximum mask value (255/65535/1.0).
// Return 0 on success.

FAST_MATH_PUSH
int mask_update_with_gradient(fits *fit) {
    if (!fit) return 1;
    if (!fit->mask || !fit->mask->data) {
        siril_log_message(_("Error: no existing mask to update\n"));
        return 1;
    }

    size_t rx = fit->rx, ry = fit->ry;
    size_t npixels = rx * ry;
    uint8_t bitpix = fit->mask->bitpix;

    if (!(bitpix == 8 || bitpix == 16 || bitpix == 32)) return 1;

    // Convert existing mask data to CV_32F for gradient computation
    cv::Mat mask_img(ry, rx, CV_32F);
    float *mask_data = mask_img.ptr<float>();

    // Copy mask data (handling different bit depths)
    switch (bitpix) {
        case 8: {
            uint8_t *src = (uint8_t*)fit->mask->data;
            for (size_t i = 0; i < npixels; i++) {
                mask_data[i] = (float)src[i] / 255.0f;
            }
            break;
        }
        case 16: {
            uint16_t *src = (uint16_t*)fit->mask->data;
            for (size_t i = 0; i < npixels; i++) {
                mask_data[i] = (float)src[i] / 65535.0f;
            }
            break;
        }
        case 32: {
            float *src = (float*)fit->mask->data;
            for (size_t i = 0; i < npixels; i++) {
                mask_data[i] = src[i];
            }
            break;
        }
    }

    // Compute gradients using Sobel operator
    cv::Mat grad_x, grad_y;
    cv::Sobel(mask_img, grad_x, CV_32F, 1, 0, 3);  // dx
    cv::Sobel(mask_img, grad_y, CV_32F, 0, 1, 3);  // dy

    // Compute gradient magnitude: sqrt(gx^2 + gy^2)
    cv::Mat grad_mag(ry, rx, CV_32F);
    float *gx_data = grad_x.ptr<float>();
    float *gy_data = grad_y.ptr<float>();
    float *mag_data = grad_mag.ptr<float>();

    float max_mag = 0.0f;
    for (size_t i = 0; i < npixels; i++) {
        mag_data[i] = sqrtf(gx_data[i] * gx_data[i] + gy_data[i] * gy_data[i]);
        if (mag_data[i] > max_mag) max_mag = mag_data[i];
    }

    // Avoid division by zero
    if (max_mag == 0.0f) {
        siril_log_message(_("Warning: gradient magnitude is zero everywhere\n"));
        max_mag = 1.0f;
    }

    // Normalize and update the existing mask
    switch (bitpix) {
        case 8: {
            uint8_t *m = (uint8_t*)fit->mask->data;
            for (size_t i = 0; i < npixels; i++) {
                float normalized = mag_data[i] / max_mag;
                m[i] = (uint8_t)roundf(normalized * 255.0f);
            }
            break;
        }
        case 16: {
            uint16_t *m = (uint16_t*)fit->mask->data;
            for (size_t i = 0; i < npixels; i++) {
                float normalized = mag_data[i] / max_mag;
                m[i] = (uint16_t)roundf(normalized * 65535.0f);
            }
            break;
        }
        case 32: {
            float *m = (float*)fit->mask->data;
            for (size_t i = 0; i < npixels; i++) {
                m[i] = mag_data[i] / max_mag;
            }
            break;
        }
    }

    set_mask_active(fit, TRUE);
    show_or_hide_mask_tab();
    siril_log_message(_("Mask updated with its gradient (max gradient: %.3f)\n"), max_mag);

    return 0;
}
FAST_MATH_POP
