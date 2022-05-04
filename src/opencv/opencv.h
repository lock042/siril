#ifndef SIRIL_OPENCV_H_
#define SIRIL_OPENCV_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "registration/registration.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "gui/progress_and_log.h"

WORD *fits_to_bgrbgr_ushort(fits *image);
float *fits_to_bgrbgr_float(fits *image);

int cvResizeGaussian(fits *, int, int, int);

void cvResizeArray(double *, double *, int, int, int, int);

int cvRotateImage(fits *, point, double, int, int);

int cvAffineTransformation(fits *image, pointf *refpoints, pointf *curpoints, int nb_points,
		gboolean upscale2x, int interpolation, Homography *Hom);

unsigned char *cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *H, transformation_type type);


int cvTransformImage(fits *image, unsigned int width, unsigned int height, Homography Hom, gboolean upscale2x, int interpolation);

int cvUnsharpFilter(fits* image, double sigma, double amount);

int cvClahe(fits *image, double clip_limit, int size);

void cvRotateImageRefPoint(fits *image, point center, double angle, int cropped, point refpointin, point *refpointout);

void cvGetEye(Homography *H);

void cvTransfPoint(double *x, double *y, Homography Href, Homography Himg);

#ifdef __cplusplus
}
#endif

#endif
