#ifndef SIRIL_OPENCV_H_
#define SIRIL_OPENCV_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdint.h>
#include "registration/registration.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "gui/progress_and_log.h"

#ifdef __cplusplus
extern "C" {
#endif

WORD *fits_to_bgrbgr_ushort(fits *image);
float *fits_to_bgrbgr_float(fits *image);

int cvResizeGaussian(fits *, int, int, int, gboolean);

void cvResizeArray(double *, double *, int, int, int, int);

int cvRotateImage(fits *image, int angle); // only for fast rotations

unsigned char *cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *H, transformation_type type, float offset);


int cvTransformImage(fits *image, unsigned int width, unsigned int height, Homography Hom, gboolean upscale2x, int interpolation, gboolean clamp);

int cvUnsharpFilter(fits* image, double sigma, double amount);

int cvClahe(fits *image, double clip_limit, int size);

void cvTransformImageRefPoint(Homography Hom, point refpointin, point *refpointout);

void cvGetEye(Homography *H);

void cvTransfPoint(double *x, double *y, Homography Href, Homography Himg);

void cvTransfH(Homography Href, Homography Himg, Homography *Hres);

double cvCalculRigidTransform(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *Hom);

void cvMultH(Homography H1, Homography H2, Homography *Hout);
void cvInvertH(Homography *Hom);
void cvApplyFlips(Homography *Hom, int source_ry, int target_ry);

void cvGetMatrixReframe(double x, double y, int w, int h, double angle, Homography *Hom);
void cvGetMatrixResize(double cxin, double cyin, double cxout, double cyout, double scale, Homography *Hom);
void cvGetBoundingRectSize(fits *image, point center, double angle, int *w, int *h);

#ifdef __cplusplus
}
#endif

#endif
