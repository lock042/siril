#ifndef SIRIL_OPENCV_H_
#define SIRIL_OPENCV_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdint.h>
#include "registration/registration.h"
#include "registration/distorsion.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "gui/progress_and_log.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"

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


int cvTransformImage(fits *image, unsigned int width, unsigned int height, Homography Hom, float scale, int interpolation, gboolean clamp, disto_data *disto);

void cvDownscaleBlendMask(int rx, int ry, int out_rx, int out_ry, uint8_t *maskin, float *maskout);
void cvUpscaleBlendMask(int rx, int ry, int out_rx, int out_ry, float *maskin, float *maskout);

int cvUnsharpFilter(fits* image, double sigma, double amount);

int cvBilateralFilter(fits* image, double d, double sigma_col, double sigma_spatial);

int cvGuidedFilter(fits* image, fits *guide, double r, double eps);

int cvClahe(fits *image, double clip_limit, int size);

void cvTransformImageRefPoint(Homography Hom, point refpointin, point *refpointout);

void cvGetEye(Homography *H);

void cvTransfPoint(double *x, double *y, Homography Href, Homography Himg, double scale);

void cvTransfH(Homography *Href, Homography *Himg, Homography *Hres);

double cvCalculRigidTransform(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *Hom);

void cvMultH(Homography H1, Homography H2, Homography *Hout);
void cvInvertH(Homography *Hom);
void cvApplyFlips(Homography *Hom, int source_ry, int target_ry);
void cvPrepareDrizzleH(Homography *Hom, double scale, int source_rx, int source_ry, int target_rx, int target_ry);
void cvdisplay2ocv(Homography *Hom);

void cvGetMatrixReframe(double x, double y, int w, int h, double angle, Homography *Hom);
void cvGetBoundingRectSize(fits *image, point center, double angle, int *w, int *h);

// TODO: create and move to cvMosaic.h
gboolean cvRotMat3(double angles[3], rotation_type rottype[3], gboolean W2C, Homography *Hom);
void cvRelRot(Homography *Ref, Homography *R, Homography *Rout);
void cvcalcH_fromKKR(Homography *Kref, Homography *K, Homography *R, Homography *H);
int cvCalcH_from_corners(double *x_img, double *y_img, double *x_ref, double *y_ref, Homography *Hom);

int mask_apply_gaussian_blur(fits *fit, float radius);
void set_poly_in_mask(UserPolygon *poly, fits *fit, gboolean state);
int mask_feather(fits *fit, float feather_dist, feather_mode mode);

#ifdef __cplusplus
}
#endif

#endif
