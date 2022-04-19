#include "core/siril.h"
#include "opencv/kombat/kombat.h"
#include "opencv/opencv.h"

#include "registration/registration.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"

#include "gui/progress_and_log.h"
#include "core/siril_log.h"

#include <fitsio.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
/****************************************************************************************\
*                                       KOMBAT                     *
 \****************************************************************************************/

/* TODO: A local, KOMBAT specific version of image_to_Mat()
     (formerly from opencv.cpp). Versions might be merged? */
 int image_to_Mat(fits *image, Mat& in)
 {
	 if (image->naxes[2] != 1 && image->naxes[2] != 3) {
		return -1;
	}
	if (image->type == DATA_USHORT) {
		if (image->naxes[2] == 1) {
			in  = Mat(image->ry, image->rx, CV_16UC1, image->data);
		}
		else if (image->naxes[2] == 3) {
			WORD *bgr_u = fits_to_bgrbgr_ushort(image);
			if (!bgr_u) return -1;
			free(image->data);
			image->data = NULL;
			memset(image->pdata, 0, sizeof image->pdata);
			in = Mat(image->ry, image->rx, CV_16UC3, bgr_u);
		}
	}
	else if (image->type == DATA_FLOAT) {
		if (image->naxes[2] == 1) {
			in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
		}
		else if (image->naxes[2] == 3) {
			float *bgr_f = fits_to_bgrbgr_float(image);
			if (!bgr_f) return -1;
			free(image->fdata);
			image->fdata = NULL;
			memset(image->fpdata, 0, sizeof image->fpdata);
			in = Mat(image->ry, image->rx, CV_32FC3, bgr_f);
		}
	}
	else return -1;
	return 0;
}

void image_crop_params(reg_kombat *templ_in_ref, int sel_w, int sel_h, int full_w, int full_h, Rect& roi)
{
	/* we accept pixels to be moved of 25% of the image size;
	    smaller values give a great speed boost, but may prevent alignment if
		captures are noisy */
	float PC = 0.25f;

	roi.x = (int) (templ_in_ref->dx - full_w*PC/2);
	if (roi.x < 0) roi.x = 0;
	if (roi.x >= full_w) roi.x = full_w - 1;

	roi.width = sel_w + (int) (full_w*PC);
	if (roi.width+roi.x>=full_w) roi.width = full_w - roi.x - 1;

	roi.y = (int) (templ_in_ref->dy - full_h*PC/2);
	if (roi.y < 0) roi.y = 0;
	if (roi.y >= full_h) roi.y = full_h - 1;

	roi.height = sel_h + (int) (full_h*PC);
	if (roi.height+roi.y>=full_h) roi.height = full_h - roi.y - 1;
}

typedef struct {
	Mat ref_mat;
	Mat ref_mat_t;
	Mat result;
	Rect crop_rect;
	int ready;
} kombat_cache;

int kombat_find_template(int idx, struct registration_args *args, fits *templ, fits *image,  reg_kombat *reg_param, reg_kombat *ref_align, void **vcache)
{
	kombat_cache *cache;
	Mat im;
	int ret = 0;

	Mat im_t;
	Mat ref_mat_t;

	if (!vcache) {
		cache = new kombat_cache();
		cache->ready = 0;
	} else {
		if (!*vcache) {
			(*vcache) = (void*) new kombat_cache();
			cache = (kombat_cache*) *vcache;
			cache->ready = 0;
		} else
			cache = (kombat_cache*) *vcache;
	}

	if (!cache->ready) {
		ret = image_to_Mat(templ, cache->ref_mat);
		cache->ref_mat.convertTo(cache->ref_mat_t, CV_8U);
		int results_w;
    	int results_h;
		/* if we have spoted template in reference frame, we will crop current image
			around template's spot to analyse less data and speed up process */
		if (ref_align) {
			image_crop_params(ref_align,
											cache->ref_mat_t.cols, cache->ref_mat_t.rows,
											image->rx, image->ry, cache->crop_rect);
			results_w = cache->crop_rect.width - cache->ref_mat.cols + 1;
    	    results_h = cache->crop_rect.height - cache->ref_mat.rows + 1;
		} else {
			results_w =  image->rx - cache->ref_mat_t.cols + 1;
			results_h =  image->ry - cache->ref_mat_t.rows + 1;
		}
		cache->result.create( results_h, results_w, CV_8UC1 );
		cache->ready = 1;
	}

	if (ret) {
		if (!vcache) delete(cache);
		return ret;
	}

	ret = image_to_Mat(image, im);

	if (ret) {
		if (!vcache) delete(cache);
		return ret;
	}

	/* if we can crop, we do it. Only a small part of image will thus have to be analysed */
	if (ref_align)
		im = im(cache->crop_rect);

	// TODO: image_to_Mat should produce expected result...
	im.convertTo(im_t, CV_8U);

	matchTemplate( im_t,  cache->ref_mat_t, cache->result, TM_CCORR_NORMED );

    /* few variables to analyse result just computed */
    double ignored1;
    double score;
    Point ignored2;
    Point matchLoc;

    /* score would be the matching score [0 .. 1.0] of pattern in img
         matchLoc.x and matchLoc.y store coordinates of pattern within img. */
    minMaxLoc( cache->result, &ignored1, &score, &ignored2, &matchLoc, Mat() );

	reg_param->dx = matchLoc.x;
	reg_param->dy = matchLoc.y;

	if (ref_align) {
		reg_param->dx += cache->crop_rect.x;
		reg_param->dy += cache->crop_rect.y;
	}

	ret = (score<0.70);

	if (!vcache) delete(cache);
	return(ret);
}

void kombat_done(void **vcache)
{
	if (vcache) {
			kombat_cache *cache = (kombat_cache*)  (*vcache);
			delete(cache);
			(*vcache) = NULL;
	}
}
