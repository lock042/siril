#include "core/siril.h"
#include "opencv/kombat/kombat.h"

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

/* TODO: Duplicated, code stolen from opencv.cpp */
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

typedef struct {
	Mat ref_mat;
	Mat ref_mat_t;
	Mat result;
	int ready;
} kombat_cache;

int kombat_find_template(int idx, struct registration_args *args, fits *templ, fits *image,  reg_kombat *reg_param, void **vcache)
{
	int layer = args->layer;
	// Mat ref_mat;	
	kombat_cache *cache;
	Mat im;
	// Mat result;
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
		int results_w = image->rx - cache->ref_mat.cols + 1;
    	int results_h = image->ry - cache->ref_mat.rows + 1;
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