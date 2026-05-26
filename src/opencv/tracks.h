#ifndef _OPENCV_TRACKS_H_
#define _OPENCV_TRACKS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "core/siril.h"

struct track {
	pointi start;	// display coordinates
	pointi end;
	float angle;	// with the X axis

	float mag;	// 99.9 means undefined
	float mag_err;
	gboolean mag_is_absolute;
	gboolean mag_is_accurate;
	float snr;
};

int cvHoughLines(fits *image, int layer, float threshvalue, int minlen, struct track **tracks);

int create_streak_masks(int rx, int ry, int sx, int sy, int ex, int ey, float fwhm, BYTE **aperture_mask, BYTE **background_mask);

#ifdef __cplusplus
}
#endif

#endif
