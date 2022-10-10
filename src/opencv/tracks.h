#ifndef _OPENCV_TRACKS_H_
#define _OPENCV_TRACKS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "core/siril.h"

struct track {
	pointi start;
	pointi end;
	float angle;	// with the X axis
};

int cvHoughLines(fits *image, int layer, float threshvalue, int minlen, struct track **tracks);


#ifdef __cplusplus
}
#endif

#endif
