#ifndef _TRACKING_H
#define _TRACKING_H

#include "algos/PSF.h"

struct linetrack_conf {
	fits *fit;
	int layer;
	float threshold;	// pixel value in image
	int minlen;
	psf_star **targets;
	int nb_targets;
};

gpointer tracking_worker(gpointer ptr);

#endif

