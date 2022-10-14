#ifndef _TRACKING_H
#define _TRACKING_H

#include "algos/PSF.h"

struct linetrack_conf {
	fits *fit;
	int layer;
	float threshold;	// pixel value in image
	int minlen;
	psf_star **fixed_targets;
	int nb_fixed_targets;
	gboolean display_lines;
	gboolean use_idle;
};

gpointer tracking_worker(gpointer ptr);

#endif

