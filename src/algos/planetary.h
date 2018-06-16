#ifndef _PLANETARY_H
#define _PLANETARY_H

#include "core/siril.h"
#include "core/processing.h"

struct mpr_args {
	sequence *seq;
	int layer;
	seq_image_filter filtering_criterion;
	double filtering_parameter;
};

gpointer the_multipoint_registration(gpointer p);

#endif
