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

char *get_reference_image_name(sequence *seq, int layer);
int refimage_is_set();
const fits *get_refimage();
const char *get_refimage_filename();
void update_refimage_on_layer_change(sequence *seq, int layer);

gpointer the_multipoint_registration(gpointer p);
int point_is_inside_zone(int px, int py, stacking_zone *zone);

#endif
