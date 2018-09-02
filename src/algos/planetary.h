#ifndef _PLANETARY_H
#define _PLANETARY_H

#include "core/siril.h"
#include "core/processing.h"

struct mpregdata {
	float x, y;
};

struct mpr_args {
	/* sequence configuration fields */
	sequence *seq;
	int layer;
	seq_image_filter filtering_criterion;
	double filtering_parameter;

	/* configuration for stacking */
	int nb_closest_AP;	// max number of closest AP to use
	double max_distance;	// discard AP farther than this, in pixels
	double own_distance_f;	// factor to the half-side of zone that gives
       				// the distance at which other zones are ignored
	
	char *output_filename;	// the file name for stacking result
	gboolean output_overwrite;

	/* internal use, set by registration for stacking */
	struct mpregdata **regdata;// regdata[image][zone]
	unsigned long *sum[3];	// the new image's channels
};

char *get_reference_image_name(sequence *seq, int layer);
int refimage_is_set();
const fits *get_refimage();
const char *get_refimage_filename();
void update_refimage_on_layer_change(sequence *seq, int layer);

gpointer the_multipoint_processing(gpointer ptr);

int point_is_inside_zone(int px, int py, stacking_zone *zone);

#endif
