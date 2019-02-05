#ifndef _PLANETARY_H
#define _PLANETARY_H

#include "core/siril.h"
#include "core/processing.h"

/* This struct conveys data from the start of the multi-point registration to
 * the end of the stacking. */
struct mpr_args {
	/* sequence configuration fields */
	sequence *seq;
	int layer;
	seq_image_filter filtering_criterion;
	double filtering_parameter;	// the best quality
	double filtering_percent;	// the percentage of images to keep
	int using_homography;		// using translation or homography

	/* configuration for stacking */
	int nb_closest_AP;	// max number of closest AP to use
	double max_distance;	// discard AP farther than this, in pixels
	double own_distance_f;	// factor to the half-side of zone that gives
       				// the distance at which other zones are ignored
	
	char *output_filename;	// the file name for stacking result
	gboolean output_overwrite;

	double *global_image;	// the global image buffer as normalized double
};

gpointer sequence_analysis_thread_func(gpointer p);
char *get_reference_image_name(sequence *seq, int layer);
int refimage_is_set();
const fits *get_refimage();
const char *get_refimage_filename();
void update_refimage_on_layer_change(sequence *seq, int layer);

gpointer the_multipoint_analysis(gpointer ptr);
gpointer the_multipoint_processing(gpointer ptr);

int point_is_inside_zone(int px, int py, stacking_zone *zone);

#endif
