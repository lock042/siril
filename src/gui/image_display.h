#ifndef _IMAGE_DISPLAY_H_
#define _IMAGE_DISPLAY_H_

#include "core/siril.h"

typedef enum {
	REDRAW_OVERLAY, // only overlays changed
	REDRAW_IMAGE,   // the image changed, render it and overlays
	REMAP_ALL       // the image data changed, remap and render all
} remap_type;

void initialize_image_display();

void redraw(remap_type doremap);	// redraw the image, possibly with a remap
void queue_redraw(remap_type doremap); // call redraw from another thread

double get_zoom_val();	// for image_interactions

point get_center_of_vport();
void add_image_and_label_to_cairo(cairo_t *cr, int vport);

#endif

