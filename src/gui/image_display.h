#ifndef _IMAGE_DISPLAY_H_
#define _IMAGE_DISPLAY_H_

#include "core/siril.h"

typedef enum {
	REDRAW_OVERLAY, // only overlays changed
	REDRAW_IMAGE,   // the image changed, render it and overlays
	REMAP_ALL       // the image data changed, remap and render all
} remap_type;

void check_gfit_profile_identical_to_monitor();

void allocate_hd_remap_indices();
void hd_remap_indices_cleanup();

void initialize_image_display();

void copy_roi_into_gfit();

void redraw(remap_type doremap);	// redraw the image, possibly with a remap
void queue_redraw(remap_type doremap); // call redraw from another thread
void queue_redraw_and_wait_for_it(remap_type doremap); // call redraw from another thread and wait for it

double get_zoom_val();	// for image_interactions

point get_center_of_vport();
void add_image_and_label_to_cairo(cairo_t *cr, int vport);

gboolean get_context_rotation_matrix(double rotation, cairo_matrix_t *transform, gboolean invert); //computes rotation matrix about center of com.selection

#endif

