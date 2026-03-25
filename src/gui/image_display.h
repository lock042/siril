#ifndef _IMAGE_DISPLAY_H_
#define _IMAGE_DISPLAY_H_

#include "core/siril.h"

typedef enum {
	REDRAW_OVERLAY, // only overlays changed
	REDRAW_IMAGE,   // the image changed, render it and overlays
	REMAP_ALL       // the image data changed, remap and render all
} remap_type;

/* FLIS composite cache control — defined in image_display.c.
 * flis_invalidate_composite() must be called whenever any layer property
 * changes.  flis_composite_free() must be called when a FLIS file is
 * closed (done automatically by free_image_data() in single_image.c). */
void flis_invalidate_composite(void);
void flis_composite_free(void);

/* flis_render_layers: composites a GSList of flis_layer_t* into a newly
 * allocated float-RGB fits*.  The caller must clearfits()+free() the result.
 * Used by the display pipeline, merge-down, flatten, and flat-FITS save. */
fits *flis_render_layers(GSList *layers);

/* Controls whether the mask viewport and tint show the active FLIS layer's
 * layer mask (TRUE) or its processing mask (FALSE, default). */
gboolean get_flis_show_layer_mask(void);
void     set_flis_show_layer_mask(gboolean val);

void check_gfit_profile_identical_to_monitor();

void allocate_hd_remap_indices();
void hd_remap_indices_cleanup();

void initialize_image_display();

void copy_roi_into_gfit();

void redraw(remap_type doremap);	// redraw the image, possibly with a remap
void queue_redraw(remap_type doremap); // call redraw from another thread
void queue_redraw_and_wait_for_it(remap_type doremap); // call redraw from another thread and wait for it
gboolean redraw_mask_idle(gpointer p);
void queue_redraw_mask(); // queue a redraw of the mask only

double get_zoom_val();	// for image_interactions

point get_center_of_vport();
void add_image_and_label_to_cairo(cairo_t *cr, int vport);

gboolean get_context_rotation_matrix(double rotation, cairo_matrix_t *transform, gboolean invert); //computes rotation matrix about center of com.selection

#endif

