#ifndef _IMAGE_DISPLAY_H_
#define _IMAGE_DISPLAY_H_

#include "core/siril.h"

void initialize_image_display();

void queue_redraw(int doremap);
void redraw(int vport, int remap);
double get_zoom_val();	// for image_interactions
point get_center_of_vport();
void add_label_to_cairo(GtkWidget *widget, cairo_t *cr);

#endif

