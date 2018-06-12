#ifndef _IMAGE_DISPLAY_H
#define _IMAGE_DISPLAY_H

double get_zoom_val();
void adjust_vport_size_to_image();

void initialize_remap();

void redraw(int vport, int doremap);
void queue_redraw(int doremap); // request a redraw from another thread

#endif
