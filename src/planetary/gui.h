#ifndef _PLANETARY_CALLBACKS_H
#define _PLANETARY_CALLBACKS_H

void update_zones_list();
gboolean end_mpp_analysis(gpointer arg);

void planetary_click_in_image(double x, double y);

void display_refimage_if_needed();

void activate_mpp_processing_button();
#endif

