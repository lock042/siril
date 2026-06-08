#ifndef SRC_GUI_MPP_SHIFT_VIEWER_H_
#define SRC_GUI_MPP_SHIFT_VIEWER_H_

#include <gtk/gtk.h>

/* True if the shift-viewer dialog is currently visible. The image
 * overlay (draw_mpp_aps in image_display.c) checks this to decide
 * whether to paint per-AP shift arrows on top of the AP boxes. */
gboolean mpp_shift_viewer_is_open(void);

/* Zero-based index of the frame whose shifts to display, or -1 if the
 * viewer isn't usable (no shifts cached). The overlay reads this each
 * paint. */
int mpp_shift_viewer_get_frame(void);

/* Pixel multiplier applied to the (dy, dx) shift values when drawing
 * arrows — per-AP shifts are typically sub-pixel and need scaling to
 * be visible. The overlay reads this each paint. */
double mpp_shift_viewer_get_scale(void);

/* Refresh the button's sensitivity (called from end_register_idle and
 * close_sequence_idle, like the AP editor's button). Sensitive iff a
 * cached run with populated shifts exists. */
void mpp_shift_viewer_update_button_sensitivity(void);

#endif /* SRC_GUI_MPP_SHIFT_VIEWER_H_ */
