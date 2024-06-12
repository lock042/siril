#ifndef SRC_GUI_ASINH_H_
#define SRC_GUI_ASINH_H_

int command_asinh(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clipmode);
int asinhlut(fits *fit, float beta, float offset, gboolean rgb_space);
void asinh_change_between_roi_and_image();
void apply_asinh_cancel();

#endif /* SRC_GUI_ASINH_H_ */
