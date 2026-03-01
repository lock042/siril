#ifndef SRC_GUI_ASINH_H_
#define SRC_GUI_ASINH_H_

/* Structure to hold asinh-specific parameters */
typedef struct {
	destructor destroy_fn;  // generic deallocator
	float beta;
	float offset;
	gboolean human_luminance;
	clip_mode_t clip_mode;
} asinh_params;

gchar *asinh_log_hook(gpointer p, log_hook_detail detail);
int asinh_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
int asinhlut(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clip_mode);
void asinh_change_between_roi_and_image();
void apply_asinh_cancel();

#endif /* SRC_GUI_ASINH_H_ */
