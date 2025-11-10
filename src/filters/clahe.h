#ifndef SRC_GUI_CLAHE_H_
#define SRC_GUI_CLAHE_H_

typedef struct {
	destructor destroy_fn;  // Must be first member
	double clip;
	int tileSize;
} clahe_params;

int clahe_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
void apply_clahe_cancel();

#endif /* SRC_GUI_CLAHE_H_ */
