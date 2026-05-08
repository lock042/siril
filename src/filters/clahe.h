#ifndef SRC_FILTERS_CLAHE_H_
#define SRC_FILTERS_CLAHE_H_

#include "core/siril.h"
#include "core/processing.h"

typedef struct {
	destructor destroy_fn;  /* Must be first member */
	double clip;
	int tileSize;
} clahe_params;

int clahe_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);
gchar *clahe_log_hook(gpointer p, log_hook_detail detail);
/* apply_clahe_cancel() is declared in gui/clahe.h */

#endif /* SRC_FILTERS_CLAHE_H_ */
