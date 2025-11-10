#ifndef SRC_FILTERS_BANDING_H_
#define SRC_FILTERS_BANDING_H_

#include "core/siril.h"

/* Banding data for both single image and sequence processing */
struct banding_data {
	destructor destroy_fn;  // Must be first member
	double sigma;
	double amount;
	gboolean protect_highlights;
	gboolean applyRotation;
	char *seqEntry;
	sequence *seq;
	fits *fit;
};

/* Allocator and destructor functions */
struct banding_data *new_banding_data();
void free_banding_data(void *args);

/* Image processing hooks */
int banding_single_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

void apply_banding_to_sequence(struct banding_data *banding_args);

#endif /* SRC_FILTERS_BANDING_H_ */
