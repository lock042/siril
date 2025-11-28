#ifndef SRC_FILTERS_SCNR_H_
#define SRC_FILTERS_SCNR_H_

#include "core/siril.h"

typedef enum {
	SCNR_AVERAGE_NEUTRAL,
	SCNR_MAXIMUM_NEUTRAL,
	SCNR_MAXIMUM_MASK,
	SCNR_ADDITIVE_MASK
} scnr_type;

/* scnr data from GUI */
struct scnr_data {
	destructor destroy_fn;  // Must be first member for generic_image_worker
	scnr_type type;
	double amount;
	gboolean preserve;
	gboolean verbose;
	gboolean applying;  // TRUE for final apply, FALSE for preview
};

/* Allocator and destructor functions */
struct scnr_data *new_scnr_data();
void free_scnr_data(void *args);

/* Image processing hook for generic_image_worker */
int scnr_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

/* Idle functions */
gboolean scnr_preview_idle(gpointer p);
gboolean scnr_apply_idle(gpointer p);

const char *scnr_type_to_string(scnr_type t);

#endif /* SRC_FILTERS_SCNR_H_ */
