#ifndef SRC_FILTERS_STARNET_H_
#define SRC_FILTERS_STARNET_H_

#include "core/siril.h"
#include "core/processing.h"

typedef struct starnet_data {
	destructor destroy_fn;  // Must be first member
	gboolean force_ser;
	fits *starnet_fit;  // reference only, not freed
	fits *starmask_fit;
	struct timeval t_start;
	gchar *stride;
	gboolean linear;
	gboolean customstride;
	gboolean upscale;
	gboolean starmask;
	gboolean follow_on;
	gboolean too_small;
	int imgnumber;
	struct multi_output_data *multi_args;  // reference only, not freed
} starnet_data;

typedef struct remixargs {
	fits *fit1;
	fits *fit2;
} remixargs;

/* Allocator and destructor */
starnet_data *new_starnet_args();
void free_starnet_args(void *ptr);

/* Image processing hook for single images */
int starnet_single_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

/* Idle function for single images */
gboolean starnet_single_image_idle(gpointer p);

starnet_version starnet_executablecheck(gchar* executable);
gpointer do_starnet(gpointer p);  // Legacy function for backward compatibility
void apply_starnet_to_sequence(struct multi_output_data *multi_args);

#endif /* SRC_FILTERS_STARNET_H_ */
