#ifndef SRC_FILTERS_STARNET_H_
#define SRC_FILTERS_STARNET_H_

#include "core/siril.h"
#include "core/processing.h"

typedef struct starnet_data {
	gboolean force_ser;
	fits *starnet_fit;
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
	struct multi_output_data *multi_args;
} starnet_data;

typedef struct remixargs {
	fits *fit1;
	fits *fit2;
} remixargs;

starnet_version starnet_executablecheck(gchar* executable);
gpointer do_starnet(gpointer p);
void apply_starnet_to_sequence(struct multi_output_data *multi_args);
void free_starnet_args(starnet_data *args);

#endif /* SRC_FILTERS_STARNET_H_ */
