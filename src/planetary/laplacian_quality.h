#ifndef _LAPLACIAN_QUALITY_H_
#define _LAPLACIAN_QUALITY_H_

#include "core/siril.h"
#include "core/processing.h"

struct lapl_data {
	gboolean use_caching;
	int kernel_size;
	regdata *current_regdata;
	struct planetary_cache *cache;

	gboolean for_zones;
	int nb_zones;
	double *max;
};

int lapl_prepare_hook(struct generic_seq_args *args);
int lapl_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_);
int lapl_finalize_hook(struct generic_seq_args *args);

int laplace_quality_for_zones(sequence *seq, gboolean use_caching, int kernel_size);

#endif
