#ifndef _PCACHING_H_
#define _PCACHING_H_

#include "core/siril.h"

struct planetary_cache {
	int kernel_size;

	gboolean use_cached;	// read pre-computed data from cache
	// sequences used to load data from cache (seq using SER)
	sequence *seq_gaussian;
	sequence *seq_laplacian;

	gboolean cache_data;	// cache computed data for later reuse
	// sequences used to store data to cache (SER)
	struct ser_struct *ser_gaussian;
	struct ser_struct *ser_laplacian;
};

void init_caching(const char *seqname, struct planetary_cache *args, int kernel_size);
void finalize_caching(struct planetary_cache *args);

WORD * get_gaussian_data_for_image(int index, fits *fit, struct planetary_cache *args);
WORD * get_laplacian_data_for_image(int index, WORD *gaussian_data,
		int width, int height, struct planetary_cache *args);

WORD * get_gaussian_data_for_image_in_seq(sequence *seq, int index, struct planetary_cache *args);
#endif
