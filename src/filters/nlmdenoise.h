#ifndef SRC_FILTERS_NLMDENOISE_H_
#define SRC_FILTERS_NLMDENOISE_H_

#include "core/siril.h"

typedef struct NLMDenoise_data {
	fits *image;
	float h_lum;
	float h_AB;
} NLMDenoise_data;

gpointer nlmdenoise(gpointer p);

#endif /* SRC_FILTERS_NLMDENOISE_H_ */
