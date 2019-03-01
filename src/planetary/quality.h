#ifndef _PLANETARY_QUALITY_H_
#define _PLANETARY_QUALITY_H_

#include "core/siril.h"	// for BYTE
#include "planetary.h"

int the_laplace_multipoint_quality_analysis(struct mpr_args *args);
int the_siril_multipoint_quality_analysis(struct mpr_args *args);

void create_index_of_best_zones(BYTE ***best_zones, int nb_best, int nb_images);
BYTE * sort_zones_quality(int zone_idx, int nb_best, int nb_seq_images);

#endif
