#if !defined(MATCH_H)
#define MATCH_H

#include "core/siril.h"
#include "atpmatch.h"


#define NB_OF_MATCHING_TRY 3


int new_star_match(psf_star **s1, psf_star **s2, int n, int nobj_override,
		double s_min, double s_max, Homography *H, gboolean for_astrometry, gboolean save_photometric_data,
		pcc_star **photometric_data, int *nb_photometric_stars,
		transformation_type type, s_star **out_list_A, s_star **out_list_B);

#endif   /* MATCH_H */
