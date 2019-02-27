#ifndef _PLANETARY_REGISTRATION_H_
#define _PLANETARY_REGISTRATION_H_

/* This is the common header for registration files for planetary mode */

#include "planetary.h"

int the_multipoint_dft_registration(struct mpr_args *args);
int the_multipoint_ecc_registration(struct mpr_args *args);
int the_multipoint_steepest_descent_registration(struct mpr_args *args);

#endif
