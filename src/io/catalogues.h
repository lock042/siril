#ifndef _CATALOGUES_H
#define _CATALOGUES_H

#include "core/siril.h"

int get_stars_from_local_nomad(double ra, double dec, double radius, fits *fit, float max_mag, pcc_star **stars, int *nb_stars);

#endif
