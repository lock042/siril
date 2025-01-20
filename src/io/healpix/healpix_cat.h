#ifndef HEALPIX_CATALOGUES_H
#define HEALPIX_CATALOGUES_H

#ifdef __cplusplus
extern "C" {
#endif

int get_raw_stars_from_local_gaia_astro_catalogue(double ra, double dec, double radius, double limitmag, gboolean phot, deepStarData **stars, uint32_t *nb_stars);

#ifdef __cplusplus
}
#endif

#endif
