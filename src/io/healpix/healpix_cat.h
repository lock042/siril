#ifndef HEALPIX_CATALOGUES_H
#define HEALPIX_CATALOGUES_H

#ifdef __cplusplus
extern "C" {
#endif

int local_gaia_xpsamp_available();
int get_raw_stars_from_local_gaia_astro_catalogue(double ra,
												  double dec,
												  double radius,
												  double limitmag,
												  gboolean phot,
												  deepStarData **stars,
												  uint32_t *nb_stars);

int get_raw_stars_from_local_gaia_xpsampled_catalogue(double ra,
													  double dec,
													  double radius,
													  double limitmag,
													  SourceEntryXPsamp **stars,
													  uint32_t *nb_stars);

#ifdef __cplusplus
}
#endif

#endif
