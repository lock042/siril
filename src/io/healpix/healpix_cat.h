#ifndef HEALPIX_CATALOGUES_H
#define HEALPIX_CATALOGUES_H

#include "io/healpix/xp_continuous.h"  /* SourceEntryXPcts + xpcts_to_xpsampled */

#ifdef __cplusplus
extern "C" {
#endif

/* Photometric-catalogue kinds we can find in the user's gaia photo dir.
 * Values 2/3 mirror the cat_type byte at offset 50 of the on-disk header
 * (Siril HEALpix Catalog Format 1.0.0, https://zenodo.org/records/14697486).
 */
typedef enum {
    LOCAL_GAIA_PHOTO_NONE   =  0,  /* no chunks found / dir unreadable */
    LOCAL_GAIA_PHOTO_XPSAMP =  2,
    LOCAL_GAIA_PHOTO_XPCTS  =  3,
    LOCAL_GAIA_PHOTO_MIXED  = -1,  /* dir contains both kinds          */
    LOCAL_GAIA_PHOTO_BAD    = -2,  /* found chunks but headers bad     */
} local_gaia_photo_kind;

/* Probe the user's gaia photometric catalogue directory and report which
 * kind of catalogue it holds. Filename pattern narrows the candidates;
 * the cat_type byte in the first chunk's header confirms. */
local_gaia_photo_kind detect_local_gaia_photo_kind(const char *dir);

/* Backwards-compatible: returns non-zero iff *any* photometric chunk file
 * (xpsamp OR xpcts) is found in catalogue_paths[5]. */
int local_gaia_xpsamp_available();

/* New: distinct probes if a caller specifically wants one kind. */
int local_gaia_xpcts_available();

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

int get_raw_stars_from_remote_gaia_xpsampled_catalogue(double ra,
													  double dec,
													  double radius,
													  double limitmag,
													  SourceEntryXPsamp **stars,
													  uint32_t *nb_stars);

/* xp_continuous variants. The local one reads chunked xpcts catalogue files
 * from catalogue_paths[5]. The remote one fetches via HTTP RANGE from the
 * online xpcts catalogue mirror (if/when published). */
int get_raw_stars_from_local_gaia_xpcontinuous_catalogue(double ra,
														 double dec,
														 double radius,
														 double limitmag,
														 SourceEntryXPcts **stars,
														 uint32_t *nb_stars);

int get_raw_stars_from_remote_gaia_xpcontinuous_catalogue(double ra,
														  double dec,
														  double radius,
														  double limitmag,
														  SourceEntryXPcts **stars,
														  uint32_t *nb_stars);


#ifdef __cplusplus
}
#endif

#endif
