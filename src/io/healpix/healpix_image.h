#ifndef SRC_IO_HEALPIX_HEALPIX_IMAGE_H_
#define SRC_IO_HEALPIX_HEALPIX_IMAGE_H_

#ifdef __cplusplus
extern "C" {
#endif

struct ffit;
typedef struct ffit fits;

/* Reports, via siril_log_message, the level 1 (Nside=2) and level 8 (Nside=256)
 * NESTED HEALPix pixels that overlap with the footprint of the plate-solved
 * image. Level 8 pixels are grouped under their level 1 parent.
 *
 * Returns 0 on success, non-zero on failure. The fit must be plate-solved. */
int log_image_healpixels(fits *fit);

#ifdef __cplusplus
}
#endif

#endif /* SRC_IO_HEALPIX_HEALPIX_IMAGE_H_ */
