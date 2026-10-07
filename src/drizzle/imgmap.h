/*
 * Minimal pixel-map type used by both the STScI drizzle backend and
 * mpp_drizzle. Extracted from cdrizzleutil.h so it can be included from
 * C++ — the rest of cdrizzleutil.h pulls in driz_portability.h, which
 * does `#define private` for legacy C compatibility, and that wrecks
 * any subsequent inclusion of <string>, <unique_ptr>, etc.
 */
#ifndef SRC_DRIZZLE_IMGMAP_H_
#define SRC_DRIZZLE_IMGMAP_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _imgmap_t {
    float* xmap;
    float* ymap;
    int rx;
    int ry;
} imgmap_t;

#ifdef __cplusplus
}
#endif

#endif /* SRC_DRIZZLE_IMGMAP_H_ */
