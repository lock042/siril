#ifndef SRC_CALL_BM3D_H_
#define SRC_CALL_BM3D_H_

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif

EXTERNC int do_bm3d(fits *fit, float modulation, int da3d);
#ifdef __cplusplus
}
#endif

#endif /* SRC_CALL_BM3D_H_ */
