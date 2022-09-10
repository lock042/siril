#ifndef SRC_CALL_BM3D_H_
#define SRC_CALL_BM3D_H_

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif

typedef struct bm3d_args {
	fits *fit;
	float modulation;
} bm3d_args;

EXTERNC int do_bm3d(fits *fit, float modulation);

EXTERNC gpointer run_bm3d_on_fit(gpointer p);

#ifdef __cplusplus
}
#endif

#endif /* SRC_CALL_BM3D_H_ */
