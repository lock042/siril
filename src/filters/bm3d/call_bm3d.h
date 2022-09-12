#ifndef SRC_CALL_BM3D_H_
#define SRC_CALL_BM3D_H_

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif

typedef struct denoise_args {
	fits *fit;
	float modulation;
	int da3d;
} denoise_args;

void bgrbgr_float_to_fits(fits *image, float *bgrbgr, float modulation);
void bgrbgr_float_to_word_fits(fits *image, float *bgrbgr, float modulation);
float *fits_to_bgrbgr_wordtofloat(fits *image);

EXTERNC int do_bm3d(fits *fit, float modulation, int da3d);
#ifdef __cplusplus
}
#endif

#endif /* SRC_CALL_BM3D_H_ */
