#ifndef SRC_CALL_NLBAYES_H_
#define SRC_CALL_NLBAYES_H_

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif

typedef struct denoise_args {
	fits *fit;
	float modulation;
	int da3d;
	int sos;
	float rho;
} denoise_args;

void bgrbgr_float_to_fits(fits *image, float *bgrbgr, float modulation);
void bgrbgr_float_to_word_fits(fits *image, float *bgrbgr, float modulation);
float *fits_to_bgrbgr_wordtofloat(fits *image);

EXTERNC int do_nlbayes(fits *fit, float modulation, unsigned sos, int da3d, float rho);
#ifdef __cplusplus
}
#endif

#endif /* SRC_CALL_NLBAYES_H_ */
