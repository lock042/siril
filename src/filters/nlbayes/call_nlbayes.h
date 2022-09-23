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
	gboolean do_anscombe;
} denoise_args;

EXTERNC int do_nlbayes(fits *fit, const float modulation, unsigned sos, int da3d, const float rho, const gboolean do_anscombe);
#ifdef __cplusplus
}
#endif

#endif /* SRC_CALL_NLBAYES_H_ */
