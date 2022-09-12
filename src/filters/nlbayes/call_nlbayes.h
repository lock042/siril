#ifndef SRC_CALL_NLBAYES_H_
#define SRC_CALL_NLBAYES_H_

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif

EXTERNC int do_nlbayes(fits *fit, float modulation, int da3d);

#ifdef __cplusplus
}
#endif

#endif /* SRC_CALL_NLBAYES_H_ */
