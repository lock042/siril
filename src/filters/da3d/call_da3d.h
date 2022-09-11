#ifndef SRC_CALL_DA3D_H_
#define SRC_CALL_DA3D_H_

//#ifdef __cplusplus
//#define EXTERNC extern "C" {
//#else
//#define EXTERNC
//#endif

/*EXTERNC*/ int call_da3d(float *in, float *gd, float *out, unsigned width, unsigned height, unsigned nchans, float sigma);

//#ifdef __cplusplus
//}
//#endif

#endif /* SRC_CALL_DA3D_H_ */
