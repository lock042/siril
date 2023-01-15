// Helpful wrapper functions to call Siril functions without causing clashes with some of the img_expr functionality

#ifndef CPLUSPLUS_HELPER_FUNCTIONS_H
#define CPLUSPLUS_HELPER_FUNCTIONS_H

#ifdef __cplusplus
#define EXTERNC1 extern "C"
#else
#define EXTERNC1
#endif


EXTERNC1 void updateprogress(const char *text, double percent);
EXTERNC1 void sirillog(const char* text);
EXTERNC1 int is_thread_stopped();
EXTERNC1 int updatenoise(float *array, int nx, int ny, int nchans, double *noise);

#ifdef __cplusplus
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN int cppmaxthreads;
EXTERN unsigned cppfftwflags;
EXTERN double cppfftwtimelimit;
EXTERN int cppfftwmultithreaded;

#endif // CPLUSPLUS_HELPER_FUNCTIONS_H
