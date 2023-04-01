#ifndef _GHT_H_
#define _GHT_H_

#define STRETCH_PAYNE_NORMAL 0
#define STRETCH_PAYNE_INVERSE 1
#define STRETCH_ASINH 2
#define STRETCH_INVASINH 3
#define STRETCH_LINEAR 4

#define COL_INDEP 0
#define COL_HUMANLUM 1
#define COL_EVENLUM 2
#define COL_SAT 3

#include <math.h>
#include "core/siril.h"

typedef struct ght_params {
	float B, D, LP, SP, HP, BP;
	int stretchtype, payne_colourstretchmodel;
	gboolean do_red, do_green, do_blue;
} ght_params;

typedef struct ght_compute_params {
	float qlp, q0, qwp, q1, q, b1, a1, a2, b2, c2, d2, e2, a3, b3, c3, d3, e3, a4, b4, LPT, SPT, HPT;
} ght_compute_params;

int GHTsetup(ght_compute_params* compute_params, float B, float D, float LP, float SP, float HP, int stretchtype);

float GHT(float in, float B, float D, float LP, float SP, float HP, float BP, int stretchtype, ght_compute_params *compute_params);
float GHTp(float in, ght_params *params, ght_compute_params *compute_params);
void apply_linked_ght_to_fits(fits *from, fits *to, ght_params *params, gboolean multithreaded);
void apply_sat_ght_to_fits(fits *from, fits *to, ght_params *params, gboolean multithreaded);
void apply_linked_ght_to_fbuf_indep(float* in, float* out, size_t layersize, size_t nchans, ght_params *params, gboolean multithreaded);

#endif
