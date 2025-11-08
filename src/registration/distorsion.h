#ifndef REGISTRATION_DISTORSION_H_
#define REGISTRATION_DISTORSION_H_

#include "core/siril.h"
#include "io/path_parse.h"

#define MAX_DISTO_SIZE 7 // need to duplicate MAX_DISTO_SIZE here because of circular refs with opencv

typedef enum {
	DISTO_NONE, // none defined
	DISTO_D2S,  // computed for each image dst->src (regular interpolation)
	DISTO_S2D,  // computed for each image src->dst (drizzle interpolation)
	DISTO_MAP_D2S,  // computed from the ref image dst->src (regular interpolation)
	DISTO_MAP_S2D  // computed from the ref image dst->src (drizzle interpolation)
} disto_type;

typedef struct {
	disto_type dtype;
	double A[MAX_DISTO_SIZE][MAX_DISTO_SIZE];
	double B[MAX_DISTO_SIZE][MAX_DISTO_SIZE];
	double AP[MAX_DISTO_SIZE][MAX_DISTO_SIZE];
	double BP[MAX_DISTO_SIZE][MAX_DISTO_SIZE];
	int order;
	double xref, yref;
	float *xmap, *ymap;
} disto_data;

int disto_correct_stars(psf_star **stars, disto_data *disto);
int init_disto_map(int rx, int ry, disto_data *disto);
void map_undistortion_D2S(disto_data *disto, int rx, int ry, float *xmap, float *ymap);
void map_undistortion_S2D(disto_data *disto, int rx, int ry, float *xmap, float *ymap);

gboolean validate_disto_params(fits *reffit, const gchar *text, disto_source index, gchar **msg1, gchar **msg2);
disto_data *init_disto_data(disto_params *distoparam, sequence *seq, struct wcsprm *WCSDATA, gboolean drizzle, int *status);
void free_disto_args(disto_data *disto);
void copy_disto(disto_data *disto_in, disto_data *disto_out);

#ifdef __cplusplus
extern "C" {
#endif
void prepare_H_with_disto_4remap(double *H, int rx_in, int ry_in, int rx_out, int ry_out, disto_data *disto, float *xmap, float *ymap);
#ifdef __cplusplus
}
#endif

#endif