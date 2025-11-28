#ifndef _COLORS_H_
#define _COLORS_H_
#include <math.h>

typedef enum {
	EXTRACT_RGB,
	EXTRACT_HSL,
	EXTRACT_HSV,
	EXTRACT_CIELAB,
} channel_extract_type;

struct extract_channels_data {
	fits *fit;
	char *channel[3];
	channel_extract_type type;
	const char* str_type;
};

typedef float ccm[3][3]; // Color Conversion Matrix

struct ccm_data {
	destructor destroy_fn;
	ccm matrix;
	float power;
	fits *fit;
	sequence *seq;
	char *seqEntry;
};

void rgb_to_hsl_float_sat(float, float, float, float, float *, float *, float *);
void hsl_to_rgb_float_sat(float, float, float, float *, float *, float *);
void rgb_to_hsl(double, double, double, double *, double *, double *);
void hsl_to_rgb(double, double, double, double *, double *, double *);
void rgb_to_hslf(float r, float g, float b, float *h, float *s, float *l);
void hsl_to_rgbf(float h, float s, float l, float *r, float *g, float *b);
void rgbw_to_hslw(uint16_t r, uint16_t g, uint16_t b, uint16_t *h, uint16_t *s, uint16_t *l);
void hslw_to_rgbw(uint16_t h, uint16_t s, uint16_t l, uint16_t *r, uint16_t *g, uint16_t *b);
void rgb_to_hsv(double, double, double, double *, double *, double *);
void rgb_to_hsvf(float, float, float, float *, float *, float *);
void hsv_to_rgb(double, double, double, double *, double *, double *);
void rgb_to_xyz(double, double, double, double *, double *, double *);
void rgb_to_xyzf(float r, float g, float b, float *x, float *y, float *z);
void xyz_to_LAB(double, double, double, double *, double *, double *);
void xyz_to_LABf(float x, float y, float z, float *L, float *a, float *b);
void LAB_to_xyz(double, double, double, double *, double *, double *);
void LAB_to_xyzf(float L, float a, float b, float *x, float *y, float *z);
void xyz_to_rgb(double, double, double, double *, double *, double *);
void xyz_to_rgbf(float x, float y, float z, float *r, float *g, float *b);
void linrgb_to_xyz(double r, double g, double b, double *x, double *y, double *z, gboolean scale);
void xyz_to_linrgb(double x, double y, double z, double *r, double *g, double *b, gboolean scale);
void linrgb_to_xyzf(float r, float g, float b, float *x, float *y, float *z, gboolean scale);
void xyz_to_linrgbf(float x, float y, float z, float *r, float *g, float *b, gboolean scale);
void rgb_to_yuvf(float red, float green, float blue, float *y, float *u, float *v);
void yuv_to_rgbf(float y, float u, float v, float *red, float *green, float *blue);

double BV_to_T(double BV);

float x1931(float w);
float y1931(float w);
float z1931(float w);
float x1964(float w);
float y1964(float w);
float z1964(float w);
cmsCIExyY wl_to_xyY(double wl);

int pos_to_neg(fits *fit);

int equalize_cfa_fit_with_coeffs(fits *fit, float coeff1, float coeff2, const char *cfa_string);

gpointer extract_channels(gpointer p);

void background_neutralize(fits* fit, rectangle black_selection);
void get_coeff_for_wb(fits *fit, rectangle white, rectangle black,
		double kw[], double bg[], double norm, double low, double high);
int calibrate(fits *fit, int layer, double kw, double bg, double norm);
int ccm_calc(fits *fit, ccm matrix, float power);
void apply_ccm_to_sequence(struct ccm_data *ccm_args);
void free_ccm_data(void *ptr);
struct ccm_data *new_ccm_data();
int ccm_process_with_worker(ccm matrix, float power);
int ccm_single_image_hook(struct generic_img_args *args, fits *fit, int nb_threads);

#endif
