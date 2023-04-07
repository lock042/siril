#ifndef _COLORS_H_
#define _COLORS_H_

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

void rgb_to_hsl_float_sat(float, float, float, float, float *, float *, float *);
void hsl_to_rgb_float_sat(float, float, float, float *, float *, float *);
void rgb_to_hsl(double, double, double, double *, double *, double *);
void rgb_to_hslf(float r, float g, float b, float *h, float *s, float *l);
void hsl_to_rgb(double, double, double, double *, double *, double *);
void hsl_to_rgbf(float h, float s, float l, float * r, float * g, float * b);
void rgb_to_hsv(double, double, double, double *, double *, double *);
void hsv_to_rgb(double, double, double, double *, double *, double *);
void rgb_to_xyz(double, double, double, double *, double *, double *);
void rgb_to_xyzf(float r, float g, float b, float *x, float *y, float *z);
void xyz_to_LAB(double, double, double, double *, double *, double *);
void xyz_to_LABf(float x, float y, float z, float *L, float *a, float *b);
void LAB_to_xyz(double, double, double, double *, double *, double *);
void LAB_to_xyzf(float L, float a, float b, float *x, float *y, float *z);
void xyz_to_rgb(double, double, double, double *, double *, double *);
void xyz_to_rgbf(float x, float y, float z, float *r, float *g, float *b);
double BV_to_T(double BV);

int pos_to_neg(fits *fit);

int equalize_cfa_fit_with_coeffs(fits *fit, float coeff1, float coeff2, const char *cfa_string);

gpointer extract_channels(gpointer p);

void background_neutralize(fits* fit, rectangle black_selection);
void get_coeff_for_wb(fits *fit, rectangle white, rectangle black,
		double kw[], double bg[], double norm, double low, double high);
int calibrate(fits *fit, int layer, double kw, double bg, double norm);


#endif
