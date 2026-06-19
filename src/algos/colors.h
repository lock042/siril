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

/* ----------------------------------------------------------------------------
 * Branchless double HSL <-> RGB (h,s,l on [0,1]). Same exclusive-mask hue
 * dispatch and trailing selects as the float cores below; matches the original
 * branchy double implementation (l clamped to 0 for non-positive luminance).
 * ------------------------------------------------------------------------- */
static inline SIRIL_VECTORIZABLE
void rgb_to_hsl(double r, double g, double b, double *h, double *s, double *l) {
	double v = (r > g) ? r : g; v = (v > b) ? v : b;
	double m = (r < g) ? r : g; m = (m < b) ? m : b;
	double vm = v - m;
	double L  = (m + v) / 2.0;

	double inv = 1.0 / vm;                  /* +inf when vm==0; discarded below */
	double r2 = (v - r) * inv, g2 = (v - g) * inv, b2 = (v - b) * inv;
	double sden = (L <= 0.5) ? (v + m) : (2.0 - v - m);
	double S = vm / sden;

	int mr = (r == v);
	int mg = (g == v) & (mr ^ 1);
	int mb = (mr | mg) ^ 1;
	double H = mr * (g == m ? 5.0 + b2 : 1.0 - g2)
	         + mg * (b == m ? 1.0 + r2 : 3.0 - b2)
	         + mb * (r == m ? 3.0 + g2 : 5.0 - r2);
	H /= 6.0;

	int keep = (L > 0.0) & (vm > 0.0);
	*l = (L > 0.0) ? L : 0.0;
	*s = keep ? S : 0.0;
	*h = keep ? H : 0.0;
}

static inline SIRIL_VECTORIZABLE
void hsl_to_rgb(double h, double sl, double l, double *r, double *g, double *b) {
	h -= floor(h);                          /* range-reduce to [0,1) */
	double v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);
	double m = l + l - v;
	double sv = (v - m) / v;                /* v<=0 lane discarded by select below */
	double h6 = h * 6.0;
	int sextant = (int) h6;                 /* 0..5 */
	double fract = h6 - (double) sextant;
	double vsf = v * sv * fract;
	double mid1 = m + vsf, mid2 = v - vsf;

	double R = ((sextant==0)|(sextant==5)) ? v : ((sextant==2)|(sextant==3)) ? m : (sextant==4) ? mid1 : mid2;
	double G = ((sextant==1)|(sextant==2)) ? v : ((sextant==4)|(sextant==5)) ? m : (sextant==0) ? mid1 : mid2;
	double B = ((sextant==3)|(sextant==4)) ? v : ((sextant==0)|(sextant==1)) ? m : (sextant==2) ? mid1 : mid2;

	int dark = (v <= 0.0);
	*r = dark ? 0.0 : R;
	*g = dark ? 0.0 : G;
	*b = dark ? 0.0 : B;
}

/* Saturation-tool variants: h on [0,6] (NOT divided by 6), and a black-floor
 * `low`: pixels with (m+v) below 2*low return l=0 so the caller skips them
 * (their h,s are then unused). Otherwise identical to the cores above. */
static inline SIRIL_VECTORIZABLE
void rgb_to_hsl_float_sat(float r, float g, float b, float low,
                          float *h, float *s, float *l) {
	float v = (r > g) ? r : g; v = (v > b) ? v : b;
	float m = (r < g) ? r : g; m = (m < b) ? m : b;
	float vm = v - m;
	float L  = 0.5f * (m + v);

	float inv = 1.0f / vm;
	float r2 = (v - r) * inv, g2 = (v - g) * inv, b2 = (v - b) * inv;
	float sden = (L <= 0.5f) ? (v + m) : (2.0f - v - m);
	float S = vm / sden;

	int mr = (r == v);
	int mg = (g == v) & (mr ^ 1);
	int mb = (mr | mg) ^ 1;
	float H = mr * (g == m ? 5.0f + b2 : 1.0f - g2)
	        + mg * (b == m ? 1.0f + r2 : 3.0f - b2)
	        + mb * (r == m ? 3.0f + g2 : 5.0f - r2);   /* h on [0,6], no /6 */

	int below = (m + v < low + low);
	int keep  = (vm > 0.0f) & (below ^ 1);
	*l = below ? 0.0f : L;
	*s = keep ? S : 0.0f;
	*h = keep ? H : 0.0f;
}

static inline SIRIL_VECTORIZABLE
void hsl_to_rgb_float_sat(float h, float sl, float l, float *r, float *g, float *b) {
	h = (h >= 6.0f) ? h - 6.0f : h;         /* h already on [0,6] */
	float v = (l <= 0.5f) ? (l * (1.f + sl)) : (l + sl - l * sl);
	float m = l + l - v;
	float sv = (v - m) / v;
	int sextant = (int) h;                  /* 0..5 */
	float fract = h - (float) sextant;
	float vsf = v * sv * fract;
	float mid1 = m + vsf, mid2 = v - vsf;

	float R = ((sextant==0)|(sextant==5)) ? v : ((sextant==2)|(sextant==3)) ? m : (sextant==4) ? mid1 : mid2;
	float G = ((sextant==1)|(sextant==2)) ? v : ((sextant==4)|(sextant==5)) ? m : (sextant==0) ? mid1 : mid2;
	float B = ((sextant==3)|(sextant==4)) ? v : ((sextant==0)|(sextant==1)) ? m : (sextant==2) ? mid1 : mid2;

	int dark = (v <= 0.0f);
	*r = dark ? 0.0f : R;
	*g = dark ? 0.0f : G;
	*b = dark ? 0.0f : B;
}

/* ----------------------------------------------------------------------------
 * Branchless, vectorisable HSL <-> RGB float cores (moved here from colors.c so
 * they inline at the call site). They are branch-free (selects only) so that,
 * under SIRIL_VECTORIZABLE plus "#pragma omp simd", the compiler emits packed
 * code. The hue dispatch uses *exclusive* masks (priority r > g > b): this fixes
 * a bug in the previous independent-mask version where a two-channel tie for the
 * maximum (e.g. a clipped (1,1,0) highlight) reported a hue 60 degrees off.
 * Behaviour is otherwise identical to the old out-of-line float versions.
 * ------------------------------------------------------------------------- */
static inline SIRIL_VECTORIZABLE
void rgb_to_hslf(float r, float g, float b, float *h, float *s, float *l) {
	float v = (r > g) ? r : g; v = (v > b) ? v : b;
	float m = (r < g) ? r : g; m = (m < b) ? m : b;
	float vm = v - m;
	float L  = 0.5f * (m + v);

	float inv = 1.0f / vm;                   /* +inf when vm==0; discarded below */
	float r2 = (v - r) * inv;
	float g2 = (v - g) * inv;
	float b2 = (v - b) * inv;

	float sden = (L <= 0.5f) ? (v + m) : (2.0f - v - m);
	float S = vm / sden;

	int mr = (r == v);
	int mg = (g == v) & (mr ^ 1);            /* exclusive masks => r>g>b priority */
	int mb = (mr | mg) ^ 1;
	float H = mr * (g == m ? 5.0f + b2 : 1.0f - g2)
	        + mg * (b == m ? 1.0f + r2 : 3.0f - b2)
	        + mb * (r == m ? 3.0f + g2 : 5.0f - r2);
	H *= (1.0f / 6.0f);

	int chroma = (vm > 0.0f);
	*l = L;
	*s = chroma ? S : 0.0f;
	*h = chroma ? H : 0.0f;
}

static inline SIRIL_VECTORIZABLE
void hsl_to_rgbf(float h, float s, float l, float *r, float *g, float *b) {
	h -= floorf(h);                          /* range-reduce to [0,1) (-> roundps) */
	float v = (l <= 0.5f) ? (l * (1.0f + s)) : (l + s - l * s);
	float m = l + l - v;
	float sv = (v - m) / v;                  /* v<=0 lane discarded by select below */
	float h6 = h * 6.0f;
	int sextant = (int) h6;                  /* 0..5 */
	float fract = h6 - (float) sextant;
	float vsf = v * sv * fract;
	float mid1 = m + vsf;
	float mid2 = v - vsf;

	float R = ((sextant==0)|(sextant==5)) ? v : ((sextant==2)|(sextant==3)) ? m : (sextant==4) ? mid1 : mid2;
	float G = ((sextant==1)|(sextant==2)) ? v : ((sextant==4)|(sextant==5)) ? m : (sextant==0) ? mid1 : mid2;
	float B = ((sextant==3)|(sextant==4)) ? v : ((sextant==0)|(sextant==1)) ? m : (sextant==2) ? mid1 : mid2;

	int dark = (v <= 0.0f);
	*r = dark ? 0.0f : R;
	*g = dark ? 0.0f : G;
	*b = dark ? 0.0f : B;
}

void rgbw_to_hslw(uint16_t r, uint16_t g, uint16_t b, uint16_t *h, uint16_t *s, uint16_t *l);
void hslw_to_rgbw(uint16_t h, uint16_t s, uint16_t l, uint16_t *r, uint16_t *g, uint16_t *b);

/* ----------------------------------------------------------------------------
 * Branchless HSV <-> RGB (h,s,v on [0,1], h=0 for grey). Exclusive-mask hue
 * dispatch (priority cmax==r>g>b) + trailing selects, same tie fix as HSL.
 * ------------------------------------------------------------------------- */
static inline SIRIL_VECTORIZABLE
void rgb_to_hsv(double r, double g, double b, double *h, double *s, double *v) {
	double cmax = (r > g) ? r : g; cmax = (cmax > b) ? cmax : b;
	double cmin = (r < g) ? r : g; cmin = (cmin < b) ? cmin : b;
	double delta = cmax - cmin;
	double inv = 1.0 / delta;               /* +inf when delta==0; discarded below */
	int mr = (cmax == r);
	int mg = (cmax == g) & (mr ^ 1);
	int mb = (mr | mg) ^ 1;
	double H = mr * ((g - b) * inv)
	         + mg * (((b - r) * inv) + 2.0)
	         + mb * (((r - g) * inv) + 4.0);
	H /= 6.0;
	H = (H < 0.0) ? H + 1.0 : H;
	int chroma = (delta > 0.0);
	*v = cmax;
	*s = chroma ? (delta / cmax) : 0.0;
	*h = chroma ? H : 0.0;
}

static inline SIRIL_VECTORIZABLE
void rgb_to_hsvf(float r, float g, float b, float *h, float *s, float *v) {
	float cmax = (r > g) ? r : g; cmax = (cmax > b) ? cmax : b;
	float cmin = (r < g) ? r : g; cmin = (cmin < b) ? cmin : b;
	float delta = cmax - cmin;
	float inv = 1.0f / delta;
	int mr = (cmax == r);
	int mg = (cmax == g) & (mr ^ 1);
	int mb = (mr | mg) ^ 1;
	float H = mr * ((g - b) * inv)
	        + mg * (((b - r) * inv) + 2.0f)
	        + mb * (((r - g) * inv) + 4.0f);
	H /= 6.0f;
	H = (H < 0.0f) ? H + 1.0f : H;
	int chroma = (delta > 0.0f);
	*v = cmax;
	*s = chroma ? (delta / cmax) : 0.0f;
	*h = chroma ? H : 0.0f;
}

static inline SIRIL_VECTORIZABLE
void hsv_to_rgb(double h, double s, double v, double *r, double *g, double *b) {
	h -= floor(h);                          /* range-reduce to [0,1) */
	h *= 6.0;
	int i = (int) h;                        /* 0..5 */
	double f = h - (double) i;
	double p = v * (1.0 - s);
	double q = v * (1.0 - (s * f));
	double t = v * (1.0 - (s * (1.0 - f)));
	*r = ((i==0)|(i==5)) ? v : ((i==2)|(i==3)) ? p : (i==1) ? q : t;
	*g = ((i==1)|(i==2)) ? v : ((i==4)|(i==5)) ? p : (i==0) ? t : q;
	*b = ((i==3)|(i==4)) ? v : ((i==0)|(i==1)) ? p : (i==2) ? t : q;
}

static inline SIRIL_VECTORIZABLE
void hsv_to_rgbf(float h, float s, float v, float *r, float *g, float *b) {
	h -= floorf(h);
	h *= 6.0f;
	int i = (int) h;
	float f = h - (float) i;
	float p = v * (1.0f - s);
	float q = v * (1.0f - (s * f));
	float t = v * (1.0f - (s * (1.0f - f)));
	*r = ((i==0)|(i==5)) ? v : ((i==2)|(i==3)) ? p : (i==1) ? q : t;
	*g = ((i==1)|(i==2)) ? v : ((i==4)|(i==5)) ? p : (i==0) ? t : q;
	*b = ((i==3)|(i==4)) ? v : ((i==0)|(i==1)) ? p : (i==2) ? t : q;
}
/* ----------------------------------------------------------------------------
 * RGB <-> XYZ <-> CIE L*a*b* and linear-RGB <-> XYZ and RGB <-> YUV.
 *
 * These transforms have no data-dependent control flow (only the sRGB gamma /
 * cube-root ternaries, which are selects, and a loop-invariant scale flag), so
 * they are moved here as static inline SIRIL_VECTORIZABLE. The matrix-only
 * conversions (linrgb<->xyz, yuv) vectorise fully; the gamma/Lab ones still
 * carry a pow()/cbrt-style call but inline and branch-free at the call site.
 * Only the float variants are kept - the double ones had no remaining callers.
 * ------------------------------------------------------------------------- */

static inline SIRIL_VECTORIZABLE
void rgb_to_xyzf(float r, float g, float b, float *x, float *y, float *z) {
	r = (r <= 0.04045f) ? r / 12.92f : powf(((r + 0.055f) / 1.055f), 2.4f);
	g = (g <= 0.04045f) ? g / 12.92f : powf(((g + 0.055f) / 1.055f), 2.4f);
	b = (b <= 0.04045f) ? b / 12.92f : powf(((b + 0.055f) / 1.055f), 2.4f);
	r *= 100.f; g *= 100.f; b *= 100.f;
	*x = 0.4124564f * r + 0.3575761f * g + 0.1804375f * b;
	*y = 0.2126729f * r + 0.7151522f * g + 0.0721750f * b;
	*z = 0.0193339f * r + 0.1191920f * g + 0.9503041f * b;
}

/* linrgb<->xyz: a scale argument is provided because it is sometimes convenient
 * to have the Y channel output in the range [0..1] instead of [0..100]. Only the
 * float variants are kept (the double ones had no remaining callers). */
static inline SIRIL_VECTORIZABLE
void linrgb_to_xyzf(float r, float g, float b, float *x, float *y, float *z, gboolean scale) {
	if (scale) { r *= 100.f; g *= 100.f; b *= 100.f; }
	*x = 0.4124564f * r + 0.3575761f * g + 0.1804375f * b;
	*y = 0.2126729f * r + 0.7151522f * g + 0.0721750f * b;
	*z = 0.0193339f * r + 0.1191920f * g + 0.9503041f * b;
}

static inline SIRIL_VECTORIZABLE
void xyz_to_linrgbf(float x, float y, float z, float *r, float *g, float *b, gboolean scale) {
	if (scale) { x /= 100.f; y /= 100.f; z /= 100.f; }
	*r =  3.2404542f * x - 1.5371385f * y - 0.4985314f * z;
	*g = -0.9692660f * x + 1.8760108f * y + 0.0415560f * z;
	*b =  0.0556434f * x - 0.2040259f * y + 1.0572252f * z;
}

static inline SIRIL_VECTORIZABLE
void xyz_to_rgbf(float x, float y, float z, float *r, float *g, float *b) {
	x /= 100.0f; y /= 100.0f; z /= 100.0f;
	*r =  3.2404542f * x - 1.5371385f * y - 0.4985314f * z;
	*g = -0.9692660f * x + 1.8760108f * y + 0.0415560f * z;
	*b =  0.0556434f * x - 0.2040259f * y + 1.0572252f * z;
	*r = (*r > 0.0031308f) ? 1.055f * (powf(*r, (1.f / 2.4f))) - 0.055f : 12.92f * (*r);
	*g = (*g > 0.0031308f) ? 1.055f * (powf(*g, (1.f / 2.4f))) - 0.055f : 12.92f * (*g);
	*b = (*b > 0.0031308f) ? 1.055f * (powf(*b, (1.f / 2.4f))) - 0.055f : 12.92f * (*b);
}

static inline SIRIL_VECTORIZABLE
void xyz_to_LABf(float x, float y, float z, float *L, float *a, float *b) {
	x /= 95.047f; y /= 100.000f; z /= 108.883f;
	x = (x > 0.008856452f) ? powf(x, 1.f / 3.0f) : (7.787037037f * x) + (16.f / 116.f);
	y = (y > 0.008856452f) ? powf(y, 1.f / 3.0f) : (7.787037037f * y) + (16.f / 116.f);
	z = (z > 0.008856452f) ? powf(z, 1.f / 3.0f) : (7.787037037f * z) + (16.f / 116.f);
	*L = (116.0f * y) - 16.0f;
	*a = 500.0f * (x - y);
	*b = 200.0f * (y - z);
}

static inline SIRIL_VECTORIZABLE
void LAB_to_xyzf(float L, float a, float b, float *x, float *y, float *z) {
	*y = (L + 16.0f) / 116.0f;
	*x = a / 500.0f + (*y);
	*z = *y - b / 200.0f;
	float x3 = (*x) * (*x) * (*x), y3 = (*y) * (*y) * (*y), z3 = (*z) * (*z) * (*z);
	*x = (x3 > 0.008856452f) ? x3 : (*x - 16.f / 116.f) / 7.787037037f;
	*y = (y3 > 0.008856452f) ? y3 : (*y - 16.f / 116.f) / 7.787037037f;
	*z = (z3 > 0.008856452f) ? z3 : (*z - 16.f / 116.f) / 7.787037037f;
	*x *= 95.047f; *y *= 100.000f; *z *= 108.883f;
}

/* Orthonormal YUV (used by the colour-preserving deconvolution path). The
 * sqrt constants fold to compile-time literals. Pure matrix => vectorises. */
static inline SIRIL_VECTORIZABLE
void rgb_to_yuvf(float red, float green, float blue, float *y, float *u, float *v) {
	const float a = 1.f / sqrtf(3.f);
	const float b = 1.f / sqrtf(2.f);
	const float c = 2.f * a * sqrtf(2.f);
	*y = a * (red + green + blue);
	*u = b * (red - blue);
	*v = c * (0.25f * red - 0.5f * green + 0.25f * blue);
}

static inline SIRIL_VECTORIZABLE
void yuv_to_rgbf(float y, float u, float v, float *red, float *green, float *blue) {
	const float a = 1.f / sqrtf(3.f);
	const float b = 1.f / sqrtf(2.f);
	const float c = a / b;
	*red   = (a * y) + (b * u) + (c * 0.5f * v);
	*green = (a * y) - (c * v);
	*blue  = (a * y) - (b * u) + (c * 0.5f * v);
}

double BV_to_T(double BV);
double T_to_BV(double T);

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
gchar *ccm_log_hook(gpointer p, log_hook_detail detail);
#endif
