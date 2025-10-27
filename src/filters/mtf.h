#ifndef _MTF_H_
#define _MTF_H_

#include "core/siril.h"

struct mtf_params {
	float midtones, shadows, highlights;
	gboolean do_red, do_green, do_blue;
};

/* Auto-stretch parameters */
#define AS_DEFAULT_SHADOWS_CLIPPING -2.80f /* Shadows clipping point measured in sigma units from the main histogram peak */
#define AS_DEFAULT_TARGET_BACKGROUND 0.25f /* final "luminance" of the image for autostretch in the [0,1] range */

float MTF_pseudoinverse(float y, struct mtf_params);
float MTF(float x, float m, float lo, float hi);
float MTFp(float x, struct mtf_params params);

void apply_linked_pseudoinverse_mtf_to_fits(fits *from, fits *to, struct mtf_params params, gboolean multithreaded);
void apply_unlinked_pseudoinverse_mtf_to_fits(fits *from, fits *to, struct mtf_params *params, gboolean multithreaded);
void apply_linked_mtf_to_fits(fits *from, fits *to, struct mtf_params params, gboolean multithreaded);
int find_linked_midtones_balance(fits *fit, float shadows_clipping, float target_bg, struct mtf_params *result);
int find_linked_midtones_balance_default(fits *fit, struct mtf_params *result); // with default args

void apply_unlinked_mtf_to_fits(fits *from, fits *to, struct mtf_params *params);
int find_unlinked_midtones_balance(fits *fit, float shadows_clipping, float target_bg, struct mtf_params *results);
int find_unlinked_midtones_balance_default(fits *fit, struct mtf_params *results);

#endif
