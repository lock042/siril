#ifndef FFT_H
#define FFT_H

#include "core/processing.h"

/* fft data from GUI */
struct fft_data {
	destructor destroy_fn;  /* Must be first member */
	fits *fit;
	char *type;
	gchar *modulus, *phase;
	int type_order;
	int retval;
};

gpointer fourier_transform(gpointer p);
void free_fft_data(void *p);
int fft_image_hook(struct generic_img_args *args, fits *fit, int threads);
gchar *fft_log_hook(gpointer p, log_hook_detail detail);
gboolean fft_idle(gpointer p);
#endif
