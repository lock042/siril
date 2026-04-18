#ifndef SRC_GUI_WAVELETS_H_
#define SRC_GUI_WAVELETS_H_

#include "core/siril.h"
#include "core/processing.h"

/* wavelets filter data from GUI */
struct wavelets_filter_data {
	fits *fit;
	int Nbr_Plan;
	int Type;
};

/* Data struct for wavelet reconstruction (wrecons command and GUI OK path) */
struct wrecons_data {
	destructor destroy_fn;  /* Must be first member */
	float coef[7];
	int nb_chan;
};

/* Data struct for wavelet decomposition (wavelet command and GUI compute path) */
struct wavelet_transform_data {
	int Nbr_Plan;
	int Type_Transform;
};

void apply_wavelets_cancel();
int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer);
gpointer extract_plans(gpointer p);
void free_wrecons_data(void *p);
int wrecons_image_hook(struct generic_img_args *args, fits *fit, int threads);
gchar *wrecons_log_hook(gpointer p, log_hook_detail detail);
gpointer wavelet_transform_worker(gpointer p);

#endif /* SRC_GUI_WAVELETS_H_ */
