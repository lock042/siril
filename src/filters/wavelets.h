#ifndef SRC_FILTERS_WAVELETS_H_
#define SRC_FILTERS_WAVELETS_H_

#include "core/siril.h"
#include "core/processing.h"
#include "algos/wavelet_denoise.h"

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
	struct denoise_params denoise; /* per-scale denoising (disabled by default) */
	/* ROI preview: when for_roi is set the full image is reconstructed and the
	 * (roi_x, roi_y, full_rx*full_ry-relative) selection window is copied into
	 * the smaller fit the hook receives. */
	gboolean for_roi;
	int roi_x, roi_y;
	int full_rx, full_ry;
};

/* Data for the higher-level à trous transform (atrous / seqatrous commands and
 * the GUI "Apply to sequence" path): a full decompose + per-scale denoise +
 * weighted reconstruct in one shot, with a fixed set of levels and options. */
struct atrous_data {
	destructor destroy_fn;         /* Must be first member */
	sequence *seq;                 /* sequence to process (seqatrous); NULL otherwise */
	fits *fit;                     /* single-image target (atrous); NULL in seq mode */
	char *seqEntry;                /* output sequence prefix (default "wv_") */
	int nbr_plan;                  /* number of wavelet layers */
	int type;                      /* TO_PAVE_LINEAR or TO_PAVE_BSPLINE */
	gboolean anscombe;             /* decompose/reconstruct in the Anscombe VST domain */
	float coef[7];                 /* per-layer reconstruction weights (1 = neutral) */
	struct denoise_params denoise; /* per-scale denoising options */
};

/* Data struct for wavelet decomposition (wavelet command and GUI compute path) */
struct wavelet_transform_data {
	int Nbr_Plan;
	int Type_Transform;
	gboolean anscombe; /* decompose in the Anscombe VST domain */
	/* Optional completion idle run on the main thread once the decomposition
	 * has finished (GUI re-enables its widgets here). When set, it takes
	 * ownership of this struct and must free it; when NULL the generic idle is
	 * used and the worker frees the struct itself. */
	GSourceFunc idle;
};

int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer);
/* One-shot à trous decompose + denoise + weighted reconstruct of a single image
 * in place. `id` makes the per-channel temporary transform files unique so
 * parallel callers (and the GUI's live decomposition) never collide. Returns 0
 * on success. */
int atrous_transform_image(fits *fit, const struct atrous_data *args, int id);
/* Image hook for the single-image atrous command (used with generic_image_worker). */
int atrous_image_hook(struct generic_img_args *args, fits *fit, int threads);
gchar *atrous_log_hook(gpointer p, log_hook_detail detail);
/* Launch atrous processing over the sequence in args->seq (prefix args->seqEntry),
 * one frame at a time. Takes ownership of args. */
void apply_atrous_to_sequence(struct atrous_data *args);
void free_atrous_data(void *p);
gpointer extract_plans(gpointer p);
void free_wrecons_data(void *p);
int wrecons_image_hook(struct generic_img_args *args, fits *fit, int threads);
gchar *wrecons_log_hook(gpointer p, log_hook_detail detail);
gpointer wavelet_transform_worker(gpointer p);
/* apply_wavelets_cancel declared in gui/wavelets.h */

#endif /* SRC_FILTERS_WAVELETS_H_ */
