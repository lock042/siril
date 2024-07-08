#ifndef FILTERS_GRAXPERT_H
#define FILTERS_GRAXPERT_H

typedef enum {
	GRAXPERT_BG,
	GRAXPERT_DENOISE,
	GRAXPERT_GUI
} graxpert_operation;

typedef enum {
	GRAXPERT_BG_AI,
	GRAXPERT_BG_RBF,
	GRAXPERT_BG_KRIGING,
	GRAXPERT_BG_SPLINE,
} graxpert_bg_algo;

typedef enum {
	GRAXPERT_SUBTRACTION,
	GRAXPERT_DIVISION
} graxpert_bg_mode;

typedef enum {
	GRAXPERT_THIN_PLATE,
	GRAXPERT_QUINTIC,
	GRAXPERT_CUBIC,
	GRAXPERT_LINEAR
} graxpert_rbf_kernel;

typedef enum {
	STRETCH_OPTION_NONE,
	STRETCH_OPTION_10_BG_3_SIGMA,
	STRETCH_OPTION_15_BG_3_SIGMA,
	STRETCH_OPTION_20_BG_3_SIGMA,
	STRETCH_OPTION_30_BG_2_SIGMA
} graxpert_stretch;

typedef struct _graxpert_data {
	fits *fit;
	sequence *seq;
	GSList *bg_samples;
	graxpert_operation operation;
	double bg_smoothing;
	graxpert_bg_algo bg_algo;
	graxpert_bg_mode bg_mode;
	graxpert_stretch stretch_option;
	graxpert_rbf_kernel kernel;
	int sample_size; // size of background samples
	int spline_order;
	double bg_tol_option; // BGE sample tolerance
	gboolean keep_bg;
	double denoise_strength;
	gboolean use_gpu;
	int ai_batch_size;
	int bg_pts_option; // points per row
	gchar *path;
	gchar *configfile;
	cmsHPROFILE backup_icc;
	gboolean previewing;
} graxpert_data;

graxpert_data *new_graxpert_data();
void free_graxpert_data(graxpert_data *p);
gpointer do_graxpert (gpointer p);
void apply_graxpert_to_sequence(graxpert_data *args);

#endif
