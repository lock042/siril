#ifndef FILTERS_GRAXPERT_H
#define FILTERS_GRAXPERT_H

typedef enum {
	GRAXPERT_BG,
	GRAXPERT_DENOISE,
	GRAXPERT_GUI
} graxpert_operation;

typedef enum {
	GRAXPERT_BG_UNSET,
	GRAXPERT_BG_RBF,
	GRAXPERT_BG_KRIGING,
	GRAXPERT_BG_AI
} graxpert_bg_algo;

typedef enum {
	GRAXPERT_MODE_UNSET,
	GRAXPERT_SUBTRACTION,
	GRAXPERT_DIVISON
} graxpert_bg_mode;


typedef struct _graxpert_data {
	graxpert_operation operation;
	double bg_smoothing;
	graxpert_bg_algo bg_algo;
	graxpert_bg_mode bg_mode;
	gboolean keep_bg;
	double denoise_strength;
	gboolean use_gpu;
	gchar *path;
	cmsHPROFILE backup_icc;
} graxpert_data;

void free_graxpert_data(graxpert_data *p);
gpointer do_graxpert (gpointer p);

#endif
