#ifndef _HIST_H_
#define _HIST_H_

#include "filters/mtf.h"
#include "filters/ght.h"

#define NO_STRETCH_SET_YET 0
#define HISTO_STRETCH 1
#define GHT_STRETCH 2

struct mtf_data {
	fits *fit;
	sequence *seq;
	struct mtf_params params;
	const gchar *seqEntry;
};

struct ght_data {
	fits *fit;
	sequence *seq;
	struct ght_params params_ght;
	struct ght_compute_params compute_params;
	const gchar *seqEntry;
};

typedef enum {
	SCALE_LOW,
	SCALE_MID,
	SCALE_HI
} ScaleType;

gsl_histogram* computeHisto(fits*, int);
gsl_histogram* computeHisto_Selection(fits*, int, rectangle *);
void compute_histo_for_gfit();
void invalidate_gfit_histogram();
void update_gfit_histogram_if_needed();
void apply_histo_cancel();
void toggle_histogram_window_visibility();

void on_histoMidEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histoShadEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histoHighEntry_changed(GtkEditable *editable, gpointer user_data);

void apply_mtf_to_sequence(struct mtf_data *mtf_args);
void apply_ght_to_sequence(struct ght_data *ght_args);

#endif
