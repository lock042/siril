#ifndef _DECONVOLUTION_GUI_H
#define _DECONVOLUTION_GUI_H

#include "filters/deconvolution/deconvolution.h"
#include "core/siril.h"

typedef struct deconvolution_sequence_data {
	sequence *seq;
	char* seqEntry;
	estk_data *deconv_data;
	gboolean from_command;
} deconvolution_sequence_data;

void reset_conv_args(estk_data* args);

int load_kernel(gchar* filename);
int save_kernel(gchar* filename);
void on_bdeconv_savekernel_clicked(GtkButton *button, gpointer user_data);
gpointer deconvolve(gpointer p);
gpointer estimate_only(gpointer p);
gpointer deconvolve_sequence_command(gpointer p, sequence* seqname);

void apply_deconvolve_to_sequence(struct deconvolution_sequence_data *seqdata);

#endif
