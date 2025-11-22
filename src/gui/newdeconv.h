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
gboolean reset_conv_controls(gpointer user_data);
void reset_conv_controls_and_args();
void reset_conv_kernel();
gboolean set_kernel_size_in_gui(gpointer user_data);
gboolean DrawPSF(gpointer user_data);
void calculate_parameters();
int load_kernel(gchar* filename, estk_data *args);
int save_kernel(gchar* filename, estk_data *args);
orientation_t get_imageorientation(fits *fit);

void on_bdeconv_savekernel_clicked(GtkButton *button, gpointer user_data);
gboolean deconvolve_img_idle(gpointer p);
gboolean estimate_img_idle(gpointer p);

gboolean deconvolve_idle(gpointer arg);
gpointer deconvolve(gpointer p);
gboolean estimate_idle(gpointer arg);
gpointer estimate_only(gpointer p);
gpointer deconvolve_sequence_command(gpointer p, sequence* seqname);

void apply_deconvolve_to_sequence(struct deconvolution_sequence_data *seqdata);

void close_deconv();

#endif
