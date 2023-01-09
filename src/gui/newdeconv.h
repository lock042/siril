#ifndef _DECONVOLUTION_GUI_H
#define _DECONVOLUTION_GUI_H

#include "filters/deconvolution/deconvolution.h"
#include "core/siril.h"

typedef struct deconvolution_sequence_data {
	sequence *seq;
	const gchar* seqEntry;
	estk_data *deconv_data;
} deconvolution_sequence_data;

void reset_conv_args(estk_data* args);

gpointer deconvolve(gpointer p);
gpointer estimate_only(gpointer p);
gpointer deconvolve_sequence_command(gpointer p, sequence* seqname);

void apply_deconvolve_to_sequence(struct deconvolution_sequence_data *seqdata);

#endif
