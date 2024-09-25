#ifndef _EXTRACTION_H
#define _EXTRACTION_H

#include <glib.h>
#include "core/siril.h"
#include "core/processing.h"
#include "io/ser.h"
#include "io/fits_sequence.h"

struct simple_extract_data {
	sequence *seq;
	char *seqEntry;
	extraction_scaling scaling;
};

void update_sampling_information(fits *fit, float factor);
void update_filter_information(fits *fit, char *filter, gboolean append);

int extractHa_ushort(fits *in, fits *Ha, sensor_pattern pattern, extraction_scaling scaling);
int extractHa_float(fits *in, fits *Ha, sensor_pattern pattern, extraction_scaling scaling);
void apply_extractHa_to_sequence(struct simple_extract_data *multi_args);

int extractGreen_ushort(fits *in, fits *green, sensor_pattern pattern);
int extractGreen_float(fits *in, fits *green, sensor_pattern pattern);
void apply_extractGreen_to_sequence(struct simple_extract_data *multi_args);

int extractHaOIII_ushort(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern, extraction_scaling scaling, int threads);
int extractHaOIII_float(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern, extraction_scaling scaling, int threads);
void apply_extractHaOIII_to_sequence(struct multi_output_data *multi_args);

int split_cfa_ushort(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3);
int split_cfa_float(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3);
void apply_split_cfa_to_sequence(struct multi_output_data *multi_args);

#endif
