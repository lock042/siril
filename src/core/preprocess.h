#ifndef _PREPRO_H_
#define _PREPRO_H_

#include "siril.h"
#include "filters/cosmetic_correction.h"
#include "core/processing.h"

/* preprocessing data from GUI */
struct preprocessing_data {
	gboolean use_bias, use_dark, use_flat;
	fits *bias, *dark, *flat;
	float bias_level;	// the synthetic bias level [0, 1]
	gboolean use_dark_optim;
	gboolean autolevel;	// auto-evaluate flat normalization
	float normalisation;	// the flat normalization level
	gboolean equalize_cfa;	// convert the master flat to gray
	gboolean fix_xtrans;

	/* input and output file or sequence */
	gboolean is_sequence;
	sequence *seq;
	sequence_type output_seqtype;
	gboolean allow_32bit_output;
	const char *ppprefix;	// prefix for output files
	gboolean debayer;	// debayer at the end
	GSList *history;	// generic history to add to the FITS output

	/* cosmetic correction */
	gboolean use_cosmetic_correction, cc_from_dark;
	GFile *bad_pixel_map_file;	// if !cc_from_dark
	double sigma[2];
	long icold, ihot;	// number of cold and hot pixels to correct
	deviant_pixel *dev;	// the runtime list of deviant pixels (icold + ihot long)
	gboolean is_cfa;	// when replacing pixels, don't use direct neighbours

	gboolean ignore_exclusion; // if true, images marked as excluded will not
							// be preprocessed
	int retval;
};

int preprocess_single_image(struct preprocessing_data *args);
int preprocess_given_image(char *file, struct preprocessing_data *args);
int evaluateoffsetlevel(const char* expression, fits *fit);
void start_sequence_preprocessing(struct preprocessing_data *prepro);
gboolean check_for_cosme_file_sanity(GFile *file);

/* used in live stacking */
int prepro_prepare_hook(struct generic_seq_args *args);
int prepro_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads);
void clear_preprocessing_data(struct preprocessing_data *prepro);

#endif
