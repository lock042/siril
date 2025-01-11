#ifndef _SEQ_FILTERING_H_
#define _SEQ_FILTERING_H_

#include "core/siril.h"

#define MAX_FILTERS 7

struct stacking_args;
struct stacking_configuration;

typedef int (*seq_image_filter)(sequence *seq, int img_index, double param);

struct filtering_tuple {
	seq_image_filter filter;
	double param;
};

int seq_filter_all(sequence *seq, int nb_img, double any);
int seq_filter_included(sequence *seq, int nb_img, double any);
int seq_filter_fwhm(sequence *seq, int nb_img, double max_fwhm);
int seq_filter_weighted_fwhm(sequence *seq, int nb_img, double max_fwhm);
int seq_filter_quality(sequence *seq, int nb_img, double max_quality);
int seq_filter_roundness(sequence *seq, int nb_img, double min_rnd);
int seq_filter_background(sequence *seq, int nb_img, double max_bkg);
int seq_filter_nbstars(sequence *seq, int nb_img, double min_nbstars);

int compute_nb_filtered_images(sequence *seq, seq_image_filter filtering_criterion, double filtering_parameter);

seq_image_filter create_filter_prefixed_nonexisting_output(const char *prefix);

seq_image_filter create_multiple_filter(seq_image_filter filter1, double fparam1, ...);
seq_image_filter create_multiple_filter_from_list(struct filtering_tuple *filters);

/* parsing filter data and building the actual filter */

struct seq_filter_config {
	float f_fwhm, f_fwhm_p, f_wfwhm, f_wfwhm_p, f_round, f_round_p, f_quality, f_quality_p, f_bkg, f_bkg_p, f_nbstars, f_nbstars_p; // on if >0
	gboolean f_fwhm_k, f_wfwhm_k, f_round_k, f_quality_k, f_bkg_k, f_nbstars_k;
	gboolean filter_included;
};

int convert_parsed_filter_to_filter(struct seq_filter_config *arg, sequence *seq, seq_image_filter *criterion, double *param);
int setup_filtered_data(struct stacking_args *args);

int stack_fill_list_of_unfiltered_images(struct stacking_args *args);
double compute_highest_accepted_fwhm(sequence *seq, double criterion, gboolean is_ksigma);
double compute_highest_accepted_weighted_fwhm(sequence *seq, double criterion, gboolean is_ksigma);
double compute_lowest_accepted_quality(sequence *seq, double criterion, gboolean is_ksigma);
double compute_lowest_accepted_roundness(sequence *seq, double criterion, gboolean is_ksigma);
double compute_highest_accepted_background(sequence *seq, double criterion, gboolean is_ksigma);
double compute_lowest_accepted_nbstars(sequence *seq, double criterion, gboolean is_ksigma);


gchar *describe_filter(sequence *seq, seq_image_filter filtering_criterion, double filtering_parameter);


#endif
