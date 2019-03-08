#ifndef STACKING_H_
#define STACKING_H_

#include "core/processing.h"

/* the stacking method */
struct stacking_args;
typedef int (*stack_method)(struct stacking_args *args);

typedef struct normalization_coeff norm_coeff;


/* TYPE OF SIGMA CLIPPING */
typedef enum {
	NO_REJEC,
	PERCENTILE,
	SIGMA,
	SIGMEDIAN,
	WINSORIZED,
	LINEARFIT,
} rejection;

/* TYPE OF NORMALIZATION */
typedef enum {
	NO_NORM,
	ADDITIVE,
	MULTIPLICATIVE,
	ADDITIVE_SCALING,
	MULTIPLICATIVE_SCALING,
} normalization;

/* identical to the combo box items */
typedef enum {
	ALL_IMAGES,
	SELECTED_IMAGES,
	BEST_PSF_IMAGES,
	BEST_ROUND_IMAGES,
	BEST_QUALITY_IMAGES
} stackType;

struct normalization_coeff {
	double *offset;
	double *mul;
	double *scale;
};

struct stacking_args {
	stack_method method;
	sequence *seq;
	seq_image_filter filtering_criterion;
	double filtering_parameter;
	int nb_images_to_stack; // calculated from the above, for display purposes
	int *image_indices;	// conversion between selected image indices and sequence image indices
	char description[100];	// description of the filtering
	char *output_filename;	// used in the idle function only
	gboolean output_overwrite;	// used in the idle function only
	struct timeval t_start;
	int retval;
	int max_number_of_rows;	/* number of rows that can be processed simultaneously,
				   function of max memory, image size and nb_images_to_stack */
	double sig[2];		/* low and high sigma rejection */
	rejection type_of_rejection;	/* type of rejection */
	normalization normalize;	/* type of normalization */
	norm_coeff coeff;		/* normalization data */
	gboolean force_norm;		/* TRUE = force normalization */
	gboolean norm_to_16;		/* normalize final image to 16bits */
	int reglayer;		/* layer used for registration data */
};

void initialize_stacking_methods();

int stack_get_max_number_of_rows(sequence *seq, int nb_images_to_stack);
void stack_fill_list_of_unfiltered_images(struct stacking_args *args);
double compute_highest_accepted_fwhm(double percent);
double compute_lowest_accepted_quality(double percent);
int compute_nb_filtered_images(sequence *seq, seq_image_filter filtering_criterion,
		double filtering_parameter, int reglayer);
int compute_nb_filtered_images_stack(struct stacking_args *stack_args);
double compute_lowest_accepted_roundness(double percent);

int stack_median(struct stacking_args *args);
int stack_mean_with_rejection(struct stacking_args *args);
int stack_addmax(struct stacking_args *args);
int stack_addmin(struct stacking_args *args);

int upscale_sequence(struct stacking_args *args);

gboolean end_stacking(gpointer p);
void clean_end_stacking(struct stacking_args *args);
void update_stack_interface(gboolean dont_change_stack_type);

int stack_filter_all(sequence *seq, int layer, int nb_img, double any);
int stack_filter_included(sequence *seq, int layer, int nb_img, double any);
int stack_filter_included_and_registered(sequence *seq, int layer, int nb_img, double any);
int stack_filter_fwhm(sequence *seq, int layer, int nb_img, double max_fwhm);
int stack_filter_roundness(sequence *seq, int layer, int nb_img, double min_rnd);
int stack_filter_quality(sequence *seq, int layer, int nb_img, double max_quality);

	/* normalization functions, normalize.c */
int do_normalization(struct stacking_args *args);

#endif
