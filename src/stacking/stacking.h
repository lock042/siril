#ifndef STACKING_H_
#define STACKING_H_

#include "core/processing.h"
#include "core/sequence_filtering.h"

//#define STACK_DEBUG

/* the stacking method */
struct stacking_args;
typedef int (*stack_method)(struct stacking_args *args);

typedef struct normalization_coeff norm_coeff;

enum {
	ST_ALLOC_ERROR = -10,
	ST_CANCEL = -9,
	ST_SEQUENCE_ERROR = -2,
	ST_GENERIC_ERROR = -1,
	ST_OK = 0
};

typedef enum {
	STACK_SUM,
	STACK_MEAN,
	STACK_MEDIAN,
	STACK_MAX,
	STACK_MIN,
} stackMethod;

/* identical to the combo box items */
typedef enum {
	ALL_IMAGES,
	SELECTED_IMAGES,
	BEST_PSF_IMAGES,
	BEST_WPSF_IMAGES,
	BEST_ROUND_IMAGES,
	BEST_BKG_IMAGES,
	BEST_NBSTARS_IMAGES,
	BEST_QUALITY_IMAGES
} stackType;

typedef enum {
	NO_WEIGHT,
	NBSTARS_WEIGHT,
	WFWHM_WEIGHT,
	NOISE_WEIGHT,
	NBSTACK_WEIGHT
} weighingType;

struct normalization_coeff {
	double *offset;
	double *mul;
	double *scale;
	double *poffset[3];
	double *pmul[3];
	double *pscale[3];
};

/* Warning: editing this will require changing the long initializer in stacking.c */
struct stacking_args {
	stack_method method;
	sequence *seq;
	int ref_image; // takes precedence over seq->reference_image which may not be applicable
	seq_image_filter filtering_criterion;
	double filtering_parameter;
	int nb_images_to_stack; 	/* calculated from the above */
	int *image_indices;		/* mapping between selected and sequence image indices */
	char *description;		/* description of the filtering */

	const char *output_filename;	/* used in the idle function only */
	gchar *output_parsed_filename;	/* used in the idle function only */
	gboolean output_overwrite;	/* used in the idle function only */

	gboolean lite_norm;		/* enable lightweight (med,mad) normalization */
	normalization normalize;	/* type of normalization */
	norm_coeff coeff;		/* normalization data */
	gboolean force_norm;		/* TRUE = force normalization */
	gboolean output_norm;		/* normalize final image to the [0, 1] range */
	gboolean use_32bit_output;	/* output to 32 bit float */
	int reglayer;			/* layer used for registration data */
	gboolean equalizeRGB;		/* enable RGB equalization through normalization */

	rejection type_of_rejection;	/* type of rejection */
	float sig[2];			/* low and high sigma rejection or GESTD parameters */
	float *critical_value;		/* index of critical_values for GESTD */
	gboolean create_rejmaps;	/* main activation flag for the rejection maps creation */
	gboolean merge_lowhigh_rejmaps;	/* create only one map */
	fits *rejmap_low, *rejmap_high;	/* rejection maps */

	gboolean apply_noise_weights;	/* enable weights */
	gboolean apply_nbstack_weights;	/* enable weights */
	gboolean apply_wfwhm_weights;	/* enable weights */
	gboolean apply_nbstars_weights;	/* enable weights */
	double *weights; 		/* computed weights for each (layer, image)*/

	float (*sd_calculator)(const WORD *, const int); // internal, for ushort
	float (*mad_calculator)(const WORD *, const size_t, const double, threading_type) ; // internal, for ushort

	struct timeval t_start;
	int retval;
	fits result;
};

/* configuration from the command line */
struct stacking_configuration {
	struct timeval t_start;
	gchar *seqfile;
	gchar *result_file;
	stack_method method;
	rejection type_of_rejection;
	double sig[2];
	gboolean create_rejmaps;
	gboolean merge_lowhigh_rejmaps;
	gboolean force_no_norm;
	gboolean output_norm;
	gboolean equalizeRGB;
	gboolean lite_norm;
	normalization norm;
	int number_of_loaded_sequences;
	struct seq_filter_config filters;
	gboolean apply_noise_weights;
	gboolean apply_nbstack_weights;
	gboolean apply_wfwhm_weights;
	gboolean apply_nbstars_weights;
};

typedef struct {
	GDateTime *date_obs;
	gdouble exposure;
} DateEvent;


/*
 * ESD test statistic data.
 */
struct outliers {
	float x;
	int i;
	int out;
};


void initialize_stacking_methods();
gboolean evaluate_stacking_should_output_32bits(const stack_method method,
		sequence *seq, int nb_img_to_stack, gchar **err);

int stack_median(struct stacking_args *args);
int stack_mean_with_rejection(struct stacking_args *args);
int stack_addmax(struct stacking_args *args);
int stack_addmin(struct stacking_args *args);
gboolean stack_regdata_is_valid(struct stacking_args *args);
void main_stack(struct stacking_args *args);
void clean_end_stacking(struct stacking_args *args);

void get_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter);
void update_stack_interface(gboolean dont_change_stack_type);

void describe_stack_for_history(struct stacking_args *args, GSList **hist, gboolean for_rejmap, gboolean low_rejmap);

	/* normalization functions, normalize.c */

int do_normalization(struct stacking_args *args);
int *compute_thread_distribution(int nb_workers, int max);


	/* median and mean functions */

int check_G_values(float Gs, float Gc);
void confirm_outliers(struct outliers *out, int N, double median, int *rejected, int rej[2]);

struct _image_block {
	long channel, start_row, end_row, height; // long matches naxes type
};

int stack_compute_parallel_blocks(struct _image_block **blocksptr, long max_number_of_rows,
		const long naxes[3], int nb_threads, long *largest_block_height, int *nb_blocks);

/* pool of memory blocks for parallel processing */
struct _data_block {
	void *tmp;	// the actual single buffer for all others below
	void **pix;	// buffer for a block on all images
	void *stack;	// the reordered stack for one pixel in all images
	int *rejected;	// 0 if pixel ok, 1 or -1 if rejected
	void *o_stack;	// original unordered stack
	void *w_stack;	// stack for the winsorized rejection
	float *xf, *yf, m_x, m_dx2;// data for the linear fit rejection
	int layer;	// to identify layer for normalization
};

int find_refimage_in_indices(const int *indices, int nb, int ref);

	/* up-scaling functions */

int upscale_sequence(struct stacking_args *args);
void remove_tmp_drizzle_files(struct stacking_args *args);


	/* rejection_float.c */

int apply_rejection_float(struct _data_block *data, int nb_frames, struct stacking_args *args, int crej[2]);

	/* keeping metadata */
void free_list_date(gpointer data);
DateEvent* new_date_item(GDateTime *dt, gdouble exposure);
void compute_date_time_keywords(GList *list_date, fits *fit);

#endif
