#ifndef _REGISTRATION_H_
#define _REGISTRATION_H_

#include "core/siril.h"
#include "algos/PSF.h"
#include "core/processing.h"
#include "algos/star_finder.h"
#include "registration/distorsion.h"

#define NUMBER_OF_METHODS 7

struct registration_args;
typedef int (*registration_function)(struct registration_args *);

typedef enum {
	REQUIRES_NO_SELECTION,	// selection is not used
	REQUIRES_ANY_SELECTION,		// selection can be of any size and shape
	REQUIRES_SQUARED_SELECTION	// selection needs to be square-shaped
} selection_type;

typedef enum {
	REGTYPE_DEEPSKY, REGTYPE_PLANETARY, REGTYPE_APPLY
} registration_type;

typedef enum {
	REG_UNDEF = -1,
	REG_GLOBAL,
	REG_3STARS,
	REG_DFT,
	REG_KOMBAT,
	REG_COMET,
	REG_APPLY,
	REG_2PASS
} regmethod_index;

typedef enum {
	PLANETARY_FULLDISK, PLANETARY_SURFACE
} planetary_type;


typedef enum {
	UNDEFINED_TRANSFORMATION = -3,
	NULL_TRANSFORMATION = -2,
	IDENTITY_TRANSFORMATION = -1,
	SHIFT_TRANSFORMATION,
	SIMILARITY_TRANSFORMATION,
	AFFINE_TRANSFORMATION,
	HOMOGRAPHY_TRANSFORMATION
} transformation_type;

typedef enum {
	FRAMING_CURRENT,
	FRAMING_MAX,
	FRAMING_MIN,
	FRAMING_COG,
} framing_type;

/* used to register a registration method */
struct registration_method {
	const char *name;
	registration_function method_ptr;
	selection_type sel;
	registration_type type;
};

/* used to draw framing in image_display */
typedef struct {
	point pt[4];
} regframe;

/* used to define rotation matrices*/
typedef enum {
	ROTX,
	ROTY,
	ROTZ
} rotation_type;

/* same as rectangle but avoids conflicts with rectangle defined in opencv namespace */
typedef struct {
	int x, y, w, h;
} framing_roi;

struct astrometric_args{
	int nb;
	int refindex;
	Homography *Rs;
	Homography *Ks;
	disto_data *disto;
	float scale;
	framing_roi roi;
};

/**** star alignment (global and 3-star) registration ****/

struct star_align_data {
	struct registration_args *regargs;
	regdata *current_regdata;
	struct astrometric_args *astargs;
	psf_star **refstars;
	int fitted_stars;
	BYTE *success;
	point ref;
};

typedef struct {
	Homography Htransf;
	framing_roi roi_out;
} framing_data;

/* arguments passed to registration functions */
struct registration_args {
	registration_function func;	// the registration function
	sequence *seq;			// the sequence to register
	gboolean bayer;			// whether we are dealing with a Bayer pattern
	int reference_image;		// reference image index
	struct seq_filter_config filters; // parsed image filters (.filter_included always used)
	int layer;			// layer of images on which the registration is computed
	int retval;			// retval of func
	gboolean run_in_thread;		// true if the registration was run in a thread
	gboolean follow_star;		// follow star position between frames
	gboolean matchSelection;	// Match stars found in the seleciton of reference image
	rectangle selection;		// the selection rectangle
	struct starfinder_data *sfargs;		// star finder configuration for global/2pass
	int min_pairs;			// Minimum number of star pairs for success
	int max_stars_candidates;	// Max candidates after psf fitting for global reg
	transformation_type type;	// Use affine transform  or homography
	float percent_moved;		// for KOMBAT algorithm
	pointf velocity;			// for comet algorithm
	gboolean two_pass;		// use the two-pass computation to find a good ref image
	seq_image_filter filtering_criterion; // the filter, (seqapplyreg only)
	double filtering_parameter;	// and its parameter (seqapplyreg only)
	gboolean no_starlist;		// disable star list creation (2pass only)
	float output_scale;		// scaling factor
	disto_source undistort;		// type of undistorsion to apply
	disto_params distoparam;		// type of undistorsion to apply
	disto_data *disto;			// undistorsion information
	struct driz_args_t *driz;	// drizzle-specific data
	framing_type framing;		// used by seqapplyreg to determine framing
	framing_data framingd;	// precomputed framing data

	opencv_interpolation interpolation; // type of rotation interpolation
	gboolean clamp;				// should Bicubic and Lanczos4 interpolation be clamped?
	double clamping_factor;		// used to set amount of interpolation clamping

	/* data for generated sequence, for star alignment/mosaic registration */
	gboolean no_output;		// write transformation to .seq
	int new_total;                  // remaining images after registration
	imgdata *imgparam;		// imgparam for the new sequence
	regdata *regparam;		// regparam for the new sequence
	char *prefix;		// prefix of the created sequence if any
	gboolean load_new_sequence;	// load the new sequence if success
	gchar *new_seq_name;
};

struct registration_method *new_reg_method(const char *name, registration_function f,
		selection_type s, registration_type t); // for compositing
int register_shift_dft(struct registration_args *args); // REG_DFT
int register_shift_fwhm(struct registration_args *args);
int register_star_alignment(struct registration_args *args); // REG_GLOBAL
int register_multi_step_global(struct registration_args *regargs); // REG_2PASS
int register_comet(struct registration_args *regargs); // REG_COMET
int register_3stars(struct registration_args *regargs); // REG_3STARS
int register_apply_reg(struct registration_args *regargs); // REG_APPLY
int register_kombat(struct registration_args *args); // REG_KOMBAT
int register_manual(struct registration_args *regargs); // defined in compositing/compositing.c
int register_astrometric(struct registration_args *regargs);

// 3 stars
void reset_3stars();
int _3stars_get_number_selected_stars();
gboolean _3stars_check_selection();

// comet
pointf compute_velocity(GDateTime *t1, GDateTime *t2, point d1, point d2);
int get_comet_shift(GDateTime *ref, GDateTime *img, pointf px_per_hour, pointf *reg);

pointf get_velocity();
void update_reg_interface(gboolean dont_change_reg_radio);
void compute_fitting_selection(rectangle *area, int hsteps, int vsteps, int preserve_square);
void get_the_registration_area(struct registration_args *regargs, const struct registration_method *method); // for compositing
gpointer register_thread_func(gpointer p);

/** getter */
int get_registration_layer(const sequence *seq);
int seq_has_any_regdata(const sequence *seq); // same as get_registration_layer but does not rely on GUI for com.seq
gboolean seq_has_any_distorsion(const sequence *seq);

regdata *apply_reg_get_current_regdata(struct registration_args *regargs);
regdata *star_align_get_current_regdata(struct registration_args *regargs);
int star_align_prepare_results(struct generic_seq_args *args);
int star_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads);
int star_align_finalize_hook(struct generic_seq_args *args);
void create_output_sequence_for_global_star(struct registration_args *args, int refindex);
int star_match_and_checks(psf_star **ref_stars, psf_star **stars, int nb_ref_stars, int nb_stars, struct registration_args *regargs, int filenum, Homography *H);
int apply_reg_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads);

const char *describe_transformation_type(transformation_type type);

void selection_H_transform(rectangle *selection, Homography Href, Homography Himg);
void guess_transform_from_seq(sequence *seq, int layer,
		transformation_type *min, transformation_type *max, gboolean excludenull);
transformation_type guess_transform_from_H(Homography H);
gboolean check_before_applyreg(struct registration_args *regargs);
gboolean layer_has_registration(const sequence *seq, int layer);
gboolean layer_has_usable_registration(sequence *seq, int layer);
gboolean layer_has_distorsion(const sequence *seq, int layer);
int get_first_selected(sequence *seq);
regdata *apply_reg_get_current_regdata(struct registration_args *regargs);

void translation_from_H(Homography H, double *dx, double *dy);
Homography H_from_translation(double dx, double dy);
void SetNullH(Homography *H);
int shift_fit_from_reg(fits *fit, Homography H);
void compute_Hmax(Homography *Himg, Homography *Href, int src_rx_in, int src_ry_in, double scale, Homography *H, Homography *Hshift, int *dst_rx_out, int *dst_ry_out);

int compute_Hs_from_astrometry(sequence *seq, struct wcsprm *WCSDATA, framing_type framing, int layer);

int minidx(const float *arr, const gboolean *mask, int nb, float *val);
void free_astrometric_args(struct astrometric_args *astargs);
int initialize_drizzle_params(struct generic_seq_args *args, struct registration_args *regargs);
#endif
