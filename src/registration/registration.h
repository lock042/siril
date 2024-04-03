#ifndef _REGISTRATION_H_
#define _REGISTRATION_H_

#include "core/siril.h"
#include "algos/PSF.h"
#include "core/processing.h"

#define NUMBER_OF_METHODS 8
#define MAX_DISTO_SIZE 7 // need to duplicate MAX_SIP_SIZE here because of circular refs with opencv

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
	PLANETARY_FULLDISK, PLANETARY_SURFACE
} planetary_type;

typedef enum {
	REG_PAGE_GLOBAL,
	REG_PAGE_COMET,
	REG_PAGE_3_STARS,
	REG_PAGE_KOMBAT,
	REG_PAGE_APPLYREG,
	REG_PAGE_MISC
} reg_notebook_page;

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

/* arguments passed to registration functions */
struct registration_args {
	registration_function func;	// the registration function
	sequence *seq;			// the sequence to register
	int reference_image;		// reference image index
	struct seq_filter_config filters; // parsed image filters (.filter_included always used)
	int layer;			// layer of images on which the registration is computed
	int retval;			// retval of func
	gboolean run_in_thread;		// true if the registration was run in a thread
	gboolean follow_star;		// follow star position between frames
	gboolean matchSelection;	// Match stars found in the seleciton of reference image
	rectangle selection;		// the selection rectangle
	gboolean x2upscale;		// apply an x2 upscale for pseudo drizzle
	gboolean cumul;			// cumul reg data with previous one
	int min_pairs;			// Minimum number of star pairs for success
	int max_stars_candidates;	// Max candidates after psf fitting for global reg
	transformation_type type;	// Use affine transform  or homography
	float percent_moved;		// for KOMBAT algorithm
	gboolean two_pass;		// use the two-pass computation to find a good ref image
	seq_image_filter filtering_criterion; // the filter, (seqapplyreg only)
	double filtering_parameter;	// and its parameter (seqapplyreg only)
	gboolean no_starlist;		// disable star list creation (2pass only)
	float astrometric_scale;		// scaling factor (for mosaic only)
	gboolean undistort;		// apply undistorsion with SIP data

	/* data for generated sequence, for star alignment/mosaic registration */
	gboolean no_output;		// write transformation to .seq
	int new_total;                  // remaining images after registration
	imgdata *imgparam;		// imgparam for the new sequence
	regdata *regparam;		// regparam for the new sequence
	char *prefix;		// prefix of the created sequence if any
	gboolean load_new_sequence;	// load the new sequence if success
	gchar *new_seq_name;
	opencv_interpolation interpolation; // type of rotation interpolation
	framing_type framing;		// used by seqapplyreg to determine framing
	gboolean clamp;				// should Bicubic and Lanczos4 interpolation be clamped?
	double clamping_factor;		// used to set amount of interpolation clamping
	opencv_projector projector; // used by mosaic registration
};

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
} astrometric_roi;

typedef struct {
	double AP[MAX_DISTO_SIZE][MAX_DISTO_SIZE];
	double BP[MAX_DISTO_SIZE][MAX_DISTO_SIZE];
	int order;
	double xref, yref;
} disto_data;

struct astrometric_args{
	int nb;
	int refindex;
	Homography *Rs;
	Homography *Ks;
	disto_data *disto;
	float scale;
	astrometric_roi roi;
};
struct registration_method *new_reg_method(const char *name, registration_function f,
		selection_type s, registration_type t); // for compositing
void initialize_registration_methods();
struct registration_method * get_selected_registration_method();
int register_shift_dft(struct registration_args *args);
int register_shift_fwhm(struct registration_args *args);
int register_star_alignment(struct registration_args *args);
int register_multi_step_global(struct registration_args *regargs);
int register_comet(struct registration_args *regargs);
int register_3stars(struct registration_args *regargs);
int register_apply_reg(struct registration_args *regargs);
int register_kombat(struct registration_args *args);
int register_manual(struct registration_args *regargs); // defined in compositing/compositing.c
int register_astrometric(struct registration_args *regargs);

void reset_3stars();
int _3stars_check_registration_ready();
gboolean _3stars_check_selection();

pointf get_velocity();
void update_reg_interface(gboolean dont_change_reg_radio);
void compute_fitting_selection(rectangle *area, int hsteps, int vsteps, int preserve_square);
void get_the_registration_area(struct registration_args *reg_args, const struct registration_method *method); // for compositing
gpointer register_thread_func(gpointer p);

/** getter */
int get_registration_layer(const sequence *seq);
int seq_has_any_regdata(const sequence *seq); // same as get_registration_layer but does not rely on GUI for com.seq

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

regdata *star_align_get_current_regdata(struct registration_args *regargs);
int star_align_prepare_results(struct generic_seq_args *args);
int star_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_, int threads);
int star_align_finalize_hook(struct generic_seq_args *args);
int star_match_and_checks(psf_star **ref_stars, psf_star **stars, int nb_ref_stars, int nb_stars, struct registration_args *regargs, int filenum, Homography *H);

const char *describe_transformation_type(transformation_type type);

void selection_H_transform(rectangle *selection, Homography Href, Homography Himg);
void guess_transform_from_seq(sequence *seq, int layer,
		transformation_type *min, transformation_type *max, gboolean excludenull);
transformation_type guess_transform_from_H(Homography H);
gboolean check_before_applyreg(struct registration_args *regargs);
gboolean layer_has_registration(const sequence *seq, int layer);
gboolean layer_has_usable_registration(sequence *seq, int layer);
int get_first_selected(sequence *seq);
regdata *apply_reg_get_current_regdata(struct registration_args *regargs);

void translation_from_H(Homography H, double *dx, double *dy);
Homography H_from_translation(double dx, double dy);
void SetNullH(Homography *H);
int shift_fit_from_reg(fits *fit, Homography H);

int minidx(const float *arr, const gboolean *mask, int nb, float *val);

void free_astrometric_args(struct astrometric_args *astargs);

#endif
