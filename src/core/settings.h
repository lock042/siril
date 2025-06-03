#ifndef _SIRIL_SETTING_H
#define _SIRIL_SETTING_H

#include <glib.h>
#include <lcms2.h>

/********************* S E T T I N G S *********************
 * Settings are variables used at runtime but kept across several lunches of
 * siril, but serializing their values to a file, historically called the init
 * file in siril, previously managed with libconfig, now with glib's key-value
 * file.
 * The structures below store these values for in-memory runtime access. They
 * are loaded from the init file at start-up and saved to it when required.
 */


/* first, types used here but that cannot be defined in siril.h because of the loop dependency */
typedef struct {
	int x, y, w, h;
} rectangle;

enum {
	WCS_FORMALISM_2,
	WCS_FORMALISM_1
};

/* TYPE OF NORMALIZATION */
typedef enum {
	NO_NORM,
	ADDITIVE,
	MULTIPLICATIVE,
	ADDITIVE_SCALING,
	MULTIPLICATIVE_SCALING,
} normalization;

/* TYPE OF SIGMA CLIPPING */
typedef enum {
	NO_REJEC,
	PERCENTILE,
	SIGMA,
	MAD,
	SIGMEDIAN,
	WINSORIZED,
	LINEARFIT,
	GESDT
} rejection;

typedef enum {
	BAYER_FILTER_RGGB,
	BAYER_FILTER_BGGR,
	BAYER_FILTER_GBRG,
	BAYER_FILTER_GRBG,
	XTRANS_FILTER_1,
	XTRANS_FILTER_2,
	XTRANS_FILTER_3,
	XTRANS_FILTER_4,
	BAYER_FILTER_NONE = -1		//case where pattern is undefined or untested
} sensor_pattern;
#define BAYER_FILTER_MIN BAYER_FILTER_RGGB
#define BAYER_FILTER_MAX BAYER_FILTER_GRBG

typedef enum {
	BAYER_BILINEAR,
	BAYER_VNG,
	BAYER_AHD,
	BAYER_AMAZE,
	BAYER_DCB,
	BAYER_HPHD,
	BAYER_IGV,
	BAYER_LMMSE,
	BAYER_RCD,
	BAYER_SUPER_PIXEL,
	XTRANS
} interpolation_method;

typedef enum {
	PSF_GAUSSIAN,
	PSF_MOFFAT_BFREE,
	PSF_MOFFAT_BFIXED
} starprofile;

typedef enum {
	NIL = 0,
	V2 = 1,
	V1MONO = 2,
	V1RGB = 4,
	TORCH = 8
} starnet_version;

typedef enum {
	MMB_ZOOM_FIT,
	MMB_ZOOM_100,
	MMB_ZOOM_TOGGLE
} mmb_action_t;

typedef enum {
	CMF_1931_2DEG = 0,
	CMF_1964_10DEG = 1
} cmf_pref;

typedef enum {
	ROW_ORDER_HEADER_TOPDOWN,
	ROW_ORDER_HEADER_BOTTOMUP,
	ROW_ORDER_FORCE_TOPDOWN,
	ROW_ORDER_FORCE_BOTTOMUP
} row_order_t;

/***********************************************************************************************/


/* This structure is used for all the elements in the box libraw_settings.
 * Don't forget to update conversion.c:initialize_libraw_settings() data when
 * modifying the glade settings */
struct libraw_config {
	double mul[3], bright;					// Color & brightness adjustement mul[0] = red, mul[1] = green = 1, mul[2] = blue
	int auto_mul, use_camera_wb, use_auto_wb;		// White Balance parameters
	int user_qual;						// Index of the Matrix interpolation set in dcraw, 0: bilinear, 1: VNG, 2: PPG, 3: AHD
	int user_black;						// black point correction
	double gamm[3];						// Gamma correction
};

struct astrometry_config {
	// common to siril and asnet solvers
	gboolean update_default_scale;	// update default focal length and pixel size from the result
	int percent_scale_range;	// percent below and above the expected sampling to allow
	int sip_correction_order;	// degrees of the polynomial correction
	double radius_degrees;		// radius around the target coordinates (degrees)
	int max_seconds_run;		// maximum seconds of CPU time to try solving
	/* for local astrometry.net  only */
	gboolean keep_xyls_files;	// do not delete .xyls FITS tables
	gboolean keep_wcs_files;	// do not delete .wcs result files
	gboolean show_asnet_output;	// show solve-field output in main log
	gchar* default_obscode;	// default observatory code
};


/* This structure is used for storing all parameters used in photometry module */
struct phot_config {
	double gain;		// A/D converter gain in electrons per ADU
	double inner;		// Inner radius of the annulus used to measure local background.
	double outer;		// Outer radius of the annulus used to measure local background.
	double auto_inner_factor;// factor for automatic inner radius computation from FWHM
	double auto_outer_factor;// factor for automatic outer radius computation from FWHM
	double aperture;	// flux aperture
	double auto_aperture_factor;	// ratio between the aperture and the half-FWHM, used with the dynamic aperture option
	gboolean force_radius;	// force the aperture radius value
	double minval, maxval;	// consider pixels outside this range as invalid for photometry
				// minval and maxval are stored as int, but adapted to image type
				// when used, so normalized to 1 for float, hence the double type
	int discard_var_catalogues; 	// encodes the catalogues to be used to discard the variable stars from the comparison stars list
				// consider this integer in its binary form:
				// b0 (LSB) sets (or not) the VSX catologue
				// b1 sets (or not) the GCVS catologue
				// b2 sets (or not) the Varisum catologue
};

struct analysis_config {
	int mosaic_panel;
	int mosaic_window;
};

struct debayer_config {
	gboolean open_debayer;			// debayer images being opened
	gboolean use_bayer_header;		// use the pattern given in the file header
	sensor_pattern bayer_pattern;		// user-defined Bayer pattern
	interpolation_method bayer_inter;	// interpolation method for non-libraw debayer
	row_order_t orientation;			// orientation preference
	int xbayeroff, ybayeroff;		// x and y Bayer offsets
	int xtrans_passes;			// number of passes for X-Trans debayer
};

// GUI data backup
struct stack_config {
	int method;		// 0=sum, 1=median, 2=average, 3=pixel max, 4=pixel min
	int normalisation_method;
	int rej_method;
	int weighting_method;
	double sigma_low, sigma_high;
	double linear_low, linear_high;
	double percentile_low, percentile_high;
};

struct comp_config {
	gboolean fits_enabled;
	int fits_method;		// 0=Rice, 1=GZIP1, 2=GZIP2, 3=Hcompress
	double fits_quantization;	// quantization factor for floating point compression
	double fits_hcompress_scale;	// scale factor for Hcompress compression
};

struct editor_config {
	gboolean highlight_syntax;
	gboolean highlight_bracketmatch;
	gboolean rmargin;
	int rmargin_pos;
	gboolean show_linenums;
	gboolean show_linemarks;
	gboolean highlight_currentline;
	gboolean autoindent;
	gboolean indentontab;
	gboolean smartbs;
	gboolean smarthomeend;
	gboolean showspaces;
	gboolean shownewlines;
	gboolean minimap;
};

typedef enum {
	ROI_MANUAL,
	ROI_AUTO
} roi_mode_t;

typedef struct configurable_colors_param {
	gchar* color_bkg_samples;
	gchar* color_std_annotations;
	gchar* color_dso_annotations;
	gchar* color_sso_annotations;
	gchar* color_tmp_annotations;
} configurable_colors;

struct mouse_config {
	GSList* mouse_actions_array;
	GSList* scroll_actions_array;
};

struct gui_config {
	gchar *first_start;		// use to display information at first use
	gboolean silent_quit;
	gboolean silent_linear;
	gboolean remember_windows;	// restore windows at their previous location
	rectangle main_w_pos;
	gint pan_position;
	gboolean is_extended;
	gboolean is_maximized;

	gint combo_theme;	// index of the combobox theme
	gdouble font_scale;	// font scale
	gboolean icon_symbolic;	// icon style

	GSList *script_path;	// script path directories
	gboolean warn_scripts_run; // show the notice when starting a script
							// updated with new warning wording to force redisplay

	gboolean show_thumbnails; // show thumbnails in open dialog box
	gint thumbnail_size;

	int position_compass;	// compass position, can be moved
	gboolean catalog[11];	// 8 system catalogs and 2 user catalogs for annotations and 1
				// short-lived catalogue for "who's in the field" annotations
				// see also cat in annotation_catalogues.c

	gint selection_guides;	// number of elements of the grid guides
				// (2 for a simple cross, 3 for the 3 thirds rule, etc.)

	gboolean show_deciasec;	// show tenths of arcseconds in the displayed hover coordinates

	// single registration GUI variable
	int reg_settings;	// selected registration method
	int reg_interpolation; // selected interpolation method
	gboolean reg_clamping; // clamping for Lanczos/Cubic
	GSList *pm_presets; // list of pixel math presets
	int default_rendering_mode; // Default view STF to use at startup
	int display_histogram_mode; // Default histogram view to use at startup
	roi_mode_t roi_mode; // Whether to set the ROI manually or auto from selection
	gboolean enable_roi_warning; // Whether to notify when a ROI-enabled dialog starts
	configurable_colors config_colors; // This used to configure some colors in Siril
	mmb_action_t mmb_action; // Defines middle mouse button double click behaviour
	double mouse_speed_limit; // Defines a mximum step size for GDK_SCROLL_SMOOTH actions
	struct mouse_config mouse_cfg; // String representation of mouse & scroll actions
	struct editor_config editor_cfg; // Configuration for the script editor
};

// TODO: is any of the following used for something else than providing the default GUI value?
struct prepro_config {
	gboolean cfa;		// type of sensor for cosmetic correction in preprocessing
	gboolean equalize_cfa;	// if flat will be equalized in preprocessing
	gboolean fix_xtrans;	// Use to fix xtrans darks and biases with the AF square
	rectangle xtrans_af;	// if no xtrans model found, use these values
	rectangle xtrans_sample;// if no xtrans model found, use these values
	gchar *bias_lib;
	gboolean use_bias_lib;
	gchar *dark_lib;
	gboolean use_dark_lib;
	gchar *flat_lib;
	gboolean use_flat_lib;
	gchar *disto_lib;
	gboolean use_disto_lib;
	gchar *stack_default;
	gboolean use_stack_default;
};

typedef struct {
	int radius;
	double sigma;
	double roundness;
	double focal_length;
	double pixel_size_x;
	int convergence;
	gboolean relax_checks;
	starprofile profile;
	double min_beta;
	double min_A, max_A;
	double max_r;
} star_finder_params;

typedef struct fftw_params {
	double timelimit;
	int strategy;
	gboolean multithreaded;
	gchar* wisdom_file;
	int fft_cutoff;
} fftw_params;

typedef enum {
	TYPE_SRGB,
	TYPE_REC2020,
	TYPE_CUSTOM
} working_gamut_type;

typedef enum {
	EXPORT_SRGB,
	EXPORT_WORKING,
	EXPORT_IMAGE_ICC
} export_icc_type;

typedef enum {
	ICC_ASSIGN_NEVER = 0,
	ICC_ASSIGN_ON_LOAD = 1,
	ICC_ASSIGN_ON_STACK = 2,
	ICC_ASSIGN_ON_STRETCH = 4,
	ICC_ASSIGN_ON_COMPOSITION = 8
} icc_assign_type;

typedef enum {
	ICC_NEVER_AUTOCONVERT,
	ICC_ASK_TO_CONVERT,
	ICC_ALWAYS_AUTOCONVERT
} icc_autoconvert_type;

typedef struct icc_params {
	cmsUInt32Number rendering_intent;
	cmsUInt32Number proofing_intent;
	cmsUInt32Number export_intent;
	cmsUInt32Number processing_intent;
	working_gamut_type working_gamut;
	gchar *icc_path_monitor;
	gchar *icc_path_soft_proof;
	gboolean custom_monitor_profile_active;
	gboolean soft_proofing_profile_active;
	gchar* custom_icc_trc;
	gchar* custom_icc_gray;
	export_icc_type export_8bit_method;
	export_icc_type export_16bit_method;
	gboolean default_to_srgb;
	gboolean rendering_bpc;
	gboolean proofing_bpc;
	icc_autoconvert_type autoconversion;
	icc_assign_type autoassignment;
	gboolean pedantic_linear;
	cmf_pref cmf;
} icc_params;

struct spcc_favourites {
	gboolean use_spcc_repository;
	gboolean auto_spcc_update; // automatically update spcc repository at startup
	gchar *redpref;
	gchar *greenpref;
	gchar *bluepref;
	gchar *lpfpref;
	gchar *oscfilterpref;
	gchar *monosensorpref;
	gchar *oscsensorpref;
	gboolean is_mono;
	gboolean is_dslr;
	gboolean nb_mode;
	double red_wl;
	double green_wl;
	double blue_wl;
	double red_bw;
	double green_bw;
	double blue_bw;
};

/**
 * This is the preference structure.
 * WARNING!!
 * If you update something in this structure you absolutely need to
 * update pref_init and the table in settings.c
 * All strings are allocated g_char *, to be freed with g_free().
 */
struct pref_struct {
	gchar *wd;		// saved working directory, not impacted by scripts like com.wd
	gchar *ext;		// FITS extension used in SIRIL
	gboolean force_16bit;	// don't use 32 bits for pixel depth
	gboolean allow_heterogeneous_fitseq; // allow images in FITS cubes to have different sizes
	gboolean fits_save_icc;		// save ICC profiles in FITS

	enum { RATIO, AMOUNT } mem_mode; // mode of memory management
	double memory_ratio;		// ratio of available memory to use for stacking (and others)
	double memory_amount;		// amount of memory in GB to use for stacking (and others)

	int hd_bitdepth; // Default bit depth for HD AutoStretch

	gboolean script_check_requires;	// check the requires command in scripts
	gboolean pipe_check_requires;	// check the requires command in pipes

	gboolean check_update;	// check update at startup

	gchar *lang;		// string value of the combobox lang

	gchar *swap_dir;	// swap directory

	gboolean binning_update;// update pixel size of binned images

	int wcs_formalism;	// formalism used in FITS header
	gchar *catalogue_paths[6]; // local star catalogues for plate solving and PCC

	gboolean rgb_aladin;	// Add CTYPE3='RGB' in the FITS header
	gboolean use_checksum;  // Verify checksum in FITS header
	gchar *copyright;	// User copyright when saving image as TIFF

	gchar *starnet_exe;	// Location of starnet++ executable
	gchar *starnet_weights;	// Location of StarNet weights file (optional, Torch based StarNet only)
	gchar *asnet_dir;	// Location of solve-field or asnet-ansvr installation on Windows
	gchar *graxpert_path; // Location of GraXpert executable

	star_finder_params starfinder_conf;
	struct prepro_config prepro;
	struct gui_config gui;
	struct debayer_config debayer;
	struct phot_config phot_set;
	struct astrometry_config astrometry;
	struct analysis_config analysis;
	struct stack_config stack;
	struct comp_config comp;
	struct spcc_favourites spcc;
	fftw_params fftw_conf;
	int max_slice_size; // Used when processing img_t in slices to limit the wisdom required
	icc_params icc;
	GList *selected_scripts;
	gboolean use_scripts_repository;
	gboolean auto_script_update; // automatically update scripts repository at startup
};

typedef struct pref_struct preferences;
/**
 * End of preference structure. Read above if modified.
 */

enum settings_type {
	STYPE_BOOL,
	STYPE_INT,
	STYPE_DOUBLE,
	STYPE_STR,
	STYPE_STRDIR,	// not native, for path checks
	STYPE_STRLIST	// not native, for script paths
	// long and ulong also possible
};

struct range_int_s {
	int min, max;
};

struct range_double_s {
	double min, max;
};

struct settings_access {
	const char *group;	// group name
	const char *key;	// key name
	enum settings_type type;// type of the option
	const char *desc;	// short description
	void *data;		// pointer to the data in com.pref
	union {			// allowed values range
		struct range_int_s range_int;
		struct range_double_s range_double;
	};
};

struct settings_access *get_all_settings();
struct settings_access *get_key_settings(const char *group, const char *key);

gchar* get_settings_key(const char *group, const char *key, gboolean with_details);
int print_settings_key(const char *group, const char *key, gboolean with_details);
int print_all_settings(gboolean with_details);

void free_preferences(preferences *pref);	// TODO check if they're used
void initialize_default_settings();
void set_wisdom_file();

void update_gain_from_gfit();

#endif
