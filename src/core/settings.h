#ifndef _SIRIL_SETTING_H
#define _SIRIL_SETTING_H

#include <glib.h>

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
	XTRANS_FILTER,
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

/* This structure is used for storing all parameters used in photometry module */
struct phot_config {
	double gain;		// A/D converter gain in electrons per ADU
	double inner;		// Inner radius of the annulus used to measure local background.
	double outer;		// Outer radius of the annulus used to measure local background.
	double aperture;	// flux aperture
	gboolean force_radius;	// force the aperture radius value
	double minval, maxval;	// consider pixels outside this range as invalid for photometry
				// minval and maxval are stored as int, but adapted to image type
				// when used, so normalized to 1 for float, hence the double type
};

struct analysis_config {
	int mosaic_panel;
};

struct debayer_config {
	gboolean open_debayer;			// debayer images being opened
	gboolean use_bayer_header;		// use the pattern given in the file header
	sensor_pattern bayer_pattern;		// user-defined Bayer pattern
	interpolation_method bayer_inter;	// interpolation method for non-libraw debayer
	gboolean top_down;			// debayer top-down orientation
	int xbayeroff, ybayeroff;		// x and y Bayer offsets
};

// GUI data backup
struct stack_config {
	int method;		// 0=sum, 1=median, 2=average, 3=pixel max, 4=pixel min
	int normalisation_method;
	int rej_method;
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

struct gui_config {
	gboolean first_start;		// use to display information at first use
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
	gboolean warn_script_run; // show the notice when starting a script

	gboolean show_thumbnails; // show thumbnails in open dialog box
	gint thumbnail_size;

	int position_compass;	// compass position, can be moved
	gboolean catalog[7];	// 6 catalogs and 1 user catalog for annotations

	gint selection_guides;	// number of elements of the grid guides
				// (2 for a simple cross, 3 for the 3 thirds rule, etc.)

	// single registration GUI variable
	int reg_settings;	// selected registration method
	int reg_interpolation; // selected interpolation method
	GSList *pm_presets; // list of pixel math presets
	int default_rendering_mode; // Default view STF to use at startup
	int display_histogram_mode; // Default histogram view to use at startup
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
	gchar *bias_synth;
	gboolean use_bias_synth;
	gchar *dark_lib;
	gboolean use_dark_lib;
	gchar *flat_lib;
	gboolean use_flat_lib;
};

typedef struct {
	int radius;
	int adj_radius;
	gboolean adjust;
	double sigma;
	double roundness;
	double focal_length;
	double pixel_size_x;
	int convergence;
	gboolean relax_checks;
} star_finder_params;

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

	enum { RATIO, AMOUNT } mem_mode; // mode of memory management
	double memory_ratio;		// ratio of available memory to use for stacking (and others)
	double memory_amount;		// amount of memory in GB to use for stacking (and others)

	int hd_bitdepth; // Default bit depth for HD AutoStretch

	gboolean script_check_requires;	// check the requires command in scripts
	gboolean pipe_check_requires;	// check the requires command in pipes

	gboolean check_update;	// check update at startup

	gchar *lang;		// string value of the combobox lang

	gchar *swap_dir;	// swap directory

	// TODO: do we actually need these two?
	gdouble focal;		// focal length saved in config file
	gdouble pitch;		// pixel pitch saved in config file

	int wcs_formalism;	// formalism used in FITS header
	gchar *catalogue_paths[4]; // local star catalogues for plate solving and PCC

	gboolean rgb_aladin;	// Add CTYPE3='RGB' in the FITS header
	gchar *copyright;	// User copyright when saving image as TIFF

	gchar *starnet_dir;	// Location of starnet++ installation (requires v2.0.2 or greater)

	star_finder_params starfinder_conf;
	struct prepro_config prepro;
	struct gui_config gui;
	struct debayer_config debayer;
	struct phot_config phot_set;
	struct analysis_config analysis;
	struct stack_config stack;
	struct comp_config comp;
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

int print_settings_key(const char *group, const char *key, gboolean with_details);
int print_all_settings(gboolean with_details);

void free_preferences(preferences *pref);	// TODO check if they're used
void initialize_default_settings();

void update_gain_from_gfit();

#endif
