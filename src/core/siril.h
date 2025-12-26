#ifndef SIRIL_H
#define SIRIL_H
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <glib.h>
#include <glib/gstdio.h>
#include <glib/gprintf.h>
#include <gtk/gtk.h>
#include <gsl/gsl_histogram.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <libintl.h>

#include <fitsio.h>	// fitsfile

#include "core/settings.h"

#define _(String) gettext (String)
#define gettext_noop(String) String
#define N_(String) gettext_noop (String)

#ifdef SIRIL_OUTPUT_DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define START_TIMER struct timeval t_start, t_end ; gettimeofday(&t_start, NULL)
#define END_TIMER gettimeofday(&t_end, NULL); show_time(t_start, t_end)

#if defined (HAVE_FFTW3F_OMP) || defined (HAVE_FFTW3F_THREADS)
#define HAVE_FFTW3F_MULTITHREAD
#endif

/* https://stackoverflow.com/questions/1644868/define-macro-for-debug-printing-in-c */
#define siril_debug_print(fmt, ...) \
	do { if (DEBUG_TEST) fprintf(stdout, fmt, ##__VA_ARGS__); } while (0)

#define PRINT_ALLOC_ERR fprintf(stderr, "Out of memory in %s (%s:%d) - aborting\n", __func__, __FILE__, __LINE__)
#define PRINT_ANOTHER_THREAD_RUNNING siril_log_message(_("Another task is already in progress, ignoring new request.\n"))

#ifndef __cplusplus
// excluding from C++ because it conflicts with stl
#undef max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#undef min
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


#define SWAPD(a,b) { double temp = (a); (a) = (b); (b) = temp; }

#define SQR(x) ((x)*(x))
#endif
#define RADCONV (((3600.0 * 180.0) / M_PI) / 1.0E3)

// Used for sanity checking reported sizes from image files.
// Not a guarantee that the file will fit in memory
#define MAX_IMAGE_DIM 100000

#define USHRT_MAX_DOUBLE ((double)USHRT_MAX)
#define SHRT_MAX_DOUBLE ((double)SHRT_MAX)
#define USHRT_MAX_SINGLE ((float)USHRT_MAX)
#define SHRT_MAX_SINGLE ((float)SHRT_MAX)
#define UCHAR_MAX_DOUBLE ((double)UCHAR_MAX)
#define UCHAR_MAX_SINGLE ((float)UCHAR_MAX)
#define INV_USHRT_MAX_SINGLE .000015259022f	 // 1/65535
#define INV_UCHAR_MAX_SINGLE .0039215686f	 // 1/255
#define INV_USHRT_MAX_DOUBLE .000015259022	 // 1/65535
#define INV_UCHAR_MAX_DOUBLE .0039215686	 // 1/255

#define BYTES_IN_A_MB 1048576	// 1024*1024

#define SEQUENCE_DEFAULT_INCLUDE TRUE	// select images by default

typedef unsigned char BYTE;	// default type for image display data
typedef unsigned short WORD;	// default type for internal image data

#define MAX_REF_STARS 20	// uchar
#define MAX_SEQPSF MAX_REF_STARS// max number of stars for which seqpsf can be run
				// if we use seqpsf for the light curve, we must allow as many as MAX_REF_STARS

#define CMD_HISTORY_SIZE 50	// size of the command line history

#define ZOOM_MAX	128
#define ZOOM_MIN	0.03125
#define ZOOM_IN		1.5
#define ZOOM_OUT	1.0 / ZOOM_IN
#define ZOOM_NONE	1.0
#define ZOOM_FIT	-1.0	// or any value < 0
#define ZOOM_DEFAULT	ZOOM_FIT

/* Some statistic constants */
#define AVGDEV_NORM 1.2533
#define MAD_NORM 1.4826
#define BWMV_NORM 0.9901
#define PBMV_NORM 0.9709
#define SN_NORM 1.1926
#define QN_NORM 2.2191

typedef struct _SirilDialogEntry SirilDialogEntry;

/* used for open and savedialog */
typedef GtkWidget SirilWidget;

#if (defined _WIN32) || (defined(OS_OSX))
#define SIRIL_EOL "\r\n"
#else
#define SIRIL_EOL "\n"
#endif

typedef enum {
	RESPONSE_CANCEL = 1,
	RESPONSE_CLIP,
	RESPONSE_RESCALE_CLIPNEG,
	RESPONSE_RESCALE_ALL
} OverrangeResponse;

typedef enum {
	TYPEUNDEF = (1 << 1),
	TYPEFITS = (1 << 2),
	TYPETIFF = (1 << 3),
	TYPEBMP = (1 << 4),
	TYPEPNG = (1 << 5),
	TYPEJPG = (1 << 6),
	TYPEHEIF = (1 << 7),
	TYPEPNM = (1 << 8),
	TYPEPIC = (1 << 9),
	TYPERAW = (1 << 10),
	TYPEXISF = (1 << 11),
	TYPEJXL = (1 << 12),
	TYPEAVIF = (1 << 13),
	TYPEAVI = (1 << 20),
	TYPESER = (1 << 21)
} image_type;

/* indices of the image data layers */
#define BW_LAYER 	0
#define RLAYER		0
#define GLAYER		1
#define BLAYER		2
#define RGB_LAYER 	3

/* indices of the viewports (graphical elements) */
#define BW_VPORT 	0
#define RED_VPORT 	0
#define GREEN_VPORT 	1
#define BLUE_VPORT 	2
#define RGB_VPORT 	3
#define MAXGRAYVPORT 	3	// 3 gray vports supported only (R, G, B)
#define MAXCOLORVPORT	1	// 1 color vport supported only (RGB)
#define MAXVPORT 	MAXGRAYVPORT + MAXCOLORVPORT

/* defines for copyfits actions */
#define CP_INIT		0x01	// initialize data array with 0s
#define CP_ALLOC	0x02	// reallocs data array
#define CP_COPYA	0x04	// copy data array content
#define CP_FORMAT	0x08	// copy metadata
#define CP_EXPAND	0x20	// expands a one-channel to a three channels

#define PREVIEW_NB 2

/* special values for com.seq.current, the currently loaded image of the
 * sequence. Negative values can be used to indicate that a specific image is
 * loaded while there is a sequence or not in com.seq */
#define RESULT_IMAGE -1		// used as current image index in a sequence
				// when the result of a processing is displayed
#define UNRELATED_IMAGE -2	// image loaded while a sequence was loaded too
#define SCALED_IMAGE -3		// the single image has a different size than
				// the loaded sequence

#define MAX_STARS 200000	// maximum length of com.stars
#define MAX_STARS_FITTED 2000	// maximum number of stars fitted for registration
#define MIN_STARS_FITTED 100	// minimum number of stars fitted for registration
#define DEF_BOX_RADIUS 5	// default radius of the box in starfinder_conf

#define INDEX_MAX 65535		// maximum index for images

typedef struct sequ sequence;
typedef struct ffit fits;
typedef struct guiinf guiinfo;
typedef struct cominf cominfo;
typedef struct historic_struct historic;
typedef struct fwhm_struct psf_star;
typedef struct photometry_struct photometry;
typedef struct tilt_struct sensor_tilt;

typedef struct {
	double x, y;
} point;

typedef struct {
	float x, y;
} pointf;

typedef struct {
	int x, y;
} pointi;

/* global structures */

typedef enum {
	BOOL_NOT_SET = -1,
	BOOL_FALSE = 0,
	BOOL_TRUE = 1,
} super_bool;

/* same order as in the combo box 'comboExport' */
typedef enum {
	EXPORT_FITS,
	EXPORT_FITSEQ,
	EXPORT_TIFF,
	EXPORT_SER,
	EXPORT_AVI,
	EXPORT_MP4,
	EXPORT_MP4_H265,
	EXPORT_WEBM_VP9
} export_format;

typedef enum {
	LINEAR_DISPLAY,
	LOG_DISPLAY,
	SQRT_DISPLAY,
	SQUARED_DISPLAY,
	ASINH_DISPLAY,
	STF_DISPLAY,
	HISTEQ_DISPLAY
} display_mode;
#define DISPLAY_MODE_MAX HISTEQ_DISPLAY

typedef enum {
	NORMAL_COLOR,
	RAINBOW_COLOR
} color_map;

typedef enum {
	CLIP,
	RESCALE,
	RGBBLEND,
	RESCALEGLOBAL
} clip_mode_t;

typedef enum {
	MIPSLOHI,
	MINMAX,
	USER
} sliders_mode;

typedef enum {
	OPEN_IMAGE_ERROR = -1,
	OPEN_IMAGE_OK = 0,
	OPEN_IMAGE_CANCEL = 10,
} open_image_status;

typedef enum {
	OPENCV_NEAREST = 0,
	OPENCV_LINEAR = 1,
	OPENCV_CUBIC = 2,
	OPENCV_AREA = 3,
	OPENCV_LANCZOS4 = 4,
	OPENCV_NONE = 5 // this one will use the pixel-wise shift transform w/o opencv
} opencv_interpolation;


typedef enum {
	SEQ_REGULAR = 0,
	SEQ_SER = 1,
	SEQ_FITSEQ = 2,
#ifdef HAVE_FFMS2
	SEQ_AVI = 3,
#endif
	SEQ_INTERNAL = 4
} sequence_type;

typedef enum {
	MULTI_THREADED = 0,	// equivalent to com.max_threads
	SINGLE_THREADED = 1
	/* and higher values are the number of threads */
} threading_type;

typedef enum {
	RICE_COMP,
	GZIP1_COMP,
	GZIP2_COMP,
	HCOMPRESS_COMP
} compression_mode;

typedef enum {
	EXT_NONE,
	EXT_STARNET,
	EXT_ASNET,
	EXT_GRAXPERT,
	EXT_PYTHON,
	INT_PROC_THREAD
} external_program;

typedef enum {
	SCALING_NONE,
	SCALING_HA_UP,
	SCALING_OIII_DOWN
} extraction_scaling;

typedef enum {
	DISTO_UNDEF, // No distortion
	DISTO_IMAGE, // Distortion from current image
	DISTO_FILE,  // Distortion from given file
	DISTO_MASTER, // Distortion from master files
	DISTO_FILES, // Distortion stored in each file (true only from seq platesolve, even with no distortion, it will be checked upon reloading)
	DISTO_FILE_COMET // special for cometary alignment, to be detected by apply reg. Enables to
} disto_source;

// defined in src/io/pythonmodule.h
typedef struct _Connection Connection;

/* image data, exists once for each image */
typedef struct {
	int filenum;		/* real file index in the sequence, i.e. for mars9.fit = 9 */
	gboolean incl;		/* selected in the sequence, included for future processings? */
	GDateTime *date_obs;	/* date of the observation, processed and copied from the header */
	double airmass;     /* airmass of the image, used in photometry */
	int rx, ry;
} imgdata;

/* this structure is used to characterize the statistics of the image */
typedef struct {
	long total,	// number of pixels
	     ngoodpix;	// number of non-zero pixels
	double mean, median, sigma, avgDev, mad, sqrtbwmv,
	       location, scale, min, max, normValue, bgnoise;
	gint _nb_refs;	// reference counting for data management
} imstats;

/* this structure is used to characterize the statistics of the overalps in a sequence of images */
typedef struct {
	int i, j;
	size_t Nij;
	rectangle areai, areaj;
	float medij, medji;
	float madij, madji;
	float locij, locji;
	float scaij, scaji;
} overlap_stats_t;

typedef struct {
	double h00, h01, h02;
	double h10, h11, h12;
	double h20, h21, h22;
	int pair_matched;
	int Inliers;
} Homography;

/* registration data, exists once for each image and each layer */
typedef struct {
	psf_star *fwhm_data;	// used in PSF/FWHM registration, not saved
	float fwhm;		// copy of fwhm->fwhmx, used as quality indicator, saved data
	float weighted_fwhm;	// used to exclude spurious images.
	float roundness;	// fwhm->fwhmy / fwhm->fwhmx, 0 when uninit, ]0, 1] when set
	double quality;
	float background_lvl;
	int number_of_stars;
	Homography H;
} regdata;

// to be stored in the seq file
typedef struct {
	disto_source index; // disto_source enum
	gchar *filename; // filename if DISTO_FILE
	pointf velocity; // shift velocity if DISTO_FILE_COMET
} disto_params;

/* see explanation about sequence and single image management in io/sequence.c */

struct sequ {
	char *seqname;		// name of the sequence, as in name.seq
	int number;		// number of images in the sequence
	int selnum;		// number of selected images in the sequence
	int fixed;		// fixed length of image index in filename (like %3d)
	int nb_layers;		// number of layers embedded in each image file, -1 if unknown
	unsigned int rx;	// first image width (or ref if set)
	unsigned int ry;	// first image height (or ref if set)
	gboolean is_variable;	// sequence has images of different sizes (imgparam->r[xy])
	gboolean is_drizzle; 	// sequence is a drizzle sequence, weights files are stored in ./drizzletmp
	int bitpix;		// image pixel format, from fits
	int reference_image;	// reference image for registration
	imgdata *imgparam;	// a structure for each image of the sequence
	regdata **regparam;	// *regparam[nb_layers], may be null if nb_layers is unknown
	imstats ***stats;	// statistics of the images for each layer, may be null too
	overlap_stats_t **ostats;	// statistics of the overlaps for each layer, may be null too
	/* in the case of a CFA sequence, depending on the opening mode, we cannot store
	 * and use everything that was in the seqfile, so we back them up here */
	regdata **regparam_bkp;	// *regparam[3], null if nothing to back up
	imstats ***stats_bkp;	// statistics of the images for 3 layers, may be null too
	disto_params *distoparam;	// the distortion parameters used for the registration if any, one per layer
	char* cached_ext;	// FIT sequences only: cache the ext to speed up filename search

	/* beg and end are used prior to imgparam allocation, hence their usefulness */
	int beg;		// imgparam[0]->filenum
	int end;		// imgparam[number-1]->filenum
	double exposure;	// exposure of frames (we assume they are all identical)
	gboolean fz;
	sequence_type type;
	struct ser_struct *ser_file;
	gboolean cfa_opened_monochrome;	// in case the CFA SER was opened in monochrome mode
	struct fits_sequence *fitseq_file; // FITS sequence data structure
#ifdef HAVE_FFMS2
	struct film_struct *film_file;
	const char *ext;	// extension of video, NULL if not video
#endif
	fits **internal_fits;	// for INTERNAL sequences: images references. Length: number
	fitsfile **fptr;	// file descriptors for open-mode operations
#ifdef _OPENMP
	omp_lock_t *fd_lock;	// locks for open-mode threaded operations
#endif

	int current;		// file number currently loaded in gfit (or displayed)

	/* registration previsualisation and manual alignment data */
	int previewX[PREVIEW_NB], previewY[PREVIEW_NB];	// center, -1 is uninitialized value
	int previewW[PREVIEW_NB], previewH[PREVIEW_NB];	// 0 is uninitialized value

	gboolean needs_saving;	// a dirty flag for the sequence, avoid saving it too often
	gboolean reg_invalidated; // a flag to detect if regframe can be plotted

	psf_star **photometry[MAX_SEQPSF+1];// psf for multiple stars for all images
	int reference_star;	// reference star for apparent magnitude (index of photometry)
	double reference_mag;	// reference magnitude for the reference star
	double photometry_colors[MAX_SEQPSF][3]; // colors for each photometry curve
};

/* this struct is used to manage data associated with a single image loaded, outside a sequence */
typedef struct {
	char *filename;		// the name of the file
	gboolean fileexist;// flag of existing file
	char *comment;		// comment on how the file got there (user load, result...)
	int nb_layers;		// number of layers embedded in each image file
	fits *fit;		// the fits is still gfit, but a reference doesn't hurt
} single;

typedef struct {
	char objctra[FLEN_VALUE];
	char objctdec[FLEN_VALUE];
	double ra;
	double dec;
	gboolean pltsolvd;
	char pltsolvd_comment[FLEN_COMMENT];
} wcs_info;

typedef struct {
	double norm[3];			// Normalization value
	char type[FLEN_VALUE];		// spectrum, phase
	char ord[FLEN_VALUE];		// regular, centered
} dft_info;

typedef enum { DATA_USHORT, DATA_FLOAT, DATA_UNSUPPORTED } data_type;

#define DEFAULT_DOUBLE_VALUE -999.0
#define DEFAULT_FLOAT_VALUE -999.f
#define DEFAULT_INT_VALUE -INT_MAX
#define DEFAULT_UINT_VALUE INT_MAX
#define DEFAULT_USHORT_VALUE 0

typedef struct {
	/* data obtained from the FITS file */
	double bscale, bzero;
	WORD lo;		// MIPS-LO key in FITS file, "Lower visualization cutoff"
	WORD hi;		// MIPS-HI key in FITS file, "Upper visualization cutoff"
	float flo, fhi; // Same but float format.
	char program[FLEN_VALUE];           // Software that created this HDU
	char filename[FLEN_VALUE];           // Original Filename
	double data_max;	// used to check if 32b float is in the [0, 1] range
	double data_min;	// used to check if 32b float is in the [0, 1] range
	double pixel_size_x, pixel_size_y;	// XPIXSZ and YPIXSZ keys
	unsigned int binning_x, binning_y;	// XBINNING and YBINNING keys
	char row_order[FLEN_VALUE];
	GDateTime *date, *date_obs;		// creation and acquisition UTC dates
	double mjd_obs;					// date-obs in Julian
	double expstart, expend;		// Julian dates
	char filter[FLEN_VALUE];		// FILTER key
	char image_type[FLEN_VALUE];		// IMAGETYP key
	char object[FLEN_VALUE];		// OBJECT key
	char instrume[FLEN_VALUE];		// INSTRUME key
	char telescop[FLEN_VALUE];		// TELESCOP key
	char observer[FLEN_VALUE];		// OBSERVER key
	double centalt, centaz;
	double sitelat, sitelong;		// SITE LAT and LONG as double
	char sitelat_str[FLEN_VALUE];		// SITE LATITUDE key as string
	char sitelong_str[FLEN_VALUE];		// SITE LONGITUDE key as string
	double siteelev;			// SITE ELEVATION key as double
	char bayer_pattern[FLEN_VALUE];		// BAYERPAT key Bayer Pattern if available
	int bayer_xoffset, bayer_yoffset;
	double airmass;                   // relative optical path length through atmosphere.
	/* data obtained from FITS or RAW files */
	double focal_length, flength, iso_speed, exposure, aperture, ccd_temp, set_temp;
	double livetime;		// sum of exposures (s)
	guint stackcnt;			// number of stacked frame
	double cvf;			// Conversion factor (e-/adu)
	int key_gain, key_offset;	// Gain, Offset values read in camera headers.
	char focname[FLEN_VALUE]; // focuser name
	int focuspos, focussz;
	double foctemp;

	/* Plate Solving data */
	wcs_info wcsdata;		// data from the header
	struct wcsprm *wcslib;		// struct of the lib

	/* data used in the Fourier space */
	dft_info dft;
} fkeywords;

struct ffit {
	unsigned int rx;	// image width	(naxes[0])
	unsigned int ry;	// image height	(naxes[1])
	int bitpix;		// current bitpix of loaded data
	int orig_bitpix;	// original bitpix of the file
	/* bitpix can take the following values:
	 * BYTE_IMG	(8-bit byte pixels, 0 - 255)
	 * SHORT_IMG	(16 bit signed integer pixels)
	 * USHORT_IMG	(16 bit unsigned integer pixels)	(used by Siril, a bit unusual)
	 * LONG_IMG	(32-bit integer pixels)
	 * FLOAT_IMG	(32-bit floating point pixels)
	 * DOUBLE_IMG	(64-bit floating point pixels)
	 * http://heasarc.nasa.gov/docs/software/fitsio/quick/node9.html
	 */
	int naxis;		// number of dimensions of the image
	long naxes[3];		// size of each dimension
	/* naxes[0] is rx, naxes[1] is ry
	 * Then, for gray images, naxes[2] is unused in FITS but set to 1, and naxis is 2.
	 * For RGB images, naxes[2] is 3 and naxis is 3.
	 * */

	fkeywords keywords; // keywords structure
	gboolean checksum; // flag to save checksum

	char *header;		// entire header of the FITS file. NULL for non-FITS file.
	gchar *unknown_keys; // list of unknown keys

	/* data computed or set by Siril */
	imstats **stats;	// stats of fit for each layer, null if naxes[2] is unknown
	double mini, maxi;	// min and max of the stats->max[3]
	float neg_ratio;	// ratio of pixels with a negative value on total number of pixels

	fitsfile *fptr;		// file descriptor. Only used for file read and write.

	data_type type;		// use of data or fdata is managed by this
	WORD *data;		// 16-bit image data (depending on image type)
	WORD *pdata[3];		// pointers on data, per layer data access (RGB)
	float *fdata;		// same with float
	float *fpdata[3];	// same with float

	gboolean top_down;	// image data is stored top-down, normally false for FITS, true for SER
	gboolean debayer_checked; // whether bayer pattern has already been checked and adjusted or not. This is set true for SER upon opening, not for other formats
	unsigned int orig_ry; // original ry of the image (only set when reading partial)
	int x_offset, y_offset; // x and y offset of partial read wrt to original image
	gboolean focalkey, pixelkey; // flag to define if pixel and focal lengths were read from prefs or from the header keys

	GSList *history;	// Former HISTORY comments of FITS file

	/* ICC Color Management data */
	gboolean color_managed; // Whether color management applies to this FITS
	cmsHPROFILE icc_profile; // ICC color management profile
};

typedef enum {
	SPCC_RED = 1 << RLAYER,
	SPCC_GREEN = 1 << GLAYER,
	SPCC_BLUE = 1 << BLAYER,
	SPCC_CLEAR = SPCC_RED | SPCC_GREEN | SPCC_BLUE,
	SPCC_INVIS = 1 << 7
} spcc_channel;

/* Filter spectral responses are defined by unevenly spaced frequency samples
 * and accompanying spectral responses corresponding to the sampling points. */
typedef struct _spcc_object {
	gchar *model;
	gchar *name;
	gchar *filepath;
	gchar *comment; // optional comment
	gboolean is_dslr; // optional flag to indicate if a DSLR, defaults to FALSE
	int index; // index in the JSON file
	int type;
	int quality;
	spcc_channel channel;
	gchar *manufacturer;
	gchar *source;
	int version;
	gboolean arrays_loaded;
	double *x;  // Wavelength array
	double *y;  // Quantity array
	int n; // Number of points in x and y
} spcc_object;

typedef struct _osc_sensor {
	spcc_object channel[3];
} osc_sensor;

struct spcc_data_store {
	GList *mono_sensors;
	GList *osc_sensors;
	GList *mono_filters[4]; // R, G, B, Lum
	GList *osc_lpf;
	GList *osc_filters;
	GList *wb_ref;
};

/* xpsampdata provides a fixed size struct matched to hold 2nm-spaced data between
 * 336-1020nm for use with Gaia DR3 xp_sampled data products. */
#define XPSAMPLED_LEN 343
typedef struct _xpsampdata {
	const double *x;
	double y[XPSAMPLED_LEN];
} xpsampled;

typedef enum {
	CUT_MONO,
	CUT_COLOR
} cut_mode;

typedef struct {
	guint major_version;
	guint minor_version;
	guint micro_version;
	guint patched_version;
	gboolean beta_version;
	gboolean rc_version;
} version_number;

typedef struct cut_struct {
	point cut_start;			// point marking start of cut line
	point cut_end;			// point dragged while selecting the cut line
	point cut_wn1;				// point for wavenumber 1 for spectroscopic cut
	point cut_wn2;				// point for wavenumber 2 for spectroscopic cut
	gboolean cut_measure;		// Whether or not to measure cuts
	double wavenumber1;
	double wavenumber2;
	gboolean plot_as_wavenumber; // If true, plot as wavenumber; if false, plot as wavelength
	gboolean tri;
	gboolean cfa;
	cut_mode mode;
	int width;
	int step;
	int bg_poly_order;	// polynomial order for background removal
	gboolean plot_spectro_bg;
	gboolean display_graph;
	gboolean save_png_too; // Only takes effect if display_graph == TRUE - ignored otherwise
	char* filename;
	gchar* title; // this is the working copy of the title
	gchar* user_title; // this is the title set by the user, may include brackets
	gboolean title_has_sequence_numbers;
	fits* fit;
	sequence* seq;
	int imgnumber;
	gboolean save_dat;
	gboolean pref_as;
	int vport;
} cut_struct;

// Used in a GSList so we can choose what to stop if multiple children are running
typedef struct _child_info {
	GPid childpid; // Platform-agnostic way to store the child PID
	external_program program; // type of program, eg EXT_STARNET, EXT_ASNET etc
	gchar *name; // argv[0], i.e. the executable name
	GDateTime *datetime; // start time of program - distinguishes between multiple children with the same argv0
} child_info;

struct historic_struct {
	char *filename;
	char history[FLEN_VALUE];
	int rx, ry, nchans;
	data_type type;
	wcs_info wcsdata;
	struct wcsprm *wcslib;
	double focal_length;
	cmsHPROFILE icc_profile;
	gboolean spcc_applied;
};

/* the structure storing information for each layer to be composed
 * (one layer = one source image) and one associated colour */
typedef struct {
	/* widgets data */
	GtkButton *remove_button;
	GtkDrawingArea *color_w;		// the simulated color chooser
	GtkButton *chooser_button;		// the file choosers
	gchar *selected_filename;		// selected filename
	GtkLabel *label;				// the labels
	GtkSpinButton *spinbutton_x;	// the X spin button
	GtkSpinButton *spinbutton_y;	// the Y spin button
	GtkSpinButton *spinbutton_r;	// the rotation spin button
	GtkToggleButton *centerbutton;	// the button to set the center
	double spinbutton_x_value;
	double spinbutton_y_value;
	double spinbutton_r_value;
	/* useful data */
	GdkRGBA color;					// real color of the layer in the image colorspace
	GdkRGBA saturated_color;		// saturated color of the layer in the image colorspace
	GdkRGBA display_color;			// color of the layer in the display colorspace
	fits the_fit;					// the fits for layers
	point center;
} layer;

/* The rendering of the main image is cached. As it can be much larger than the
 * widget in which it's displayed, it can take a lot of time to transform it
 * for rendering. Unfortunately, rendering is requested on each update of a
 * label, so every time the pointer moves above the image.
 * With the caching, the rendering is done once in a cache surface
 * (disp_surface) and this surface is just displayed on subsequent calls if the
 * image, the transform or the widget size has not changed.
 */
struct image_view {
	GtkWidget *drawarea;

	guchar *buf;	// display buffer (image mapped to 3 times 8-bit)
	int full_surface_stride;
	int full_surface_height;
	int view_width;	// drawing area size
	int view_height;

	cairo_surface_t *full_surface;
	cairo_surface_t *disp_surface;	// the cache
};

typedef struct roi {
	fits fit;
	rectangle selection;
	gboolean active;
	gboolean operation_supports_roi;
} roi_t;

typedef struct draw_data {
	cairo_t *cr;	// the context to draw to
	int vport;	// the viewport index to draw
	double zoom;	// the current zoom value
	gboolean neg_view;	// negative view
	cairo_filter_t filter;	// the type of image filtering to use
	guint image_width, image_height;	// image size
	guint window_width, window_height;	// drawing area size
} draw_data_t;

struct gui_icc {
	cmsHPROFILE monitor;
	cmsHPROFILE soft_proof;
	cmsHTRANSFORM proofing_transform;
	cmsUInt32Number proofing_flags;
	gboolean same_primaries;
	gboolean profile_changed;
	guint sh_r;
	guint sh_g;
	guint sh_b;
	guint sh_rgb;
	gboolean iso12646;
};

/* The global data structure of siril gui */
struct guiinf {
	GtkBuilder *builder;		// the builder of the glade interface

	/*********** rendering of gfit, the currently loaded image ***********/
	struct image_view view[MAXVPORT];
	int cvport;			// current viewport, index in the list vport above

	cairo_matrix_t display_matrix;	// matrix used for image rendering (convert image to display coordinates)
	cairo_matrix_t image_matrix;	// inverse of display_matrix (convert display to image coordinates)
	double zoom_value;		// 1.0 is normal zoom, use get_zoom_val() to access it
	point display_offset;		// image display offset

	gboolean translating;		// the image is being translated

	gboolean show_excluded;		// show excluded images in sequences

	int selected_star;		// current selected star in the GtkListStore

	gboolean show_wcs_grid;		// show the equatorial grid over the image
	gboolean show_wcs_disto;		// show the distortions (if present) include in the wcs solution

	psf_star *qphot;		// quick photometry result, highlight a star

	point measure_start;	// quick alt-drag measurement
	point measure_end;

	GSList *user_polygons;	// user defined polygons for the overlay

	void (*draw_extra)(draw_data_t *dd);

	/* List of all scripts from the repository */
	GSList* repo_scripts; // the list of selected scripts is in com.pref
	/* gboolean to confirm the script repository has been opened without error */
	gboolean script_repo_available;
	gboolean spcc_repo_available;

	/*********** Color mapping **********/
	WORD lo, hi;			// the values of the cutoff sliders
	gboolean cut_over;		// display values over hi as negative
	sliders_mode sliders;		// lo/hi, minmax, user
	display_mode rendering_mode;	// pixel value scaling, defaults to LINEAR_DISPLAY or default_rendering_mode if set in preferences
	gboolean unlink_channels;	// only for autostretch
	BYTE remap_index[3][USHRT_MAX + 1];	// abstracted here so it can be used for previews and is easier to change the bit depth
	BYTE *hd_remap_index[3];	// HD remap indexes for the high precision LUTs.
	guint hd_remap_max;		// the maximum index value to use for the HD LUT. Default is 2^22
	gboolean use_hd_remap;		// use high definition LUT for auto-stretch
	struct gui_icc icc;

	/* selection rectangle for registration, FWHM, PSF, coords in com.selection */
	gboolean drawing;		// true if the rectangle is being set (clicked motion)
	cut_struct cut;				// Struct to hold data relating to intensity
								// profile cuts
	pointi start;			// where the mouse was originally clicked to
	pointi origin;			// where the selection was originally located
	gboolean freezeX, freezeY;	// locked axis during modification of a selection
	double ratio;			// enforced ratio of the selection (default is 0: none)
	double rotation;		// selection rotation for dynamic crop

	/* alignment preview data */
	cairo_surface_t *preview_surface[PREVIEW_NB];
	GtkWidget *preview_area[PREVIEW_NB];
	guchar *refimage_regbuffer;	// the graybuf[registration_layer] of the reference image
	cairo_surface_t *refimage_surface;

	int file_ext_filter;		// file extension filter for open/save dialogs

	/* history of the command line. This is a circular buffer (cmd_history)
	 * of size cmd_hist_size, position to be written is cmd_hist_current and
	 * position being browser for display of the history is cmd_hist_display.
	 */
	char **cmd_history;		// the history of the command line
	int cmd_hist_size;		// allocated size
	int cmd_hist_current;		// current command index
	int cmd_hist_display;		// displayed command index
	layer* comp_layer_centering;	// pointer to the layer to center in RGB compositing tool

	roi_t roi; // Region of interest struct
	GSList *mouse_actions;
	GSList *scroll_actions;
	gboolean drawing_polygon; // whether drawing a polygon or not
	GSList *drawing_polypoints; // list of points drawn in MOUSE_ACTION_DRAW_POLY mode
	GdkRGBA poly_ink; // Color and alpha for drawing polygons
	gboolean poly_fill; // Whether to fill drawn polygons
};

struct common_icc {
	cmsHPROFILE mono_linear;
	cmsHPROFILE working_standard;
	cmsHPROFILE mono_standard;
	cmsHPROFILE srgb_profile;
	cmsHPROFILE srgb_out;
	cmsHPROFILE working_out;
	cmsHPROFILE mono_out;
	cmsContext context_single;
	cmsContext context_threaded;
	/* This sets BPC for rendering. It is in com. because it
	 * is also used as the processing intent for non-rendering
	 * color transforms. */
	cmsUInt32Number rendering_flags;
	gboolean srgb_hint; // Hint that FITS has been stretched, for use with file/open
};

/* The global data structure of siril core */
struct cominf {
	GResource *resource; // resources
	gchar *wd;			// current working directory, where images and sequences are

	preferences pref;		// some variables are stored in settings
	gchar *initfile;		// the path of the init file (stores settings)

	sequence seq;			// currently loaded sequence
	single *uniq;			// currently loaded image, if outside sequence

	gsl_histogram *layers_hist[MAXVPORT]; // current image's histograms
	gsl_histogram *sat_hist;
					      // TODO: move in ffit?

	gboolean headless;		// pure console, no GUI
	gboolean script;		// script being executed, always TRUE when headless is
	gboolean python_script;	// python script being executed
	gboolean python_command;	// python is running a Siril command
	GThread *python_init_thread; // python initialization thread, used to monitor startup completion
	GThread *thread;		// the thread for processing
	GMutex mutex;			// a mutex we use for this thread
	GMutex env_mutex;		// a mutex used for updating environment vars (g_setenv is not threadsafe)
	GThread *python_thread;	// the thread for the python interpreter
	char python_magic[9];	// magic number for the python interpreter, used to check .pyc compatibility
	gboolean run_thread;		// the main thread loop condition
	gboolean python_claims_thread;	// prevent other things acquiring the processing thread while a python script has it
	gboolean stop_script;		// abort script execution, not just a command
	GThread *script_thread;		// reads a script and executes its commands
	gboolean script_thread_exited;	// boolean set by the script thread when it exits
	GThread *update_scripts_thread;	// thread used to update the scripts repository, so the GUI can wait it

	int max_images;			// max number of image threads used for parallel execution
	int max_thread;			// max total number of threads used for parallel execution
	int fftw_max_thread;	// max number of threads for FFTW execution

	rectangle selection;		// coordinates of the selection rectangle

	psf_star **stars;		// list of stars detected in the current image
	gboolean star_is_seqdata;	// the only star in stars belongs to seq, don't free it
	double magOffset;		// offset to reduce the real magnitude, single image

	/* history of operations, for the FITS header and the undo feature */
	historic *history;		// the history of all operations on the current image
	int hist_size;			// allocated size
	int hist_current;		// current index
	int hist_display;		// displayed index

	/* all fields below are used by some specific features as a temporary storage */
	GSList *grad_samples;		// list of samples for the background extraction

	GSList *found_object;		// list of objects found in the image from catalogues

	struct spcc_data_store spcc_data; // library of SPCC filters, sensors

	sensor_tilt *tilt;		// computed tilt information

	float* kernel;			// float* to hold kernel for new deconvolution process
	unsigned kernelsize;		// Holds size of kernel (kernel is square kernelsize * kernelsize)
	unsigned kernelchannels;	// Holds number of channels for the kernel
	struct common_icc icc;		// Holds common ICC color profile data
	version_number python_version; // Holds the python version number
	GSList *children;		// List of children; children->data is of type child_info
	gchar *spcc_remote_catalogue;	// Which catalogue to use for SPCC
};

#ifndef MAIN
extern guiinfo gui;
extern cominfo com;		// the main data struct
extern fits gfit;		// currently loaded image
#endif

#endif /*SIRIL */
