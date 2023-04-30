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
#include "core/atomic.h"

#define _(String) gettext (String)
#define gettext_noop(String) String
#define N_(String) gettext_noop (String)

#ifdef SIRIL_OUTPUT_DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#if defined (HAVE_FFTW3F_OMP) || defined (HAVE_FFTW3F_THREADS)
#define HAVE_FFTW3F_MULTITHREAD
#endif

#define GLADE_FILE "siril3.glade"

/* https://stackoverflow.com/questions/1644868/define-macro-for-debug-printing-in-c */
#define siril_debug_print(fmt, ...) \
	do { if (DEBUG_TEST) fprintf(stdout, fmt, ##__VA_ARGS__); } while (0)

#define PRINT_ALLOC_ERR fprintf(stderr, "Out of memory in %s (%s:%d) - aborting\n", __func__, __FILE__, __LINE__)
#define PRINT_ANOTHER_THREAD_RUNNING siril_log_message(_("Another task is already in progress, ignoring new request.\n"))

#ifndef RT_INCLUDE
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
#define SIGMA_PER_FWHM 2.35482
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
	TYPEAVI = (1 << 20),
	TYPESER = (1 << 21),
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

#define MAX_STARS 200000		// maximum length of com.stars
#define MAX_STARS_FITTED 2000	// maximum number of stars fitted for registration
#define MIN_STARS_FITTED 100	// minimum number of stars fitted for registration
#define DEF_BOX_RADIUS 5 // default radius of the box in starfinder_conf

#define INDEX_MAX 65535		// maximum index for images

typedef struct sequ sequence;
typedef struct ffit fits;
typedef struct guiinf guiinfo;
typedef struct cominf cominfo;
typedef struct historic_struct historic;
typedef struct fwhm_struct psf_star;
typedef struct tilt_struct sensor_tilt;

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
	SEQ_REGULAR, SEQ_SER, SEQ_FITSEQ,
#ifdef HAVE_FFMS2
	SEQ_AVI,
#endif
	SEQ_INTERNAL
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
	EXT_GNUPLOT // not used for now
} external_program;

/* image data, exists once for each image */
typedef struct {
	int filenum;		/* real file index in the sequence, i.e. for mars9.fit = 9 */
	gboolean incl;		/* selected in the sequence, included for future processings? */
	GDateTime *date_obs;	/* date of the observation, processed and copied from the header */
	int rx, ry;
} imgdata;

/* this structure is used to characterize the statistics of the image */
typedef struct {
	long total,	// number of pixels
	     ngoodpix;	// number of non-zero pixels
	double mean, median, sigma, avgDev, mad, sqrtbwmv,
	       location, scale, min, max, normValue, bgnoise;

	atomic_int* _nb_refs;	// reference counting for data management
} imstats;

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
	int bitpix;		// image pixel format, from fits
	int reference_image;	// reference image for registration
	imgdata *imgparam;	// a structure for each image of the sequence
	regdata **regparam;	// *regparam[nb_layers], may be null if nb_layers is unknown
	imstats ***stats;	// statistics of the images for each layer, may be null too
	/* in the case of a CFA sequence, depending on the opening mode, we cannot store
	 * and use everything that was in the seqfile, so we back them up here */
	regdata **regparam_bkp;	// *regparam[3], null if nothing to back up
	imstats ***stats_bkp;	// statistics of the images for 3 layers, may be null too

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

	double upscale_at_stacking;// up-scale factor during stacking (see #215)

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
	double equinox;
	double crpix[2];
	double crval[2];
	double cdelt[2];
	double pc[2][2];
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

	/* data obtained from the FITS file */
	char *header;		// entire header of the FITS file. NULL for non-FITS file.
	WORD lo;		// MIPS-LO key in FITS file, "Lower visualization cutoff"
	WORD hi;		// MIPS-HI key in FITS file, "Upper visualization cutoff"
	double data_max;	// used to check if 32b float is in the [0, 1] range
	double data_min;	// used to check if 32b float is in the [0, 1] range
	float pixel_size_x, pixel_size_y;	// XPIXSZ and YPIXSZ keys
	unsigned int binning_x, binning_y;	// XBINNING and YBINNING keys
	char row_order[FLEN_VALUE];
	GDateTime *date, *date_obs;
	double expstart, expend;
	char filter[FLEN_VALUE];		// FILTER key
	char image_type[FLEN_VALUE];		// IMAGETYP key
	char object[FLEN_VALUE];		// OBJECT key
	char instrume[FLEN_VALUE];		// INSTRUME key
	char telescop[FLEN_VALUE];		// TELESCOP key
	char observer[FLEN_VALUE];		// OBSERVER key
	double sitelat;				// SITE LATITUDE key
	double sitelong;			// SITE LONGITUDE key
	double siteelev;			// SITE LONGITUDE key
	char bayer_pattern[FLEN_VALUE];		// BAYERPAT key Bayer Pattern if available
	int bayer_xoffset, bayer_yoffset;
	double airmass;                   // relative optical path length through atmosphere.
	/* data obtained from FITS or RAW files */
	double focal_length, iso_speed, exposure, aperture, ccd_temp;
	double livetime;		// total exposure
	guint stackcnt;			// number of stacked frame
	double cvf;			// Conversion factor (e-/adu)
	int key_gain, key_offset;	// Gain, Offset values read in camera headers.

	/* Plate Solving data */
	wcs_info wcsdata;		// data from the header
#ifdef HAVE_WCSLIB
	struct wcsprm *wcslib;		// struct of the lib
#endif

	/* data used in the Fourier space */
	dft_info dft;

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

	GSList *history;	// Former HISTORY comments of FITS file
};

typedef struct {
	double x, y;
} point;

typedef struct {
	float x, y;
} pointf;

typedef struct {
	int x, y;
} pointi;

struct historic_struct {
	char *filename;
	char history[FLEN_VALUE];
	int rx, ry;
	data_type type;
	wcs_info wcsdata;
	double focal_length;
};

typedef struct _GNUPLOT_CTRL_ {
    /** Pipe to gnuplot process */
    FILE* gnucmd ;
	FILE* gnumon ;

    /** Number of currently active plots */
    int nplots ;
	/** Current plot window **/
	gboolean replot; // Add additional plots to the current one in the current window
    /** Current plotting style */
    char pstyle[32] ;

    /** Pointer to table of names of temporary files */
    char** tmp_filename_tbl ;
    /** Number of temporary files */
    int ntmp ;
	GThread* thread;
	int child_fd_stdin;
	int child_fd_stderr;
	GPid child_pid;
	gboolean running;
} gnuplot_ctrl ;

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

typedef struct draw_data {
	cairo_t *cr;	// the context to draw to
	int vport;	// the viewport index to draw
	double zoom;	// the current zoom value
	gboolean neg_view;	// negative view
	cairo_filter_t filter;	// the type of image filtering to use
	guint image_width, image_height;	// image size
	guint window_width, window_height;	// drawing area size
} draw_data_t;

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

	psf_star *qphot;		// quick photometry result, highlight a star

	void (*draw_extra)(draw_data_t *dd);

	/*********** Color mapping **********/
	WORD lo, hi;			// the values of the cutoff sliders
	gboolean cut_over;		// display values over hi as negative
	sliders_mode sliders;		// lo/hi, minmax, user
	display_mode rendering_mode;	// pixel value scaling, defaults to LINEAR_DISPLAY or default_rendering_mode if set in preferences
	gboolean unlink_channels;	// only for autostretch
	BYTE remap_index[3][USHRT_MAX];	// abstracted here so it can be used for previews and is easier to change the bit depth
	BYTE *hd_remap_index[3];	// HD remap indexes for the high precision LUTs.
	guint hd_remap_max;		// the maximum index value to use for the HD LUT. Default is 2^22
	gboolean use_hd_remap;		// use high definition LUT for auto-stretch

	/* selection rectangle for registration, FWHM, PSF, coords in com.selection */
	gboolean drawing;		// true if the rectangle is being set (clicked motion)
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
};

/* The global data structure of siril core */
struct cominf {
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
	GThread *thread;		// the thread for processing
	GMutex mutex;			// a mutex we use for this thread
	gboolean run_thread;		// the main thread loop condition
	gboolean stop_script;		// abort script execution, not just a command
	GThread *script_thread;		// reads a script and executes its commands
	gboolean script_thread_exited;	// boolean set by the script thread when it exits

	int max_images;			// max number of image threads used for parallel execution
	int max_thread;			// max total number of threads used for parallel execution

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

	sensor_tilt *tilt;		// computed tilt information

	external_program child_is_running;	// external_program id to check if there is a child process running

	float* kernel;			// float* to hold kernel for new deconvolution process
	unsigned kernelsize;		// Holds size of kernel (kernel is square kernelsize * kernelsize)
	unsigned kernelchannels;	// Holds number of channels for the kernel
#ifdef _WIN32
	void* childhandle;		// For Windows, handle of a child process
#else
	pid_t childpid;			// For other OSes, PID of a child process
#endif
	gnuplot_ctrl **gnuplot_handles; // list of gnuplot handles
	int num_gnuplot_handles; // how many gnuplot handles are in the list

};

#ifndef MAIN
extern guiinfo gui;
extern cominfo com;		// the main data struct
extern fits gfit;		// currently loaded image
#endif

#endif /*SIRIL */
