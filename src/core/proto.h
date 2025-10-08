#ifndef PROTO_H_
#define PROTO_H_
#include <stdint.h>
#include <math.h> // for fabs() in test_double_eq
#include <sys/time.h>
#include "core/siril.h"
#ifdef HAVE_LIBTIFF
#define uint64 uint64_hack_
#define int64 int64_hack_
#include <tiffio.h>
#undef uint64
#undef int64
#endif

#ifdef __cplusplus
extern "C" {
#endif

/****************** image_formats_internal.h ******************/
/* BMP */
int readbmp(const char*, fits*);
int savebmp(const char*, fits*);

/* PNM */
int import_pnm_to_fits(const char *filename, fits *fit);
int saveNetPBM(const char *name, fits *fit);

/* PIC */
struct pic_struct {
	unsigned long magic;
	unsigned short width;
	unsigned short height;
	unsigned short bin[6];
	unsigned short nbplane;
	unsigned short hi;
	unsigned short lo;
	char *date;
	char *time;

	// internal stuff
	FILE *file;
};
int readpic(const char *name, fits *fit);

/****************** image_formats_libraries.h ******************/
#ifdef HAVE_LIBTIFF
gboolean get_tiff_compression();
int readtif(const char *name, fits *fit, gboolean force_float, gboolean verbose);
void get_tif_data_from_ui(fits *fit, gchar **description, gchar **copyright);
int savetif(const char *name, fits *fit, uint16_t bitspersample,
		const char *description, const char *copyright,
		gboolean tiff_compression, gboolean embeded_icc, gboolean verbose);
#endif

#ifdef HAVE_LIBXISF
int readxisf(const char *name, fits *fit, gboolean force_float);
#endif

#ifdef HAVE_LIBJPEG
int readjpg(const char*, fits*);
int savejpg(const char*, fits*, int quality, gboolean verbose);
#endif

#ifdef HAVE_LIBPNG
int readpng(const char*, fits*);
int savepng(const char *filename, fits *fit, uint32_t bytes_per_sample,
		gboolean is_colour);
#endif

#ifdef HAVE_LIBRAW
int open_raw_files(const char*, fits*);
#endif

#ifdef HAVE_LIBHEIF
int readheif(const char* name, fits *fit, gboolean interactive);
int saveheifavif(const char* name, fits *fit, int quality, gboolean lossless, gboolean is_av1f, int max_bitdepth);
#endif

#ifdef HAVE_LIBJXL
int readjxl(const char* name, fits *fit);
int savejxl(const char* name, fits* fit, int effort, double quality, gboolean force_8bit);
#endif
/****************** utils.h ******************/

WORD *float_buffer_to_ushort(const float *buffer, size_t ndata);
signed short *float_buffer_to_short(const float *buffer, size_t ndata);
signed short *ushort_buffer_to_short(const WORD *buffer, size_t ndata);
float *uchar_buffer_to_float(BYTE *buffer, size_t ndata);
float *ushort_buffer_to_float(WORD *buffer, size_t ndata);
float *ushort8_buffer_to_float(WORD *buffer, size_t ndata);
const char *channel_number_to_name(int channel);
int get_extension_index(const char *filename);
image_type get_type_from_filename(const gchar *filename);
char* remove_ext_from_filename(const char *basename);
char *remove_all_ext_from_filename(const char *filename);
gchar *replace_ext(const char *path, const char *new_ext);
gboolean string_is_a_path(const char *file);
int is_readable_file(const char *filename);
int is_symlink_file(const char *filename);
gint siril_mkdir_with_parents(const gchar* pathname, gint mode);
gboolean is_forbiden_in_filename(gchar c);
gboolean file_name_has_invalid_chars(const char *name);
void replace_invalid_chars(char *name, char repl);
gchar* replace_wide_char(const gchar *str);
int stat_file(const char *filename2, image_type *type, char **realname);
const char* get_filename_ext(const char *filename);

int siril_change_dir(const char *dir, gchar **err);
int update_sequences_list(const char *sequence_name_to_select);
void expand_home_in_filename(char *filename, int size);
double get_normalized_value(fits*);
void swap_param(double*, double*);
gchar* str_append(char **data, const char *newdata);
char *format_basename(char *root, gboolean can_free);
float compute_slope(WORD *lo, WORD *hi);
gchar *siril_get_file_info(const gchar *filename, GdkPixbuf *pixbuf);
gchar *siril_truncate_str(gchar *str, gint size);
gchar **glist_to_array(GList *list, int *arg_count);
gchar* url_cleanup(const gchar *uri_string);
void remove_spaces_from_str(gchar *s);
gboolean string_has_space(const gchar *str);
void remove_trailing_eol(char *str);
gboolean string_is_a_number(const char *str);
#if !GLIB_CHECK_VERSION(2,68,0)
guint g_string_replace(GString *string, const gchar *find, const gchar *replace,
		guint limit);
#endif
char *str_replace(char *orig, const char *rep, char *with);
void replace_spaces_from_str(gchar *s, char c);
void replace_char_from_str(gchar *s, gchar in, gchar out);
gchar *build_string_from_words(char **words);
void append_elements_to_array(char **array, char **elements);
const gchar *get_com_ext(gboolean fz);
gchar *siril_any_to_utf8 (const gchar  *str, gssize len, const gchar *warning_format, ...);

int siril_to_display(double fx, double fy, double *dx, double *dy, int ry);
int display_to_siril(double dx, double dy, double *fx, double *fy, int ry);
int fits_to_display(double fx, double fy, double *dx, double *dy, int ry);
gchar *siril_file_chooser_get_filename(GtkFileChooser *chooser);
GSList *siril_file_chooser_get_filenames(GtkFileChooser *chooser);
int interleave(fits *fit, int max_bitdepth, void **interleaved_buffer, int *bit_depth, gboolean force_even);
int count_lines_in_textfile(const gchar *filename);
void copy_filename(const char *filename, char *truncated_filename, size_t max_length);
gboolean is_string_numeric(const gchar *str);
const gchar* find_first_numeric(const gchar *string);
const gchar* find_first_nonnumeric(const gchar *string);
int count_pattern_occurence(const gchar *string, const gchar *pattern);
guint gui_function(GSourceFunc idle_function, gpointer data);
gchar *find_file_in_directory(gchar *basename, const gchar *path);
gchar *find_file_recursively(gchar *basename, const gchar *top_path);
char *strdupnullok(char *data);
gchar* remove_extension_from_path(const gchar* filepath);
gboolean delete_directory(const gchar *dir_path, GError **error);
gchar *posix_path_separators(const gchar *path);

/****************** quantize.h ***************/
int siril_fits_img_stats_ushort(WORD *array, long nx, long ny,
		long *ngoodpix, WORD *minvalue, WORD *maxvalue,
		double *mean, double *sigma, double *noise1, double *noise2,
		double *noise3, double *noise5, threading_type threads, int *status);

int siril_fits_img_stats_float(float *array, long nx, long ny,
		long *ngoodpix, float *minvalue, float *maxvalue,
		double *mean, double *sigma, double *noise1, double *noise2,
		double *noise3, double *noise5, threading_type threads, int *status);

/****************** siril.h ******************/

int threshlo(fits *fit, WORD level);
int threshhi(fits *fit, WORD level);
int nozero(fits *fit, WORD level);
int gaussian_blur_RT(fits *fit, double sigma, int threads);
int unsharp(fits*, double sigma, double mult, gboolean verbose);
float entropy(fits *fit, int layer, rectangle *area, const imstats *opt_stats);
int loglut(fits *fit);
int ddp(fits *a, float lev, float coef, float sig);
int visu(fits *fit, int low, int high);
int fill(fits *fit, int level, const rectangle *arearg);
int off(fits *a, float level);
double background(fits* fit, int reqlayer, rectangle *selection, threading_type threads);
void compute_grey_flat(fits *fit);

/****************** seqfile.h ******************/
sequence* readseqfile(const char *name);
int writeseqfile(sequence *seq);
gboolean existseq(const char *name);
int buildseqfile(sequence *seq, int force_recompute);

/****************** statistics_list.h ******************/
void computeStat();

/****************** chelperfuncs.h *********************/
float bilinear(float *x, int w, int h, float i, float j);
float bilinear_ushort(WORD *x, int w, int h, float i, float j);

/******** static inlines to aid vectorization **********/

/**
 * Round double value to an integer (branchless)
 * @param x value to round
 * @return an integer
 */
static inline int round_to_int(double x) {
	x = (x > (double)INT_MAX - 0.5) ? (double)INT_MAX - 0.5 : x;
	x = (x < (double)INT_MIN + 0.5) ? (double)INT_MIN + 0.5 : x;
	double offset = (x >= 0.0) ? 0.5 : -0.5;
	return (int)(x + offset);
}

/**
 * Round float value to an integer (branchless)
 * @param x value to round
 * @return an integer
 */
static inline int roundf_to_int(float x) {
	x = (x > (float)INT_MAX - 0.5f) ? (float)INT_MAX - 0.5f : x;
	x = (x < (float)INT_MIN + 0.5f) ? (float)INT_MIN + 0.5f : x;
	float offset = (x >= 0.0f) ? 0.5f : -0.5f;
	return (int)(x + offset);
}

/**
 * Round double value to a WORD (branchless)
 * @param x value to round
 * @return a WORD
 */
static inline WORD round_to_WORD(double x) {
	x = x + 0.5;
	x = (x > USHRT_MAX_DOUBLE) ? USHRT_MAX_DOUBLE : x;
	x = (x < 0.0) ? 0.0 : x;
	return (WORD)x;
}

/**
 * Round double value to a BYTE (branchless)
 * @param x value to round
 * @return a BYTE
 */
static inline BYTE round_to_BYTE(double x) {
	x = x + 0.5;
	x = (x > UCHAR_MAX_DOUBLE) ? UCHAR_MAX_DOUBLE : x;
	x = (x < 0.0) ? 0.0 : x;
	return (BYTE)x;
}

/**
 * Round float value to a BYTE (branchless)
 * @param f value to round
 * @return a truncated and rounded BYTE
 */
static inline BYTE roundf_to_BYTE(float f) {
	f = f + 0.5f;
	f = (f > (float)UCHAR_MAX) ? (float)UCHAR_MAX : f;
	f = (f < 0.0f) ? 0.0f : f;
	return (BYTE)f;
}

/**
 * Round float value to a short (branchless)
 * @param f value to round
 * @return a truncated and rounded short
 */
static inline signed short roundf_to_short(float f) {
	f = f + 0.5f;
	f = (f > (float)SHRT_MAX) ? (float)SHRT_MAX : f;
	f = (f < (float)SHRT_MIN) ? (float)SHRT_MIN : f;
	return (signed short)f;
}

/**
 * Scale float value to a maximum value up to 2^32-1
 * and return as guint32 (branchless)
 * @param f value to scale
 * @param max float range [0f..1f] scales to guint32 range [0..max]
 * @return a guint32
 */
static inline guint float_to_max_range(float f, guint max) {
	f = f * max + 0.5f;
	f = (f > (float)max) ? (float)max : f;
	f = (f < 0.0f) ? 0.0f : f;
	return (guint)f;
}

/**
 * Compute a ceiling factor (branchless)
 * @param x the number to test
 * @param factor the factor
 * @return x if it is a factor of factor or the next factor
 */
static inline int round_to_ceiling_multiple(int x, int factor) {
	int remainder = x % factor;
	int needs_rounding = (remainder != 0);
	return x + (factor - remainder) * needs_rounding;
}

/**
 * convert double value to a BYTE (branchless)
 * @param x value to convert
 * @return a BYTE
 */
static inline BYTE conv_to_BYTE(double x) {
	x = (x / USHRT_MAX_DOUBLE) * UCHAR_MAX_DOUBLE;
	return (BYTE)x;
}

/**
 * truncate a 64 bit unsigned int to a 32 bit signed int (branchless)
 * @param x value to truncate
 * @return an int
 */
static inline int truncate_to_int32(uint64_t x) {
	uint64_t clamped = (x < (uint64_t)INT_MAX) ? x : (uint64_t)INT_MAX;
	return (int)clamped;
}

/**
 * Truncate int to WORD (branchless)
 */
static inline WORD truncate_to_WORD(int x) {
	x = (x < 0) ? 0 : x;
	x = (x > USHRT_MAX) ? USHRT_MAX : x;
	return (WORD)x;
}

/**
 * Scale WORD value to standard float range [0,1]
 * and return as float (branchless)
 */
static inline float ushort_to_float_range(WORD x) {
	return (float) x * INV_USHRT_MAX_SINGLE;
}

/**
 * Helper for roundf_to_WORD (used by float_to_ushort_range)
 */
static inline WORD roundf_to_WORD(float f) {
	f = f + 0.5f;
	f = (f > (float)USHRT_MAX) ? (float)USHRT_MAX : f;
	f = (f < 0.0f) ? 0.0f : f;
	return (WORD)f;
}

/**
 * Truncate WORD to BYTE (branchless)
 */
static inline BYTE truncate_to_BYTE(WORD x) {
	x = (x > UCHAR_MAX) ? UCHAR_MAX : x;
	return (BYTE)x;
}

/**
 * Round WORD to BYTE with scaling (branchless)
 */
static inline BYTE roundw_to_BYTE(WORD input) {
	input = (input >= USHRT_MAX - 127) ? (USHRT_MAX - 127) : input;
	return (uint8_t)((input + 128) >> 8);
}

/**
 * Clamp an integer value in the interval given by [low, high] (branchless)
 * @param val value to be checked
 * @param low low value of the interval
 * @param high high value of the interval
 * @return a new value set in the [low, high] interval
 */
static inline int set_int_in_interval(int val, int low, int high) {
	val = (val < low) ? low : val;
	val = (val > high) ? high : val;
	return val;
}

/**
 * Clamp a float value in the interval given by [low, high] (branchless)
 * @param val value to be checked
 * @param low low value of the interval
 * @param high high value of the interval
 * @return a new value set in the [low, high] interval
 */
static inline float set_float_in_interval(float val, float low, float high) {
	val = (val < low) ? low : val;
	val = (val > high) ? high : val;
	return val;
}

/**
 * Clamp a double value in the interval given by [low, high] (branchless)
 * @param val value to be checked
 * @param low low value of the interval
 * @param high high value of the interval
 * @return a new value set in the [low, high] interval
 */
static inline double set_double_in_interval(double val, double low, double high) {
	val = (val < low) ? low : val;
	val = (val > high) ? high : val;
	return val;
}

/**
 * convert an unsigned char value to siril's representation of float values [0, 1]
 * @param w value to convert
 * @return the float equivalent
 */
static inline float uchar_to_float_range(BYTE w) {
	return (float)w * INV_UCHAR_MAX_SINGLE;
}

/**
 * convert an double value from the unsigned short range to siril's representation
 * of float values [0, 1]
 * @param d value to convert
 * @return the float equivalent
 */
static inline float double_ushort_to_float_range(double d) {
	return (float)d * INV_USHRT_MAX_SINGLE;
}

/**
 * convert a siril float [0, 1] to an unsigned short
 * @param f value to convert
 * @return the unsigned short equivalent
 */
static inline WORD float_to_ushort_range(float f) {
	return roundf_to_WORD(f * USHRT_MAX_SINGLE);
}

/**
 * convert a siril float [0, 1] to a signed short
 * @param f value to convert
 * @return the signed short equivalent
 * (-SHRT_MAX - 1)
 */
static inline signed short float_to_short_range(float f) {
	return roundf_to_short((f * USHRT_MAX_SINGLE) - SHRT_MAX_SINGLE - 1.0f);
}

/**
 * convert a siril float [0, 1] to an unsigned char
 * @param f value to convert
 * @return the unsigned char equivalent
 */
static inline BYTE float_to_uchar_range(float f) {
	return roundf_to_BYTE(f * UCHAR_MAX_SINGLE);
}

/**
 * convert the pixel value of an image to a float [0, 1] normalized using bitpix
 * value depending on btpix (branchless)
 * @param fit the image the data is from
 * @return a float [0, 1] value for the given integer value
 */
static inline float ushort_to_float_bitpix(const fits *fit, const WORD value) {
	const float fval = (float)value;
	int is_byte = (fit->orig_bitpix == BYTE_IMG);
	float scale = is_byte ? INV_UCHAR_MAX_SINGLE : INV_USHRT_MAX_SINGLE;
	return fval * scale;
}

/**
 * Test equality between two double number
 * @param a
 * @param b
 * @param epsilon
 * @return
 */
static inline gboolean test_double_eq(double a, double b, double epsilon) {
	return (fabs(a - b) <= epsilon);
}

/**
 * change endianness of a 16 bit unsigned int
 * @param x value to convert
 * @return byte-swapped value
 */
static inline uint16_t change_endianness16(uint16_t x) {
    return (x >> 8) | (x << 8);
}

/**
 * convert a 16 bit unsigned int in CPU byte order to little endian
 * @param x value to convert
 * @return little endian value
 */
static inline uint16_t cpu_to_le16(uint16_t x) {
#ifdef __BIG_ENDIAN__
    return change_endianness16(x);
#else
    return x;
#endif
}

/**
 * convert a 16 bit unsigned int in CPU byte order to big endian
 * @param x value to convert
 * @return big endian value
 */
static inline uint16_t cpu_to_be16(uint16_t x) {
#ifdef __BIG_ENDIAN__
    return x;
#else
    return change_endianness16(x);
#endif
}

/**
 * convert a 16 bit unsigned int from little endian to CPU byte order
 * @param x little endian value to convert
 * @return value
 */
static inline uint16_t le16_to_cpu(uint16_t x) {
    return cpu_to_le16(x);
}

/**
 * convert a 16 bit unsigned int from big endian to CPU byte order
 * @param x big endian value to convert
 * @return value
 */
static inline uint16_t be16_to_cpu(uint16_t x) {
    return cpu_to_be16(x);
}

static inline uint32_t be24_to_cpu(const BYTE x[3]) {
#ifdef __BIG_ENDIAN__
	uint32_t r = ((x[2] << 16) | (x[1] << 8) | x[0]);
#else
	uint32_t r = ((x[0] << 16) | (x[1] << 8) | x[2]);
#endif
	return r;
}

/**
 * change endianness of a 32 bit unsigned int
 * @param x value to convert
 * @return byte-swapped value
 */
static inline uint32_t change_endianness32(uint32_t x) {
    return (x >> 24) | ((x & 0xFF0000) >> 8) | ((x & 0xFF00) << 8) | (x << 24);
}

/**
 * convert a 32 bit unsigned int in CPU byte order to little endian
 * @param x value to convert
 * @return little endian value
 */
static inline uint32_t cpu_to_le32(uint32_t x) {
#ifdef __BIG_ENDIAN__
    return change_endianness32(x);
#else
    return x;
#endif
}

/**
 * convert a 32 bit unsigned int in CPU byte order to big endian
 * @param x value to convert
 * @return big endian value
 */
static inline uint32_t cpu_to_be32(uint32_t x) {
#ifdef __BIG_ENDIAN__
    return x;
#else
    return change_endianness32(x);
#endif
}

/**
 * convert a 32 bit unsigned int from little endian to CPU byte order
 * @param x little endian value to convert
 * @return value
 */
static inline uint32_t le32_to_cpu(uint32_t x) {
    return cpu_to_le32(x);
}

/**
 * convert a 32 bit unsigned int from big endian to CPU byte order
 * @param x big endian value to convert
 * @return value
 */
static inline uint32_t be32_to_cpu(uint32_t x) {
    return cpu_to_be32(x);
}

/**
 * change endianness of a 64 bit unsigned int
 * @param x value to convert
 * @return byte-swapped value
 */
static inline uint64_t change_endianness64(uint64_t x) {
    return
        (x >> 56)
        | ((x & 0xFF000000000000) >> 40)
        | ((x & 0xFF0000000000) >> 24)
        | ((x & 0xFF00000000) >> 8)
        | ((x & 0xFF000000) << 8)
        | ((x & 0xFF0000) << 24)
        | ((x & 0xFF00) << 40)
        | (x << 56);
}

/**
 * convert a 64 bit unsigned int in CPU byte order to little endian
 * @param x value to convert
 * @return little endian value
 */
static inline uint64_t cpu_to_le64(uint64_t x) {
#ifdef __BIG_ENDIAN__
	return change_endianness64(x);
#else
	return x;
#endif
}

/**
 * convert a 64 bit unsigned int in CPU byte order to big endian
 * @param x value to convert
 * @return big endian value
 */
static inline uint64_t cpu_to_be64(uint64_t x) {
#ifdef __BIG_ENDIAN__
	return x;
#else
	return change_endianness64(x);
#endif
}

/**
 * convert a 64 bit unsigned int from little endian to CPU byte order
 * @param x little endian value to convert
 * @return value
 */
static inline uint64_t le64_to_cpu(uint64_t x) {
	return cpu_to_le64(x);
}

/**
 * convert a 64 bit unsigned int from big endian to CPU byte order
 * @param x big endian value to convert
 * @return value
 */
static inline uint64_t be64_to_cpu(uint64_t x) {
	return cpu_to_be64(x);
}

/**
 * Test if fit has 3 channels
 * @param fit input FITS image
 * @return TRUE if fit image has 3 channels
 */
static inline gboolean isrgb(const fits *fit) {
	return (fit->naxis == 3);
}

#ifdef __cplusplus
}
#endif

#endif
