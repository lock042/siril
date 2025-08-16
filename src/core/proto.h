#ifndef PROTO_H_
#define PROTO_H_
#include <stdint.h>
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
int round_to_int(double x);
int roundf_to_int(float x);
WORD round_to_WORD(double x);
BYTE round_to_BYTE(double x);
BYTE roundf_to_BYTE(float f);
WORD roundf_to_WORD(float f);
signed short roundf_to_short(float f);
guint float_to_max_range(float f, guint max);
int round_to_ceiling_multiple(int x, int factor);
BYTE conv_to_BYTE(double x);
int truncate_to_int32(uint64_t x);
WORD truncate_to_WORD(int x);
BYTE truncate_to_BYTE(WORD x);
int set_int_in_interval(int val, int low, int high);
float set_float_in_interval(float val, float low, float high);
double set_double_in_interval(double val, double low, double high);
float ushort_to_float_range(WORD w);
float uchar_to_float_range(BYTE w);
float double_ushort_to_float_range(double d);
WORD float_to_ushort_range(float f);
signed short float_to_short_range(float f);
BYTE float_to_uchar_range(float f);
float ushort_to_float_bitpix(const fits *fit, const WORD value);
WORD *float_buffer_to_ushort(const float *buffer, size_t ndata);
signed short *float_buffer_to_short(const float *buffer, size_t ndata);
signed short *ushort_buffer_to_short(const WORD *buffer, size_t ndata);
float *uchar_buffer_to_float(BYTE *buffer, size_t ndata);
float *ushort_buffer_to_float(WORD *buffer, size_t ndata);
float *ushort8_buffer_to_float(WORD *buffer, size_t ndata);
gboolean test_double_eq(double a, double b, double epsilon);
uint16_t change_endianness16(uint16_t x);
uint16_t cpu_to_le16(uint16_t x);
uint16_t cpu_to_be16(uint16_t x);
uint16_t le16_to_cpu(uint16_t x);
uint16_t be16_to_cpu(uint16_t x);
uint32_t change_endianness32(uint32_t x);
uint32_t cpu_to_le32(uint32_t x);
uint32_t cpu_to_be32(uint32_t x);
uint32_t le32_to_cpu(uint32_t x);
uint32_t be32_to_cpu(uint32_t x);
uint32_t be24_to_cpu(const BYTE x[3]);
uint64_t change_endianness64(uint64_t x);
uint64_t cpu_to_le64(uint64_t x);
uint64_t cpu_to_be64(uint64_t x);
uint64_t le64_to_cpu(uint64_t x);
uint64_t be64_to_cpu(uint64_t x);
gboolean isrgb(const fits *fit);
const char *channel_number_to_name(int channel);
int get_extension_index(const char *filename);
image_type get_type_from_filename(const gchar *filename);
char* remove_ext_from_filename(const char *basename);
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

#ifdef __cplusplus
}
#endif

#endif
