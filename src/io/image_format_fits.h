#ifndef _IMAGE_FORMAT_FITS_H
#define _IMAGE_FORMAT_FITS_H

#include <fitsio.h>
#include "core/siril.h"

#define FITS_DOUBLE_BLOC_SIZE 2 * IOBUFLEN // 2 * 2880, the size of a double FITS block, used to allocate bigger chunk and avoid reallocating

/****************** image_format_fits.h ******************/
void read_fits_header(fits *fit);
char *copy_header(fits *fit);

typedef struct {
	char *key;
	char *value;
} header_record;
GSList *read_header_keyvals_strings(fitsfile *fptr);
cmsHPROFILE read_icc_profile_from_fptr(fitsfile *fptr);
int read_icc_profile_from_fits(fits *fit);
int write_icc_profile_to_fits(fits *fit);
int write_icc_profile_to_fptr(fitsfile *fptr, cmsHPROFILE icc_profile);
data_type get_data_type(int bitpix);
void fit_get_photometry_data(fits *fit);
int fit_stats(fitsfile *fptr, float *mini, float *maxi);
int readfits(const char *filename, fits *fit, char *realname, gboolean force_float);
void get_date_data_from_fitsfile(fitsfile *fptr, GDateTime **dt, double *exposure, double *livetime, unsigned int *stack_count);
int import_metadata_from_fitsfile(fitsfile *fptr, fits *to);
void clearfits(fits*);
void clearfits_header(fits*);
int readfits_partial(const char *filename, int layer, fits *fit,
		const rectangle *area, gboolean read_date);
int readfits_partial_all_layers(const char *filename, fits *fit, const rectangle *area);
int read_fits_metadata(fits *fit);
int read_fits_metadata_from_path(const char *filename, fits *fit);
int read_fits_metadata_from_path_first_HDU(const char *filename, fits *fit);
void flip_buffer(int bitpix, void *buffer, const rectangle *area);
int read_opened_fits_partial(sequence *seq, int layer, int index, void *buffer,
		const rectangle *area);
int siril_fits_compress(fits *f);
int save_opened_fits(fits *f);
gchar *set_right_extension(const char *name);
int savefits(const char *name, fits *f);
int copyfits(fits *from, fits *to, unsigned char oper, int layer);
void copy_fits_metadata(fits *from, fits *to);
int copy_fits_from_file(const char *source, const char *destination);
int save1fits16(const char *filename, fits *fit, int layer);
int save1fits32(const char *filename, fits *fit, int layer);
int siril_fits_open_diskfile_img(fitsfile **fptr, const char *filename, int iomode, int *status);
GDateTime *get_date_from_fits(const gchar *filename);

void rgb24bit_to_fits48bit(unsigned char *rgbbuf, fits *fit, gboolean inverted);
void rgb8bit_to_fits16bit(const unsigned char *graybuf, fits *fit);
void rgb48bit_to_fits48bit(const WORD *rgbbuf, fits *fit, gboolean inverted,
		gboolean change_endian);

void fits_flip_top_to_bottom(fits *fit);
void extract_region_from_fits(fits *from, int layer, fits *to,
		const rectangle *area);
int new_fit_image(fits **fit, int width, int height, int nblayer, data_type type);
int new_fit_image_with_data(fits **fit, int width, int height, int nblayer, data_type type, void *data);
void fit_replace_buffer(fits *fit, void *newbuf, data_type newtype);
int fits_change_depth(fits *fit, int layers);
int extract_fits(fits *from, fits *to, int channel, gboolean to_float);
void keep_only_first_channel(fits *fit);
void fit_debayer_buffer(fits *fit, void *newbuf);

GdkPixbuf* get_thumbnail_from_fits(char *filename, gchar **descr);

// internal read of FITS file, for FITS images and FITS sequences
void manage_bitpix(fitsfile *fptr, int *bitpix, int *orig_bitpix);
int read_fits_with_convert(fits* fit, const char* filename, gboolean force_float);
int internal_read_partial_fits(fitsfile *fptr, unsigned int ry,
		int bitpix, void *dest, int layer, const rectangle *area);
int siril_fits_create_diskfile(fitsfile **fptr, const char *filename, int *status);
void update_fits_header(fits *fit);
void save_fits_header(fits *fit);
void report_fits_error(int status);

int check_fits_params(fitsfile *fptr, int *oldbitpix, int *oldnaxis, long *oldnaxes, gboolean relax_dimcheck);
int check_loaded_fits_params(fits *ref, ...);

void merge_fits_headers_to_result2(fits *result, fits **f, gboolean do_sum);
void merge_fits_headers_to_result(fits *result, gboolean do_sum, fits *f1,...);
int get_xpsampled(double* xps, const gchar *filename, int i);
void process_keyword_string_value(const char *input, char *output, gboolean condition);
int updateFITSKeyword(fits *fit, const gchar *key, const gchar *newkey, const gchar *value, const gchar *comment, gboolean verbose, gboolean isfitseq);
int associate_header_to_memfile(const char *header, fitsfile *fptr);
int fits_parse_header_str(fits *fit, const char *header);
int fits_swap_image_data(fits *a, fits *b);

int save_wcs_fits(fits *f, const gchar *filename);
int save_mask_fits(int rx, int ry, float *buffer, const gchar *name);
int read_mask_fits_area(const gchar *name, rectangle *area, int ry, float *mask);
int read_drizz_fits_area(const gchar *name, int layer, rectangle *area, int ry, float *drizz);

void interpolate_nongreen(fits *fit);
const char* get_cfa_from_pattern(sensor_pattern pattern);

#endif
