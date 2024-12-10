#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include <stdint.h>
#include "../core/siril.h"
#include "../core/processing.h"
#include "../algos/PSF.h"

gboolean populate_seqcombo(gpointer user_data);
int	read_single_sequence(char *realname, image_type imagetype);
char *normalize_seqname(char *name, gboolean add_underscore);
int	check_seq();
int	seq_check_basic_data(sequence *seq, gboolean load_ref_into_gfit);
gboolean set_seq(gpointer user_data);
char *	seq_get_image_filename(sequence *seq, int index, char *name_buf);
int	seq_read_frame(sequence *seq, int index, fits *dest, gboolean force_float, int thread_id);
int seq_read_frame_metadata(sequence *seq, int index, fits *dest);
int	seq_read_frame_part(sequence *seq, int layer, int index, fits *dest, const rectangle *area, gboolean do_photometry, int thread_id);
int	seq_load_image(sequence *seq, int index, gboolean load_it);
gboolean seq_load_image_in_thread(gpointer user_data);
int64_t seq_compute_size(sequence *seq, int nb_frames, data_type type);
gboolean check_if_seq_exist(gchar *name, gboolean name_is_base);
int	seq_open_image(sequence *seq, int index);
void	seq_close_image(sequence *seq, int index);
int	seq_opened_read_region(sequence *seq, int layer, int index, void *buffer, const rectangle *area, int thread_id);
void	set_fwhm_star_as_star_list(sequence *seq);
char *	fit_sequence_get_image_filename(sequence *seq, int index, char *name_buffer, gboolean add_fits_ext);
char *	fit_sequence_get_image_filename_prefixed(sequence *seq, const char *prefix, int index);
char *	get_possible_image_filename(sequence *seq, int image_number, char *name_buffer);
int	get_index_and_basename(const char *filename, char **basename, int *index, int *fixed, const gchar *com_ext);
void	remove_prefixed_sequence_files(sequence *seq, const char *prefix);
void	remove_prefixed_star_files(sequence *seq, const char *prefix);
void	initialize_sequence(sequence *seq, gboolean is_zeroed);
void	free_sequence(sequence *seq, gboolean free_seq_too);
void	free_photometry_set(sequence *seq, int set);
void	close_sequence(int loading_another);
gboolean check_seq_is_comseq(sequence *seq);
gboolean check_seq_is_variable(sequence *seq);
gboolean sequence_is_loaded();
gboolean check_cachefile_date(sequence *seq, int index, const gchar *star_filename) ;

typedef enum {
	ORIGINAL_FRAME,
	FOLLOW_STAR_FRAME,
	REGISTERED_FRAME
} framing_mode;

struct seqpsf_args {
	gboolean for_photometry;
	super_bool allow_use_as_regdata;
	framing_mode framing;
	char bayer_pattern[FLEN_VALUE];

	/* The seqpsf result for each image, list of seqpsf_data */
	GSList *list;
};

struct seqpsf_data {
	int image_index;
	psf_star *psf;
	double exposure;
};

int	sequence_find_refimage(sequence *seq);
void check_or_allocate_regparam(sequence *seq, int layer);
void set_shifts(sequence *seq, int frame, int layer, double shiftx, double shifty, gboolean data_is_top_down);
void cum_shifts(regdata *regparam, int frame, double shiftx, double shifty);
gboolean test_regdata_is_valid_and_shift(sequence *seq, int reglayer);
sequence *create_internal_sequence(int size);
void	internal_sequence_set(sequence *seq, int index, fits *fit);
int	internal_sequence_find_index(sequence *seq, const fits *fit);
fits	*internal_sequence_get(sequence *seq, int index);
gboolean sequence_is_rgb(sequence *seq);
gboolean	enforce_area_in_image(rectangle *area, sequence *seq, int index);

int seqpsf(sequence *seq, int layer, gboolean for_registration, gboolean regall,
		framing_mode framing, gboolean run_in_thread, gboolean no_GUI);
int seqpsf_image_hook(struct generic_seq_args *args, int out_index, int index, fits *fit, rectangle *area, int threads);
void free_reference_image();

/* in export.c now */
void	update_export_crop_label();

size_t get_max_seq_dimension(sequence *seq, int *rx, int *ry);
int compute_nb_images_fit_memory(sequence *seq, double factor, gboolean force_float, unsigned int *MB_per_orig_image, unsigned int *MB_per_scaled_image, unsigned int *max_mem_MB);

int compute_nb_images_fit_memory_from_fit(fits *fit, double factor, gboolean force_float, unsigned int *MB_per_orig_image, unsigned int *MB_per_scaled_image, unsigned int *max_mem_MB);

void fix_selnum(sequence *seq, gboolean warn);

gboolean sequence_ref_has_wcs(sequence *seq);
struct wcsprm *get_wcs_ref(sequence *seq);

gboolean sequence_drifts(sequence *seq, int reglayer, int threshold);

void clean_sequence(sequence *seq, gboolean cleanreg, gboolean cleanstat, gboolean cleansel);

#endif
