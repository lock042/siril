#ifndef SINGLE_H_
#define SINGLE_H_

void close_single_image();
void free_image_data();
int create_uniq_from_gfit(char *filename, gboolean exists);
int read_single_image(const char* filename, fits *dest, char **realname_out, gboolean allow_sequences, gboolean *is_sequence, gboolean allow_dialogs, gboolean force_float, gboolean no_debayer);
/* Set/query whether the shared image readers should report decode progress on
 * the progress bar.  Active only during an interactive single-image open. */
void set_read_progress_active(gboolean active);
gboolean read_progress_active(void);
gboolean end_open_single_image(gpointer arg);
int open_single_image(const char* filename);
/* Guard-free open body, for callers already owning the processing thread. */
int open_single_image_internal(const char* filename);
/* Threaded single-image open: the worker decodes into a private fits
 * (read_new_single_image), the main thread swaps it into gfit under the writer
 * lock (install_new_single_image).  Sequences are routed to the synchronous
 * path via single_image_path_is_sequence(). */
gboolean single_image_path_is_sequence(const char *filename);
fits *read_new_single_image(const char *filename, char **realname_out, int *retval);
int install_new_single_image(fits *newfit, char *realname);
gboolean open_single_image_from_gfit(gpointer user_data);
gboolean update_single_image_from_gfit(gpointer user_data);
int image_find_minmax(fits *fit);
double fit_get_max(fits *fit, int layer);
double fit_get_min(fits *fit, int layer);
void init_layers_hi_and_lo_values(sliders_mode force_minmax);
void unique_free_preprocessing_data(single *uniq);
int single_image_is_loaded();

void adjust_cutoff_from_updated_gfit();		// was level_adjust(), deprecated too
gboolean end_gfit_operation(gpointer data);
void gfit_modified_update_gui();		// to be called after all gfit modifications
void notify_gfit_data_modified();

gboolean enforce_area_in_fits(fits *fit, rectangle *area);

#endif
