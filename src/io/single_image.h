#ifndef SINGLE_H_
#define SINGLE_H_

void close_single_image();
void free_image_data();
int create_uniq_from_gfit(char *filename, gboolean exists);
int read_single_image(const char* filename, fits *dest, char **realname_out, gboolean allow_sequences, gboolean *is_sequence, gboolean allow_dialogs, gboolean force_float);
gboolean end_open_single_image(gpointer arg);
int open_single_image(const char* filename);
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
void notify_gfit_modified();		// to be called after all gfit modifications

#endif
