#ifndef SIRIL_FILE_CHOOSER_BUTTON_H
#define SIRIL_FILE_CHOOSER_BUTTON_H

#include <gtk/gtk.h>

G_BEGIN_DECLS

#define SIRIL_TYPE_FILE_CHOOSER_BUTTON siril_file_chooser_button_get_type()
G_DECLARE_FINAL_TYPE(SirilFileChooserButton, siril_file_chooser_button, SIRIL, FILE_CHOOSER_BUTTON, GtkButton)

SirilFileChooserButton *siril_file_chooser_button_new(const gchar *title, GtkFileChooserAction action);
void siril_file_chooser_button_set_width_chars(SirilFileChooserButton *button, gint width_chars);
void siril_file_chooser_set_current_folder(SirilFileChooserButton *button, const gchar *folder_path);

G_END_DECLS

#endif /* SIRIL_FILE_CHOOSER_BUTTON_H */

/*
 * Functions still needed:
 *
 * siril_filechooser_button_unselect_all
 * siril_file_chooser_button_set_current_folder
 * siril_file_chooser_button_get_current_folder
 * siril_file_chooser_button_add_filter
 * siril_file_chooser_button_set_filter
 * siril_file_chooser_button_set_select_multiple
 * siril_file_chooser_button_set_current_name
 * siril_file_chooser_button_set_do_overwrite_confirmation
 * siril_file_chooser_button_set_filename
 * siril_file_chooser_button_get_filename
 *
