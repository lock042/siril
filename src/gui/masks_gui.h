#ifndef _MASK_FROM_IMAGE_CALLBACKS_H_
#define _MASK_FROM_IMAGE_CALLBACKS_H_

#include <gtk/gtk.h>

/* Dialog show handler */
void on_mask_from_image_dialog_show(GtkWidget *widget, gpointer user_data);

/* Button click handlers */
void on_mask_from_image_close_clicked(GtkButton *button, gpointer user_data);
void on_mask_from_image_apply_clicked(GtkButton *button, gpointer user_data);

/* ComboBox change handlers */
void on_combo_mask_from_image_type_changed(GtkComboBox *combo, gpointer user_data);
void on_combo_mask_luminance_type_changed(GtkComboBox *combo, gpointer user_data);

#endif /* _MASK_FROM_IMAGE_CALLBACKS_H_ */
