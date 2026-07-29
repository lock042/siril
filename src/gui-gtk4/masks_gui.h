#ifndef _MASK_FROM_IMAGE_CALLBACKS_H_
#define _MASK_FROM_IMAGE_CALLBACKS_H_

#include <gtk/gtk.h>

/* Dialog show handler */
void on_mask_from_image_dialog_show(GtkWidget *widget, gpointer user_data);

/* Button click handlers */
void on_mask_from_image_close_clicked(GtkButton *button, gpointer user_data);
void on_mask_from_image_apply_clicked(GtkButton *button, gpointer user_data);

/* ComboBox change handlers */
void on_combo_mask_from_image_type_changed(GObject *obj, GParamSpec *pspec, gpointer user_data);
void on_combo_mask_luminance_type_changed(GObject *obj, GParamSpec *pspec, gpointer user_data);

/* Color picker */
void mask_color_handle_image_click(int x, int y);

/* Idle function: sync the mask-enable toggle to state (for gui_iface impl). */
gboolean set_mask_active_idle(gpointer p);

/* Re-read the layer stack into every mask dialog's "Target" dropdown.
 * Called from the Layers panel's refresh, because those dropdowns were built
 * when the dialog was shown and go stale the moment a layer is added,
 * duplicated, removed or renamed — a dialog left open across a duplicate
 * offered the old stack, so the mask landed on the wrong layer.  The current
 * selection is preserved by item id when the layer is still there. */
void masks_gui_refresh_target_combos(void);

#endif /* _MASK_FROM_IMAGE_CALLBACKS_H_ */
