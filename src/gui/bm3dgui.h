#ifndef SRC_GUI_BM3D_H_
#define SRC_GUI_BM3D_H_

void on_bm3d_dialog_show(GtkWidget *widget, gpointer user_data);
void on_bm3d_cancel_clicked(GtkButton *button, gpointer user_data);

void on_toggle_bm3d_pre_median_filter_toggled(GtkCheckButton *button, gpointer user_data);
void on_spin_bm3d_modulation_value_changed(GtkSpinButton *button, gpointer user_data);

void on_bm3d_apply_clicked(GtkButton *button, gpointer user_data);

#endif /* SRC_GUI_BM3D_H_ */
