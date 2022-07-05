#ifndef SRC_FILTERS_STARNETGUI_H_
#define SRC_FILTERS_STARNETGUI_H_

void on_starnet_cancel_clicked(GtkButton *button, gpointer user_data);
void on_starnet_execute_clicked(GtkButton *button, gpointer user_data);
void on_spin_starnet_stride_changed(GtkSpinButton *button, gpointer user_data);
void on_starnet_stretch_toggled(GtkToggleButton *button, gpointer user_data);
void on_starnet_upsample_toggled(GtkToggleButton *button, gpointer user_data);
void on_starnet_starmask_toggled(GtkToggleButton *button, gpointer user_data);


#endif /* SRC_FILTERS_STARNET_H_ */
