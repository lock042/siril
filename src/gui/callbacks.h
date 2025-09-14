#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <sys/time.h>
#include "core/siril.h"	// for sliders_mode

void handle_owner_change(GtkClipboard *clipboard, GdkEvent *event, gpointer data);
void on_press_seq_field();
gboolean launch_clipboard_survey(gpointer user_data);

// Region of Interest processing
typedef void (*ROICallback)();
int populate_roi();
gpointer on_set_roi();
gpointer on_clear_roi();
void add_roi_callback(ROICallback func);
void remove_roi_callback(ROICallback func);
void update_roi_config();

gboolean is_gui_ready();
void lock_roi_mutex();
void unlock_roi_mutex();
void roi_supported(gboolean state);
void initialize_all_GUI(gchar *files);
void siril_set_theme(int active);
void load_prefered_theme(gint theme);
void set_cutoff_sliders_max_values();		// was set_upper_minmax
void set_cutoff_sliders_values();		// was set_ranges
void set_sliders_value_to_gfit();
void set_accel_map(const gchar * const *accelmap);
void initialize_display_mode();
void set_display_mode();
void set_unlink_channels(gboolean unlinked);
void adjust_exclude(int n, gboolean changed);
void adjust_sellabel();
gpointer update_seq_gui_idle_thread_func(gpointer data);
gboolean set_GUI_CWD(gpointer user_data);
void set_icon_entry(GtkEntry *entry, gchar *string);
gboolean update_MenuItem(gpointer user_data);
void sliders_mode_set_state(sliders_mode);
display_mode get_display_mode_from_menu();
int copy_rendering_settings();

void clear_sampling_setting_box();
void set_GUI_CAMERA();

int match_drawing_area_widget(const GtkWidget *drawing_area, gboolean allow_rgb);
void update_display_selection();
void update_display_fwhm();
void display_filename();
gboolean set_precision_switch(gpointer user_data);
void set_layers_for_assign();
int set_layers_for_registration();
void show_dialog(const char *text, const char *title, const char *icon);
void show_txt_and_data_dialog(const char *text, const char *data, const char *title, const char *icon);
void show_data_dialog(char *text, char *title, gchar *parent, gchar *extra_button);
GtkWindow *siril_get_active_window();
void initialize_FITS_name_entries();

void adjust_vport_size_to_image();
void set_output_filename_to_sequence_name();
gboolean close_tab(gpointer user_data);
void activate_tab(int vport);
gboolean init_right_tab(gpointer user_data);

void update_prepro_interface(gboolean allow_debayer);

void on_treeview_selection_convert_changed(GtkTreeSelection *treeselection, gpointer user_data);
void update_statusbar_convert();

gboolean update_spinCPU(gpointer user_data);

gpointer update_scripts(gpointer user_data);
gpointer update_spcc(gpointer user_data);
gpointer initialize_scripts(gpointer user_data);
gpointer initialize_spcc(gpointer user_data);

gboolean save_main_window_state(gpointer user_data);
gboolean load_main_window_state(gpointer user_data);
GPid show_child_process_selection_dialog(GSList *children);
gboolean set_seq_browser_active(gpointer user_data);
gboolean siril_quit(void);

/* for image_display */
void set_viewer_mode_widgets_sensitive(gboolean sensitive);

int seq_qphot(sequence *seq, int layer);

/*****************************************************************************
*      P U B L I C      C A L L B A C K      F U N C T I O N S               *
 ****************************************************************************/
void setup_stretch_sliders();
void on_radiobutton_minmax_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_radiobutton_hilo_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_radiobutton_user_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_max_entry_changed(GtkEditable *editable, gpointer user_data);
void on_min_entry_changed(GtkEditable *editable, gpointer user_data);

void on_seqproc_entry_changed (GtkComboBox *widget,	gpointer user_data);
void on_excludebutton_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_ref_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data);

void on_spin_w_changed(GtkSpinButton *spinbutton, gpointer user_data);

void on_check_button_pref_bias_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_check_button_pref_bias_bis_toggled(GtkToggleButton *togglebutton, gpointer user_data);

void on_focal_entry_changed(GtkEditable *editable, gpointer user_data);
void on_pitchX_entry_changed(GtkEditable *editable, gpointer user_data);
void on_pitchY_entry_changed(GtkEditable *editable, gpointer user_data);
void on_combobinning_changed(GtkComboBox *box, gpointer user_data);

#endif
