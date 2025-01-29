#ifndef SRC_GUI_SCRIPT_MENU_H_
#define SRC_GUI_SCRIPT_MENU_H_

#define SCRIPT_EXT "ssf"
#define PYSCRIPT_EXT "py"
#define PYCSCRIPT_EXT "pyc"

GSList *get_list_from_preferences_dialog();
GSList *set_list_to_preferences_dialog(GSList *list);
int initialize_script_menu(gboolean verbose);
int refresh_script_menu(gboolean verbose);
int refresh_scripts(gboolean update_list, gchar **error);
void siril_get_on_script_pages();
void script_widgets_enable(gboolean status);
gboolean script_widgets_idle(gpointer user_data);
gboolean accept_script_warning_dialog();
gboolean test_last_subdir(const gchar *path, const gchar *expected_subdir);

#endif /* SRC_GUI_SCRIPT_MENU_H_ */
