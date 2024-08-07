#ifndef SRC_GUI_SCRIPT_MENU_H_
#define SRC_GUI_SCRIPT_MENU_H_

#define SCRIPT_EXT ".ssf"

GSList *get_list_from_preferences_dialog();
GSList *set_list_to_preferences_dialog(GSList *list);
int initialize_script_menu(gboolean verbose);
int refresh_script_menu(gboolean verbose);
int refresh_scripts(gboolean update_list, gchar **error);
void siril_get_on_script_pages();
void script_widgets_enable(gboolean status);
gboolean script_widgets_idle(gpointer user_data);

#endif /* SRC_GUI_SCRIPT_MENU_H_ */
