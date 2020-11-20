#ifndef SRC_GUI_SCRIPT_MENU_H_
#define SRC_GUI_SCRIPT_MENU_H_

GSList *get_list_from_preferences();
int initialize_script_menu(gboolean UpdateScriptPath);
int refresh_scripts(gchar **error);
void siril_get_on_script_pages();

#endif /* SRC_GUI_SCRIPT_MENU_H_ */
