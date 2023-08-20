#ifndef _GITSCRIPTS_H_
#define _GITSCRIPTS_H_

int update_gitscripts(gboolean sync);
void fill_script_repo_list(gboolean as_idle);
void on_pref_use_gitscripts_toggled(GtkToggleButton *button, gpointer user_data);


#endif

