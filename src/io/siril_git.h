#ifndef _GITSCRIPTS_H_
#define _GITSCRIPTS_H_

#ifdef HAVE_LIBGIT2
int auto_update_gitscripts(gboolean sync);
int auto_update_gitspcc(gboolean sync);
void on_disable_gitscripts();
void fill_script_repo_tree(gboolean as_idle);
int reset_repository(const gchar *local_path);
gboolean is_scripts_repo_cloned();
gboolean is_spcc_repo_cloned();
gboolean fill_spcc_widgets_in_thread(gpointer user_data);
void gui_repo_scripts_mutex_lock();
void gui_repo_scripts_mutex_unlock();
gpointer update_repo_scripts_list_and_menu_in_thread();
gchar *get_script_content_string_from_file_revision(const char *filepath,
													int file_revisions_to_backtrack,
													size_t *content_size,
													gchar **commit_message,
													size_t *message_size);
#else
void hide_git_widgets();
#endif

#endif

