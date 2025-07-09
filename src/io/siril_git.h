#ifndef _GITSCRIPTS_H_
#define _GITSCRIPTS_H_

#ifdef HAVE_LIBGIT2
int auto_update_gitscripts(gboolean sync);
int auto_update_gitspcc(gboolean sync);
void on_disable_gitscripts();
void fill_script_repo_tree(gboolean as_idle);
int reset_repository(const gchar *local_path);
int preview_scripts_update(GString **git_pending_commit_buffer);
int preview_spcc_update(GString **git_pending_commit_buffer);
gboolean is_scripts_repo_cloned();
gboolean is_spcc_repo_cloned();
gboolean fill_spcc_widgets_in_thread(gpointer user_data);
void gui_repo_scripts_mutex_lock();
void gui_repo_scripts_mutex_unlock();
#else
void hide_git_widgets();
#endif

#endif

