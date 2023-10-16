#ifndef _GITSCRIPTS_H_
#define _GITSCRIPTS_H_

#ifdef HAVE_LIBGIT2
int auto_update_gitscripts(gboolean sync);
void on_disable_gitscripts();
void fill_script_repo_list(gboolean as_idle);
#else
void hide_git_widgets();
#endif

#endif

