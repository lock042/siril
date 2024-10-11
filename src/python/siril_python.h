#ifndef _SIRIL_PYTHON_H
#define _SIRIL_PYTHON_H

gpointer init_python(gpointer user_data);
void finalize_python(void);

gboolean run_python_script_from_file(gpointer p);
gboolean run_python_script_from_mem(gpointer p);
void run_python_script_in_python_thread(const char *script, gboolean from_file);

#endif
