#ifndef _SIRIL_PYTHON_H
#define _SIRIL_PYTHON_H

void init_python(void);
void finalize_python(void);

gpointer run_python_script_from_file(gpointer p);
gpointer run_python_script_from_mem(gpointer p);

#endif
