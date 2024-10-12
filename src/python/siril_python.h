#ifndef _SIRIL_PYTHON_H
#define _SIRIL_PYTHON_H

#include <Python.h>

typedef struct {
	PyObject_HEAD
	fits *fit;
	int should_free_data;  // Flag to indicate if the data itself should be freed
	// note: if should_free_data is 0, the struct will not be freed either
	int should_free;  // Flag to indicate if the fits pointer itself should be freed
} PyFits;

// Object type definitions
PyTypeObject PyFitsType;
PyTypeObject PySeqType;
PyTypeObject PyImgDataType;

// Common helper functions
PyObject* gdatetime_to_pydatetime(GDateTime *gdt);

// Interpreter initialization and finalization functions
gpointer init_python(gpointer user_data);
void finalize_python(void);

void run_python_script_in_python_thread(const char *script, gboolean from_file);

#endif
