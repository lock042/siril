#ifndef _SIRIL_PYTHON_H
#define _SIRIL_PYTHON_H

#include <Python.h>

// Object type definitions
typedef struct {
	PyObject_HEAD
	fits *fit;
	int should_free_data;  // Flag to indicate if the data itself should be freed
	// note: if should_free_data is 0, the struct will not be freed either
	int should_free;  // Flag to indicate if the fits pointer itself should be freed
} PyFits;

// PyImgDataObject
typedef struct {
    PyObject_HEAD
    imgdata *img;
} PyImgDataObject;

// PySeqObject
typedef struct {
    PyObject_HEAD
    struct sequ *seq;
} PySeqObject;

// PyImStatsObject
typedef struct {
    PyObject_HEAD
    imstats *stats;
} PyImStatsObject;

// PyHomographyObject
typedef struct {
    PyObject_HEAD
    Homography *homography;
} PyHomographyObject;

// PyRegDataObject
typedef struct {
    PyObject_HEAD
    regdata *reg;
} PyRegDataObject;

// PyFWHMObject
typedef struct {
    PyObject_HEAD
    struct fwhm_struct *fwhm;
} PyFWHMObject;

PyTypeObject PyFitsType;
PyTypeObject PySeqType;
PyTypeObject PyImgDataType;
PyTypeObject PyImStatsType;
PyTypeObject PyHomographyType;
PyTypeObject PyRegDataType;
PyTypeObject PyFWHMType;

// Common helper functions
PyObject* gdatetime_to_pydatetime(GDateTime *gdt);

// Interpreter initialization and finalization functions
gpointer init_python(gpointer user_data);
void finalize_python(void);

void run_python_script_in_python_thread(const char *script, gboolean from_file);

#endif
