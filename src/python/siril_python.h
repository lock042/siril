#ifndef _SIRIL_PYTHON_H
#define _SIRIL_PYTHON_H

#include <Python.h>

// Object type definitions
typedef struct {
	PyObject_HEAD
	fits *fit;
	gchar* filename;
	int should_free_data;  // Flag to indicate if the data itself should be freed
	// note: if should_free_data is 0, the struct will not be freed either
	int should_free;  // Flag to indicate if the fits pointer itself should be freed
	Py_ssize_t ref_count;
} PyFits;

// PySeqObject
typedef struct {
	PyObject_HEAD
	struct sequ *seq;
	int should_free;
	Py_ssize_t ref_count;
} PySeqObject;

// PyImgDataObject
typedef struct {
	PyObject_HEAD
	imgdata *img;
	int should_free;
	PySeqObject *seq;
} PyImgDataObject;

// PyImStatsObject
typedef struct {
	PyObject_HEAD
	imstats *stats;
	int should_free;
	PyObject *parent;  // This can be either PySequenceObject or PyFitsObject
	char parent_type;  // 'S' for Sequence, 'F' for Fits
} PyImStatsObject;

// PyHomographyObject
typedef struct {
	PyObject_HEAD
	Homography *homography;
	int should_free;
} PyHomographyObject;

// PyRegDataObject
typedef struct {
	PyObject_HEAD
	regdata *reg;
	int should_free;
	Py_ssize_t ref_count;
	PySeqObject *seq;
} PyRegDataObject;

// PyFWHMObject
typedef struct {
	PyObject_HEAD
	struct fwhm_struct *fwhm;
	int should_free;
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
