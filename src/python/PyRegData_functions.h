#ifndef PYREGDATA_FUNCTIONS_H
#define PYREGDATA_FUNCTIONS_H

#include "python/siril_python.h"

void PyRegData_dealloc(PyRegDataObject *self);
PyObject *PyRegData_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int PyRegData_init(PyRegDataObject *self, PyObject *args, PyObject *kwds);
PyObject *PyRegData_FromExisting(regdata *data, PyObject *seq);

PyObject* PyRegData_get_fwhm(PyRegDataObject *self, void *closure);
PyObject* PyRegData_get_wfwhm(PyRegDataObject *self, void *closure);
PyObject* PyRegData_get_weighted_fwhm(PyRegDataObject *self, void *closure);
PyObject* PyRegData_get_roundness(PyRegDataObject *self, void *closure);
PyObject* PyRegData_get_quality(PyRegDataObject *self, void *closure);
PyObject* PyRegData_get_bg(PyRegDataObject *self, void *closure);
PyObject* PyRegData_get_number_of_stars(PyRegDataObject *self, void *closure);
PyObject* PyRegData_get_homography(PyRegDataObject *self, PyObject *args);

#endif
