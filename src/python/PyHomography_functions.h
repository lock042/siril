#ifndef PYHOMOGRAPHY_FUNCTIONS_H
#define PYHOMOGRAPHY_FUNCTIONS_H

#include "python/siril_python.h"

void PyHomography_dealloc(PyHomographyObject *self);
PyObject *PyHomography_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int PyHomography_init(PyHomographyObject *self, PyObject *args, PyObject *kwds);
PyObject *PyHomography_FromExisting(Homography *homography);
PyBufferProcs PyHomography_as_buffer;

PyObject* PyHomography_get_h00(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_h01(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_h02(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_h10(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_h11(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_h12(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_h20(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_h21(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_h22(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_pair_matched(PyHomographyObject *self, void *closure);
PyObject* PyHomography_get_Inliers(PyHomographyObject *self, void *closure);

#endif
