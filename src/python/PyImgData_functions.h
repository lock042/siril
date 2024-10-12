#ifndef PYIMGDATA_FUNCTIONS_H
#define PYIMGDATA_FUNCTIONS_H

#include "python/siril_python.h"

PyObject *PyImgData_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int PyImgData_init(PyImgDataObject *self, PyObject *args, PyObject *kwds);
void PyImgData_dealloc(PyImgDataObject *self);
PyObject *PyImgData_FromExisting(imgdata *data, PyObject *seq);

PyObject* PyImgData_get_filenum(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_incl(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_date_obs(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_airmass(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_rx(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_ry(PyImgDataObject *self, void *closure);

#endif
