#ifndef PYIMGDATA_FUNCTIONS_H
#define PYIMGDATA_FUNCTIONS_H

#include "python/siril_python.h"

PyObject* PyImgData_get_filenum(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_incl(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_date_obs(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_airmass(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_rx(PyImgDataObject *self, void *closure);
PyObject* PyImgData_get_ry(PyImgDataObject *self, void *closure);

#endif
