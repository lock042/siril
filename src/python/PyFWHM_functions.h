#ifndef PYFWHM_FUNCTIONS_H
#define PYFWHM_FUNCTIONS_H

#include "python/siril_python.h"

void PyFWHM_dealloc(PyFWHMObject *self);
PyObject *PyFWHM_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int PyFWHM_init(PyFWHMObject *self, PyObject *args, PyObject *kwds);
PyObject *PyFWHM_FromExisting(struct fwhm_struct *fwhm);

PyObject* PyFWHM_get_star_name(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_R(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_B(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_A(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_x0(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_y0(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_sx(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_sy(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_fwhmx(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_fwhmy(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_fwhmx_arcsec(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_fwhmy_arcsec(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_angle(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_rmse(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_sat(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_has_saturated(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_xpos(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_ypos(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_mag(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_Bmag(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_SNR(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_layer(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_ra(PyFWHMObject *self, void *closure);
PyObject* PyFWHM_get_dec(PyFWHMObject *self, void *closure);

#endif
