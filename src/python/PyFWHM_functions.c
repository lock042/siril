/*
* This file is part of Siril, an astronomy image processor.
* Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
* Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
* Reference site is https://siril.org
*
* Siril is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Siril is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Siril. If not, see <http://www.gnu.org/licenses/>.
*/

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "core/siril.h"
#include "core/settings.h"
#include "algos/PSF.h"
#include "python/siril_python.h"

PyObject *PyFWHM_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    PyFWHMObject *self;
    self = (PyFWHMObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->fwhm = NULL;
        self->should_free = 0;
    }
    return (PyObject *)self;
}

int PyFWHM_init(PyFWHMObject *self, PyObject *args, PyObject *kwds) {
    self->fwhm = calloc(1, sizeof(struct fwhm_struct));
    if (self->fwhm == NULL) {
        PyErr_SetString(PyExc_MemoryError, N_("Failed to allocate memory for fwhm_struct"));
        return -1;
    }
    self->should_free = 1;
    return 0;
}

PyObject *PyFWHM_FromExisting(struct fwhm_struct *fwhm) {
    PyFWHMObject *obj = (PyFWHMObject *)PyFWHM_new(&PyFWHMType, NULL, NULL);
    if (obj != NULL) {
        obj->fwhm = fwhm;
        obj->should_free = 0;  // We don't own this data
    }
    return (PyObject *)obj;
}

void PyFWHM_dealloc(PyFWHMObject *self) {
    if (self->fwhm != NULL) {
        if (self->should_free) {
            free(self->fwhm->star_name);
            free(self->fwhm->units);
            // Free any other dynamically allocated members here
        }
        free(self->fwhm);
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

// PyFWHMType methods and getters
PyObject* PyFWHM_get_star_name(PyFWHMObject *self, void *closure) {
	return PyUnicode_FromString(self->fwhm->star_name);
}

PyObject* PyFWHM_get_profile_type(PyFWHMObject *self, void *closure) {
	return PyUnicode_FromString(self->fwhm->profile == PSF_GAUSSIAN ? "Gaussian" : "Moffat");
}

PyObject* PyFWHM_get_R(PyFWHMObject *self, void *closure) {
	return PyLong_FromLong(self->fwhm->R);
}

PyObject* PyFWHM_get_B(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->B);
}

PyObject* PyFWHM_get_A(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->A);
}

PyObject* PyFWHM_get_x0(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->x0);
}

PyObject* PyFWHM_get_y0(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->y0);
}

PyObject* PyFWHM_get_sx(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->sx);
}

PyObject* PyFWHM_get_sy(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->sy);
}

PyObject* PyFWHM_get_fwhmx(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->fwhmx);
}

PyObject* PyFWHM_get_fwhmy(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->fwhmy);
}

PyObject* PyFWHM_get_fwhmx_arcsec(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->fwhmx_arcsec);
}

PyObject* PyFWHM_get_fwhmy_arcsec(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->fwhmy_arcsec);
}

PyObject* PyFWHM_get_angle(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->angle);
}

PyObject* PyFWHM_get_rmse(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->rmse);
}

PyObject* PyFWHM_get_sat(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->sat);
}

PyObject* PyFWHM_has_saturated(PyFWHMObject *self, void *closure) {
	return PyBool_FromLong(self->fwhm->has_saturated ? 1 : 0);
}

PyObject* PyFWHM_get_xpos(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->xpos);
}

PyObject* PyFWHM_get_beta(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->beta);
}

PyObject* PyFWHM_get_ypos(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->ypos);
}

PyObject* PyFWHM_get_mag(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->mag);
}

PyObject* PyFWHM_get_Bmag(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->Bmag);
}

PyObject* PyFWHM_get_SNR(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->SNR);
}

PyObject* PyFWHM_get_layer(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->layer);
}

PyObject* PyFWHM_get_ra(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->ra);
}

PyObject* PyFWHM_get_dec(PyFWHMObject *self, void *closure) {
	return PyFloat_FromDouble(self->fwhm->dec);
}
