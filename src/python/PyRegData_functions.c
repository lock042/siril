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
#include "python/siril_python.h"
#include "python/PyHomography_functions.h"

PyObject *PyRegData_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	PyRegDataObject *self;
	self = (PyRegDataObject *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->reg = NULL;
		self->should_free = 0;
	}
	return (PyObject *)self;
}

int PyRegData_init(PyRegDataObject *self, PyObject *args, PyObject *kwds) {
	self->reg = (regdata *)calloc(1, sizeof(regdata));
	if (self->reg == NULL) {
		PyErr_SetString(PyExc_MemoryError, N_("Failed to allocate memory for regdata"));
		return -1;
	}
	self->should_free = 1;
	return 0;
}

void PyRegData_dealloc(PyRegDataObject *self) {
	if (self->reg != NULL && self->should_free) {
		if (self->reg->fwhm_data != NULL) {
			free(self->reg->fwhm_data);
		}
		free(self->reg);
	}
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *PyRegData_FromExisting(regdata *data, PyObject *seq) {
	PyRegDataObject *obj = (PyRegDataObject *)PyRegData_new(&PyRegDataType, NULL, NULL);
	if (obj != NULL) {
		obj->reg = data;
		obj->should_free = 0;
        Py_INCREF(seq);
        obj->seq = (PySeqObject *)seq;
	}
	return (PyObject *)obj;
}

// PyRegDataType methods and getters
PyObject* PyRegData_get_fwhm(PyRegDataObject *self, void *closure) {
	return PyFloat_FromDouble(self->reg->fwhm);
}

PyObject* PyRegData_get_wfwhm(PyRegDataObject *self, void *closure) {
	return PyFloat_FromDouble(self->reg->fwhm);
}

PyObject* PyRegData_get_roundness(PyRegDataObject *self, void *closure) {
	return PyFloat_FromDouble(self->reg->fwhm);
}

PyObject* PyRegData_get_quality(PyRegDataObject *self, void *closure) {
	return PyFloat_FromDouble(self->reg->fwhm);
}

PyObject* PyRegData_get_bg(PyRegDataObject *self, void *closure) {
	return PyFloat_FromDouble(self->reg->fwhm);
}

PyObject* PyRegData_get_number_of_stars(PyRegDataObject *self, void *closure) {
	return PyFloat_FromDouble(self->reg->fwhm);
}

PyObject* PyRegData_get_homography(PyRegDataObject *self, PyObject *Py_UNUSED(args))
{
    if (self->reg == NULL) {
        PyErr_SetString(PyExc_AttributeError, N_("RegData object has no data"));
        return NULL;
    }

    // Create a new PyHomographyObject using the from_existing method
    PyObject *homography_obj = PyHomography_FromExisting(&(self->reg->H));

    if (homography_obj == NULL) {
        PyErr_SetString(PyExc_RuntimeError, N_("Failed to create Homography object"));
        return NULL;
    }

    return homography_obj;
}
