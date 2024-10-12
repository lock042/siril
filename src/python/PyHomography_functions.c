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

PyObject *PyHomography_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    PyHomographyObject *self;
    self = (PyHomographyObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->homography = NULL;
        self->should_free = 0;
    }
    return (PyObject *)self;
}

int PyHomography_init(PyHomographyObject *self, PyObject *args, PyObject *kwds) {
    self->homography = (Homography *)calloc(1, sizeof(Homography));
    if (self->homography == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for Homography");
        return -1;
    }
    self->should_free = 1;
    return 0;
}

void PyHomography_dealloc(PyHomographyObject *self) {
    if (self->homography != NULL && self->should_free) {
        free(self->homography);
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *PyHomography_FromExisting(Homography *homography) {
    PyHomographyObject *obj = (PyHomographyObject *)PyHomography_new(&PyHomographyType, NULL, NULL);
    if (obj != NULL) {
        obj->homography = homography;
        obj->should_free = 0;  // We don't own this data
    }
    return (PyObject *)obj;
}

static int PyHomography_getbuffer(PyObject *obj, Py_buffer *view, int flags) {
	PyHomographyObject *self = (PyHomographyObject *)obj;
	Homography *H = self->homography;

	if (!H) {
		PyErr_SetString(PyExc_ValueError, "Homography object is NULL");
		return -1;
	}

	// Determine dimensionality
	view->ndim = 2;
	view->readonly = 1;  // Assume writable for now

	// Allocate memory for shape and strides
	Py_ssize_t *shape = PyMem_Malloc(2 * sizeof(Py_ssize_t));
	Py_ssize_t *strides = PyMem_Malloc(2 * sizeof(Py_ssize_t));

	if (!shape || !strides) {
		PyErr_NoMemory();
		PyMem_Free(shape);
		PyMem_Free(strides);
		return -1;
	}

	// Set shape
	shape[0] = shape[1] = 3;
	view->shape = shape;

	// Set buffer and itemsize based on data type
	view->buf = self->homography;
	view->itemsize = sizeof(double);
	view->format = "d";  // double

	// Calculate total length
	view->len = view->itemsize * 9;

	strides[0] = strides[1] = 3 * view->itemsize;
	view->strides = strides;

	view->suboffsets = NULL;
	view->internal = NULL;

	return 0;
}

static void PyHomography_releasebuffer(PyObject *obj, Py_buffer *view) {
	// Free the shape and strides arrays we allocated in getbuffer
	PyMem_Free(view->shape);
	PyMem_Free(view->strides);

	// Clear the Py_buffer struct without freeing the actual data
	// which is still owned by the PyHomographyObject
	memset(view, 0, sizeof(Py_buffer));
}

PyBufferProcs PyHomography_as_buffer = {
	(getbufferproc)PyHomography_getbuffer,
	(releasebufferproc)PyHomography_releasebuffer  // Standard release function
};

// PyHomographyType methods and getters
PyObject* PyHomography_get_h00(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h00);
}

PyObject* PyHomography_get_h01(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h01);
}

PyObject* PyHomography_get_h02(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h02);
}

PyObject* PyHomography_get_h10(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h10);
}

PyObject* PyHomography_get_h11(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h11);
}

PyObject* PyHomography_get_h12(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h12);
}

PyObject* PyHomography_get_h20(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h20);
}

PyObject* PyHomography_get_h21(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h21);
}

PyObject* PyHomography_get_h22(PyHomographyObject *self, void *closure) {
	return PyFloat_FromDouble(self->homography->h22);
}

PyObject* PyHomography_get_pair_matched(PyHomographyObject *self, void *closure) {
	return PyLong_FromLong(self->homography->pair_matched);
}

PyObject* PyHomography_get_Inliers(PyHomographyObject *self, void *closure) {
	return PyLong_FromLong(self->homography->Inliers);
}
