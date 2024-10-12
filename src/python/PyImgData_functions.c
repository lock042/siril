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

PyObject *PyImgData_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	PyImgDataObject *self;
	self = (PyImgDataObject *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->img = NULL;
		self->should_free = 0;
	}
	return (PyObject *)self;
}

int PyImgData_init(PyImgDataObject *self, PyObject *args, PyObject *kwds) {
	self->img = (imgdata *)calloc(1, sizeof(imgdata));
	if (self->img == NULL) {
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for imgdata");
		return -1;
	}
	self->should_free = 1;
	return 0;
}

void PyImgData_dealloc(PyImgDataObject *self) {
	if (self->img != NULL && self->should_free) {
		if (self->img->date_obs != NULL) {
			g_date_time_unref(self->img->date_obs);
		}
		free(self->img);
	}
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *PyImgData_FromExisting(imgdata *data, PyObject *seq) {
	PyImgDataObject *obj = (PyImgDataObject *)PyImgData_new(&PyImgDataType, NULL, NULL);
	if (obj != NULL) {
		obj->img = data;
		obj->should_free = 0;
		Py_INCREF(seq);
		obj->seq = (PySeqObject *)seq;
	}
	return (PyObject *)obj;
}

PyObject* PyImgData_get_filenum(PyImgDataObject *self, void *closure) {
	return PyLong_FromLong(self->img->filenum);
}

PyObject* PyImgData_get_incl(PyImgDataObject *self, void *closure) {
	return PyBool_FromLong(self->img->incl);
}

PyObject* PyImgData_get_date_obs(PyImgDataObject *self, void *closure) {
	// Assuming you have a function to convert GDateTime to PyDateTime
	return gdatetime_to_pydatetime(self->img->date_obs);
}

PyObject* PyImgData_get_airmass(PyImgDataObject *self, void *closure) {
	return PyFloat_FromDouble(self->img->airmass);
}

PyObject* PyImgData_get_rx(PyImgDataObject *self, void *closure) {
	return PyLong_FromLong(self->img->rx);
}

PyObject* PyImgData_get_ry(PyImgDataObject *self, void *closure) {
	return PyLong_FromLong(self->img->ry);
}
