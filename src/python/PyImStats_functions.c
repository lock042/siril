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

PyObject *PyImStats_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    PyImStatsObject *self;
    self = (PyImStatsObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->stats = NULL;
        self->should_free = 0;
        self->parent = NULL;
        self->parent_type = '\0';
    }
    return (PyObject *)self;
}

int PyImStats_init(PyImStatsObject *self, PyObject *args, PyObject *kwds) {
	self->stats = (imstats *)calloc(1, sizeof(imstats));
	if (self->stats == NULL) {
		PyErr_SetString(PyExc_MemoryError, N_("Failed to allocate memory for imstats"));
		return -1;
	}
	self->should_free = 1;
	return 0;
}

void PyImStats_dealloc(PyImStatsObject *self) {
    if (self->parent != NULL) {
        if (self->parent_type == 'S') {
            PySeqObject *seq = (PySeqObject *)self->parent;
            seq->ref_count--;
            if (seq->ref_count == 0) {
                Py_DECREF(self->parent);
            }
        } else if (self->parent_type == 'F') {
            // Assuming PyFitsObject has a similar ref_count mechanism
            PyFits *fits = (PyFits *)self->parent;
            fits->ref_count--;
            if (fits->ref_count == 0) {
                Py_DECREF(self->parent);
            }
        }
    }
    if (self->stats != NULL && self->should_free) {
        free(self->stats);
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *PyImStats_FromExisting(imstats *stats, PyObject *parent, char parent_type) {
    PyImStatsObject *obj = (PyImStatsObject *)PyImStats_new(&PyImStatsType, NULL, NULL);
    if (obj != NULL) {
        obj->stats = stats;
        obj->should_free = 0;
        Py_INCREF(parent);
        obj->parent = parent;
        obj->parent_type = parent_type;
    }
    return (PyObject *)obj;
}

int PyImStats_traverse(PyImStatsObject *self, visitproc visit, void *arg) {
    Py_VISIT(self->parent);
    return 0;
}

int PyImStats_clear(PyImStatsObject *self) {
    Py_CLEAR(self->parent);
    return 0;
}

// PyImStatsType methods and getters
PyObject* PyImStats_get_total(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyLong_FromLong(self->stats->total);
}

PyObject* PyImStats_get_ngoodpix(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyLong_FromLong(self->stats->ngoodpix);
}

PyObject* PyImStats_get_mean(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->mean);
}

PyObject* PyImStats_get_median(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->median);
}

PyObject* PyImStats_get_sigma(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->sigma);
}

PyObject* PyImStats_get_avgDev(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->avgDev);
}

PyObject* PyImStats_get_mad(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->mad);
}

PyObject* PyImStats_get_sqrtbwmv(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->sqrtbwmv);
}

PyObject* PyImStats_get_location(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->location);
}

PyObject* PyImStats_get_scale(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->scale);
}

PyObject* PyImStats_get_min(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->min);
}

PyObject* PyImStats_get_max(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->max);
}

PyObject* PyImStats_get_normValue(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->normValue);
}

PyObject* PyImStats_get_bgnoise(PyImStatsObject *self, void *closure) {
	if (self->stats == NULL) Py_RETURN_NONE;
	return PyFloat_FromDouble(self->stats->bgnoise);
}
