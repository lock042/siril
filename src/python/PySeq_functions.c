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
#include "python/PyImgData_functions.h"
#include "python/PyImStats_functions.h"
#include "python/PyRegData_functions.h"

PyObject *PySeq_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	PySeqObject *self;
	self = (PySeqObject *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->seq = NULL;
		self->should_free = 0;
		self->ref_count = 1;  // Initialize ref_count to 1
	}
	return (PyObject *)self;
}

int PySeq_init(PySeqObject *self, PyObject *args, PyObject *kwds) {
	self->seq = (struct sequ *)calloc(1, sizeof(struct sequ));
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_MemoryError, N_("Failed to allocate memory for sequence"));
		return -1;
	}
	self->should_free = 1;
	self->ref_count = 1;  // Initialize ref_count to 1
	return 0;
}

void PySeq_dealloc(PySeqObject *self) {
	if (self->seq != NULL && self->should_free) {
		free(self->seq->seqname);
		free(self->seq->imgparam);
		if (self->seq->regparam) {
			for (int i = 0; i < self->seq->nb_layers; i++) {
				free(self->seq->regparam[i]);
			}
			free(self->seq->regparam);
		}
		if (self->seq->stats) {
			for (int i = 0; i < self->seq->number; i++) {
				if (self->seq->stats[i]) {
					for (int j = 0; j < self->seq->nb_layers; j++) {
						free(self->seq->stats[i][j]);
					}
					free(self->seq->stats[i]);
				}
			}
			free(self->seq->stats);
		}
		// Free other dynamically allocated members here
		free(self->seq);
	}
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *PySeq_FromExisting(struct sequ *seq) {
	PySeqObject *obj = (PySeqObject *)PySeq_new(&PySeqType, NULL, NULL);
	if (obj != NULL) {
		obj->seq = seq;
		obj->should_free = 0;
		obj->ref_count = 1;  // Initialize ref_count to 1
	}
	return (PyObject *)obj;
}

int PySeq_traverse(PySeqObject *self, visitproc visit, void *arg) {
	// We need to visit all PyObject members that could be part of a reference cycle
	// In this case, we don't have any such objects in PySequenceObject itself
	return 0;
}

int PySeq_clear(PySeqObject *self) {
	// Clear all references that could be part of a reference cycle
	// Again, in this case, we don't have any such references
	return 0;
}

// Getter functions for PySeqType
PyObject* PySeq_get_seqname(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	if (self->seq->seqname == NULL) {
		Py_RETURN_NONE;
	}
	return PyUnicode_FromString(self->seq->seqname);
}

PyObject* PySeq_get_number(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->seq->number);
}

PyObject* PySeq_get_selnum(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->seq->selnum);
}

PyObject* PySeq_get_fixed(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->seq->fixed);
}

PyObject* PySeq_get_nb_layers(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->seq->nb_layers);
}

PyObject* PySeq_get_rx(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->seq->rx);
}

PyObject* PySeq_get_ry(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->seq->ry);
}

PyObject* PySeq_get_is_variable(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyBool_FromLong(self->seq->is_variable);
}

PyObject* PySeq_get_bitpix(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->seq->bitpix);
}

PyObject* PySeq_get_reference_image(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->seq->reference_image);
}

PyObject* PySeq_get_type(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->seq->type);
}

PyObject* PySeq_get_current(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->seq->current);
}

PyObject* PySeq_get_needs_saving(PySeqObject *self, void *closure) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	return PyBool_FromLong(self->seq->needs_saving);
}

// Method to get imgdata for a specified frame in the sequence
PyObject *PySeq_get_imgdata(PySeqObject *self, PyObject *args) {
	int frame;
	if (!PyArg_ParseTuple(args, "i", &frame))
		return NULL;

	if (frame < 0 || frame >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, "Frame index out of range");
		return NULL;
	}

	if (!self->seq->imgparam) {
		Py_RETURN_NONE;
	}

	PyObject *imgdata = PyImgData_FromExisting(&self->seq->imgparam[frame], (PyObject *)self);
	if (imgdata != NULL) {
		self->ref_count++;
	}
	return imgdata;
}

PyObject* PySeq_get_regdata(PySeqObject *self, PyObject *args) {
	int frame, channel;
	if (!PyArg_ParseTuple(args, "ii", &frame, &channel))
		return NULL;

	if (frame < 0 || frame >= self->seq->number || channel < 0 || channel >= self->seq->nb_layers) {
		PyErr_SetString(PyExc_IndexError, "Frame or channel index out of range");
		return NULL;
	}

	if (!self->seq->regparam || !self->seq->regparam[channel]) {
		Py_RETURN_NONE;
	}

	PyObject *regdata = PyRegData_FromExisting(&self->seq->regparam[channel][frame], (PyObject *)self);
	if (regdata != NULL) {
		self->ref_count++;
	}
	return regdata;
}

PyObject *PySeq_get_imstats(PySeqObject *self, PyObject *args) {
	int frame, channel;
	if (!PyArg_ParseTuple(args, "ii", &frame, &channel))
		return NULL;

	if (frame < 0 || frame >= self->seq->number || channel < 0 || channel >= self->seq->nb_layers) {
		PyErr_SetString(PyExc_IndexError, N_("Frame or channel index out of range"));
		return NULL;
	}

	if (!self->seq->stats || !self->seq->stats[frame] || !self->seq->stats[frame][channel]) {
		Py_RETURN_NONE;
	}

	PyObject *stats = PyImStats_FromExisting(self->seq->stats[frame][channel], (PyObject *)self, 'S');
	if (stats != NULL) {
		self->ref_count++;
	}
	return stats;
}

// Method to access gfit
PyObject *PySeq_comseq(PyObject *cls, PyObject *args) {
	PyFits *self = (PySeq *)PySeqType.tp_alloc(&PySeqType, 0);
	if (self != NULL) {
		self->seq = &com.seq;
		self->should_free = 0;  // com.seq is statically allocated, don't free it
	}
	return (PyObject *)self;
}
