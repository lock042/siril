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

#include "io/sequence.h"

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
		free_sequence(self->seq, TRUE);
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
	// Parse the frame argument
	if (!PyArg_ParseTuple(args, "i", &frame)) {
		return NULL; // Return NULL if the frame is not passed correctly
	}
	// Bounds checking: Ensure 0 < frame <= self->seq->number
	if (frame < 0 || frame >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, _("Frame index out of range"));
		return NULL;
	}
	// Create a new PyImgDataObject
	PyImgDataObject *py_imgdata = (PyImgDataObject *)PyObject_New(PyImgDataObject, &PyImgDataType);
	if (py_imgdata == NULL) {
		return PyErr_NoMemory(); // Return memory error if allocation fails
	}
	if (self->seq->imgparam == NULL) {
		PyErr_SetString(PyExc_IndexError, _("imgdata unallocated"));
		return NULL;
	}
	// Copy the imgparam data from the selected frame
	py_imgdata->img = (imgdata *)calloc(1, sizeof(imgdata));
	if (py_imgdata->img == NULL) {
		Py_DECREF(py_imgdata);
		return PyErr_NoMemory(); // Return memory error if the allocation fails
	}
	// Perform the copy from self->seq->imgparam[frame] to the new PyImgDataObject
	memcpy(py_imgdata->img, &(self->seq->imgparam[frame]), sizeof(imgdata));
	// Set additional fields (if necessary)
	py_imgdata->should_free = 1;  // The object owns the img copy and should free it later
	py_imgdata->seq = NULL;    // No parent is set for this object
	// Return the new PyImgDataObject
	return (PyObject *)py_imgdata;
}

PyObject *PySeq_get_regdata(PySeqObject *self, PyObject *args) {
	int frame, channel;
	// Parse the frame and channel arguments
	if (!PyArg_ParseTuple(args, "ii", &frame, &channel)) {
		return NULL; // Return NULL if frame and channel are not passed correctly
	}
	// Bounds checking: Ensure 0 <= frame < self->seq->number and 0 <= channel < self->seq->nb_layers
	if (frame < 0 || frame >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, _("Frame index out of range"));
		return NULL;
	}
	if (channel < 0 || channel >= self->seq->nb_layers) {
		PyErr_SetString(PyExc_IndexError, _("Channel index out of range"));
		return NULL;
	}
	if (self->seq->regparam == NULL || self->seq->regparam[frame] == NULL) {
		PyErr_SetString(PyExc_IndexError, _("regdata unallocated"));
		return NULL;
	}
	// Create a new PyRegDataObject
	PyRegDataObject *py_regdata = (PyRegDataObject *)PyObject_New(PyRegDataObject, &PyRegDataType);
	if (py_regdata == NULL) {
		return PyErr_NoMemory(); // Return memory error if allocation fails
	}
	// Copy the regparam data from the selected frame and channel
	py_regdata->reg = (regdata *)calloc(1, sizeof(regdata));
	if (py_regdata->reg == NULL) {
		Py_DECREF(py_regdata);
		return PyErr_NoMemory(); // Return memory error if the allocation fails
	}
	// Perform the copy from self->seq->regparam[frame][channel] to the new PyRegDataObject
	memcpy(py_regdata->reg, &(self->seq->regparam[frame][channel]), sizeof(regdata));
	// don't copy fwhm_data, it's only used during registration and not saved

	// Set additional fields (if necessary)
	py_regdata->should_free = 1;  // The object owns the reg copy and should free it later
	py_regdata->seq = NULL;    // No parent is set for this object
	py_regdata->ref_count = 1; // initialize the refcount at 1
	// Return the new PyRegDataObject
	return (PyObject *)py_regdata;
}

PyObject *PySeq_get_imstats(PySeqObject *self, PyObject *args) {
	int frame, channel;
	// Parse the frame and channel arguments
	if (!PyArg_ParseTuple(args, "ii", &frame, &channel)) {
		return NULL; // Return NULL if frame and channel are not passed correctly
	}
	// Bounds checking: Ensure 0 <= frame < self->seq->number and 0 <= channel < self->seq->nb_layers
	if (frame < 0 || frame >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, _("Frame index out of range"));
		return NULL;
	}
	if (channel < 0 || channel >= self->seq->nb_layers) {
		PyErr_SetString(PyExc_IndexError, _("Channel index out of range"));
		return NULL;
	}
	// Check if stats exist for the given frame and channel
	if (self->seq->stats == NULL || self->seq->stats[frame] == NULL || self->seq->stats[frame][channel] == NULL) {
		PyErr_SetString(PyExc_ValueError, _("No stats available for the requested frame and channel"));
		return NULL;
	}
	// Create a new PyImStatsObject
	PyImStatsObject *py_imstats = (PyImStatsObject *)PyObject_New(PyImStatsObject, &PyImStatsType);
	if (py_imstats == NULL) {
		return PyErr_NoMemory(); // Return memory error if allocation fails
	}
	// Copy the imstats from the selected frame and channel
	py_imstats->stats = (imstats *)malloc(sizeof(imstats));
	if (py_imstats->stats == NULL) {
		Py_DECREF(py_imstats);
		return PyErr_NoMemory(); // Return memory error if the allocation fails
	}
	// Perform the copy from self->seq->stats[frame][channel] to the new PyImStatsObject
	memcpy(py_imstats->stats, self->seq->stats[frame][channel], sizeof(imstats));
	// Set additional fields (if necessary)
	py_imstats->should_free = 1;  // The object owns the stats copy and should free it later
	py_imstats->parent = NULL;    // No parent is set for this object
	py_imstats->parent_type = '\0';  // No parent type, so set to null character
	// Return the new PyImStatsObject
	return (PyObject *)py_imstats;
}

PyObject *PySeq_get_fits(PySeqObject *self, PyObject *args) {
	int n;
	// Parse the 'n' argument (the frame index)
	if (!PyArg_ParseTuple(args, "i", &n)) {
		return NULL; // Return NULL if 'n' is not passed correctly
	}
	// Bounds checking: Ensure 0 <= n < self->seq->nb_frames
	if (n < 0 || n >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, "Frame index out of range");
		return NULL;
	}
	// Create a new PyFits object
	PyFits *py_fits = (PyFits *)PyObject_New(PyFits, &PyFitsType);
	if (py_fits == NULL) {
		return PyErr_NoMemory(); // Return memory error if allocation fails
	}
	// Set additional fields for the PyFits object
	py_fits->should_free_data = 1; // The object owns the fits data and should free it later
	py_fits->should_free = 1;      // The PyFits object is responsible for freeing the fits structure
	// Allocate memory for the fits structure
	py_fits->fit = (fits *)calloc(1, sizeof(fits));
	if (py_fits->fit == NULL) {
		Py_DECREF(py_fits); // Free the allocated PyFits object
		return PyErr_NoMemory(); // Return memory error if allocation fails
	}
	// Call seq_read_frame() to read the frame data
	int result = seq_read_frame(self->seq, n, py_fits->fit, !com.pref.force_16bit, -1);
	if (result != 0) {
		Py_DECREF(py_fits);  // Free the allocated PyFits object
		PyErr_SetString(PyExc_RuntimeError, "Failed to read frame data");
		return NULL;
	}
	// Return the newly created PyFits object
	return (PyObject *)py_fits;
}

// Method to access com.seq
PyObject *PySeq_comseq(PyObject *cls, PyObject *args) {
	if (!sequence_is_loaded()) {
		PyErr_SetString(PyExc_ValueError, _("Main sequence is not loaded"));
		return NULL;
	}

	PySeqObject *self = (PySeqObject *)PySeqType.tp_alloc(&PySeqType, 0);
	if (self != NULL) {
		self->seq = &com.seq;
		self->should_free = 0;  // com.seq is statically allocated, don't free it
	}
	return (PyObject *)self;
}
