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
PyObject* PySeq_get_imgdata(PySeqObject *self, PyObject *args) {
	if (self->seq == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("sequence is not initialized"));
		return NULL;
	}
	int frame;
	if (!PyArg_ParseTuple(args, "i", &frame))
		return NULL;

	if (frame < 0 || frame >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, _("Frame index out of range"));
		return NULL;
	}
	if (self->seq->imgparam == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("imgdata is not allocated"));
		return NULL;
	}

	PyImgDataObject *imgdata_obj = (PyImgDataObject*)PyObject_New(PyImgDataObject, &PyImgDataType);
	if (!imgdata_obj)
		return NULL;

	imgdata_obj->img = &(self->seq->imgparam[frame]);
	return (PyObject*)imgdata_obj;
}

PyObject* PySeq_get_regdata(PySeqObject *self, PyObject *args) {
	int frame, layer;
	if (!PyArg_ParseTuple(args, "ii", &frame, &layer))
		return NULL;

	if (frame < 0 || frame >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, "Frame index out of range");
		return NULL;
	}

	if (layer < 0 || layer >= self->seq->nb_layers) {
		PyErr_SetString(PyExc_IndexError, "Channel index out of range");
		return NULL;
	}

	PyRegDataObject *regdata_obj = (PyRegDataObject*)PyObject_New(PyRegDataObject, &PyRegDataType);
	if (!regdata_obj)
		return NULL;

	regdata_obj->reg = &self->seq->regparam[layer][frame];
	return (PyObject*)regdata_obj;
}

PyObject* PySeq_get_homography(PySeqObject *self, PyObject *args) {
	int frame, layer;
	if (!PyArg_ParseTuple(args, "ii", &frame, &layer))
		return NULL;

	if (frame < 0 || frame >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, "Frame index out of range");
		return NULL;
	}

	PyHomographyObject *homography_obj = (PyHomographyObject*)PyObject_New(PyHomographyObject, &PyHomographyType);
	if (!homography_obj)
		return NULL;

	homography_obj->homography = &(self->seq->regparam[layer][frame].H);
	return (PyObject*)homography_obj;
}

PyObject* PySeq_get_imstats(PySeqObject *self, PyObject *args) {
	int frame, channel;
	if (!PyArg_ParseTuple(args, "ii", &frame, &channel))
		return NULL;

	if (frame < 0 || frame >= self->seq->number) {
		PyErr_SetString(PyExc_IndexError, "Frame index out of range");
		return NULL;
	}

	if (channel < 0 || channel >= self->seq->nb_layers) {
		PyErr_SetString(PyExc_IndexError, "Channel index out of range");
		return NULL;
	}

	PyImStatsObject *imstats_obj = (PyImStatsObject*)PyObject_New(PyImStatsObject, &PyImStatsType);
	if (!imstats_obj)
		return NULL;

	imstats_obj->stats = self->seq->stats[frame][channel];
	return (PyObject*)imstats_obj;
}
