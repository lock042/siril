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
	return PyUnicode_FromString(self->seq->seqname);
}

PyObject* PySeq_get_number(PySeqObject *self, void *closure) {
	return PyLong_FromLong(self->seq->number);
}

PyObject* PySeq_get_selnum(PySeqObject *self, void *closure) {
	return PyLong_FromLong(self->seq->selnum);
}

PyObject* PySeq_get_fixed(PySeqObject *self, void *closure) {
	return PyLong_FromLong(self->seq->fixed);
}

PyObject* PySeq_get_nb_layers(PySeqObject *self, void *closure) {
	return PyLong_FromLong(self->seq->nb_layers);
}

PyObject* PySeq_get_rx(PySeqObject *self, void *closure) {
	return PyLong_FromUnsignedLong(self->seq->rx);
}

PyObject* PySeq_get_ry(PySeqObject *self, void *closure) {
	return PyLong_FromUnsignedLong(self->seq->ry);
}

PyObject* PySeq_get_is_variable(PySeqObject *self, void *closure) {
	return PyBool_FromLong(self->seq->is_variable);
}

PyObject* PySeq_get_bitpix(PySeqObject *self, void *closure) {
	return PyLong_FromLong(self->seq->bitpix);
}

PyObject* PySeq_get_reference_image(PySeqObject *self, void *closure) {
	return PyLong_FromLong(self->seq->reference_image);
}

PyObject* PySeq_get_type(PySeqObject *self, void *closure) {
	return PyLong_FromLong(self->seq->type);
}

PyObject* PySeq_get_current(PySeqObject *self, void *closure) {
	return PyLong_FromLong(self->seq->current);
}

PyObject* PySeq_get_needs_saving(PySeqObject *self, void *closure) {
	return PyBool_FromLong(self->seq->needs_saving);
}
