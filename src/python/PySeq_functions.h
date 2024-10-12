#ifndef PYFITS_FUNCTIONS_H
#define PYFITS_FUNCTIONS_H

#include "python/siril_python.h"

PyObject* PySeq_get_seqname(PySeqObject *self, void *closure);
PyObject* PySeq_get_number(PySeqObject *self, void *closure);
PyObject* PySeq_get_selnum(PySeqObject *self, void *closure);
PyObject* PySeq_get_fixed(PySeqObject *self, void *closure);
PyObject* PySeq_get_nb_layers(PySeqObject *self, void *closure);
PyObject* PySeq_get_rx(PySeqObject *self, void *closure);
PyObject* PySeq_get_ry(PySeqObject *self, void *closure);
PyObject* PySeq_get_is_variable(PySeqObject *self, void *closure);
PyObject* PySeq_get_bitpix(PySeqObject *self, void *closure);
PyObject* PySeq_get_reference_image(PySeqObject *self, void *closure);
PyObject* PySeq_get_type(PySeqObject *self, void *closure);
PyObject* PySeq_get_current(PySeqObject *self, void *closure);
PyObject* PySeq_get_needs_saving(PySeqObject *self, void *closure);

#endif
