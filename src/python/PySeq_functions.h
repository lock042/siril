#ifndef PYSEQ_FUNCTIONS_H
#define PYSEQ_FUNCTIONS_H

#include "python/siril_python.h"

PyObject *PySeq_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int PySeq_init(PySeqObject *self, PyObject *args, PyObject *kwds);
void PySeq_dealloc(PySeqObject *self);
PyObject *PySeq_FromExisting(struct sequ *seq);
int PySeq_traverse(PySeqObject *self, visitproc visit, void *arg);
int PySeq_clear(PySeqObject *self);

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
PyObject* PySeq_get_imgdata(PySeqObject *self, PyObject *args);
PyObject* PySeq_get_regdata(PySeqObject *self, PyObject *args);
PyObject* PySeq_get_homography(PySeqObject *self, PyObject *args);
PyObject* PySeq_get_imstats(PySeqObject *self, PyObject *args);
PyObject *PySeq_comseq(PyObject *cls, PyObject *args);
#endif