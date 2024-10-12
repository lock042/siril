#ifndef PYIMSTATS_FUNCTIONS_H
#define PYIMSTATS_FUNCTIONS_H

#include "python/siril_python.h"

void PyImStats_dealloc(PyImStatsObject *self);
PyObject *PyImStats_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int PyImStats_init(PyImStatsObject *self, PyObject *args, PyObject *kwds);
PyObject *PyImStats_FromExisting(imstats *stats, PyObject *parent, char parent_type);
int PyImStats_traverse(PyImStatsObject *self, visitproc visit, void *arg);
int PyImStats_clear(PyImStatsObject *self);

PyObject* PyImStats_get_total(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_ngoodpix(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_mean(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_median(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_sigma(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_avgDev(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_mad(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_sqrtbwmv(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_location(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_scale(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_min(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_max(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_normValue(PyImStatsObject *self, void *closure);
PyObject* PyImStats_get_bgnoise(PyImStatsObject *self, void *closure);

#endif
