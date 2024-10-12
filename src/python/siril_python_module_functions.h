#ifndef SIRIL_PYTHON_MODULE_H
#define SIRIL_PYTHON_MODULE_H

PyObject* py_gui_block(PyObject* self, PyObject* args);
PyObject* py_gui_unblock(PyObject* self, PyObject* args);
PyObject* siril_get_continue(PyObject* self, PyObject* args);
PyObject *siril_get_config_item(PyObject *self, PyObject *args);
PyObject* siril_get_filename(PyObject *self, PyObject *args);
PyObject* siril_log_message_wrapper(PyObject *self, PyObject *args);
PyObject* siril_notify_gfit_modified(PyObject* self, PyObject* args);
PyObject* siril_processcommand(PyObject* self, PyObject* args);
PyObject* siril_progress_wrapper(PyObject* self, PyObject* args);
PyObject* siril_get_wd(PyObject *self, PyObject *args);

#endif
