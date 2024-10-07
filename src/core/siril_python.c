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
#include <structmember.h>

#include "core/siril.h"
#include "core/command_line_processor.h"
#include "gui/script_menu.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"

extern PyTypeObject PyFitsType;

// Python object to wrap our fits structure - note, this is the Siril-specific fits
// structure and is NOT interchangeable with astropy.io.fits - however a fits file
// can be saved and opened by importing astropy into a script and using it directly
// if astropy processing is desired. The result would have to be saved and then loaded
// back using siril.process_command("load", f"{filename}")

typedef struct {
	PyObject_HEAD
	fits *fit;
	int should_free;  // Flag to indicate if the fits pointer itself should be freed
} PyFits;

// Deallocation function for PyFits
static void PyFits_dealloc(PyFits *self) {
	if (self->fit != NULL) {
		clearfits(self->fit);  // Always clear the internal structures
		if (self->should_free) {
			free(self->fit);  // Only free the pointer if it was dynamically allocated
		}
	}
	Py_TYPE(self)->tp_free((PyObject *)self);
}

// New function for PyFits
static PyObject *PyFits_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	PyFits *self;
	self = (PyFits *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->fit = NULL;
		self->should_free = 0;
	}
	return (PyObject *)self;
}

// Initialize function for PyFits
static int PyFits_init(PyFits *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"width", "height", "nblayer", "type", NULL};
	int width = 0, height = 0, nblayer = 1;
	data_type type = DATA_USHORT;  // Assuming DATA_USHORT is the default

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiii", kwlist,
									&width, &height, &nblayer, &type))
		return -1;

	if (width > 0 && height > 0) {
		// Create a new fits structure
		self->fit = malloc(sizeof(fits));
		if (self->fit == NULL) {
			PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for fits");
			return -1;
		}
		if (new_fit_image(&self->fit, width, height, nblayer, type) != 0) {
			free(self->fit);
			self->fit = NULL;
			PyErr_SetString(PyExc_RuntimeError, "Failed to create new fits image");
			return -1;
		}
		self->should_free = 1;  // This fits was dynamically allocated
	}

	return 0;
}

// Method to get rx
static PyObject *PyFits_get_rx(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, "fit is not initialized");
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->rx);
}

// Method to set rx
static int PyFits_set_rx(PyFits *self, PyObject *value, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, "fit is not initialized");
		return -1;
	}
	if (value == NULL) {
		PyErr_SetString(PyExc_TypeError, "Cannot delete rx");
		return -1;
	}
	if (!PyLong_Check(value)) {
		PyErr_SetString(PyExc_TypeError, "rx must be an integer");
		return -1;
	}
	self->fit->rx = PyLong_AsUnsignedLong(value);
	return 0;
}

// TODO: the rx getter and setter is a proof of concept. Need to decide whether
// such low level manipulation of fits is desirable (in which case add more
// getters and setters) or whether the python interface should remain at a
// higher level.

// Define getter and setter for rx
static PyGetSetDef PyFits_getsetters[] = {
	{"rx", (getter)PyFits_get_rx, (setter)PyFits_set_rx, "image width", NULL},
	{NULL}
};

// Method to access gfit
static PyObject *PyFits_gfit(PyObject *cls, PyObject *args) {
	PyFits *self = (PyFits *)PyFitsType.tp_alloc(&PyFitsType, 0);
	if (self != NULL) {
		self->fit = &gfit;
		self->should_free = 0;  // gfit is statically allocated, don't free it
	}
	return (PyObject *)self;
}

// Define methods for PyFits
static PyMethodDef PyFits_methods[] = {
	{"gfit", (PyCFunction)PyFits_gfit, METH_CLASS | METH_NOARGS, "Get the global fit object"},
	{NULL}
};

// Define the PyFits type
PyTypeObject PyFitsType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.Fits",
	.tp_doc = "Fits object",
	.tp_basicsize = sizeof(PyFits),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_new = PyFits_new,
	.tp_init = (initproc)PyFits_init,
	.tp_dealloc = (destructor)PyFits_dealloc,
	.tp_methods = PyFits_methods,
	.tp_getset = PyFits_getsetters,
};

static PyObject* siril_processcommand(PyObject* self, PyObject* args) {
	Py_ssize_t num_args = PyTuple_Size(args);
	if (num_args < 1) {
		PyErr_SetString(PyExc_TypeError, "At least one argument is required");
		return NULL;
	}

	// Process the first argument (command)
	PyObject* command_obj = PyTuple_GetItem(args, 0);  // Borrowed reference
	if (!command_obj) {
		PyErr_SetString(PyExc_TypeError, "Failed to get the command argument");
		return NULL;
	}

	gchar* command = NULL;
	if (PyUnicode_Check(command_obj)) {
		command = g_strdup(PyUnicode_AsUTF8(command_obj));
	} else if (PyBytes_Check(command_obj)) {
		command = g_strdup(PyBytes_AsString(command_obj));
	} else {
		PyErr_SetString(PyExc_TypeError, "Command must be a string");
		return NULL;
	}

	GString* full_command = g_string_new(command);
	g_free(command);

	// Process additional arguments if any
	for (Py_ssize_t i = 1; i < num_args; i++) {
		PyObject* arg_obj = PyTuple_GetItem(args, i);  // Borrowed reference

		if (arg_obj != NULL) {
			gchar* arg_str = NULL;
			if (PyUnicode_Check(arg_obj)) {
				arg_str = g_strdup(PyUnicode_AsUTF8(arg_obj));
			} else if (PyBytes_Check(arg_obj)) {
				arg_str = g_strdup(PyBytes_AsString(arg_obj));
			} else {
				// Convert non-string objects to string
				PyObject* str_item = PyObject_Str(arg_obj);
				if (str_item != NULL) {
					arg_str = g_strdup(PyUnicode_AsUTF8(str_item));
					Py_DECREF(str_item);
				}
			}

			if (arg_str != NULL) {
				g_string_append_printf(full_command, " %s", arg_str);
				g_free(arg_str);
			} else {
				PyErr_SetString(PyExc_ValueError, "Failed to process argument");
				g_string_free(full_command, TRUE);
				return NULL;
			}
		}
	}

	int result = processcommand(full_command->str);
	g_string_free(full_command, TRUE);

	return PyLong_FromLong(result);
}


// Method to call siril_log_message
static PyObject *siril_log_message_wrapper(PyObject *self, PyObject *args) {
	const char *message;
	if (!PyArg_ParseTuple(args, "s", &message))
		return NULL;

	// Create a new string with a newline character - it feels more python-like not to require a \n in the string
	size_t message_length = strlen(message);
	char *message_with_newline = malloc(message_length + 2); // +1 for '\n' +1 for '\0'

	if (message_with_newline == NULL) {
		return NULL; // Handle memory allocation failure
	}

	// Copy the original message and add a newline
	strcpy(message_with_newline, message);
	message_with_newline[message_length] = '\n'; // Add newline
	message_with_newline[message_length + 1] = '\0'; // Null-terminate the string

	siril_log_message(message_with_newline); // Call the logging function

	free(message_with_newline); // Free the allocated memory
	Py_RETURN_NONE;
}

// Function to return com.wd
static PyObject *siril_get_wd(PyObject *self, PyObject *args) {
	if (com.wd == NULL) {
		Py_RETURN_NONE;  // If com.wd is NULL, return None
	}
	return PyUnicode_FromString(com.wd);
}

// Function to return com.wd
static PyObject *siril_get_filename(PyObject *self, PyObject *args) {
	if (single_image_is_loaded()) {
		if (com.uniq == NULL) {
			return PyUnicode_FromString(_("unsaved_file.fit"));
		}
		if (com.uniq->filename == NULL) {
			return PyUnicode_FromString(_("unsaved_file.fit"));
		}
		return PyUnicode_FromString(com.uniq->filename);
	} else if (sequence_is_loaded() && com.seq.type == SEQ_REGULAR) {
		char filename[256];
		fit_sequence_get_image_filename(&com.seq, com.seq.current, filename, TRUE);
		gchar* path = g_strdup_printf("%s/%s", com.wd, filename);
		PyObject* retval = PyUnicode_FromString(path);
		g_free(path);
		return retval;
	} else {
		// We shouldn't try to handle single-file sequences in this way
		Py_RETURN_NONE;
	}

}

// Define methods for the module
static PyMethodDef SirilMethods[] = {
	{"processcommand", siril_processcommand, METH_VARARGS, "Execute a Siril command"},
	{"log_message", siril_log_message_wrapper, METH_VARARGS, "Log a message"},
	{"filename", (PyCFunction)siril_get_filename, METH_NOARGS, "Get the current image filename"},
	{"wd", (PyCFunction)siril_get_wd, METH_NOARGS, "Get the current working directory"},
	{NULL, NULL, 0, NULL}  /* Sentinel */
};

// Module definition
static struct PyModuleDef sirilmodule = {
	PyModuleDef_HEAD_INIT,
	"siril",   /* name of module */
	NULL, /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
				or -1 if the module keeps state in global variables. */
	SirilMethods
};

// Module initialization function
PyMODINIT_FUNC PyInit_siril(void) {
	PyObject *m;

	if (PyType_Ready(&PyFitsType) < 0)
		return NULL;

	m = PyModule_Create(&sirilmodule);
	if (m == NULL)
		return NULL;

	Py_INCREF(&PyFitsType);
	if (PyModule_AddObject(m, "Fits", (PyObject *) &PyFitsType) < 0) {
		Py_DECREF(&PyFitsType);
		Py_DECREF(m);
		return NULL;
	}

	return m;
}

// Function to initialize Python interpreter and load our module
void init_python(void) {
	PyImport_AppendInittab("siril", PyInit_siril);
	Py_Initialize();
	PyImport_ImportModule("siril");
	PyEval_SaveThread();  // Save the current thread state and release the GIL
	siril_log_message(_("Python scripting module initialized.\n"));
}

// Function to run a Python script from a file
gpointer run_python_script_from_file(gpointer p) {
	const char *script_path = (const char*) p;
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();  // Acquire the GIL
	FILE *fp = fopen(script_path, "r");
	int retval = -1;
	if (fp) {
		retval = PyRun_SimpleFile(fp, script_path);
		fclose(fp);
	} else {
		fprintf(stderr, "Failed to open script file: %s\n", script_path);
	}
	PyGILState_Release(gstate);  // Release the GIL
	g_idle_add(script_widgets_idle, NULL);
	return GINT_TO_POINTER(retval);
}

// Function to run a Python script from memory
gpointer run_python_script_from_mem(gpointer p) {
	const char *script_contents = (const char*) p;
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();  // Acquire the GIL

	int retval = PyRun_SimpleString(script_contents);
	PyGILState_Release(gstate);  // Release the GIL
	// use case for this not fully defined yet (python scribble pad?)
	// but probably doesn't need the script_widgets_idle() call
	return GINT_TO_POINTER(retval);
}

// Function to finalize Python interpreter
void finalize_python(void) {
	Py_Finalize();
	siril_log_message(_("Python scripting module cleaned up.\n"));
}
