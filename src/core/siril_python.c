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
			PyErr_SetString(PyExc_MemoryError, _("Failed to allocate memory for fits"));
			return -1;
		}
		if (new_fit_image(&self->fit, width, height, nblayer, type) != 0) {
			free(self->fit);
			self->fit = NULL;
			PyErr_SetString(PyExc_RuntimeError, _("Failed to create new fits image"));
			return -1;
		}
		self->should_free = 1;  // This fits was dynamically allocated
	}

	return 0;
}

// Typed memoryview getter (returns a memoryview with the correct type, wrapping existing planar data)
static PyObject* PyFits_get_pixel_data(PyFits *self, PyObject *Py_UNUSED(ignored)) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}

	void *data_ptr;
	Py_ssize_t buffer_size;
	const char *format;
	int itemsize;

	// Determine data type and size based on fit->type
	if (self->fit->type == DATA_FLOAT) {
		data_ptr = self->fit->fdata;
		itemsize = sizeof(float);
		format = "f";  // 32-bit float format
	} else if (self->fit->type == DATA_USHORT) {
		data_ptr = self->fit->data;
		itemsize = sizeof(WORD);
		format = "H";  // 16-bit unsigned int format
	} else {
		PyErr_SetString(PyExc_TypeError, _("Unsupported FITS data type"));
		return NULL;
	}

	int ndim = self->fit->naxis;
	Py_ssize_t width = self->fit->naxes[0];
	Py_ssize_t height = self->fit->naxes[1];
	Py_ssize_t depth = (ndim > 2) ? self->fit->naxes[2] : 1;

	buffer_size = width * height * depth * itemsize;

	// Allocate memory for shape and strides
	Py_ssize_t *shape = PyMem_Malloc(ndim * sizeof(Py_ssize_t));
	Py_ssize_t *strides = PyMem_Malloc(ndim * sizeof(Py_ssize_t));
	if (shape == NULL || strides == NULL) {
		PyMem_Free(shape);
		PyMem_Free(strides);
		PyErr_NoMemory();
		return NULL;
	}

	// Set shape and strides for planar data
	shape[0] = depth;
	shape[1] = height;
	shape[2] = width;
	strides[0] = height * width * itemsize;  // stride between planes
	strides[1] = width * itemsize;           // stride between rows
	strides[2] = itemsize;                   // stride between elements in a row

	// Create a new Py_buffer and fill it with the correct information
	Py_buffer view;
	if (PyBuffer_FillInfo(&view, (PyObject*)self, data_ptr, buffer_size, 0, PyBUF_CONTIG) == -1) {
		PyMem_Free(shape);
		PyMem_Free(strides);
		return NULL;
	}

	// Manually set the format, itemsize, shape, and strides
	view.format = PyMem_Malloc(strlen(format) + 1);
	if (view.format == NULL) {
		PyMem_Free(shape);
		PyMem_Free(strides);
		PyErr_NoMemory();
		return NULL;
	}
	strcpy((char *)view.format, format);

	view.itemsize = itemsize;
	view.ndim = ndim;
	view.shape = shape;
	view.strides = strides;

	// Create and return a memoryview wrapping the existing data
	PyObject *memoryview = PyMemoryView_FromBuffer(&view);
	if (memoryview == NULL) {
		PyMem_Free((void *)view.format);
		PyMem_Free(shape);
		PyMem_Free(strides);
	}
	return memoryview;
}

static int PyFits_set_pixel_data_typed_with_size(PyFits *self, PyObject *value, int rx, int ry, int nchans) {
	// Ensure the input is a memoryview object
	if (!PyMemoryView_Check(value)) {
		PyErr_SetString(PyExc_TypeError, _("Expected a memoryview object"));
		return -1;
	}

	Py_buffer view;
	if (PyObject_GetBuffer(value, &view, PyBUF_CONTIG_RO) != 0) {
		return -1;
	}

	// Validate the input dimensions: rx * ry * nchans
	Py_ssize_t expected_size = rx * ry * nchans;
	if (view.len / view.itemsize != expected_size) {
		PyErr_SetString(PyExc_ValueError, _("Input data size does not match the specified dimensions (rx * ry * nchans)"));
		PyBuffer_Release(&view);
		return -1;
	}

	// Free existing data if it is not NULL
	if (self->fit->data != NULL) {
		free(self->fit->data);
		self->fit->data = NULL;
	}
	if (self->fit->fdata != NULL) {
		free(self->fit->fdata);
		self->fit->fdata = NULL;
	}

	// Set the correct type and allocate new memory based on the buffer's item size
	if (view.itemsize == sizeof(float)) {
		// 32-bit float data
		self->fit->type = DATA_FLOAT;
		self->fit->fdata = (float *)malloc(view.len);
		if (self->fit->fdata == NULL) {
			PyErr_SetString(PyExc_MemoryError, _("Failed to allocate memory for FITS float data"));
			PyBuffer_Release(&view);
			return -1;
		}
		// Copy the data from the memoryview into fit->fdata
		memcpy(self->fit->fdata, view.buf, view.len);
	} else if (view.itemsize == sizeof(WORD)) {
		// 16-bit unsigned integer (WORD) data
		self->fit->type = DATA_USHORT;
		self->fit->data = (WORD *)malloc(view.len);
		if (self->fit->data == NULL) {
			PyErr_SetString(PyExc_MemoryError, _("Failed to allocate memory for FITS WORD data"));
			PyBuffer_Release(&view);
			return -1;
		}
		// Copy the data from the memoryview into fit->data
		memcpy(self->fit->data, view.buf, view.len);
	} else {
		PyErr_SetString(PyExc_TypeError, _("Unsupported memoryview type: expected 16-bit or 32-bit"));
		PyBuffer_Release(&view);
		return -1;
	}

	// Update FITS dimensions
	self->fit->rx = rx;
	self->fit->ry = ry;
	self->fit->naxes[0] = rx;
	self->fit->naxes[1] = ry;
	self->fit->naxes[2] = nchans;
	self->fit->naxis = nchans == 1 ? 2 : 3;

	// Release the buffer
	PyBuffer_Release(&view);
	return 0;
}

static int PyFits_set_pixel_data(PyFits *self, PyObject *args, PyObject *kwds) {
	PyObject *value;
	int rx = 0, ry = 0, nchans = 0;
	static char *kwlist[] = {"buffer", "rx", "ry", "nchans", NULL};

	// Parse the arguments: memoryview buffer, rx, ry, nchans
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|iii", kwlist, &value, &rx, &ry, &nchans)) {
		return -1;
	}

	// If rx, ry, and nchans are not provided, use the current values
	if (rx == 0 || ry == 0 || nchans == 0) {
		rx = self->fit->rx;
		ry = self->fit->ry;
		nchans = self->fit->naxes[2];
	}

	// Call the generalized setter with the determined parameters
	return PyFits_set_pixel_data_typed_with_size(self, value, rx, ry, nchans);
}

// Method to get rx
static PyObject *PyFits_get_rx(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->rx);
}

// Method to get ry
static PyObject *PyFits_get_ry(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->ry);
}

// Method to get naxes[2] (number of channels)
static PyObject *PyFits_get_nchans(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->naxes[2]);
}

// Method to get mini (min pixel value across all channels)
static PyObject *PyFits_get_mini(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->mini);  // Changed to return a float
}

// Method to get negratio
static PyObject *PyFits_get_neg_ratio(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->neg_ratio);  // Changed to return a float
}

// Method to get maxi (max pixel value across all channels)
static PyObject *PyFits_get_maxi(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->maxi);
}

// Method to get top_down (boolean indicating top-down orientation)
static PyObject *PyFits_get_top_down(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	// Return Python boolean: 1 for TRUE, 0 for FALSE
	return PyBool_FromLong(self->fit->top_down ? 1 : 0);
}

// Method to get row_order (indicates row order of the FITS data)
static PyObject *PyFits_get_row_order(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}

	// Assuming fits->keywords.row_order is a string (char*)
	const char *row_order = self->fit->keywords.row_order;

	// Convert TOP_DOWN and BOTTOM_UP to more Pythonic format
	if (strcmp(row_order, "TOP-DOWN") == 0) {
		return PyUnicode_FromString("top_down");
	} else if (strcmp(row_order, "BOTTOM-UP") == 0) {
		return PyUnicode_FromString("bottom_up");
	} else {
		// If row_order contains some other value, return it as is
		return PyUnicode_FromString(row_order);
	}
}

// Method to get bit depth
static PyObject *PyFits_get_bitdepth(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->type == DATA_FLOAT ? 32 : 16);
}

static int PyFits_set_bitdepth(PyFits *self, PyObject *value, void *closure) {
	int layers;

	// Ensure the PyFits instance is not NULL
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return -1; // Indicate failure
	}

	// Parse the input to get the layers (bit depth)
	if (!PyLong_Check(value)) {
		PyErr_SetString(PyExc_TypeError, _("Bit depth must be an integer"));
		return -1; // Indicate failure
	}

	layers = PyLong_AsLong(value);

	// Call the internal function to change the bit depth
	int result = fits_change_depth(self->fit, layers);

	// Check the result of the function
	if (result != 0) {
		PyErr_SetString(PyExc_RuntimeError, _("Failed to change bit depth"));
		return -1; // Indicate failure
	}

	return 0; // Indicate success
}

// Define getters and setters
static PyGetSetDef PyFits_getsetters[] = {
	// rx, ry, nchans have getters only: the fits size should not be changed using these properties
	{"rx", (getter)PyFits_get_rx, NULL, N_("image width"), NULL},
	{"ry", (getter)PyFits_get_ry, NULL, N_("image height"), NULL},
	{"nchans", (getter)PyFits_get_nchans, NULL, N_("image channel depth"), NULL},
	// these don't have setters either, as they are computed values
	{"mini", (getter)PyFits_get_mini, NULL, N_("min value across all channels"), NULL},
	{"maxi", (getter)PyFits_get_maxi, NULL, N_("max value across all channels"), NULL},
	{"neg_ratio", (getter)PyFits_get_neg_ratio, NULL, N_("ratio of negative to nonnegative pixels"), NULL},
	{"top_down", (getter)PyFits_get_top_down, NULL, N_("native roworder of the original file format"), NULL},
	{"roworder", (getter)PyFits_get_row_order, NULL, N_("row order of the sensor that produced the image"), NULL},
	// these properties have both getters and setters
	{"data", (getter)PyFits_get_pixel_data, (setter)PyFits_set_pixel_data, N_("Pixel data buffer"), NULL},
	{"bitdepth", (getter)PyFits_get_bitdepth,(setter)PyFits_set_bitdepth, N_("image bitdepth"), NULL},
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
	{"gfit", (PyCFunction)PyFits_gfit, METH_CLASS | METH_NOARGS, N_("Get the global fits object")},
	{NULL}
};

// Define the PyFits type
PyTypeObject PyFitsType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.Fits",
	.tp_doc = N_("Siril fits object"),
	.tp_basicsize = sizeof(PyFits),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_new = PyFits_new,
	.tp_init = (initproc)PyFits_init,
	.tp_dealloc = (destructor)PyFits_dealloc,
	.tp_methods = PyFits_methods,
	.tp_getset = PyFits_getsetters,
};

/*************************************************
* Functions providing control over the Siril UI *
************************************************/

static PyObject* py_gui_block(PyObject* self, PyObject* args) {
	if (!g_main_context_iteration(NULL, FALSE)) {
		siril_log_color_message(_("Warning: siril.gui_block() must not be called except from a script's GTK main loop.\n"), "red");
		Py_RETURN_NONE;
	}
	script_widgets_enable(FALSE);  // Disable main control window GUI elements
	Py_RETURN_NONE;
}

static PyObject* py_gui_unblock(PyObject* self, PyObject* args) {
	if (!g_main_context_iteration(NULL, FALSE)) {
		siril_log_color_message(_("Warning: siril.gui_unblock() must not be called except from a script's GTK main loop.\n"), "red");
		Py_RETURN_NONE;
	}
	script_widgets_enable(TRUE);   // Enable main control window GUI elements
	Py_RETURN_NONE;
}

static PyObject* PyNotifyGfitModified(PyObject* self, PyObject* args) {
    // Call the C function
    notify_gfit_modified();
    Py_RETURN_NONE;
}

/*************************************************************
* Functions providing control over Siril command processing *
*************************************************************/

static PyObject* siril_processcommand(PyObject* self, PyObject* args) {
	Py_ssize_t num_args = PyTuple_Size(args);
	if (num_args < 1) {
		PyErr_SetString(PyExc_TypeError, _("At least one argument is required"));
		return NULL;
	}

	// Process the first argument (command)
	PyObject* command_obj = PyTuple_GetItem(args, 0);  // Borrowed reference
	if (!command_obj) {
		PyErr_SetString(PyExc_TypeError, _("Failed to get the command argument"));
		return NULL;
	}

	gchar* command = NULL;
	if (PyUnicode_Check(command_obj)) {
		command = g_strdup(PyUnicode_AsUTF8(command_obj));
	} else if (PyBytes_Check(command_obj)) {
		command = g_strdup(PyBytes_AsString(command_obj));
	} else {
		PyErr_SetString(PyExc_TypeError, _("Command must be a string"));
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
				PyErr_SetString(PyExc_ValueError, _("Failed to process argument"));
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

/*****************************************************
* Functions providing access to important variables *
****************************************************/

// Function to return com.wd
static PyObject *siril_get_wd(PyObject *self, PyObject *args) {
	if (com.wd == NULL) {
		Py_RETURN_NONE;  // If com.wd is NULL, return None
	}
	return PyUnicode_FromString(com.wd);
}

// Function to return the current image filename
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
	{"filename", (PyCFunction)siril_get_filename, METH_NOARGS, N_("Get the current image filename")},
	{"gui_block", py_gui_block, METH_NOARGS, N_("Block the GUI by disabling script widgets")},
	{"gui_unblock", py_gui_unblock, METH_NOARGS, N_("Unblock the GUI by enabling script widgets")},
	{"log", siril_log_message_wrapper, METH_VARARGS, N_("Log a message")},
    {"notify_gfit_modified", (PyCFunction)PyNotifyGfitModified, METH_NOARGS, N_("Notify that the main image has been modified.")},
	{"cmd", siril_processcommand, METH_VARARGS, N_("Execute a Siril command")},
	{"wd", (PyCFunction)siril_get_wd, METH_NOARGS, N_("Get the current working directory")},
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
	const char *script_path = (const char*) p; // must not be freed, it is owned by the list of script menu items
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();  // Acquire the GIL

	FILE *fp = g_fopen(script_path, "r");
	int retval = -1;
	if (fp) {
		// Create a Python file object from the C file pointer
		PyObject *py_file = PyFile_FromFd(fileno(fp), script_path, "r", -1, NULL, NULL, NULL, 1);
		if (py_file != NULL) {
			// Run the script and catch any exceptions
			PyObject *main_module = PyImport_AddModule("__main__");
			PyObject *globals = PyModule_GetDict(main_module);
			PyObject *result = PyRun_File(fp, script_path, Py_file_input, globals, globals);

			if (result == NULL) {
				// An exception occurred
				PyObject *type, *value, *traceback;
				PyErr_Fetch(&type, &value, &traceback);

				// Convert the error to a string
				PyObject *str_exc_value = PyObject_Str(value);
				const char* err_msg = PyUnicode_AsUTF8(str_exc_value);

				// Log the error
				siril_log_color_message(_("Error in Python script: %s\n"), "red", err_msg);

				Py_XDECREF(str_exc_value);
				Py_XDECREF(type);
				Py_XDECREF(value);
				Py_XDECREF(traceback);

				retval = -1;
			} else {
				// Script completed successfully
				Py_DECREF(result);
				retval = 0;
			}

			Py_DECREF(py_file);
		} else {
			siril_log_color_message(_("Failed to create Python file object from: %s\n"), "red", script_path);
			retval = -1;
		}
		fclose(fp);
	} else {
		siril_log_color_message(_("Failed to open script file: %s\n"), "red", script_path);
		retval = -1;
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

	int retval = 0;
	PyObject *main_module = PyImport_AddModule("__main__");
	if (main_module == NULL) {
		siril_log_message(_("Failed to get __main__ module\n"));
		retval = -1;
	} else {
		PyObject *globals = PyModule_GetDict(main_module);
		PyObject *result = PyRun_String(script_contents, Py_file_input, globals, globals);

		if (result == NULL) {
			// An exception occurred
			PyObject *type, *value, *traceback;
			PyErr_Fetch(&type, &value, &traceback);

			// Convert the error to a string
			PyObject *str_exc_value = PyObject_Str(value);
			const char* err_msg = PyUnicode_AsUTF8(str_exc_value);

			// Log the error
			siril_log_message(_("Error in Python script (memory): %s\n"), err_msg);

			Py_XDECREF(str_exc_value);
			Py_XDECREF(type);
			Py_XDECREF(value);
			Py_XDECREF(traceback);

			retval = -1;
		} else {
			// Script completed successfully
			Py_DECREF(result);
			retval = 0;
		}
	}

	PyGILState_Release(gstate);  // Release the GIL
	// Note: script_widgets_idle() call is omitted as per the original comment
	return GINT_TO_POINTER(retval);
}

// Function to finalize Python interpreter
void finalize_python(void) {
	Py_Finalize();
	siril_log_message(_("Python scripting module cleaned up.\n"));
}
