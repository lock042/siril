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
#include "core/command_line_processor.h"
#include "gui/script_menu.h"
#include "gui/progress_and_log.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "algos/statistics.h"

/******************
******************
**              **
** siril module **
**              **
******************
*****************/

static PyObject *get_config_item_as_pyobject(char *input) {
	/* parsing a single variable command */
	int sep, len = strlen(input);
	for (sep = 1; sep < len; sep++)
		if (input[sep] == '.')
			break;
	if (sep == len) {
		siril_log_message(_("syntax: group.key=value\n"));
		Py_RETURN_NONE;
	}
	input[sep] = '\0';
	gchar *str = get_settings_key(input, input+sep+1, FALSE);
	PyObject *retval = PyUnicode_FromString(str);
	g_free(str);
	return retval;
}

PyObject *siril_get_config_item(PyObject *self, PyObject *args) {
	char *input_str;

	// Parse the Python argument as a string
	if (!PyArg_ParseTuple(args, "s", &input_str)) {
		PyErr_SetString(PyExc_TypeError, N_("Expected a string argument"));
		return NULL;
	}

	// Check if the string is NULL or empty
	if (input_str == NULL || strlen(input_str) == 0) {
		PyErr_SetString(PyExc_ValueError, N_("Input string cannot be NULL or empty"));
		return NULL;
	}

	// Call the python_get function with the input string
	PyObject *result = get_config_item_as_pyobject(input_str);

	// Check if python_get returned a valid PyObject
	if (result == NULL) {
		// If NULL was returned, raise an error in Python
		PyErr_SetString(PyExc_RuntimeError, N_("Failed to get configuration item"));
		return NULL;
	}

	// Return the result from python_get (which is already a PyObject)
	return result;
}

/*************************************************
* Functions providing control over the Siril UI *
************************************************/

PyObject* py_gui_block(PyObject* self, PyObject* args) {
	if (!g_main_context_iteration(NULL, FALSE)) {
		siril_log_color_message(_("Warning: siril.gui_block() must not be called except from a script's GTK main loop.\n"), "red");
		Py_RETURN_NONE;
	}
	script_widgets_enable(FALSE);  // Disable main control window GUI elements
	Py_RETURN_NONE;
}

PyObject* py_gui_unblock(PyObject* self, PyObject* args) {
	if (!g_main_context_iteration(NULL, FALSE)) {
		siril_log_color_message(_("Warning: siril.gui_unblock() must not be called except from a script's GTK main loop.\n"), "red");
		Py_RETURN_NONE;
	}
	script_widgets_enable(TRUE);   // Enable main control window GUI elements
	Py_RETURN_NONE;
}

PyObject* siril_notify_gfit_modified(PyObject* self, PyObject* args) {
	siril_debug_print("end of gfit operation\n");
	notify_gfit_modified();
	Py_RETURN_NONE;
}

// Method to call siril_log_message
PyObject *siril_log_message_wrapper(PyObject *self, PyObject *args) {
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

PyObject* siril_progress_wrapper(PyObject* self, PyObject* args) {
	const char* text;
	double progress;

	if (!PyArg_ParseTuple(args, "sd", &text, &progress)) {
		return NULL;
	}

	// Call the set_progress_bar_data function
	set_progress_bar_data(text, progress);

	// Return None (Python's equivalent of void)
	Py_RETURN_NONE;
}

/************************************************************
* Functions providing control over Siril command processing *
************************************************************/

PyObject* siril_processcommand(PyObject* self, PyObject* args) {
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
	siril_debug_print("starting command\n");
	int result = processcommand(full_command->str, TRUE);
	siril_debug_print("ended command\n");
	g_string_free(full_command, TRUE);

	return PyLong_FromLong(result);
}

/*****************************************************
* Functions providing access to important variables *
****************************************************/

// Function to return com.wd
PyObject *siril_get_wd(PyObject *self, PyObject *args) {
	if (com.wd == NULL) {
		Py_RETURN_NONE;  // If com.wd is NULL, return None
	}
	return PyUnicode_FromString(com.wd);
}

// Function to return the current image filename
PyObject *siril_get_filename(PyObject *self, PyObject *args) {
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

// Property-style getter "continue" for long running scripts to check
// periodically to see if the user has requested them to stop
PyObject* siril_get_continue(PyObject* self, PyObject* args) {
	return PyBool_FromLong(com.stop_script);
}

PyObject* siril_pipinstall(PyObject* self, PyObject* args) {
	const char* module_name;
	if (!PyArg_ParseTuple(args, "s", &module_name)) {
		return NULL;
	}

	// Import pip
	PyObject* pip_module = PyImport_ImportModule("pip");
	if (pip_module == NULL) {
		PyErr_SetString(PyExc_ImportError, "Failed to import pip. Make sure it's installed in the virtual environment.");
		return NULL;
	}

	// Prepare arguments for pip install
	PyObject* pip_args = Py_BuildValue("[sss]", "install", "--upgrade", module_name);
	if (pip_args == NULL) {
		Py_DECREF(pip_module);
		return NULL;
	}

	// Get pip.main function
	PyObject* pip_main = PyObject_GetAttrString(pip_module, "main");
	if (pip_main == NULL) {
		Py_DECREF(pip_module);
		Py_DECREF(pip_args);
		PyErr_SetString(PyExc_AttributeError, "Failed to get pip.main function.");
		return NULL;
	}

	// Call pip.main(["install", "--upgrade", module_name])
	PyObject* result = PyObject_CallFunctionObjArgs(pip_main, pip_args, NULL);

	// Clean up
	Py_DECREF(pip_module);
	Py_DECREF(pip_main);
	Py_DECREF(pip_args);

	if (result == NULL) {
		return NULL;  // Python exception already set
	}

	int exit_code = PyLong_AsLong(result);
	Py_DECREF(result);

	if (exit_code != 0) {
		PyErr_Format(PyExc_RuntimeError, "pip install failed with exit code %d", exit_code);
		return NULL;
	}

	Py_RETURN_NONE;
}
