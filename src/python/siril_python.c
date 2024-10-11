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

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "python/siril_python_functions.h"
#include "gui/script_menu.h"
#include "core/siril_log.h"

// Define getters and setters
static PyGetSetDef PyFits_getsetters[] = {
	// rx, ry, nchans have getters only: the fits size should not be changed using these properties
	{"rx", (getter)PyFits_get_rx, NULL, N_("image width"), NULL},
	{"ry", (getter)PyFits_get_ry, NULL, N_("image height"), NULL},
	{"nchans", (getter)PyFits_get_nchans, NULL, N_("image channel depth"), NULL},
	// These don't have setters either, as they are computed values
	{"mini", (getter)PyFits_get_mini, NULL, N_("min value across all channels"), NULL},
	{"maxi", (getter)PyFits_get_maxi, NULL, N_("max value across all channels"), NULL},
	{"neg_ratio", (getter)PyFits_get_neg_ratio, NULL, N_("ratio of negative to nonnegative pixels"), NULL},
	{"top_down", (getter)PyFits_get_top_down, NULL, N_("native roworder of the original file format"), NULL},
	{"roworder", (getter)PyFits_get_row_order, NULL, N_("row order of the sensor that produced the image"), NULL},
	// These could have setters, but the Siril python module only provides direct read-only access in Siril 1.4
	// Modification of images must be carried out using Siril commands and the siril.cmd() function
	{"bitdepth", (getter)PyFits_get_bitdepth,NULL, N_("image bitdepth"), NULL},
	{"bscale", (getter)PyFits_get_bscale, NULL, N_("bscale value"), NULL},
	{"bzero", (getter)PyFits_get_bzero, NULL, N_("bzero value"), NULL},
	{"lo", (getter)PyFits_get_lo, NULL, N_("Lower visualization cutoff (WORD or float as appropriate to the bitdepth)"), NULL},
	{"hi", (getter)PyFits_get_hi, NULL, N_("Upper visualization cutoff (WORD or float as appropriate to the bitdepth)"), NULL},
	{"program", (getter)PyFits_get_program, NULL, N_("Software that created this HDU"), NULL},
	{"filename", (getter)PyFits_get_filename, NULL, N_("Original Filename"), NULL},
	{"data_max", (getter)PyFits_get_data_max, NULL, N_("Maximum data value"), NULL},
	{"data_min", (getter)PyFits_get_data_min, NULL, N_("Minimum data value"), NULL},
	{"pixel_size_x", (getter)PyFits_get_pixel_size_x, NULL, N_("Pixel size in X direction"), NULL},
	{"pixel_size_y", (getter)PyFits_get_pixel_size_y, NULL, N_("Pixel size in Y direction"), NULL},
	{"binning_x", (getter)PyFits_get_binning_x, NULL, N_("Binning in X direction"), NULL},
	{"binning_y", (getter)PyFits_get_binning_y, NULL, N_("Binning in Y direction"), NULL},
	{"row_order", (getter)PyFits_get_row_order, NULL, N_("Row order"), NULL},
	{"date", (getter)PyFits_get_date, NULL, N_("Creation date (UTC)"), NULL},
	{"date_obs", (getter)PyFits_get_date_obs, NULL, N_("Observation date (UTC)"), NULL},
	{"expstart", (getter)PyFits_get_expstart, NULL, N_("Exposure start time (Julian date)"), NULL},
	{"expend", (getter)PyFits_get_expend, NULL, N_("Exposure end time (Julian date)"), NULL},
	{"filter", (getter)PyFits_get_filter, NULL, N_("Filter used"), NULL},
	{"image_type", (getter)PyFits_get_image_type, NULL, N_("Image type"), NULL},
	{"object", (getter)PyFits_get_object, NULL, N_("Object name"), NULL},
	{"instrume", (getter)PyFits_get_instrume, NULL, N_("Instrument name"), NULL},
	{"telescop", (getter)PyFits_get_telescop, NULL, N_("Telescope name"), NULL},
	{"observer", (getter)PyFits_get_observer, NULL, N_("Observer name"), NULL},
	{"centalt", (getter)PyFits_get_centalt, NULL, N_("Center altitude"), NULL},
	{"centaz", (getter)PyFits_get_centaz, NULL, N_("Center azimuth"), NULL},
	{"sitelat", (getter)PyFits_get_sitelat, NULL, N_("Site latitude"), NULL},
	{"sitelong", (getter)PyFits_get_sitelong, NULL, N_("Site longitude"), NULL},
	{"sitelat_str", (getter)PyFits_get_sitelat_str, NULL, N_("Site latitude (string)"), NULL},
	{"sitelong_str", (getter)PyFits_get_sitelong_str, NULL, N_("Site longitude (string)"), NULL},
	{"siteelev", (getter)PyFits_get_siteelev, NULL, N_("Site elevation"), NULL},
	{"bayer_pattern", (getter)PyFits_get_bayer_pattern, NULL, N_("Bayer pattern"), NULL},
	{"bayer_xoffset", (getter)PyFits_get_bayer_xoffset, NULL, N_("Bayer X offset"), NULL},
	{"bayer_yoffset", (getter)PyFits_get_bayer_yoffset, NULL, N_("Bayer Y offset"), NULL},
	{"airmass", (getter)PyFits_get_airmass, NULL, N_("Relative optical path length through atmosphere"), NULL},
	{"focal_length", (getter)PyFits_get_focal_length, NULL, N_("Focal length"), NULL},
	{"flength", (getter)PyFits_get_flength, NULL, N_("Focal length (alternative)"), NULL},
	{"iso_speed", (getter)PyFits_get_iso_speed, NULL, N_("ISO speed"), NULL},
	{"exposure", (getter)PyFits_get_exposure, NULL, N_("Exposure time"), NULL},
	{"aperture", (getter)PyFits_get_aperture, NULL, N_("Aperture"), NULL},
	{"ccd_temp", (getter)PyFits_get_ccd_temp, NULL, N_("CCD temperature"), NULL},
	{"set_temp", (getter)PyFits_get_set_temp, NULL, N_("Set temperature"), NULL},
	{"livetime", (getter)PyFits_get_livetime, NULL, N_("Sum of exposures"), NULL},
	{"stackcnt", (getter)PyFits_get_stackcnt, NULL, N_("Number of stacked frames"), NULL},
	{"cvf", (getter)PyFits_get_cvf, NULL, N_("Conversion factor (e-/ADU)"), NULL},
	{"key_gain", (getter)PyFits_get_key_gain, NULL, N_("Gain value from camera headers"), NULL},
	{"key_offset", (getter)PyFits_get_key_offset, NULL, N_("Offset value from camera headers"), NULL},
	{"focname", (getter)PyFits_get_focname, NULL, N_("Focuser name"), NULL},
	{"focuspos", (getter)PyFits_get_focuspos, NULL, N_("Focuser position"), NULL},
	{"focussz", (getter)PyFits_get_focussz, NULL, N_("Focuser size"), NULL},
	{"foctemp", (getter)PyFits_get_foctemp, NULL, N_("Focuser temperature"), NULL},
	{"header", (getter)PyFits_get_header, NULL, N_("FITS header"), NULL},
	{"unknown_keys", (getter)PyFits_get_unknown_keys, NULL, N_("Unknown keys"), NULL},
	{"icc_profile", (getter)PyFits_get_icc_profile, NULL, N_("ICC profile (as PyBytes)"), NULL},
	{"history", (getter)PyFits_get_history, NULL, N_("History (as a list of strings)"), NULL},
};

// Define methods for PyFits
static PyMethodDef PyFits_methods[] = {
	{"get_config_item", (PyCFunction)PyFits_get_config_item, METH_CLASS | METH_VARARGS, N_("Get a config item")},
	{"image", (PyCFunction)PyFits_gfit, METH_CLASS | METH_NOARGS, N_("Get the Siril main image as a PyFits object")},
	{"get_total", (PyCFunction)PyFits_get_total, METH_VARARGS, N_("Return the total pixel count for the specified channel.")},
	{"get_ngoodpix", (PyCFunction)PyFits_get_ngoodpix, METH_VARARGS, N_("Return the number of good pixels for the specified channel.")},
	{"get_mean", (PyCFunction)PyFits_get_mean, METH_VARARGS, N_("Return the mean pixel value for the specified channel.")},
	{"get_median", (PyCFunction)PyFits_get_median, METH_VARARGS, N_("Return the median pixel value for the specified channel.")},
	{"get_sigma", (PyCFunction)PyFits_get_sigma, METH_VARARGS, N_("Return the sigma (standard deviation) for the specified channel.")},
	{"get_avgdev", (PyCFunction)PyFits_get_avgdev, METH_VARARGS, N_("Return the average deviation for the specified channel.")},
	{"get_mad", (PyCFunction)PyFits_get_mad, METH_VARARGS, N_("Return the median absolute deviation (MAD) for the specified channel.")},
	{"get_sqrtbwmv", (PyCFunction)PyFits_get_sqrtbwmv, METH_VARARGS, N_("Return the square root of the biweight midvariance for the specified channel.")},
	{"get_location", (PyCFunction)PyFits_get_location, METH_VARARGS, N_("Return the location value for the specified channel.")},
	{"get_scale", (PyCFunction)PyFits_get_scale, METH_VARARGS, N_("Return the scale value for the specified channel.")},
	{"get_min", (PyCFunction)PyFits_get_min, METH_VARARGS, N_("Return the minimum pixel value for the specified channel.")},
	{"get_max", (PyCFunction)PyFits_get_max, METH_VARARGS, N_("Return the maximum pixel value for the specified channel.")},
	{"get_normvalue", (PyCFunction)PyFits_get_normvalue, METH_VARARGS, N_("Return the normalization value for the specified channel.")},
	{"get_bgnoise", (PyCFunction)PyFits_get_bgnoise, METH_VARARGS, N_("Return the background noise value for the specified channel.")},
	{NULL}
};

// Define the PyFits type
PyTypeObject PyFitsType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.fits",
	.tp_doc = N_("Siril fits object"),
	.tp_basicsize = sizeof(PyFits),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_new = PyFits_new,
	.tp_init = (initproc)PyFits_init,
	.tp_dealloc = (destructor)PyFits_dealloc,
    .tp_as_buffer = &PyFits_as_buffer,
	.tp_methods = PyFits_methods,
	.tp_getset = PyFits_getsetters,
};

// Define methods for the module
static PyMethodDef SirilMethods[] = {
	{"filename", (PyCFunction)siril_get_filename, METH_NOARGS, N_("Get the current image filename")},
	{"gui_block", py_gui_block, METH_NOARGS, N_("Block the GUI by disabling widgets except for Stop")},
	{"gui_unblock", py_gui_unblock, METH_NOARGS, N_("Unblock the GUI by enabling all widgets")},
	{"log", siril_log_message_wrapper, METH_VARARGS, N_("Log a message")},
	{"notify_image_modified", (PyCFunction)PyNotifyGfitModified, METH_NOARGS, N_("Notify that the main image has been modified.")},
	{"cmd", siril_processcommand, METH_VARARGS, N_("Execute a Siril command")},
	{"wd", (PyCFunction)siril_get_wd, METH_NOARGS, N_("Get the current working directory")},
	{NULL, NULL, 0, NULL}  /* Sentinel */
};

// Module definition
static struct PyModuleDef sirilmodule = {
    PyModuleDef_HEAD_INIT,
    "siril",
    N_("Siril Python API"),
    -1,
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
	if (PyModule_AddObject(m, "fits", (PyObject *) &PyFitsType) < 0) {
		Py_DECREF(&PyFitsType);
		Py_DECREF(m);
		return NULL;
	}

	return m;
}

// Functions to do with initializing and finalizing the interpreter

static gboolean check_or_create_python_venv(const char *venv_dir, gboolean *already_active) {
    const char *current_venv = g_getenv("VIRTUAL_ENV");
	*already_active = FALSE;
    if (current_venv) {
        siril_log_message(_("A virtual environment is already active: %s. Using the active virtual environment.\n"), current_venv);
		*already_active = TRUE;
        return TRUE;
    }

    char *venv_python = NULL;

    #ifdef _WIN32
    venv_python = g_build_filename(venv_dir, "Scripts", "python.exe", NULL);
    #else
    venv_python = g_build_filename(venv_dir, "bin", "python", NULL);
    #endif

    gboolean venv_exists = g_file_test(venv_python, G_FILE_TEST_IS_EXECUTABLE);
    g_free(venv_python);

    if (venv_exists) {
        siril_debug_print("A virtual environment already exists: %s.\n", current_venv);
        return TRUE;
    }

    // Create a new venv
    char *python_path = NULL;
    #ifdef _WIN32
    python_path = g_find_program_in_path("python.exe");
    #else
    python_path = g_find_program_in_path("python3");
    if (!python_path) {
        python_path = g_find_program_in_path("python");
    }
    #endif

    if (!python_path) {
        siril_log_message(_("Python not found. Cannot create virtual environment.\n"));
        return FALSE;
    }

    char *argv[] = {python_path, "-m", "venv", (char *)venv_dir, NULL};
    gchar *stdout_output = NULL, *stderr_output = NULL;
    GError *error = NULL;
    gint exit_status;

    gboolean success = g_spawn_sync(
        NULL, argv, NULL, G_SPAWN_SEARCH_PATH, NULL, NULL,
        &stdout_output, &stderr_output, &exit_status, &error
    );

    g_free(python_path);

    if (!success) {
        siril_log_message(_("Failed to create Python virtual environment: %s"), error->message);
        g_clear_error(&error);
        g_free(stdout_output);
        g_free(stderr_output);
        return FALSE;
    }

    if (exit_status != 0) {
        siril_log_message(_("Python virtual environment creation failed with exit code %d: %s"), exit_status, stderr_output);
        g_free(stdout_output);
        g_free(stderr_output);
        return FALSE;
    }
	siril_log_message(_("Created Python virtual environment: %s\n"), venv_dir);

    g_free(stdout_output);
    g_free(stderr_output);

    return TRUE;
}

static gboolean activate_python_venv(const char *venv_dir) {
    // Set environment variables to activate the venv
    g_setenv("VIRTUAL_ENV", venv_dir, TRUE);

    char *new_path = NULL;
    #ifdef _WIN32
    char *scripts_dir = g_build_filename(venv_dir, "Scripts", NULL);
    #else
    char *scripts_dir = g_build_filename(venv_dir, "bin", NULL);
    #endif

    const char *old_path = g_getenv("PATH");
    if (old_path) {
        new_path = g_strdup_printf("%s%c%s", scripts_dir, G_SEARCHPATH_SEPARATOR, old_path);
    } else {
        new_path = g_strdup(scripts_dir);
    }
    g_setenv("PATH", new_path, TRUE);

    g_free(scripts_dir);
    g_free(new_path);

    // Unset PYTHONHOME if it's set
    g_unsetenv("PYTHONHOME");

    // Update sys.prefix and sys.exec_prefix
    PyObject *sys_module = PyImport_ImportModule("sys");
    if (sys_module) {
        PyObject *py_venv_dir = PyUnicode_DecodeFSDefault(venv_dir);
        if (py_venv_dir) {
            PyObject_SetAttrString(sys_module, "prefix", py_venv_dir);
            PyObject_SetAttrString(sys_module, "exec_prefix", py_venv_dir);
            Py_DECREF(py_venv_dir);
        } else {
            PyErr_Print();
        }
        Py_DECREF(sys_module);
    } else {
        PyErr_Print();
    }

    // Update sys.path
    if (PyRun_SimpleString("import sys; sys.path = [p for p in sys.path if 'site-packages' not in p]") != 0) {
        PyErr_Print();
    }
    if (PyRun_SimpleString("import site; site.main()") != 0) {
        PyErr_Print();
    }

    siril_log_message(_("Activated Python virtual environment: %s\n"), venv_dir);
    return TRUE;
}

// Function to initialize Python interpreter and load our module
gpointer init_python(void *user_data) {
	gchar* venv_dir = g_build_filename(g_get_user_data_dir(), "siril", "venv", NULL);
	gboolean already_active;
	gboolean venv_created = check_or_create_python_venv(venv_dir, &already_active);
	PyImport_AppendInittab("siril", PyInit_siril);
	Py_Initialize();
	if (venv_created && !already_active)
		activate_python_venv(venv_dir);
	g_free(venv_dir);
	PyImport_ImportModule("siril");
	PyEval_SaveThread();  // Save the current thread state and release the GIL
	// Update the global thread reference
	com.python_thread = g_thread_self();
	siril_log_message(_("Python scripting module initialized.\n"));
	// Create a GMainContext and GMainLoop for this thread
	com.python_context = g_main_context_new();
	com.python_loop = g_main_loop_new(com.python_context, FALSE);
	// Enter the main loop (this will keep the thread alive and waiting for tasks)
	g_main_loop_run(com.python_loop);
	// Finalize Python interpreter when done
	Py_Finalize();
	return NULL;
}

// Function to run a Python script from a file
gboolean run_python_script_from_file(gpointer p) {
	const char *script_path = (const char*) p; // must not be freed, it is owned by the list of script menu items
	PyGILState_STATE gstate;
	com.script = TRUE;
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

				retval = FALSE;
			} else {
				// Script completed successfully
				Py_DECREF(result);
				retval = FALSE;
			}
			Py_DECREF(py_file);
		} else {
			siril_log_color_message(_("Failed to create Python file object from: %s\n"), "red", script_path);
			retval = FALSE;
		}
		fclose(fp);
	} else {
		siril_log_color_message(_("Failed to open script file: %s\n"), "red", script_path);
		retval = FALSE;
	}
	PyGC_Collect(); // Force garbage collection, in case the script didn't bother
	PyGILState_Release(gstate);  // Release the GIL
	com.script = FALSE;
	g_idle_add(script_widgets_idle, NULL);
	return retval;
}

// Function to run a Python script from memory
gboolean run_python_script_from_mem(gpointer p) {
	const char *script_contents = (const char*) p;
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();  // Acquire the GIL
	com.script = TRUE;
	int retval = FALSE;
	PyObject *main_module = PyImport_AddModule("__main__");
	if (main_module == NULL) {
		siril_log_message(_("Failed to get __main__ module\n"));
		retval = FALSE;
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

			retval = FALSE;
		} else {
			// Script completed successfully
			Py_DECREF(result);
			retval = FALSE;
		}
	}
	PyGC_Collect(); // Force garbage collection, in case the script didn't bother
	PyGILState_Release(gstate);  // Release the GIL
	com.script = FALSE;
	// Note: script_widgets_idle() call is omitted as per the original comment
	return retval;
}

// Function to run Python script, delegating to the Python thread if necessary
void run_python_script_in_python_thread(const char *script, gboolean from_file) {
    // Check if we're in the Python thread
    if (g_thread_self() == com.python_thread) {
        // We're already in the Python thread, so just run the script directly
        if (from_file) {
			run_python_script_from_file((gpointer) script);
		} else {
			run_python_script_from_mem((gpointer) script);
		}
    } else {
        // Use g_main_context_invoke() to schedule the script execution
        // in the Python thread's main loop/context
		if (from_file) {
			g_main_context_invoke(com.python_context, run_python_script_from_file, (gpointer) script);
		} else {
			g_main_context_invoke(com.python_context, run_python_script_from_mem, (gpointer) script);
		}
    }
}

// Function to finalize Python interpreter
void finalize_python(void) {
	Py_Finalize();
	siril_log_message(_("Python scripting module cleaned up.\n"));
}
