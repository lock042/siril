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
#include <datetime.h>

#include "core/siril.h"
#include "python/siril_python.h"
#include "python/PyImStats_functions.h"
#include "core/command_line_processor.h"
#include "gui/script_menu.h"
#include "gui/progress_and_log.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "algos/statistics.h"

// Helper function to convert GDateTime to PyDateTime
PyObject* gdatetime_to_pydatetime(GDateTime *gdt) {
	// Ensure the datetime API is ready to use
	if (!PyDateTimeAPI) {
		PyDateTime_IMPORT;
	}

	if (!gdt) {
		siril_log_message(_("FITS does not contain a valid DateTime\n"), "red");
		Py_RETURN_NONE;
	}

	// Extract fields from GDateTime
	int year = g_date_time_get_year(gdt);
	int month = g_date_time_get_month(gdt);
	int day = g_date_time_get_day_of_month(gdt);
	int hour = g_date_time_get_hour(gdt);
	int minute = g_date_time_get_minute(gdt);
	int second = g_date_time_get_second(gdt);
	int microsecond = g_date_time_get_microsecond(gdt);

	// Create a PyDateTime object using the Python C API
	PyObject *py_datetime = PyDateTime_FromDateAndTime(
		year, month, day, hour, minute, second, microsecond
	);

	if (!py_datetime) {
		PyErr_SetString(PyExc_RuntimeError, N_("Failed to create Python datetime object"));
		return NULL;
	}

	// Return the new PyDateTime object
	return py_datetime;
}

// Deallocation function for PyFits
void PyFits_dealloc(PyFits *self) {
	if (self->fit != NULL && self->should_free_data) {
		clearfits(self->fit);  // Always clear the internal structures
		if (self->should_free) {
			free(self->fit);  // Only free the pointer if it was dynamically allocated
		}
	}
	Py_TYPE(self)->tp_free((PyObject *)self);
}

// New function for PyFits
PyObject *PyFits_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
	PyFits *self;
	self = (PyFits *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->fit = NULL;
		self->should_free_data = 0;
		self->should_free = 0;
	}
	return (PyObject *)self;
}

// Initialize function for PyFits
int PyFits_init(PyFits *self, PyObject *args, PyObject *kwds) {
	static char *kwlist[] = {"width", "height", "nblayer", "type", NULL};
	int width = 0, height = 0, nblayer = 1;
	data_type type = DATA_USHORT;  // Assuming DATA_USHORT is the default

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiii", kwlist,
									&width, &height, &nblayer, &type))
		return -1;

	if (width > 0 && height > 0) {
		// Create a new fits structure
		self->fit = calloc(1, sizeof(fits));
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
		self->should_free_data = 1;  // We generally want to clean up after ourselves (not for gfit though, maybe other exeptions)
	}

	return 0;
}

static int PyFits_getbuffer(PyObject *obj, Py_buffer *view, int flags) {
	PyFits *self = (PyFits *)obj;
	fits *fit = self->fit;

	if (!fit) {
		PyErr_SetString(PyExc_ValueError, "Fits object is NULL");
		return -1;
	}

	// Determine dimensionality
	int ndim = (fit->naxes[2] > 1) ? 3 : 2;
	view->ndim = ndim;
	view->readonly = 0;  // Assume writable for now

	// Allocate memory for shape and strides
	Py_ssize_t *shape = PyMem_Malloc(ndim * sizeof(Py_ssize_t));
	Py_ssize_t *strides = PyMem_Malloc(ndim * sizeof(Py_ssize_t));

	if (!shape || !strides) {
		PyErr_NoMemory();
		PyMem_Free(shape);
		PyMem_Free(strides);
		return -1;
	}

	// Set shape
	if (ndim == 3) {
		shape[0] = fit->naxes[2];
		shape[1] = fit->naxes[1];
		shape[2] = fit->naxes[0];
	} else {  // ndim == 2
		shape[0] = fit->naxes[1];
		shape[1] = fit->naxes[0];
	}
	view->shape = shape;

	// Set buffer and itemsize based on data type
	if (fit->type == DATA_USHORT) {
		view->buf = fit->data;
		view->itemsize = sizeof(WORD);
		view->format = "H";  // Unsigned short
	} else if (fit->type == DATA_FLOAT) {
		view->buf = fit->fdata;
		view->itemsize = sizeof(float);
		view->format = "f";  // Float
	} else {
		PyErr_SetString(PyExc_ValueError, "Unsupported data type");
		PyMem_Free(shape);
		PyMem_Free(strides);
		return -1;
	}

	// Calculate total length
	view->len = view->itemsize * fit->naxes[0] * fit->naxes[1] * fit->naxes[2];

	// Set strides for row-major order
	if (ndim == 3) {
		strides[0] = fit->naxes[1] * fit->naxes[0] * view->itemsize;
		strides[1] = fit->naxes[0] * view->itemsize;
		strides[2] = view->itemsize;
	} else {  // ndim == 2
		strides[0] = fit->naxes[0] * view->itemsize;
		strides[1] = view->itemsize;
	}
	view->strides = strides;

	view->suboffsets = NULL;
	view->internal = NULL;

	// Check if writable buffer is requested but data is read-only
	if ((flags & PyBUF_WRITABLE) && view->readonly) {
		PyErr_SetString(PyExc_BufferError, "Object is not writable");
		PyMem_Free(shape);
		PyMem_Free(strides);
		return -1;
	}

	return 0;
}

static void PyFits_releasebuffer(PyObject *obj, Py_buffer *view) {
	// Free the shape and strides arrays we allocated in getbuffer
	PyMem_Free(view->shape);
	PyMem_Free(view->strides);

	// Clear the Py_buffer struct without freeing the actual data
	// which is still owned by the PyFits (or by Siril)
	memset(view, 0, sizeof(Py_buffer));
}

PyBufferProcs PyFits_as_buffer = {
	(getbufferproc)PyFits_getbuffer,
	(releasebufferproc)PyFits_releasebuffer  // Standard release function
};

// Method to get rx
PyObject *PyFits_get_rx(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->rx);
}

// Method to get ry
PyObject *PyFits_get_ry(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->ry);
}

// Method to get naxes[2] (number of channels)
PyObject *PyFits_get_nchans(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->naxes[2]);
}

// Method to get mini (min pixel value across all channels)
PyObject *PyFits_get_mini(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->mini);  // Changed to return a float
}

// Method to get negratio
PyObject *PyFits_get_neg_ratio(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->neg_ratio);  // Changed to return a float
}

// Method to get maxi (max pixel value across all channels)
PyObject *PyFits_get_maxi(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->maxi);
}

// Method to get top_down (boolean indicating top-down orientation)
PyObject *PyFits_get_top_down(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	// Return Python boolean: 1 for TRUE, 0 for FALSE
	return PyBool_FromLong(self->fit->top_down ? 1 : 0);
}

// Method to get bit depth
PyObject *PyFits_get_bitdepth(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->type == DATA_FLOAT ? 32 : 16);
}

// Method to get bscale
PyObject *PyFits_get_bscale(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.bscale);
}

// Method to get bzero
PyObject *PyFits_get_bzero(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.bzero);
}

// Method to get lo
PyObject *PyFits_get_lo(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	if (self->fit->type == DATA_USHORT) {
		return PyLong_FromLong(self->fit->keywords.lo);
	} else if (self->fit->type == DATA_FLOAT) {
		return PyLong_FromLong(self->fit->keywords.flo);
	} else {
		PyErr_SetString(PyExc_AttributeError, _("unknown fit data type"));
		return NULL;
	}
}

// Method to get hi
PyObject *PyFits_get_hi(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	if (self->fit->type == DATA_USHORT) {
		return PyLong_FromLong(self->fit->keywords.hi);
	} else if (self->fit->type == DATA_FLOAT) {
		return PyLong_FromLong(self->fit->keywords.fhi);
	} else {
		PyErr_SetString(PyExc_AttributeError, _("unknown fit data type"));
		return NULL;
	}
}

// Method to get program
PyObject *PyFits_get_program(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.program);
}

// Method to get filename
PyObject *PyFits_get_filename(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.filename);
}

// Method to get data_max
PyObject *PyFits_get_data_max(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.data_max);
}

// Method to get data_min
PyObject *PyFits_get_data_min(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.data_min);
}

// Method to get pixel_size_x
PyObject *PyFits_get_pixel_size_x(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.pixel_size_x);
}

// Method to get pixel_size_y
PyObject *PyFits_get_pixel_size_y(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.pixel_size_y);
}

// Method to get binning_x
PyObject *PyFits_get_binning_x(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->keywords.binning_x);
}

// Method to get binning_y
PyObject *PyFits_get_binning_y(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->keywords.binning_y);
}

// Method to get row_order
PyObject *PyFits_get_row_order(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.row_order);
}

// Method to get date
PyObject *PyFits_get_date(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	// Assuming you have a way to convert GDateTime to Python datetime
	// You might need to implement this conversion
	return gdatetime_to_pydatetime(self->fit->keywords.date);
}

// Method to get date_obs
PyObject *PyFits_get_date_obs(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	// Assuming you have a way to convert GDateTime to Python datetime
	// You might need to implement this conversion
	return gdatetime_to_pydatetime(self->fit->keywords.date_obs);
}

// Method to get expstart
PyObject *PyFits_get_expstart(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.expstart);
}

// Method to get expend
PyObject *PyFits_get_expend(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.expend);
}

// Method to get filter
PyObject *PyFits_get_filter(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.filter);
}

// Method to get image_type
PyObject *PyFits_get_image_type(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.image_type);
}

// Method to get object
PyObject *PyFits_get_object(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.object);
}

// Method to get instrume
PyObject *PyFits_get_instrume(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.instrume);
}

// Method to get telescop
PyObject *PyFits_get_telescop(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.telescop);
}

// Method to get observer
PyObject *PyFits_get_observer(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.observer);
}

// Method to get centalt
PyObject *PyFits_get_centalt(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.centalt);
}

// Method to get centaz
PyObject *PyFits_get_centaz(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.centaz);
}

// Method to get sitelat
PyObject *PyFits_get_sitelat(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.sitelat);
}

// Method to get sitelong
PyObject *PyFits_get_sitelong(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.sitelong);
}

// Method to get sitelat_str
PyObject *PyFits_get_sitelat_str(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.sitelat_str);
}

// Method to get sitelong_str
PyObject *PyFits_get_sitelong_str(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.sitelong_str);
}

// Method to get siteelev
PyObject *PyFits_get_siteelev(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.siteelev);
}

// Method to get bayer_pattern
PyObject *PyFits_get_bayer_pattern(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.bayer_pattern);
}

// Method to get bayer_xoffset
PyObject *PyFits_get_bayer_xoffset(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.bayer_xoffset);
}

// Method to get bayer_yoffset
PyObject *PyFits_get_bayer_yoffset(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.bayer_yoffset);
}

// Method to get airmass
PyObject *PyFits_get_airmass(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.airmass);
}

// Method to get focal_length
PyObject *PyFits_get_focal_length(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.focal_length);
}

// Method to get flength
PyObject *PyFits_get_flength(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.flength);
}

// Method to get iso_speed
PyObject *PyFits_get_iso_speed(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.iso_speed);
}

// Method to get exposure
PyObject *PyFits_get_exposure(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.exposure);
}

// Method to get aperture
PyObject *PyFits_get_aperture(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.aperture);
}

// Method to get ccd_temp
PyObject *PyFits_get_ccd_temp(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.ccd_temp);
}

// Method to get set_temp
PyObject *PyFits_get_set_temp(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.set_temp);
}

// Method to get livetime
PyObject *PyFits_get_livetime(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.livetime);
}

// Method to get stackcnt
PyObject *PyFits_get_stackcnt(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->keywords.stackcnt);
}

// Method to get cvf
PyObject *PyFits_get_cvf(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.cvf);
}

// Method to get key_gain
PyObject *PyFits_get_key_gain(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.key_gain);
}

// Method to get key_offset
PyObject *PyFits_get_key_offset(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.key_offset);
}

// Method to get focname
PyObject *PyFits_get_focname(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.focname);
}

// Method to get focuspos
PyObject *PyFits_get_focuspos(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.focuspos);
}

// Method to get focussz
PyObject *PyFits_get_focussz(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.focussz);
}

// Method to get foctemp
PyObject *PyFits_get_foctemp(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.foctemp);
}

// Method to get header
PyObject *PyFits_get_header(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	if (self->fit->header == NULL) {
		Py_RETURN_NONE;
	}
	return PyUnicode_FromString(self->fit->header);
}

// Method to get unknown_keys
PyObject *PyFits_get_unknown_keys(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	if (self->fit->unknown_keys == NULL) {
		Py_RETURN_NONE;
	}
	return PyUnicode_FromString(self->fit->unknown_keys);
}

// Method to get ICC profile
PyObject *PyFits_get_icc_profile(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}

	cmsBool ret = FALSE;
	cmsUInt32Number length = 0;
	void *block = NULL;

	// Check if the ICC profile is available
	if (!self->fit->icc_profile) {
		Py_RETURN_NONE;
	}

	// First call to get the length of the ICC profile
	ret = cmsSaveProfileToMem(self->fit->icc_profile, NULL, &length);
	if (!ret || length == 0) {
		Py_RETURN_NONE;
	}

	// Allocate memory for the ICC profile buffer
	block = malloc(length);
	if (!block) {
		PyErr_SetString(PyExc_MemoryError, N_("Unable to allocate memory for ICC profile"));
		return NULL;
	}

	// Second call to actually save the ICC profile into the buffer
	ret = cmsSaveProfileToMem(self->fit->icc_profile, block, &length);
	if (!ret) {
		free(block);  // Free the allocated buffer on error
		PyErr_SetString(PyExc_RuntimeError, N_("Failed to save ICC profile to memory"));
		return NULL;
	}

	// Create a Python bytes object from the buffer
	PyObject *py_icc_profile = PyBytes_FromStringAndSize((const char*)block, length);

	// Free the C buffer, since PyBytes_FromStringAndSize copies the data
	free(block);

	// Return the Python bytes object
	return py_icc_profile;
}

PyObject* PyFits_get_history(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}

	GSList *history = self->fit->history;
	PyObject *py_list = PyList_New(0);
	if (py_list == NULL) {
		Py_RETURN_NONE; // There is no history. We can return None rather than throwing an exception
	}

	for (GSList *current = history; current != NULL; current = current->next) {
		char *history_item = (char *)current->data;
		PyObject *py_str = PyUnicode_FromString(history_item);
		if (py_str == NULL) {
			Py_DECREF(py_list);
			Py_RETURN_NONE;
		}
		if (PyList_Append(py_list, py_str) < 0) {
			Py_DECREF(py_str);
			Py_DECREF(py_list);
			Py_RETURN_NONE;
		}
		Py_DECREF(py_str);
	}

	return py_list;
}

// Method to access gfit
PyObject *PyFits_gfit(PyObject *cls, PyObject *args) {
	PyFits *self = (PyFits *)PyFitsType.tp_alloc(&PyFitsType, 0);
	if (self != NULL) {
		self->fit = &gfit;
		self->should_free = 0;  // gfit is statically allocated, don't free it
		self->should_free_data = 0;  // never free the data in gfit just because we don't want the view of it from python any more!
	}
	return (PyObject *)self;
}

// Stats are obtained by methods not getters as they take a channel parameter
// Stats methods are available for PyFits as well as directly for PyImStats
// because we can be more rigorous here and calculate the stats if they are
// no yet calculated: the PyImStats methods are more low-level.

// Helper function to check validity of channel index and stats
static int check_stats(PyFits *self, int n, int option) {
	if (self->fit == NULL || self->fit->stats == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("FITS data or image statistics not initialized"));
		return 0;
	}
	if (n < 0 || n >= self->fit->naxes[2]) {
		PyErr_SetString(PyExc_IndexError, _("Channel index out of range"));
		return 0;
	}
	// If the stats for the requested channel are unavailable, compute them
	if (self->fit->stats[n] == NULL) {
		statistics(NULL, -1, self->fit, n, NULL, option, MULTI_THREADED);
		if (self->fit->stats[n] == NULL) {
			PyErr_SetString(PyExc_AttributeError, _("Image statistics computation failed for this channel"));
			return 0;
		}
	}
	return 1;
}

// Getter for total
PyObject* PyFits_get_total(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyLong_FromLong(self->fit->stats[n]->total);
}

// Getter for ngoodpix
PyObject* PyFits_get_ngoodpix(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_MAIN))
		return NULL;
	// Perhaps stats were available but not computed with a high enough option, so we try recomputing
	if (self->fit->stats[n]->avgDev == NULL_STATS)
		statistics(NULL, -1, self->fit, n, NULL, STATS_MAIN, MULTI_THREADED);
	return PyLong_FromLong(self->fit->stats[n]->ngoodpix);
}

// Getter for mean
PyObject* PyFits_get_mean(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->mean);
}

// Getter for median
PyObject* PyFits_get_median(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->median);
}

// Getter for sigma
PyObject* PyFits_get_sigma(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->sigma);
}

// Getter for avgDev
PyObject* PyFits_get_avgdev(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_MAIN))
		return NULL;
	if (self->fit->stats[n]->avgDev == NULL_STATS)
		statistics(NULL, -1, self->fit, n, NULL, STATS_MAIN, MULTI_THREADED);

	return PyFloat_FromDouble(self->fit->stats[n]->avgDev);
}

// Getter for mad
PyObject* PyFits_get_mad(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_MAIN))
		return NULL;
	if (self->fit->stats[n]->avgDev == NULL_STATS)
		statistics(NULL, -1, self->fit, n, NULL, STATS_MAIN, MULTI_THREADED);
	return PyFloat_FromDouble(self->fit->stats[n]->mad);
}

// Getter for sqrtbwmv
PyObject* PyFits_get_sqrtbwmv(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_MAIN))
		return NULL;
	if (self->fit->stats[n]->avgDev == NULL_STATS)
		statistics(NULL, -1, self->fit, n, NULL, STATS_MAIN, MULTI_THREADED);
	return PyFloat_FromDouble(self->fit->stats[n]->sqrtbwmv);
}

// Getter for location
PyObject* PyFits_get_location(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_EXTRA))
		return NULL;
	if (self->fit->stats[n]->avgDev == NULL_STATS)
		statistics(NULL, -1, self->fit, n, NULL, STATS_EXTRA, MULTI_THREADED);
	return PyFloat_FromDouble(self->fit->stats[n]->location);
}

// Getter for scale
PyObject* PyFits_get_scale(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_EXTRA))
		return NULL;
	if (self->fit->stats[n]->avgDev == NULL_STATS)
		statistics(NULL, -1, self->fit, n, NULL, STATS_EXTRA, MULTI_THREADED);
	return PyFloat_FromDouble(self->fit->stats[n]->scale);
}

// Getter for min
PyObject* PyFits_get_min(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->min);
}

// Getter for max
PyObject* PyFits_get_max(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->max);
}

// Getter for normValue
PyObject* PyFits_get_normvalue(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_NORM))
		return NULL;
	if (self->fit->stats[n]->avgDev == NULL_STATS)
		statistics(NULL, -1, self->fit, n, NULL, STATS_NORM, MULTI_THREADED);
	return PyFloat_FromDouble(self->fit->stats[n]->normValue);
}

// Getter for bgnoise
PyObject* PyFits_get_bgnoise(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->bgnoise);
}

PyObject* PyFits_get_ImStats(PyFits *self, PyObject *args) {
	int channel;
	if (!PyArg_ParseTuple(args, "i", &channel))
		return NULL;

	if (channel < 0 || channel >= self->fit->naxes[2]) {
		PyErr_SetString(PyExc_IndexError, N_("Channel index out of range"));
		return NULL;
	}

	if (!self->fit->stats || !self->fit->stats[channel]) {
		Py_RETURN_NONE;
	}

	PyObject *stats = PyImStats_FromExisting(self->fit->stats[channel], (PyObject *)self, 'F');
	if (stats != NULL) {
		self->ref_count++;
	}
	return stats;
}