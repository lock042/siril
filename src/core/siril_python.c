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
#include <datetime.h>

#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/command_line_processor.h"
#include "gui/script_menu.h"
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

extern PyTypeObject PyFitsType;

// Python object to wrap our fits structure - note, this is the Siril-specific fits
// structure and is NOT interchangeable with astropy.io.fits - however a fits file
// can be saved and opened by importing astropy into a script and using it directly
// if astropy processing is desired. The result would have to be saved and then loaded
// back using siril.process_command("load", f"{filename}")

typedef struct {
	PyObject_HEAD
	fits *fit;
	int should_free_data;  // Flag to indicate if the data itself should be freed
	// note: if should_free_data is 0, the struct will not be freed either
	int should_free;  // Flag to indicate if the fits pointer itself should be freed
} PyFits;

// Deallocation function for PyFits
static void PyFits_dealloc(PyFits *self) {
	if (self->fit != NULL && self->should_free_data) {
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
		self->should_free_data = 0;
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

static PyBufferProcs PyFits_as_buffer = {
    (getbufferproc)PyFits_getbuffer,
    (releasebufferproc)PyFits_releasebuffer  // Standard release function
};

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

// Method to get bit depth
static PyObject *PyFits_get_bitdepth(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->type == DATA_FLOAT ? 32 : 16);
}

// Method to get bscale
static PyObject *PyFits_get_bscale(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.bscale);
}

// Method to get bzero
static PyObject *PyFits_get_bzero(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.bzero);
}

// Method to get lo
static PyObject *PyFits_get_lo(PyFits *self, void *closure) {
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
static PyObject *PyFits_get_hi(PyFits *self, void *closure) {
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
static PyObject *PyFits_get_program(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.program);
}

// Method to get filename
static PyObject *PyFits_get_filename(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.filename);
}

// Method to get data_max
static PyObject *PyFits_get_data_max(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.data_max);
}

// Method to get data_min
static PyObject *PyFits_get_data_min(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.data_min);
}

// Method to get pixel_size_x
static PyObject *PyFits_get_pixel_size_x(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.pixel_size_x);
}

// Method to get pixel_size_y
static PyObject *PyFits_get_pixel_size_y(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.pixel_size_y);
}

// Method to get binning_x
static PyObject *PyFits_get_binning_x(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->keywords.binning_x);
}

// Method to get binning_y
static PyObject *PyFits_get_binning_y(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->keywords.binning_y);
}

// Method to get row_order
static PyObject *PyFits_get_row_order(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.row_order);
}

// Method to get date
static PyObject *PyFits_get_date(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	// Assuming you have a way to convert GDateTime to Python datetime
	// You might need to implement this conversion
	return gdatetime_to_pydatetime(self->fit->keywords.date);
}

// Method to get date_obs
static PyObject *PyFits_get_date_obs(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	// Assuming you have a way to convert GDateTime to Python datetime
	// You might need to implement this conversion
	return gdatetime_to_pydatetime(self->fit->keywords.date_obs);
}

// Method to get expstart
static PyObject *PyFits_get_expstart(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.expstart);
}

// Method to get expend
static PyObject *PyFits_get_expend(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.expend);
}

// Method to get filter
static PyObject *PyFits_get_filter(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.filter);
}

// Method to get image_type
static PyObject *PyFits_get_image_type(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.image_type);
}

// Method to get object
static PyObject *PyFits_get_object(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.object);
}

// Method to get instrume
static PyObject *PyFits_get_instrume(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.instrume);
}

// Method to get telescop
static PyObject *PyFits_get_telescop(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.telescop);
}

// Method to get observer
static PyObject *PyFits_get_observer(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.observer);
}

// Method to get centalt
static PyObject *PyFits_get_centalt(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.centalt);
}

// Method to get centaz
static PyObject *PyFits_get_centaz(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.centaz);
}

// Method to get sitelat
static PyObject *PyFits_get_sitelat(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.sitelat);
}

// Method to get sitelong
static PyObject *PyFits_get_sitelong(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.sitelong);
}

// Method to get sitelat_str
static PyObject *PyFits_get_sitelat_str(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.sitelat_str);
}

// Method to get sitelong_str
static PyObject *PyFits_get_sitelong_str(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.sitelong_str);
}

// Method to get siteelev
static PyObject *PyFits_get_siteelev(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.siteelev);
}

// Method to get bayer_pattern
static PyObject *PyFits_get_bayer_pattern(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.bayer_pattern);
}

// Method to get bayer_xoffset
static PyObject *PyFits_get_bayer_xoffset(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.bayer_xoffset);
}

// Method to get bayer_yoffset
static PyObject *PyFits_get_bayer_yoffset(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.bayer_yoffset);
}

// Method to get airmass
static PyObject *PyFits_get_airmass(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.airmass);
}

// Method to get focal_length
static PyObject *PyFits_get_focal_length(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.focal_length);
}

// Method to get flength
static PyObject *PyFits_get_flength(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.flength);
}

// Method to get iso_speed
static PyObject *PyFits_get_iso_speed(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.iso_speed);
}

// Method to get exposure
static PyObject *PyFits_get_exposure(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.exposure);
}

// Method to get aperture
static PyObject *PyFits_get_aperture(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.aperture);
}

// Method to get ccd_temp
static PyObject *PyFits_get_ccd_temp(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.ccd_temp);
}

// Method to get set_temp
static PyObject *PyFits_get_set_temp(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.set_temp);
}

// Method to get livetime
static PyObject *PyFits_get_livetime(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.livetime);
}

// Method to get stackcnt
static PyObject *PyFits_get_stackcnt(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromUnsignedLong(self->fit->keywords.stackcnt);
}

// Method to get cvf
static PyObject *PyFits_get_cvf(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.cvf);
}

// Method to get key_gain
static PyObject *PyFits_get_key_gain(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.key_gain);
}

// Method to get key_offset
static PyObject *PyFits_get_key_offset(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.key_offset);
}

// Method to get focname
static PyObject *PyFits_get_focname(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyUnicode_FromString(self->fit->keywords.focname);
}

// Method to get focuspos
static PyObject *PyFits_get_focuspos(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.focuspos);
}

// Method to get focussz
static PyObject *PyFits_get_focussz(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyLong_FromLong(self->fit->keywords.focussz);
}

// Method to get foctemp
static PyObject *PyFits_get_foctemp(PyFits *self, void *closure) {
	if (self->fit == NULL) {
		PyErr_SetString(PyExc_AttributeError, _("fit is not initialized"));
		return NULL;
	}
	return PyFloat_FromDouble(self->fit->keywords.foctemp);
}

// Method to get header
static PyObject *PyFits_get_header(PyFits *self, void *closure) {
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
static PyObject *PyFits_get_unknown_keys(PyFits *self, void *closure) {
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
static PyObject *PyFits_get_icc_profile(PyFits *self, void *closure) {
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

static PyObject* PyFits_get_history(PyFits *self, void *closure) {
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

// Method to access gfit
static PyObject *PyFits_gfit(PyObject *cls, PyObject *args) {
	PyFits *self = (PyFits *)PyFitsType.tp_alloc(&PyFitsType, 0);
	if (self != NULL) {
		self->fit = &gfit;
		self->should_free = 0;  // gfit is statically allocated, don't free it
		self->should_free_data = 0;  // never free the data in gfit just because we don't want the view of it from python any more!
	}
	return (PyObject *)self;
}

// Stats are obtained by methods not getters as they take a channel parameter

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
static PyObject* PyFits_get_total(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyLong_FromLong(self->fit->stats[n]->total);
}

// Getter for ngoodpix
static PyObject* PyFits_get_ngoodpix(PyFits *self, PyObject *args) {
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
static PyObject* PyFits_get_mean(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->mean);
}

// Getter for median
static PyObject* PyFits_get_median(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->median);
}

// Getter for sigma
static PyObject* PyFits_get_sigma(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->sigma);
}

// Getter for avgDev
static PyObject* PyFits_get_avgdev(PyFits *self, PyObject *args) {
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
static PyObject* PyFits_get_mad(PyFits *self, PyObject *args) {
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
static PyObject* PyFits_get_sqrtbwmv(PyFits *self, PyObject *args) {
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
static PyObject* PyFits_get_location(PyFits *self, PyObject *args) {
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
static PyObject* PyFits_get_scale(PyFits *self, PyObject *args) {
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
static PyObject* PyFits_get_min(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->min);
}

// Getter for max
static PyObject* PyFits_get_max(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->max);
}

// Getter for normValue
static PyObject* PyFits_get_normvalue(PyFits *self, PyObject *args) {
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
static PyObject* PyFits_get_bgnoise(PyFits *self, PyObject *args) {
	int n;
	if (!PyArg_ParseTuple(args, "i", &n))
		return NULL;
	if (!check_stats(self, n, STATS_BASIC))
		return NULL;
	return PyFloat_FromDouble(self->fit->stats[n]->bgnoise);
}

PyObject *get_config_item_as_pyobject(char *input) {
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

static PyObject *PyFits_get_config_item(PyFits *self, PyObject *args) {
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
	siril_debug_print("end of gfit operation\n");
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
	if (PyModule_AddObject(m, "fits", (PyObject *) &PyFitsType) < 0) {
		Py_DECREF(&PyFitsType);
		Py_DECREF(m);
		return NULL;
	}

	return m;
}

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

	PyGILState_Release(gstate);  // Release the GIL
	g_idle_add(script_widgets_idle, NULL);
	return retval;
}

// Function to run a Python script from memory
gboolean run_python_script_from_mem(gpointer p) {
	const char *script_contents = (const char*) p;
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();  // Acquire the GIL

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

	PyGILState_Release(gstate);  // Release the GIL
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
