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
#include <wchar.h>
#include <datetime.h>

#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "python/PyFits_functions.h"
#include "python/PyFWHM_functions.h"
#include "python/PyHomography_functions.h"
#include "python/PyImgData_functions.h"
#include "python/PyImStats_functions.h"
#include "python/PyRegData_functions.h"
#include "python/PySeq_functions.h"
#include "python/siril_python_module_functions.h"
#include "gui/script_menu.h"
#include "core/siril_log.h"

// Helper functions

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
	{"type", (getter)PyFits_get_image_type, NULL, N_("Image datatype"), NULL},
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
	{"current", (PyCFunction)PyFits_gfit, METH_CLASS | METH_NOARGS, N_("Get the current Siril image as a PyFits object")},
	{"total", (PyCFunction)PyFits_get_total, METH_VARARGS, N_("Return the total pixel count for the specified channel.")},
	{"ngoodpix", (PyCFunction)PyFits_get_ngoodpix, METH_VARARGS, N_("Return the number of good pixels for the specified channel.")},
	{"mean", (PyCFunction)PyFits_get_mean, METH_VARARGS, N_("Return the mean pixel value for the specified channel.")},
	{"median", (PyCFunction)PyFits_get_median, METH_VARARGS, N_("Return the median pixel value for the specified channel.")},
	{"sigma", (PyCFunction)PyFits_get_sigma, METH_VARARGS, N_("Return the sigma (standard deviation) for the specified channel.")},
	{"avgdev", (PyCFunction)PyFits_get_avgdev, METH_VARARGS, N_("Return the average deviation for the specified channel.")},
	{"mad", (PyCFunction)PyFits_get_mad, METH_VARARGS, N_("Return the median absolute deviation (MAD) for the specified channel.")},
	{"sqrtbwmv", (PyCFunction)PyFits_get_sqrtbwmv, METH_VARARGS, N_("Return the square root of the biweight midvariance for the specified channel.")},
	{"location", (PyCFunction)PyFits_get_location, METH_VARARGS, N_("Return the location value for the specified channel.")},
	{"scale", (PyCFunction)PyFits_get_scale, METH_VARARGS, N_("Return the scale value for the specified channel.")},
	{"min", (PyCFunction)PyFits_get_min, METH_VARARGS, N_("Return the minimum pixel value for the specified channel.")},
	{"max", (PyCFunction)PyFits_get_max, METH_VARARGS, N_("Return the maximum pixel value for the specified channel.")},
	{"normvalue", (PyCFunction)PyFits_get_normvalue, METH_VARARGS, N_("Return the normalization value for the specified channel.")},
	{"bgnoise", (PyCFunction)PyFits_get_bgnoise, METH_VARARGS, N_("Return the background noise value for the specified channel.")},
    {"stats", (PyCFunction)PyFits_stats, METH_VARARGS, N_("Returns a copy of the statistics for a given channel in the FITS image as a PyImStats object.")},
    {"open", (PyCFunction)PyFits_open, METH_VARARGS | METH_CLASS, N_("Open a FITS file and load it into the PyFits object")},
	{"set_as_current", (PyCFunction)PyFits_move_to_gfit, METH_VARARGS | METH_KEYWORDS,
     N_("Move the PyFits data to the main Siril image. Takes an optional argument return_new: if return_new is True, returns a new PyFits object pointing to the main Siril image.")},
    {"save", (PyCFunction)PyFits_save, METH_VARARGS, N_("Save the current PyFits object to a FITS file")},
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

static PyGetSetDef PySeq_getsetters[] = {
	{"seqname", (getter)PySeq_get_seqname, NULL, N_("sequence name"), NULL},
	{"number", (getter)PySeq_get_number, NULL, N_("number of images"), NULL},
	{"selnum", (getter)PySeq_get_selnum, NULL, N_("number of selected images"), NULL},
	{"fixed", (getter)PySeq_get_fixed, NULL, N_("fixed length of image index in filename"), NULL},
	{"nb_layers", (getter)PySeq_get_nb_layers, NULL, N_("number of layers"), NULL},
	{"rx", (getter)PySeq_get_rx, NULL, N_("image width"), NULL},
	{"ry", (getter)PySeq_get_ry, NULL, N_("image height"), NULL},
	{"is_variable", (getter)PySeq_get_is_variable, NULL, N_("sequence has images of different sizes"), NULL},
	{"bitpix", (getter)PySeq_get_bitpix, NULL, N_("image pixel format"), NULL},
	{"reference_image", (getter)PySeq_get_reference_image, NULL, N_("reference image for registration"), NULL},
	{"type", (getter)PySeq_get_type, NULL, N_("sequence type"), NULL},
	{"current", (getter)PySeq_get_current, NULL, N_("current file number"), NULL},
	{"needs_saving", (getter)PySeq_get_needs_saving, NULL, N_("sequence needs saving"), NULL},
	{NULL}
};

PyMethodDef PySeq_methods[] = {
	{"current", (PyCFunction)PySeq_comseq, METH_CLASS | METH_NOARGS, N_("The current sequence loaded in Siril")},
	{"imstats", (PyCFunction)PySeq_get_imstats, METH_VARARGS, N_("Get stats for a specific frame and channel")},
	{"imgdata", (PyCFunction)PySeq_get_imgdata, METH_VARARGS, N_("Get imgparam for a specific frame")},
	{"regdata", (PyCFunction)PySeq_get_regdata, METH_VARARGS, N_("Get regparam for a specific frame and channel")},
	{"frame", (PyCFunction)PySeq_get_fits, METH_VARARGS, N_("Get the frame number specified by the parameter as a PyFits object")},
	{NULL}  /* Sentinel */
};

PyTypeObject PySeqType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.seq",
	.tp_doc = N_("Siril Sequence object"),
	.tp_basicsize = sizeof(PySeqObject),
	.tp_itemsize = 0,
	.tp_traverse = (traverseproc)PySeq_traverse,
	.tp_clear = (inquiry)PySeq_clear,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
	.tp_new = PySeq_new,
	.tp_init = (initproc)PySeq_init,
	.tp_dealloc = (destructor)PySeq_dealloc,
	.tp_getset = PySeq_getsetters,
	.tp_methods = PySeq_methods,
};

static PyGetSetDef PyImgData_getsetters[] = {
	{"filenum", (getter)PyImgData_get_filenum, NULL, N_("real file index in the sequence"), NULL},
	{"incl", (getter)PyImgData_get_incl, NULL, N_("indicator whether or not the frame is included for future processings"), NULL},
	{"date_obs", (getter)PyImgData_get_date_obs, NULL, N_("date of the observation"), NULL},
	{"airmass", (getter)PyImgData_get_airmass, NULL, N_("airmass of the image"), NULL},
	{"rx", (getter)PyImgData_get_rx, NULL, N_("image width"), NULL},
	{"ry", (getter)PyImgData_get_ry, NULL, N_("image height"), NULL},
	{NULL}  /* Sentinel */
};

PyTypeObject PyImgDataType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.imgdata",
	.tp_doc = N_("Siril ImgData object"),
	.tp_basicsize = sizeof(PyImgDataObject),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_new = PyType_GenericNew,
	.tp_getset = PyImgData_getsetters,
};

PyGetSetDef PyImStats_getsetters[] = {
	{"total", (getter)PyImStats_get_total, NULL, N_("total number of pixels"), NULL},
	{"ngoodpix", (getter)PyImStats_get_ngoodpix, NULL, N_("number of non-zero pixels"), NULL},
	{"mean", (getter)PyImStats_get_mean, NULL, N_("mean value"), NULL},
	{"median", (getter)PyImStats_get_median, NULL, N_("median value"), NULL},
	{"sigma", (getter)PyImStats_get_sigma, NULL, N_("sigma value"), NULL},
	{"avgDev", (getter)PyImStats_get_avgDev, NULL, N_("average deviation"), NULL},
	{"mad", (getter)PyImStats_get_mad, NULL, N_("median absolute deviation"), NULL},
	{"sqrtbwmv", (getter)PyImStats_get_sqrtbwmv, NULL, N_("square root of BWMV"), NULL},
	{"location", (getter)PyImStats_get_location, NULL, N_("location"), NULL},
	{"scale", (getter)PyImStats_get_scale, NULL, N_("scale"), NULL},
	{"min", (getter)PyImStats_get_min, NULL, N_("minimum value"), NULL},
	{"max", (getter)PyImStats_get_max, NULL, N_("maximum value"), NULL},
	{"normValue", (getter)PyImStats_get_normValue, NULL, N_("normalization value"), NULL},
	{"bgnoise", (getter)PyImStats_get_bgnoise, NULL, N_("background noise"), NULL},
	{NULL}
};

PyTypeObject PyImStatsType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.imstats",
	.tp_doc = N_("Siril ImStats object"),
	.tp_basicsize = sizeof(PyImStatsObject),
	.tp_itemsize = 0,
	.tp_traverse = (traverseproc)PyImStats_traverse,
	.tp_clear = (inquiry)PyImStats_clear,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
	.tp_new = PyType_GenericNew,
	.tp_getset = PyImStats_getsetters,
};

PyGetSetDef PyHomography_getsetters[] = {
	{"h00", (getter)PyHomography_get_h00, NULL, N_("h00 value"), NULL},
	{"h01", (getter)PyHomography_get_h01, NULL, N_("h01 value"), NULL},
	{"h02", (getter)PyHomography_get_h02, NULL, N_("h02 value"), NULL},
	{"h10", (getter)PyHomography_get_h10, NULL, N_("h10 value"), NULL},
	{"h11", (getter)PyHomography_get_h11, NULL, N_("h11 value"), NULL},
	{"h12", (getter)PyHomography_get_h12, NULL, N_("h12 value"), NULL},
	{"h20", (getter)PyHomography_get_h20, NULL, N_("h20 value"), NULL},
	{"h21", (getter)PyHomography_get_h21, NULL, N_("h21 value"), NULL},
	{"h22", (getter)PyHomography_get_h22, NULL, N_("h22 value"), NULL},
	{"pair_matched", (getter)PyHomography_get_pair_matched, NULL, N_("pairs matched"), NULL},
	{"inliers", (getter)PyHomography_get_Inliers, NULL, N_("inliers"), NULL},
	{NULL}
};

PyTypeObject PyHomographyType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.homography",
	.tp_doc = N_("Siril Homography object"),
	.tp_basicsize = sizeof(PyHomographyObject),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_new = PyHomography_new,
	.tp_init = (initproc)PyHomography_init,
	.tp_dealloc = (destructor)PyHomography_dealloc,
	.tp_as_buffer = &PyHomography_as_buffer,
	.tp_getset = PyHomography_getsetters,
};

PyGetSetDef PyRegData_getsetters[] = {
	{"fwhm", (getter)PyRegData_get_fwhm, NULL, N_("FWHM value"), NULL},
	{"wfwhm", (getter)PyRegData_get_wfwhm, NULL, N_("Weighted FWHM value"), NULL},
	{"roundness", (getter)PyRegData_get_roundness, NULL, N_("Roundness value"), NULL},
	{"quality", (getter)PyRegData_get_quality, NULL, N_("Quality value"), NULL},
	{"bg", (getter)PyRegData_get_bg, NULL, N_("Background level"), NULL},
	{NULL}
};

PyMethodDef PyRegData_methods[] = {
	{"get_homography", (PyCFunction)PyRegData_get_homography, METH_NOARGS, "Get Homography object"},
	{NULL}
};

PyTypeObject PyRegDataType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.regdata",
	.tp_doc = N_("Siril RegData object"),
	.tp_basicsize = sizeof(PyRegDataObject),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_new = PyType_GenericNew,
	.tp_getset = PyRegData_getsetters,
	.tp_methods = PyRegData_methods,
};

PyGetSetDef PyFWHM_getsetters[] = {
	{"star_name", (getter)PyFWHM_get_star_name, NULL, N_("starname"), NULL},
	{"R", (getter)PyFWHM_get_star_name, NULL, N_("optimized box size to enclose sufficient background pixels"), NULL},
	{"B", (getter)PyFWHM_get_star_name, NULL, N_("average sky background value"), NULL},
	{"A", (getter)PyFWHM_get_star_name, NULL, N_("amplitude"), NULL},
	{"x0", (getter)PyFWHM_get_star_name, NULL, N_("x coordinate of the PSF peak"), NULL},
	{"y0", (getter)PyFWHM_get_star_name, NULL, N_("y coordinate of the PSF peak"), NULL},
	{"sx", (getter)PyFWHM_get_star_name, NULL, N_("Size of the fitted function on the PSF x-axis"), NULL},
	{"sy", (getter)PyFWHM_get_star_name, NULL, N_("Size of the fitted function on the PSF y-axis"), NULL},
	{"fwhmx", (getter)PyFWHM_get_star_name, NULL, N_("FWHM along the PSF's x-axis in pixels"), NULL},
	{"fwhmy", (getter)PyFWHM_get_star_name, NULL, N_("FWHM along the PSF's y-axis in pixels"), NULL},
	{"fwhmx_arcsec", (getter)PyFWHM_get_star_name, NULL, N_("FWHM along the PSF x-axis in arcsec"), NULL},
	{"fwhmy_arcsec", (getter)PyFWHM_get_star_name, NULL, N_("FWHM along the PSF y-axis in arcsec"), NULL},
	{"angle", (getter)PyFWHM_get_star_name, NULL, N_("rotation angle of the PSF's x, y axes with respect to the image axes"), NULL},
	{"rmse", (getter)PyFWHM_get_star_name, NULL, N_("RMS error"), NULL},
	{"sat", (getter)PyFWHM_get_star_name, NULL, N_("level above which a star is considered saturated"), NULL},
	{"has_saturated", (getter)PyFWHM_get_star_name, NULL, N_("flags if the star is saturated"), NULL},
	{"beta", (getter)PyFWHM_get_star_name, NULL, N_("Moffat beta parameter"), NULL},
	{"xpos", (getter)PyFWHM_get_star_name, NULL, N_("x position of the star in the image"), NULL},
	{"ypos", (getter)PyFWHM_get_star_name, NULL, N_("y position of the star in the image"), NULL},
	{"mag", (getter)PyFWHM_get_star_name, NULL, N_("Magnitude"), NULL},
	{"Bmag", (getter)PyFWHM_get_star_name, NULL, N_("Blue filter magnitude"), NULL},
	{"SNR", (getter)PyFWHM_get_star_name, NULL, N_("Signal to Noise Ratio"), NULL},
	{"layer", (getter)PyFWHM_get_star_name, NULL, N_("Channel"), NULL},
	{"ra", (getter)PyFWHM_get_star_name, NULL, N_("Right Ascension"), NULL},
	{"dec", (getter)PyFWHM_get_star_name, NULL, N_("Declination"), NULL},
	{NULL}
};

PyTypeObject PyFWHMType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "siril.fwhm",
	.tp_doc = N_("Siril FWHM object"),
	.tp_basicsize = sizeof(PyFWHMObject),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_new = PyFWHM_new,
	.tp_init = (initproc)PyFWHM_init,
	.tp_dealloc = (destructor)PyFWHM_dealloc,
	.tp_getset = PyFWHM_getsetters,
};


// Function to strip control characters and "--More--" lines
char* strip_control_chars(const char* input) {
	size_t input_len = strlen(input);
	char* cleaned_output = malloc(input_len + 1);  // Allocate memory
	if (!cleaned_output) return NULL;  // Check for allocation failure

	size_t j = 0;
	int skip_line = 0;

	for (size_t i = 0; i < input_len; i++) {
		if (skip_line) {
			if (input[i] == '\n') {
				skip_line = 0;
			}
			continue;
		}

		// Check for "--More--" at the start of a line
		if (i == 0 || input[i-1] == '\n') {
			if (strncmp(input + i, "--More--", 8) == 0) {
				skip_line = 1;
				continue;
			}
		}

		if (input[i] == '\x08' && j > 0) {
			j--;  // Backtrack one character for backspace
		} else if (input[i] >= 32 || input[i] == '\n' || input[i] == '\t') {
			cleaned_output[j++] = input[i];
		}
	}

	cleaned_output[j] = '\0';
	return cleaned_output;
}

// Python method for write(), which logs output after stripping control chars
static PyObject* py_log_write(PyObject* self, PyObject* args) {
	const char* text;
	if (!PyArg_ParseTuple(args, "s", &text)) {
		return NULL;
	}
	char* cleaned_text = strip_control_chars(text);
	if (cleaned_text) {
		siril_log_message(cleaned_text);
		free(cleaned_text);  // Free the allocated memory
	}
	Py_RETURN_NONE;
}

// Python method for flush(), which does nothing but is required
static PyObject* py_log_flush(PyObject* self, PyObject* args) {
	Py_RETURN_NONE;
}

// Define the methods of the custom Python object
static PyMethodDef LogMethods[] = {
	{"write", py_log_write, METH_VARARGS, "Log output to custom logger"},
	{"flush", py_log_flush, METH_NOARGS, "Flush (no-op)"},
	{NULL, NULL, 0, NULL} // Sentinel
};

// Define the Python type for the log object
static PyTypeObject PyLog_Type = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "Log",
	.tp_basicsize = sizeof(PyObject),
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_methods = LogMethods,
};

// Function to set up custom pager with simulated input
static void setup_custom_pager() {
	PyRun_SimpleString(
		"import sys, pydoc\n"
		"class SimulatedInputPager:\n"
		"    def __init__(self, write_func):\n"
		"        self.write_func = write_func\n"
		"    def __call__(self, text):\n"
		"        lines = text.split('\\n')\n"
		"        for i, line in enumerate(lines):\n"
		"            self.write_func(line + '\\n')\n"
		"            if (i + 1) % 23 == 0:  # Simulate pressing space every 23 lines\n"
		"                self.write_func('--More--\\n')\n"
		"        self.write_func('\\n')  # Final newline\n"
		"pydoc.pager = SimulatedInputPager(sys.stdout.write)\n"
	);
}

// Initialize the custom log object and set up the environment
void init_custom_logger() {
	if (PyType_Ready(&PyLog_Type) < 0) return;

	PyObject* log_obj = PyObject_New(PyObject, &PyLog_Type);
	if (!log_obj) return;

	// Redirect sys.stdout and sys.stderr to the custom logger
	PySys_SetObject("stdout", log_obj);
	PySys_SetObject("stderr", log_obj);

	// Set up custom pager
	setup_custom_pager();

	Py_DECREF(log_obj);
}

// Define methods for the module
static PyMethodDef SirilMethods[] = {
	{"filename", (PyCFunction)siril_get_filename, METH_NOARGS, N_("Get the current image filename")},
	{"get_config_item", (PyCFunction)siril_get_config_item, METH_VARARGS, N_("Get a config item")},
	{"gui_block", py_gui_block, METH_NOARGS, N_("Block the GUI by disabling widgets except for Stop")},
	{"gui_unblock", py_gui_unblock, METH_NOARGS, N_("Unblock the GUI by enabling all widgets")},
	{"log", siril_log_message_wrapper, METH_VARARGS, N_("Log a message")},
	{"update_progress", siril_progress_wrapper, METH_VARARGS, N_("Update the progress bar")},
	{"notify_image_modified", (PyCFunction)siril_notify_gfit_modified, METH_NOARGS, N_("Notify that the main image has been modified.")},
	{"cmd", siril_processcommand, METH_VARARGS, N_("Execute a Siril command")},
	{"wd", (PyCFunction)siril_get_wd, METH_NOARGS, N_("Get the current working directory")},
	{"pipinstall", (PyCFunction)siril_pipinstall, METH_VARARGS, N_("Install a module using pip")},
	{"should_stop", (PyCFunction)siril_get_continue, METH_NOARGS, N_("Check if the Stop button has been pressed")},
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
	if (PyType_Ready(&PySeqType) < 0)
		return NULL;
	if (PyType_Ready(&PyImgDataType) < 0)
		return NULL;
	if (PyType_Ready(&PyImStatsType) < 0)
		return NULL;
	if (PyType_Ready(&PyRegDataType) < 0)
		return NULL;
	if (PyType_Ready(&PyHomographyType) < 0)
		return NULL;
	if (PyType_Ready(&PyFWHMType) < 0)
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

	Py_INCREF(&PySeqType);
	if (PyModule_AddObject(m, "seq", (PyObject *)&PySeqType) < 0) {
		Py_DECREF(&PySeqType);
		Py_DECREF(m);
		return NULL;
	}

	Py_INCREF(&PyImgDataType);
	if (PyModule_AddObject(m, "imgdata", (PyObject *)&PyImgDataType) < 0) {
		Py_DECREF(&PyImgDataType);
		Py_DECREF(m);
		return NULL;
	}

	Py_INCREF(&PyFitsType);
	if (PyModule_AddObject(m, "fwhm", (PyObject *) &PyFWHMType) < 0) {
		Py_DECREF(&PyFWHMType);
		Py_DECREF(m);
		return NULL;
	}
	Py_INCREF(&PyFitsType);
	if (PyModule_AddObject(m, "imstats", (PyObject *) &PyImStatsType) < 0) {
		Py_DECREF(&PyImStatsType);
		Py_DECREF(m);
		return NULL;
	}
	Py_INCREF(&PyFitsType);
	if (PyModule_AddObject(m, "regdata", (PyObject *) &PyRegDataType) < 0) {
		Py_DECREF(&PyRegDataType);
		Py_DECREF(m);
		return NULL;
	}
	Py_INCREF(&PyFitsType);
	if (PyModule_AddObject(m, "homography", (PyObject *) &PyHomographyType) < 0) {
		Py_DECREF(&PyHomographyType);
		Py_DECREF(m);
		return NULL;
	}
	return m;
}

// Functions to do with initializing and finalizing the interpreter

// Helper function to check if venv exists and is properly configured
static int check_venv(const char* venv_path) {
    char pip_path[1024];
    snprintf(pip_path, sizeof(pip_path), "%s/bin/pip", venv_path);

    FILE* f = fopen(pip_path, "r");
    if (f) {
        fclose(f);
        return 1;
    }
    return 0;
}

// Helper function to create and configure venv
static int setup_venv(PyConfig* config, const char* venv_path) {
    PyStatus status = Py_InitializeFromConfig(config);
    if (PyStatus_Exception(status)) {
        return 0;
    }

    // Import venv and create a new environment
    PyRun_SimpleString("import venv\n"
                      "import sys\n"
                      "import os\n"
                      "\n"
                      "def setup_venv(path):\n"
                      "    builder = venv.EnvBuilder(\n"
                      "        system_site_packages=False,\n"
                      "        clear=True,\n"
                      "        with_pip=True,\n"
                      "        upgrade_deps=True\n"
                      "    )\n"
                      "    builder.create(path)\n"
                      "\n"
                      "    # Ensure pip is up to date\n"
                      "    pip_path = os.path.join(path, 'bin', 'pip')\n"
                      "    if os.path.exists(pip_path):\n"
                      "        os.system(f'{pip_path} install --upgrade pip setuptools wheel')\n"
                      "\n"
                      "venv_path = sys.argv[1] if len(sys.argv) > 1 else None\n"
                      "if venv_path:\n"
                      "    setup_venv(venv_path)\n");

    // Clean up the temporary Python instance
    Py_Finalize();
    return 1;
}

// Function to initialize Python interpreter and load our module
gpointer init_python(void *user_data) {
	// Set up venv
	gchar* venv_path = g_build_filename(g_get_user_data_dir(), "siril", "venv", NULL);
    PyStatus status;

    // Check if venv exists and is properly configured
    if (!check_venv(venv_path)) {
        printf("Virtual environment not found or incomplete. Setting up...\n");

        // Set up initial config for venv creation
        PyConfig config_setup;
        PyConfig_InitPythonConfig(&config_setup);

        // Minimal configuration for venv setup
        config_setup.isolated = 0;  // Temporarily disable isolation for setup
        config_setup.use_environment = 1;
        config_setup.site_import = 1;

        // Configure argv for the venv setup script
        wchar_t venv_path_w[1024];
        mbstowcs(venv_path_w, venv_path, 1024);

        wchar_t* setup_argv[] = {
            L"python",
            venv_path_w
        };
        config_setup.argv.length = 2;
        config_setup.argv.items = setup_argv;

        if (!setup_venv(&config_setup, venv_path)) {
            fprintf(stderr, "Failed to create virtual environment\n");
            PyConfig_Clear(&config_setup);
            return NULL;
        }

        PyConfig_Clear(&config_setup);
    }

    // Convert venv path to wide string
    wchar_t venv_path_w[1024];
    mbstowcs(venv_path_w, venv_path, 1024);

    // First, set PYTHONHOME environment variable
    // This ensures Python can find its core modules during early initialization
    char pythonhome[1024];
    snprintf(pythonhome, sizeof(pythonhome), "PYTHONHOME=%s", venv_path);
    putenv(pythonhome);

    PyPreConfig preconfig;
    PyConfig config;

    PyPreConfig_InitPythonConfig(&preconfig);
    preconfig.isolated = 1;
    preconfig.use_environment = 1;  // Changed to 1 to use PYTHONHOME
    preconfig.configure_locale = 1;

    status = Py_PreInitialize(&preconfig);
    if (PyStatus_Exception(status)) {
        goto exception;
    }

    PyConfig_InitPythonConfig(&config);

    // Basic isolation settings
    config.isolated = 1;
    config.use_environment = 1;  // Changed to 1 to use PYTHONHOME
    config.parse_argv = 0;
//    config.safe_path = 1; // commented out as not available in Python-3.9 (appimage)
    config.user_site_directory = 0;
    config.write_bytecode = 0;
    config.site_import = 1;

    config.stdio_encoding = L"utf-8";
    config.stdio_errors = L"surrogateescape";

    // Set program name (important for initialization)
    status = PyConfig_SetString(&config, &config.program_name, venv_path_w);
    if (PyStatus_Exception(status)) {
        goto exception;
    }

    // Set executable to python in venv
    wchar_t executable_path[1024];
    swprintf(executable_path, 1024, L"%ls/bin/python", venv_path_w);
    status = PyConfig_SetString(&config, &config.executable, executable_path);
    if (PyStatus_Exception(status)) {
        goto exception;
    }

    // Let Python configure its own paths based on these settings
    status = PyConfig_Read(&config);
    if (PyStatus_Exception(status)) {
        goto exception;
    }

	PyImport_AppendInittab("siril", PyInit_siril);

    // Initialize Python with our configuration
    status = Py_InitializeFromConfig(&config);
    if (PyStatus_Exception(status)) {
        goto exception;
    }

    init_custom_logger();
	PyImport_ImportModule("siril");

    // Test the environment
    PyRun_SimpleString(
        "import sys\n"
        "import os\n"
        "print('Python version:', sys.version)\n"
        "print('Executable:', sys.executable)\n"
        "print('Prefix:', sys.prefix)\n"
        "print('Base Prefix:', sys.base_prefix)\n"
        "print('Python Path:', sys.path)\n"
        "\n"
        "try:\n"
        "    import encodings\n"
        "    print('Successfully imported encodings module')\n"
        "except ImportError as e:\n"
        "    print('Failed to import encodings:', e)\n"
        "\n"
        "try:\n"
        "    import pip\n"
        "    print('Pip version:', pip.__version__)\n"
        "except ImportError as e:\n"
        "    print('Warning: Could not import pip:', e)\n"
    );

	g_free(venv_path);
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

exception:

    if (PyStatus_Exception(status)) {
        Py_ExitStatusException(status);
    }
    PyConfig_Clear(&config);
	return NULL;
}

// Function to run a Python script from a file
gboolean run_python_script_from_file(gpointer p) {
	const char *script_path = (const char*) p; // must not be freed, it is owned by the list of script menu items
	PyGILState_STATE gstate;
	com.python_script = TRUE;
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
	com.python_script = FALSE;
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
		com.python_script = TRUE;
		PyObject *result = PyRun_String(script_contents, Py_file_input, globals, globals);
		com.python_script = FALSE;

		if (result == NULL) {
			// An exception occurred
			PyObject *type, *value, *traceback;
			PyErr_Fetch(&type, &value, &traceback);

			// Convert the error to a string
			PyObject *str_exc_value = PyObject_Str(value);
			const char* err_msg = PyUnicode_AsUTF8(str_exc_value);

			// Log the error
			siril_log_message(_("Error in Python script: %s\n"), err_msg);

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
	// Note: script_widgets_idle() call is omitted as per the original comment
	return retval;
}

// Function to run Python script, delegating to the Python thread if necessary
void run_python_script_in_python_thread(const char *script, gboolean from_file) {
	com.stop_script = FALSE;
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
