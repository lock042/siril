/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
 *
 *
 * FITS sequences are not a sequence of FITS files but a FITS file containing a
 * sequence. It simply has as many elements in the third dimension as the
 * number of images in the sequence multiplied by the number of channels per
 * image. Given its use of the third dimension, it's sometimes called FITS cube.
 */

#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "algos/photometric_cc.h"
#include "algos/spcc.h"
#include "yyjson.h"

// Uncomment the following line for verbose confirmation of loading each JSON object
// #define DEBUG_JSON

void spcc_object_free(spcc_object *data, gboolean free_struct);
void osc_sensor_free(osc_sensor *data, gboolean free_struct);
spcc_object* spcc_object_copy(spcc_object *data);
gboolean spcc_metadata_loaded = FALSE;

static int load_spcc_object_from_file(const gchar *jsonFilePath, spcc_object *data, int index, gboolean from_osc_sensor) {
	if (!jsonFilePath)
		return FALSE;

	// Ensure data is zero-filled to prevent any issues with freeing members at validation fail time
	memset(data, 0, sizeof(spcc_object));

	// Read the entire file
	yyjson_read_err err;
	yyjson_doc *doc = yyjson_read_file(jsonFilePath, 0, NULL, &err);
	if (!doc) {
		siril_log_color_message(_("Error loading SPCC JSON file %s: %s\n"), "red", jsonFilePath, err.msg);
		return 0;
	}

	// Get root array
	yyjson_val *root = yyjson_doc_get_root(doc);
	if (!yyjson_is_arr(root)) {
		siril_log_color_message(_("Error: The JSON file should contain an array of objects.\n"), "red");
		yyjson_doc_free(doc);
		return 0;
	}

	// Get the object at specified index
	size_t num_objects = yyjson_arr_size(root);
	if (index >= num_objects) {
		siril_debug_print("Error: index out of range.\n");
		yyjson_doc_free(doc);
		return 0;
	}
	yyjson_val *object = yyjson_arr_get(root, index);

	// Get type and validate
	const char *typestring = yyjson_get_str(yyjson_obj_get(object, "type"));
	if (!typestring) goto validation_error;

	if (!strcmp(typestring, "MONO_SENSOR")) {
		data->type = MONO_SENSORS;
	} else if (!strcmp(typestring, "OSC_SENSOR")) {
		if (!from_osc_sensor) {
			yyjson_doc_free(doc);
			return 2; // OSC sensors are handled by a different routine
		} else {
			data->type = OSC_SENSORS;
		}
	} else if (!strcmp(typestring, "MONO_FILTER")) {
		data->type = MONO_FILTERS;
	} else if (!strcmp(typestring, "OSC_FILTER")) {
		data->type = OSC_FILTERS;
	} else if (!strcmp(typestring, "OSC_LPF")) {
		data->type = OSC_LPFS;
	} else if (!strcmp(typestring, "WB_REF")) {
		data->type = WB_REFS;
	} else {
		goto validation_error;
	}

	// Get required string fields
	const char *model = yyjson_get_str(yyjson_obj_get(object, "model"));
	const char *name = yyjson_get_str(yyjson_obj_get(object, "name"));
	const char *manufacturer = yyjson_get_str(yyjson_obj_get(object, "manufacturer"));
	const char *source = yyjson_get_str(yyjson_obj_get(object, "dataSource"));

	if (!model || !name || !manufacturer || !source) {
		goto validation_error;
	}

	data->model = g_strdup(model);
	data->name = g_strdup(name);
	data->manufacturer = g_strdup(manufacturer);
	data->source = g_strdup(source);

	// Get optional string fields
	yyjson_val *comment = yyjson_obj_get(object, "comment");
	if (comment) {
		data->comment = g_strdup(yyjson_get_str(comment));
	}

	// Get optional boolean fields
	yyjson_val *is_dslr = yyjson_obj_get(object, "is_dslr");
	if (is_dslr) {
		data->is_dslr = yyjson_get_bool(is_dslr);
	}

	// Get required integer fields
	yyjson_val *quality = yyjson_obj_get(object, "dataQualityMarker");
	yyjson_val *version = yyjson_obj_get(object, "version");
	if (!quality || !version) {
		goto validation_error;
	}
	data->quality = yyjson_get_int(quality);
	data->version = yyjson_get_int(version);
	if (!data->quality || !data->version) {
		goto validation_error;
	}

	// Handle channel field if present
	yyjson_val *channel_val = yyjson_obj_get(object, "channel");
	if (channel_val) {
		const char *channel_string = yyjson_get_str(channel_val);
		if (!strcmp(channel_string, "RED")) {
			data->channel = SPCC_RED;
		} else if (!strcmp(channel_string, "GREEN")) {
			data->channel = SPCC_GREEN;
		} else if (!strcmp(channel_string, "BLUE")) {
			data->channel = SPCC_BLUE;
		} else if (!strcmp(channel_string, "RED GREEN") || !strcmp(channel_string, "GREEN RED")) {
			data->channel = SPCC_RED | SPCC_GREEN;
		} else if (!strcmp(channel_string, "GREEN BLUE") || !strcmp(channel_string, "BLUE GREEN")) {
			data->channel = SPCC_GREEN | SPCC_BLUE;
		} else if (!strcmp(channel_string, "RED BLUE") || !strcmp(channel_string, "BLUE RED")) {
			data->channel = SPCC_RED | SPCC_BLUE;
		} else if (!strcmp(channel_string, "ALL") || !strcmp(channel_string, "RED GREEN BLUE") ||
			!strcmp(channel_string, "BLUE GREEN RED")) {
			data->channel = SPCC_CLEAR;
			} else {
				data->channel = SPCC_INVIS;
			}
	}

	// Store filepath and index
	data->filepath = g_strdup(jsonFilePath);
	data->index = index;

	// Validate arrays
	yyjson_val *wavelengthObject = yyjson_obj_get(object, "wavelength");
	yyjson_val *wavelengthArray = yyjson_obj_get(wavelengthObject, "value");
	yyjson_val *valuesObject = yyjson_obj_get(object, "values");
	yyjson_val *valuesArray = yyjson_obj_get(valuesObject, "value");

	if (!yyjson_is_arr(wavelengthArray) || !yyjson_is_arr(valuesArray)) {
		goto validation_error;
	}

	data->n = yyjson_arr_size(wavelengthArray);
	size_t valuesLength = yyjson_arr_size(valuesArray);

	if (data->n != valuesLength) {
		siril_log_color_message(_("Error loading SPCC JSON file: arrays have not the same size (%zu != %zu)\n"),
								"red", data->n, valuesLength);
		goto validation_error;
	}

	yyjson_doc_free(doc);
	return 1;

	validation_error:
	yyjson_doc_free(doc);
	g_free(data->model);
	data->model = NULL;
	g_free(data->name);
	data->name = NULL;
	g_free(data->source);
	data->source = NULL;
	g_free(data->manufacturer);
	data->manufacturer = NULL;
	return 0;
}

static int compare_spcc_chan(const void *a, const void *b) {
	const spcc_object *obj_a = (spcc_object*) a;
	const spcc_object *obj_b = (spcc_object*) b;
	return obj_a->channel - obj_b->channel;
}

static gboolean load_osc_sensor_from_file(const gchar *jsonFilePath, osc_sensor *data) {
	gboolean retbool = TRUE;
	gboolean found_dslr = FALSE;

	if (!jsonFilePath)
		return FALSE;
	for (int i = 0 ; i < 3 ; i++) {
		if (load_spcc_object_from_file(jsonFilePath, &data->channel[i], i, TRUE) != 1) {
			retbool = FALSE;
			for (int j = 0 ; j <= i ; j++) {
				spcc_object_free(&data->channel[j], FALSE);
			}
			return FALSE;
		}
		if (data->channel[i].is_dslr) {
			found_dslr = TRUE;
		}
	}

	int chan[3];
	gboolean correct_channels = TRUE;
	// Check if channels are in the correct order and assign is_dslr consistently
	for (int i = 0 ; i < 3 ; i++) {
		chan[i] = data->channel[i].channel;
		if (chan[i] != i) {
			correct_channels = FALSE;
		}
		data->channel[i].is_dslr = found_dslr;
	}
	// Ensure the channels are in the correct order, if needed
	if (!correct_channels) {
		qsort(data->channel, 3, sizeof(spcc_object), compare_spcc_chan);
		for (int i = 0 ; i < 3 ; i++) {
			data->channel[i].index = i;
		}
	}

	return retbool;
}

static gboolean processJsonFile(const char *file_path) {
	int retval = 0;

	if (!file_path)
		return FALSE;

	// Read and parse JSON file
	yyjson_read_err err;
	yyjson_doc *doc = yyjson_read_file(file_path, 0, NULL, &err);
	if (!doc) {
		siril_log_color_message(_("Error loading SPCC JSON file %s: %s\n"), "red", file_path, err.msg);
		return FALSE;
	}

	// Get root array
	yyjson_val *root = yyjson_doc_get_root(doc);
	if (!yyjson_is_arr(root)) {
		siril_debug_print("Error: The JSON file should contain an array of objects.\n");
		yyjson_doc_free(doc);
		return FALSE;
	}

	// Get number of objects
	size_t num_objects = yyjson_arr_size(root);

	// Process each object in the array
	for (size_t index = 0; index < num_objects; index++) {
		spcc_object *data = g_new0(spcc_object, 1);
		// Ensure data is zero-filled to prevent any issues with freeing members at validation fail time
		memset(data, 0, sizeof(spcc_object));

		retval = load_spcc_object_from_file(file_path, data, (int)index, FALSE);

		if (retval == 1) {
			#ifdef DEBUG_JSON
			siril_debug_print("Read JSON object: %s\n", data->name);
			#endif
			// Place the data into the correct list based on its type
			switch (data->type) {
				case MONO_SENSORS:
					com.spcc_data.mono_sensors = g_list_append(com.spcc_data.mono_sensors, data);
					break;

				case OSC_SENSORS:
					siril_debug_print("Error, this should have been trapped and handled by load_osc_sensor_from_file!\n");
					break;

				case MONO_FILTERS:;
				gboolean added = FALSE;
				for (int chan = RLAYER; chan <= BLAYER; chan++) {
					if ((int)data->channel & (1 << chan)) {
						if (!added) {
							com.spcc_data.mono_filters[chan] = g_list_append(com.spcc_data.mono_filters[chan], data);
							added = TRUE;
						} else {
							// If a spcc object has already been added to one channel and needs adding
							// to another we MUST make a copy, otherwise we get corruption when reloading
							// the objects
							spcc_object *copy = spcc_object_copy(data);
							com.spcc_data.mono_filters[chan] = g_list_append(com.spcc_data.mono_filters[chan], copy);
						}
					}
				}
				if (!added) {
					spcc_object_free(data, TRUE); // Free if not added anywhere
				}
				break;

				case OSC_FILTERS:
					com.spcc_data.osc_filters = g_list_append(com.spcc_data.osc_filters, data);
					break;

				case OSC_LPFS:
					com.spcc_data.osc_lpf = g_list_append(com.spcc_data.osc_lpf, data);
					break;

				case WB_REFS:
					com.spcc_data.wb_ref = g_list_append(com.spcc_data.wb_ref, data);
					break;

				default:
					siril_debug_print("Unknown type: %d", data->type);
					spcc_object_free(data, TRUE);
					yyjson_doc_free(doc);
					return FALSE;
			}
		} else {
			spcc_object_free(data, TRUE);
			if (retval == 2) {
				osc_sensor *osc = g_new0(osc_sensor, 1);
				retval = load_osc_sensor_from_file(file_path, osc);
				if (retval) {
					#ifdef JSON_DEBUG
					siril_debug_print("Read JSON object: %s\n", osc->channel[0].model);
					#endif
					com.spcc_data.osc_sensors = g_list_append(com.spcc_data.osc_sensors, osc);
					yyjson_doc_free(doc);
					return TRUE;
				} else {
					siril_log_color_message(_("Error reading JSON object in file %s\n"), "red", file_path);
					osc_sensor_free(osc, TRUE);
					yyjson_doc_free(doc);
					return FALSE;
				}
			}
		}
	}

	yyjson_doc_free(doc);
	return TRUE;
}

/*********************** PUBLIC FUNCTIONS ****************************/

spcc_object* spcc_object_copy(spcc_object *data) {
	spcc_object *copy = malloc(sizeof(spcc_object));
	if (!copy) return NULL;
	memcpy(copy, data, sizeof(spcc_object));
	copy->model = g_strdup(data->model);
	copy->name = g_strdup(data->name);
	copy->filepath = g_strdup(data->filepath);
	copy->comment = g_strdup(data->comment);
	copy->manufacturer = g_strdup(data->manufacturer);
	copy->source = g_strdup(data->source);
	copy->arrays_loaded = FALSE;
	copy->x = NULL;
	copy->y = NULL;
	copy->n = 0;
	return copy;
}

// Call to free the members of a spcc_object
void spcc_object_free(spcc_object *data, gboolean free_struct) {
	if (!data)
		return;
	g_free(data->name);
	data->name = NULL;
	g_free(data->model);
	data->model = NULL;
	g_free(data->manufacturer);
	data->manufacturer = NULL;
	g_free(data->filepath);
	data->filepath = NULL;
	g_free(data->source);
	data->source = NULL;
	free(data->x);
	data->x = NULL;
	free(data->y);
	data->y = NULL;
	if (free_struct) {
		g_free(data);
		data = NULL;
	}
	return;
}

void osc_sensor_free(osc_sensor *data, gboolean free_struct) {
	if (!data)
		return;
	for (int i = 0; i < 3; i++) {
		g_free(data->channel[i].model);
		data->channel[i].model = NULL;
		g_free(data->channel[i].name);
		data->channel[i].name = NULL;
		g_free(data->channel[i].manufacturer);
		data->channel[i].manufacturer = NULL;
		g_free(data->channel[i].filepath);
		data->channel[i].filepath = NULL;
		g_free(data->channel[i].source);
		data->channel[i].source = NULL;
		free(data->channel[i].x);
		data->channel[i].x = NULL;
		free(data->channel[i].y);
		data->channel[i].y = NULL;
	}
	if (free_struct) {
		g_free(data);
		data = NULL;
	}
	return;
}

// Frees and nulls the arrays of a spcc_object. This should be called
// after use to release the memory used by the arrays, but without
// destroying the entire spcc_object
void spcc_object_free_arrays(spcc_object *data) {
	if (!data)
		return;

	free(data->x);
	data->x = NULL;
	free(data->y);
	data->y = NULL;
	data->arrays_loaded = FALSE;
}

static int compare_pair_x(const void *a, const void *b) {
	const point *aa = (point*) a;
	const point *bb = (point*) b;
	if (aa->x == bb->x)
		return 0;
	else if (aa->x > bb->x)
		return 1;
	else
		return -1;;
}

int remove_duplicate_x(point *points, int n, const gchar *filename) {
	int write_index = 0;
	gchar *basename = g_path_get_basename(filename);
	for (int i = 1; i < n; i++) {
		if (points[i].x != points[write_index].x) {
			write_index++;
			points[write_index] = points[i];
		} else {
			siril_log_color_message(_("Warning: Duplicate x value detected in JSON file %s: %.1f\n"), "salmon", basename, points[i].x);
		}
	}
	g_free(basename);
	return write_index + 1;
}

gboolean load_spcc_object_arrays(spcc_object *data) {
	if (!data || !data->filepath)
		return FALSE;
	if (data->arrays_loaded)
		return TRUE;

	// Clear existing arrays
	free(data->x);
	free(data->y);
	data->x = data->y = NULL;

	// Read JSON file
	yyjson_read_err err;
	yyjson_doc *doc = yyjson_read_file(data->filepath, 0, NULL, &err);
	if (!doc) {
		siril_log_color_message(_("Error loading SPCC JSON file %s: %s\n"), "red", data->filepath, err.msg);
		return FALSE;
	}

	// Get root and handle both array and single object cases
	yyjson_val *root = yyjson_doc_get_root(doc);
	yyjson_val *object = NULL;

	if (yyjson_is_arr(root)) {
		size_t num_objects = yyjson_arr_size(root);
		if (data->index >= num_objects) {
			if (num_objects == 0) {
				siril_log_color_message(_("Error: index % " G_GSIZE_FORMAT " out of range (no objects).\n"),
										"red", data->index);
			} else {
				siril_log_color_message(_("Error: index % " G_GSIZE_FORMAT " out of range (max: % " G_GSIZE_FORMAT ").\n"),
										"red", data->index, num_objects - 1);
			}
			goto error_cleanup;
		}
		object = yyjson_arr_get(root, data->index);
	} else if (yyjson_is_obj(root)) {
		if (data->index != 0) {
			siril_log_color_message(_("Error: index must be 0 for single object JSON.\n"), "red");
			goto error_cleanup;
		}
		object = root;
	} else {
		siril_log_color_message(_("Error: JSON root must be an array or object.\n"), "red");
		goto error_cleanup;
	}

	// Get wavelength and values objects
	yyjson_val *wavelengthObject = yyjson_obj_get(object, "wavelength");
	yyjson_val *valuesObject = yyjson_obj_get(object, "values");
	if (!yyjson_is_obj(wavelengthObject) || !yyjson_is_obj(valuesObject)) {
		siril_log_color_message(_("Error: Missing or invalid wavelength/values objects.\n"), "red");
		goto error_cleanup;
	}

	// Get wavelength units and determine scale factor
	yyjson_val *units = yyjson_obj_get(wavelengthObject, "units");
	if (!yyjson_is_str(units)) {
		siril_log_color_message(_("Error: Missing or invalid wavelength units.\n"), "red");
		goto error_cleanup;
	}

	double scalefactor;
	const char *wavelengthUnit = yyjson_get_str(units);
	if (strcmp(wavelengthUnit, "nm") == 0)        scalefactor = 1.0;
	else if (strcmp(wavelengthUnit, "micrometer") == 0) scalefactor = 1000.0;
	else if (strcmp(wavelengthUnit, "angstrom") == 0)   scalefactor = 0.1;
	else if (strcmp(wavelengthUnit, "m") == 0)          scalefactor = 1.0e9;
	else {
		siril_log_color_message(_("Warning: error in JSON file %s: unrecognised wavelength unit\n"),
								"salmon", data->filepath);
		goto error_cleanup;
	}

	// Get value range for normalization
	yyjson_val *valuerange_val = yyjson_obj_get(valuesObject, "range");
	if (!valuerange_val || !yyjson_is_num(valuerange_val)) {
		siril_log_color_message(_("Error: Missing or invalid value range.\n"), "red");
		goto error_cleanup;
	}
	double valuerange = yyjson_get_num(valuerange_val);  // handles both int and real

	// Get arrays and validate
	yyjson_val *wavelengthArray = yyjson_obj_get(wavelengthObject, "value");
	yyjson_val *valuesArray = yyjson_obj_get(valuesObject, "value");
	if (!yyjson_is_arr(wavelengthArray) || !yyjson_is_arr(valuesArray)) {
		siril_log_color_message(_("Error: Missing or invalid value arrays.\n"), "red");
		goto error_cleanup;
	}

	size_t wavelengthCount = yyjson_arr_size(wavelengthArray);
	if (wavelengthCount != yyjson_arr_size(valuesArray) ||
		wavelengthCount < 5 || wavelengthCount > 2000) {
		siril_log_color_message(_("Error: Invalid array sizes (wavelength: %zu, values: %zu).\n"),
								"red", wavelengthCount, yyjson_arr_size(valuesArray));
		goto error_cleanup;
		}

		// Set number of points to process
		data->n = (data->n == 0 || data->n > wavelengthCount) ? wavelengthCount : data->n;

	// Allocate temporary storage
	point *pairs = malloc(data->n * sizeof(point));
	if (!pairs) {
		siril_log_color_message(_("Error: Memory allocation failed.\n"), "red");
		goto error_cleanup;
	}

	// Read data points using array iteration
	size_t idx = 0;
	yyjson_val *w_val, *v_val;
	size_t w_idx, max;

	yyjson_arr_foreach(wavelengthArray, w_idx, max, w_val) {
		if (idx >= data->n) break;

		v_val = yyjson_arr_get(valuesArray, w_idx);
		if (!v_val || !yyjson_is_num(w_val) || !yyjson_is_num(v_val)) {
			siril_log_color_message(_("Error: Invalid number at index %zu.\n"), "red", idx);
			free(pairs);
			goto error_cleanup;
		}

		pairs[idx].x = yyjson_get_num(w_val) * scalefactor;
		pairs[idx].y = yyjson_get_num(v_val) / valuerange;
		idx++;
	}

	// Process and store the data
	qsort(pairs, data->n, sizeof(point), compare_pair_x);
	data->n = remove_duplicate_x(pairs, data->n, data->filepath);

	// Allocate final arrays
	data->x = malloc(data->n * sizeof(double));
	data->y = malloc(data->n * sizeof(double));
	if (!data->x || !data->y) {
		siril_log_color_message(_("Error: Memory allocation failed.\n"), "red");
		free(pairs);
		free(data->x);
		free(data->y);
		data->x = data->y = NULL;
		goto error_cleanup;
	}

	// Copy sorted data
	for (size_t i = 0; i < data->n; i++) {
		data->x[i] = pairs[i].x;
		data->y[i] = pairs[i].y;
	}

	// Handle WB_REFS type scaling
	if (data->type == WB_REFS) {
		size_t norm_ref = 0;
		while (norm_ref < data->n - 1 && data->x[norm_ref] < 550) {
			norm_ref++;
		}

		double norm = data->y[norm_ref];
		for (size_t i = 0; i < data->n; i++) {
			data->y[i] = (data->y[i] * data->x[i]) / norm;
		}
	}

	free(pairs);
	yyjson_doc_free(doc);
	data->arrays_loaded = TRUE;
	return TRUE;

	error_cleanup:
	yyjson_doc_free(doc);
	return FALSE;
}

// Call once to populate com.spcc_data with metadata for all known spcc_objects
// This doesn't populate the arrays (which take up more memory): the arrays can
// be populated for a required spcc_object with load_spcc_object_arrays()
static void processDirectory(const gchar *directory_path) {
	if (!directory_path)
		return;
	GError *error = NULL;
	GDir *dir = g_dir_open(directory_path, 0, &error);

	if (dir == NULL) {
		siril_debug_print("Unable to open directory: %s", directory_path);
		g_error_free(error);
		return;
	}

	const gchar *filename;

	while ((filename = g_dir_read_name(dir)) != NULL) {
		gchar *file_path = g_build_filename(directory_path, filename, NULL);

		if (g_file_test(file_path, G_FILE_TEST_IS_DIR)) {
			// If the current item is a directory, recursively process it (ignore ., .. and utils and .git)
			if (g_strcmp0(filename, ".") && g_strcmp0(filename, "..") && g_strcmp0(filename, ".git")) {
				processDirectory(file_path);
			}
		} else {
			// Check if the file has a .json extension
			if (g_str_has_suffix(filename, ".json") && !g_strrstr(filename, "schema")) {
				processJsonFile(file_path);
			}
		}

		g_free(file_path);
	}

	g_dir_close(dir);
}

gint compare_spcc_object_names(gconstpointer a, gconstpointer b) {
	const spcc_object *object1 = a;
	const spcc_object *object2 = b;
	return g_strcmp0(object1->name, object2->name);
}

gint compare_osc_object_models(gconstpointer a, gconstpointer b) {
	const osc_sensor *object1 = a;
	const osc_sensor *object2 = b;
	return g_strcmp0(object1->channel[0].model, object2->channel[0].model);
}

static void spcc_object_destroy(void *user_data) {
	spcc_object *object = (spcc_object*) user_data;
	spcc_object_free(object, TRUE);
}

static void osc_sensor_destroy(void *user_data) {
	osc_sensor *object = (osc_sensor*) user_data;
	osc_sensor_free(object, TRUE);
}

void load_all_spcc_metadata() {
	// Ensure any previous content in the lists is removed and freed properly
	g_list_free_full(com.spcc_data.wb_ref, (GDestroyNotify) spcc_object_destroy);
	com.spcc_data.wb_ref = NULL;
	g_list_free_full(com.spcc_data.osc_lpf, (GDestroyNotify)spcc_object_destroy);
	com.spcc_data.osc_lpf = NULL;
	g_list_free_full(com.spcc_data.osc_filters, (GDestroyNotify)spcc_object_destroy);
	com.spcc_data.osc_filters = NULL;
	g_list_free_full(com.spcc_data.mono_sensors, (GDestroyNotify)spcc_object_destroy);
	com.spcc_data.mono_sensors = NULL;
	g_list_free_full(com.spcc_data.osc_sensors, (GDestroyNotify)osc_sensor_destroy);
	com.spcc_data.osc_sensors = NULL;
	for (int i = 0 ; i < 3 ; i++) {
		g_list_free_full(com.spcc_data.mono_filters[i], (GDestroyNotify)spcc_object_destroy);
		com.spcc_data.mono_filters[i] = NULL;
	}

	const gchar *path = siril_get_spcc_repo_path();
	processDirectory(path);
	siril_debug_print("SPCC JSON metadata loaded\n");

	com.spcc_data.wb_ref = g_list_sort(com.spcc_data.wb_ref, compare_spcc_object_names);
	com.spcc_data.osc_sensors = g_list_sort(com.spcc_data.osc_sensors, compare_osc_object_models);
	com.spcc_data.osc_lpf = g_list_sort(com.spcc_data.osc_lpf, compare_spcc_object_names);
	com.spcc_data.osc_filters = g_list_sort(com.spcc_data.osc_filters, compare_spcc_object_names);
	com.spcc_data.mono_sensors = g_list_sort(com.spcc_data.mono_sensors, compare_spcc_object_names);
	for (int i = 0 ; i < 3 ; i++)
		com.spcc_data.mono_filters[i] = g_list_sort(com.spcc_data.mono_filters[i], compare_spcc_object_names);
	spcc_metadata_loaded = TRUE;
}

void load_spcc_metadata_if_needed() {
	if (!spcc_metadata_loaded)
		load_all_spcc_metadata();
}
