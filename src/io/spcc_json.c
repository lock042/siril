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
#include <json-glib/json-glib.h>

// Uncomment the following line for verbose confirmation of loading each JSON object
// #define DEBUG_JSON

void spcc_object_free(spcc_object *data, gboolean free_struct);
void osc_sensor_free(osc_sensor *data, gboolean free_struct);
spcc_object* spcc_object_copy(spcc_object *data);
gboolean spcc_metadata_loaded = FALSE;

static int load_spcc_object_from_file(const gchar *jsonFilePath, spcc_object *data, int index, gboolean from_osc_sensor) {

	if (!jsonFilePath)
		return FALSE;

	GError *error = NULL;
	JsonParser *parser;
	JsonObject *object;
	JsonNode *node;
	JsonArray *array;

	// Ensure data is zero-filled to prevent any issues with freeing members at validation fail time
	memset(data, 0, sizeof(spcc_object));

	// Create a JSON parser
	parser = json_parser_new();

	// Load JSON file
	if (!json_parser_load_from_file(parser, jsonFilePath, &error)) {
		siril_log_color_message(_("Error loading SPCC JSON file: %s\n"), "red", error->message);
		g_error_free(error);
		g_object_unref(parser);
		return 0;
	}

	// Parse JSON data
	node = json_parser_get_root(parser);
	// Ensure the root is an array
	if (!JSON_NODE_HOLDS_ARRAY(node)) {
		siril_log_color_message(_("Error: The JSON file should contain an array of objects.\n"), "red");
		g_object_unref(parser);
		return 0;
	}

	// Get the array of objects
	array = json_node_get_array(node);
	int num_objects = json_array_get_length(array);
	if (index > num_objects) {
		siril_debug_print("Error: index out of range.\n");
		g_object_unref(parser);
		return 0;
	}
	object = json_array_get_object_element(array, index);

	// Get values from JSON and store in the struct
	const gchar *typestring = json_object_get_string_member(object, "type");
	if (!strcmp(typestring, "MONO_SENSOR")) {
		data->type = MONO_SENSORS;
	} else if (!strcmp(typestring, "OSC_SENSOR")) {
		if (!from_osc_sensor) {
			g_object_unref(parser);
			return 2; // OSC sensors are handled by a different routine
		} else {
			data->type = OSC_SENSORS; // We have been called from the OSC handler
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
	data->model = g_strdup(json_object_get_string_member(object, "model"));
	if (!data->model) {
		goto validation_error;
	}
	data->name = g_strdup(json_object_get_string_member(object, "name"));
	if (!data->name) {
		goto validation_error;
	}
	if (json_object_has_member(object, "comment")) {
		data->comment = g_strdup(json_object_get_string_member(object, "comment"));
	}
	if (json_object_has_member(object, "is_dslr")) {
		data->is_dslr = json_object_get_boolean_member(object, "is_dslr");
	}
	data->quality = json_object_get_int_member(object, "dataQualityMarker");
	if (!data->quality) {
		goto validation_error;
	}
	if (json_object_has_member(object, "channel")) {
		const gchar *channel_string = json_object_get_string_member(object, "channel");
		if (!strcmp(channel_string, "RED")) {
			data->channel = SPCC_RED;
		} else if (!strcmp(channel_string, "GREEN")) {
			data->channel = SPCC_GREEN;
		} else if (!strcmp(channel_string, "BLUE")) {
			data->channel = SPCC_BLUE;
		} else if (!strcmp(channel_string, "RED GREEN") || !strcmp(channel_string, "GREEN RED")) {
			data->channel = SPCC_RED | SPCC_BLUE;
		} else if (!strcmp(channel_string, "GREEN BLUE") || !strcmp(channel_string, "BLUE GREEN")) {
			data->channel = SPCC_GREEN | SPCC_BLUE;
		} else if (!strcmp(channel_string, "RED BLUE") || !strcmp(channel_string, "BLUE RED")) { // This is unlikely ever to be used
		// but is included for completeness
			data->channel = SPCC_RED | SPCC_BLUE;
		} else if (!strcmp(channel_string, "ALL") || !strcmp(channel_string, "RED GREEN BLUE") || !strcmp(channel_string, "BLUE GREEN RED")) {
			data->channel = SPCC_CLEAR;
		} else {
			data->channel = SPCC_INVIS; // Other filters e.g. UV or IR filters
			// should not be shown in SPCC combos
		}
	}
	data->manufacturer = g_strdup(json_object_get_string_member(object, "manufacturer"));
	if (!data->manufacturer) {
		goto validation_error;
	}
	data->source = g_strdup(json_object_get_string_member(object, "dataSource"));
	if (!data->source) {
		goto validation_error;
	}
	data->version = json_object_get_int_member(object, "version");
	if (!data->version) {
		goto validation_error;
	}
	data->filepath = g_strdup(jsonFilePath);
	data->index = index;

	// Get array lengths for validation
	JsonObject *wavelengthObject = json_object_get_object_member(object, "wavelength");
	JsonArray *wavelengthArray = json_object_get_array_member(wavelengthObject, "value");
	JsonObject *valuesObject = json_object_get_object_member(object, "values");
	JsonArray *valuesArray = json_object_get_array_member(valuesObject, "value");
	data->n = json_array_get_length(wavelengthArray);
	int valuesLength = json_array_get_length(valuesArray);
	if (data->n != valuesLength) {
		siril_log_color_message(_("Error loading SPCC JSON file: arrays have not the same size (%d != %d)\n"), "red", data->n, valuesLength);
		goto validation_error;
	}

    // Cleanup
	g_object_unref(parser);
	return 1;

validation_error:
	g_object_unref(parser);
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
	GError *error = NULL;
	JsonParser *parser;
	JsonNode *node;
	JsonArray *array;

	if (!file_path)
		return FALSE;

	// Create a JSON parser
	parser = json_parser_new();

	// Load JSON file
	if (!json_parser_load_from_file(parser, file_path, &error)) {
        siril_log_color_message(_("Error loading SPCC JSON file: %s\n"), "red", error->message);
        g_error_free(error);
        g_object_unref(parser);
		return FALSE;
	}

	// Parse JSON data
	node = json_parser_get_root(parser);
	// Ensure the root is an array
	if (!JSON_NODE_HOLDS_ARRAY(node)) {
		siril_debug_print("Error: The JSON file should contain an array of objects.\n");
		g_object_unref(parser);
		return FALSE;
	}

	// Get the array of objects
	array = json_node_get_array(node);
	int num_objects = json_array_get_length(array);

	g_object_unref(parser);

	for (int index = 0 ; index < num_objects ; index++) {
		spcc_object *data = g_new0(spcc_object, 1);
		// Ensure data is zero-filled to prevent any issues with freeing members at validation fail time
		memset(data, 0, sizeof(spcc_object));

		retval = load_spcc_object_from_file(file_path, data, index, FALSE);
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
					for (int chan = RLAYER ; chan <= BLAYER ; chan++) {
						if ((int) data->channel & (1 << chan)) {
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
					return TRUE;
				} else {
					siril_log_color_message(_("Error reading JSON object in file %s\n"), "red", file_path);
					osc_sensor_free(osc, TRUE);
					return FALSE;
				}
			}
		}
	}
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

// Call to populate the arrays of a specific spcc_object
gboolean load_spcc_object_arrays(spcc_object *data) {
	if (!data || !data->filepath) // Avoid dereferencing null pointers, if the spcc_object isn't prepopulated we can't load the arrays
		return FALSE;

	if (data->arrays_loaded)
		return TRUE;

	GError *error = NULL;
	JsonParser *parser;
	JsonObject *object;
	JsonNode *node;
	JsonArray *array;
	int index = data->index;

	// Ensure arrays are not already populated, to avoid memory leaks
	if (data->x) {
		free(data->x);
		data->x = NULL;
	}
	if (data->y) {
		free(data->y);
		data->y = NULL;
	}

	// Create a JSON parser
	parser = json_parser_new();

	// Load JSON file
	if (!json_parser_load_from_file(parser, data->filepath, &error)) {
		siril_log_color_message(_("Error loading SPCC JSON file: %s\n"), "red", error->message);
		g_error_free(error);
		g_object_unref(parser);
		return FALSE;
	}

	// Parse JSON data
	node = json_parser_get_root(parser);
	// Ensure the root is an array
	if (!JSON_NODE_HOLDS_ARRAY(node)) {
		siril_log_color_message(_("Error: The JSON file should contain an array of objects.\n"), "red");
		g_object_unref(parser);
		return FALSE;
	}

	// Get the array of objects
	array = json_node_get_array(node);
	int num_objects = json_array_get_length(array);
	if (index > num_objects) {
		siril_debug_print("Error: index out of range.\n");
		g_object_unref(parser);
		return FALSE;
	}

	object = json_array_get_object_element(array, index);

    // Get 'wavelength' and 'values' arrays
	double scalefactor = 1.0;
	JsonObject *wavelengthObject = json_object_get_object_member(object, "wavelength");
	JsonArray *wavelengthArray = json_object_get_array_member(wavelengthObject, "value");
	JsonObject *valuesObject = json_object_get_object_member(object, "values");
	JsonArray *valuesArray = json_object_get_array_member(valuesObject, "value");
	gchar *wavelengthUnit = g_strdup(json_object_get_string_member(wavelengthObject, "units"));
	double valuerange = json_object_get_double_member(valuesObject, "range");
	if (!strcmp(wavelengthUnit, "nm"))
		scalefactor = 1.0;
	else if (!strcmp(wavelengthUnit, "micrometer"))
		scalefactor = 1000.0;
	else if (!strcmp(wavelengthUnit, "angstrom"))
		scalefactor = 0.1;
	else if (!strcmp(wavelengthUnit, "m"))
		scalefactor = 1.0e9;
	else {
		siril_log_color_message(_("Warning: error in JSON file %s: unrecognised wavelength unit\n"), "salmon", data->filepath);
		g_object_unref(parser);
		g_free(wavelengthUnit);
		return FALSE;
	}
	point *pairs = (point*) malloc(data->n * sizeof(point));
	for (int i = 0; i < data->n; i++) {
		pairs[i].x = json_array_get_double_element(wavelengthArray, i) * scalefactor;
		pairs[i].y = json_array_get_double_element(valuesArray, i) / valuerange;
	}
	qsort(pairs, data->n, sizeof(point), compare_pair_x);
	data->n = remove_duplicate_x(pairs, data->n, data->filepath);
	data->x = (double*) malloc(data->n * sizeof(double));
	data->y = (double*) malloc(data->n * sizeof(double));
	for (int i = 0; i < data->n; i++) {
		data->x[i] = pairs[i].x;
		data->y[i] = pairs[i].y;
	}
	if (data->type == WB_REFS) {
		int norm_ref = 0;
		while (data->x[norm_ref] < 550) {
			if (norm_ref == data->n - 1) {
				norm_ref = 0;
				break;
			}
			norm_ref++;
		}
		for (int i = 0; i < data->n; i++) {
			data->y[i] *= data->x[i];
		}
		double norm = data->y[norm_ref];
		for (int i = 0; i < data->n; i++) {
			data->y[i] /= norm;
		}
	}
	free(pairs);
	data->arrays_loaded = TRUE;

	// Cleanup
	g_object_unref(parser);
	g_free(wavelengthUnit);
	return TRUE;
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
