/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#include <json-glib/json-glib.h>

void spcc_object_free(spcc_object *data, gboolean free_struct);
void osc_sensor_free(osc_sensor *data, gboolean free_struct);

static int load_spcc_object_from_file(const gchar *jsonFilePath, spcc_object *data, int index, gboolean from_osc_sensor) {
#ifndef HAVE_JSON_GLIB
	siril_log_color_message(_("json-glib was not found at build time, cannot proceed. Install and rebuild.\n"), "red");
	return 0;
#else

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
        fprintf(stderr, "Error loading SPCC JSON file: %s\n", error->message);
        g_error_free(error);
        g_object_unref(parser);
        return 0;
    }

    // Parse JSON data
    node = json_parser_get_root(parser);
	// Ensure the root is an array
    if (!JSON_NODE_HOLDS_ARRAY(node)) {
        fprintf(stderr, "Error: The JSON file should contain an array of objects.\n");
        g_object_unref(parser);
        return 0;
    }

    // Get the array of objects
    array = json_node_get_array(node);
    int num_objects = json_array_get_length(array);
	if (index > num_objects) {
        fprintf(stderr, "Error: index out of range.\n");
        g_object_unref(parser);
        return 0;
    }
    object = json_array_get_object_element(array, index);

    // Get values from JSON and store in the struct
	const gchar *typestring = json_object_get_string_member(object, "type");
	if (!strcmp(typestring, "MONO_SENSOR")) {
		data->type = 1;
	} else if (!strcmp(typestring, "OSC_SENSOR")) {
		if (!from_osc_sensor) {
			g_object_unref(parser);
			return 2; // OSC sensors are handled by a different routine
		} else {
			data->type = 2; // We have been called from the OSC handler
		}
	} else if (!strcmp(typestring, "MONO_FILTER")) {
		data->type = 3;
	} else if (!strcmp(typestring, "OSC_FILTER")) {
		data->type = 4;
	} else if (!strcmp(typestring, "OSC_LPF")) {
		data->type = 5;
	} else if (!strcmp(typestring, "WB_REF")) {
		data->type = 6;
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
    data->quality = json_object_get_int_member(object, "dataQualityMarker");
	if (!data->quality) {
		goto validation_error;
	}
    if (json_object_has_member(object, "channel")) {
		const gchar *channel_string = json_object_get_string_member(object, "channel");
		if (!strcmp(channel_string, "RED")) {
			data->channel = 0;
		} else if (!strcmp(channel_string, "GREEN")) {
			data->channel = 1;
		} else if (!strcmp(channel_string, "BLUE")) {
			data->channel = 2;
		} else {
			// "OTHER" may be used for filters that aren't clearly R, G or B
			// we don't care about them for SPCC though, so move on
			goto validation_error;
		}
	} else {
		data->channel = -1;
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
		goto validation_error;
	}

    // Cleanup
    g_object_unref(parser);
	return 1;

validation_error:
    g_object_unref(parser);
	g_free(data->name);
    g_free(data->manufacturer);
    free(data->x);
    free(data->y);
	return 0;
#endif
}

static int compare_spcc_chan(const void *a, const void *b) {
	spcc_object *obj_a = (spcc_object*) a;
	spcc_object *obj_b = (spcc_object*) b;
	return obj_a->n - obj_b->n;
}

static gboolean load_osc_sensor_from_file(const gchar *jsonFilePath, osc_sensor *data) {
	gboolean retbool = TRUE;

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
	}
	int chan[3];
	gboolean correct_channels = TRUE;
	for (int i = 0 ; i < 3 ; i++) {
		chan[i] = data->channel[i].channel;
		if (chan[i] != i) {
			correct_channels = FALSE;
		}
	}
	// Ensure the channels are in the correct order
	if (!correct_channels) {
		qsort(data->channel, 3, sizeof(spcc_object), compare_spcc_chan);
		for (int i = 0 ; i < 3 ; i++) {
			data->channel[i].index = i;
		}
	}

	return retbool;
}

static gboolean processJsonFile(const char *file_path) {
	int retval;
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
        fprintf(stderr, "Error loading SPCC JSON file: %s\n", error->message);
        g_error_free(error);
        g_object_unref(parser);
        return FALSE;
    }

    // Parse JSON data
    node = json_parser_get_root(parser);
	// Ensure the root is an array
    if (!JSON_NODE_HOLDS_ARRAY(node)) {
        fprintf(stderr, "Error: The JSON file should contain an array of objects.\n");
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
			if (data->type == 3 && (data->channel < 0 || data->channel > 2)) {
				spcc_object_free(data, TRUE);
				return 0;
			}
			siril_debug_print("Read JSON object: %s\n", data->name);
			// Place the data into the correct list based on its type
			switch (data->type) {
				case 1:
					com.spcc_data.mono_sensors = g_list_append(com.spcc_data.mono_sensors, data);
					break;
				case 2:
					siril_debug_print("Error, this should have been trapped and handled by load_osc_sensor_from_file!\n");
					break;
				case 3:
					com.spcc_data.mono_filters[data->channel] = g_list_append(com.spcc_data.mono_filters[data->channel], data);
					break;
				case 4:
					com.spcc_data.osc_filters = g_list_append(com.spcc_data.osc_filters, data);
					break;
				case 5:
					com.spcc_data.osc_lpf = g_list_append(com.spcc_data.osc_lpf, data);
					break;
				case 6:
					com.spcc_data.wb_ref = g_list_append(com.spcc_data.wb_ref, data);
					break;
				default:
					g_warning("Unknown type: %d", data->type);
					spcc_object_free(data, TRUE);
					return 0;
			}
		}
		else if (retval == 2) {
			osc_sensor *osc = g_new0(osc_sensor, 1);
			retval = load_osc_sensor_from_file(file_path, osc);
			if (retval) {
				siril_debug_print("Read JSON object: %s\n", osc->channel[0].model);
				com.spcc_data.osc_sensors = g_list_append(com.spcc_data.osc_sensors, osc);
				return retval;
			} else {
				siril_debug_print("Error reading JSON object in file %s\n", file_path);
				return retval;
			}
		}
	}
	return retval;
}

/*********************** PUBLIC FUNCTIONS ****************************/

// Call to free the members of a spcc_object
void spcc_object_free(spcc_object *data, gboolean free_struct) {
	if (!data)
		return;
	g_free(data->name);
    g_free(data->manufacturer);
	g_free(data->filepath);
	g_free(data->source);
    free(data->x);
    free(data->y);
	if (free_struct)
		g_free(data);
	return;
}

void osc_sensor_free(osc_sensor *data, gboolean free_struct) {
	if(!data)
		return;
	for (int i = 0 ; i < 3 ; i++) {
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
	point* aa = (point*) a;
	point* bb = (point*) b;
	if (aa->x == bb->x)
    return 0;
	else if(aa->x > bb->x)
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
#ifndef HAVE_JSON_GLIB
	siril_log_color_message(_("Siril was not built with json-glib, cannot proceed.\n"), "red");
	return FALSE;
#else
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
        fprintf(stderr, "Error loading SPCC JSON file: %s\n", error->message);
        g_error_free(error);
        g_object_unref(parser);
        return FALSE;
    }

    // Parse JSON data
    node = json_parser_get_root(parser);
	// Ensure the root is an array
    if (!JSON_NODE_HOLDS_ARRAY(node)) {
        fprintf(stderr, "Error: The JSON file should contain an array of objects.\n");
        g_object_unref(parser);
        return FALSE;
    }

    // Get the array of objects
    array = json_node_get_array(node);
    int num_objects = json_array_get_length(array);
	if (index > num_objects) {
        fprintf(stderr, "Error: index out of range.\n");
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
		siril_debug_print("Error in JSON: unrecognised wavelength unit\n");
		g_object_unref(parser);
		return FALSE;
	}
	point *pairs = (point*) malloc(data->n * sizeof(point));
    for (int i = 0; i < data->n; i++) {
		pairs[i].x = json_array_get_double_element(wavelengthArray, i) * scalefactor;
		pairs[i].y = json_array_get_double_element(valuesArray, i) / valuerange;
    }
    qsort(pairs, data->n, sizeof(point), compare_pair_x);
	data->n = remove_duplicate_x(pairs, data->n, data->filepath);
    data->x = (double *)malloc(data->n * sizeof(double));
    data->y = (double *)malloc(data->n * sizeof(double));
    for (int i = 0; i < data->n; i++) {
        data->x[i] = pairs[i].x;
        data->y[i] = pairs[i].y;
    }
    free(pairs);
	data->arrays_loaded = TRUE;

	// Cleanup
    g_object_unref(parser);
	return TRUE;
#endif
}

// Call once to populate com.spcc_data with metadata for all known spcc_objects
// This doesn't populate the arrays (which take up more memory): the arrays can
// be populated for a required spcc_object with load_spcc_object_arrays()
static void processDirectory(const gchar *directory_path) {
	if (!directory_path)
		return;
	GError *error;
    GDir *dir = g_dir_open(directory_path, 0, &error);

    if (dir == NULL) {
        g_warning("Unable to open directory: %s", directory_path);
        return;
    }

    const gchar *filename;

    while ((filename = g_dir_read_name(dir)) != NULL) {
        gchar *file_path = g_build_filename(directory_path, filename, NULL);

        if (g_file_test(file_path, G_FILE_TEST_IS_DIR)) {
            // If the current item is a directory, recursively process it (ignore ., .. and .git)
            if (g_strcmp0(filename, ".") != 0 && g_strcmp0(filename, "..") && g_strcmp0(filename, ".git") != 0) {
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

void load_all_spcc_metadata() {
    const gchar *path = siril_get_spcc_repo_path();
    processDirectory(path);
	com.spcc_data.wb_ref = g_list_sort(com.spcc_data.wb_ref, compare_spcc_object_names);
	com.spcc_data.osc_sensors = g_list_sort(com.spcc_data.osc_sensors, compare_osc_object_models);
	com.spcc_data.osc_lpf = g_list_sort(com.spcc_data.osc_lpf, compare_spcc_object_names);
	com.spcc_data.osc_filters = g_list_sort(com.spcc_data.osc_filters, compare_spcc_object_names);
	com.spcc_data.mono_sensors = g_list_sort(com.spcc_data.mono_sensors, compare_spcc_object_names);
	for (int i = 0 ; i < 3 ; i++)
		com.spcc_data.mono_filters[i] = g_list_sort(com.spcc_data.mono_filters[i], compare_spcc_object_names);
}