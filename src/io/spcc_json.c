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
#include <json-glib/json-glib.h>

void spcc_object_free(spcc_object *data, gboolean free_struct);

static gboolean load_spcc_object_from_file(const gchar *jsonFilePath, spcc_object *data) {
#ifndef HAVE_JSON_GLIB
	siril_log_color_message(_("json-glib was not found at build time, cannot proceed. Install and rebuild.\n"), "red");
	return FALSE;
#else
    GError *error = NULL;
    JsonParser *parser;
    JsonObject *root;
    JsonNode *node;

	// Ensure data is zero-filled to prevent any issues with freeing members at validation fail time
	memset(data, 0, sizeof(spcc_object));

    // Create a JSON parser
    parser = json_parser_new();

    // Load JSON file
    if (!json_parser_load_from_file(parser, jsonFilePath, &error)) {
        fprintf(stderr, "Error loading SPCC JSON file: %s\n", error->message);
        g_error_free(error);
        g_object_unref(parser);
        return FALSE;
    }

    // Parse JSON data
    node = json_parser_get_root(parser);
    root = json_node_get_object(node);

    // Get values from JSON and store in the struct
    data->name = g_strdup(json_object_get_string_member(root, "name"));
	if (!data->name) {
		goto validation_error;
	}
	const gchar *typestring = json_object_get_string_member(root, "type");
	if (!strcmp(typestring, "MONO_SENSOR")) {
		data->type = 1;
	} else if (!strcmp(typestring, "OSC_SENSOR")) {
		data->type = 2;
	} else if (!strcmp(typestring, "MONO_FILTER")) {
		data->type = 3;
	} else if (!strcmp(typestring, "OSC_FILTER")) {
		data->type = 4;
	} else {
		goto validation_error;
	}
    data->quality = json_object_get_int_member(root, "dataQualityMarker");
	if (!data->quality) {
		goto validation_error;
	}
	    data->manufacturer = g_strdup(json_object_get_string_member(root, "manufacturer"));
	if (!data->manufacturer) {
		goto validation_error;
	}
    data->version = json_object_get_int_member(root, "version");
	if (!data->type) {
		goto validation_error;
	}
	data->filepath = g_strdup(jsonFilePath);

    // Cleanup
    g_object_unref(parser);
	return TRUE;

validation_error:
    g_object_unref(parser);
	g_free(data->name);
    g_free(data->manufacturer);
    free(data->x);
    free(data->y);
	return FALSE;
#endif
}

static gboolean processJsonFile(const char *file_path) {
	int retval;
	spcc_object *data = g_new0(spcc_object, 1);
    retval = load_spcc_object_from_file(file_path, data);
   if (!retval) {
	   // Place the data into the correct list based on its type
		switch (data->type) {
			case 1:
				com.spcc_data.mono_sensors = g_list_append(com.spcc_data.mono_sensors, data);
				break;
			case 2:
				com.spcc_data.osc_sensors = g_list_append(com.spcc_data.osc_sensors, data);
				break;
			case 3:
				com.spcc_data.mono_filters = g_list_append(com.spcc_data.mono_filters, data);
				break;
			case 4:
				com.spcc_data.osc_filters = g_list_append(com.spcc_data.osc_filters, data);
				break;
			default:
				g_warning("Unknown type: %d", data->type);
				spcc_object_free(data, TRUE);
				break;
		}
   }
   return retval;
}

/*********************** PUBLIC FUNCTIONS ****************************/

// Call to free the members of a spcc_object
void spcc_object_free(spcc_object *data, gboolean free_struct) {
    g_free(data->name);
    g_free(data->manufacturer);
	g_free(data->filepath);
    free(data->x);
    free(data->y);
	if (free_struct)
		free(data);
	return;
}

// Frees and nulls the arrays of a spcc_object. This should be called
// after use to release the memory used by the arrays, but without
// destroying the entire spcc_object
void spcc_object_free_arrays(spcc_object *data) {
	free(data->x);
	data->x = NULL;
	free(data->y);
	data->y = NULL;
}

// Call to populate the arrays of a specific spcc_object
gboolean load_spcc_object_arrays(spcc_object *data) {
#ifndef HAVE_JSON_GLIB
	siril_log_color_message(_("json-glib was not found at build time, cannot proceed. Install and rebuild.\n"), "red");
	return FALSE;
#else
    GError *error = NULL;
    JsonParser *parser;
    JsonObject *root;
    JsonNode *node;

	// Ensure data is zero-filled to prevent any issues with freeing members at validation fail time
	memset(data, 0, sizeof(spcc_object));

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
    root = json_node_get_object(node);
    // Get 'wavelength' and 'values' arrays
    JsonArray *wavelengthArray = json_object_get_array_member(root, "wavelength");
    JsonArray *valuesArray = json_object_get_array_member(root, "values");
    data->n = json_array_get_length(wavelengthArray);
    int valuesLength = json_array_get_length(valuesArray);
	if (data->n != valuesLength) {
		goto validation_error;
	}
	data->x = (double *)malloc(data->n * sizeof(double));
    for (int i = 0; i < data->n; i++) {
        data->x[i] = json_array_get_double_element(wavelengthArray, i);
    }
    data->y = (double *)malloc(data->n * sizeof(double));
    for (int i = 0; i < data->n; i++) {
        data->y[i] = json_array_get_double_element(valuesArray, i);
    }
    // Cleanup
    g_object_unref(parser);
	return TRUE;

validation_error:
    g_object_unref(parser);
	// This function does not free the struct members: the caller must do this
	// and also remove the struct from its GList.
	return FALSE;
#endif
}

// Call once to populate com.spcc_data with metadata for all known spcc_objects
// This doesn't populate the arrays (which take up more memory): the arrays can
// be populated for a required spcc_object with load_spcc_object_arrays()
void load_all_spcc_metadata(gchar *path) {
    GDir *dir = g_dir_open(path, 0, NULL);

    if (dir == NULL) {
        g_warning("Unable to open directory: %s", path);
        return;
    }

    const gchar *filename;

    while ((filename = g_dir_read_name(dir)) != NULL) {
        gchar *file_path = g_build_filename(path, filename, NULL);

        // Check if the file has a .json extension
        if (g_str_has_suffix(filename, ".json")) {
            processJsonFile(file_path);
        }

        g_free(file_path);
    }

    g_dir_close(dir);
}
