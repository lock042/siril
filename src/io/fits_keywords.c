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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"

#include "fits_keywords.h"

#define KEYWORD(group, key, type, comment, data) { group, key, type, comment, data, TRUE }

KeywordInfo *initialize_keywords(fits *fit) {
	KeywordInfo keyword_list[] = {
        KEYWORD( "image", "ROWORDER", KTYPE_STR, "Order of the rows in image array", &(fit->keywords.row_order)),
        KEYWORD( "setup", "INSTRUME", KTYPE_STR, "Instrument name", &(fit->keywords.instrume)),
        KEYWORD( "setup", "TELESCOP", KTYPE_STR, "Telescope used to acquire this image", &(fit->keywords.telescop)),
        KEYWORD( "setup", "OBSERVER", KTYPE_STR, "Observer name", &(fit->keywords.observer)),
        KEYWORD( "image", "BAYERPAT", KTYPE_STR, "Bayer color pattern", &(fit->keywords.bayer_pattern)),
        KEYWORD( "image", "XBAYROFF", KTYPE_INT, "X offset of Bayer array", &(fit->keywords.bayer_xoffset)),
        KEYWORD( "image", "YBAYROFF", KTYPE_INT, "Y offset of Bayer array", &(fit->keywords.bayer_yoffset)),
        KEYWORD( "image", "GAIN", KTYPE_USHORT, "Camera gain", &(fit->keywords.key_gain)),
        KEYWORD( "image", "OFFSET", KTYPE_USHORT, "Camera offset", &(fit->keywords.key_offset)),
        KEYWORD( "image", "CVF", KTYPE_DOUBLE, "Conversion factor (e-/adu)", &(fit->keywords.cvf)),
        KEYWORD( "date", "DATE-OBS", KTYPE_DATE, "YYYY-MM-DDThh:mm:ss observation start, UT", &(fit->keywords.date_obs)),
        KEYWORD( "image", "STACKCNT;NCOMBINE", KTYPE_UINT, "Stack frames", &(fit->keywords.stackcnt)),
        KEYWORD( "image", "EXPTIME;EXPOSURE", KTYPE_DOUBLE, "Exposure time [s]", &(fit->keywords.exposure)),
        KEYWORD( "image", "LIVETIME", KTYPE_DOUBLE, "Exposure time after deadtime correction", &(fit->keywords.livetime)),
        {NULL, NULL, KTYPE_BOOL, NULL, NULL, FALSE}
    };

	// Count the number of keywords in the list
	int num_keywords = 0;
	while (keyword_list[num_keywords].group != NULL) {
		num_keywords++;
	}

    // Allocate memory dynamically for the keyword array
	KeywordInfo *all_keywords = (KeywordInfo*) malloc((num_keywords + 1) * sizeof(KeywordInfo));

    // Copy keyword information from the list to the dynamic array
	for (int i = 0; i < num_keywords; i++) {
		all_keywords[i] = keyword_list[i];
	}

    // Mark the end of the list
	all_keywords[num_keywords] = keyword_list[num_keywords];

	return all_keywords;
}

int save_fits_keywords(fits *fit) {
	KeywordInfo *keys = initialize_keywords(fit);
	int status;
	gchar *str;
	GDateTime *date;
	while (keys->group) {
		if (!keys->is_used) {
			keys++;
			continue;
		}
		gchar** tokens = g_strsplit(keys->key, ";", -1);
		switch (keys->type) {
			case KTYPE_INT:
				status = 0;
				fits_update_key(fit->fptr, TINT, tokens[0], &(*((int*)keys->data)), keys->comment, &status);
				break;
			case KTYPE_UINT:
				status = 0;
				fits_update_key(fit->fptr, TUINT, tokens[0], &(*((guint*)keys->data)), keys->comment, &status);
				break;
			case KTYPE_USHORT:
				status = 0;
				fits_update_key(fit->fptr, TUSHORT, tokens[0], &(*((int*)keys->data)), keys->comment, &status);
				break;
			case KTYPE_DOUBLE:
				status = 0;
				fits_update_key(fit->fptr, TDOUBLE, tokens[0], &(*((double*)keys->data)), keys->comment, &status);
				break;
			case KTYPE_STR:
				status = 0;
				str = ((gchar*)keys->data);
				if (str && str[0] != '\0') {
					fits_update_key(fit->fptr, TSTRING, tokens[0], str, keys->comment, &status);
				}
				break;
			case KTYPE_DATE:
				status = 0;
				date = *((GDateTime**)keys->data);
				if (date) {
					if (!g_strcmp0("DATE", tokens[0])) {
						int itmp;
						char fit_date[40];
						fits_get_system_time(fit_date, &itmp, &status);
						fits_update_key(fit->fptr, TSTRING, tokens[0], fit_date, keys->comment, &status);
					} else {
						gchar *formatted_date = date_time_to_FITS_date(date);
						fits_update_key(fit->fptr, TSTRING, tokens[0], formatted_date, keys->comment, &status);
						g_free(formatted_date);
					}
				}
				break;
			default:
				siril_debug_print("Save_fits_keywords: Error. Type is not handled.\n");
		}
		keys++;
	}

	free(keys);

	return 0;
}
