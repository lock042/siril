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

#define KEYWORD(group, key, type, comment) { group, key, type, comment, NULL, TRUE }

static gboolean should_use_keyword(const fkeywords *keywords, const gchar *keyword) {
    if (g_strcmp0(keyword, "XBAYROFF") == 0) {
        return keywords->bayer_pattern[0] != '\0';
    } else if (g_strcmp0(keyword, "YBAYROFF") == 0) {
        return keywords->bayer_pattern[0] != '\0';
    }

    return TRUE;
}

struct keywords_access *get_all_keywords(fits *fit) {
    static struct keywords_access all_keywords[] = {
//    		KEYWORD( "image", "BZERO", KTYPE_DOUBLE, "Offset data range to that of unsigned short" ),
//    		KEYWORD( "image", "BSCALE", KTYPE_DOUBLE, "Default scaling factor" ),
    		KEYWORD( "image", "ROWORDER", KTYPE_STR, "Order of the rows in image array" ),
			KEYWORD( "setup", "INSTRUME", KTYPE_STR, "Instrument name" ),
			KEYWORD( "setup", "TELESCOP", KTYPE_STR, "Telescope used to acquire this image" ),
			KEYWORD( "setup", "OBSERVER", KTYPE_STR, "Observer name" ),
			KEYWORD( "image", "BAYERPAT", KTYPE_STR, "Bayer color pattern" ),
			KEYWORD( "image", "XBAYROFF", KTYPE_INT, "X offset of Bayer array" ),
			KEYWORD( "image", "YBAYROFF", KTYPE_INT, "Y offset of Bayer array" ),
			KEYWORD( "date",  "DATE-OBS", KTYPE_DATE, "YYYY-MM-DDThh:mm:ss observation start, UT" ),
			KEYWORD( "image", "STACKCNT;NCOMBINE", KTYPE_UINT, "Stack frames" ),
			KEYWORD( "image", "EXPTIME;EXPOSURE", KTYPE_DOUBLE, "Exposure time [s]" ),
			KEYWORD( "image", "LIVETIME", KTYPE_DOUBLE, "Exposure time after deadtime correction" ),
		    { NULL, NULL, KTYPE_BOOL, NULL, NULL }
    };

    /** Handle data **/
    int i = 0;
    all_keywords[i++].data = &fit->keywords.row_order;
    all_keywords[i++].data = &fit->keywords.instrume;
    all_keywords[i++].data = &fit->keywords.telescop;
    all_keywords[i++].data = &fit->keywords.observer;
    all_keywords[i++].data = &fit->keywords.bayer_pattern;
   	all_keywords[i++].data = &fit->keywords.bayer_xoffset;
   	all_keywords[i++].data = &fit->keywords.bayer_yoffset;
    all_keywords[i++].data = &fit->keywords.date_obs;
    all_keywords[i++].data = &fit->keywords.stackcnt;
    all_keywords[i++].data = &fit->keywords.exposure;
    all_keywords[i++].data = &fit->keywords.livetime;

    /** Handle conditions if needed */
    for (i = 0; all_keywords[i].group != NULL; i++) {
        all_keywords[i].is_used = should_use_keyword(&fit->keywords, all_keywords[i].key);
    }

    return all_keywords;
}

int save_fits_keywords(fits *fit) {
	struct keywords_access *keys = get_all_keywords(fit);
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
			case KTYPE_BOOL:

				break;
			case KTYPE_INT:
				status = 0;
				fits_update_key(fit->fptr, TINT, tokens[0], &(*((int*)keys->data)), keys->comment, &status);
				break;
			case KTYPE_UINT:
				status = 0;
				fits_update_key(fit->fptr, TUINT, tokens[0], &(*((guint*)keys->data)), keys->comment, &status);
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
					gchar *formatted_date = date_time_to_FITS_date(date);
					fits_update_key(fit->fptr, TSTRING, tokens[0], formatted_date, keys->comment, &status);
					g_free(formatted_date);
				}
				break;
		}
		keys++;
	}

	return 0;
}
