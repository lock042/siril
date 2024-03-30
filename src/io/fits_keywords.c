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

static fkeywords key = { 0 };

struct keywords_access all_keywords[] = {
	{ "image", "ROWORDER", KTYPE_STR, "Order of the rows in image array", &key.row_order },
	{ "setup", "INSTRUME", KTYPE_STR, "Instrument name", &key.instrume },
	{ "setup", "TELESCOP", KTYPE_STR, "Telescope used to acquire this image", &key.telescop },
	{ "date", "DATE-OBS", KTYPE_DATE, "YYYY-MM-DDThh:mm:ss observation start, UT", &key.date_obs },

	{ NULL, NULL, KTYPE_BOOL, NULL, NULL }
};

struct keywords_access *get_all_keywords(fits *fit) {
	memcpy(&key, &fit->keywords, sizeof(fkeywords));
	return all_keywords;
}


int save_fits_keywords(fits *fit) {
	struct keywords_access *keys = get_all_keywords(fit);
	int status;
	gchar *str;
	GDateTime *date;
	while (keys->group) {
		switch (keys->type) {
			case KTYPE_BOOL:

				break;
			case KTYPE_INT:

				break;
			case KTYPE_UINT:

				break;
			case KTYPE_DOUBLE:

				break;
			case KTYPE_STR:
				status = 0;
				str = ((gchar*)keys->data);
				if (str && str[0] != '\0')
					fits_update_key(fit->fptr, TSTRING, keys->key, str, keys->comment, &status);
				break;
			case KTYPE_DATE:
				status = 0;
				date = *((GDateTime**)keys->data);
				if (date) {
					gchar *formatted_date = date_time_to_FITS_date(date);
					fits_update_key(fit->fptr, TSTRING, keys->key, formatted_date, keys->comment, &status);
					g_free(formatted_date);
				}
				break;
		}
		keys++;
	}

	return 0;
}
