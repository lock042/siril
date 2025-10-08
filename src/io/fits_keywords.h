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
 */
#ifndef SRC_IO_FITS_KEYWORDS_H_
#define SRC_IO_FITS_KEYWORDS_H_

#include "core/siril.h"

enum keywords_type {
	KTYPE_BOOL,
	KTYPE_INT,
	KTYPE_UINT,
	KTYPE_USHORT,
	KTYPE_FLOAT,
	KTYPE_DOUBLE,
	KTYPE_STR,
	KTYPE_DATE
};

typedef struct KeywordInfo KeywordInfo;

typedef void (*special_handler_read_func)(fits *fit, const char *comment, KeywordInfo *info);
typedef void (*special_handler_save_func)(fits *fit,  KeywordInfo *info);

struct KeywordInfo {
    const char *group;    // group name
    const char *key;    // key name
    enum keywords_type type;    // type of the keyword
    const char *comment;    // comment
    void *data;    // pointer to the data in keyword struct
    special_handler_read_func special_handler_read;
    special_handler_save_func special_handler_save;
    gboolean is_saved;
    gboolean fixed_value;
    gboolean used;
};

struct keywords_data {
	fits *fit;
	sequence *seq;
	gchar *FITS_key;
	gchar *newkey;
	gchar *value;
	gchar *comment;
};

int save_fits_keywords(fits *fit);
int save_fits_unknown_keywords(fits *fit);
int save_history_keywords(fits *fit);
int save_wcs_keywords(fits *fit);
void save_gps_keywords(fits *fit);
void read_fits_date_obs_header(fits *fit);
int read_fits_keywords(fits *fit);
void read_qhy_gps_data(fits *fit);
int remove_all_fits_keywords(fits *fit);
void clear_Bayer_information(fits *fit);

void set_all_keywords_default(fits *fit);
void remove_keyword_in_fit_keywords(const gchar *keyword, fits *fit);

void start_sequence_keywords(sequence *seq, struct keywords_data *args);

int parse_wcs_image_dimensions(fits *fit, int *rx, int *ry);

gboolean keyword_is_protected(char *card, fits *fit);

#endif /* SRC_IO_FITS_KEYWORDS_H_ */
