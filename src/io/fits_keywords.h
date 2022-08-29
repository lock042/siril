/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
 */
#ifndef SRC_IO_FITS_KEYWORDS_H_
#define SRC_IO_FITS_KEYWORDS_H_

typedef struct kkey key;

enum key_type {
	KEY_MANDATORY,
	KEY_STANDARD,
	KEY_EXTENDED,
	KEY_COMMENT,
	KEY_HISTORY,
	KEY_SIRIL
};

struct fits_keywords {
	const char *key;	// key name
	enum key_type ktype; //
	int type; // type of the FITS keyword
	const char *desc;	// short description
	union {			// allowed values range
		char *str_value;
		double value;
	};
};

struct kkey {
	char object[FLEN_VALUE];		// OBJECT key
	char instrume[FLEN_VALUE];		// INSTRUME key
	char telescop[FLEN_VALUE];		// TELESCOP key
	char observer[FLEN_VALUE];		// OBSERVER key
};

void add_keyword(GHashTable *hash, const char *key, gpointer value);
char *get_keyword_key_value_str(GHashTable *hash, const char *key);
double get_keyword_key_value_double(GHashTable *hash, const char *key);

#endif /* SRC_IO_FITS_KEYWORDS_H_ */
