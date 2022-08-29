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

#include "core/siril.h"
#include "core/proto.h"

#include "fits_keywords.h"

struct fits_keywords all_keywords[] = {
		{ "OBJECT",   KEY_STANDARD, TSTRING, "Name of the object of interest", "" },
		{ "INSTRUME", KEY_STANDARD, TSTRING, "Instrument name", "" },
		{ "TELESCOP", KEY_STANDARD, TSTRING, "Telescope used to acquire this image", "" },
		{ "OBSERVER", KEY_STANDARD, TSTRING, "Observer name", "" },
};

GHashTable *create_keywords() {
	return g_hash_table_newx(g_str_hash, g_str_equal);
}

void add_keyword(GHashTable *hash, const char *key, gpointer value) {
	g_hash_table_insert(hash, key, value);
}

static gpointer get_keyword(GHashTable *hash, const char *key) {
	return g_hash_table_lookup (hash, key);
}

char *get_keyword_key_value_str(GHashTable *hash, const char *key) {
	struct fits_keywords one_key = get_keyword(hash, key);

	g_assert(one_key.type == TSTRING);

	return one_key.str_value;
}

double get_keyword_key_value_double(GHashTable *hash, const char *key) {
	struct fits_keywords one_key = get_keyword(hash, key);

	g_assert(one_key.type == TDOUBLE);

	return one_key.value;
}
