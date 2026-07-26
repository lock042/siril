/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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

#ifndef SRC_CORE_SIRIL_LANGUAGE_H_
#define SRC_CORE_SIRIL_LANGUAGE_H_

#include <glib.h>

typedef struct {
	char *locale;
	char *language_name;
} parsed_code;


void siril_language_parser_init();
void language_init(const gchar *language);

/* Helpers used by gui/preferences.c to populate/read the language combo.
 * Not for use outside gui/. */
GHashTable  *siril_language_get_full_list(void);
gchar       *extract_locale_from_string(const gchar *str);
gint         locale_compare(gconstpointer a, gconstpointer b);

/* Defined in gui/preferences.c; declared here so gui/callbacks.c can call. */
void  siril_language_fill_combo(const gchar *language);
gchar *get_interface_language(void);

#endif /* SRC_CORE_SIRIL_LANGUAGE_H_ */
