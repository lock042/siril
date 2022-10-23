/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2022 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_IO_PATH_PARSE_H_
#define SRC_IO_PATH_PARSE_H_

#include "core/siril.h"

typedef enum {
	PATH_PARSE_OK = 0,
	PATH_PARSE_HEADER_NULL = 1,
	PATH_PARSE_WRONG_CALL = 2,
	PATH_PARSE_WRONG_RESERVED_KEYWORD = 3,
	PATH_PARSE_LIBRARY_NOTDEFINED = 4,
	PATH_PARSE_KEY_NOT_FOUND = 10,
	PATH_PARSE_WRONG_FORMAT = 11,
	PATH_PARSE_UNSUPPORTED_FORMAT = 12,
	PATH_PARSE_WRONG_DATE = 13,
	PATH_PARSE_WRONG_WCS = 14
} path_parse_errors;

gchar *path_parse(fits *fit, gchar *expression, gboolean *success);

#endif /* SRC_IO_PATH_PARSE_H_ */